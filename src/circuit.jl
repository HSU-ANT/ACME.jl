# Copyright 2015, 2016, 2017, 2018, 2019 Martin Holters
# See accompanying license file.

export Circuit, add!, connect!, disconnect!, @circuit, composite_element

import Base: delete!

struct CircuitNLFunc{Fs}
    fs::Fs
end

@inline (cnlf::CircuitNLFunc{Fs})(q) where {Fs} = _apply_all(q, cnlf.fs...)

_apply_all(q) = (SVector{0,Float64}(), SMatrix{0,0,Float64}())
@inline function _apply_all(q, f1::F, fs...) where F
    (res1, J1) = f1(q)
    (resrem, Jrem) = _apply_all(q, fs...)
    return ([res1; resrem], dcat(J1, Jrem))
end

dcat(a::SMatrix{Ma,Na,Ta}, b::SMatrix{Mb,Nb,Tb}) where {Ma,Na,Mb,Nb,Ta,Tb} =
    [[a zeros(SMatrix{Ma,Nb,Ta})]; [zeros(SMatrix{Mb,Na,Tb}) b]]

const Net = Vector{Tuple{Symbol,Symbol}} # pairs of element designator and pin name

struct Circuit
    elements::OrderedDict{Symbol, Element}
    nets :: Vector{Net}
    net_names :: Dict{Symbol, Net}
    Circuit() = new(OrderedDict{Symbol, Element}(), [], Dict{Symbol, Net}())
end

elements(c::Circuit) = values(c.elements)

for n in [:nb; :nx; :nq; :nu; :nl; :ny; :nn]
    @eval ($n)(c::Circuit) = reduce(+, init=0, $n(elem) for elem in elements(c))
end

for mat in [:mv; :mi; :mx; :mxd; :mq; :mu; :pv; :pi; :px; :pxd; :pq]
    # blockdiag() does not work, so include an empty matrix of desired type in
    # case elements(c) is empty
    # as blockdiag for unknown number of arguments cannot be inferred properly,
    # add type-assertion
    @eval ($mat)(c::Circuit) =
         blockdiag(spzeros(Rational{BigInt}, 0, 0),
                 (convert(SparseMatrixCSC{Rational{BigInt}}, elem.$mat)
                  for elem in elements(c))...
                )::SparseMatrixCSC{Rational{BigInt},Int}
end

u0(c::Circuit) = vcat((elem.u0 for elem in elements(c))...)

function incidence(c::Circuit)
    i = sizehint!(Int[], 2nb(c))
    j = sizehint!(Int[], 2nb(c))
    v = sizehint!(Int[], 2nb(c))
    for (row, pins) in enumerate(c.nets), (elemname, pinname) in pins
        elem = c.elements[elemname]
        offset = branch_offset(c, elem)
        bps = elem.pins[pinname]
        for (branch, polarity) in bps
            push!(i, row)
            push!(j, offset + branch)
            push!(v, polarity)
        end
    end
    # ensure zeros due to short-circuited branches are removed with dropzeros!
    dropzeros!(sparse(i, j, v, length(c.nets), nb(c)))
end

function nonlinear_eq_func(c::Circuit, elem_idxs=1:length(elements(c)))
    row_offset = 0
    col_offset = 0
    funcs = Function[]
    for elem in collect(elements(c))[elem_idxs]
        if nn(elem) == 0 && nq(elem) == 0
            continue
        end
        push!(funcs,
            let q_indices=SVector{nq(elem),Int}(col_offset+1:col_offset+nq(elem)),
                nleqfunc=nonlinear_eq_func(elem)
                @inline function (q)
                    nleqfunc(q[q_indices])
                end
            end)

        row_offset += nn(elem)
        col_offset += nq(elem)
    end
    return CircuitNLFunc((funcs...,))
end

"""
    add!(c::Circuit, elem::Element)

Adds the element `elem` to the circuit `c`, creating and returning a new, unique
reference designator, leaving its pins unconnected.
"""
function add!(c::Circuit, elem::Element)
    for (k, v) in c.elements
        if v == elem
            return k
        end
    end
    designator = gensym()
    add!(c, designator, elem)
    return designator
end

add!(c::Circuit, elems::Element...) = ((add!(c, elem) for elem in elems)...,)

"""
    add!(c::Circuit, designator::Symbol, elem::Element)

Adds the element `elem` to the circuit `c` with the reference designator
`designator`, leaving its pins unconnected. If the circuit already contained
an element named `designator`, it is removed first.
"""
function add!(c::Circuit, designator::Symbol, elem::Element)
    if haskey(c.elements, designator)
        delete!(c, designator)
    end
    for pin in keys(elem.pins)
        push!(c.nets, [(designator, pin)])
    end
    c.elements[designator] = elem
end

"""
    delete!(c::Circuit, designator::Symbol)

Deletes the element named `designator` from the circuit `c` (disconnecting all
its pins).
"""
function delete!(c::Circuit, designator::Symbol)
    for net in c.nets
        filter!(elempin -> elempin[1] != designator, net)
    end
    delete!(c.elements, designator)
end

function branch_offset(c::Circuit, elem::Element)
    offset = 0
    for el in elements(c)
        el == elem && return offset
        offset += nb(el)
    end
    throw(ArgumentError("Element not found in circuit"))
end

function netfor!(c::Circuit, p::Tuple{Symbol,Symbol})
    for net in c.nets
        p ∈ net && return net
    end
    throw(ArgumentError("Unknown pin $p"))
end
netfor!(c::Circuit, p::Tuple{Symbol,Any}) = netfor!(c, (p[1], Symbol(p[2])))

function netfor!(c::Circuit, name::Symbol)
    haskey(c.net_names, name) || push!(c.nets, get!(c.net_names, name, []))
    c.net_names[name]
end

"""
    connect!(c::Circuit, pins::Union{Symbol,Tuple{Symbol,Any}}...)

Connects the given pins (or named nets) to each other in the circuit `c`. Named
nets are given as `Symbol`s, pins are given as `Tuple{Symbols,Any}`s, where the
first entry is the reference designator of an element in `c`, and the second
entry is the pin name. For convenience, the latter is automatically converted to
a `Symbol` as needed.

# Example

```jldoctest; output = false, setup = :(using ACME)
circ = Circuit()
add!(circ, :r, resistor(1e3))
add!(circ, :src, voltagesource(5))
connect!(circ, (:src, -), (:r, 2), :gnd) # connect to gnd net

# output

```
"""
function connect!(c::Circuit, pins::Union{Symbol,Tuple{Symbol,Any}}...)
    nets = unique([netfor!(c, pin) for pin in pins])
    for net in nets[2:end]
        append!(nets[1], net)
        deleteat!(c.nets, findfirst(n -> n == net, c.nets))
        for (name, named_net) in c.net_names
            if named_net == net
                c.net_names[name] = nets[1]
            end
        end
    end
end

function disconnect!(c::Circuit, pin::Tuple{Symbol,Symbol})
    net = netfor!(c, pin)
    filter!(p -> p != pin, net)
    push!(c.nets, [pin])
end
"""
    disconnect!(c::Circuit, p::Tuple{Symbol,Any})

Disconnects the given pin `p` from anything else in the circuit `c`. The pin is
given as a`Tuple{Symbols,Any}`, where the first entry is the reference
designator of an element in `c`, and the second entry is the pin name. For
convenience, the latter is automatically converted to a `Symbol` as needed. Note
that if e.g. three pin `p1`, `p2`, and `p3` are connected then
`disconnect!(c, p1)` will disconnect `p1` from `p2` and `p3`, but leave `p2` and
`p3` connected to each other.
"""
disconnect!(c::Circuit, p::Tuple{Symbol,Any}) = disconnect!(c, (p[1], Symbol(p[2])))

function topomat!(incidence::SparseMatrixCSC{T}) where {T<:Integer}
    @assert all(x -> abs(x) == 1, nonzeros(incidence))
    @assert all(sum(incidence, dims=1) .== 0)

    t = falses(size(incidence)[2]);

    row = 1;
    for col = 1:size(incidence)[2]
        rows = filter(r -> r ≥ row, findall(!iszero, incidence[:, col]))
        @assert length(rows) ≤ 2

        isempty(rows) && continue
        t[col] = true;

        if rows[1] ≠ row
            incidence[[rows[1], row], :] = incidence[[row, rows[1]], :]
        end
        if length(rows) == 2
            @assert incidence[row, col] + incidence[rows[2], col] == 0
            incidence[rows[2],:] = incidence[rows[2],:] + incidence[row,:]
        end
        if incidence[row, col] < 0
            cols = findall(!iszero, incidence[row, :])
            incidence[row,cols] = -incidence[row,cols]
        end
        rows = findall(==(1), incidence[1:row-1, col])
        incidence[rows, :] .-= incidence[row, :]'
        rows = findall(==(-1), incidence[1:row-1, col])
        incidence[rows, :] .+= incidence[row, :]'
        row += 1
    end

    ti = incidence[1:row-1, :]

    dl = ti[:, (!).(t)]
    tv = spzeros(T, size(dl, 2), size(incidence, 2))
    # findall here works around JuliaLang/julia#27013 (introdcued in 0.7.0-DEV.4983)
    tv[:, findall(t)] = -dl'
    tv[:, findall((!).(t))] = SparseMatrixCSC{T}(I, size(dl, 2), size(dl, 2))

    tv, ti
end

topomat(incidence::SparseMatrixCSC{<:Integer}) = topomat!(copy(incidence))
topomat(c::Circuit) = topomat!(incidence(c))

@doc doc"""
    @circuit begin #= ... =# end

Provides a simple domain-specific language to decribe circuits. The
`begin`/`end` block can hold element definitions of the form
`refdes = elementfunc(params)` and connection specifications of the form
`refdes1[pin1] ⟷ refdes2[pin2]`.

# Example

To create a circuit with a voltage source connected to a resistor:
```jldoctest; output = false, setup = :(using ACME), filter = r"(ACME\.)?Circuit\(.*"s
@circuit begin
    src = voltagesource(5)
    r = resistor(1000)
    src[+] ⟷ r[1]
    src[-] ⟷ r[2]
end

# output

Circuit(...)
```

Alternatively, connection specifications can be given after an element
specification, separated by commas. In that case, the `refdes` may be omitted,
defaulting to the current element.

# Example

```jldoctest; output = false, setup = :(using ACME), filter = r"(ACME\.)?Circuit\(.*"s
@circuit begin
    src = voltagesource(5)
    r = resistor(1000), src[+] ⟷ [1], src[-] ⟷ [2]
end

# output

Circuit(...)
```

Finally, a connection endpoint may simply be of the form `netname`, to connect
to a named net. (Such named nets are created as needed.)

# Example

```jldoctest; output = false, setup = :(using ACME), filter = r"(ACME\.)?Circuit\(.*"s
@circuit begin
    src = voltagesource(5), [-] ⟷ gnd
    r = resistor(1000), [1] ⟷ src[+], [2] ⟷ gnd
end

# output

Circuit(...)
```

If a net or pin specification is not just a single symbol or number, and has to
be put in quotes (e.g. `"in+"`, `"9V"`)

!!! note
    Instead of `⟷` (`\\longleftrightarrow`), one can also use `==`.
""" macro circuit(cdef)
    is_conn_spec(expr::Expr) =
        (expr.head === :call && (expr.args[1] === :(⟷) || expr.args[1] === :(↔) || expr.args[1] === :(==))) ||
        (expr.head === :comparison && all(c -> c === :(==), expr.args[2:2:end]))
    is_conn_spec(::Any) = false

    function elem_spec(expr)
        if !isa(expr, Expr) || expr.head !== :(=)
            error("invalid element specification$locinfo: $(expr)")
        end
        if !isa(expr.args[1], Symbol)
            error("invalid element identifier$locinfo: $(expr.args[1])")
        end
        if isa(expr.args[2], Expr) && expr.args[2].head === :tuple
            if isempty(expr.args[2].args)
                error("invalid element specification$locinfo: $(expr.args[2])")
            end
            elemspec = expr.args[2].args[1]
            conn_exprs = expr.args[2].args[2:end]
        else
            elemspec = expr.args[2]
            conn_exprs = []
        end
        push!(ccode.args, :(add!(circ, $(QuoteNode(expr.args[1])), $(esc(elemspec)))))
        for conn_expr in conn_exprs
            if !is_conn_spec(conn_expr)
                error("invalid connection specification$locinfo: $conn_expr")
            end
            push!(ccode.args, Expr(:call, :connect!, :circ, extractpins(conn_expr, expr.args[1])...))
        end
    end

    function extractpins(expr::Expr, default_element=nothing)
        if expr.head === :call && (expr.args[1] === :(⟷) || expr.args[1] === :(↔) || expr.args[1] === :(==))
            return vcat((extractpins(a, default_element) for a in expr.args[2:end])...)
        elseif expr.head === :comparison && all(c -> c === :(==), expr.args[2:2:end])
            return vcat((extractpins(a, default_element) for a in expr.args[1:2:end])...)
        elseif expr.head === :ref
            return [:(($(QuoteNode(expr.args[1])), $(QuoteNode(expr.args[2]))))]
        elseif expr.head === :vect && length(expr.args) == 1
            if default_element === nothing
                error("missing element$(locinfo): $expr")
            end
            return [:(($(QuoteNode(default_element)), $(QuoteNode(expr.args[1]))))]
        else
            error("invalid pin specification$(locinfo): $expr")
        end
    end

    function extractpins(netname::Symbol, default_element=nothing)
        return [QuoteNode(netname)]
    end

    extractpins(netname::String, default_element=nothing) =
        extractpins(Symbol(netname), default_element)

    if !isa(cdef, Expr) || cdef.head !== :block
        error("@circuit must be followed by a begin/end block")
    end
    ccode = Expr(:block)
    push!(ccode.args, :(circ = Circuit()))
    locinfo = ""
    for expr in cdef.args
        if isa(expr, LineNumberNode)
            locinfo = " at $(expr.file):$(expr.line)"
            continue
        end
        if !isa(expr, Expr)
            error("invalid statement in circuit definition$locinfo: $expr")
        end
        if expr.head === :line
            locinfo = " at $(expr.args[2]):$(expr.args[1])"
        elseif expr.head === :(=)
            elem_spec(expr)
        elseif is_conn_spec(expr)
            push!(ccode.args, Expr(:call, :connect!, :circ, extractpins(expr)...))
        else
            error("invalid statement in circuit definition$locinfo: $expr")
        end
    end

    push!(ccode.args, :(circ))
    return ccode
end

@doc doc"""
    composite_element(circ; pinmap=Dict(), ports)

Create a circuit element from the (sub-)circuit `circ`. The `pinmap` defines
the mapping from pin names of the element to be created to pins (or nets) in
`circ`. Optionally, `ports` can be used to explicitly specify (as a list `Pair`s
of pins) the ports to use. By default, a port will be created between one
arbitrarily chosen pin and every other pin.

# Example

The following will create an element for a 2.5V source, created using a 5V
source and a voltage divider, stabilized with a capacitor.

```jldoctest; output = false, setup = :(using ACME), filter = r"(ACME\.)?Element\(.*"s
circ = @circuit begin
   r1 = resistor(10e3)
   r2 = resistor(10e3), [1] == r1[2]
   c = capacitor(1e-6), [1] == r2[1], [2] == r2[2]
   src = voltagesource(5), [+] == r1[1], [-] == r2[2]
end

composite_element(circ, pinmap=Dict(1 => (:r2, 1), 2 => (:r2, 2)))

# output

Element(...)
```

Note that the pins still implicitly specifiy ports, so an even number must be
given, but the same pin may be given multiple times to be part of multiple
ports.
""" function composite_element(circ::Circuit; pinmap=Dict(), ports=ports_from_pinmap(pinmap))
    if ny(circ) > 0
        throw(ArgumentError("creating composite elements from circuits with outputs is not supported"))
    end
    numports = length(ports)
    # construct system matrix with norators for connection ports included
    Mᵥ = blockdiag(mv(circ), spzeros(numports, numports))
    Mᵢ = blockdiag(mi(circ), spzeros(numports, numports))
    Mₓ = [mx(circ); spzeros(numports, nx(circ))]
    Mₓ´ = [mxd(circ); spzeros(numports, nx(circ))]
    Mq = [mq(circ); spzeros(numports, nq(circ))]
    Mu = [mu(circ); spzeros(numports, nu(circ))]
    u0 = [ACME.u0(circ); spzeros(numports)]
    incid = [incidence(circ) spzeros(Int, length(circ.nets), numports)]
    for i in eachindex(ports)
        port = ports[i]
        net = netfor!(circ, pinmap[port[1]])
        row = findfirst(==(net), circ.nets)
        b = nb(circ)+i
        incid[row, b] = 1
        net = netfor!(circ, pinmap[port[2]])
        row = findfirst(==(net), circ.nets)
        incid[row, b] = -1
    end

    tv, ti = topomat!(incid)
    S = SparseMatrixCSC{Rational{BigInt}}([Mᵥ Mᵢ Mₓ Mₓ´ Mq;
        blockdiag(tv, ti) spzeros(nb(circ)+numports, 2nx(circ)+nq(circ))])
    ũ, M = gensolve(S, SparseMatrixCSC{Rational{BigInt}}([Mu u0; spzeros(nb(circ)+numports, nu(circ)+1)]))
    # [v' i' x' xd' q']' = ũ + My for arbitrary y
    # now drop all rows concerning only internal voltages and currents
    indices = vcat(consecranges([nb(circ), numports, nb(circ), numports+2nx(circ)+nq(circ)])[[2,4]]...)
    ũ = ũ[indices,:]
    M = M[indices,:]
    S̃ = SparseMatrixCSC(gensolve(M', spzeros(size(M,2), 0))[2]') # output matrices p as RHS?
    # S̃ spans nullspace of M', so that
    #    S̃*[v' i' x' xd' q']' = S̃*ũ + S̃*My = S̃*ũ
    # i.e. S̃ acts as a new system matrix
    M̃ᵥ, M̃ᵢ, M̃ₓ, M̃ₓ´, M̃q =
        matsplit(S̃, [size(S̃,1)], [numports, numports, nx(circ), nx(circ), nq(circ)])
    M̃u = S̃*ũ[:, 1:nu(circ)]
    ũ0 = S̃*ũ[:, end]
    # note that we need to flip the sign of M̃ᵢ to view the ports from the other side
    return Element(mv = M̃ᵥ, mi = -M̃ᵢ, mx = M̃ₓ, mxd = M̃ₓ´, mq = M̃q,
        mu = M̃u, u0 = ũ0,
        nonlinear_eq=nonlinear_eq_func(circ),
        ports=ports)
end

function ports_from_pinmap(pinmap)
    ks = keys(pinmap)
    refpin = first(ks)
    ks = Iterators.drop(ks, 1)
    return [refpin => pin for pin in ks]
end
