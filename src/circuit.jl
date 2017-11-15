# Copyright 2015, 2016, 2017 Martin Holters
# See accompanying license file.

export Circuit, add!, connect!, disconnect!, @circuit

import Base: delete!

const Net = Vector{Tuple{Symbol,Symbol}} # pairs of element designator and pin name

#struct Circuit
@struct Circuit begin
    elements::OrderedDict{Symbol, Element}
    nets :: Vector{Net}
    net_names :: Dict{Symbol, Net}
    Circuit() = new(OrderedDict{Symbol, Element}(), [], Dict{Symbol, Net}())
end

elements(c::Circuit) = values(c.elements)

for n in [:nb; :nx; :nq; :nu; :nl; :ny; :nn]
    @eval ($n)(c::Circuit) = sum([$n(elem) for elem in elements(c)])
end

for mat in [:mv; :mi; :mx; :mxd; :mq; :mu; :pv; :pi; :px; :pxd; :pq]
    # blkdiag() does not work, so include an empty matrix of desired type in
    # case elements(c) is empty
    # as blkdiag for unknown number of arguments cannot be inferred properly,
    # add type-assertion
    @eval ($mat)(c::Circuit) =
         blkdiag(spzeros(Rational{BigInt}, 0, 0),
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
    # ensure zeros due to short-circuited branches are removed, hence the
    # additional sparse(findnz(...))
    sparse(findnz(sparse(i,j,v))..., length(c.nets), nb(c))
end

function nonlinear_eq(c::Circuit, elem_idxs=1:length(elements(c)))
    # construct a block expression containing all element's expressions after
    # offsetting their indexes into q, J and res

    row_offset = 0
    col_offset = 0
    nl_expr = Expr(:block)
    for elem in collect(elements(c))[elem_idxs]
        index_offsets = Dict( :q => (col_offset,),
                              :J => (row_offset, col_offset),
                              :res => (row_offset,) )

        function offset_indexes(expr::Expr)
            ret = Expr(expr.head)
            ret.typ = expr.typ
            if expr.head == :ref && haskey(index_offsets, expr.args[1])
                push!(ret.args, expr.args[1])
                offsets = index_offsets[expr.args[1]]
                length(expr.args) == length(offsets) + 1 ||
                    throw(ArgumentError(string(expr.args[1], " must be indexed with exactly ", length(offsets), " index(es)")))
                for i in 1:length(offsets)
                    push!(ret.args,
                          :($(offsets[i]) + $(offset_indexes(expr.args[i+1]))))
                end
            else
                append!(ret.args, offset_indexes.(expr.args))
            end
            ret
        end

        function offset_indexes(s::Symbol)
            haskey(index_offsets, s) && throw(ArgumentError(string(s, " used without indexing expression")))
            s
        end

        offset_indexes(x::Any) = x

        # wrap into a let to keep variables local
        push!(nl_expr.args, :( let; $(offset_indexes(elem.nonlinear_eq)) end))

        row_offset += nn(elem)
        col_offset += nq(elem)
    end
    nl_expr
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

function netfor!(c::Circuit, p::Pin)
    Base.depwarn("pin specification $p is deprecated, use (refdes, pinname) instead", :netfor!)
    element = p[1]
    designator = add!(c, element)
    local pinname
    for (pname, pbps) in element.pins
        if pbps == p[2]
            pinname = pname
            break
        end
    end
    return netfor!(c, (designator, pinname))
end

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

```julia
circ = Circuit()
add!(circ, :r, resistor(1e3))
add!(circ, :src, voltages(5))
connect!(circ, (:src, -), (:r, 2), :gnd) # connect to gnd net
```
"""
function connect!(c::Circuit, pins::Union{Pin,Symbol,Tuple{Symbol,Any}}...)
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

function disconnect!(c::Circuit, pin::Pin)
    Base.depwarn("disconnect!(::Circuit, $p) is deprecated, use (refdes, pinname) to specify the the pin", :disconnect!)
    element = pin[1]
    designator = add!(c, element)
    local pinname
    for (pname, pbps) in element.pins
        if pbps == pin[2]
            pinname = pname
            break
        end
    end
    disconnect!(c, (designator, pinname))
end

@pfunction topomat!(incidence::SparseMatrixCSC{T}) [T<:Integer] begin
    @assert all(x -> abs(x) == 1, nonzeros(incidence))
    @assert all(sum(incidence, 1) .== 0)

    t = falses(size(incidence)[2]);

    row = 1;
    for col = 1:size(incidence)[2]
        rows = filter(r -> r ≥ row, find(!iszero, incidence[:, col]))
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
            cols = find(!iszero, incidence[row, :])
            incidence[row,cols] = -incidence[row,cols]
        end
        rows = find(incidence[1:row-1, col] .== 1)
        incidence[rows, :] .-= incidence[row, :]'
        rows = find(incidence[1:row-1, col] .== -1)
        incidence[rows, :] .+= incidence[row, :]'
        row += 1
    end

    ti = incidence[1:row-1, :]

    dl = ti[:, (!).(t)]
    tv = spzeros(T, size(dl, 2), size(incidence, 2))
    tv[:, t] = -dl.'
    tv[:, (!).(t)] = speye(T, size(dl, 2))

    tv, ti
end

@pfunction topomat(incidence::SparseMatrixCSC{T}) [T<:Integer] begin
    topomat!(copy(incidence))
 end
topomat(c::Circuit) = topomat!(incidence(c))

"""
    @circuit begin #= ... =# end

Provides a simple domain-specific language to decribe circuits. The
`begin`/`end` block can hold element definitions of the form
`refdes = elementfunc(params)` and connection specifications of the form
`refdes1[pin1] ⟷ refdes2[pin2]`.

# Example

To create a circuit with a voltage source connected to a resistor:
```julia
@circuit begin
    src = voltagesource(5)
    r = resistor(1000)
    src[+] ⟷ r[1]
    src[-] ⟷ r[2]
end
```

Alternatively, connection specifications can be given after an element
specification, separated by commas. In that case, the `refdes` may be omitted,
defaulting to the current element.

# Example

```julia
@circuit begin
    src = voltagesource(5)
    r = resistor(1000), src[+] ⟷ [1], src[-] ⟷ [2]
end
```

Finally, a connection endpoint may simply be of the form `netname`, to connect
to a named net. (Such named nets are created as needed.)

# Example

```julia
@circuit begin
    src = voltagesource(5), [-] ⟷ gnd
    r = resistor(1000), [1] ⟷ src[+], [2] ⟷ gnd
end
```
!!! note
    Instead of `⟷` (`\\longleftrightarrow`), one can also use `==`.
"""
macro circuit(cdef)
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
