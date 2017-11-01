# Copyright 2015, 2016, 2017 Martin Holters
# See accompanying license file.

export Circuit, add!, connect!, disconnect!

import Base: delete!

const Net = Vector{Tuple{Symbol,Symbol}} # pairs of element designator and pin name

#struct Circuit
@struct Circuit begin
    elements :: Vector{Element}
    element_names :: Dict{Symbol, Int} # map name to offset in elements
    nets :: Vector{Net}
    net_names :: Dict{Symbol, Net}
    Circuit() = new([], Dict{Symbol, Int}(), [], Dict{Symbol, Net}())
end

for n in [:nb; :nx; :nq; :nu; :nl; :ny; :nn]
    @eval ($n)(c::Circuit) = sum([$n(elem) for elem in c.elements])
end

for mat in [:mv; :mi; :mx; :mxd; :mq; :mu; :pv; :pi; :px; :pxd; :pq]
    # blkdiag() does not work, so include an empty matrix of desired type in
    # case c.elements is empty
    # as blkdiag for unknown number of arguments cannot be inferred properly,
    # add type-assertion
    @eval ($mat)(c::Circuit) =
         blkdiag(spzeros(Rational{BigInt}, 0, 0),
                 [convert(SparseMatrixCSC{Rational{BigInt}}, elem.$mat)
                  for elem in c.elements]...
                )::SparseMatrixCSC{Rational{BigInt},Int}
end

u0(c::Circuit) = vcat([elem.u0 for elem in c.elements]...)

function incidence(c::Circuit)
    i = sizehint!(Int[], 2nb(c))
    j = sizehint!(Int[], 2nb(c))
    v = sizehint!(Int[], 2nb(c))
    for (row, pins) in enumerate(c.nets), (elemname, pinname) in pins
        elem = c.elements[c.element_names[elemname]]
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

function nonlinear_eq(c::Circuit, elem_idxs=1:length(c.elements))
    # construct a block expression containing all element's expressions after
    # offsetting their indexes into q, J and res

    row_offset = 0
    col_offset = 0
    nl_expr = Expr(:block)
    for elem in c.elements[elem_idxs]
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
                append!(ret.args, map(offset_indexes, expr.args))
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

function add!(c::Circuit, elem::Element)
    for (k, v) in c.element_names
        if c.elements[v] == elem
            return k
        end
    end
    designator = gensym()
    add!(c, designator, elem)
    return designator
end

add!(c::Circuit, elems::Element...) = ([add!(c, elem) for elem in elems]...,)

function add!(c::Circuit, designator::Symbol, elem::Element)
    idx = findfirst(e -> e == elem, c.elements)
    if idx == 0
        push!(c.elements, elem)
        for pin in keys(elem.pins)
            push!(c.nets, [(designator, pin)])
        end
        idx = length(c.elements)
    end
    c.element_names[designator] = idx
end

function delete!(c::Circuit, designator::Symbol)
    idx = c.element_names[designator]
    elem = c.elements[idx]
    for net in c.nets
        filter!(elempin -> elempin[1] != designator, net)
    end
    delete!(c.element_names, designator)
    deleteat!(c.elements, idx)
    for (des, i) in c.element_names
        if i > idx
            c.element_names[des] = i - 1
        end
    end
end

function branch_offset(c::Circuit, elem::Element)
    offset = 0
    for el in c.elements
        el == elem && return offset
        offset += nb(el)
    end
    throw(ArgumentError("Element not found in circuit"))
end

function netfor!(c::Circuit, p::Pin)
    element = p[1]
    designator = add!(c, element)
    local pinname
    for (pname, pbps) in element.pins
        if pbps == p[2]
            pinname = pname
            break
        end
    end
    for net in c.nets
        (designator, pinname) ∈ net && return net
    end
end

function netfor!(c::Circuit, name::Symbol)
    haskey(c.net_names, name) || push!(c.nets, get!(c.net_names, name, []))
    c.net_names[name]
end

function connect!(c::Circuit, pins::Union{Pin,Symbol}...)
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

function disconnect!(c::Circuit, pin::Pin)
    element = pin[1]
    net = netfor!(c, pin)
    designator = add!(c, element)
    local pinname
    for (pname, pbps) in element.pins
        if pbps == pin[2]
            pinname = pname
            break
        end
    end
    filter!(p -> p != (designator, pinname), net)
    push!(c.nets, [(designator, pinname)])
end

# lines marked with !SV avoid creation of SparseVector by indexing with Ranges
# instead of Ints; a better way for cross-julia-version compatibilty would be
# nice; maybe Compat helps in the future...
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
        incidence[rows, :] .-= incidence[row:row, :] # !SV
        rows = find(incidence[1:row-1, col] .== -1)
        incidence[rows, :] .+= incidence[row:row, :] # !SV
        row += 1
    end

    ti = incidence[1:row-1, :]

    dl = ti[:, broadcast(!, t)]
    tv = spzeros(T, size(dl, 2), size(incidence, 2))
    tv[:, t] = -dl.'
    tv[:, broadcast(!, t)] = speye(T, size(dl, 2))

    tv, ti
end

@pfunction topomat(incidence::SparseMatrixCSC{T}) [T<:Integer] begin
    topomat!(copy(incidence))
 end
topomat(c::Circuit) = topomat!(incidence(c))
