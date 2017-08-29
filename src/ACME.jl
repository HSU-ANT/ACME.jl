# Copyright 2015, 2016, 2017 Martin Holters
# See accompanying license file.

__precompile__()

module ACME

export Circuit, add!, connect!, DiscreteModel, run!, steadystate, steadystate!,
    linearize, ModelRunner

using ProgressMeter
using Compat
using IterTools

import Base.getindex

if VERSION ≥ v"0.6.0"
    macro pfunction(sig, params, body)
        esc(Expr(:function, Expr(:where, sig, params.args...), body))
    end
else
    macro pfunction(sig, params, body)
        ts = copy(params.args)
        if VERSION < v"0.5.0"
            for i in eachindex(ts)
                if isa(ts[i], Expr) && ts[i].head === :comparison && ts[i].args[2] === :<:
                    ts[i] = Expr(:<:, ts[i].args[1], ts[i].args[3])
                end
            end
        end
        esc(Expr(:function,
                 Expr(:call, Expr(:curly, sig.args[1], ts...),
                      sig.args[2:end]...),
                 body))
    end
end

macro expandafter(mc)
    args = copy(mc.args)
    for i in eachindex(args)
        if isa(args[i], Expr) && args[i].head === :macrocall
            args[i] = macroexpand(Compat.@__MODULE__, args[i])
        end
    end
    esc(Expr(:macrocall, args...))
end

if VERSION ≥ v"0.6.0"
    @eval macro $(:struct)(head, body)
        Expr(parse("struct Foo end").head, false, esc(head), Expr(:block, [esc(a) for a in body.args]...))
    end
    macro mutable_struct(head, body)
        Expr(parse("struct Foo end").head, true, esc(head), Expr(:block, [esc(a) for a in body.args]...))
    end
else
    @eval macro $(:struct)(head, body)
        Expr(:type, false, esc(head), Expr(:block, [esc(a) for a in body.args]...))
    end
    macro mutable_struct(head, body)
        Expr(:type, true, esc(head), Expr(:block, [esc(a) for a in body.args]...))
    end
end

include("kdtree.jl")
include("solvers.jl")


#mutable struct Element
@mutable_struct Element begin
  mv :: SparseMatrixCSC{Real,Int}
  mi :: SparseMatrixCSC{Real,Int}
  mx :: SparseMatrixCSC{Real,Int}
  mxd :: SparseMatrixCSC{Real,Int}
  mq :: SparseMatrixCSC{Real,Int}
  mu :: SparseMatrixCSC{Real,Int}
  u0 :: SparseMatrixCSC{Real,Int}
  pv :: SparseMatrixCSC{Real,Int}
  pi :: SparseMatrixCSC{Real,Int}
  px :: SparseMatrixCSC{Real,Int}
  pxd :: SparseMatrixCSC{Real,Int}
  pq :: SparseMatrixCSC{Real,Int}
  nonlinear_eq :: Expr
  pins :: Dict{Symbol, Vector{Tuple{Int, Int}}}

  function Element(;args...)
    sizes = Dict{Symbol,Int}(:n0 => 1)

    function update_sizes(mat, syms)
      for (sym, s) in zip(syms, size(mat))
        if !haskey(sizes, sym)
          sizes[sym] = s
        elseif sizes[sym] ≠ s
          error("Inconsistent sizes for ", sym)
        end
      end
    end

    function make_pin_dict(syms)
      dict = Dict{Symbol,Vector{Tuple{Int, Int}}}()
      for i in 1:length(syms)
        branch = (i+1) ÷ 2
        polarity = 2(i % 2) - 1
        push!(get!(dict, Symbol(syms[i]), []), (branch, polarity))
      end
      dict
    end
    make_pin_dict(dict::Dict) = dict

    mat_dims =
        Dict( :mv => (:nl,:nb), :mi => (:nl,:nb), :mx => (:nl,:nx),
              :mxd => (:nl,:nx), :mq => (:nl,:nq), :mu => (:nl,:nu),
              :u0 => (:nl, :n0),
              :pv => (:ny,:nb), :pi => (:ny,:nb), :px => (:ny,:nx),
              :pxd => (:ny,:nx), :pq => (:ny,:nq) )

    elem = new()
    for (key, val) in args
      if haskey(mat_dims, key)
        val = convert(SparseMatrixCSC{Real}, sparse(hcat(val))) # turn val into a sparse matrix whatever it is
        update_sizes(val, mat_dims[key])
      elseif key == :pins
        val = make_pin_dict(val)
      end
      setfield!(elem, key, val)
    end
    for (m, ns) in mat_dims
      if !isdefined(elem, m)
        setfield!(elem, m, spzeros(Real, get(sizes, ns[1], 0), get(sizes, ns[2], 0)))
      end
    end
    if !isdefined(elem, :nonlinear_eq)
      elem.nonlinear_eq = Expr(:block)
    end
    if !isdefined(elem, :pins)
      elem.pins = make_pin_dict(1:2nb(elem))
    end
    elem
  end
end

for (n,m) in Dict(:nb => :mv, :nx => :mx, :nq => :mq, :nu => :mu)
  @eval ($n)(e::Element) = size(e.$m, 2)
end
nl(e::Element) = size(e.mv, 1)
ny(e::Element) = size(e.pv, 1)
nn(e::Element) = nb(e) + nx(e) + nq(e) - nl(e)

# a Pin combines an element with a branch/polarity list
const Pin = Tuple{Element, Vector{Tuple{Int,Int}}}

# allow elem[:pin] notation to get an elements pin
getindex(e::Element, p) = (e, e.pins[Symbol(p)])

include("elements.jl")

const Net = Vector{Tuple{Int,Int}} # each net is a list of branch/polarity pairs

#struct Circuit
@struct Circuit begin
    elements :: Vector{Element}
    nets :: Vector{Net}
    net_names :: Dict{Symbol, Net}
    Circuit() = new([], [], Dict{Symbol, Net}())
end

for n in [:nb; :nx; :nq; :nu; :nl; :ny; :nn]
    @eval ($n)(c::Circuit) = sum([$n(elem) for elem in c.elements])
end

for mat in [:mv; :mi; :mx; :mxd; :mq; :mu; :pv; :pi; :px; :pxd; :pq]
    # blkdiag() does not work, so include an empty matrix of desired type in
    # case c.elements is empty
    # as blkdiag for unknown numbner of arguments cannot be inferred properly,
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
    for (row, pins) in enumerate(c.nets), (branch, polarity) in pins
        push!(i, row)
        push!(j, branch)
        push!(v, polarity)
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
    elem ∈ c.elements && return
    b_offset = nb(c)
    push!(c.elements, elem)
    for branch_pols in values(elem.pins)
        push!(c.nets, [(b_offset + b, pol) for (b, pol) in branch_pols])
    end
    nothing
end

add!(c::Circuit, elems::Element...) = for elem in elems add!(c, elem) end

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
    add!(c, element)
    b_offset = branch_offset(c, element)
    for (branch, pol) in p[2], net in c.nets
        (branch + b_offset, pol) ∈ net && return net
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
        deleteat!(c.nets, findfirst(c.nets, net))
        for (name, named_net) in c.net_names
            if named_net == net
                c.net_names[name] = nets[1]
            end
        end
    end
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
        rows = filter(r -> r ≥ row, find(incidence[:, col]))
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
            cols = find(incidence[row, :])
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

#mutable struct DiscreteModel{Solvers}
@mutable_struct DiscreteModel{Solvers} begin
    a::Matrix{Float64}
    b::Matrix{Float64}
    c::Matrix{Float64}
    x0::Vector{Float64}
    pexps::Vector{Matrix{Float64}}
    dqs::Vector{Matrix{Float64}}
    eqs::Vector{Matrix{Float64}}
    fqprevs::Vector{Matrix{Float64}}
    fqs::Vector{Matrix{Float64}}
    q0s::Vector{Vector{Float64}}
    dy::Matrix{Float64}
    ey::Matrix{Float64}
    fy::Matrix{Float64}
    y0::Vector{Float64}

    nonlinear_eqs::Vector{Expr}

    solvers::Solvers
    x::Vector{Float64}

    @expandafter @compat @pfunction (::Type{DiscreteModel{Solver}})(circ::Circuit,
            t::Float64) [Solver] begin
        Base.depwarn("DiscreteModel{Solver}(circ, t) is deprecated, use DiscreteModel(circ, t, Solver) instead.",
                     :DiscreteModel)
        DiscreteModel(circ, t, Solver)
    end

    @expandafter @compat @pfunction (::Type{DiscreteModel{Solvers}})(mats::Dict{Symbol},
            nonlinear_eqs::Vector{Expr}, solvers::Solvers) [Solvers] begin
        model = new{Solvers}()

        for mat in (:a, :b, :c, :pexps, :dqs, :eqs, :fqprevs, :fqs, :dy, :ey, :fy, :x0, :q0s, :y0)
            setfield!(model, mat, convert(fieldtype(typeof(model), mat), mats[mat]))
        end

        model.nonlinear_eqs = nonlinear_eqs
        model.solvers = solvers
        model.x = zeros(nx(model))
        return model
    end
end

@pfunction DiscreteModel(circ::Circuit, t::Real, ::Type{Solver}=HomotopySolver{CachingSolver{SimpleSolver}};
                         decompose_nonlinearity=true) [Solver] begin
    mats = model_matrices(circ, t)

    nns = Int[nn(e) for e in circ.elements]
    nqs = Int[nq(e) for e in circ.elements]
    if decompose_nonlinearity
        nl_elems = nldecompose!(mats, nns, nqs)
    else
        nl_elems = Vector{Int}[filter(e -> nn(circ.elements[e]) > 0, eachindex(circ.elements))]
    end

    model_nns = Int[sum(nns[nles]) for nles in nl_elems]
    model_qidxs = [vcat(consecranges(nqs)[nles]...) for nles in nl_elems]
    split_nl_model_matrices!(mats, model_qidxs, model_nns)

    reduce_pdims!(mats)

    model_nonlinear_eqs = [quote
        #copy!(q, pfull + fq * z)
        copy!(q, pfull)
        BLAS.gemv!('N',1.,fq,z,1.,q)
        let J=Jq
            $(nonlinear_eq(circ, nles))
        end
        #copy!(J, Jq*model.fq)
        BLAS.gemm!('N', 'N', 1., Jq, fq, 0., J)
    end for nles in nl_elems]

    model_nps = map(dq -> size(dq, 1), mats[:dqs])
    model_nqs = map(pexp -> size(pexp, 1), mats[:pexps])

    @assert nn(circ) == sum(model_nns)

    q0s = map(m -> convert(Array{Float64}, m), mats[:q0s])
    fqs = map(m -> convert(Array{Float64}, m), mats[:fqs])
    fqprev_fulls = map(m -> convert(Array{Float64}, m), mats[:fqprev_fulls])

    init_zs = [zeros(nn) for nn in model_nns]
    for idx in eachindex(model_nonlinear_eqs)
        q = q0s[idx] + fqprev_fulls[idx] * vcat(init_zs...)
        init_zs[idx] = initial_solution(model_nonlinear_eqs[idx], q, fqs[idx])
    end

    while any(np -> np == 0, model_nps)
        const_idxs = find(np -> np == 0, model_nps)
        const_zidxs = vcat(consecranges(model_nns)[const_idxs]...)
        varying_zidxs = filter(idx -> !(idx in const_zidxs), 1:sum(model_nns))
        for idx in eachindex(mats[:q0s])
            mats[:q0s][idx] += mats[:fqprev_fulls][idx][:,const_zidxs] * vcat(init_zs[const_idxs]...)
            mats[:fqprev_fulls][idx] = mats[:fqprev_fulls][idx][:,varying_zidxs]
        end
        mats[:x0] += mats[:c][:,const_zidxs] * vcat(init_zs[const_idxs]...)
        mats[:y0] += mats[:fy][:,const_zidxs] * vcat(init_zs[const_idxs]...)
        deleteat!(mats[:q0s], const_idxs)
        deleteat!(mats[:dq_fulls], const_idxs)
        deleteat!(mats[:eq_fulls], const_idxs)
        deleteat!(mats[:fqs], const_idxs)
        deleteat!(mats[:fqprev_fulls], const_idxs)
        deleteat!(init_zs, const_idxs)
        deleteat!(model_nns, const_idxs)
        deleteat!(model_nqs, const_idxs)
        deleteat!(model_nonlinear_eqs, const_idxs)
        deleteat!(nl_elems, const_idxs)
        mats[:fy] = mats[:fy][:,varying_zidxs]
        mats[:c] = mats[:c][:,varying_zidxs]
        reduce_pdims!(mats)
        model_nps = map(dq -> size(dq, 1), mats[:dqs])
    end

    q0s = map(m -> convert(Array{Float64}, m), mats[:q0s])
    fqs = map(m -> convert(Array{Float64}, m), mats[:fqs])
    fqprev_fulls = map(m -> convert(Array{Float64}, m), mats[:fqprev_fulls])
    pexps = map(m -> convert(Array{Float64}, m), mats[:pexps])

    nonlinear_eq_funcs = [eval(quote
        if VERSION < v"0.5.0-dev+2396"
            # wrap up in named function because anonymous functions are slow
            # in old Julia versions
            function $(gensym())(res, J, scratch, z)
                pfull=scratch[1]
                Jq=scratch[2]
                q=$(zeros(nq))
                fq=$fq
                $(nonlinear_eq)
                return nothing
            end
        else
            (res, J, scratch, z) ->
                let pfull=scratch[1], Jq=scratch[2], q=$(zeros(nq)), fq=$fq
                    $(nonlinear_eq)
                    return nothing
                end
        end
    end) for (nonlinear_eq, fq, nq) in zip(model_nonlinear_eqs, fqs, model_nqs)]
    nonlinear_eq_set_ps = [eval(quote
        if VERSION < v"0.5.0-dev+2396"
            # wrap up in named function because anonymous functions are slow
            # in old Julia versions
            function $(gensym())(scratch, p)
                #copy!(pfull, q0 + pexp * p)
                pfull = scratch[1]
                copy!(pfull, $q0)
                BLAS.gemv!('N', 1., $pexp, p, 1., pfull)
                return nothing
            end
        else
            (scratch, p) ->
                begin
                    pfull = scratch[1]
                    #copy!(pfull, q0 + pexp * p)
                    copy!(pfull, $q0)
                    BLAS.gemv!('N', 1., $pexp, p, 1., pfull)
                    return nothing
                end
        end
    end) for (pexp, q0) in zip(pexps, q0s)]
    nonlinear_eq_calc_Jps = [eval(quote
        if VERSION < v"0.5.0-dev+2396"
            # wrap up in named function because anonymous functions are slow
            # in old Julia versions
            function $(gensym())(scratch, Jp)
                Jq = scratch[2]
                #copy!(Jp, Jq*pexp)
                BLAS.gemm!('N', 'N', 1., Jq, $pexp, 0., Jp)
                return nothing
            end
        else
            (scratch, Jp) ->
                begin
                    Jq = scratch[2]
                    #copy!(Jp, Jq*pexp)
                    BLAS.gemm!('N', 'N', 1., Jq, $pexp, 0., Jp)
                    return nothing
                end
        end
    end) for pexp in pexps]
    solvers = ([eval(:($Solver(ParametricNonLinEq($nonlinear_eq_funcs[$idx],
                                          $nonlinear_eq_set_ps[$idx],
                                          $nonlinear_eq_calc_Jps[$idx],
                                          (zeros($model_nqs[$idx]), zeros($model_nns[$idx], $model_nqs[$idx])),
                                          $model_nns[$idx], $model_nps[$idx]),
                       zeros($model_nps[$idx]), $init_zs[$idx])))
                for idx in eachindex(model_nonlinear_eqs)]...)
    return DiscreteModel{typeof(solvers)}(mats, model_nonlinear_eqs, solvers)
end

function model_matrices(circ::Circuit, t::Rational{BigInt})
    lhs = convert(SparseMatrixCSC{Rational{BigInt},Int},
                  sparse([mv(circ) mi(circ) mxd(circ)//t+mx(circ)//2 mq(circ);
                   blkdiag(topomat(circ)...) spzeros(nb(circ), nx(circ) + nq(circ))]))
    rhs = convert(SparseMatrixCSC{Rational{BigInt},Int},
                  sparse([u0(circ) mu(circ) mxd(circ)//t-mx(circ)//2;
                          spzeros(nb(circ), 1+nu(circ)+nx(circ))]))
    x, f = map(full, gensolve(lhs, rhs))

    rowsizes = [nb(circ); nb(circ); nx(circ); nq(circ)]
    res = Dict{Symbol,Array}(zip([:fv; :fi; :c; :fq], matsplit(f, rowsizes)))

    nullspace = gensolve(sparse(res[:fq]::Matrix{Rational{BigInt}}),
                         spzeros(Rational{BigInt}, size(res[:fq],1), 0))[2]
    indeterminates = f * nullspace

    if normsquared(res[:c] * nullspace) > 1e-20
        warn("State update depends on indeterminate quantity")
    end
    while size(nullspace, 2) > 0
        i, j = ind2sub(size(nullspace), indmax(map(abs, nullspace)))
        nullspace = nullspace[[1:i-1; i+1:end], [1:j-1; j+1:end]]
        f = f[:, [1:j-1; j+1:end]]
        for k in [:fv; :fi; :c; :fq]
            res[k] = res[k][:, [1:j-1; j+1:end]]
        end
    end

    merge!(res, Dict(zip([:v0 :ev :dv; :i0 :ei :di; :x0 :b :a; :q0 :eq_full :dq_full],
                         matsplit(x, rowsizes, [1; nu(circ); nx(circ)]))))
    for v in (:v0, :i0, :x0, :q0)
        res[v] = squeeze(res[v], 2)
    end

    p = [pv(circ) pi(circ) px(circ)//2+pxd(circ)//t pq(circ)]
    if normsquared(p * indeterminates) > 1e-20
        warn("Model output depends on indeterminate quantity")
    end
    res[:dy] = p * x[:,2+nu(circ):end] + px(circ)//2-pxd(circ)//t
    #          p * [dv; di; a;  dq_full] + px(circ)//2-pxd(circ)//t
    res[:ey] = p * x[:,2:1+nu(circ)] # p * [ev; ei; b;  eq_full]
    res[:fy] = p * f                 # p * [fv; fi; c;  fq]
    res[:y0] = p * vec(x[:,1])       # p * [v0; i0; x0; q0]

    return res
end

model_matrices(circ::Circuit, t) = model_matrices(circ, Rational{BigInt}(t))

function tryextract(fq, numcols)
    a = eye(eltype(fq), size(fq,2))
    if numcols ≥ size(fq,2)
        return Nullable(a)
    end
    for colcnt in 1:numcols
        # determine element with maximum absolute value in unprocessed columns
        # to use as pivot
        i, j = ind2sub(size(fq), indmax(map(abs, fq[:,colcnt:end])))
        j += colcnt-1

        # swap pivot to first (unprocessed) column
        fq[:,[colcnt, j]] = fq[:,[j, colcnt]]
        a[:,[colcnt, j]] = a[:,[j, colcnt]]

        # elimnate remaining columns in i-th row and perform equivalent
        # transformation to a
        jj = colcnt+1:size(fq,2)
        a[:,jj] -= a[:,colcnt] * (fq[[i],jj] / fq[i,colcnt])
        fq[:,jj] -= fq[:,colcnt] * (fq[[i],jj] / fq[i,colcnt])
        # ignore i-th row in following processing steps
        fq = fq[[1:i-1; i+1:end],:]

        if countnz(fq[:,colcnt+1:end]) == 0
            return Nullable(a)
        end
    end
    return Nullable{typeof(a)}()
end

function nldecompose!(mats, nns, nqs)
    fq = mats[:fq]
    a = eye(eltype(fq), size(fq,2))
    sub_ranges = consecranges(nqs)
    extracted_subs = Vector{Int}[]
    rem_cols = 1:size(fq, 2)
    rem_nles = IntSet(filter!(e -> nqs[e] > 0, collect(eachindex(nqs))))

    while !isempty(rem_nles)
        for sz in 1:length(rem_nles), sub in subsets(collect(rem_nles), sz)
            nn_sub = sum(nns[sub])
            maybe_a = tryextract(fq[[sub_ranges[sub]...;],rem_cols], nn_sub)
            if !isnull(maybe_a)
                fq[:,rem_cols] = fq[:,rem_cols] * get(maybe_a)
                a[:,rem_cols] = a[:,rem_cols] * get(maybe_a)
                rem_cols = first(rem_cols)+nn_sub:size(fq, 2)
                push!(extracted_subs, sub)
                for nle in sub
                    delete!(rem_nles, nle)
                end
                break
            end
        end
    end

    mats[:c] = mats[:c] * a
    # mats[:fq] is updated as part of the loop
    mats[:fy] = mats[:fy] * a
    return extracted_subs
end


function split_nl_model_matrices!(mats, model_qidxs, model_nns)
    mats[:dq_fulls] = Matrix[mats[:dq_full][qidxs,:] for qidxs in model_qidxs]
    mats[:eq_fulls] = Matrix[mats[:eq_full][qidxs,:] for qidxs in model_qidxs]
    let fqsplit = vcat([matsplit(mats[:fq][qidxs,:], [length(qidxs)], model_nns) for qidxs in model_qidxs]...)
        mats[:fqs] = Matrix[fqsplit[i,i] for i in 1:length(model_qidxs)]
        mats[:fqprev_fulls] = Matrix[[fqsplit[i, 1:i-1]... zeros(eltype(mats[:fq]), length(model_qidxs[i]), sum(model_nns[i:end]))]
                                     for i in 1:length(model_qidxs)]
    end
    mats[:q0s] = Vector[mats[:q0][qidxs] for qidxs in model_qidxs]
end

function reduce_pdims!(mats::Dict)
    subcount = length(mats[:dq_fulls])
    mats[:dqs] = Vector{Matrix}(subcount)
    mats[:eqs] = Vector{Matrix}(subcount)
    mats[:fqprevs] = Vector{Matrix}(subcount)
    mats[:pexps] = Vector{Matrix}(subcount)
    for idx in 1:subcount
        # decompose [dq_full eq_full] into pexp*[dq eq] with [dq eq] having minimum
        # number of rows
        pexp, dqeq = rank_factorize(sparse([mats[:dq_fulls][idx] mats[:eq_fulls][idx] mats[:fqprev_fulls][idx]]))
        mats[:pexps][idx] = pexp
        colsizes = [size(mats[m][idx], 2) for m in [:dq_fulls, :eq_fulls, :fqprev_fulls]]
        mats[:dqs][idx], mats[:eqs][idx], mats[:fqprevs][idx] = matsplit(dqeq, [size(dqeq, 1)], colsizes)

        # project pexp onto the orthogonal complement of the column space of Fq
        fq = mats[:fqs][idx]
        fq_pinv = gensolve(sparse(fq'*fq), fq')[1]
        pexp = pexp - fq*fq_pinv*pexp
        # if the new pexp has lower rank, update
        pexp, f = rank_factorize(sparse(pexp))
        if size(pexp, 2) < size(mats[:pexps][idx], 2)
            mats[:pexps][idx] = pexp
            mats[:dqs][idx] = f * mats[:dqs][idx]
            mats[:eqs][idx] = f * mats[:eqs][idx]
            mats[:fqprevs][idx] = f * mats[:fqprevs][idx]
        end
    end
end

function initial_solution(nleq, q0, fq)
    # determine an initial solution with a homotopy solver that may vary q0
    # between 0 and the true q0 -> q0 takes the role of p
    nq, nn = size(fq)
    return eval(quote
        init_nl_eq_func = (res, J, scratch, z) ->
            let pfull=scratch[1], Jp=scratch[2], q=$(zeros(nq)),
                Jq=$(zeros(nn, nq)), fq=$(fq)
                $(nleq)
                return nothing
            end
        init_nleq = ParametricNonLinEq(init_nl_eq_func, $nn, $nq)
        init_solver = HomotopySolver{SimpleSolver}(init_nleq, zeros($nq), zeros($nn))
        init_z = solve(init_solver, $q0)
        if !hasconverged(init_solver)
            error("Failed to find initial solution")
        end
        return init_z
    end)
end

nx(model::DiscreteModel) = length(model.x0)
nq(model::DiscreteModel, subidx) = length(model.q0s[subidx])
np(model::DiscreteModel, subidx) = size(model.dqs[subidx], 1)
nu(model::DiscreteModel) = size(model.b, 2)
ny(model::DiscreteModel) = length(model.y0)
nn(model::DiscreteModel, subidx) = size(model.fqs[subidx], 2)
nn(model::DiscreteModel) = sum([size(fq, 2) for fq in model.fqs])

function steadystate(model::DiscreteModel, u=zeros(nu(model)))
    IA_LU = lufact(eye(nx(model))-model.a)
    steady_z = zeros(nn(model))
    zoff = 1
    for idx in 1:length(model.solvers)
        zoff_last = zoff+nn(model,idx)-1
        steady_q0 = model.q0s[idx] + model.pexps[idx]*((model.dqs[idx]/IA_LU*model.b + model.eqs[idx])*u + (model.dqs[idx]/IA_LU*model.c + model.fqprevs[idx])*steady_z) +
            model.pexps[idx]*model.dqs[idx]/IA_LU*model.x0
        steady_z[zoff:zoff_last] = eval(quote
            steady_nl_eq_func = (res, J, scratch, z) ->
                let pfull=scratch[1], Jp=scratch[2],
                    q=$(zeros(nq(model, idx))), Jq=$(zeros(nn(model, idx), nq(model, idx))),
                    fq=$(model.pexps[idx]*model.dqs[idx]/IA_LU*model.c[:,zoff:zoff_last] + model.fqs[idx])
                    $(model.nonlinear_eqs[idx])
                    return nothing
                end
            steady_nleq = ParametricNonLinEq(steady_nl_eq_func, nn($model, $idx), nq($model, $idx))
            steady_solver = HomotopySolver{SimpleSolver}(steady_nleq, zeros(nq($model, $idx)),
                                                         zeros(nn($model, $idx)))
            set_resabstol!(steady_solver, 1e-15)
            steady_z = solve(steady_solver, $steady_q0)
            if !hasconverged(steady_solver)
                error("Failed to find steady state solution")
            end
            return steady_z
        end)
        zoff += nn(model,idx)
    end
    return IA_LU\(model.b*u + model.c*steady_z + model.x0)
end

function steadystate!(model::DiscreteModel, u=zeros(nu(model)))
    x_steady = steadystate(model, u)
    copy!(model.x, x_steady)
    return x_steady
end

function linearize(model::DiscreteModel, usteady::AbstractVector{Float64}=zeros(nu(model)))
    xsteady = steadystate(model, usteady)
    zranges = Vector{UnitRange{Int64}}(length(model.solvers))
    dzdps = Vector{Matrix{Float64}}(length(model.solvers))
    dqlins = Vector{Matrix{Float64}}(length(model.solvers))
    eqlins = Vector{Matrix{Float64}}(length(model.solvers))
    zsteady = zeros(nn(model))
    zoff = 1
    x0 = copy(model.x0)
    a = copy(model.a)
    b = copy(model.b)
    c = copy(model.c)
    y0 = copy(model.y0)
    dy = copy(model.dy)
    ey = copy(model.ey)
    fy = copy(model.fy)

    for idx in 1:length(model.solvers)
        psteady = model.dqs[idx] * xsteady + model.eqs[idx] * usteady +
                  model.fqprevs[idx] * zsteady
        zsub, dzdps[idx] =
            Compat.invokelatest(linearize, model.solvers[idx], psteady)
        copy!(zsteady, zoff, zsub, 1, length(zsub))

        zranges[idx] = zoff:zoff+length(zsub)-1
        fqdzdps = [model.fqprevs[idx][:,zranges[n]] * dzdps[n] for n in 1:idx-1]
        dqlins[idx] = reduce(+, model.dqs[idx], fqdzdps .* dqlins[1:idx-1])
        eqlins[idx] = reduce(+, model.eqs[idx], fqdzdps .* eqlins[1:idx-1])

        x0 += model.c[:,zranges[idx]] * (zsub - dzdps[idx]*psteady)
        a += model.c[:,zranges[idx]] * dzdps[idx] * dqlins[idx]
        b += model.c[:,zranges[idx]] * dzdps[idx] * eqlins[idx]

        y0 += model.fy[:,zranges[idx]] * (zsub - dzdps[idx]*psteady)
        dy += model.fy[:,zranges[idx]] * dzdps[idx] * dqlins[idx]
        ey += model.fy[:,zranges[idx]] * dzdps[idx] * eqlins[idx]

        zoff += length(zsub)
    end

    mats = Dict(:a => a, :b => b, :c => zeros(nx(model), 0),
        :pexps => Matrix{Float64}[], :dqs => Matrix{Float64}[],
        :eqs => Matrix{Float64}[], :fqprevs => Matrix{Float64}[],
        :fqs => Matrix{Float64}[], :q0s => Vector{Float64}[],
        :dy => dy, :ey => ey, :fy => zeros(ny(model), 0), :x0 => x0, :y0 => y0)
    return DiscreteModel{Tuple{}}(mats, Expr[], ())
end

"""
    run!(model::DiscreteModel, u::AbstractMatrix{Float64}; showprogress=true)

Run the given `model` by feeding it the input `u` which must be a matrix with
one row for each of the circuit's inputs and one column for each time step to
simulate. Likewise, the returned output will be a matrix with one row for each
of the circuit's outputs and one column for each simulated time step. The order
of the rows will correspond to the order in which the respective input and
output elements were added to the `Circuit`. To simulate a circuit without
inputs, a matrix with zero rows may be passed. The internal state of the model
(e.g. capacitor charges) is preserved accross calls to `run!`.

By default `run!` will show a progress bar to report its progress. This can be
disabled by passing `showprogress=false`.
"""
run!(model::DiscreteModel, u::AbstractMatrix{Float64}; showprogress=true) =
    return run!(ModelRunner(model, showprogress), u)

#struct ModelRunner{Model<:DiscreteModel,ShowProgress}
@struct ModelRunner{Model<:DiscreteModel,ShowProgress} begin
    model::Model
    ucur::Vector{Float64}
    ps::Vector{Vector{Float64}}
    ycur::Vector{Float64}
    xnew::Vector{Float64}
    z::Vector{Float64}
    @expandafter @compat @pfunction (::Type{ModelRunner{Model,ShowProgress}})(
            model::Model) [Model<:DiscreteModel,ShowProgress] begin
        ucur = Array{Float64,1}(nu(model))
        ps = Vector{Float64}[Vector{Float64}(np(model, idx)) for idx in 1:length(model.solvers)]
        ycur = Array{Float64,1}(ny(model))
        xnew = Array{Float64,1}(nx(model))
        z = Array{Float64,1}(nn(model))
        return new{Model,ShowProgress}(model, ucur, ps, ycur, xnew, z)
    end
end

@pfunction ModelRunner(model::Model) [Model<:DiscreteModel] begin
     ModelRunner{Model,true}(model)
 end
@pfunction ModelRunner(model::Model, ::Val{ShowProgress}) [Model<:DiscreteModel,ShowProgress] begin
    ModelRunner{Model,ShowProgress}(model)
end

"""
    ModelRunner(model::DiscreteModel, showprogress::Bool = true)

Construct a `ModelRunner` instance for running `model`. The `ModelRunner`
instance pre-allocates some memory required for model execution. Hence, when
running the same model for multiple small input data blocks, some overhead can
be saved by explicitly using a `ModelRunner`.

By default `run!` for the constructed `ModelRunner` will show a progress bar to
report its progress. This can be disabled by passing `false` as second
parameter.
"""
@pfunction ModelRunner(model::Model, showprogress::Bool) [Model<:DiscreteModel] begin
    ModelRunner{Model,showprogress}(model)
end

"""
    run!(runner::ModelRunner, u::AbstractMatrix{Float64})

Run the given `runner` by feeding it the input `u` which must be a matrix with
one row for each of the circuit's inputs and one column for each time step to
simulate. Likewise, the returned output will be a matrix with one row for each
of the circuit's outputs and one column for each simulated time step. The order
of the rows will correspond to the order in which the respective input and
output elements were added to the `Circuit`. To simulate a circuit without
inputs, a matrix with zero rows may be passed. The internal state of the
underlying `DiscreteModel` (e.g. capacitor charges) is preserved accross calls
to `run!`.
"""
function run!(runner::ModelRunner, u::AbstractMatrix{Float64})
    y = Array{Float64,2}(ny(runner.model), size(u, 2))
    run!(runner, y, u)
    return y
end

function checkiosizes(runner::ModelRunner, u::AbstractMatrix{Float64}, y::AbstractMatrix{Float64})
    if size(u, 1) ≠ nu(runner.model)
        throw(DimensionMismatch("input matrix has $(size(u,1)) rows, but model has $(nu(runner.model)) inputs"))
    end
    if size(y, 1) ≠ ny(runner.model)
        throw(DimensionMismatch("output matrix has $(size(y,1)) rows, but model has $(ny(runner.model)) outputs"))
    end
    if size(u, 2) ≠ size(y, 2)
        throw(DimensionMismatch("input matrix has $(size(u,2)) columns, output matrix has $(size(y,2)) columns"))
    end
end

"""
    run!(runner::ModelRunner, y::AbstractMatrix{Float64}, u::AbstractMatrix{Float64})

Run the given `runner` by feeding it the input `u` and storing the output
in `y`. The input `u` must be a matrix with one row for each of the circuit's
inputs and one column for each time step to simulate. Likewise, the output `y`
must be a matrix with one row for each of the circuit's outputs and one column
for each simulated time step. The order of the rows will correspond to the order
in which the respective input and output elements were added to the `Circuit`.
To simulate a circuit without inputs, a matrix with zero rows may be passed.
The internal state of the  underlying `DiscreteModel` (e.g. capacitor charges)
is preserved accross calls to `run!`.
"""
@pfunction run!(runner::ModelRunner{Model,true}, y::AbstractMatrix{Float64},
                u::AbstractMatrix{Float64}) [Model<:DiscreteModel] begin
    checkiosizes(runner, u, y)
    @showprogress "Running model: " for n = 1:size(u, 2)
        step!(runner, y, u, n)
    end
end

@pfunction run!(runner::ModelRunner{Model,false}, y::AbstractMatrix{Float64},
                u::AbstractMatrix{Float64}) [Model<:DiscreteModel] begin
    checkiosizes(runner, u, y)
    for n = 1:size(u, 2)
        step!(runner, y, u, n)
    end
end

function step!(runner::ModelRunner, y::AbstractMatrix{Float64}, u::AbstractMatrix{Float64}, n)
    model = runner.model
    ucur = runner.ucur
    ycur = runner.ycur
    xnew = runner.xnew
    z = runner.z
    copy!(ucur, 1, u, (n-1)*nu(model)+1, nu(model))
    zoff = 1
    fill!(z, 0.0)
    for idx in 1:length(model.solvers)
        p = runner.ps[idx]
        # copy!(p, model.dqs[idx] * model.x + model.eqs[idx] * u[:,n]) + model.fqprevs[idx] * z
        if size(model.dqs[idx], 2) == 0
            fill!(p, 0.0)
        else
            BLAS.gemv!('N', 1., model.dqs[idx], model.x, 0., p)
        end
        BLAS.gemv!('N', 1., model.eqs[idx], ucur, 1., p)
        if idx > 1
            BLAS.gemv!('N', 1., model.fqprevs[idx], z, 1., p)
        end
        zsub = solve(model.solvers[idx], p)
        if !hasconverged(model.solvers[idx])
            if all(isfinite, zsub)
                warn("Failed to converge while solving non-linear equation.")
            else
                error("Failed to converge while solving non-linear equation, got non-finite result.")
            end
        end
        copy!(z, zoff, zsub, 1, length(zsub))
        zoff += length(zsub)
    end
    #y[:,n] = model.dy * model.x + model.ey * u[:,n] + model.fy * z + model.y0
    copy!(ycur, model.y0)
    BLAS.gemv!('N', 1., model.dy, model.x, 1., ycur)
    BLAS.gemv!('N', 1., model.ey, ucur, 1., ycur)
    BLAS.gemv!('N', 1., model.fy, z, 1., ycur)
    #y[:,n] = ycur
    copy!(y, (n-1)*ny(model)+1, ycur, 1, ny(model))
    #model.x = model.a * model.x + model.b * u[:,n] + model.c * z + model.x0
    copy!(xnew, model.x0)
    BLAS.gemv!('N', 1., model.a, model.x, 1., xnew)
    BLAS.gemv!('N', 1., model.b, ucur, 1.,xnew)
    BLAS.gemv!('N', 1., model.c, z, 1., xnew)
    copy!(model.x, xnew)
end

# lines marked with !SV avoid creation of SparseVector by indexing with Ranges
# instead of Ints; a better way for cross-julia-version compatibilty would be
# nice; maybe Compat helps in the future...
function gensolve(a::SparseMatrixCSC, b, x, h, thresh=0.1)
    m = size(a)[1]
    t = sortperm(vec(sum(spones(a),2))) # row indexes in ascending order of nnz
    tol = 3 * max(eps(float(eltype(a))), eps(float(eltype(h)))) * size(a, 2)
    for i in 1:m
        ait = a[t[i:i],:] # ait is a row of the a matrix # !SV
        s = ait * h;
        inz, jnz, nz_vals = findnz(s)
        nz_abs_vals = map(abs, nz_vals)
        max_abs_val = reduce(max, zero(eltype(s)), nz_abs_vals)
        if max_abs_val ≤ tol # cosidered numerical zero
            continue
        end
        jat = jnz[nz_abs_vals .≥ thresh*max_abs_val] # cols above threshold
        j = jat[indmin(sum(spones(h[:,jat])))]
        q = h[:,j:j] # !SV
        # ait*q only has a single element!
        x = x + convert(typeof(x), q * ((b[t[i:i],:] - ait*x) * (1 / (ait*q)[1]))); # !SV
        if size(h)[2] > 1
            h = h[:,[1:j-1;j+1:end]] - convert(typeof(h), q * s[1:1,[1:j-1;j+1:end]]*(1/s[1,j])) # !SV
        else
            h = similar(h, eltype(h), (size(h)[1], 0))
        end
    end
    return x, h
end

gensolve(a, b, thresh=0.1) =
    gensolve(a, b, spzeros(promote_type(eltype(a), eltype(b)), size(a)[2], size(b)[2]), speye(eltype(a), size(a)[2]), thresh)

function rank_factorize(a::SparseMatrixCSC)
    f = a
    nullspace = gensolve(a', spzeros(size(a, 2), 0))[2]
    c = eye(eltype(a), size(a, 1))
    while size(nullspace, 2) > 0
        i, j = ind2sub(size(nullspace), indmax(map(abs, nullspace)))
        c -= c[:, i] * nullspace[:, j]' / nullspace[i, j]
        c = c[:, [1:i-1; i+1:end]]
        nullspace -= nullspace[:, j] * vec(nullspace[i, :])' / nullspace[i, j]
        nullspace = nullspace[[1:i-1; i+1:end], [1:j-1; j+1:end]]
        f = f[[1:i-1; i+1:end], :]
    end
    return c, f
end

consecranges(lengths) = isempty(lengths) ? [] : map(range, cumsum([1; lengths[1:end-1]]), lengths)

matsplit(v::AbstractVector, rowsizes) = [v[rs] for rs in consecranges(rowsizes)]
matsplit(m::AbstractMatrix, rowsizes, colsizes=[size(m,2)]) =
    [m[rs, cs] for rs in consecranges(rowsizes), cs in consecranges(colsizes)]

    if VERSION < v"0.5.0"
        # this is deprecated starting in Julia 0.6
        normsquared(x) = sumabs2(x)
    else
        # prior to Julia 0.5, this fails for empty x and is slow
        normsquared(x) = sum(abs2, x)
    end

end # module
