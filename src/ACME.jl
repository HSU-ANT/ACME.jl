# Copyright 2015, 2016, 2017, 2018, 2019, 2020, 2021 Martin Holters
# See accompanying license file.

module ACME

export DiscreteModel, run!, steadystate, steadystate!, linearize, ModelRunner

using IterTools: subsets
using LinearAlgebra: BLAS, I, axpy!, lu, rmul!
using Markdown: @doc_str
using OrderedCollections: OrderedDict
using ProgressMeter: @showprogress
using SparseArrays: SparseMatrixCSC, blockdiag, dropzeros!, findnz,
    nonzeros, sparse, spzeros
using StaticArrays: @SMatrix, @SVector, SMatrix, SVector
using Statistics: mean, var

include("kdtree.jl")
include("solvers.jl")

let mat_dims = Dict(
        :mv => (:nl, :nb), :mi => (:nl, :nb), :mx => (:nl, :nx),
        :mxd => (:nl, :nx), :mq => (:nl, :nq), :mu => (:nl, :nu),
        :u0 => (:nl, :n0),
        :pv => (:ny, :nb), :pi => (:ny, :nb), :px => (:ny, :nx),
        :pxd => (:ny, :nx), :pq => (:ny, :nq),
    )

    @eval function prepare_element_matrices(; $((Expr(:kw, m, nothing) for m ∈ keys(mat_dims))...))
        matrices = Dict{Symbol,SparseMatrixCSC{Real,Int}}()
        $(
            (
                quote
                    if $(mat_name) !== nothing
                        matrices[$(QuoteNode(mat_name))] = hcat($(mat_name))::AbstractMatrix
                    end
                end
                for mat_name ∈ keys(mat_dims)
            )...
        )
        sizes = Dict{Symbol,Int}(:n0 => 1)
        for (mat_name, mat) ∈ matrices
            for (sym, s) ∈ zip($(mat_dims)[mat_name], size(mat))
                if get!(sizes, sym, s) ≠ s
                    throw(ArgumentError("Inconsistent sizes for $(sym)"))
                end
            end
        end
        for (m, ns) ∈ $(mat_dims)
            if !haskey(matrices, m)
                matrices[m] = spzeros(Real, get!(sizes, ns[1], 0), get!(sizes, ns[2], 0))
            end
        end
        return matrices, sizes
    end
end

struct Element
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
  nonlinear_eq
  pins :: Dict{Symbol, Vector{Tuple{Int, Int}}}

    function Element(;nonlinear_eq=nothing, ports=nothing, pins=nothing, mat_args...)
        matrices, sizes = prepare_element_matrices(; mat_args...)
        if nonlinear_eq === nothing
            nonlinear_eq = (q) -> (SVector{0,Float64}(), SMatrix{0,0,Float64}())
        end

        if ports !== nothing
            pins = Dict{Symbol,Vector{Tuple{Int, Int}}}()
            for branch in 1:length(ports)
                push!(get!(pins, Symbol(ports[branch][1]), Tuple{Int, Int}[]), (branch, 1))
                push!(get!(pins, Symbol(ports[branch][2]), Tuple{Int, Int}[]), (branch, -1))
            end
        end
        if pins === nothing
            pins = Dict(Symbol(i) => [((i+1) ÷ 2, 2(i % 2) - 1)] for i in 1:2sizes[:nb])
        end
        return new(
            matrices[:mv], matrices[:mi], matrices[:mx], matrices[:mxd],
            matrices[:mq], matrices[:mu], matrices[:u0],
            matrices[:pv], matrices[:pi], matrices[:px], matrices[:pxd], matrices[:pq],
            nonlinear_eq,
            pins,
        )
    end
end

for (n,m) in Dict(:nb => :mv, :nx => :mx, :nq => :mq, :nu => :mu)
  @eval ($n)(e::Element) = size(e.$m, 2)
end
nl(e::Element) = size(e.mv, 1)
ny(e::Element) = size(e.pv, 1)
nn(e::Element) = nb(e) + nx(e) + nq(e) - nl(e)

nonlinear_eq_func(e::Element) = e.nonlinear_eq

include("elements.jl")

include("circuit.jl")

mutable struct DiscreteModel{Solvers}
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

    nonlinear_eq_funcs::Vector

    solvers::Solvers
    x::Vector{Float64}

    function DiscreteModel(mats::Dict{Symbol}, nonlinear_eq_funcs::Vector,
            solvers::Solvers) where {Solvers}
        model = new{Solvers}()

        for mat in (:a, :b, :c, :pexps, :dqs, :eqs, :fqprevs, :fqs, :dy, :ey, :fy, :x0, :q0s, :y0)
            setfield!(model, mat, convert(fieldtype(typeof(model), mat), mats[mat]))
        end

        model.nonlinear_eq_funcs = nonlinear_eq_funcs
        model.solvers = solvers
        model.x = zeros(nx(model))
        return model
    end
end

function DiscreteModel(circ::Circuit, t::Real, ::Type{Solver}=HomotopySolver{CachingSolver{SimpleSolver}};
                       decompose_nonlinearity=true) where {Solver}
    mats = model_matrices(circ, t)

    nns = Int[nn(e) for e in elements(circ)]
    nqs = Int[nq(e) for e in elements(circ)]
    if decompose_nonlinearity
        nl_elems = nldecompose!(mats, nns, nqs)
    else
        nl_elems = Vector{Int}[findall(nn -> nn > 0, nns)]
    end

    model_nns = Int[sum(nns[nles]) for nles in nl_elems]
    model_qidxs = [vcat(consecranges(nqs)[nles]...) for nles in nl_elems]
    split_nl_model_matrices!(mats, model_qidxs, model_nns)

    reduce_pdims!(mats)

    model_nps = size.(mats[:dqs], 1)
    model_nqs = size.(mats[:pexps], 1)

    @assert nn(circ) == sum(model_nns)

    q0s = Vector{Float64}.(mats[:q0s])
    fqs = Matrix{Float64}.(mats[:fqs])
    fqprev_fulls = Matrix{Float64}.(mats[:fqprev_fulls])

    model_nonlinear_eq_funcs = [
        let q = zeros(nq), circ_nl_func = nonlinear_eq_func(circ, nles)
            @inline function(res, J, pfull, Jq, fq, z)
                #copyto!(q, pfull + fq * z)
                copyto!(q, pfull)
                BLAS.gemv!('N', 1., fq, z, 1., q)
                res´, Jq´ = circ_nl_func(q)
                res .= res´
                Jq .= Jq´
                #copyto!(J, Jq*model.fq)
                BLAS.gemm!('N', 'N', 1., Jq, fq, 0., J)
                return nothing
            end
        end for (nles, nq) in zip(nl_elems, model_nqs)]

    nonlinear_eq_funcs = [
        @inline function (res, J, scratch, z)
            nleq(res, J, scratch[1], scratch[2], fq, z)
        end for (nleq, fq) in zip(model_nonlinear_eq_funcs, fqs)]

    init_zs = [zeros(nn) for nn in model_nns]
    for idx in eachindex(nonlinear_eq_funcs)
        q = q0s[idx] + fqprev_fulls[idx] * vcat(init_zs...)
        init_zs[idx] = initial_solution(nonlinear_eq_funcs[idx], q, model_nns[idx])
    end

    while any(np -> np == 0, model_nps)
        const_idxs = findall(iszero, model_nps)
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
        deleteat!(model_nonlinear_eq_funcs, const_idxs)
        deleteat!(nonlinear_eq_funcs, const_idxs)
        deleteat!(nl_elems, const_idxs)
        mats[:fy] = mats[:fy][:,varying_zidxs]
        mats[:c] = mats[:c][:,varying_zidxs]
        reduce_pdims!(mats)
        model_nps = size.(mats[:dqs], 1)
    end

    q0s = Array{Float64}.(mats[:q0s])
    fqs = Array{Float64}.(mats[:fqs])
    fqprev_fulls = Array{Float64}.(mats[:fqprev_fulls])
    pexps = Array{Float64}.(mats[:pexps])

    nonlinear_eq_set_ps = [
        function(scratch, p)
            pfull = scratch[1]
            #copyto!(pfull, q0 + pexp * p)
            copyto!(pfull, q0)
            BLAS.gemv!('N', 1., pexp, p, 1., pfull)
            return nothing
        end
        for (pexp, q0) in zip(pexps, q0s)]
    nonlinear_eq_calc_Jps = [
        function (scratch, Jp)
            Jq = scratch[2]
            #copyto!(Jp, Jq*pexp)
            BLAS.gemm!('N', 'N', 1., Jq, pexp, 0., Jp)
            return nothing
        end
        for pexp in pexps]
    solvers = ((Solver(
                                  ParametricNonLinEq(nonlinear_eq_funcs[idx],
                                      nonlinear_eq_set_ps[idx],
                                      nonlinear_eq_calc_Jps[idx],
                                      (zeros(model_nqs[idx]), zeros(model_nns[idx], model_nqs[idx])),
                                      model_nns[idx], model_nps[idx]),
                                  zeros(model_nps[idx]), init_zs[idx])
                for idx in eachindex(nonlinear_eq_funcs))...,)
    return DiscreteModel(mats, model_nonlinear_eq_funcs, solvers)
end

function model_matrices(circ::Circuit, t::Rational{BigInt})
    lhs = convert(SparseMatrixCSC{Rational{BigInt},Int},
                  [mv(circ) mi(circ) mxd(circ)//t+mx(circ)//2 mq(circ);
                   blockdiag(topomat(circ)...) spzeros(nb(circ), nx(circ) + nq(circ))])
    rhs = convert(SparseMatrixCSC{Rational{BigInt},Int},
                  [u0(circ) mu(circ) mxd(circ)//t-mx(circ)//2;
                          spzeros(nb(circ), 1+nu(circ)+nx(circ))])
    x, f = Matrix.(gensolve(lhs, rhs))

    rowsizes = [nb(circ); nb(circ); nx(circ); nq(circ)]
    res = Dict{Symbol,Array}(zip([:fv; :fi; :c; :fq], matsplit(f, rowsizes)))

    nullspace = gensolve(sparse(res[:fq]::Matrix{Rational{BigInt}}),
                         spzeros(Rational{BigInt}, size(res[:fq],1), 0))[2]
    indeterminates = f * nullspace

    if sum(abs2, res[:c] * nullspace) > 1e-20
        @warn "State update depends on indeterminate quantity"
    end
    while size(nullspace, 2) > 0
        i, j = argmax(abs.(nullspace)).I
        nullspace = nullspace[[1:i-1; i+1:end], [1:j-1; j+1:end]]
        f = f[:, [1:i-1; i+1:end]]
        for k in [:fv; :fi; :c; :fq]
            res[k] = res[k][:, [1:i-1; i+1:end]]
        end
    end

    merge!(res, Dict(zip([:v0 :ev :dv; :i0 :ei :di; :x0 :b :a; :q0 :eq_full :dq_full],
                         matsplit(x, rowsizes, [1; nu(circ); nx(circ)]))))
    for v in (:v0, :i0, :x0, :q0)
        res[v] = dropdims(res[v], dims=2)
    end

    p = [pv(circ) pi(circ) px(circ)//2+pxd(circ)//t pq(circ)]
    if sum(abs2, p * indeterminates) > 1e-20
        @warn "Model output depends on indeterminate quantity"
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
    a = Matrix{eltype(fq)}(I, size(fq, 2), size(fq, 2))
    if numcols ≥ size(fq,2)
        return a
    end
    for colcnt in 1:numcols
        # determine element with maximum absolute value in unprocessed columns
        # to use as pivot
        i, j = argmax(abs.(fq[:,colcnt:end])).I
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

        if all(iszero, fq[:,colcnt+1:end])
            return a
        end
    end
    return nothing
end

function nldecompose!(mats, nns, nqs)
    fq = mats[:fq]
    a = Matrix{eltype(fq)}(I, size(fq, 2), size(fq, 2))
    sub_ranges = consecranges(nqs)
    extracted_subs = Vector{Int}[]
    rem_cols = 1:size(fq, 2)
    rem_nles = BitSet(filter!(e -> nqs[e] > 0, collect(eachindex(nqs))))

    while !isempty(rem_nles)
        for sz in 1:length(rem_nles), sub in subsets(collect(rem_nles), sz)
            nn_sub = sum(nns[sub])
            a_update = tryextract(fq[[sub_ranges[sub]...;],rem_cols], nn_sub)
            if a_update !== nothing
                fq[:,rem_cols] = fq[:,rem_cols] * a_update
                a[:,rem_cols] = a[:,rem_cols] * a_update
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
    let fqsplit = vcat((matsplit(mats[:fq][qidxs,:], [length(qidxs)], model_nns) for qidxs in model_qidxs)...)
        mats[:fqs] = Matrix[fqsplit[i,i] for i in 1:length(model_qidxs)]
        mats[:fqprev_fulls] = Matrix[[fqsplit[i, 1:i-1]... zeros(eltype(mats[:fq]), length(model_qidxs[i]), sum(model_nns[i:end]))]
                                     for i in 1:length(model_qidxs)]
    end
    mats[:q0s] = Vector[mats[:q0][qidxs] for qidxs in model_qidxs]
end

function reduce_pdims!(mats::Dict)
    subcount = length(mats[:dq_fulls])
    mats[:dqs] = Vector{Matrix}(undef, subcount)
    mats[:eqs] = Vector{Matrix}(undef, subcount)
    mats[:fqprevs] = Vector{Matrix}(undef, subcount)
    mats[:pexps] = Vector{Matrix}(undef, subcount)
    offset = 0
    for idx in 1:subcount
        # decompose [dq_full eq_full] into pexp*[dq eq] with [dq eq] having minimum
        # number of rows
        pexp, dqeq = rank_factorize(sparse([mats[:dq_fulls][idx] mats[:eq_fulls][idx] mats[:fqprev_fulls][idx]]))
        mats[:pexps][idx] = pexp
        colsizes = [size(mats[m][idx], 2) for m in [:dq_fulls, :eq_fulls, :fqprev_fulls]]
        mats[:dqs][idx], mats[:eqs][idx], mats[:fqprevs][idx] = matsplit(dqeq, [size(dqeq, 1)], colsizes)

        # project pexp onto the orthogonal complement of the column space of Fq
        fq = mats[:fqs][idx]
        nn = size(fq, 2)
        fq_pinv = gensolve(sparse(fq'*fq), fq')[1]
        pexp = pexp - fq*fq_pinv*pexp
        # if the new pexp has lower rank, update
        pexp, f = rank_factorize(sparse(pexp))
        if size(pexp, 2) < size(mats[:pexps][idx], 2)
            cols = offset .+ (1:nn)
            mats[:a] = mats[:a] - mats[:c][:,cols]*fq_pinv*mats[:pexps][idx]*mats[:dqs][idx]
            mats[:b] = mats[:b] - mats[:c][:,cols]*fq_pinv*mats[:pexps][idx]*mats[:eqs][idx]
            mats[:dy] = mats[:dy] - mats[:fy][:,cols]*fq_pinv*mats[:pexps][idx]*mats[:dqs][idx]
            mats[:ey] = mats[:ey] - mats[:fy][:,cols]*fq_pinv*mats[:pexps][idx]*mats[:eqs][idx]
            for idx2 in (idx+1):subcount
                mats[:dq_fulls][idx2] = mats[:dq_fulls][idx2] - mats[:fqprev_fulls][idx2][:,cols]*fq_pinv*mats[:pexps][idx]*mats[:dqs][idx]
                mats[:eq_fulls][idx2] = mats[:eq_fulls][idx2] - mats[:fqprev_fulls][idx2][:,cols]*fq_pinv*mats[:pexps][idx]*mats[:eqs][idx]
                mats[:fqprev_fulls][idx2][:,1:offset] = mats[:fqprev_fulls][idx2][:,1:offset] - mats[:fqprev_fulls][idx2][:,cols]*fq_pinv*mats[:pexps][idx]*mats[:fqprevs][idx][:,1:offset]
            end
            mats[:pexps][idx] = pexp
            mats[:dqs][idx] = f * mats[:dqs][idx]
            mats[:eqs][idx] = f * mats[:eqs][idx]
            mats[:fqprevs][idx] = f * mats[:fqprevs][idx]
            mats[:dq_fulls][idx] = pexp * mats[:dqs][idx]
            mats[:eq_fulls][idx] = pexp * mats[:eqs][idx]
            mats[:fqprev_fulls][idx] = pexp * mats[:fqprevs][idx]
        end
        offset += nn
    end
end

function initial_solution(init_nl_eq_func::Function, q0, nn)
    # determine an initial solution with a homotopy solver that may vary q0
    # between 0 and the true q0 -> q0 takes the role of p
    nq = length(q0)
    init_nleq = ParametricNonLinEq(init_nl_eq_func, nn, nq)
    init_solver = HomotopySolver{SimpleSolver}(init_nleq, zeros(nq), zeros(nn))
    init_z = solve(init_solver, q0)
    if !hasconverged(init_solver)
        error("Failed to find initial solution")
    end
    return init_z
end

nx(model::DiscreteModel) = length(model.x0)
nq(model::DiscreteModel, subidx) = length(model.q0s[subidx])
np(model::DiscreteModel, subidx) = size(model.dqs[subidx], 1)
nu(model::DiscreteModel) = size(model.b, 2)
ny(model::DiscreteModel) = length(model.y0)
nn(model::DiscreteModel, subidx) = size(model.fqs[subidx], 2)
nn(model::DiscreteModel) = reduce(+, init=0, size(fq, 2) for fq in model.fqs)

function steadystate(model::DiscreteModel, u=zeros(nu(model)))
    IA_LU = lu(I-model.a)
    steady_z = zeros(nn(model))
    zoff = 1
    for idx in 1:length(model.solvers)
        zoff_last = zoff+nn(model,idx)-1
        steady_q0 = model.q0s[idx] + model.pexps[idx]*((model.dqs[idx]/IA_LU*model.b + model.eqs[idx])*u + (model.dqs[idx]/IA_LU*model.c + model.fqprevs[idx])*steady_z) +
            model.pexps[idx]*model.dqs[idx]/IA_LU*model.x0
        fq = model.pexps[idx]*model.dqs[idx]/IA_LU*model.c[:,zoff:zoff_last] + model.fqs[idx]
        nleq = model.nonlinear_eq_funcs[idx]
        steady_nl_eq_func =
            (res, J, scratch, z) -> nleq(res, J, scratch[1], scratch[2], fq, z)
        steady_nleq = ParametricNonLinEq(steady_nl_eq_func, nn(model, idx), nq(model, idx))
        steady_solver = HomotopySolver{SimpleSolver}(steady_nleq, zeros(nq(model, idx)),
                                                     zeros(nn(model, idx)))
        set_resabstol!(steady_solver, 1e-15)
        steady_z[zoff:zoff_last] = solve(steady_solver, steady_q0)
        if !hasconverged(steady_solver)
            error("Failed to find steady state solution")
        end
        zoff += nn(model,idx)
    end
    return IA_LU\(model.b*u + model.c*steady_z + model.x0)
end

function steadystate!(model::DiscreteModel, u=zeros(nu(model)))
    x_steady = steadystate(model, u)
    copyto!(model.x, x_steady)
    return x_steady
end

function linearize(model::DiscreteModel, usteady::AbstractVector{Float64}=zeros(nu(model)))
    xsteady = steadystate(model, usteady)
    zranges = Vector{UnitRange{Int64}}(undef, length(model.solvers))
    dzdps = Vector{Matrix{Float64}}(undef, length(model.solvers))
    dqlins = Vector{Matrix{Float64}}(undef, length(model.solvers))
    eqlins = Vector{Matrix{Float64}}(undef, length(model.solvers))
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
        zsub, dzdps[idx] = linearize(model.solvers[idx], psteady)
        copyto!(zsteady, zoff, zsub, 1, length(zsub))

        zranges[idx] = zoff:zoff+length(zsub)-1
        fqdzdps = [model.fqprevs[idx][:,zranges[n]] * dzdps[n] for n in 1:idx-1]
        dqlins[idx] = reduce(+, init=model.dqs[idx], fqdzdps .* dqlins[1:idx-1])
        eqlins[idx] = reduce(+, init=model.eqs[idx], fqdzdps .* eqlins[1:idx-1])

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
    return DiscreteModel(mats, [], ())
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

struct ModelRunner{Model<:DiscreteModel,ShowProgress}
    model::Model
    ucur::Vector{Float64}
    ps::Vector{Vector{Float64}}
    ycur::Vector{Float64}
    xnew::Vector{Float64}
    z::Vector{Float64}
    function ModelRunner{Model,ShowProgress}(model::Model) where {Model<:DiscreteModel,ShowProgress}
        ucur = Vector{Float64}(undef, nu(model))
        ps = Vector{Float64}[Vector{Float64}(undef, np(model, idx)) for idx in 1:length(model.solvers)]
        ycur = Vector{Float64}(undef, ny(model))
        xnew = Vector{Float64}(undef, nx(model))
        z = Vector{Float64}(undef, nn(model))
        return new{Model,ShowProgress}(model, ucur, ps, ycur, xnew, z)
    end
end

ModelRunner(model::Model) where {Model<:DiscreteModel} = ModelRunner{Model,true}(model)
ModelRunner(model::Model, ::Val{ShowProgress}) where {Model<:DiscreteModel,ShowProgress} =
    ModelRunner{Model,ShowProgress}(model)

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
ModelRunner(model::Model, showprogress::Bool) where {Model<:DiscreteModel} =
    ModelRunner{Model,showprogress}(model)

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
    y = Matrix{Float64}(undef, ny(runner.model), size(u, 2))
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
function run!(runner::ModelRunner{<:DiscreteModel,true}, y::AbstractMatrix{Float64},
              u::AbstractMatrix{Float64})
    checkiosizes(runner, u, y)
    @showprogress "Running model: " for n = 1:size(u, 2)
        step!(runner, y, u, n)
    end
end

function run!(runner::ModelRunner{<:DiscreteModel,false}, y::AbstractMatrix{Float64},
              u::AbstractMatrix{Float64})
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
    copyto!(ucur, 1, u, (n-1)*nu(model)+1, nu(model))
    zoff = 1
    fill!(z, 0.0)
    for idx in 1:length(model.solvers)
        p = runner.ps[idx]
        # copyto!(p, model.dqs[idx] * model.x + model.eqs[idx] * u[:,n]) + model.fqprevs[idx] * z
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
                @warn "Failed to converge while solving non-linear equation."
            else
                error("Failed to converge while solving non-linear equation, got non-finite result.")
            end
        end
        copyto!(z, zoff, zsub, 1, length(zsub))
        zoff += length(zsub)
    end
    #y[:,n] = model.dy * model.x + model.ey * u[:,n] + model.fy * z + model.y0
    if !isempty(ycur)
        copyto!(ycur, model.y0)
        BLAS.gemv!('N', 1., model.dy, model.x, 1., ycur)
        BLAS.gemv!('N', 1., model.ey, ucur, 1., ycur)
        BLAS.gemv!('N', 1., model.fy, z, 1., ycur)
        #y[:,n] = ycur
        copyto!(y, (n-1)*ny(model)+1, ycur, 1, ny(model))
    end
    #model.x = model.a * model.x + model.b * u[:,n] + model.c * z + model.x0
    if !isempty(xnew)
        copyto!(xnew, model.x0)
        BLAS.gemv!('N', 1., model.a, model.x, 1., xnew)
        BLAS.gemv!('N', 1., model.b, ucur, 1.,xnew)
        BLAS.gemv!('N', 1., model.c, z, 1., xnew)
        copyto!(model.x, xnew)
    end
end

function gensolve(a, b, x, h, thresh=0.1)
    m = size(a)[1]
    if m == 0
        return x, h
    end
    t = sortperm(vec(mapslices(ait -> count(!iszero, ait), a, dims=2))) # row indexes in ascending order of nnz
    tol = 3 * max(eps(float(eltype(a))), eps(float(eltype(h)))) * size(a, 2)
    for i in 1:m
        ait = a[t[i],:]' # ait is a row of the a matrix
        s = ait * h;
        jnz, nz_vals = findnz(s')
        nz_abs_vals = abs.(nz_vals)
        max_abs_val = reduce(max, init=zero(eltype(s)), nz_abs_vals)
        if max_abs_val ≤ tol # cosidered numerical zero
            continue
        end
        jat = jnz[nz_abs_vals .≥ thresh*max_abs_val] # cols above threshold
        j = jat[argmin(vec(mapslices(hj -> count(!iszero, hj), h[:,jat], dims=1)))]
        q = h[:,j]
        x = x + convert(typeof(x), q * ((b[t[i],:]' - ait*x) * (1 / (ait*q))))
        if size(h)[2] > 1
            h = h[:,[1:j-1;j+1:end]] - convert(typeof(h), q * s[1,[1:j-1;j+1:end]]'*(1/s[1,j]))
        else
            h = similar(h, eltype(h), (size(h)[1], 0))
        end
    end
    return x, h
end

gensolve(a, b, thresh=0.1) =
    gensolve(a, b, spzeros(promote_type(eltype(a), eltype(b)), size(a, 2), size(b, 2)), SparseMatrixCSC{eltype(a)}(I, size(a,2), size(a,2)), thresh)

function rank_factorize(a::SparseMatrixCSC)
    f = a
    nullspace = gensolve(a', spzeros(eltype(a), size(a, 2), 0))[2]
    c = Matrix{eltype(a)}(I, size(a, 1), size(a, 1))
    while size(nullspace, 2) > 0
        i, j = argmax(abs.(nullspace)).I
        c -= c[:, i] * nullspace[:, j]' / nullspace[i, j]
        c = c[:, [1:i-1; i+1:end]]
        nullspace -= nullspace[:, j] * vec(nullspace[i, :])' / nullspace[i, j]
        nullspace = nullspace[[1:i-1; i+1:end], [1:j-1; j+1:end]]
        f = f[[1:i-1; i+1:end], :]
    end
    return c, f
end

consecranges(lengths) = map((l, e) -> (e-l+1):e, lengths, cumsum(lengths))

matsplit(v::AbstractVector, rowsizes) = [v[rs] for rs in consecranges(rowsizes)]
matsplit(m::AbstractMatrix, rowsizes, colsizes=[size(m,2)]) =
    [m[rs, cs] for rs in consecranges(rowsizes), cs in consecranges(colsizes)]

precompile(resistor, (Float64,))
precompile(resistor, (Int,))
precompile(potentiometer, (Float64,))
precompile(potentiometer, (Int,))
precompile(potentiometer, (Float64, Float64))
precompile(potentiometer, (Int, Float64))
precompile(capacitor, (Float64,))
precompile(inductor, (Float64,))
precompile(inductor, (Type{Val{:JA}},))
precompile(transformer, (Float64, Float64))
precompile(transformer, (Type{Val{:JA}},))
precompile(voltagesource, ())
precompile(voltagesource, (Float64,))
precompile(voltagesource, (Int,))
precompile(currentsource, ())
precompile(currentsource, (Float64,))
precompile(currentsource, (Int,))
precompile(voltageprobe, ())
precompile(currentprobe, ())
precompile(diode, ())
precompile(bjt, (Symbol,))
precompile(mosfet, (Symbol,))
precompile(opamp, ())
precompile(opamp, (Type{Val{:macak}}, Float64, Float64, Float64))

end # module
