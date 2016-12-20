# Copyright 2015, 2016 Martin Holters
# See accompanying license file.

__precompile__()

module ACME

export Circuit, add!, connect!, DiscreteModel, run!, steadystate, steadystate!

using ProgressMeter
using Compat

import Base.getindex

include("kdtree.jl")
include("solvers.jl")


type Element
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

    const mat_dims =
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
typealias Pin Tuple{Element, Vector{Tuple{Int,Int}}}

# allow elem[:pin] notation to get an elements pin
getindex(e::Element, p) = (e, e.pins[Symbol(p)])

include("elements.jl")

typealias Net Vector{Tuple{Int,Int}} # each net is a list of branch/polarity pairs

type Circuit
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
    @eval ($mat)(c::Circuit) = blkdiag(spzeros(Real, 0, 0),
                                       [elem.$mat for elem in c.elements]...)
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

function nonlinear_eq(c::Circuit)
    # construct a block expression containing all element's expressions after
    # offsetting their indexes into q, J and res

    row_offset = 0
    col_offset = 0
    nl_expr = Expr(:block)
    for elem in c.elements
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
    local net
    for (branch, pol) in p[2], net in c.nets
        (branch + b_offset, pol) ∈ net && break
    end
    net
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
function topomat!{T<:Integer}(incidence::SparseMatrixCSC{T})
    @assert all(x -> abs(x) == 1, nonzeros(incidence))
    @assert all(sum(incidence, 1) .== 0)

    t = falses(size(incidence)[2]);

    row = 1;
    for col = 1:size(incidence)[2]
        rows = filter(r -> r ≥ row, find(incidence[:,col:col])) # !SV
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
            cols = find(incidence[row:row,:]) # !SV
            incidence[row,cols] = -incidence[row,cols]
        end
        rows = find(incidence[1:row-1,col:col] .== 1) # !SV
        incidence[rows,:] = broadcast(-, incidence[rows, :], incidence[row:row,:]) # !SV
        rows = find(incidence[1:row-1,col:col] .== -1) # !SV
        incidence[rows,:] = broadcast(+, incidence[rows, :], incidence[row:row,:]) # !SV
        row += 1
    end

    if row > 1
        ti = incidence[1:row-1, :]
    else
        ti = spzeros(T, 0, size(incidence)[2])
    end

    if all(t)
        dl = spzeros(T, row-1, 0)
        tv = spzeros(T, 0, size(incidence)[2])
    else
        dl = ti[:, ~t]
        tv = spzeros(T, size(dl)[2], size(incidence)[2])
        if ~all(~t)
            tv[:,find(t)] = -dl' # with julia 0.3.2, sparse([1 -1]).' -> [1, 0], hence use of ' (conjugate transpose)
        end
        tv[:,find(~t)] = speye(T,size(dl)[2])
    end

    tv, ti
end

topomat{T<:Integer}(incidence::SparseMatrixCSC{T}) = topomat!(copy(incidence))
topomat(c::Circuit) = topomat!(incidence(c))

type DiscreteModel{Solver}
    a::Matrix{Float64}
    b::Matrix{Float64}
    c::Matrix{Float64}
    x0::Vector{Float64}
    pexp::Matrix{Float64}
    dq::Matrix{Float64}
    eq::Matrix{Float64}
    fq::Matrix{Float64}
    q0::Vector{Float64}
    dy::Matrix{Float64}
    ey::Matrix{Float64}
    fy::Matrix{Float64}
    y0::Vector{Float64}

    nonlinear_eq :: Expr

    solver::Solver
    x::Vector{Float64}

    function DiscreteModel(circ::Circuit, t::Float64)
        Base.depwarn("DiscreteModel{Solver}(circ, t) is deprecated, use DiscreteModel(circ, t, Solver) instead.",
                     :DiscreteModel)
        DiscreteModel(circ, t, Solver)
    end

    function DiscreteModel(mats::Dict{Symbol,Array{Float64}}, nonlinear_eq::Expr, solver::Solver)
        model = new()

        for mat in (:a, :b, :c, :pexp, :dq, :eq, :fq, :dy, :ey, :fy, :x0, :q0, :y0)
            setfield!(model, mat, mats[mat])
        end

        model.nonlinear_eq = nonlinear_eq
        model.solver = solver
        model.x = zeros(nx(model))
        return model
    end
end

function DiscreteModel{Solver}(circ::Circuit, t::Float64, ::Type{Solver}=HomotopySolver{CachingSolver{SimpleSolver}})
    mats = model_matrices(circ, t)
    reduce_pdims!(mats)

    model_nonlinear_eq = quote
        #copy!(q, pfull + fq * z)
        copy!(q, pfull)
        BLAS.gemv!('N',1.,fq,z,1.,q)
        let J=Jq
            $(nonlinear_eq(circ))
        end
        #copy!(J, Jq*model.fq)
        BLAS.gemm!('N', 'N', 1., Jq, fq, 0., J)
    end

    model_nq = length(mats[:q0])
    model_nn = size(mats[:fq],2)
    model_np = size(mats[:dq],1)

    @assert nn(circ) == model_nn

    init_z = initial_solution(model_nonlinear_eq, mats[:q0], mats[:fq])
    nonlinear_eq_func = eval(quote
        if VERSION < v"0.5.0-dev+2396"
            # wrap up in named function because anonymous functions are slow
            # in old Julia versions
            function $(gensym())(res, J, scratch, z)
                pfull=scratch[1]
                Jq=scratch[2]
                q=$(zeros(model_nq))
                fq=$(mats[:fq])
                $(model_nonlinear_eq)
                return nothing
            end
        else
            (res, J, scratch, z) ->
                let pfull=scratch[1], Jq=scratch[2], q=$(zeros(model_nq)),
                    fq=$(mats[:fq])
                    $(model_nonlinear_eq)
                    return nothing
                end
        end
    end)
    nonlinear_eq_set_p = eval(quote
        if VERSION < v"0.5.0-dev+2396"
            # wrap up in named function because anonymous functions are slow
            # in old Julia versions
            function $(gensym())(scratch, p)
                #copy!(pfull, q0 + pexp * p)
                pfull = scratch[1]
                copy!(pfull, $(mats[:q0]))
                BLAS.gemv!('N', 1., $(mats[:pexp]), p, 1., pfull)
                return nothing
            end
        else
            (scratch, p) ->
                begin
                    pfull = scratch[1]
                    #copy!(pfull, q0 + pexp * p)
                    copy!(pfull, $(mats[:q0]))
                    BLAS.gemv!('N', 1., $(mats[:pexp]), p, 1., pfull)
                    return nothing
                end
        end
    end)
    nonlinear_eq_calc_Jp = eval(quote
        if VERSION < v"0.5.0-dev+2396"
            # wrap up in named function because anonymous functions are slow
            # in old Julia versions
            function $(gensym())(scratch, Jp)
                Jq = scratch[2]
                #copy!(Jp, Jq*pexp)
                BLAS.gemm!('N', 'N', 1., Jq, $(mats[:pexp]), 0., Jp)
                return nothing
            end
        else
            (scratch, Jp) ->
                begin
                    Jq = scratch[2]
                    #copy!(Jp, Jq*pexp)
                    BLAS.gemm!('N', 'N', 1., Jq, $(mats[:pexp]), 0., Jp)
                    return nothing
                end
        end
    end)
    solver = Solver(ParametricNonLinEq(nonlinear_eq_func, nonlinear_eq_set_p,
                                       nonlinear_eq_calc_Jp,
                                       (zeros(model_nq), zeros(model_nn, model_nq)),
                                       model_nn, model_np),
                    zeros(model_np), init_z)
    return DiscreteModel{typeof(solver)}(mats, model_nonlinear_eq, solver)
end

function model_matrices(circ::Circuit, t)
    x, f = map(full,
               gensolve([mv(circ) mi(circ) 1/t*mxd(circ)+0.5*mx(circ) mq(circ);
                         blkdiag(topomat(circ)...) spzeros(nb(circ), nx(circ) + nq(circ))],
                        [u0(circ) mu(circ) 1/t*mxd(circ)-0.5*mx(circ);
                         spzeros(nb(circ), 1+nu(circ)+nx(circ))]))

    rowsizes = [nb(circ); nb(circ); nx(circ); nq(circ)]
    res = Dict{Symbol,Array{Float64}}(zip([:fv; :fi; :c; :fq], matsplit(f, rowsizes)))

    nullspace = gensolve(sparse(res[:fq]), spzeros(size(res[:fq],1), 0))[2]
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

    # This would choose a particular solution such that the rows corresponding
    # to q are column-wise orthogonal to the column space of fq (and hence have
    # a column space of minimal dimension). However, this destroys all sparsity
    # in x and leads to numerical difficulties in actually finding a rank
    # factorization of [dq eq], while no case has been found so far where it
    # actually reduces the rank. Hence, it is disabled for now, but a warning is
    # produced if it could be helpful (see reduce_pdims! below).
    #x = x - f*pinv(res[:fq])*x[end-nq(circ)+1:end,:]

    merge!(res, Dict(zip([:v0 :ev :dv; :i0 :ei :di; :x0 :b :a; :q0 :eq_full :dq_full],
                         matsplit(x, rowsizes, [1; nu(circ); nx(circ)]))))
    for v in (:v0, :i0, :x0, :q0)
        res[v] = squeeze(res[v], 2)
    end

    p = [pv(circ) pi(circ) 0.5*px(circ)+1/t*pxd(circ) pq(circ)]
    if normsquared(p * indeterminates) > 1e-20
        warn("Model output depends on indeterminate quantity")
    end
    res[:dy] = p * x[:,2+nu(circ):end] + 0.5*px(circ)-1/t*pxd(circ)
    #          p * [dv; di; a;  dq_full] + 0.5*px(circ)-1/t*pxd(circ)
    res[:ey] = p * x[:,2:1+nu(circ)] # p * [ev; ei; b;  eq_full]
    res[:fy] = p * f                 # p * [fv; fi; c;  fq]
    res[:y0] = p * vec(x[:,1])       # p * [v0; i0; x0; q0]

    return res
end

function reduce_pdims!(mats::Dict)
    dqeq_full = [mats[:dq_full] mats[:eq_full]]
    # decompose [dq_full eq_full] into pexp*[dq eq] with [dq eq] having minimum
    # number of rows
    nullspace = gensolve(sparse(dqeq_full'), spzeros(size(dqeq_full, 2), 0))[2]
    mats[:dq] = copy(mats[:dq_full])
    mats[:eq] = copy(mats[:eq_full])
    pexp = eye(size(dqeq_full, 1))
    while size(nullspace, 2) > 0
        i, j = ind2sub(size(nullspace), indmax(map(abs, nullspace)))
        pexp -= pexp[:, i] * nullspace[:, j]' / nullspace[i, j]
        pexp = pexp[:, [1:i-1; i+1:end]]

        nullspace -= nullspace[:, j] * nullspace[i:i, :] / nullspace[i, j]
        nullspace = nullspace[[1:i-1; i+1:end], [1:j-1; j+1:end]]

        mats[:dq] = mats[:dq][[1:i-1; i+1:end], :]
        mats[:eq] = mats[:eq][[1:i-1; i+1:end], :]
    end
    mats[:pexp] = pexp

    if rank(dqeq_full - mats[:fq]*pinv(mats[:fq])*dqeq_full) < size(pexp, 2)
        warn("Dimension of p could be further reduced by projecting onto the orthogonal complement of the column space of Fq. However, this has not been implemented due to numerical difficulties.")
    end
end

function initial_solution(nleq, q0, fq)
    # determine an initial solution with a homotopy solver that may vary q0
    # between 0 and the true q0 -> q0 takes the role of p
    nq, nn = size(fq)
    init_nl_eq_func = eval(quote
        (res, J, scratch, z) ->
            let pfull=scratch[1], Jp=scratch[2], q=$(zeros(nq)),
                Jq=$(zeros(nn, nq)), fq=$(fq)
                $(nleq)
                return nothing
            end
    end)
    init_nleq = ParametricNonLinEq(init_nl_eq_func, nn, nq)
    init_solver = HomotopySolver{SimpleSolver}(init_nleq, zeros(nq), zeros(nn))
    init_z = solve(init_solver, q0)
    if !hasconverged(init_solver)
        error("Failed to find initial solution")
    end
    return init_z
end

nx(model::DiscreteModel) = length(model.x0)
nq(model::DiscreteModel) = length(model.q0)
np(model::DiscreteModel) = size(model.dq, 1)
nu(model::DiscreteModel) = size(model.eq, 2)
ny(model::DiscreteModel) = length(model.y0)
nn(model::DiscreteModel) = size(model.fq, 2)

function steadystate(model::DiscreteModel, u=zeros(nu(model)))
    IA_LU = lufact(eye(nx(model))-model.a)
    steady_nl_eq_func = eval(quote
        (res, J, scratch, z) ->
            let pfull=scratch[1], Jp=scratch[2],
                q=$(zeros(nq(model))), Jq=$(zeros(nn(model), nq(model))),
                fq=$(model.pexp*model.dq/IA_LU*model.c + model.fq)
                $(model.nonlinear_eq)
                return nothing
            end
    end)
    steady_nleq = ParametricNonLinEq(steady_nl_eq_func, nn(model), nq(model))
    steady_solver = HomotopySolver{SimpleSolver}(steady_nleq, zeros(nq(model)),
                                                 zeros(nn(model)))
    set_resabs2tol!(steady_solver, 1e-30)
    steady_q0 = model.q0 + model.pexp*(model.dq/IA_LU*model.b + model.eq)*u +
                model.pexp*model.dq/IA_LU*model.x0
    steady_z = solve(steady_solver, steady_q0)
    if !hasconverged(steady_solver)
        error("Failed to find steady state solution")
    end
    return IA_LU\(model.b*u + model.c*steady_z + model.x0)
end

function steadystate!(model::DiscreteModel, u=zeros(nu(model)))
    x_steady = steadystate(model, u)
    copy!(model.x, x_steady)
    return x_steady
end

function run!(model::DiscreteModel, u::AbstractMatrix{Float64})
    if size(u, 1) ≠ nu(model)
        throw(DimensionMismatch("input matrix has $(size(u,1)) rows, but model requires $(nu(model)) inputs"))
    end
    y = Array(Float64, ny(model), size(u)[2])
    ucur = Array(Float64, nu(model))
    p = Array(Float64, np(model))
    ycur = Array(Float64, ny(model))
    xnew = Array(Float64, nx(model))
    @showprogress 1 "Running model: " for n = 1:size(u)[2]
        # copy!(p, model.dq * model.x + model.eq * u[:,n])
        copy!(ucur, 1, u, (n-1)*nu(model)+1, nu(model))
        BLAS.gemv!('N', 1., model.dq, model.x, 0., p)
        BLAS.gemv!('N', 1., model.eq, ucur, 1., p)
        z = solve(model.solver, p)
        if ~hasconverged(model.solver)
            if all(isfinite, z)
                warn("Failed to converge while solving non-linear equation.")
            else
                error("Failed to converge while solving non-linear equation, got non-finite result.")
            end
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
    return y
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
    gensolve(a, b, spzeros(size(a)[2], size(b)[2]), speye(size(a)[2]), thresh)

consecranges(lengths) = map(range, cumsum([1; lengths[1:end-1]]), lengths)

matsplit(m, rowsizes, colsizes=[size(m)[2]]) =
    [m[rs, cs] for rs in consecranges(rowsizes), cs in consecranges(colsizes)]

    if VERSION < v"0.5.0"
        # this is deprecated starting in Julia 0.6
        normsquared(x) = sumabs2(x)
    else
        # prior to Julia 0.5, this fails for empty x and is slow
        normsquared(x) = sum(abs2, x)
    end

if VERSION >= v"0.6.0-dev.1553"
    # workaround for JuliaLang/julia#19595
    import Base: +, -
    +(x1::SparseMatrixCSC{Real}, x2::SparseMatrixCSC{Real}) = broadcast!(+, similar(x1), x1, x2)
    -(x1::SparseMatrixCSC{Real}, x2::SparseMatrixCSC{Real}) = broadcast!(-, similar(x1), x1, x2)
end

end # module
