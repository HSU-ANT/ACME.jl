# Copyright 2015, 2016, 2017 Martin Holters
# See accompanying license file.

__precompile__()

module ACME

export Circuit, add!, connect!, DiscreteModel, run!, steadystate, steadystate!,
    ModelRunner

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
const Pin = Tuple{Element, Vector{Tuple{Int,Int}}}

# allow elem[:pin] notation to get an elements pin
getindex(e::Element, p) = (e, e.pins[Symbol(p)])

include("elements.jl")

const Net = Vector{Tuple{Int,Int}} # each net is a list of branch/polarity pairs

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
    @eval ($mat)(c::Circuit) =
         blkdiag(spzeros(Rational{BigInt}, 0, 0),
                 [convert(SparseMatrixCSC{Rational{BigInt}}, elem.$mat) for elem in c.elements]...)
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

    @compat function (::Type{DiscreteModel{Solver}}){Solver}(circ::Circuit, t::Float64)
        Base.depwarn("DiscreteModel{Solver}(circ, t) is deprecated, use DiscreteModel(circ, t, Solver) instead.",
                     :DiscreteModel)
        DiscreteModel(circ, t, Solver)
    end

    @compat function (::Type{DiscreteModel{Solver}}){Solver}(mats::Dict{Symbol}, nonlinear_eq::Expr, solver::Solver)
        model = new{Solver}()

        for mat in (:a, :b, :c, :pexp, :dq, :eq, :fq, :dy, :ey, :fy, :x0, :q0, :y0)
            setfield!(model, mat, convert(fieldtype(typeof(model), mat), mats[mat]))
        end

        model.nonlinear_eq = nonlinear_eq
        model.solver = solver
        model.x = zeros(nx(model))
        return model
    end
end

function DiscreteModel{Solver}(circ::Circuit, t::Real, ::Type{Solver}=HomotopySolver{CachingSolver{SimpleSolver}})
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

    q0 = convert(Vector{Float64}, mats[:q0])
    pexp = convert(Matrix{Float64}, mats[:pexp])
    fq = convert(Matrix{Float64}, mats[:fq])

    init_z = initial_solution(model_nonlinear_eq, q0, fq)
    nonlinear_eq_func = eval(quote
        if VERSION < v"0.5.0-dev+2396"
            # wrap up in named function because anonymous functions are slow
            # in old Julia versions
            function $(gensym())(res, J, scratch, z)
                pfull=scratch[1]
                Jq=scratch[2]
                q=$(zeros(model_nq))
                fq=$fq
                $(model_nonlinear_eq)
                return nothing
            end
        else
            (res, J, scratch, z) ->
                let pfull=scratch[1], Jq=scratch[2], q=$(zeros(model_nq)),
                    fq=$fq
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
    end)
    nonlinear_eq_calc_Jp = eval(quote
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
    end)
    solver = eval(:($Solver(ParametricNonLinEq($nonlinear_eq_func, $nonlinear_eq_set_p,
                                       $nonlinear_eq_calc_Jp,
                                       (zeros($model_nq), zeros($model_nn, $model_nq)),
                                       $model_nn, $model_np),
                    zeros($model_np), $init_z)))
    return DiscreteModel{typeof(solver)}(mats, model_nonlinear_eq, solver)
end

function model_matrices(circ::Circuit, t::Rational{BigInt})
    lhs = convert(SparseMatrixCSC{Rational{BigInt}},
                  sparse([mv(circ) mi(circ) mxd(circ)//t+mx(circ)//2 mq(circ);
                   blkdiag(topomat(circ)...) spzeros(nb(circ), nx(circ) + nq(circ))]))
    rhs = convert(SparseMatrixCSC{Rational{BigInt}},
                  sparse([u0(circ) mu(circ) mxd(circ)//t-mx(circ)//2;
                          spzeros(nb(circ), 1+nu(circ)+nx(circ))]))
    x, f = map(full, gensolve(lhs, rhs))

    rowsizes = [nb(circ); nb(circ); nx(circ); nq(circ)]
    res = Dict{Symbol,Array}(zip([:fv; :fi; :c; :fq], matsplit(f, rowsizes)))

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

function reduce_pdims!(mats::Dict)
    # decompose [dq_full eq_full] into pexp*[dq eq] with [dq eq] having minimum
    # number of rows
    pexp, dqeq = rank_factorize(sparse([mats[:dq_full] mats[:eq_full]]))
    mats[:pexp] = pexp
    colsizes = [size(mats[m], 2) for m in [:dq_full, :eq_full]]
    mats[:dq], mats[:eq] = matsplit(dqeq, [size(dqeq, 1)], colsizes)

    # project pexp onto the orthogonal complement of the column space of Fq
    fq_pinv = gensolve(sparse(mats[:fq]'*mats[:fq]), mats[:fq]')[1]
    pexp = pexp - mats[:fq]*fq_pinv*pexp
    # if the new pexp has lower rank, update
    pexp, f = rank_factorize(sparse(pexp))
    if size(pexp, 2) < size(mats[:pexp], 2)
        mats[:pexp] = pexp
        mats[:dq] = f * mats[:dq]
        mats[:eq] = f * mats[:eq]
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
nq(model::DiscreteModel) = length(model.q0)
np(model::DiscreteModel) = size(model.dq, 1)
nu(model::DiscreteModel) = size(model.eq, 2)
ny(model::DiscreteModel) = length(model.y0)
nn(model::DiscreteModel) = size(model.fq, 2)

function steadystate(model::DiscreteModel, u=zeros(nu(model)))
    IA_LU = lufact(eye(nx(model))-model.a)
    steady_q0 = model.q0 + model.pexp*(model.dq/IA_LU*model.b + model.eq)*u +
    model.pexp*model.dq/IA_LU*model.x0
    steady_z = eval(quote
        steady_nl_eq_func = (res, J, scratch, z) ->
            let pfull=scratch[1], Jp=scratch[2],
                q=$(zeros(nq(model))), Jq=$(zeros(nn(model), nq(model))),
                fq=$(model.pexp*model.dq/IA_LU*model.c + model.fq)
                $(model.nonlinear_eq)
                return nothing
            end
        steady_nleq = ParametricNonLinEq(steady_nl_eq_func, nn($model), nq($model))
        steady_solver = HomotopySolver{SimpleSolver}(steady_nleq, zeros(nq($model)),
                                                     zeros(nn($model)))
        set_resabstol!(steady_solver, 1e-15)
        steady_z = solve(steady_solver, $steady_q0)
        if !hasconverged(steady_solver)
            error("Failed to find steady state solution")
        end
        return steady_z
    end)
    return IA_LU\(model.b*u + model.c*steady_z + model.x0)
end

function steadystate!(model::DiscreteModel, u=zeros(nu(model)))
    x_steady = steadystate(model, u)
    copy!(model.x, x_steady)
    return x_steady
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

immutable ModelRunner{Model<:DiscreteModel,ShowProgress}
    model::Model
    ucur::Vector{Float64}
    p::Vector{Float64}
    ycur::Vector{Float64}
    xnew::Vector{Float64}
    @compat function (::Type{ModelRunner{Model,ShowProgress}}){Model<:DiscreteModel,ShowProgress}(model::Model)
        ucur = Array{Float64,1}(nu(model))
        p = Array{Float64,1}(np(model))
        ycur = Array{Float64,1}(ny(model))
        xnew = Array{Float64,1}(nx(model))
        return new{Model,ShowProgress}(model, ucur, p, ycur, xnew)
    end
end

ModelRunner{Model<:DiscreteModel}(model::Model) = ModelRunner{Model,true}(model)
ModelRunner{Model<:DiscreteModel,ShowProgress}(model::Model, ::Val{ShowProgress}) =
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
ModelRunner{Model<:DiscreteModel}(model::Model, showprogress::Bool) =
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
function run!{Model<:DiscreteModel}(runner::ModelRunner{Model,true},
                                    y::AbstractMatrix{Float64},
                                    u::AbstractMatrix{Float64})
    checkiosizes(runner, u, y)
    @showprogress "Running model: " for n = 1:size(u, 2)
        step!(runner, y, u, n)
    end
end

function run!{Model<:DiscreteModel}(runner::ModelRunner{Model,false},
                                    y::AbstractMatrix{Float64},
                                    u::AbstractMatrix{Float64})
    checkiosizes(runner, u, y)
    for n = 1:size(u, 2)
        step!(runner, y, u, n)
    end
end

function step!(runner::ModelRunner, y::AbstractMatrix{Float64}, u::AbstractMatrix{Float64}, n)
    model = runner.model
    ucur = runner.ucur
    p = runner.p
    ycur = runner.ycur
    xnew = runner.xnew
    # copy!(p, model.dq * model.x + model.eq * u[:,n])
    copy!(ucur, 1, u, (n-1)*nu(model)+1, nu(model))
    if size(model.dq, 2) == 0
        fill!(p, 0.0)
    else
        BLAS.gemv!('N', 1., model.dq, model.x, 0., p)
    end
    BLAS.gemv!('N', 1., model.eq, ucur, 1., p)
    z = solve(model.solver, p)
    if !hasconverged(model.solver)
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

end # module
