# Copyright 2015, 2016 Martin Holters
# See accompanying license file.

VERSION >= v"0.4.0-dev+6521" && __precompile__()

module ACME

export Circuit, add!, connect!, DiscreteModel, run!, steadystate, steadystate!

using ProgressMeter
using Compat

import Base.getindex

include("kdtree.jl")
include("solvers.jl")


type Element
  mv :: SparseMatrixCSC{Number,Int}
  mi :: SparseMatrixCSC{Number,Int}
  mx :: SparseMatrixCSC{Number,Int}
  mxd :: SparseMatrixCSC{Number,Int}
  mq :: SparseMatrixCSC{Number,Int}
  mu :: SparseMatrixCSC{Number,Int}
  u0 :: SparseMatrixCSC{Number,Int}
  pv :: SparseMatrixCSC{Number,Int}
  pi :: SparseMatrixCSC{Number,Int}
  px :: SparseMatrixCSC{Number,Int}
  pxd :: SparseMatrixCSC{Number,Int}
  pq :: SparseMatrixCSC{Number,Int}
  nonlinear_eq :: Expr
  pins :: Dict{Symbol, Vector{@compat Tuple{Int, Int}}}

  function Element(;args...)
    sizes = @compat Dict{Symbol,Int}(:n0 => 1)

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
      dict = @compat Dict{Symbol,Vector{@compat Tuple{Int, Int}}}()
      for i in 1:length(syms)
        branch = div(i+1, 2)
        polarity = 2mod(i, 2) - 1
        push!(get!(dict, symbol(syms[i]), []), (branch, polarity))
      end
      dict
    end
    make_pin_dict(dict::Dict) = dict

    const mat_dims =
        @compat  Dict( :mv => (:nl,:nb), :mi => (:nl,:nb), :mx => (:nl,:nx),
                       :mxd => (:nl,:nx), :mq => (:nl,:nq), :mu => (:nl,:nu),
                       :u0 => (:nl, :n0),
                       :pv => (:ny,:nb), :pi => (:ny,:nb), :px => (:ny,:nx),
                       :pxd => (:ny,:nx), :pq => (:ny,:nq) )

    elem = new()
    for (key, val) in args
      if haskey(mat_dims, key)
        val = sparse(hcat(val)) # turn val into a sparse matrix whatever it is
        update_sizes(val, mat_dims[key])
      elseif key == :pins
        val = make_pin_dict(val)
      end
      elem.(key) = val
    end
    for (m, ns) in mat_dims
      if !isdefined(elem, m)
        elem.(m) = spzeros(Int, get(sizes, ns[1], 0), get(sizes, ns[2], 0))
      end
    end
    if !isdefined(elem, :nonlinear_eq)
      elem.nonlinear_eq = Expr(:block)
    end
    if !isdefined(elem, :pins)
      elem.pins = make_pin_dict(map(string,1:2nb(elem)))
    end
    elem
  end
end

for (n,m) in @compat Dict(:nb => :mv, :nx => :mx, :nq => :mq, :nu => :mu)
  @eval ($n)(e::Element) = size(e.$m)[2]
end
nl(e::Element) = size(e.mv)[1]
ny(e::Element) = size(e.pv)[1]
nn(e::Element) = nb(e) + nx(e) + nq(e) - nl(e)

# a Pin combines an element with a branch/polarity list
typealias Pin @compat Tuple{Element, Vector{@compat Tuple{Int,Int}}}

# allow elem[:pin] notation to get an elements pin
getindex(e::Element, p::Symbol) = (e, e.pins[p])
getindex(e::Element, p::AbstractString) = getindex(e, symbol(p))
getindex(e::Element, p::Int) = getindex(e, string(p))

include("elements.jl")

typealias Net Vector{@compat Tuple{Int,Int}} # each net is a list of branch/polarity pairs

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
    @eval ($mat)(c::Circuit) = blkdiag(spzeros(Number, 0, 0),
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
        index_offsets = @compat Dict( :q => (col_offset,),
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

add!(c::Circuit, elems::Element...) = for elem in es add!(c, elem) end

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
    @assert isdefined(net)
    net
end

function netfor!(c::Circuit, name::Symbol)
    haskey(c.net_names, name) || push!(c.nets, get!(c.net_names, name, []))
    c.net_names[name]
end

function connect!(c::Circuit, pins::(@compat Union{Pin,Symbol})...)
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
    @assert all(abs(nonzeros(incidence)) .== 1)
    @assert all(sum(incidence, 1) .== 0)

    t = falses(size(incidence)[2]);

    row = 1;
    for col = 1:size(incidence)[2]
        rows = filter(r -> r ≥ row, find(incidence[:,col:col])) # !SV
        @assert length(rows) ≤ 2

        isempty(rows) && continue
        t[col] = true;

        if rows[1] ≠ row
            swaprows!(incidence, rows[1], row)
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
        model = new()

        mats = model_matrices(circ, t)
        for mat in [:a, :b, :c, :pexp, :dq, :eq, :fq, :dy, :ey, :fy]
            model.(mat)=full(mats[mat])
        end
        for vec in [:x0, :q0, :y0]
            model.(vec)=squeeze(full(mats[vec]), tuple((2:ndims(mats[vec]))...))
        end

        model.x = zeros(nx(model))

        model.nonlinear_eq = quote
            #copy!(q, q0 + pexp * p + fq * z)
            copy!(q, q0)
            BLAS.gemv!('N',1.,pexp,p,1.,q)
            BLAS.gemv!('N',1.,fq,z,1.,q)
            let J=Jq
                $(nonlinear_eq(circ))
            end
            #copy!(J, Jq*model.fq)
            BLAS.gemm!('N', 'N', 1., Jq, fq, 0., J)
            #copy!(Jp, Jq*pexp)
            BLAS.gemm!('N', 'N', 1., Jq, pexp, 0., Jp)
        end

        init_z = initial_solution(model.nonlinear_eq, model.q0, model.fq)
        nonlinear_eq_func = eval(quote
            if VERSION < v"0.5.0-dev+2396"
                # wrap up in named function because anonymous functions are slow
                # in old Julia versions
                function $(gensym())(res, J, Jp, p, z)
                    q0=$(model.q0)
                    pexp=$(model.pexp)
                    q=$(zeros(nq(model)))
                    Jq=$(zeros(nn(model), nq(model)))
                    fq=$(model.fq)
                    $(model.nonlinear_eq)
                    return nothing
                end
            else
                (res, J, Jp, p, z) ->
                    let q0=$(model.q0), pexp=$(model.pexp),
                        q=$(zeros(nq(model))),
                        Jq=$(zeros(nn(model), nq(model))), fq=$(model.fq)
                        $(model.nonlinear_eq)
                        return nothing
                    end
            end
        end)
        model.solver =
            Solver(ParametricNonLinEq(nonlinear_eq_func, nn(model), np(model)),
                   zeros(np(model)), init_z)
        return model
    end
end

DiscreteModel(circ::Circuit, t::Float64) =
    DiscreteModel{HomotopySolver{CachingSolver{SimpleSolver}}}(circ, t)

function model_matrices(circ::Circuit, t)
    x, f = gensolve([mv(circ) mi(circ) 1/t*mxd(circ)+0.5*mx(circ) mq(circ);
                     blkdiag(topomat(circ)...) spzeros(nb(circ), nx(circ) + nq(circ))],
                    [u0(circ) mu(circ) 1/t*mxd(circ)-0.5*mx(circ);
                     spzeros(nb(circ), 1+nu(circ)+nx(circ))])

    rowsizes = [nb(circ); nb(circ); nx(circ); nq(circ)]
    res = Dict{Symbol,AbstractArray}(zip([:fv; :fi; :c; :fq], matsplit(f, rowsizes)))

    # choose particular solution such that the rows corresponding to q are
    # column-wise orthogonal to the column space of fq (and hence have a column
    # space of minimal dimension)
    fq = res[:fq]
    x = x - full(f)/(full(fq)'*fq)*fq'*x[end-nq(circ)+1:end,:]

    merge!(res, Dict(zip([:v0 :ev :dv; :i0 :ei :di; :x0 :b :a; :q0 :eq_full :dq_full],
                         matsplit(x, rowsizes, [1; nu(circ); nx(circ)]))))

    dq_full = res[:dq_full]
    eq_full = res[:eq_full]
    if size(dq_full)[1] > 0
        # decompose [dq_full eq_full] into pexp*[dq eq] with [dq eq] having minimum
        # number of rows based on the QR factorization
        pexp, r, piv = qr([dq_full eq_full], Val{true})
        ref_err = vecnorm([dq_full eq_full][:,piv] - pexp*r)
        local rowcount
        for rowcount=size(r)[1]:-1:0
            err = vecnorm([dq_full eq_full][:,piv] - pexp[:,1:rowcount]*r[1:rowcount,:])
            if err > ref_err
                rowcount += 1
                break
            end
        end
        res[:dq], res[:eq] =
            matsplit(r[1:rowcount,sortperm(piv)], [rowcount], [nx(circ); nu(circ)])
        res[:pexp] = pexp[:,1:rowcount]
    else
        res[:dq] = zeros(0, nx(circ))
        res[:eq] = zeros(0, nu(circ))
        res[:pexp] = zeros(0, 0)
    end

    p = [pv(circ) pi(circ) 0.5*px(circ)+1/t*pxd(circ) pq(circ)]
    res[:dy] = p * x[:,2+nu(circ):end] + 0.5*px(circ)-1/t*pxd(circ)
    #          p * [dv; di; a;  dq_full] + 0.5*px(circ)-1/t*pxd(circ)
    res[:ey] = p * x[:,2:1+nu(circ)] # p * [ev; ei; b;  eq_full]
    res[:fy] = p * f                 # p * [fv; fi; c;  fq]
    res[:y0] = p * x[:,1]            # p * [v0; i0; x0; q0]

    return res
end

function initial_solution(nleq, q0, fq)
    # determine an initial solution with a homotopy solver that may vary q0
    # between 0 and the true q0 -> q0 takes the role of p
    nq, nn = size(fq)
    init_nl_eq_func = eval(quote
        (res, J, Jp, p, z) ->
            let q0=$(zeros(nq)), pexp=$(eye(nq)), q=$(zeros(nq)),
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
np(model::DiscreteModel) = size(model.dq)[1]
nu(model::DiscreteModel) = size(model.eq, 2)
ny(model::DiscreteModel) = length(model.y0)
nn(model::DiscreteModel) = size(model.fq)[2]

function steadystate(model::DiscreteModel, u=zeros(nu(model)))
    IA_LU = lufact(eye(nx(model))-model.a)
    steady_nl_eq_func = eval(quote
        (res, J, Jp, p, z) ->
            let q0=$(zeros(nq(model))), pexp=$(eye(nq(model))),
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
    y = Array(Float64, ny(model), size(u)[2])
    ucur = Array(Float64, nu(model))
    p = Array(Float64, np(model))
    ycur = Array(Float64, ny(model))
    xnew = Array(Float64, nx(model))
    @showprogress 1 "Running model: " for n = 1:size(u)[2]
        copy!(ucur, u[:,n])
        # copy!(p, model.dq * model.x + model.eq * u[:,n])
        BLAS.gemv!('N', 1., model.dq, model.x, 0., p)
        BLAS.gemv!('N', 1., model.eq, ucur, 1., p)
        z = solve(model.solver, p)
        if ~hasconverged(model.solver)
            if all(isfinite(z))
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
        y[:,n] = ycur
        #model.x = model.a * model.x + model.b * u[:,n] + model.c * z + model.x0
        copy!(xnew, model.x0)
        BLAS.gemv!('N', 1., model.a, model.x, 1., xnew)
        BLAS.gemv!('N', 1., model.b, ucur, 1.,xnew)
        BLAS.gemv!('N', 1., model.c, z, 1., xnew)
        copy!(model.x, xnew)
    end
    return y
end

function swaprows!(a::SparseMatrixCSC, row1, row2)
    # This sometimes gives a wrong result with julia 0.3.2:
    #    a[[row1, row2],:] = a[[row2, row1],:]
    i, j, v = findnz(a)
    row1idx = i.==row1
    row2idx = i.==row2
    i[row1idx] = row2
    i[row2idx] = row1
    a[:,:] = sparse(i, j, v, size(a)...)
end

# lines marked with !SV avoid creation of SparseVector by indexing with Ranges
# instead of Ints; a better way for cross-julia-version compatibilty would be
# nice; maybe Compat helps in the future...
function gensolve(a::SparseMatrixCSC, b, x, h, thresh=0.1)
    m = size(a)[1]
    t = sortperm(vec(sum(spones(a),2))) # row indexes in ascending order of nnz
    for i in 1:m
        ait = a[t[i:i],:] # ait is a row of the a matrix # !SV
        s = ait * h;
        if nnz(s) == 0
            continue
        end
        inz, jnz, nz_vals = findnz(s)
        nz_abs_vals = abs(nz_vals)
        jat = jnz[nz_abs_vals .≥ thresh*maximum(nz_abs_vals)] # cols above threshold
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

if VERSION < v"0.4.0-dev+2833"
    # work around a bug when indexing with empty ranges
    function matsplit(m, rowsizes, colsizes=[size(m)[2]])
        fixed_getindex(A, rs, cs) = isempty(rs) || isempty(cs) ?
            similar(A, length(rs), length(cs)) : return A[rs, cs]
        return [fixed_getindex(m, rs, cs) for rs in consecranges(rowsizes),
                                              cs in consecranges(colsizes)]
    end
else
    matsplit(m, rowsizes, colsizes=[size(m)[2]]) =
        [m[rs, cs] for rs in consecranges(rowsizes), cs in consecranges(colsizes)]
end

end # module
