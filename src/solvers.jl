# Copyright 2016, 2017 Martin Holters
# See accompanying license file.

export SimpleSolver, HomotopySolver, CachingSolver
import Base.copy!

#struct ParametricNonLinEq{F_eval<:Function,F_setp<:Function,F_calcjp<:Function,Scratch}
eval(Expr(:type, false, :(ParametricNonLinEq{F_eval<:Function,F_setp<:Function,F_calcjp<:Function,Scratch}), quote
    func::F_eval
    set_p::F_setp
    calc_Jp::F_calcjp
    res::Vector{Float64}
    Jp::Matrix{Float64}
    J::Matrix{Float64}
    scratch::Scratch
    @compat function (::Type{ParametricNonLinEq{F_eval,F_setp,F_calcjp,Scratch}}){
            F_eval<:Function,F_setp<:Function,F_calcjp<:Function,Scratch
        }(
            func::F_eval, set_p::F_setp, calc_Jp::F_calcjp, scratch::Scratch,
            nn::Integer, np::Integer
        )
        res = zeros(nn)
        Jp = zeros(nn, np)
        J = zeros(nn, nn)
        return new{F_eval,F_setp,F_calcjp,Scratch}(func, set_p, calc_Jp, res, Jp, J, scratch)
    end
end))
ParametricNonLinEq{F_eval<:Function,F_setp<:Function,F_calcjp<:Function,
                   Scratch}(func::F_eval, set_p::F_setp, calc_Jp::F_calcjp,
                            scratch::Scratch, nn::Integer, np::Integer) =
    ParametricNonLinEq{F_eval,F_setp,F_calcjp,Scratch}(func, set_p, calc_Jp, scratch, nn, np)
ParametricNonLinEq(func::Function, nn::Integer, np::Integer) =
    ParametricNonLinEq(func, default_set_p, default_calc_Jp,
                       (zeros(np), zeros(nn, np)), nn, np)

default_set_p(scratch, p) = (copy!(scratch[1], p); nothing)
default_calc_Jp(scratch, Jp) = (copy!(Jp, scratch[2]); nothing)

nn(nleq::ParametricNonLinEq) = length(nleq.res)
np(nleq::ParametricNonLinEq) = size(nleq.Jp, 2)

set_p!(nleq::ParametricNonLinEq, p) = nleq.set_p(nleq.scratch, p)
calc_Jp!(nleq::ParametricNonLinEq) = nleq.calc_Jp(nleq.scratch, nleq.Jp)
evaluate!(nleq::ParametricNonLinEq, z) =
    nleq.func(nleq.res, nleq.J, nleq.scratch, z)

#struct LinearSolver
eval(Expr(:type, false, :LinearSolver, quote
    factors::Matrix{Float64}
    ipiv::Vector{Base.LinAlg.BlasInt}
    info::typeof(Ref{Base.LinAlg.BlasInt}(0))
    function LinearSolver(n::Int)
        new(zeros(n, n), zeros(Base.LinAlg.BlasInt, n), Ref{Base.LinAlg.BlasInt}(0))
    end
end))

function setlhs!(solver::LinearSolver, A::Matrix{Float64})
    m, n = size(solver.factors)
    if (m, n) ≠ size(A)
        throw(DimensionMismatch("matrix has size $(size(A)), but must have size $(size(solver.factors))"))
    end
    copy!(solver.factors, A)
    lda  = max(1, m)
    ccall((Compat.@blasfunc(dgetrf_), Base.LinAlg.LAPACK.liblapack), Void,
          (Ptr{Base.LinAlg.BlasInt}, Ptr{Base.LinAlg.BlasInt}, Ptr{Float64},
           Ptr{Base.LinAlg.BlasInt}, Ptr{Base.LinAlg.BlasInt}, Ptr{Base.LinAlg.BlasInt}),
          &m, &n, solver.factors, &lda, solver.ipiv, solver.info)
    return solver.info[] == 0
end

function solve!(solver::LinearSolver, x::Vector{Float64}, b::Vector{Float64})
    n = size(solver.factors, 2)
    if n ≠ length(x)
        throw(DimensionMismatch("x has length $(length(x)), but needs $n"))
    end
    if x !== b
        if n ≠ length(b)
            throw(DimensionMismatch("b has length $(length(b)), but needs $n"))
        end
        copy!(x, b)
    end
    Base.LinAlg.chkstride1(solver.factors, x, solver.ipiv)
    ccall((Compat.@blasfunc(dgetrs_), Base.LinAlg.LAPACK.liblapack), Void,
          (Ptr{UInt8}, Ptr{Base.LinAlg.BlasInt}, Ptr{Base.LinAlg.BlasInt}, Ptr{Float64}, Ptr{Base.LinAlg.BlasInt},
           Ptr{Base.LinAlg.BlasInt}, Ptr{Float64}, Ptr{Base.LinAlg.BlasInt}, Ptr{Base.LinAlg.BlasInt}),
          &'N', &n, &1, solver.factors, &max(1,stride(solver.factors,2)), solver.ipiv, x, &max(1,stride(x,2)), solver.info)
    if solver.info[] ≠ 0
        throw(Base.LinAlg.LAPACKException(solver.info[]))
    end
    return nothing
end

function copy!(dest::LinearSolver, src::LinearSolver)
    copy!(dest.factors, src.factors)
    copy!(dest.ipiv, src.ipiv)
    dest.info[] = src.info[]
end

#mutable struct SimpleSolver{NLEQ<:ParametricNonLinEq}
eval(Expr(:type, true, :(SimpleSolver{NLEQ<:ParametricNonLinEq}), quote
    nleq::NLEQ
    z::Vector{Float64}
    linsolver::LinearSolver
    last_z::Vector{Float64}
    last_p::Vector{Float64}
    last_Jp::Matrix{Float64}
    last_linsolver::LinearSolver
    iters::Int
    resmaxabs::Float64
    tol::Float64
    tmp_nn::Vector{Float64}
    tmp_np::Vector{Float64}
    @compat function (::Type{SimpleSolver{NLEQ}}){NLEQ<:ParametricNonLinEq}(
            nleq::NLEQ, initial_p::Vector{Float64}, initial_z::Vector{Float64})
        z = zeros(nn(nleq))
        linsolver = LinearSolver(nn(nleq))
        last_z = zeros(nn(nleq))
        last_p = zeros(np(nleq))
        last_Jp = zeros(nn(nleq), np(nleq))
        last_linsolver = LinearSolver(nn(nleq))
        tmp_nn = zeros(nn(nleq))
        tmp_np = zeros(np(nleq))
        solver = new{NLEQ}(nleq, z, linsolver, last_z, last_p, last_Jp,
            last_linsolver, 0, 0.0, 1e-10, tmp_nn, tmp_np)
        set_extrapolation_origin(solver, initial_p, initial_z)
        return solver
    end
end))
@doc """
    SimpleSolver

The `SimpleSolver` is the simplest available solver. It uses Newton iteration
which features fast local convergence, but makes no guarantees about global
convergence. The initial solution of the iteration is obtained by extrapolating
the last solution found (or another solution provided externally) using the
available Jacobians. Due to the missing global convergence, the `SimpleSolver`
is rarely useful as such.
""" SimpleSolver

SimpleSolver{NLEQ<:ParametricNonLinEq}(nleq::NLEQ, initial_p::Vector{Float64},
                                       initial_z::Vector{Float64}) =
    SimpleSolver{NLEQ}(nleq, initial_p, initial_z)

set_resabstol!(solver::SimpleSolver, tol) = solver.tol = tol

function set_extrapolation_origin(solver::SimpleSolver, p, z)
    set_p!(solver.nleq, p)
    evaluate!(solver.nleq, z)
    setlhs!(solver.linsolver, solver.nleq.J)
    calc_Jp!(solver.nleq)
    set_extrapolation_origin(solver, p, z, solver.nleq.Jp, solver.linsolver)
end

function set_extrapolation_origin(solver::SimpleSolver, p, z, Jp, linsolver)
    copy!(solver.last_linsolver, linsolver)
    copy!(solver.last_Jp, Jp)
    copy!(solver.last_p, p)
    copy!(solver.last_z, z)
end

get_extrapolation_origin(solver::SimpleSolver) = solver.last_p, solver.last_z

hasconverged(solver::SimpleSolver) = solver.resmaxabs < solver.tol

needediterations(solver::SimpleSolver) = solver.iters

function solve(solver::SimpleSolver, p::AbstractVector{Float64}, maxiter=500)
    set_p!(solver.nleq, p)
    #solver.z = solver.last_z - solver.last_J\(solver.last_Jp * (p-solver.last_p))
    copy!(solver.tmp_np, p)
    BLAS.axpy!(-1.0, solver.last_p, solver.tmp_np)
    BLAS.gemv!('N', 1.,solver.last_Jp, solver.tmp_np, 0., solver.tmp_nn)
    solve!(solver.last_linsolver, solver.tmp_nn, solver.tmp_nn)
    copy!(solver.z, solver.last_z)
    BLAS.axpy!(-1.0, solver.tmp_nn, solver.z)

    for solver.iters=1:maxiter
        evaluate!(solver.nleq, solver.z)
        solver.resmaxabs = isempty(solver.nleq.res) ? 0.0 : maximum(abs, solver.nleq.res)
        if !isfinite(solver.resmaxabs) || !all(isfinite, solver.nleq.J)
            return solver.z
        end
        if !setlhs!(solver.linsolver, solver.nleq.J) # J was singular
            return solver.z
        end
        hasconverged(solver) && break
        #solver.z -= solver.nleq.J\solver.nleq.res
        solve!(solver.linsolver, solver.tmp_nn, solver.nleq.res)
        BLAS.axpy!(-1.0, solver.tmp_nn, solver.z)
    end
    if hasconverged(solver)
        calc_Jp!(solver.nleq)
        set_extrapolation_origin(solver, p, solver.z, solver.nleq.Jp, solver.linsolver)
    end
    return solver.z
end

#mutable struct HomotopySolver{BaseSolver}
eval(Expr(:type, true, :(HomotopySolver{BaseSolver}), quote
    basesolver::BaseSolver
    start_p::Vector{Float64}
    pa::Vector{Float64}
    iters::Int
    @compat function (::Type{HomotopySolver{BaseSolver}}){BaseSolver}(
            basesolver::BaseSolver, np::Integer)
        return new{BaseSolver}(basesolver, zeros(np), zeros(np), 0)
    end
    @compat function (::Type{HomotopySolver{BaseSolver}}){BaseSolver}(
            nleq::ParametricNonLinEq, initial_p::Vector{Float64},
            initial_z::Vector{Float64})
        basesolver = BaseSolver(nleq, initial_p, initial_z)
        return HomotopySolver{typeof(basesolver)}(basesolver, np(nleq))
    end
end))
@doc """
    HomotopySolver{BaseSolver}

The `HomotopySolver` extends an existing solver (provided as the type parameter)
by applying homotopy to (at least theoretically) ensure global convergence. It
can be combined with the `SimpleSolver` as `HomotopySolver{SimpleSolver}` to
obtain a useful Newton homtopy solver with generally good convergence
properties.
""" -> HomotopySolver

set_resabstol!(solver::HomotopySolver, tol) =
    set_resabstol!(solver.basesolver, tol)

set_extrapolation_origin(solver::HomotopySolver, p, z) =
    set_extrapolation_origin(solver.basesolver, p, z)

function solve(solver::HomotopySolver, p)
    z = solve(solver.basesolver, p)
    solver.iters = needediterations(solver.basesolver)
    if !hasconverged(solver)
        a = 0.5
        best_a = 0.0
        copy!(solver.start_p, get_extrapolation_origin(solver.basesolver)[1])
        while best_a < 1
            # copy!(solver.pa = (1-a) * solver.start_p + a * p)
            copy!(solver.pa, solver.start_p)
            LinAlg.scale!(1-a, solver.pa)
            LinAlg.axpy!(a, p, solver.pa)
            z = solve(solver.basesolver, solver.pa)
            solver.iters += needediterations(solver.basesolver)
            if hasconverged(solver)
                best_a = a
                a = 1.0
            else
                new_a = (a + best_a) / 2
                if !(best_a < new_a < a)
                    # no floating point value inbetween best_a and a
                    break
                end
                a = new_a
            end
        end
    end
    return z
end

hasconverged(solver::HomotopySolver) = hasconverged(solver.basesolver)
needediterations(solver::HomotopySolver) = solver.iters

#mutable struct CachingSolver{BaseSolver}
eval(Expr(:type, true, :(CachingSolver{BaseSolver}), quote
    basesolver::BaseSolver
    ps_tree::KDTree{Vector{Float64}, Matrix{Float64}}
    zs::Matrix{Float64}
    num_ps::Int
    new_count::Int
    new_count_limit::Int
    alts::Alts{Float64}
    @compat function (::Type{CachingSolver{BaseSolver}}){BaseSolver}(basesolver::BaseSolver,
            initial_p::Vector{Float64}, initial_z::Vector{Float64}, nn::Integer)
         ps_tree = KDTree(hcat(initial_p))
         zs = reshape(copy(initial_z), nn, 1)
         alts = Alts(initial_p)
         return new{BaseSolver}(basesolver, ps_tree, zs, 1, 0, 2, alts)
    end
    @compat function (::Type{CachingSolver{BaseSolver}}){BaseSolver}(nleq::ParametricNonLinEq,
            initial_p::Vector{Float64}, initial_z::Vector{Float64})
        basesolver = BaseSolver(nleq, initial_p, initial_z)
        return CachingSolver{typeof(basesolver)}(basesolver, initial_p, initial_z, nn(nleq))
    end
end))
@doc """
    CachingSolver{BaseSolver}

The `CachingSolver` extends an existing solver (provided as the type parameter)
by storing found solutions in a k-d tree to use as initial solutions in the
future. Whenever the underlying solver needs more than a preset number of
iterations (defaults to five), the solution will be stored. Storing new
solutions is a relatively expensive operation, so until the stored solutions
suffice to ensure convergence in few iterations throughout, use of a
`CachingSolver` may actually slow things down.

See [M. Holters, U. Zölzer, "A k-d Tree Based Solution Cache for the Non-linear
Equation of Circuit Simulations"](http://www.eurasip.org/Proceedings/Eusipco/Eusipco2016/papers/1570255150.pdf)
for a more detailed discussion.
""" -> CachingSolver

set_resabstol!(solver::CachingSolver, tol) =
    set_resabstol!(solver.basesolver, tol)

hasconverged(solver::CachingSolver) = hasconverged(solver.basesolver)
needediterations(solver::CachingSolver) = needediterations(solver.basesolver)

function solve(solver::CachingSolver, p)
    origin_p = get_extrapolation_origin(solver.basesolver)[1]
    best_diff = 0.0
    for i in eachindex(origin_p)
        best_diff += abs2(p[i] - origin_p[i])
    end
    idx = 0
    for i in (solver.num_ps-solver.new_count+1):solver.num_ps
        diff = 0.
        for j in 1:size(solver.ps_tree.ps, 1)
            diff += abs2(solver.ps_tree.ps[j,i] - p[j])
        end
        if diff < best_diff
            best_diff = diff
            idx = i
        end
    end

    init!(solver.alts, best_diff, idx)
    idx = indnearest(solver.ps_tree, p, solver.alts)

    if idx ≠ 0
        set_extrapolation_origin(solver.basesolver,
                                 solver.ps_tree.ps[:,idx], solver.zs[:,idx])
    end

    z = solve(solver.basesolver, p)
    if needediterations(solver.basesolver) > 5 && hasconverged(solver.basesolver)
        solver.num_ps += 1
        if solver.num_ps > size(solver.ps_tree.ps, 2)
            solver.ps_tree.ps =
                copy!(zeros(size(solver.ps_tree.ps, 1), 2solver.num_ps),
                      solver.ps_tree.ps)
            solver.zs =
                copy!(zeros(size(solver.zs, 1), 2solver.num_ps), solver.zs)
        end
        solver.ps_tree.ps[:,solver.num_ps] = p
        solver.zs[:,solver.num_ps] = z
        solver.new_count += 1
    end
    if solver.new_count > 0
        solver.new_count_limit -= 1
    end
    if solver.new_count > solver.new_count_limit
        solver.ps_tree = KDTree(solver.ps_tree.ps, solver.num_ps)
        solver.new_count = 0
        solver.new_count_limit = 2size(solver.ps_tree.ps, 2)
    end
    return z
end

get_extrapolation_origin(solver::CachingSolver) =
    get_extrapolation_origin(solver.basesolver)


function set_resabs2tol!(solver, tol)
    Base.depwarn(string("set_resabs2tol!(solver, tol) is deprecated, use set_resabstol!(solver, sqrt(tol)) instead."),
                 :set_resabs2tol!)
    set_resabstol!(solver, sqrt(tol))
end
