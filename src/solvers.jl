# Copyright 2016, 2017 Martin Holters
# See accompanying license file.

export SimpleSolver, HomotopySolver, CachingSolver
import Base.copy!

#struct ParametricNonLinEq{F_eval<:Function,F_setp<:Function,F_calcjp<:Function,Scratch}
@struct ParametricNonLinEq{F_eval<:Function,F_setp<:Function,F_calcjp<:Function,Scratch} begin
    func::F_eval
    set_p::F_setp
    calc_Jp::F_calcjp
    res::Vector{Float64}
    Jp::Matrix{Float64}
    J::Matrix{Float64}
    scratch::Scratch
    @pfunction (::Type{ParametricNonLinEq{F_eval,F_setp,F_calcjp,Scratch}})(
            func::F_eval, set_p::F_setp, calc_Jp::F_calcjp, scratch::Scratch,
            nn::Integer, np::Integer
        ) [F_eval<:Function,F_setp<:Function,F_calcjp<:Function,Scratch] begin
        res = zeros(nn)
        Jp = zeros(nn, np)
        J = zeros(nn, nn)
        return new{F_eval,F_setp,F_calcjp,Scratch}(func, set_p, calc_Jp, res, Jp, J, scratch)
    end
end
@pfunction ParametricNonLinEq(func::F_eval, set_p::F_setp, calc_Jp::F_calcjp,
        scratch::Scratch, nn::Integer, np::Integer) [F_eval<:Function,
        F_setp<:Function,F_calcjp<:Function,Scratch] begin
    ParametricNonLinEq{F_eval,F_setp,F_calcjp,Scratch}(func, set_p, calc_Jp, scratch, nn, np)
end
ParametricNonLinEq(func::Function, nn::Integer, np::Integer) =
    ParametricNonLinEq(func, default_set_p, default_calc_Jp,
                       (zeros(np), zeros(nn, np)), nn, np)

default_set_p(scratch, p) = (copyto!(scratch[1], p); nothing)
default_calc_Jp(scratch, Jp) = (copyto!(Jp, scratch[2]); nothing)

nn(nleq::ParametricNonLinEq) = length(nleq.res)
np(nleq::ParametricNonLinEq) = size(nleq.Jp, 2)

set_p!(nleq::ParametricNonLinEq, p) = nleq.set_p(nleq.scratch, p)
calc_Jp!(nleq::ParametricNonLinEq) = nleq.calc_Jp(nleq.scratch, nleq.Jp)
evaluate!(nleq::ParametricNonLinEq, z) =
    nleq.func(nleq.res, nleq.J, nleq.scratch, z)

#struct LinearSolver
@struct LinearSolver begin
    factors::Matrix{Float64}
    ipiv::Vector{Int}
    function LinearSolver(n::Int)
        new(zeros(n, n), zeros(Int, n))
    end
end

function setlhs!(solver::LinearSolver, A::Matrix{Float64})
    m, n = size(solver.factors)
    if (m, n) ≠ size(A)
        throw(DimensionMismatch("matrix has size $(size(A)), but must have size $(size(solver.factors))"))
    end
    copyto!(solver.factors, A)

    # based on Julia's generic_lufact!, but storing inverses on the diagonal;
    # faster than calling out to dgetrf for sizes up to about 60×60
    factors = solver.factors
    minmn = min(m,n)
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            amax = 0.0
            for i = k:m
                absi = abs(factors[i,k])
                if absi > amax
                    kp = i
                    amax = absi
                end
            end
            solver.ipiv[k] = kp
            if factors[kp,k] != 0.0
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = factors[k,i]
                        factors[k,i] = factors[kp,i]
                        factors[kp,i] = tmp
                    end
                end
                # Scale first column
                fkkinv = factors[k,k] = inv(factors[k,k])
                for i = k+1:m
                    factors[i,k] *= fkkinv
                end
            else
                return false
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    factors[i,j] -= factors[i,k]*factors[k,j]
                end
            end
        end
    end
    return true
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
        copyto!(x, b)
    end

    # native Julia implementation seems to be faster than dgetrs up to about
    # n=45 (and not slower up to about n=70)
    @inbounds for i in 1:n
        x[i], x[solver.ipiv[i]] = x[solver.ipiv[i]], x[i]
    end
    # taken from Julia's naivesub!(::UnitLowerTriangular, ...)
    @inbounds for j in 1:n
        xj = x[j]
        for i in j+1:n
            x[i] -= solver.factors[i,j] * xj
        end
    end
    # based on Julia's naivesub!(::UpperTriangular, ...), but with factors[j,j]
    # holding inverses
    @inbounds for j in n:-1:1
        xj = x[j] = solver.factors[j,j] * x[j]
        for i in 1:j-1
            x[i] -= solver.factors[i,j] * xj
        end
    end

    return nothing
end

function copy!(dest::LinearSolver, src::LinearSolver)
    copyto!(dest.factors, src.factors)
    copyto!(dest.ipiv, src.ipiv)
end

"""
    SimpleSolver

The `SimpleSolver` is the simplest available solver. It uses Newton iteration
which features fast local convergence, but makes no guarantees about global
convergence. The initial solution of the iteration is obtained by extrapolating
the last solution found (or another solution provided externally) using the
available Jacobians. Due to the missing global convergence, the `SimpleSolver`
is rarely useful as such.
"""
#mutable struct SimpleSolver{NLEQ<:ParametricNonLinEq}
@mutable_struct SimpleSolver{NLEQ<:ParametricNonLinEq} begin
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
    @pfunction (::Type{SimpleSolver{NLEQ}})(
            nleq::NLEQ, initial_p::Vector{Float64}, initial_z::Vector{Float64}) [NLEQ<:ParametricNonLinEq] begin
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
end

@pfunction SimpleSolver(nleq::NLEQ, initial_p::Vector{Float64},
                        initial_z::Vector{Float64}) [NLEQ<:ParametricNonLinEq] begin
    SimpleSolver{NLEQ}(nleq, initial_p, initial_z)
end

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
    copyto!(solver.last_Jp, Jp)
    copyto!(solver.last_p, p)
    copyto!(solver.last_z, z)
end

get_extrapolation_origin(solver::SimpleSolver) = solver.last_p, solver.last_z

get_extrapolation_jacobian(solver::SimpleSolver) =
    -solver.nleq.J \ solver.nleq.Jp

hasconverged(solver::SimpleSolver) = solver.resmaxabs < solver.tol

needediterations(solver::SimpleSolver) = solver.iters

function solve(solver::SimpleSolver, p::AbstractVector{Float64}, maxiter=500)
    set_p!(solver.nleq, p)
    #solver.z = solver.last_z - solver.last_J\(solver.last_Jp * (p-solver.last_p))
    copyto!(solver.tmp_np, p)
    BLAS.axpy!(-1.0, solver.last_p, solver.tmp_np)
    BLAS.gemv!('N', 1.,solver.last_Jp, solver.tmp_np, 0., solver.tmp_nn)
    solve!(solver.last_linsolver, solver.tmp_nn, solver.tmp_nn)
    copyto!(solver.z, solver.last_z)
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

"""
    HomotopySolver{BaseSolver}

The `HomotopySolver` extends an existing solver (provided as the type parameter)
by applying homotopy to (at least theoretically) ensure global convergence. It
can be combined with the `SimpleSolver` as `HomotopySolver{SimpleSolver}` to
obtain a useful Newton homtopy solver with generally good convergence
properties.
"""
#mutable struct HomotopySolver{BaseSolver}
@mutable_struct HomotopySolver{BaseSolver} begin
    basesolver::BaseSolver
    start_p::Vector{Float64}
    pa::Vector{Float64}
    iters::Int
    @pfunction (::Type{HomotopySolver{BaseSolver}})(
            basesolver::BaseSolver, np::Integer) [BaseSolver] begin
        return new{BaseSolver}(basesolver, zeros(np), zeros(np), 0)
    end
    @pfunction (::Type{HomotopySolver{BaseSolver}})(
            nleq::ParametricNonLinEq, initial_p::Vector{Float64},
            initial_z::Vector{Float64}) [BaseSolver] begin
        basesolver = BaseSolver(nleq, initial_p, initial_z)
        return HomotopySolver{typeof(basesolver)}(basesolver, np(nleq))
    end
end

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
        copyto!(solver.start_p, get_extrapolation_origin(solver.basesolver)[1])
        while best_a < 1
            # copyto!(solver.pa, (1-a) * solver.start_p + a * p)
            copyto!(solver.pa, solver.start_p)
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

get_extrapolation_jacobian(solver::HomotopySolver) =
    get_extrapolation_jacobian(solver.basesolver)

"""
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
"""
#mutable struct CachingSolver{BaseSolver}
@mutable_struct CachingSolver{BaseSolver} begin
    basesolver::BaseSolver
    ps_tree::KDTree{Vector{Float64}, Matrix{Float64}}
    zs::Matrix{Float64}
    num_ps::Int
    new_count::Int
    new_count_limit::Int
    alts::Alts{Float64}
    @pfunction (::Type{CachingSolver{BaseSolver}})(basesolver::BaseSolver,
            initial_p::Vector{Float64}, initial_z::Vector{Float64}, nn::Integer) [BaseSolver] begin
         ps_tree = KDTree(hcat(initial_p))
         zs = reshape(copy(initial_z), nn, 1)
         alts = Alts(initial_p)
         return new{BaseSolver}(basesolver, ps_tree, zs, 1, 0, 2, alts)
    end
    @pfunction (::Type{CachingSolver{BaseSolver}})(nleq::ParametricNonLinEq,
            initial_p::Vector{Float64}, initial_z::Vector{Float64}) [BaseSolver] begin
        basesolver = BaseSolver(nleq, initial_p, initial_z)
        return CachingSolver{typeof(basesolver)}(basesolver, initial_p, initial_z, nn(nleq))
    end
end

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
                copyto!(zeros(size(solver.ps_tree.ps, 1), 2solver.num_ps),
                      solver.ps_tree.ps)
            solver.zs =
                copyto!(zeros(size(solver.zs, 1), 2solver.num_ps), solver.zs)
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

set_extrapolation_origin(solver::CachingSolver, p, z) =
    set_extrapolation_origin(solver.basesolver, p, z)

get_extrapolation_jacobian(solver::CachingSolver) =
    get_extrapolation_jacobian(solver.basesolver)


function set_resabs2tol!(solver, tol)
    Base.depwarn(string("set_resabs2tol!(solver, tol) is deprecated, use set_resabstol!(solver, sqrt(tol)) instead."),
                 :set_resabs2tol!)
    set_resabstol!(solver, sqrt(tol))
end

function linearize(solver, p::AbstractVector{Float64})
    z = solve(solver, p)
    set_extrapolation_origin(solver, p, z)
    if !hasconverged(solver)
        throw(ArgumentError("Cannot linearize because no solution found at p=$p"))
    end
    return z, get_extrapolation_jacobian(solver)
end
