# Copyright 2016 Martin Holters
# See accompanying license file.

type HomotopySolver{BaseSolver}
    basesolver::BaseSolver
    start_p::Vector{Float64}
    iters::Int
    function HomotopySolver(model::DiscreteModel)
        basesolver = BaseSolver(model)
        return new(basesolver, zeros(np(model)), 0)
    end
end

function solve(solver::HomotopySolver, p)
    z = solve(solver.basesolver, p)
    solver.iters = needediterations(solver.basesolver)
    if ~hasconverged(solver)
        a = 0.5
        best_a = 0.0
        copy!(solver.start_p, solver.basesolver.last_p)
        while best_a < 1 && a > 0
            pa = (1-a) * solver.start_p + a * p
            z = solve(solver.basesolver, pa)
            if hasconverged(solver)
                best_a = a
                a = 1.0
            else
                a = (a + best_a) / 2
            end
        end
    end
    return z
end

hasconverged(solver::HomotopySolver) = hasconverged(solver.basesolver)
needediterations(solver::HomotopySolver) = solver.iters


type CachingSolver{BaseSolver}
    basesolver::BaseSolver
    ps_tree::KDTree{Vector{Float64}, Matrix{Float64}}
    zs::Matrix{Float64}
    new_count::Int
    new_count_limit::Int
    iters::Int
    function CachingSolver(model::DiscreteModel)
        basesolver = BaseSolver(model)
        p = zeros(np(model))
        z = solve(basesolver, p, 2500)
        if ~hasconverged(basesolver)
            error("Failed to find initial solution.")
        end
        ps_tree = KDTree(zeros(np(model), 1))
        zs = reshape(copy(z), nn(model), 1)
        return new(basesolver, ps_tree, zs, 0, 2, 0)
    end
end

hasconverged(solver::CachingSolver) = hasconverged(solver.basesolver)
needediterations(solver::CachingSolver) = solver.iters

function solve(solver::CachingSolver, p, recurse=true)
    best_diff = Inf
    idx = 0
    num_ps = size(solver.ps_tree.ps, 2)
    for i in (num_ps-solver.new_count+1):num_ps
        diff = 0.
        for j in 1:size(solver.ps_tree.ps, 1)
            diff += abs2(solver.ps_tree.ps[j,i] - p[j])
        end
        if diff < best_diff
            best_diff = diff
            idx = i
        end
    end

    idx = indnearest(solver.ps_tree, p,
                     Alts([AltEntry(1, zeros(p), 0.0)], best_diff, idx))[1]

    set_extrapolation_origin(solver.basesolver,
                             solver.ps_tree.ps[:,idx], solver.zs[:,idx])

    if recurse
        solver.iters = 0
    end

    z = solve(solver.basesolver, p, 2500)
    solver.iters += needediterations(solver.basesolver)
    if needediterations(solver.basesolver) > 5 || ~hasconverged(solver.basesolver)
        if recurse && ~hasconverged(solver)
            a = 0.5
            best_a = 0.0
            while best_a < 1 && a > 0
                pa = (1-a) * solver.ps_tree.ps[:,idx] + a * p
                z = solve(solver, pa, false)
                if hasconverged(solver)
                    best_a = a
                    a = 1.0
                else
                    a = (a + best_a) / 2
                end
            end
        end
        if hasconverged(solver)
            solver.ps_tree.ps = [solver.ps_tree.ps p]
            solver.zs = [solver.zs z]
            solver.new_count += 1
        end
    end
    if solver.new_count > 0
        solver.new_count_limit -= 1
    end
    if solver.new_count > solver.new_count_limit
        solver.ps_tree = KDTree(solver.ps_tree.ps)
        solver.new_count = 0
        solver.new_count_limit = 2size(solver.ps_tree.ps, 2)
    end
    return z
end
