# Copyright 2016 Martin Holters
# See accompanying license file.

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
        zs = reshape(z, nn(model), 1)
        return new(basesolver, ps_tree, zs, 0, 50, 0)
    end
end

hasconverged(solver::CachingSolver) = hasconverged(solver.basesolver)
needediterations(solver::CachingSolver) = solver.iters

function solve(solver::CachingSolver, p, recurse=true)
    idx = indnearest(solver.ps_tree, p)[1]

    best_new_diff = Inf
    best_new_idx = 0
    num_ps = size(solver.ps_tree.ps, 2)
    for i in (num_ps-solver.new_count+1):num_ps
        diff = 0.
        for j in 1:size(solver.ps_tree.ps, 1)
            diff += (solver.ps_tree.ps[j,i] - p[j])^2
        end
        if diff < best_new_diff
            best_new_diff = diff
            best_new_idx = i
        end
    end

    if best_new_diff < sumabs2(solver.ps_tree.ps[:,idx] - p)
        idx = best_new_idx
    end
    set_extrapolation_origin(solver.basesolver,
                             solver.ps_tree.ps[:,idx], solver.zs[:,idx])

    if recurse
        solver.iters = 0
    end
    z = solve(solver.basesolver, p, 5)
    solver.iters += needediterations(solver.basesolver)
    if ~hasconverged(solver.basesolver)
        z = solve(solver.basesolver, p, 2500)
        solver.iters += needediterations(solver.basesolver)
        if recurse && ~hasconverged(solver)
            a = 0.5
            best_a = 0
            while best_a < 1 && a > 0
                pa = (1-a) * solver.ps_tree.ps[:,idx] + a * p
                z = solve(solver, pa, false)
                if hasconverged(solver)
                    best_a = a
                    a = 1
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
        if solver.new_count > 0
            solver.new_count_limit -= 1
        end
        if solver.new_count > solver.new_count_limit
            solver.ps_tree = KDTree(solver.ps_tree.ps)
            solver.new_count = 0
            solver.new_count_limit = 50
        end
    end
    return z
end
