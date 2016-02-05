# Copyright 2016 Martin Holters
# See accompanying license file.

type CachingSolver{BaseSolver}
    basesolver::BaseSolver
    ps_tree::KDTree{Vector{Float64}, Matrix{Float64}}
    zs::Matrix{Float64}
    new_ps::Matrix{Float64}
    new_zs::Matrix{Float64}
    new_count_limit::Int
    function CachingSolver(model::DiscreteModel)
        basesolver = BaseSolver(model)
        p = zeros(np(model))
        z = solve(basesolver, p, 2500)
        if ~hasconverged(basesolver)
            error("Falied to find initial solution.")
        end
        ps_tree = KDTree(zeros(np(model), 1))
        zs = reshape(z, nn(model), 1)
        new_ps = zeros(np(model), 0)
        new_zs = zeros(nn(model), 0)
        return new(basesolver, ps_tree, zs, new_ps, new_zs, 50)
    end
end

hasconverged(solver::CachingSolver) = hasconverged(solver.basesolver)

function solve(solver::CachingSolver, p, recurse=true)
    idx = indnearest(solver.ps_tree, p)[1]

    best_new_diff = Inf
    best_new_idx = 0
    for i in 1:size(solver.new_ps, 2)
        diff = 0.
        for j in 1:size(solver.new_ps, 1)
            diff += (solver.new_ps[j,i] - p[j])^2
        end
        if diff < best_new_diff
            best_new_diff = diff
            best_new_idx = i
        end
    end

    if sumabs2(solver.ps_tree.ps[:,idx] - p) < best_new_diff
        set_extrapolation_origin(solver.basesolver,
                                 solver.ps_tree.ps[:,idx], solver.zs[:,idx])
    else
        set_extrapolation_origin(solver.basesolver,
                                 solver.new_ps[:,best_new_idx],
                                 solver.new_zs[:,best_new_idx])
    end

    z = solve(solver.basesolver, p, 2)
    if ~hasconverged(solver.basesolver)
        z = solve(solver.basesolver, p, 2500)
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
            solver.new_ps = [solver.new_ps p]
            solver.new_zs = [solver.new_zs z]
        end
        if ~isempty(solver.new_ps)
            solver.new_count_limit -= 1
        end
        if size(solver.new_ps, 2) > solver.new_count_limit
            solver.ps_tree = KDTree([solver.ps_tree.ps solver.new_ps])
            solver.zs = [solver.zs solver.new_zs]
            solver.new_ps = zeros(size(solver.new_ps, 1), 0)
            solver.new_zs = zeros(size(solver.new_zs, 1), 0)
            solver.new_count_limit = 50
        end
    end
    return z
end
