# Copyright 2016 Martin Holters
# See accompanying license file.

type CachingSolver{BaseSolver}
    basesolver::BaseSolver
    ps_tree::KDTree{Vector{Float64}, Matrix{Float64}}
    zs::Matrix{Float64}
    function CachingSolver(model::DiscreteModel)
        basesolver = BaseSolver(model)
        p = zeros(np(model))
        z = solve(basesolver, p, 2500)
        if ~hasconverged(basesolver)
            error("Falied to find initial solution.")
        end
        ps_tree = KDTree(zeros(np(model), 1))
        zs = reshape(z, nn(model), 1)
        return new(basesolver, ps_tree, zs)
    end
end

hasconverged(solver::CachingSolver) = hasconverged(solver.basesolver)

function solve(solver::CachingSolver, p, recurse=true)
    idx = indnearest(solver.ps_tree, p)[1]
    set_extrapolation_origin(solver.basesolver,
                             solver.ps_tree.ps[:,idx], solver.zs[:,idx])
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
            solver.ps_tree = KDTree([solver.ps_tree.ps p])
            solver.zs = [solver.zs z]
        end
    end
    return z
end
