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

function solve(solver::CachingSolver, p)
    idx = indnearest(solver.ps_tree, p)[1]
    set_extrapolation_origin(solver.basesolver,
                             solver.ps_tree.ps[:,idx], solver.zs[:,idx])
    z = solve(solver.basesolver, p, 5)
    if ~hasconverged(solver.basesolver)
        z = solve(solver.basesolver, p, 2500)
        if hasconverged(solver.basesolver)
            solver.ps_tree = KDTree([solver.ps_tree.ps p])
            solver.zs = [solver.zs z]
        end
    end
    return z
end
