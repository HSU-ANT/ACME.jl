# Copyright 2016 Martin Holters
# See accompanying license file.

type CachingSolver{BaseSolver}
    basesolver::BaseSolver
    ps_tree::KDTree{Vector{Float64}, Matrix{Float64}}
    zs::Matrix{Float64}
    function CachingSolver(model::DiscreteModel)
        basesolver = BaseSolver(model)
        ps_tree = KDTree(@compat Array{Float64}(np(model),0))
        return new(basesolver, ps_tree, @compat Array{Float64}(nn(model),0))
    end
end

hasconverged(solver::CachingSolver) = hasconverged(solver.basesolver)

function solve(solver::CachingSolver, p)
    if size(solver.zs)[2] > 0
        idx = indnearest(solver.ps_tree, p)[1]
        set_extrapolation_origin(solver.basesolver,
                                 solver.ps_tree.ps[:,idx], solver.zs[:,idx])
    end
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
