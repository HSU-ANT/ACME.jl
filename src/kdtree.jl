# Copyright 2016 Martin Holters
# See accompanying license file.

type KDTree{E<:Number,Tcv<:AbstractVector{E},Tp<:AbstractMatrix{E}}
    cut_dim::Vector{Int}
    cut_val::Tcv
    ps_idx::Vector{Int}
    ps::Tp
end

function KDTree(p::AbstractMatrix)
    depth = ceil(Int, log2(size(p)[2]))
    min_idx = zeros(Int, 2^depth-1)
    max_idx = zeros(Int, 2^depth-1)
    cut_idx = zeros(Int, 2^depth-1)
    cut_dim = zeros(Int, 2^depth-1)
    cut_val = zeros(eltype(p), 2^depth-1)

    if depth == 0
        return KDTree{eltype(p),typeof(cut_val),typeof(p)}(cut_dim, cut_val, [1], p)
    end

    dim = indmax(var(p,2))
    p_idx = sortperm(vec(p[dim,:]))

    min_idx[1] = 1;
    max_idx[1] = size(p)[2]
    cut_idx[1] = div(size(p)[2]+1, 2)
    cut_dim[1] = dim
    cut_val[1] = p[dim, p_idx[cut_idx[1]]]

    for n in 2:2^depth-1
        parent_n = div(n, 2)
        if mod(n, 2) == 0
            min_idx[n] = min_idx[parent_n]
            max_idx[n] = cut_idx[parent_n]
        else
            min_idx[n] = cut_idx[parent_n]+1
            max_idx[n] = max_idx[parent_n]
        end
        dim = indmax(var(p[:,p_idx[min_idx[n]:max_idx[n]]],2))
        idx = sortperm(vec(p[dim,p_idx[min_idx[n]:max_idx[n]]]))
        p_idx[min_idx[n]:max_idx[n]] = p_idx[idx + min_idx[n] - 1]
        cut_idx[n] = div(min_idx[n]+max_idx[n], 2)
        cut_dim[n] = dim
        cut_val[n] = p[dim, p_idx[cut_idx[n]]]
    end

    p_idx_final = zeros(Int, 1, 2^depth)
    for n in 1:2^depth
        parent_n = div(n+2^depth-1, 2);
        if mod(n,2) == 1
            p_idx_final[n] = p_idx[min_idx[parent_n]]
        else
            p_idx_final[n] = p_idx[max_idx[parent_n]]
        end
    end

    return KDTree{eltype(p),typeof(cut_val),typeof(p)}(cut_dim, cut_val, vec(p_idx_final), p)
end

type Alts{T}
    idx::Vector{Int}
    delta::Matrix{T}
    best_dist::T
end

Alts(T, k) = Alts([1], zeros(T, k, 1), Inf)

indnearest(tree::KDTree, p::AbstractVector, alt = Alts(eltype(p), length(p))) =
    indnearest(tree, p, typemax(Int), alt)

function indnearest(tree::KDTree, p::AbstractVector, max_leaves::Int,
                    alt = Alts(eltype(p), length(p)))
    depth = round(Int, log2(length(tree.cut_dim)+1))

    l = 0
    p_idx = 0
    while l < max_leaves && ~isempty(alt.idx)
        min_idx = indmin(sum(alt.delta.^2, 1))
        idx = alt.idx[min_idx]
        start_depth = floor(Int, log2(idx)) + 1
        new_alt_idx = zeros(Int, depth - start_depth + 1)
        new_alt_delta = repmat(alt.delta[:,min_idx], 1, depth - start_depth + 1)
        deleteat!(alt.idx, min_idx)
        alt.delta = [alt.delta[:,1:min_idx-1] alt.delta[:,min_idx+1:end]]
        for d in 1:depth - start_depth + 1
            dim = tree.cut_dim[idx]
            new_alt_delta[dim,d] = p[dim] - tree.cut_val[idx]
            if p[dim] ≤ tree.cut_val[idx]
                new_alt_idx[d] = 2idx + 1
                idx *= 2
            else
                new_alt_idx[d] = 2idx
                idx = 2idx + 1
            end
        end
        idx -= 2^depth - 1
        append!(alt.idx, new_alt_idx)
        alt.delta = [alt.delta new_alt_delta]
        p_idx = tree.ps_idx[idx]
        alt.best_dist = min(alt.best_dist, sum((p - tree.ps[:,p_idx]).^2))
        alt_keep = vec(sum(alt.delta.^2, 1) .< alt.best_dist)
        alt.idx = alt.idx[alt_keep]
        alt.delta = alt.delta[:,alt_keep]
        l += 1
    end

    return p_idx, alt
end
