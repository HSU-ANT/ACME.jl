# Copyright 2016 Martin Holters
# See accompanying license file.

import Base.isless

type KDTree{Tcv<:AbstractVector,Tp<:AbstractMatrix}
    cut_dim::Vector{Int}
    cut_val::Tcv
    ps_idx::Vector{Int}
    ps::Tp
end

function KDTree(p::AbstractMatrix)
    function calc_cut_idx(min_idx, max_idx)
        N = max_idx - min_idx + 1
        N2 = 2^floor(Int, log2(N-1))
        if 3*N2/2 ≤ N
            return min_idx+N2-1
        else
            return min_idx+N-(N2/2)-1
        end
    end

    if size(p)[2] == 0
        return KDTree{Vector{eltype(p)},typeof(p)}([], [], [], p)
    end

    min_idx = zeros(Int, size(p,2)-1)
    max_idx = zeros(Int, size(p,2)-1)
    cut_idx = zeros(Int, size(p,2)-1)
    cut_dim = zeros(Int, size(p,2)-1)
    cut_val = zeros(eltype(p), size(p,2)-1)

    if size(p,2) == 1
        return KDTree{typeof(cut_val),typeof(p)}(cut_dim, cut_val, [1], p)
    end

    dim = indmax(var(p,2))
    p_idx = sortperm(vec(p[dim,:]))

    min_idx[1] = 1
    max_idx[1] = size(p)[2]
    cut_idx[1] = calc_cut_idx(min_idx[1], max_idx[1])
    cut_dim[1] = dim
    cut_val[1] = mean(p[dim, p_idx[cut_idx[1]:cut_idx[1]+1]])

    for n in 2:size(p,2)-1
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
        cut_idx[n] = calc_cut_idx(min_idx[n], max_idx[n])
        cut_dim[n] = dim
        cut_val[n] = mean(p[dim, p_idx[cut_idx[n]:cut_idx[n]+1]])
    end

    p_idx_final = zeros(Int, 1, size(p,2))
    for n in 1:size(p,2)
        parent_n = div(n+size(p,2)-1, 2);
        if mod(n+size(p,2),2) == 1
            p_idx_final[n] = p_idx[min_idx[parent_n]]
        else
            p_idx_final[n] = p_idx[max_idx[parent_n]]
        end
    end

    return KDTree{typeof(cut_val),typeof(p)}(cut_dim, cut_val, vec(p_idx_final), p)
end

type AltEntry{T}
    idx::Int
    delta::Vector{T}
    delta_norm::T
end

isless(e1::AltEntry, e2::AltEntry) = isless(e1.delta_norm, e2.delta_norm)

type Alts{T}
    entries::Vector{AltEntry{T}}
    best_dist::T
    best_pidx::Int
end

Alts{T}(p::Vector{T}) =
    Alts([AltEntry(1, zeros(T, length(p)), zero(T))], typemax(T), 0)

function siftup!(alts::Alts, i)
    entries = alts.entries
    parent = div(i, 2)
    while i > 1 && entries[i] < entries[parent]
        t = entries[i]
        entries[i] = entries[parent]
        entries[parent] = t
        i = parent
        parent = div(i, 2)
    end
end

function siftdown!(alts::Alts, i)
    entries = alts.entries
    N = length(entries)
    while true
        min = i
        if 2i <= N && entries[2i] < entries[min]
            min = 2i
        end
        if 2i+1 <= N && entries[2i+1] < entries[min]
            min = 2i+1
        end
        if min == i
            break
        end
        t = entries[i]
        entries[i] = entries[min]
        entries[min] = t
        i = min
    end
end

peek(alts::Alts) = minimum(alts.entries)

function dequeue!(alts::Alts)
    e = alts.entries[1]
    last_idx = length(alts.entries)
    alts.entries[1] = alts.entries[last_idx]
    deleteat!(alts.entries, last_idx)
    siftdown!(alts, 1)
    return e
end

function enqueue!(alts::Alts, entry::AltEntry)
    if entry.delta_norm < alts.best_dist
        push!(alts.entries, entry)
        siftup!(alts, length(alts.entries))
    end
end

function update_best_dist!(alts, dist, p_idx)
    if dist < alts.best_dist
        alts.best_dist = dist
        alts.best_pidx = p_idx
        i = length(alts.entries)
        while i > 0
            if alts.entries[i].delta_norm .≥ alts.best_dist
                last_idx = length(alts.entries)
                if last_idx ≠ i
                    alts.entries[i] = alts.entries[last_idx]
                end
                deleteat!(alts.entries, last_idx)
                if last_idx ≠ i
                    if i==1 || alts.entries[i] > alts.entries[div(i,2)]
                        siftdown!(alts, i)
                    else
                        siftup!(alts, i)
                    end
                end
            end
            i -= 1
        end
    end
end

indnearest(tree::KDTree, p::AbstractVector, alt = Alts(p)) =
    indnearest(tree, p, typemax(Int), alt)

function indnearest(tree::KDTree, p::AbstractVector, max_leaves::Int,
                    alt = Alts(p))
    l = 0
    p_idx = 0
    while l < max_leaves && ~isempty(alt.entries)
        best_alt = dequeue!(alt)
        idx = best_alt.idx
        delta = best_alt.delta
        delta_norm = best_alt.delta_norm

        while idx ≤ length(tree.cut_dim)
            dim = tree.cut_dim[idx]
            new_alt_delta_norm = delta_norm - delta[dim]^2 + (p[dim] - tree.cut_val[idx])^2
            if new_alt_delta_norm < alt.best_dist
                new_alt_delta = copy(delta)
                new_alt_delta[dim] = p[dim] - tree.cut_val[idx]
                new_alt_idx = p[dim] ≤ tree.cut_val[idx] ? 2idx+1 : 2idx;
                enqueue!(alt, AltEntry(new_alt_idx, new_alt_delta, new_alt_delta_norm))
            end
            if p[dim] ≤ tree.cut_val[idx]
                idx *= 2
            else
                idx = 2idx + 1
            end
        end
        idx -= length(tree.cut_dim)

        p_idx = tree.ps_idx[idx]
        dist = 0.
        for i in 1:length(p)
            dist += (p[i] - tree.ps[i, p_idx])^2
        end
        update_best_dist!(alt, dist, p_idx)

        l += 1
    end

    return alt.best_pidx, alt
end
