# Copyright 2016, 2017, 2018, 2019 Martin Holters
# See accompanying license file.

import Base.deleteat!
import Base.isless
import Base.isempty
using Statistics: mean, var

mutable struct KDTree{Tcv<:AbstractVector,Tp<:AbstractMatrix}
    cut_dim::Vector{Int}
    cut_val::Tcv
    ps_idx::Vector{Int}
    ps::Tp
end

function KDTree(p::AbstractMatrix, Np=size(p,2))
    function calc_cut_idx(min_idx, max_idx)
        N = max_idx - min_idx + 1
        N2 = 2^floor(Int, log2(N-1))
        if 3(N2÷2) ≤ N
            return min_idx+N2-1
        else
            return min_idx + N - N2÷2 - 1
        end
    end

    if Np == 0
        return KDTree{Vector{eltype(p)},typeof(p)}([], [], [], p)
    end

    min_idx = zeros(Int, Np-1)
    max_idx = zeros(Int, Np-1)
    cut_idx = zeros(Int, Np-1)
    cut_dim = zeros(Int, Np-1)
    cut_val = zeros(eltype(p), Np-1)

    if Np == 1
        return KDTree{typeof(cut_val),typeof(p)}(cut_dim, cut_val, [1], p)
    end

    dim = argmax(vec(var(p[:,1:Np], dims=2)))
    p_idx = sortperm(vec(p[dim,:]))

    min_idx[1] = 1
    max_idx[1] = Np
    cut_idx[1] = calc_cut_idx(min_idx[1], max_idx[1])
    cut_dim[1] = dim
    cut_val[1] = mean(p[dim, p_idx[cut_idx[1]:cut_idx[1]+1]])

    for n in 2:Np-1
        parent_n = n ÷ 2
        if n % 2 == 0
            min_idx[n] = min_idx[parent_n]
            max_idx[n] = cut_idx[parent_n]
        else
            min_idx[n] = cut_idx[parent_n]+1
            max_idx[n] = max_idx[parent_n]
        end
        dim = argmax(vec(var(p[:,p_idx[min_idx[n]:max_idx[n]]], dims=2)))
        idx = sortperm(vec(p[dim,p_idx[min_idx[n]:max_idx[n]]]))
        p_idx[min_idx[n]:max_idx[n]] = p_idx[idx .+ min_idx[n] .- 1]
        cut_idx[n] = calc_cut_idx(min_idx[n], max_idx[n])
        cut_dim[n] = dim
        cut_val[n] = mean(p[dim, p_idx[cut_idx[n]:cut_idx[n]+1]])
    end

    p_idx_final = zeros(Int, Np)
    for n in 1:Np
        parent_n = (n+Np-1) ÷ 2;
        if (n+Np) % 2 == 1
            p_idx_final[n] = p_idx[min_idx[parent_n]]
        else
            p_idx_final[n] = p_idx[max_idx[parent_n]]
        end
    end

    return KDTree{typeof(cut_val),typeof(p)}(cut_dim, cut_val, p_idx_final, p)
end

mutable struct AltEntry{T}
    idx::Int
    delta::Vector{T}
    delta_norm::T
end

isless(e1::AltEntry, e2::AltEntry) = isless(e1.delta_norm, e2.delta_norm)

mutable struct Alts{T}
    entries::Vector{AltEntry{T}}
    best_dist::T
    best_pidx::Int
    number_valid::Int
end

Alts(p::Vector{T}) where {T} =
    Alts([AltEntry(1, zeros(T, length(p)), zero(T))], typemax(T), 0, 1)

function init!(alts::Alts{T}, best_dist, best_pidx) where {T}
    alts.number_valid = 1
    alts.entries[1].idx = 1
    fill!(alts.entries[1].delta, zero(T))
    alts.entries[1].delta_norm = zero(T)
    alts.best_dist = best_dist
    alts.best_pidx = best_pidx
end

function siftup!(alts::Alts, i)
    entries = alts.entries
    if i > alts.number_valid
        throw(BoundsError())
    end
    parent = i ÷ 2
    @inbounds while i > 1 && entries[i] < entries[parent]
        entries[i], entries[parent] = entries[parent], entries[i]
        i = parent
        parent = i ÷ 2
    end
end

function siftdown!(alts::Alts, i)
    if i < 1
        throw(BoundsError())
    end
    entries = alts.entries
    N = alts.number_valid
    @inbounds while true
        min = i
        if 2i ≤ N && entries[2i] < entries[min]
            min = 2i
        end
        if 2i+1 ≤ N && entries[2i+1] < entries[min]
            min = 2i+1
        end
        if min == i
            break
        end
        entries[i], entries[min] = entries[min], entries[i]
        i = min
    end
end

isempty(alts::Alts) = alts.number_valid == 0
peek(alts::Alts) = alts.entries[1]

function deleteat!(alts::Alts, i::Integer)
    alts.entries[i], alts.entries[alts.number_valid] = alts.entries[alts.number_valid], alts.entries[i]
    alts.number_valid -= 1
    if i ≤ alts.number_valid
        if i==1 || alts.entries[i] > alts.entries[i ÷ 2]
            siftdown!(alts, i)
        else
            siftup!(alts, i)
        end
    end
end

function dequeue!(alts::Alts)
    e = alts.entries[1]
    deleteat!(alts, 1)
    return e
end

function enqueue!(alts::Alts{T}, new_idx::Int, ref_delta::Vector{T},
                  delta_update_dim::Int, delta_update_val::T,
                  new_delta_norm::T) where {T}
    if alts.number_valid == length(alts.entries)
        delta = copy(ref_delta)
        delta[delta_update_dim] = delta_update_val
        push!(alts.entries, AltEntry(new_idx, delta, new_delta_norm))
    else
        alts.entries[alts.number_valid+1].idx = new_idx
        copyto!(alts.entries[alts.number_valid+1].delta, ref_delta)
        alts.entries[alts.number_valid+1].delta[delta_update_dim] = delta_update_val
        alts.entries[alts.number_valid+1].delta_norm = new_delta_norm
    end
    if alts.entries[alts.number_valid+1].delta_norm < alts.best_dist
        alts.number_valid += 1
        siftup!(alts, alts.number_valid)
    end
end

function update_best_dist!(alts, dist, p_idx)
    if dist < alts.best_dist
        alts.best_dist = dist
        alts.best_pidx = p_idx
        for i in alts.number_valid:-1:1
            if alts.entries[i].delta_norm ≥ alts.best_dist
                deleteat!(alts, i)
            end
        end
    end
end

indnearest(tree::KDTree, p::AbstractVector, alt = Alts(p)) =
    indnearest(tree, p, typemax(Int), alt)

function indnearest(tree::KDTree, p::AbstractVector, max_leaves::Int,
                    alt = Alts(p))
    if length(p) ≠ size(tree.ps, 1)
        throw(DimensionMismatch())
    end
    l = 0
    p_idx = 0
    while l < max_leaves && !isempty(alt)
        best_alt = dequeue!(alt)
        idx = best_alt.idx
        delta = best_alt.delta
        delta_norm = best_alt.delta_norm

        while idx ≤ length(tree.cut_dim)
            dim = tree.cut_dim[idx]
            new_alt_delta_norm = delta_norm - delta[dim]^2 + (p[dim] - tree.cut_val[idx])^2
            if new_alt_delta_norm < alt.best_dist
                new_alt_idx = p[dim] ≤ tree.cut_val[idx] ? 2idx+1 : 2idx;
                enqueue!(alt, new_alt_idx, delta, dim, p[dim] - tree.cut_val[idx], new_alt_delta_norm)
            end
            if p[dim] ≤ tree.cut_val[idx]
                idx *= 2
            else
                idx = 2idx + 1
            end
        end
        idx -= length(tree.cut_dim)

        p_idx = tree.ps_idx[idx]
        if p_idx < 1 || p_idx > size(tree.ps, 2)
            throw(BoundsError)
        end
        dist = 0.
        @simd for i in 1:length(p)
            @inbounds dist += (p[i] - tree.ps[i, p_idx])^2
        end
        update_best_dist!(alt, dist, p_idx)

        l += 1
    end

    return alt.best_pidx
end
