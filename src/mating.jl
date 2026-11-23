# Mating policies (docs/src/roadmap.md, Phase 3): pairing the slots of the
# selected sires and dams by minimum coancestry, by index (assortative or
# disassortative), or at random under a limit per pair. The minimum-cost
# pairings solve a transportation problem exactly.

"""
    transport(supply, demand, cost; cap = typemax(Int)) -> Matrix{Int}

The integer flows `F` (`p × q`) minimising `sum(cost .* F)` subject to row
sums `supply`, column sums `demand` (equal totals) and `0 ≤ F ≤ cap` in
every cell: a transportation problem, solved exactly by successive
shortest paths with potentials (Dijkstra on reduced costs; `cost ≥ 0`).
The priority queue orders by distance then node, so ties go to the
lowest-numbered sire or dam and the result is deterministic. Throws an
`ArgumentError` when no flow meets `cap`.
"""
function transport(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer}, C::AbstractMatrix{<:Real};
    cap::Integer = typemax(Int))
    p, q = length(a), length(b)
    size(C) == (p, q) || throw(DimensionMismatch("cost must be $p × $q"))
    cap >= 1 || throw(ArgumentError("cap must be at least 1"))
    sum(a) == sum(b) || throw(ArgumentError("supplies and demands differ"))
    all(>=(0), C) || throw(ArgumentError("costs must be non-negative"))
    N = p + q + 2
    src, snk = p + q + 1, p + q + 2
    to, rcap, cost, nxt = Int[], Int[], Float64[], Int[]
    head = zeros(Int, N)
    function edge!(u, v, c, w)
        push!(to, v); push!(rcap, c); push!(cost, w); push!(nxt, head[u]); head[u] = length(to)
        push!(to, u); push!(rcap, 0); push!(cost, -w); push!(nxt, head[v]); head[v] = length(to)
    end
    for s = 1:p
        edge!(src, s, a[s], 0.0)
    end
    cell = zeros(Int, p, q)
    for s = 1:p, d = 1:q
        edge!(s, p + d, min(cap, a[s], b[d]), Float64(C[s, d]))
        cell[s, d] = length(to) - 1
    end
    for d = 1:q
        edge!(p + d, snk, b[d], 0.0)
    end
    π = zeros(N)
    dist = fill(Inf, N)
    pe = zeros(Int, N)
    flow, need = 0, sum(a)
    while flow < need
        fill!(dist, Inf); fill!(pe, 0)
        dist[src] = 0.0
        h = BinaryMinHeap{Tuple{Float64,Int}}()
        push!(h, (0.0, src))
        while !isempty(h)
            du, u = pop!(h)
            du > dist[u] && continue
            e = head[u]
            while e != 0
                v = to[e]
                if rcap[e] > 0
                    nd = du + cost[e] + π[u] - π[v]
                    if nd < dist[v] - 1e-12
                        dist[v] = nd; pe[v] = e; push!(h, (nd, v))
                    end
                end
                e = nxt[e]
            end
        end
        isfinite(dist[snk]) || throw(ArgumentError("no mating plan meets the limit per pair"))
        for v = 1:N
            isfinite(dist[v]) && (π[v] += dist[v])
        end
        δ, v = need - flow, snk
        while v != src
            e = pe[v]; δ = min(δ, rcap[e]); v = to[isodd(e) ? e + 1 : e - 1]
        end
        v = snk
        while v != src
            e = pe[v]; rcap[e] -= δ; r = isodd(e) ? e + 1 : e - 1; rcap[r] += δ; v = to[r]
        end
        flow += δ
    end
    return [rcap[cell[s, d]+1] for s = 1:p, d = 1:q]     # flow = capacity of the reverse edge
end
