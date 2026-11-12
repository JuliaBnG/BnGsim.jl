"""
    pack_loci(xy::AbstractMatrix{Bool}, loci) -> Matrix{UInt64}

Pack the haplotypes at `loci` (rows of `xy`) into 64-bit words, one column
per locus, so that allele and haplotype counts between two loci are popcounts.
"""
function pack_loci(xy::AbstractMatrix{Bool}, loci::AbstractVector{<:Integer})
    nh = size(xy, 2)
    nw = cld(nh, 64)
    W = zeros(UInt64, nw, length(loci))
    Threads.@threads for k in eachindex(loci)
        l = loci[k]
        @inbounds for h = 1:nh
            xy[l, h] && (W[(h-1)>>>6+1, k] |= UInt64(1) << ((h - 1) & 63))
        end
    end
    return W
end

"""
    r2(W, i, j, nh) -> Float64

Squared correlation of alleles between packed loci `i` and `j` over `nh`
phased haplotypes (Hill and Robertson 1968). `NaN` if either is monomorphic.
"""
@inline function r2(W::Matrix{UInt64}, i::Integer, j::Integer, nh::Integer)
    a = b = ab = 0
    @inbounds for w in axes(W, 1)
        x, y = W[w, i], W[w, j]
        a += count_ones(x)
        b += count_ones(y)
        ab += count_ones(x & y)
    end
    pa, pb = a / nh, b / nh
    D = ab / nh - pa * pb
    den = pa * (1 - pa) * pb * (1 - pb)
    return den > 0 ? D^2 / den : NaN
end

"""
    ld_decay(xy, lmp; minmaf = 0.05, edges = [5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6],
             npairs = 20_000, M = 1e8, rng = Random.default_rng()) -> DataFrame

Haplotype r² between loci with sample MAF ≥ `minmaf`, in distance bins
`[edges[k], edges[k+1])` bp. For each bin, `npairs` pairs on the same
chromosome are drawn by picking a random anchor locus and a random partner
in the bin downstream of it.

Columns: `lo`, `hi`, `d` (mean distance of the sampled pairs), `n`, `r2`
(mean), `r2adj = r2 - 1/nh` (sample-size correction for `nh` phased
haplotypes), and the Sved (1971) LD-based effective size
`ne = (1/r2adj - 1) / (4c)`, `c = d/M` Morgan, which refers to about
`t = 1/(2c)` generations ago (Hayes et al. 2003). `ne` is `NaN` when
`r2adj ≤ 0`.
"""
function ld_decay(
    xy::AbstractMatrix{Bool},
    lmp::DataFrame;
    minmaf::Real = 0.05,
    edges::AbstractVector{<:Real} = [5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6],
    npairs::Integer = 20_000,
    M::Real = 1e8,
    rng::AbstractRNG = Random.default_rng(),
)
    nh = size(xy, 2)
    frq = Float64.(lmp.frq)
    panel = findall(f -> min(f, 1 - f) ≥ minmaf, frq)
    W = pack_loci(xy, panel)
    pos = Int.(lmp.pos[panel])
    chr = lmp.chr[panel]
    np = length(panel)

    out = DataFrame(lo = Float64[], hi = Float64[], d = Float64[], n = Int[],
        r2 = Float64[], r2adj = Float64[], t = Float64[], ne = Float64[])
    for k = 1:(length(edges)-1)
        lo, hi = edges[k], edges[k+1]
        sr = sd = 0.0
        n = tries = 0
        while n < npairs && tries < 50npairs
            tries += 1
            i = rand(rng, 1:np)
            a = searchsortedfirst(pos, pos[i] + lo)
            b = searchsortedlast(pos, pos[i] + hi - 1)
            a > b && continue
            j = rand(rng, a:b)
            chr[j] == chr[i] || continue
            r = r2(W, i, j, nh)
            isnan(r) && continue
            sr += r
            sd += pos[j] - pos[i]
            n += 1
        end
        n == 0 && continue
        r̄, d̄ = sr / n, sd / n
        radj = r̄ - 1 / nh
        c = d̄ / M
        push!(out, (lo, hi, d̄, n, r̄, radj, 1 / (2c), radj > 0 ? (1 / radj - 1) / (4c) : NaN))
    end
    return out
end

"""
    vanraden_g(xy, loci) -> Matrix{Float64}

VanRaden (2008) method-1 genomic relationships among the individuals of `xy`
(columns `2i-1`, `2i`) at `loci`, centred on the allele frequencies of the
same individuals: `G = ZZ' / 2Σpq`.
"""
function vanraden_g(xy::AbstractMatrix{Bool}, loci::AbstractVector{<:Integer})
    nid = size(xy, 2) ÷ 2
    Z = Matrix{Float64}(undef, nid, length(loci))
    Threads.@threads for k in eachindex(loci)
        l = loci[k]
        @inbounds for i = 1:nid
            Z[i, k] = xy[l, 2i-1] + xy[l, 2i]
        end
    end
    p = vec(mean(Z; dims = 1)) ./ 2
    Z .-= 2p'
    return (Z * Z') ./ (2sum(p .* (1 .- p)))
end

"""
    founder_summary(xy, lmp; nchip = 20_000, chip_maf = 0.2, nqtl = 10_000, qtl_maf = 0.0,
                    presample = 3, rare = 0.05, ld = (;), rng = Random.default_rng())
        -> (scalars::NamedTuple, sfs::DataFrame, ld::DataFrame, chip, qtl)

Founder-level outcomes of plan §4.3 for one replicate of one arm. `xy` holds
the founders' phased haplotypes at the loci segregating among them, `lmp`
their map with sample frequencies (`founder_haplotypes`).

Chip and QTL loci are drawn as in the WCGALP pipeline (xyBnG `founder_pop`):
`presample × nchip` chip candidates with sample MAF > `chip_maf` are set
aside first, and the QTL are drawn from the loci with MAF > `qtl_maf` that
are left, so the QTL are depleted of common loci when the founders have
few of them. The chip is `nchip` of the candidates and `nqtl` QTL are kept.
(The pipeline makes its final cut after the founder expansion, among the
candidates that still meet the MAF condition; here it is made among the
founders.) `presample = 1` removes only the chip, which is what this
function did before 2026-11-11 and does not describe the pipeline. Scalars:

- `nid`, `nseg`: founders, segregating sites;
- `sum2pq`: Σ2pq over all segregating sites (expected heterozygosity × sites);
- `n_maf02`: sites with MAF ≥ 0.2, the pool for the chip;
- `chip_r2_adj`: mean r² between adjacent chip markers on a chromosome,
  and `chip_gap`, their mean spacing (bp);
- `qtl_maf_mean`, `qtl_rare`: mean QTL MAF and the fraction with MAF < `rare`;
- `qtl_rare_va`: share of Σ2pq over QTL carried by QTL with MAF < `rare`,
  i.e. the expected share of V_A from rare QTL under i.i.d. effects;
- `qtl_maf02`: fraction of QTL with MAF ≥ 0.2;
- `g_diag`, `g_offsd`: mean diagonal and SD of off-diagonals of the
  founder-centred VanRaden G at the chip;
- `me`: `1 / var(G off-diagonal)`, the effective number of independent
  chromosome segments implied by realized relationships (Goddard 2009).

`sfs` is the folded MAF spectrum in 0.05 bins; `ld` is `ld_decay` with
keyword arguments `ld`.
"""
function founder_summary(
    xy::AbstractMatrix{Bool},
    lmp::DataFrame;
    nchip::Integer = 20_000,
    chip_maf::Real = 0.2,
    nqtl::Integer = 10_000,
    qtl_maf::Real = 0.0,
    presample::Integer = 3,
    rare::Real = 0.05,
    ld::NamedTuple = (;),
    rng::AbstractRNG = Random.default_rng(),
)
    nlc, nh = size(xy)
    nlc == nrow(lmp) || throw(DimensionMismatch("xy has $nlc loci, lmp $(nrow(lmp))"))
    p = Float64.(lmp.frq)
    maf = min.(p, 1 .- p)
    h = 2 .* p .* (1 .- p)

    presample ≥ 1 || throw(ArgumentError("presample must be at least 1"))
    pool = findall(>(chip_maf), maf)
    length(pool) ≥ presample * nchip ||
        error("only $(length(pool)) loci with MAF > $chip_maf for $presample × $nchip chip candidates")
    cand = sample(rng, pool, presample * nchip; replace = false)   # in random order
    chip = sort(cand[1:nchip])
    pool = setdiff(findall(>(qtl_maf), maf), cand)
    length(pool) ≥ nqtl || error("only $(length(pool)) loci with MAF > $qtl_maf for $nqtl QTL")
    qtl = sort(sample(rng, pool, nqtl; replace = false))

    W = pack_loci(xy, chip)
    sr, nr, sg = 0.0, 0, 0.0
    for k = 1:(nchip-1)
        lmp.chr[chip[k]] == lmp.chr[chip[k+1]] || continue
        sr += r2(W, k, k + 1, nh)
        sg += Int(lmp.pos[chip[k+1]]) - Int(lmp.pos[chip[k]])
        nr += 1
    end

    qm, qh = maf[qtl], h[qtl]
    isr = qm .< rare
    G = vanraden_g(xy, chip)
    off = [G[i, j] for j in axes(G, 2) for i = 1:(j-1)]

    scalars = (
        nid = nh ÷ 2,
        nseg = nlc,
        sum2pq = sum(h),
        n_maf02 = count(≥(0.2), maf),
        chip_r2_adj = sr / nr,
        chip_gap = sg / nr,
        qtl_maf_mean = mean(qm),
        qtl_rare = mean(isr),
        qtl_rare_va = sum(qh[isr]) / sum(qh),
        qtl_maf02 = mean(qm .≥ 0.2),
        g_diag = mean(diag(G)),
        g_offsd = std(off),
        me = 1 / var(off),
    )
    edges = 0:0.05:0.5
    sfs = DataFrame(lo = edges[1:end-1], hi = edges[2:end],
        n = [count(x -> edges[k] < x ≤ edges[k+1], maf) for k = 1:(length(edges)-1)])
    return (scalars = scalars, sfs = sfs, ld = ld_decay(xy, lmp; rng = rng, ld...), chip = chip, qtl = qtl)
end
