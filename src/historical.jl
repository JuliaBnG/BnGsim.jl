"""
    founder_haplotypes(res::FisherWrightResult, nid; rng = Random.default_rng())
        -> (xy::BitMatrix, lmp::DataFrame)

Sample `nid` individuals without replacement from the final generation of
`res` and return their phased haplotypes at the loci segregating among them.
Columns `2i-1` and `2i` are the two haplotypes of founder `i`. Alleles are
coded ancestral = 0, derived = 1 (no flipping), and `lmp.frq` is the
derived-allele frequency in the sample. `lmp.pos` is the genome-wide
coordinate; `lmp.chr` the chromosome.
"""
function founder_haplotypes(res::FisherWrightResult, nid::Integer; rng::AbstractRNG = Random.default_rng())
    np = length(res.active_haplotypes) ÷ 2
    0 < nid ≤ np || throw(ArgumentError("nid = $nid not in 1:$np"))
    ids = nid == np ? collect(1:np) : sort(sample(rng, 1:np, nid; replace = false))
    muts = Vector{Vector{UInt32}}(undef, 2nid)
    for (k, i) in enumerate(ids)
        muts[2k-1] = res.active_haplotypes[2i-1]
        muts[2k] = res.active_haplotypes[2i]
    end
    return dense_haplotypes(muts, res.chromosome_ends; rng=copy(rng))   # the caller's stream does not advance
end

"""
    founder_haplotypes(res::FisherWrightResult, nid, loci::TraitLoci; rng) -> (xy, lmp, rows, effects)

Sample `nid` founders from `res` as the method above (same draws from
`rng`), with the alleles of the historical trait of `loci` (FisherWright's
`selected_haplotypes`) merged into the neutral ones: a neutral mutation
at the coordinate of a trait locus is dropped, as the site is a trait
locus. Also returns the rows of `lmp` that are trait loci segregating
among the founders and the effect of their derived allele (the first
trait of `loci`; additive, `h = 0.5`).
"""
function founder_haplotypes(res::FisherWrightResult, nid::Integer, loci::TraitLoci;
    rng::AbstractRNG = Random.default_rng())
    all(==(0.5), loci.h) || throw(ArgumentError("the breeding loop handles additive trait loci only (h = 0.5)"))
    np = length(res.active_haplotypes) ÷ 2
    0 < nid ≤ np || throw(ArgumentError("nid = $nid not in 1:$np"))
    ids = nid == np ? collect(1:np) : sort(sample(rng, 1:np, nid; replace = false))
    tp = loci.pos
    muts = Vector{Vector{UInt32}}(undef, 2nid)
    for (k, i) in enumerate(ids), h in (2i - 1, 2i)
        neutral = filter(p -> !insorted(p, tp), res.active_haplotypes[h])
        muts[2k - (h == 2i - 1)] = sort!(vcat(neutral, res.selected_haplotypes[h]))
    end
    xy, lmp = dense_haplotypes(muts, res.chromosome_ends; rng = copy(rng))
    rows = findall(p -> insorted(p, tp), lmp.pos)
    effects = [loci.a[searchsortedfirst(tp, lmp.pos[r]), 1] for r in rows]
    return xy, lmp, rows, effects
end

"""
    founder_haplotypes(fxy::AbstractString, fmp::AbstractString, nid, cbp;
                       rng = Random.default_rng()) -> (xy::BitMatrix, lmp::DataFrame)

Sample `nid` founders from a historical population stored by xyBnG (the
WCGALP base, arm F1): `fxy` is an `xy` file holding a locus-major BitMatrix
(8-byte header `xy-`, then `nrow`, `ncol` as `Int64`, then the bits) and
`fmp` a serialized map `DataFrame` with genome-wide `pos`.

FisherWright v0.1.x wrote some loci out of order (a few blocks near the end
of the genome). Here rows are sorted by `pos` and `chr` is recomputed from the
cumulative chromosome ends `cbp`, so the output follows the same
conventions as the `FisherWrightResult` method: loci segregating among the
sampled founders, sorted, `frq` the sample frequency of allele 1.
"""
function founder_haplotypes(
    fxy::AbstractString,
    fmp::AbstractString,
    nid::Integer,
    cbp::AbstractVector{<:Integer};
    rng::AbstractRNG = Random.default_rng(),
)
    lmp = deserialize(fmp)
    xy = open(fxy) do io
        hdr = read(io, 8)
        Char.(hdr[1:3]) == ['x', 'y', '-'] && hdr[4] == 0 && hdr[5] == 0 && hdr[6] == 0 ||
            error("$fxy is not a locus-major BitArray xy file")
        nr, nc = read(io, Int), read(io, Int)
        nr == nrow(lmp) || throw(DimensionMismatch("$fxy has $nr loci, $fmp $(nrow(lmp))"))
        np = nc ÷ 2
        0 < nid ≤ np || throw(ArgumentError("nid = $nid not in 1:$np"))
        bits = Mmap.mmap(io, BitMatrix, (nr, nc))
        ids = sort(sample(rng, 1:np, nid; replace = false))
        bits[:, vec([2ids' .- 1; 2ids'])]
    end
    frq = vec(sum(xy; dims = 2)) ./ size(xy, 2)
    keep = findall(f -> 0 < f < 1, frq)
    ord = keep[sortperm(lmp.pos[keep])]
    pos = UInt32.(lmp.pos[ord])
    out = DataFrame(
        chr = Int8.(searchsortedfirst.(Ref(cbp), pos)),
        pos = pos,
        ref = lmp.ref[ord],
        alt = lmp.alt[ord],
        frq = Float32.(frq[ord]),
    )
    return xy[ord, :], out
end
