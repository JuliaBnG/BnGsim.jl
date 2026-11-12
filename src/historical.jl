"""
    struct HistoricalSpec

A historical (pre-founder) population run with `FisherWright.fisher_wright`.

# Fields
- `name::String`: arm label, e.g. `"F2"`.
- `ne::Int`: diploid population size, constant over the run.
- `nt::Int`: number of discrete, non-overlapping generations.
- `mr::Float64`: mutations per `1e8` bp per meiosis (μ = `mr × 1e-8`).
- `chromosome::Vector{UInt32}`: autosome lengths in bp.
- `M::Float64`: bp per Morgan (`1e8` is 1 cM/Mb).

Model, as implemented by FisherWright with its neutral defaults: diploid, dioecious Wright-Fisher
population with random union of gametes; infinite-sites mutation; Poisson
crossovers without interference; independent assortment of chromosomes;
the genome starts monomorphic, so `nt` determines how close the run is to
mutation-drift equilibrium (about `10ne` generations for the frequency
spectrum, much less for LD at short distances).
"""
struct HistoricalSpec
    name::String
    ne::Int
    nt::Int
    mr::Float64
    chromosome::Vector{UInt32}
    M::Float64
end

HistoricalSpec(name::AbstractString, sp::Species, nt::Integer, mr::Real) =
    HistoricalSpec(String(name), Int(sp.nid), Int(nt), Float64(mr), copy(sp.chromosome), Float64(sp.M))

"""
    historical_population(spec::HistoricalSpec; seed, verbose = false, progress_interval = 1000)
        -> FisherWrightResult

Run the historical population described by `spec`. The global RNG is seeded
with `seed`; FisherWright draws from task-local RNGs inside `Threads.@threads`,
so the result is reproducible only at the same thread count.
"""
function historical_population(
    spec::HistoricalSpec;
    seed::Integer,
    verbose::Bool = false,
    progress_interval::Integer = 1000,
)
    Random.seed!(seed)
    return fisher_wright(
        spec.ne,
        spec.nt,
        spec.chromosome,
        spec.mr;
        M = spec.M,
        result = true,
        verbose = verbose,
        progress_interval = progress_interval,
    )
end

"""
    continue_population(res::FisherWrightResult, ns, mr; M = 1e8, mut_base = 1e8,
                        rng = Random.default_rng()) -> FisherWrightResult

Continue a FisherWright population for `length(ns)` further generations with
census size `ns[g]` in generation `g`, e.g. a decline from the historical
size to a breed-like size before founders are sampled (arm F5). The model is
the same as in `fisher_wright` (dioecious WF, Poisson crossovers at `1/M` per
bp, infinite-sites mutation at `mr / mut_base` per bp per meiosis), except
that the population size may change every generation.

Runs on one thread with the explicit `rng`, so it is reproducible regardless
of the thread count. Mutations that fix during the continuation stay in the
haplotypes; `muts2bitarray` drops them as monomorphic.

Uses FisherWright internals (`random_mate!`, `cobp!`, `recombine`,
`merge_sorted`); `test/runtests.jl` guards that contract.
"""
function continue_population(
    res::FisherWrightResult,
    ns::AbstractVector{<:Integer},
    mr::Real;
    M::Real = 1e8,
    mut_base::Real = 1e8,
    rng::AbstractRNG = Random.default_rng(),
)
    all(>(1), ns) || throw(ArgumentError("every generation needs at least 2 individuals"))
    cbp = res.chromosome_ends
    chr = diff([0; Int.(cbp)])
    rmap = uniform_recombination_map(chr; M = M)
    p_mut = Poisson(Float64(cbp[end]) / mut_base * mr)
    span = UInt32(1):cbp[end]

    prt = res.active_haplotypes # read only; offspring are new vectors
    co = UInt32[]
    newm = UInt32[]
    for n in ns
        np = length(prt) ÷ 2
        pm = Matrix{Int}(undef, n, 2)
        FisherWright.random_mate!(pm, Vector{Bool}(undef, np), Int[], Int[]; rng = rng)
        off = [UInt32[] for _ = 1:2n]
        for i = 1:n
            s, d = pm[i, 1], pm[i, 2]
            FisherWright.recombine(prt[2s-1], prt[2s], off[2i-1],
                FisherWright.cobp!(co, rmap; rng = rng); rng = rng)
            FisherWright.recombine(prt[2d-1], prt[2d], off[2i],
                FisherWright.cobp!(co, rmap; rng = rng); rng = rng)
        end
        for h in eachindex(off)
            nm = rand(rng, p_mut)
            nm == 0 && continue
            resize!(newm, nm)
            for k = 1:nm
                newm[k] = rand(rng, span)
            end
            sort!(newm)
            off[h] = FisherWright.merge_sorted(off[h], newm)
        end
        prt = off
    end
    return FisherWrightResult(prt, copy(cbp), copy(res.substitutions))
end

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
    return muts2bitarray(muts, res.chromosome_ends)
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

"""
    read_packed_haplotypes(path, cbp; rng = Random.default_rng())
        -> (xy::BitMatrix, lmp::DataFrame)

Read phased haplotypes written by an external simulator (e.g. the msprime
driver of the founder-design paper) in the `BNGH1` format: 8-byte magic
`"BNGH1\\0\\0\\0"`, `Int64` `nloci` and `nhap`, `UInt32` genome-wide
positions, then each locus as `cld(nhap, 8)` bytes with haplotype `h` in bit
`(h-1) % 8` of byte `(h-1) ÷ 8` (numpy `packbits(bitorder = "little")`).
1 is the derived allele; haplotypes `2i-1`, `2i` form individual `i`.

Returns the same layout as the `FisherWrightResult` method: loci segregating
in the sample, sorted by position, `chr` from the cumulative chromosome ends
`cbp`, `frq` the sample frequency of allele 1, and random `ref`/`alt` bases.
"""
function read_packed_haplotypes(
    path::AbstractString,
    cbp::AbstractVector{<:Integer};
    rng::AbstractRNG = Random.default_rng(),
)
    xy, pos = open(path) do io
        read(io, 8) == UInt8['B', 'N', 'G', 'H', '1', 0, 0, 0] || error("$path is not a BNGH1 file")
        nlc, nhp = read(io, Int), read(io, Int)
        (nhp > 0 && iseven(nhp)) || error("nhap = $nhp must be positive and even")
        pos = Vector{UInt32}(undef, nlc)
        read!(io, pos)
        nb = cld(nhp, 8)
        packed = Matrix{UInt8}(undef, nb, nlc)
        read!(io, packed)
        xy = falses(nlc, nhp)
        for h = 1:nhp
            b, s = ((h - 1) >>> 3) + 1, (h - 1) & 7
            col = view(xy, :, h)
            @inbounds for l = 1:nlc
                col[l] = (packed[b, l] >>> s) & 0x01 == 0x01
            end
        end
        xy, pos
    end
    frq = vec(sum(xy; dims = 2)) ./ size(xy, 2)
    keep = findall(f -> 0 < f < 1, frq)
    ord = keep[sortperm(pos[keep])]
    allunique(pos[ord]) || error("$path has duplicate positions")
    bases = ('A', 'C', 'G', 'T')
    ref = [rand(rng, bases) for _ in ord]
    alt = [rand(rng, filter(!=(r), bases)) for r in ref]
    lmp = DataFrame(
        chr = Int8.(searchsortedfirst.(Ref(cbp), pos[ord])),
        pos = pos[ord],
        ref = ref,
        alt = alt,
        frq = Float32.(frq[ord]),
    )
    return xy[ord, :], lmp
end

"""
    write_xy(path, xy::BitMatrix)

Write a locus-major `BitMatrix` in xyBnG's `xy` format (8-byte header
`xy-` with zero type bytes, `Int64` `nrow` and `ncol`, then the bits), the
format read by xyBnG's pipeline and by `founder_haplotypes(fxy, fmp, …)`.
"""
function write_xy(path::AbstractString, xy::BitMatrix)
    open(path, "w") do io
        write(io, Int8['x', 'y', '-', 0, 0, 0, 0, 0])
        write(io, Int[size(xy)...])
        write(io, xy)
    end
    return path
end
