# Keyed meiosis on packed haplotypes (roadmap §4.4).
#
# Model, per gamete:
# - every chromosome assorts independently: the gamete starts the
#   chromosome on either parental haplotype with probability ½;
# - crossovers on chromosome c are a Poisson process with mean equal to its
#   genetic length L_c in Morgans (`BnGStructs.map_length`), placed
#   uniformly on the physical length (no interference; Haldane map, with a
#   recombination rate constant along each chromosome);
# - the gamete switches parental haplotype at every crossover.
# Loci are the rows of a `Haplotype`, at the genome-wide positions of a
# `VariantMap`; a crossover at base x falls after the last locus at or
# before x.

"""
    MeiosisMap(species, map::VariantMap)

The genome layout used by [`make_progeny!`](@ref), built once: for every
chromosome its physical start and end (genome-wide bp), its genetic length
in Morgans (from [`BnGStructs.map_length`](@ref)), and the range of loci of
`map` on it. The map is checked against the species with
`BnGStructs.check_coordinates`.
"""
struct MeiosisMap
    pos::Vector{UInt32}          # genome-wide position of every locus
    first::Vector{UInt32}        # first base of each chromosome (genome-wide)
    last::Vector{UInt32}         # last base of each chromosome
    morgan::Vector{Float64}      # genetic length of each chromosome
    emorgan::Vector{Float64}     # exp(-morgan), for the Poisson draws
    loci::Vector{UnitRange{Int}} # loci on each chromosome (may be empty)
end

function MeiosisMap(sp::Species, vm::VariantMap)
    check_coordinates(vm, sp)
    ends = cbp(sp)
    nchr = length(ends)
    first = [c == 1 ? UInt32(1) : ends[c-1] + UInt32(1) for c = 1:nchr]
    loci = Vector{UnitRange{Int}}(undef, nchr)
    i = 1
    n = length(vm)
    for c = 1:nchr
        j = i
        while j <= n && vm.chr[j] == c
            j += 1
        end
        loci[c] = i:(j-1)
        i = j
    end
    L = Float64.(map_length(sp))
    MeiosisMap(copy(vm.pos), first, copy(ends), L, exp.(-L), loci)
end

Base.length(mm::MeiosisMap) = length(mm.pos)

function Base.show(io::IO, mm::MeiosisMap)
    print(io, "MeiosisMap with ", length(mm.pos), " loci on ", length(mm.first),
        " chromosomes, ", round(sum(mm.morgan); digits = 3), " M")
end

# Poisson(λ) by inversion for the small means of one chromosome; its draws
# depend only on `rand(rng)`, so they do not change with Distributions.jl.
@inline function _poisson(rng::AbstractRNG, λ::Float64, eλ::Float64)
    λ > 30 && return rand(rng, Poisson(λ))
    k = 0
    p = eλ
    s = p
    u = rand(rng)
    while u > s
        k += 1
        p *= λ / k
        s += p
        p == 0 && break
    end
    return k
end

# Copy loci a:b of haplotype column `scol` of `src` into column `dcol` of
# `dst`; both `BitMatrix` with `nr` (padded) rows.
@inline function _copy_segment!(dst::BitMatrix, dcol::Int, src::BitMatrix, scol::Int, nr::Int, a::Int, b::Int)
    b < a && return nothing
    copyto!(dst, (dcol - 1) * nr + a, src, (scol - 1) * nr + a, b - a + 1)
    return nothing
end

# One gamete of the parent in columns 2p-1, 2p of `src`, written to column
# `dcol` of `dst`. `xo` is a work buffer. If `rec` is not `nothing`, the
# crossover positions are appended to it and the starting phases (`true` =
# the parent's second haplotype) are written to `phases`.
function _gamete!(dst::BitMatrix, dcol::Int, src::BitMatrix, p::Int, nr::Int,
    mm::MeiosisMap, rng::AbstractRNG, xo::Vector{UInt32}, rec, phases)
    pos = mm.pos
    for c in eachindex(mm.first)
        phase = rand(rng, Bool)
        k = _poisson(rng, mm.morgan[c], mm.emorgan[c])
        resize!(xo, k)
        for t = 1:k
            xo[t] = rand(rng, mm.first[c]:mm.last[c])
        end
        k > 1 && sort!(xo; alg = InsertionSort)
        if phases !== nothing
            phases[c] = phase
            append!(rec, xo)
        end
        r = mm.loci[c]
        isempty(r) && continue
        a = first(r)
        for x in xo
            j = searchsortedlast(view(pos, a:last(r)), x) + a - 1   # last locus ≤ x
            _copy_segment!(dst, dcol, src, 2p - 1 + phase, nr, a, j)
            a = max(a, j + 1)
            phase = !phase
        end
        _copy_segment!(dst, dcol, src, 2p - 1 + phase, nr, a, last(r))
    end
    return nothing
end

"""
    MeiosisRecord

Crossovers of the gametes made by one call of [`make_progeny!`](@ref).
Gamete `g` is column `g` of the offspring `Haplotype` (`2i-1` from the
sire, `2i` from the dam of mating `i`).

- `parent[g]`: the parent (column of the mating plan) that made it;
- `phase[c, g]`: `true` if chromosome `c` starts on the parent's second
  haplotype;
- `crossovers(rec, g)`: genome-wide positions of its crossovers, in
  increasing order within each chromosome.

With the parental records this determines the founder origin of every
segment, the input of segment-level IBD.
"""
struct MeiosisRecord
    parent::Vector{Int}
    phase::BitMatrix
    ptr::Vector{Int}
    pos::Vector{UInt32}
end

crossovers(r::MeiosisRecord, g::Integer) = view(r.pos, r.ptr[g]:(r.ptr[g+1]-1))
ncrossovers(r::MeiosisRecord, g::Integer) = r.ptr[g+1] - r.ptr[g]

"""
    make_progeny!(offspring::Haplotype, parents::Haplotype, plan, mm::MeiosisMap, key::SimKey;
                  record = false) -> offspring, or (offspring, MeiosisRecord)

Make one offspring per row of the mating plan `plan` (`n × 2`, sire and dam
as columns of `parents`): row `i` writes haplotypes `2i-1` (from the sire)
and `2i` (from the dam) of `offspring`, which must have at least `2n`
haplotypes and the loci of `parents`.

The draws of gamete `(i, k)` come from `keyed_rng(SimKey(key; domain =
Domain.MEIOSIS, entity = i, subcounter = k - 1))`, `k = 1` for the sire and
`2` for the dam. The result is a function of `key`, the plan and the
parents only: it does not depend on the number of threads or on the order
in which rows are processed. Several offspring of one pair are separate
rows. The caller sets `seed`, `replicate`, `generation` and `stage` of
`key`; its `domain`, `entity` and `subcounter` are overwritten.

See [`MeiosisMap`](@ref) for the meiosis model.
"""
function make_progeny!(off::Haplotype, par::Haplotype, plan::AbstractMatrix{<:Integer},
    mm::MeiosisMap, key::SimKey; record::Bool = false)
    n = size(plan, 1)
    size(plan, 2) == 2 || throw(DimensionMismatch("the mating plan needs two columns, sire and dam"))
    off.nhp >= 2n || throw(DimensionMismatch("offspring has $(off.nhp ÷ 2) individuals; the plan has $n rows"))
    off.nlc == par.nlc == length(mm) ||
        throw(DimensionMismatch("loci: offspring $(off.nlc), parents $(par.nlc), map $(length(mm))"))
    np = par.nhp ÷ 2
    all(p -> 1 <= p <= np, plan) || throw(ArgumentError("parents in the plan must lie in 1:$np"))
    dst, src = packed(off), packed(par)
    nr = size(src, 1)
    nchr = length(mm.first)
    recs = record ? [UInt32[] for _ = 1:2n] : nothing
    phase = record ? zeros(Bool, nchr, 2n) : nothing   # bytes: columns can be written from different tasks
    chunks = collect(Iterators.partition(1:n, max(1, cld(n, 4Threads.nthreads()))))
    @sync for rows in chunks
        Threads.@spawn begin
            xo = UInt32[]
            ph = record ? falses(nchr) : nothing
            for i in rows, k = 1:2
                g = 2i - 2 + k
                rng = keyed_rng(SimKey(key; domain = Domain.MEIOSIS, entity = i, subcounter = k - 1))
                _gamete!(dst, g, src, Int(plan[i, k]), nr, mm, rng, xo,
                    record ? recs[g] : nothing, ph)
                record && (phase[:, g] = ph)
            end
        end
    end
    record || return off
    ptr = cumsum([1; length.(recs)])
    rec = MeiosisRecord([Int(plan[(g+1)÷2, 2-g%2]) for g = 1:2n], BitMatrix(phase), ptr,
        reduce(vcat, recs; init = UInt32[]))
    return off, rec
end
