# Keyed meiosis on packed haplotypes (roadmap §4.4).
#
# Model, per gamete:
# - every chromosome assorts independently: the gamete starts the
#   chromosome on either parental haplotype with probability ½;
# - crossovers on chromosome c are a Poisson process with mean equal to its
#   genetic length L_c in Morgans (`BnGStructs.map_length`), placed
#   uniformly on the physical length (no interference; Haldane map, with a
#   recombination rate constant along each chromosome); or, with
#   interference ν ≠ 1, the gamma model (below); the length is the sire's
#   or the dam's with sex-specific maps;
# - the gamete switches parental haplotype at every crossover.
# Loci are the rows of a `Haplotype`, at the genome-wide positions of a
# `VariantMap`; a crossover at base x falls after the last locus at or
# before x.

"""
    MeiosisMap(species, map::VariantMap)

    MeiosisMap(species, map::VariantMap; male = nothing, female = nothing, interference = 1.0)

The genome layout used by [`make_progeny!`](@ref), built once: for every
chromosome its physical start and end (genome-wide bp), its genetic length
in Morgans (from `BnGStructs.map_length`), and the range of loci of
`map` on it. The map is checked against the species with
`BnGStructs.check_coordinates`.

Sex-specific maps: `male` and `female` are the genetic lengths of the
chromosomes in the meioses of sires and dams (default: those of the
species); a sire's gamete uses the male map, a dam's the female.

Crossover interference (`interference = ν`): the chiasmata of the
four-strand bundle form a stationary gamma renewal process along the
chromosome, two per Morgan, with gaps ``\\mathrm{Gamma}(\\nu, 1/(2\\nu))``
(the gamma model; McPeek and Speed 1995, Broman and Weber 2000), and each
chiasma is in the gamete with probability ½ (no chromatid interference).
The first chiasma lies at ``UX``, ``X`` the length-biased gap covering
the start of the chromosome, ``\\mathrm{Gamma}(\\nu + 1, 1/(2\\nu))``, and
``U`` uniform, so the process is stationary. A gamete then has on
average ``L`` crossovers on a chromosome of ``L`` Morgans. `ν = 1` is
Haldane's map (no interference), drawn as a Poisson count as before;
`ν > 1` is positive interference (human data fit `ν ≈ 4.3`). Positions are
linear in the physical length.
"""
struct MeiosisMap
    pos::Vector{UInt32}          # genome-wide position of every locus
    first::Vector{UInt32}        # first base of each chromosome (genome-wide)
    last::Vector{UInt32}         # last base of each chromosome
    morgan::Vector{Float64}      # genetic length of each chromosome (sires)
    emorgan::Vector{Float64}     # exp(-morgan), for the Poisson draws
    loci::Vector{UnitRange{Int}} # loci on each chromosome (may be empty)
    fmorgan::Vector{Float64}     # genetic length of each chromosome (dams)
    efmorgan::Vector{Float64}
    ν::Float64                   # interference (gamma shape; 1: none)
end

function MeiosisMap(sp::Species, vm::VariantMap; male = nothing, female = nothing, interference::Real = 1.0)
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
    Lm = male === nothing ? L : Float64.(collect(male))
    Lf = female === nothing ? L : Float64.(collect(female))
    length(Lm) == length(Lf) == nchr && all(x -> x >= 0 && isfinite(x), [Lm; Lf]) ||
        throw(ArgumentError("male and female need one finite, non-negative length per chromosome ($nchr)"))
    isfinite(interference) && interference > 0 || throw(ArgumentError("interference must be positive"))
    MeiosisMap(copy(vm.pos), first, copy(ends), Lm, exp.(-Lm), loci, Lf, exp.(-Lf), Float64(interference))
end

Base.length(mm::MeiosisMap) = length(mm.pos)

function Base.show(io::IO, mm::MeiosisMap)
    print(io, "MeiosisMap with ", length(mm.pos), " loci on ", length(mm.first),
        " chromosomes, ", round(sum(mm.morgan); digits = 3), " M",
        mm.fmorgan == mm.morgan ? "" : " (sires), $(round(sum(mm.fmorgan); digits = 3)) M (dams)",
        mm.ν == 1 ? "" : ", interference ν = $(mm.ν)")
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
    mm::MeiosisMap, rng::AbstractRNG, xo::Vector{UInt32}, rec, phases; dam::Bool = false)
    pos = mm.pos
    for c in eachindex(mm.first)
        phase = rand(rng, Bool)
        L = dam ? mm.fmorgan[c] : mm.morgan[c]
        if mm.ν == 1
            k = _poisson(rng, L, dam ? mm.efmorgan[c] : mm.emorgan[c])
            resize!(xo, k)
            for t = 1:k
                xo[t] = keyed_int(rng, mm.first[c], mm.last[c])
            end
            k > 1 && sort!(xo; alg = InsertionSort)
        else
            _gamma_crossovers!(xo, rng, L, mm.ν, mm.first[c], mm.last[c])
        end
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

# Crossovers of one gamete on a chromosome of `L` Morgans, bases a:b,
# under the gamma model with shape ν (see `MeiosisMap`), sorted: chiasmata
# at a stationary gamma renewal process of rate 2 per Morgan, each kept
# with probability ½, and their genetic positions mapped linearly to bases.
function _gamma_crossovers!(xo::Vector{UInt32}, rng::AbstractRNG, L::Float64, ν::Float64, a::UInt32, b::UInt32)
    empty!(xo)
    L > 0 || return xo
    θ = 1 / (2ν)
    t = rand(rng) * rand(rng, Gamma(ν + 1, θ))        # first chiasma: uniform in the length-biased gap
    span = Float64(b - a + 1)
    while t < L
        if rand(rng, Bool)
            push!(xo, a + UInt32(min(floor(Int, t / L * span), Int(b - a))))
        end
        t += rand(rng, Gamma(ν, θ))
    end
    return xo
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
                    record ? recs[g] : nothing, ph; dam = k == 2)
                record && (phase[:, g] = ph)
            end
        end
    end
    record || return off
    ptr = cumsum([1; length.(recs)])
    rec = MeiosisRecord([Int(plan[(g+1)÷2, 2-g%2]) for g = 1:2n], BitMatrix(phase), ptr,
        isempty(recs) ? UInt32[] : reduce(vcat, recs))   # no `init`: it defeats the fast vcat
    return off, rec
end
