# The breeding loop (roadmap Phase 4, BnGsim 0.4–0.6): additive traits,
# evaluation through Breeding.jl, selection on an index and keyed mating,
# in discrete or overlapping generations.
#
# Model of one time step t (one birth cohort per step):
# - every animal has one phenotype per trait, y = μ + g + e, drawn at
#   birth, with residuals correlated across traits; a record is kept only
#   in the sex that expresses its trait (`trait.sex`) and enters
#   evaluations from age `trait.age` (in time steps) on;
# - all records available at t and the full pedigree enter the evaluation
#   (single- or multi-trait animal model), with the true or, optionally,
#   REML-re-estimated covariance components;
# - candidates are ranked on the index I = w'â of their EBVs; sires and
#   dams are chosen among the males and females of age at most
#   `max_sire_age` and `max_dam_age` (age 0: born at t) and given numbers
#   of offspring; matings pair the sire and dam slots at random;
# - the offspring form cohort t + 1; between t and t + 1 animals survive
#   with an age-specific probability and, optionally, unselected adults
#   are culled. Nothing is removed from the population: dead or culled
#   animals keep their records, pedigree and haplotypes, but are no longer
#   candidates and make no later records.
# Every random draw is keyed by a `SimKey` (domains MATING, PHENOTYPE,
# SELECTION, SURVIVAL), so a run is a function of its key and inputs.

# Stage codes within domains; fixed values, never schedule positions.
const _STAGE_SEX = 0x0001
const _STAGE_PAIRING = 0x0002
const _STAGE_RANDOM_SELECTION = 0x0003
const _STAGE_HERD = 0x0004              # Domain.MATING: a founder's herd
const _STAGE_CG = 0x0005                # Domain.PHENOTYPE: herd × step effects
const _STAGE_PAIR_COST = 0x0006         # Domain.MATING: random costs under a pair limit
const _STAGE_PE = 0x0007                # Domain.PHENOTYPE: permanent environment, at birth
const _STAGE_RECORD = 0x0008            # Domain.PHENOTYPE: a repeated record, per step
const _STAGE_MISSING = 0x0009           # Domain.PHENOTYPE: a record missing at random
const _STAGE_HERD_MISSING = 0x000a      # Domain.PHENOTYPE: herds not recording a trait

# The repeated records of a population (traits with `records > 1`), one
# entry per record: animal, trait, time step, value.
struct RepeatedRecords
    id::Vector{Int32}
    trait::Vector{Int32}
    step::Vector{Int32}
    y::Vector{Float64}
end
RepeatedRecords() = RepeatedRecords(Int32[], Int32[], Int32[], Float64[])
Base.:(==)(a::RepeatedRecords, b::RepeatedRecords) = a.id == b.id && a.trait == b.trait && a.step == b.step && a.y == b.y
Base.isequal(a::RepeatedRecords, b::RepeatedRecords) =
    isequal(a.id, b.id) && isequal(a.trait, b.trait) && isequal(a.step, b.step) && isequal(a.y, b.y)
Base.hash(a::RepeatedRecords, h::UInt) = hash((a.id, a.trait, a.step, a.y), hash(:RepeatedRecords, h))

"""
    Population

Animals of a breeding program, one row per animal in order of birth.
Animal ids are `1:n` and equal the row; nothing is ever removed, so the
equality holds (invariant I4 holds trivially). Tables that leave the
population carry the id explicitly ([`pedigree`](@ref)), and relationship
matrices on a subset are indexed through a [`PedigreeView`](@ref).

Fields:

- `sire`, `dam` (`0` = unknown), `male`, `generation` (the birth cohort,
  i.e. the time step of birth; founders are 0; with discrete generations
  it is the generation number), `exit` (the first time step at which the
  animal is dead or culled; `typemax(Int32)` while it lives);
- `tbv`, `phenotype`, `ebv`: `n × k` matrices, one column per trait, of
  genetic values, phenotypes (`NaN` in the sex that does not express the
  trait) and the last evaluation;
- `blocks`, `first`, `slot`: the haplotypes, one `Haplotype` per cohort
  (`nothing` once no genome of the cohort is kept), with `first[t + 1]`
  the id of the first animal of cohort `t` and `slot[i]` the position of
  animal `i` in its cohort's block (`0` once its genome is dropped, see
  `retain` in [`Scheme`](@ref)); read genomes with [`genomes`](@ref);
- `qtl` (one `TraitQTL` per trait), `traits`, `Lr` (lower-triangular
  factor of the true residual covariance matrix, ``R = L_rL_r'``), `key`;
- `G0`, `R0`: the genetic and residual covariance matrices of the next
  evaluation (the true ones of the founders unless re-estimated);
- `meioses`: the crossover records of each cohort (`meioses[t]` for
  cohort `t ≥ 1`, parents as population ids) when the scheme records
  them, `nothing` otherwise; see [`founder_segments`](@ref).

Build one with [`founder_population`](@ref), extend it with
[`step!`](@ref), and save or restore it with [`save_checkpoint`](@ref)
and [`load_checkpoint`](@ref).
"""
mutable struct Population
    sire::Vector{Int32}
    dam::Vector{Int32}
    male::BitVector
    generation::Vector{Int32}
    exit::Vector{Int32}
    tbv::Matrix{Float64}
    phenotype::Matrix{Float64}
    ebv::Matrix{Float64}
    blocks::Vector{Union{Nothing,Haplotype}}
    first::Vector{Int}
    slot::Vector{Int32}
    qtl::Vector{TraitQTL}
    traits::Vector{Trait}
    Lr::Matrix{Float64}
    key::SimKey
    G0::Matrix{Float64}
    R0::Matrix{Float64}
    meioses::Vector{Union{Nothing,MeiosisRecord}}
    herd::Vector{Int32}
    σcg::Vector{Float64}
    nrec::Vector{Int}
    σpe::Vector{Float64}
    pe::Matrix{Float64}
    PE0::Matrix{Float64}
    records::RepeatedRecords
    Lpe::Matrix{Float64}
    pmiss::Vector{Float64}
    hmiss::Vector{Float64}
    nobs::Matrix{Int32}
    parity::Vector{Int32}
end

Base.length(pop::Population) = length(pop.sire)
ntraits(pop::Population) = length(pop.traits)
ngenerations(pop::Population) = length(pop.blocks)
generation_ids(pop::Population, g::Integer) = pop.first[g+1]:(g + 2 <= length(pop.first) ? pop.first[g+2] - 1 : length(pop))

function Base.show(io::IO, pop::Population)
    print(io, "Population of ", length(pop), " animals in ", ngenerations(pop), " cohorts (",
        join(("\"$(t.name)\" (h² = $(t.h²))" for t in pop.traits), ", "), ")")
end

"""
    founder_population(hap::Haplotype, qtl::TraitQTL, trait::Trait, key::SimKey) -> Population
    founder_population(hap::Haplotype, qtl::MultiTraitQTL, traits, key::SimKey; Re = I) -> Population

Cohort 0 of a breeding program: the individuals of `hap`, unrelated
founders, with sexes drawn at random (`Domain.MATING`), genetic values
from `qtl` and one phenotype per trait (`Domain.PHENOTYPE`). Each trait's
residual standard deviation is fixed from the trait as
``\\sigma_e = \\sigma_a\\sqrt{(1 - h^2)/h^2}``, with ``\\sigma_a`` the trait's
`σₐ`, so ``h^2`` is the heritability in the founders; `Re` is the
residual correlation matrix of several traits. `key` carries the seed
and replicate; its other fields are set per draw.

The covariance components of the evaluations are those of the founders:
``G_0 = D\\,C\\,D``, ``D`` the diagonal of the traits' `σₐ` and ``C`` the
correlation matrix of the founders' genetic values (exactly the `R` of a
multi-trait [`eqtl`](@ref)), and ``R_0 = D_e R_e D_e``.

Contemporary groups: with `herds = H`, each founder is placed in one of
`H` herds at random (a keyed draw, `Domain.MATING`), and every offspring
is born in its dam's herd. A record of trait `j` belongs to the group of
the animal's herd and the time step of the record (birth plus the trait's
`age`); each group and trait has an effect ``N(0, \\sigma^2_{cg,j})``, a
keyed draw per herd and step (`Domain.PHENOTYPE`), added to the
phenotype, with ``\\sigma^2_{cg,j} = c^2_j V_{p,j}``, ``V_p = \\sigma^2_a/h^2``
(the variance within groups), from `c2` (a number for all traits, or one
per trait; default `0`, no effect). `h^2` stays the heritability within
groups. Effects of different traits are
independent. When some `c2 > 0`, evaluations fit the groups as fixed
effects in place of the overall mean (see [`evaluate!`](@ref)).

Repeated records: with `records = m > 1` for a trait (a number, or one
per trait), an animal gets a record of it at every time step it is alive
from the trait's `age` on, at most `m`, ``y = \\mu + a + pe + cg + e``: a
permanent-environment effect ``pe`` drawn at birth (keyed,
`Domain.PHENOTYPE`) and a temporary residual ``e`` drawn per record,
with ``\\sigma^2_{pe} = (r - h^2)V_p`` and ``\\sigma^2_e = (1 - r)V_p`` for
the `repeatability` ``r`` (default ``h^2``: no permanent environment; one
number, or one per trait with `nothing` for ``h^2``). With several
traits, the permanent environments are correlated by `Rpe` (among the
traits with ``r > h^2``) and the temporary residuals of the records one
animal gets at one step by `Re`, drawn jointly from the animal's stream
of that step; single-record traits are then recorded step by step too.
The records are made at the start of each [`step!`](@ref)
(`pop.records`); `pop.phenotype` holds the mean of an animal's records
so far (`NaN` before the first). `pop.R0` is then the covariance of the
temporary residuals and `pop.PE0` that of the permanent environments
used by the evaluations.

Missing records: `missing` (a number, or one per trait) is the
probability that a record is missing at random, a keyed draw per animal,
trait and time step of the record; `unrecorded_herds` the probability
that a herd never records a trait, a keyed draw per herd and trait. A
missing record is not made (`NaN` as phenotype, no entry in
`pop.records`); the values of the other records do not change, as the
masks have their own keyed streams.

A trait's `sex` (`0` females, `1` males, `2` both) limits its records to
one sex; its `age`, a whole number of time steps, is the age at which a
record becomes available to evaluations (`0`: at birth). Phenotypes are
drawn for every animal and trait either way, from the animal's own keyed
stream, so these settings change which records are seen, not the values
drawn.
"""
function founder_population(hap::Haplotype, qtl::TraitQTL, trait::Trait, key::SimKey; herds::Integer = 1,
    c2 = 0.0, records::Integer = 1, repeatability = nothing, missing = 0.0, unrecorded_herds = 0.0)
    _founders(hap, [qtl], Trait[trait], key, ones(1, 1), herds, c2, records, repeatability, ones(1, 1), missing,
        unrecorded_herds)
end

function founder_population(hap::Haplotype, qtl::MultiTraitQTL, traits::AbstractVector{<:Trait}, key::SimKey;
    Re = I, herds::Integer = 1, c2 = 0.0, records = 1, repeatability = nothing, Rpe = I, missing = 0.0,
    unrecorded_herds = 0.0)
    k = length(traits)
    qtl.trait_names == [t.name for t in traits] ||
        throw(ArgumentError("the traits must be those of the QTL, in order: $(qtl.trait_names)"))
    isempty(qtl.dominance) || all(iszero, qtl.dominance) ||
        throw(ArgumentError("the breeding loop handles additive traits only"))
    Rm = Re isa UniformScaling ? Matrix{Float64}(Re(k)) : Matrix{Float64}(Re)
    size(Rm) == (k, k) && issymmetric(Rm) && all(Rm[j, j] == 1 for j = 1:k) && isposdef(Rm) ||
        throw(ArgumentError("Re must be a $k × $k positive definite correlation matrix"))
    qs = [TraitQTL(qtl.trait_names[j], qtl.loci, qtl.additive[:, j]; intercept = qtl.intercept[j]) for j = 1:k]
    Pm = Rpe isa UniformScaling ? Matrix{Float64}(Rpe(k)) : Matrix{Float64}(Rpe)
    size(Pm) == (k, k) && issymmetric(Pm) && all(Pm[j, j] == 1 for j = 1:k) ||
        throw(ArgumentError("Rpe must be a $k × $k correlation matrix"))
    _founders(hap, qs, Vector{Trait}(traits), key, Rm, herds, c2, records, repeatability, Pm, missing,
        unrecorded_herds)
end

function _founders(hap, qs, traits, key, Re, herds, c2, records, repeatability, Rpe, missing, unrecorded_herds)
    k = length(traits)
    per(x, what) = (v = x isa Real ? fill(Float64(x), k) : Float64.(collect(x));
        length(v) == k && all(p -> 0 <= p < 1, v) || throw(ArgumentError("$what must lie in [0, 1), one number or one per trait"));
        v)
    pmiss, hmiss = per(missing, "missing"), per(unrecorded_herds, "unrecorded_herds")
    h²s = [t.h² for t in traits]
    nrec = records isa Integer ? fill(Int(records), k) : Int.(collect(records))
    length(nrec) == k && all(>=(1), nrec) || throw(ArgumentError("records must be at least 1, one number or one per trait"))
    rv = repeatability === nothing ? copy(h²s) : repeatability isa Real ? fill(Float64(repeatability), k) :
         [x === nothing || isnan(x) ? h²s[j] : Float64(x) for (j, x) in enumerate(repeatability)]
    length(rv) == k || throw(ArgumentError("repeatability must be one number or one per trait"))
    for j = 1:k
        nrec[j] == 1 && rv[j] != h²s[j] && throw(ArgumentError("a repeatability other than h² needs records > 1"))
        h²s[j] <= rv[j] < 1 || throw(ArgumentError("the repeatability must lie in [h², 1)"))
    end
    herds >= 1 || throw(ArgumentError("herds must be at least 1"))
    c2v = c2 isa Real ? fill(Float64(c2), k) : Float64.(collect(c2))
    length(c2v) == k && all(x -> 0 <= x < 1, c2v) ||
        throw(ArgumentError("c2 must be a number in [0, 1) or one per trait"))
    any(_has_dominance, qs) && throw(ArgumentError("the breeding loop handles additive traits only"))
    all(t -> isinteger(t.age), traits) ||
        throw(ArgumentError("a trait's age must be a whole number of time steps"))
    allunique(t.name for t in traits) || throw(ArgumentError("trait names must be unique"))
    n = hap.nhp ÷ 2
    n >= 2 || throw(ArgumentError("need at least two founders"))
    G = reduce(hcat, [tbv(hap, q) for q in qs])
    σa = [t.σₐ for t in traits]
    Vp = [t.σₐ^2 / t.h² for t in traits]
    # with repeated records, σe is the temporary residual: Vp(1 - r)
    σe = [nrec[j] > 1 ? sqrt(Vp[j] * (1 - rv[j])) : t.σₐ * sqrt((1 - t.h²) / t.h²) for (j, t) in enumerate(traits)]
    σpe = [nrec[j] > 1 ? sqrt(Vp[j] * (rv[j] - h²s[j])) : 0.0 for j = 1:k]
    pej = findall(>(0), σpe)
    Lpe = zeros(k, k)
    if !isempty(pej)
        F = cholesky(Symmetric(Rpe[pej, pej]); check = false)
        issuccess(F) || throw(ArgumentError("Rpe is not positive definite among the traits with a permanent environment"))
        Lpe[pej, pej] = Diagonal(σpe[pej]) * Matrix(F.L)
    end
    C = length(traits) == 1 ? ones(1, 1) : cor(G)
    Lr = Diagonal(σe) * Matrix(cholesky(Symmetric(Re)).L)
    G0 = σa .* C .* σa'
    R0 = σe .* Re .* σe'
    pop = Population(zeros(Int32, n), zeros(Int32, n), falses(n), zeros(Int32, n), fill(typemax(Int32), n), G,
        fill(NaN, size(G)), fill(NaN, size(G)), Union{Nothing,Haplotype}[hap], [1], Int32.(1:n), qs, traits, Lr, key, G0, R0,
        Union{Nothing,MeiosisRecord}[], ones(Int32, n), sqrt.(c2v .* Vp),
        nrec, σpe, zeros(n, k), Lpe * Lpe', RepeatedRecords(), Lpe, pmiss, hmiss, zeros(Int32, n, k), zeros(Int32, n))
    if herds > 1
        for i = 1:n
            pop.herd[i] = keyed_int(keyed_rng(SimKey(key; generation = 0, domain = Domain.MATING,
                stage = _STAGE_HERD, entity = i)), Int32(1), Int32(herds))
        end
    end
    _sexes!(pop, 1:n, 0)
    _phenotypes!(pop, 1:n, 0)
    return pop
end

_has_dominance(q::TraitQTL) = any(!iszero, q.dominance)

# the current time step: the cohort of the youngest animals
current_step(pop::Population) = ngenerations(pop) - 1

"""
    pedigree(pop::Population[, ids]) -> DataFrame

The animals `ids` of `pop` (default: all), one row each, with the id
explicit, so that a subset, a reordering or a concatenation of tables
keeps the identity of every animal:

- `id`, `sire`, `dam` (`0` = unknown; parents are ids, and may be absent
  from a subset), `male`, `generation` (the birth cohort), `exit` (the
  first time step at which the animal is dead or culled; `missing` while
  it lives), `herd`;
- per trait `<name>`: `tbv_<name>`, the true genetic value;
  `phenotype_<name>`, the phenotype (`missing` in the sex that does not
  express the trait; it is drawn at birth, and is available to an
  evaluation only from the trait's `age` on, see [`observed`](@ref);
  with repeated records the mean of the animal's records so far, which
  are in `pop.records`);
  and `ebv_<name>`, the last evaluation (`missing` before one, and for
  animals outside an evaluation window, see [`Scheme`](@ref)).

Inside `pop`, ids equal rows; use the `id` column, not the row number,
once a table leaves it. [`PedigreeView`](@ref) builds the dense rows of
relationship matrices from these ids.
"""
function pedigree(pop::Population, ids = 1:length(pop))
    ids = collect(Int, ids)
    all(i -> 1 <= i <= length(pop), ids) || throw(ArgumentError("ids must lie in 1:$(length(pop))"))
    nan(x) = isnan(x) ? missing : x
    df = DataFrame(
        id = ids,
        sire = Int.(pop.sire[ids]),
        dam = Int.(pop.dam[ids]),
        male = pop.male[ids],
        generation = Int.(pop.generation[ids]),
        exit = [x == typemax(Int32) ? missing : Int(x) for x in pop.exit[ids]],
        herd = Int.(pop.herd[ids]),
    )
    for (j, t) in enumerate(pop.traits)
        df[!, "tbv_" * t.name] = pop.tbv[ids, j]
        df[!, "phenotype_" * t.name] = nan.(pop.phenotype[ids, j])
        df[!, "ebv_" * t.name] = nan.(pop.ebv[ids, j])
    end
    return df
end

"""
    alive(pop::Population, t = current_step(pop)) -> Vector{Int}

Ids of the animals born at or before time step `t` that are neither dead
nor culled at `t`.
"""
alive(pop::Population, t::Integer = current_step(pop)) =
    filter(i -> pop.generation[i] <= t < pop.exit[i], 1:length(pop))

_expresses(t::Trait, male::Bool) = t.sex == 2 || (t.sex == 1) == male

"""
    recorded(pop::Population, t = current_step(pop); trait = 1) -> Vector{Int}

Ids of the animals whose phenotype of trait number `trait` is available
to an evaluation at time step `t`, i.e. expressed in their sex and of age
at least the trait's `age`. See also [`observed`](@ref).
"""
recorded(pop::Population, t::Integer = current_step(pop); trait::Integer = 1) =
    findall(view(observed(pop, t), :, trait))

"""
    observed(pop::Population, t = current_step(pop)) -> BitMatrix

The `n × k` pattern of the records available to an evaluation at time
step `t`: entry `(i, j)` is `true` if animal `i` expresses trait `j` (by
sex), is at least the trait's `age` old, and was alive at that age (an
animal that died or was culled earlier has no record).
"""
function observed(pop::Population, t::Integer = current_step(pop))
    obs = falses(length(pop), ntraits(pop))
    for j = 1:ntraits(pop), i = 1:length(pop)
        a = Int(pop.traits[j].age)
        obs[i, j] = !isnan(pop.phenotype[i, j]) && t - pop.generation[i] >= a &&
                    pop.generation[i] + a < pop.exit[i]
    end
    return obs
end

function _sexes!(pop::Population, ids, g)
    for i in ids
        pop.male[i] = rand(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.MATING,
            stage = _STAGE_SEX, entity = i)), Bool)
    end
    any(pop.male[ids]) && !all(pop.male[ids]) ||
        throw(ArgumentError("generation $g has one sex only; increase its size"))
end

# Residuals e = L_r z, z standard normal, drawn one by one from the
# animal's stream (trait 1 first): the first trait's residual is the same
# whatever the number of traits.
function _phenotypes!(pop::Population, ids, g)
    k = ntraits(pop)
    z = zeros(k)
    for i in ids
        rng = keyed_rng(SimKey(pop.key; generation = g, domain = Domain.PHENOTYPE, entity = i))
        for j = 1:k
            z[j] = randn(rng)
        end
        if _record_mode(pop)                   # records come with the steps
            if any(>(0), pop.σpe)
                prng = keyed_rng(SimKey(pop.key; generation = g, domain = Domain.PHENOTYPE, stage = _STAGE_PE,
                    entity = i))
                zp = [randn(prng) for _ = 1:k]
                pop.pe[i, :] = pop.Lpe * zp
            end
            pop.phenotype[i, :] .= NaN
            continue
        end
        for j = 1:k
            e = 0.0
            for l = 1:j
                e += pop.Lr[j, l] * z[l]
            end
            t = pop.traits[j]
            pop.phenotype[i, j] = _expresses(t, pop.male[i]) && !_missing(pop, i, j, g + Int(t.age)) ?
                t.μ + pop.tbv[i, j] + e + _cg(pop, i, j) : NaN
        end
    end
end

# Whether animal i's record of trait j made at step s is missing: its herd
# does not record the trait (one keyed draw per herd and trait), or the
# record is missing at random (one keyed draw per animal, step and trait).
function _missing(pop::Population, i, j, s)
    if pop.hmiss[j] > 0
        rng = keyed_rng(SimKey(pop.key; generation = 0, domain = Domain.PHENOTYPE, stage = _STAGE_HERD_MISSING,
            entity = pop.herd[i]))
        u = 0.0
        for _ = 1:j
            u = rand(rng)
        end
        u < pop.hmiss[j] && return true
    end
    pop.pmiss[j] > 0 || return false
    rng = keyed_rng(SimKey(pop.key; generation = s, domain = Domain.PHENOTYPE, stage = _STAGE_MISSING, entity = i))
    u = 0.0
    for _ = 1:j
        u = rand(rng)
    end
    return u < pop.pmiss[j]
end

# Whether records are made step by step (some trait has repeated records):
# then every trait's records are, single ones included.
_record_mode(pop::Population) = any(>(1), pop.nrec)

# The records of time step g: an animal alive at g gets a record of each
# trait it expresses, is at least the trait's age old and has fewer than
# its `records` of (one per step from that age on). The temporary
# residuals of one animal and step are drawn jointly, e = L_r z with z
# standard normal from the animal's stream of that step, so records of
# several traits made together have the residual covariance R0 (a
# subvector of e has the covariance of the matching submatrix). The
# phenotype becomes the mean of the animal's records.
function _record!(pop::Population, g)
    _record_mode(pop) || return pop
    k = ntraits(pop)
    z = zeros(k)
    for i = 1:length(pop)
        pop.generation[i] <= g < pop.exit[i] || continue
        nths = ntuple(j -> (t = pop.traits[j]; _expresses(t, pop.male[i]) ?
            g - pop.generation[i] - Int(t.age) + 1 : 0), k)
        any(j -> 1 <= nths[j] <= pop.nrec[j], 1:k) || continue
        rng = keyed_rng(SimKey(pop.key; generation = g, domain = Domain.PHENOTYPE, stage = _STAGE_RECORD, entity = i))
        for j = 1:k
            z[j] = randn(rng)
        end
        for j = 1:k
            nth = nths[j]
            1 <= nth <= pop.nrec[j] || continue
            e = 0.0
            for l = 1:j
                e += pop.Lr[j, l] * z[l]
            end
            _missing(pop, i, j, g) && continue
            y = pop.traits[j].μ + pop.tbv[i, j] + pop.pe[i, j] + _cg(pop, i, j, g) + e
            push!(pop.records.id, i); push!(pop.records.trait, j); push!(pop.records.step, g); push!(pop.records.y, y)
            c = pop.nobs[i, j] += 1
            m = pop.phenotype[i, j]
            pop.phenotype[i, j] = c == 1 ? y : m + (y - m) / c
        end
    end
    return pop
end

# The effect of the contemporary group of animal i's record of trait j:
# its herd and the step of the record. All traits' effects of a group come
# from one keyed stream, trait 1 first.
function _cg(pop::Population, i, j, t = pop.generation[i] + Int(pop.traits[j].age))
    pop.σcg[j] > 0 || return 0.0
    rng = keyed_rng(SimKey(pop.key; generation = t, domain = Domain.PHENOTYPE, stage = _STAGE_CG,
        entity = pop.herd[i]))
    z = 0.0
    for _ = 1:j
        z = randn(rng)
    end
    return pop.σcg[j] * z
end

# The contemporary group of each record of trait j among `ids`: herd and
# step of the record, numbered 1, 2, … in order of first appearance.
function _groups(pop::Population, ids, j; steps = pop.generation[ids] .+ Int(pop.traits[j].age))
    keys = [(pop.herd[i], Int32(s)) for (i, s) in zip(ids, steps)]
    lv = Dict{eltype(keys),Int}()
    return [get!(lv, x, length(lv) + 1) for x in keys]
end

_fits_groups(pop::Population) = any(>(0), pop.σcg)

# Fixed effects of trait j's records `rec` (rows of the view `v`): the
# contemporary groups when they are fitted (a sparse incidence matrix),
# else the overall mean.
function _fixed(pop::Population, v, rec, j; nrow = length(rec), rows = 1:length(rec))
    _fits_groups(pop) || (X = zeros(nrow, 1); X[rows, 1] .= 1; return X)
    g = _groups(pop, v.ids[rec], j)
    lv = zeros(Int, nrow)
    lv[rows] = g
    return incidence(lv, maximum(g; init = 0))
end

"""
    Scheme(; noff, nsire, ndam, evaluation = :pblup, selection = :truncation,
           K = NaN, λ = 0.0, κ = 0.0, ub = 1.0, chip = Int32[], blend = 0.05,
           max_sire_age = 0, max_dam_age = 0, reml = false, weights = Float64[],
           survival = Float64[], cull_unselected = false, genotyped_cohorts = 0,
           ocs = :fixed, ΔF = NaN, record_meiosis = false, report_reliability = false,
           mating = :random, pair_limit = 0, min_sire_age = 0, min_dam_age = 0,
           max_sire_offspring = 0, max_dam_offspring = 0, max_parity = 0)

A selection and mating scheme, applied every time step:

- `noff`: offspring per step (one per mating);
- `nsire`, `ndam`: parents selected by truncation, and the number of each
  sex used by random selection;
- `evaluation`: `:tbv` (true values, an oracle), `:phenotype` (own
  records), `:pblup` (animal model with the pedigree relationship matrix),
  `:gblup` (animal model with ``G^* = (1-w)G + wA``, `w = blend`, `G`
  VanRaden at the loci `chip` with base allele frequencies from the
  founders), `:ssgblup` (single-step, ``H^{-1}`` of Legarra et al. 2009
  and Christensen & Lund 2010 from `RelationshipMatrices.hinv`, with the
  animals of the last `genotyped_cohorts` cohorts genotyped at the loci
  `chip` (`0`: all animals), ``G`` as for `:gblup` and blended with
  ``A_{22}`` by `blend`; with every animal genotyped it equals `:gblup`).
  With several traits the animal model is multi-trait, with an intercept
  per trait and residuals correlated within animals;
- `weights`: the index weights ``w`` (one per trait; empty: all `1`).
  Candidates are ranked on ``I = w'\\hat a``, which with economic values
  as weights and multi-trait BLUP EBVs is the optimal index for the
  aggregate genotype ``H = w'g`` (Henderson 1963). With `:phenotype` it is
  ``w'y``, a fixed-weight index, not the optimal phenotypic index;
- `selection`: `:truncation` (best `nsire` males and `ndam` females on
  the index, equal numbers of offspring), `:ocs` (optimum contributions;
  see `ocs` below), `:rocs`
  (`Breeding.rocs` with `λ`, `κ`, the prediction error covariance of the
  candidates' index as `Ω`, and an upper bound `ub` on each parent's
  contribution, in `rocs` units where each sex contributes ½;
  `ub = 0.5/nsire` with `λ = κ = 0` is truncation selection), or
  `:random` (`nsire` and `ndam` chosen at random, a neutral reference);
- `max_sire_age`, `max_dam_age`: the oldest age, in time steps, at which an
  animal can be selected as a parent (age 0: born in the current step).
  `0` for both gives discrete generations; larger values give overlapping
  generations, with all eligible ages competing on one criterion. With
  `:ocs` and `:rocs` the contributions are those of the candidates to the
  next cohort, without the age-class weighting of Meuwissen & Sonesson
  (1998) or Bijma & Woolliams (2000); the constraint then does not bound
  the long-term rate of inbreeding exactly;
- `reml`: if `true`, the covariance components of `:pblup` and `:gblup`
  are re-estimated by AI-REML (Breeding.jl) from all records before every
  evaluation, starting from the previous estimates, once every trait has
  a record on an animal with a known parent (before that, as with the
  founders' records alone, genetic and residual components are not
  separable and the true values are used); otherwise the true values of
  the founders are used. With all records and the full pedigree, REML
  estimates the base-population components even under selection
  (Sorensen & Kennedy 1984; infinitesimal model). Several traits are
  estimated with AI steps on the Cholesky factors (`param = :cholesky` of
  `Breeding.reml`), which keeps the estimates positive definite; early
  estimates may lie on the boundary (e.g. a genetic correlation near ±1
  after one cohort);
- `ocs`: the constraint of `:ocs`. `:fixed`: `Breeding.TM1997` with
  ``c'Ac \\le K`` at every step (average coancestry ``K/2``). The others
  hold the rate of inbreeding at the target `ΔF` per generation
  (see [`ocs_contributions`](@ref)): `:cohort` constrains the mean kinship
  of the newborn cohort; `:age_year` the mean kinship of the population
  at the next step, age classes weighted by their remaining expected
  contributions (Meuwissen & Sonesson 1998; optiSel), with the per-step
  increase set by the estimated generation interval; `:age_generation`
  the same with the generation interval of the contributions being chosen
  (a per-generation constraint, Sonesson et al. 2000). With discrete
  generations the three coincide;
- `mating`: how the slots of the selected sires and dams (one per
  offspring) are paired, see [`mating_plan`](@ref): `:random`,
  `:min_coancestry` (the pairing of least total coancestry),
  `:assortative` or `:disassortative` (by the index);
- `min_sire_age`, `min_dam_age`: the youngest age at which an animal can
  be a parent (time steps; `0`: from birth), at most the maximum ages;
  the founders, the adult base population, are parents from step 0;
- `max_sire_offspring`, `max_dam_offspring`: at most this many offspring
  per parent and time step (`0`: no limit), e.g. a litter size, or
  natural service against AI. Offspring numbers are capped and the
  rest spread over the other parents in proportion to their
  contributions; an `ArgumentError` when the selected parents cannot
  have `noff` offspring under the limits;
- `max_parity`: at most this many time steps with offspring per dam
  (`0`: no limit; `pop.parity` counts them);
- `pair_limit`: at most this many offspring per sire–dam pair (`0`: no
  limit; `1`: no full sibs, a factorial-type design), with `:random` or
  `:min_coancestry`;
- `record_meiosis`: keep the crossover records of every cohort in
  `pop.meioses`, for realized IBD from founder segments
  ([`founder_segments`](@ref)); about 8 bytes per crossover;
- `report_reliability`: report the mean reliability of the candidates of
  each step from the prediction error variances of the evaluation
  (`:pblup`, `:gblup`, `:ssgblup`; see [`step!`](@ref)). It costs a
  solve of the mixed-model equations per candidate and trait, so it is
  off by default (the column is then `NaN`), and the resource estimate
  does not include it (on 5500 animals and 10 steps of PBLUP it took
  0.3 s more than the 0.4 s of the run);
- `survival`: `survival[a + 1]` is the probability that an animal alive
  at age `a` survives to age `a + 1` (the last entry for all older ages;
  empty: no deaths). Deaths are drawn between steps (`Domain.SURVIVAL`),
  independently of genotype and records (missing at random), so they
  bias neither BLUP nor REML;
- `cull_unselected`: if `true`, eligible animals of age 1 or more that
  are not selected as parents at a step are culled. Young candidates
  (age 0) are kept for the next step. Culling depends on the EBV, so the
  animals left are a selected sample; BLUP and REML with all records and
  the pedigree account for that;
- `window`: if positive, `:pblup`, `:gblup` and `:ssgblup` use the records
  of the animals born in the last `window` time steps only (the current
  cohort and the `window - 1` before it), on the pedigree of these animals
  and their ancestors traced for `ancestor_generations` generations (a
  [`PedigreeView`](@ref)); the other animals get no EBV (`NaN`). It must
  exceed the maximum ages, so that every candidate is evaluated. `0`
  (default): all records and the whole pedigree. Under selection,
  discarding the records of earlier cohorts and the parents beyond the
  traced depth biases BLUP (Henderson 1975), so a window trades accuracy
  for the size of the equations; REML with `reml = true` is then
  estimated from the window's records too;
- `ancestor_generations`: generations of ancestors traced above the
  animals of the window (default: all; then the relationships among the
  window's animals are exact, and only records are discarded);
- `retain`: `:all` (default) keeps every genome. `:active` drops, after
  each step, the genomes no later step can use: those of animals that
  can no longer be parents (dead, culled, or past the maximum age of
  their sex) and are not genotyped by `:ssgblup` (outside the last
  `genotyped_cohorts` cohorts). The founders are kept (the base allele
  frequencies). Pedigree and records are kept for every animal; memory
  then grows with the animals that are still candidates, not with all
  animals born. Not with `:gblup` or `:ssgblup` with every animal
  genotyped, which need every genome. A dropped genome cannot be read
  later ([`genomes`](@ref) throws), so summaries of past cohorts that
  need genomes are taken during the run.

Dead and culled animals are not candidates, and a record taken at an age
the animal did not reach is missing.
"""
Base.@kwdef struct Scheme
    noff::Int
    nsire::Int
    ndam::Int
    evaluation::Symbol = :pblup
    selection::Symbol = :truncation
    K::Float64 = NaN
    λ::Float64 = 0.0
    κ::Float64 = 0.0
    ub::Float64 = 1.0
    chip::Vector{Int32} = Int32[]
    blend::Float64 = 0.05
    max_sire_age::Int = 0
    max_dam_age::Int = 0
    reml::Bool = false
    weights::Vector{Float64} = Float64[]
    survival::Vector{Float64} = Float64[]
    cull_unselected::Bool = false
    genotyped_cohorts::Int = 0
    ocs::Symbol = :fixed
    ΔF::Float64 = NaN
    record_meiosis::Bool = false
    window::Int = 0
    ancestor_generations::Int = typemax(Int)
    retain::Symbol = :all
    report_reliability::Bool = false
    mating::Symbol = :random
    pair_limit::Int = 0
    min_sire_age::Int = 0
    min_dam_age::Int = 0
    max_sire_offspring::Int = 0
    max_dam_offspring::Int = 0
    max_parity::Int = 0
    function Scheme(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, chip, blend,
        max_sire_age, max_dam_age, reml, weights, survival, cull_unselected, genotyped_cohorts, ocs, ΔF,
        record_meiosis, window, ancestor_generations, retain, report_reliability, mating, pair_limit,
        min_sire_age, min_dam_age, max_sire_offspring, max_dam_offspring, max_parity)
        0 <= min_sire_age <= max_sire_age && 0 <= min_dam_age <= max_dam_age ||
            throw(ArgumentError("minimum ages must lie in 0:maximum age"))
        max_sire_offspring >= 0 && max_dam_offspring >= 0 && max_parity >= 0 ||
            throw(ArgumentError("offspring and parity limits must be non-negative"))
        mating in (:random, :min_coancestry, :assortative, :disassortative) ||
            throw(ArgumentError("mating must be :random, :min_coancestry, :assortative or :disassortative"))
        pair_limit >= 0 || throw(ArgumentError("pair_limit must be non-negative"))
        pair_limit > 0 && mating in (:assortative, :disassortative) &&
            throw(ArgumentError("pair_limit applies to :random and :min_coancestry mating"))
        noff > 1 && nsire > 0 && ndam > 0 || throw(ArgumentError("noff > 1, nsire > 0 and ndam > 0 are required"))
        max_sire_age >= 0 && max_dam_age >= 0 || throw(ArgumentError("maximum ages must be non-negative"))
        reml && !(evaluation in (:pblup, :gblup, :ssgblup)) &&
            throw(ArgumentError("reml = true needs evaluation :pblup, :gblup or :ssgblup"))
        evaluation in (:tbv, :phenotype, :pblup, :gblup, :ssgblup) ||
            throw(ArgumentError("evaluation must be :tbv, :phenotype, :pblup, :gblup or :ssgblup"))
        all(p -> 0 <= p <= 1, survival) || throw(ArgumentError("survival probabilities must lie in [0, 1]"))
        genotyped_cohorts >= 0 || throw(ArgumentError("genotyped_cohorts must be non-negative"))
        selection in (:truncation, :ocs, :rocs, :random) ||
            throw(ArgumentError("selection must be :truncation, :ocs, :rocs or :random"))
        ocs in (:fixed, :cohort, :age_year, :age_generation) ||
            throw(ArgumentError("ocs must be :fixed, :cohort, :age_year or :age_generation"))
        selection == :ocs && ocs == :fixed && !(K > 0) && throw(ArgumentError(":ocs needs a constraint K > 0"))
        selection == :ocs && ocs != :fixed && !(0 < ΔF < 1) &&
            throw(ArgumentError("ocs = :$ocs needs a target 0 < ΔF < 1 per generation"))
        evaluation in (:gblup, :ssgblup) && isempty(chip) &&
            throw(ArgumentError(":$evaluation needs chip loci"))
        0 <= blend < 1 || throw(ArgumentError("blend must lie in [0, 1)"))
        0 < ub <= 1 || throw(ArgumentError("ub must lie in (0, 1]"))
        all(isfinite, weights) && (isempty(weights) || any(!iszero, weights)) ||
            throw(ArgumentError("weights must be finite and not all zero"))
        retain in (:all, :active) || throw(ArgumentError("retain must be :all or :active"))
        retain == :active && (evaluation == :gblup || evaluation == :ssgblup && genotyped_cohorts == 0) &&
            throw(ArgumentError("retain = :active needs the genomes of a bounded set of animals; :$evaluation " *
                                "with every animal genotyped uses all of them (use :ssgblup with genotyped_cohorts > 0)"))
        window >= 0 && ancestor_generations >= 0 ||
            throw(ArgumentError("window and ancestor_generations must be non-negative"))
        window == 0 || window > max(max_sire_age, max_dam_age) || throw(ArgumentError(
            "window = $window leaves candidates of age $(max(max_sire_age, max_dam_age)) unevaluated; " *
            "it must exceed the maximum ages"))
        new(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, Int32.(chip), blend,
            max_sire_age, max_dam_age, reml, Float64.(weights), Float64.(survival), cull_unselected,
            genotyped_cohorts, ocs, ΔF, record_meiosis, window, ancestor_generations, retain,
            report_reliability, mating, pair_limit, min_sire_age, min_dam_age, max_sire_offspring,
            max_dam_offspring, max_parity)
    end
end

"""
    index_weights(pop::Population, scheme::Scheme) -> Vector{Float64}

The index weights of `scheme` for the traits of `pop` (all `1` when the
scheme gives none).
"""
function index_weights(pop::Population, sc::Scheme)
    isempty(sc.weights) && return ones(ntraits(pop))
    length(sc.weights) == ntraits(pop) ||
        throw(ArgumentError("the scheme has $(length(sc.weights)) weights for $(ntraits(pop)) traits"))
    return sc.weights
end

# ---- evaluation --------------------------------------------------------

"""
    genomes(pop::Population[, ids]) -> Haplotype

The phased genomes of the animals `ids` (default: all), two haplotypes
per animal in the order of `ids`. An `ArgumentError` if the genome of one
of them was dropped (`retain = :active` in [`Scheme`](@ref)); see
[`retained`](@ref).
"""
genomes(pop::Population, ids = 1:length(pop)) = _haplotypes(pop, ids)

"""
    retained(pop::Population) -> Vector{Int}

The ids of the animals whose genomes are kept.
"""
retained(pop::Population) = findall(>(0), pop.slot)

function _haplotypes(pop::Population, ids)
    nlc = pop.blocks[1].nlc
    out = Haplotype(nlc, 2length(ids))
    for (k, i) in enumerate(ids)
        g = pop.generation[i] + 1
        j = Int(pop.slot[i])
        j > 0 || throw(ArgumentError("the genome of animal $i was dropped (retain = :active)"))
        for h = 0:1
            copyto!(packed(out), (2k - 2 + h) * size(packed(out), 1) + 1, packed(pop.blocks[g]),
                (2j - 2 + h) * size(packed(pop.blocks[g]), 1) + 1, size(packed(out), 1))
        end
    end
    return out
end

"""
    genotyped(pop::Population, scheme::Scheme, t = current_step(pop)) -> Vector{Int}

The animals genotyped for `:gblup` (all) and `:ssgblup` (those born in
the last `scheme.genotyped_cohorts` cohorts at step `t`, or all when it
is `0`).
"""
function genotyped(pop::Population, sc::Scheme, t::Integer = current_step(pop))
    (sc.evaluation == :gblup || sc.genotyped_cohorts == 0) && return collect(1:length(pop))
    filter(i -> t - pop.generation[i] < sc.genotyped_cohorts, 1:length(pop))
end

pedigree_view(pop::Population, ids; generations::Integer = typemax(Int)) =
    pedigree_view(pop.sire, pop.dam, ids; generations = generations)

"""
    evaluation_view(pop::Population, scheme::Scheme, t = current_step(pop)) -> PedigreeView

The pedigree on which `scheme` evaluates at step `t`: every animal when
`scheme.window` is `0`, otherwise the animals born in the last `window`
steps and their ancestors traced for `scheme.ancestor_generations`
generations. The animal effects of the equations returned by
[`evaluate!`](@ref) are in the rows of this view.
"""
function evaluation_view(pop::Population, sc::Scheme, t::Integer = current_step(pop))
    sc.window == 0 && return PedigreeView(pop.sire, pop.dam, 1:length(pop))
    ids = filter(i -> t - pop.generation[i] < sc.window, 1:length(pop))
    return pedigree_view(pop, ids; generations = sc.ancestor_generations)
end

function _kinv(pop::Population, sc::Scheme, v::PedigreeView)
    sc.evaluation == :pblup && return ainv(v.sire, v.dam)
    ids = intersect(genotyped(pop, sc), v.ids)            # both ascending
    H = _haplotypes(pop, ids)
    base = pop.blocks[1]                                  # the founders
    p = zeros(base.nlc)                                   # grm wants one per locus
    p[sc.chip] = vec(sum(dosage(base, sc.chip); dims = 2)) ./ base.nhp
    G = grm(H; p = p, loci = sc.chip)
    sc.evaluation == :ssgblup && return hinv(v.sire, v.dam, G, id_to_row(v, ids); delta = sc.blend)
    sc.blend > 0 && (G = (1 - sc.blend) .* G .+ sc.blend .* nrm(v.sire, v.dam))
    inv(cholesky(Symmetric(G)))
end

"""
    evaluate!(pop::Population, scheme::Scheme) -> Union{Nothing, MME}

Fill `pop.ebv` for every animal and trait by the scheme's evaluation at
the current time step. For `:pblup`, `:gblup` and `:ssgblup` the mixed-model
equations over the records available (see [`observed`](@ref)) are built
with Breeding.jl and solved directly: ``y = 1\\mu + Za + e`` with
``\\sigma^2_e`` scaling for one trait; for several traits, records stacked
trait-major with an intercept per recorded trait,
``\\operatorname{var}(a) = G_0 \\otimes K`` and ``R_0`` among the records
of one animal (missing records handled exactly). The covariance
components are `pop.G0` and `pop.R0`, first re-estimated by REML if
`scheme.reml`. The equations are returned for prediction error variances;
their animal effects are in the rows of [`evaluation_view`](@ref), which
with `scheme.window > 0` covers the recent cohorts and their ancestors
only (the EBV of the other animals is `NaN`).
`:phenotype` needs every trait recorded at birth in both sexes.
"""
function evaluate!(pop::Population, sc::Scheme)
    n, k = length(pop), ntraits(pop)
    if sc.evaluation == :tbv
        pop.ebv .= pop.tbv
        return nothing
    elseif sc.evaluation == :phenotype
        all(t -> t.sex == 2 && t.age == 0, pop.traits) ||
            throw(ArgumentError(":phenotype needs every trait recorded at birth in both sexes"))
        all(iszero, pop.pmiss) && all(iszero, pop.hmiss) ||
            throw(ArgumentError(":phenotype needs every record; missing records leave candidates without one"))
        pop.ebv .= pop.phenotype
        return nothing
    end
    v = evaluation_view(pop, sc)
    t = current_step(pop)
    obs = observed(pop)[v.ids, :]
    # records of the window's animals only, not of the ancestors traced
    sc.window > 0 && (obs[[t - pop.generation[i] >= sc.window for i in v.ids], :] .= false)
    any(obs) || throw(ArgumentError("no records are available at time step $t"))
    P = pop.phenotype[v.ids, :]
    n = length(v)
    Kinv = _kinv(pop, sc, v)
    # Genetic and residual components are not separable from records of
    # unrelated animals only (the founders at step 0): REML waits for every
    # trait to be recorded on an animal with a known parent.
    known = (v.sire .> 0) .| (v.dam .> 0)
    do_reml = sc.reml && all(j -> any(view(obs, :, j) .& known), 1:k)
    if k == 1 && pop.nrec[1] > 1
        m = _evaluate_repeated!(pop, sc, v, Kinv, t, obs, do_reml)
    elseif _record_mode(pop)
        m = _evaluate_records!(pop, sc, v, Kinv, t, do_reml)
    elseif k == 1
        rec = findall(view(obs, :, 1))
        y = P[rec, 1]
        X = _fixed(pop, v, rec, 1)
        effects = [RandomEffect(incidence(rec, n), Kinv, pop.G0[1, 1])]
        if do_reml
            r = reml(y, X, effects; σ²e = pop.R0[1, 1], tol = 1e-6, maxiter = 100)
            _converged(r, pop)
            pop.G0 .= r.G0[1]
            pop.R0[1, 1] = r.σ²e
            effects = [RandomEffect(incidence(rec, n), Kinv, pop.G0[1, 1])]
        end
        m = mme(y, X, effects; σ²e = pop.R0[1, 1])
    else
        y, X, Zs = _mt_design(P, obs, j -> (rec = findall(view(obs, :, j));
            _fixed(pop, v, rec, j; nrow = size(obs, 1), rows = rec)))
        if do_reml
            for j = 1:k, l = 1:j-1
                any(view(obs, :, j) .& view(obs, :, l)) || throw(ArgumentError(
                    "traits $(pop.traits[l].name) and $(pop.traits[j].name) are never recorded on " *
                    "the same animal: their residual covariance is not estimable; use reml = false"))
            end
            # on the Cholesky scale: early estimates often lie on the boundary
            # (G0 singular after one cohort), where natural-scale AI steps
            # are rejected and the EM fall-back crawls
            r = reml(y, X, [RandomEffect(Zs, Kinv, pop.G0)];
                residual = MultiTrait(obs, pop.R0), param = :cholesky, tol = 1e-6, maxiter = 100)
            _converged(r, pop)
            pop.G0 .= r.G0[1]
            pop.R0 .= r.R0[1]
        end
        m = mme(y, X, [RandomEffect(Zs, Kinv, pop.G0)]; Rinv = mt_rinv(pop.R0, obs))
    end
    x = solve_mme(m)
    a = first(m.blocks[2])
    pop.ebv .= NaN
    pop.ebv[v.ids, :] = reshape(x[a:a+k*n-1], n, k)
    return m
end

# Repeatability model of one trait: y = Xb + Za + Zpe + e over the
# repeated records of the view's animals made by step t (the window's
# animals only), Var(a) = σ²a K, Var(pe) = σ²pe I, Var(e) = σ²e I. The
# permanent environment is fitted when the simulation has one (σpe > 0).
function _evaluate_repeated!(pop::Population, sc::Scheme, v, Kinv, t, obs, do_reml)
    R = pop.records
    keep = [k for k in eachindex(R.id) if R.trait[k] == 1 && R.step[k] <= t && v.row[R.id[k]] > 0 &&
            obs[v.row[R.id[k]], 1]]
    ids = Int.(R.id[keep])
    rows = v.row[ids]
    n = length(v)
    y = R.y[keep]
    X = if _fits_groups(pop)
        g = _groups(pop, ids, 1; steps = R.step[keep])
        incidence(g, maximum(g; init = 0))
    else
        ones(length(y), 1)
    end
    Z = incidence(rows, n)
    fitpe = pop.σpe[1] > 0
    effs() = fitpe ? [RandomEffect(Z, Kinv, pop.G0[1, 1]), RandomEffect(Z, I, pop.PE0[1, 1])] :
             [RandomEffect(Z, Kinv, pop.G0[1, 1])]
    if do_reml
        r = reml(y, X, effs(); σ²e = pop.R0[1, 1], tol = 1e-6, maxiter = 100)
        _converged(r, pop)
        pop.G0 .= r.G0[1]
        fitpe && (pop.PE0[1, 1] = r.G0[2][1, 1])
        pop.R0[1, 1] = r.σ²e
    end
    return mme(y, X, effs(); σ²e = pop.R0[1, 1])
end

# Several traits with records made step by step: one row per record unit
# (animal and step) with the traits recorded in it, stacked trait-major;
# y = Xb + Za + Z_pe pe + e with Var(a) = G0 ⊗ K, Var(pe) = P0 ⊗ I over the
# traits with a permanent environment, and R0 among the records of one
# unit (missing records handled exactly). Records of the view's animals
# made by step t, of the window's animals only.
function _evaluate_records!(pop::Population, sc::Scheme, v, Kinv, t, do_reml)
    k = ntraits(pop)
    R = pop.records
    keep = [q for q in eachindex(R.id) if R.step[q] <= t && v.row[R.id[q]] > 0 &&
            (sc.window == 0 || t - pop.generation[R.id[q]] < sc.window)]
    unit = Dict{Tuple{Int32,Int32},Int}()
    for q in keep
        get!(unit, (R.id[q], R.step[q]), length(unit) + 1)
    end
    U = length(unit)
    uid, ust = zeros(Int32, U), zeros(Int32, U)
    for ((i, s), u) in unit
        uid[u], ust[u] = i, s
    end
    Y = fill(NaN, U, k)
    for q in keep
        Y[unit[(R.id[q], R.step[q])], R.trait[q]] = R.y[q]
    end
    obs = .!isnan.(Y)
    n = length(v)
    y, _ = mt_stack(Y)
    function fixed(j)
        _fits_groups(pop) || return ones(U, 1)
        rows = findall(view(obs, :, j))
        isempty(rows) && return ones(U, 1)
        g = _groups(pop, Int.(uid[rows]), j; steps = ust[rows])
        lv = zeros(Int, U)
        lv[rows] = g
        return incidence(lv, maximum(g))
    end
    ot = findall(j -> any(view(obs, :, j)), 1:k)
    X = reduce(hcat, mt_blocks(obs, [fixed(j) for j = 1:k])[ot])
    Zs = mt_blocks(obs, [incidence(v.row[uid], n) for _ = 1:k])
    pej = findall(>(0), pop.σpe)
    effs() = isempty(pej) ? [RandomEffect(Zs, Kinv, pop.G0)] :
             [RandomEffect(Zs, Kinv, pop.G0), RandomEffect(Zs[pej], I, pop.PE0[pej, pej])]
    if do_reml
        for j = 1:k, l = 1:j-1
            any(view(obs, :, j) .& view(obs, :, l)) || throw(ArgumentError(
                "traits $(pop.traits[l].name) and $(pop.traits[j].name) are never recorded at the same step " *
                "on the same animal: their residual covariance is not estimable; use reml = false"))
        end
        r = reml(y, X, effs(); residual = MultiTrait(obs, pop.R0), param = :cholesky, tol = 1e-6, maxiter = 100)
        _converged(r, pop)
        pop.G0 .= r.G0[1]
        isempty(pej) || (pop.PE0[pej, pej] = r.G0[2])
        pop.R0 .= r.R0[1]
    end
    return mme(y, X, effs(); Rinv = mt_rinv(pop.R0, obs))
end

function _converged(r, pop)
    r.iterations < 100 ||
        throw(ErrorException("REML did not converge in 100 iterates at time step $(current_step(pop))"))
end

# Multi-trait design over the records `obs` (n × k), stacked trait-major.
# Every trait keeps its genetic component, so a trait without records yet
# (e.g. recorded late in life) is predicted through G0; only the
# intercepts of traits with records enter X.
# `fixed(j)` is trait j's fixed-effect design over all n rows (only the
# rows with records are used); the default is the overall mean.
function _mt_design(P, obs, fixed = j -> ones(size(obs, 1), 1))
    n, k = size(obs)
    ot = findall(j -> any(view(obs, :, j)), 1:k)
    y, _ = mt_stack(ifelse.(obs, P, NaN))
    X = reduce(hcat, mt_blocks(obs, [fixed(j) for j = 1:k])[ot])
    Zs = mt_blocks(obs, [incidence(1:n, n) for _ = 1:k])
    return y, X, Zs
end

# ---- selection ---------------------------------------------------------

# Integer numbers of offspring for parents with contributions `w` (summing
# to 1 within the group): largest remainder, ties by ascending id.
function _allocate(ids, w, n; cap::Integer = 0)
    cap > 0 && return _allocate_capped(ids, w, n, cap)
    raw = w ./ sum(w) .* n
    k = floor.(Int, raw)
    rest = n - sum(k)
    order = sortperm(collect(zip(-(raw .- k), ids)))
    for t in order[1:rest]
        k[t] += 1
    end
    return k
end

# The allocation above with at most `cap` per parent: shares in proportion
# to `w`, those above `cap` set to `cap` and the rest spread over the others
# in proportion to their `w` (repeated until none exceeds it), then rounded
# by largest remainder without passing `cap`, ties by ascending id.
function _allocate_capped(ids, w, n, cap)
    length(ids) * cap >= n || throw(ArgumentError("$(length(ids)) parents with at most $cap offspring each " *
        "cannot have $n offspring; select more parents or raise the limit"))
    raw = zeros(length(ids))
    free = trues(length(ids))
    rem = Float64(n)
    while true
        tot = sum(w[free])
        tot > 0 || throw(ArgumentError("the parents with a contribution cannot have $n offspring under " *
            "the limit of $cap each"))
        raw[free] = w[free] ./ tot .* rem
        over = free .& (raw .> cap)
        any(over) || break
        raw[over] .= cap
        free .&= .!over
        rem = n - sum(raw[.!free])
    end
    k = floor.(Int, raw)
    rest = n - sum(k)
    order = filter(t -> k[t] < cap, sortperm(collect(zip(-(raw .- k), ids))))
    length(order) >= rest || throw(ArgumentError("cannot allocate $n offspring under the limit of $cap each"))
    for t in order[1:rest]
        k[t] += 1
    end
    return k
end

# Indices of the `m` best by (criterion descending, id ascending).
_best(ids, crit, m) = ids[sortperm(collect(zip(-crit, ids)))[1:min(m, length(ids))]]

# Prediction error covariance of the index w'â among the candidates:
# Σ_jl w_j w_l Γ_jl, Γ_jl the block of C⁻¹ (times the scale) of the
# genetic components j and l.
function _pev_index(m::MME, cand, w, n)
    k = length(w)
    k == 1 && return w[1]^2 .* pev_candidates(m, 1, cand)
    a = first(m.blocks[2]) - 1
    rows = [a + (j - 1) * n + c for j = 1:k for c in cand]
    F = cholesky(Symmetric(m.lhs))
    E = zeros(size(m.lhs, 1), length(rows))
    for (c, r) in enumerate(rows)
        E[r, c] = 1
    end
    Γ = (F\E)[rows, :] .* m.scale
    nc = length(cand)
    Ω = zeros(nc, nc)
    for j = 1:k, l = 1:k
        Ω .+= w[j] * w[l] .* view(Γ, (j-1)*nc+1:j*nc, (l-1)*nc+1:l*nc)
    end
    return (Ω .+ Ω') ./ 2
end

# living males and females of an age to be parents at step g (between the
# minimum and maximum ages; the founders, the adult base population, from
# step 0), dams below the maximum parity
function _eligible(pop::Population, sc::Scheme, g)
    young(a) = pop.first[max(0, g - a)+1]:length(pop)
    old(i, a) = pop.generation[i] == 0 || g - pop.generation[i] >= a
    males = filter(i -> pop.male[i] && pop.exit[i] > g && old(i, sc.min_sire_age), young(sc.max_sire_age))
    females = filter(i -> !pop.male[i] && pop.exit[i] > g && old(i, sc.min_dam_age) &&
        (sc.max_parity == 0 || pop.parity[i] < sc.max_parity), young(sc.max_dam_age))
    return males, females
end

"""
    select_parents(pop::Population, scheme::Scheme, m) -> (sires, nsire_off, dams, ndam_off, info)

Parents and their numbers of offspring (each list summing to
`scheme.noff`), among the living males of age at most
`scheme.max_sire_age` and females of age at most `scheme.max_dam_age` at
the current time step (the youngest cohort only, with the default ages
0), ranked on the index
``w'\\hat a`` (see [`Scheme`](@ref)). `m` is the `MME` returned by
[`evaluate!`](@ref) (`nothing` for `:tbv` and `:phenotype`); `:rocs` takes
the prediction error covariance of the candidates' index from it. Every
ranking breaks ties by ascending id (invariant I1a). `info` is empty
except for `:ocs`, where it holds the bound, the constrained kinship
reached, the generation interval used and whether the bound could be met
(`ocs_bound`, `ocs_kinship`, `ocs_L`, `ocs_feasible`). With `ocs = :fixed`
the bound is `K/2` and the kinship ``c'Ac/2`` of `Breeding.TM1997`'s
contributions, on the same kinship scale, and `ocs_L` is `NaN`; when `K`
cannot be met, TM1997 returns the minimum-coancestry contributions and
`ocs_feasible` is `false`.
"""
function select_parents(pop::Population, sc::Scheme, m)
    g = current_step(pop)
    w = index_weights(pop, sc)
    index = pop.ebv * w
    males, females = _eligible(pop, sc, g)
    (isempty(males) || isempty(females)) && throw(ArgumentError(
        "no living $(isempty(males) ? "males" : "females") are eligible at step $g; " *
        "raise the maximum ages or the survival, or lower the minimum ages or the parity limit"))
    cand = sort!([males; females])
    if sc.selection in (:truncation, :random)
        crit = sc.selection == :truncation ? index : _random_criterion(pop, g, cand)
        s = _best(males, crit[males], sc.nsire)
        d = _best(females, crit[females], sc.ndam)
        s = sort(s)
        d = sort(d)
        return s, _allocate(s, ones(length(s)), sc.noff; cap = sc.max_sire_offspring), d,
            _allocate(d, ones(length(d)), sc.noff; cap = sc.max_dam_offspring), NamedTuple()
    end
    A = nrm(pop.sire, pop.dam, cand)
    u = index[cand]
    sex = Int.(pop.male[cand])
    info = NamedTuple()
    c = if sc.selection == :ocs && sc.ocs == :fixed
        r = TM1997(u, sex, A, sc.K)
        # c'Ac with ½ per sex is twice the mean kinship of the contributions
        info = (ocs_bound = sc.K / 2, ocs_kinship = r.cAc / 2, ocs_L = NaN, ocs_feasible = feasible(r))
        2 .* r
    elseif sc.selection == :ocs
        cc, info = _dynamic_ocs(pop, sc, g, cand, u)
        2 .* cc
    else
        Ω = if m === nothing
            zeros(length(cand), length(cand))
        else
            v = evaluation_view(pop, sc)
            _pev_index(m, id_to_row(v, cand), w, length(v))
        end
        2 .* rocs(u, A; Ω = Ω, sex = sex, λ = sc.λ, κ = sc.κ, ub = sc.ub).w
    end
    ms = findall(k -> pop.male[cand[k]] && c[k] > 1e-9, eachindex(cand))
    fs = findall(k -> !pop.male[cand[k]] && c[k] > 1e-9, eachindex(cand))
    s, d = cand[ms], cand[fs]
    ns = _allocate(s, c[ms], sc.noff; cap = sc.max_sire_offspring)
    nd = _allocate(d, c[fs], sc.noff; cap = sc.max_dam_offspring)
    keep_s, keep_d = ns .> 0, nd .> 0
    return s[keep_s], ns[keep_s], d[keep_d], nd[keep_d], info
end

# A keyed uniform per candidate (entity = id), NaN for the other animals:
# drawn for the candidates only, as their values do not depend on which
# other animals exist.
function _random_criterion(pop::Population, g, cand)
    crit = fill(NaN, length(pop))
    for i in cand
        crit[i] = rand(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.SELECTION,
            stage = _STAGE_RANDOM_SELECTION, entity = i)))
    end
    return crit
end

"""
    mating_plan(pop, sires, nsire_off, dams, ndam_off; mating = :random, pair_limit = 0, index = nothing)
        -> Matrix{Int}

Pair the slots of the sires `sires` (`nsire_off[k]` offspring each) and
the dams `dams` (`ndam_off[k]`): row `i` is mating `i`, as ids. Every
parent keeps its number of offspring.

- `:random`: a keyed shuffle of the dam slots (`Domain.MATING`). With
  `pair_limit > 0`, the pairing of least total cost under the limit, with
  a keyed uniform cost per sire–dam pair: random, but not uniform over
  the pairings that meet the limit.
- `:min_coancestry`: the pairing of least total coancestry
  ``\\sum \\theta_{sd}`` (from the pedigree, ``\\theta = A/2``), exactly
  ([`transport`](@ref)); the mean coancestry of the matings is the
  expected inbreeding of the offspring.
- `:assortative`: sire slots by descending `index`, paired with dam slots
  by descending `index` (best with best); `:disassortative`: with dam
  slots by ascending `index`. Ties by ascending id.

`pair_limit > 0` allows at most that many offspring per sire–dam pair
(`:random`, `:min_coancestry`); an `ArgumentError` when no pairing meets
it. Rows of the pairings by cost are ordered by sire, then dam.

The numbers of offspring come from [`select_parents`](@ref), so a policy
applies to any selection: with `:ocs` or `:rocs`, the optimised
contributions fix each parent's offspring and `mating` pairs them (two
stages, contributions then mate allocation; not a joint optimisation of
both). Minimum-coancestry mating lowers the inbreeding of the offspring
below the coancestry of the contributions in the short term; the
asymptotic rate of inbreeding is set by the contributions.
"""
function mating_plan(pop::Population, s, ns, d, nd; mating::Symbol = :random, pair_limit::Integer = 0,
    index = nothing)
    sum(ns) == sum(nd) || throw(ArgumentError("sire and dam slots differ: $(sum(ns)) and $(sum(nd))"))
    g = current_step(pop)
    if mating == :random && pair_limit == 0
        sl = reduce(vcat, [fill(Int(s[k]), ns[k]) for k in eachindex(s)])
        dl = reduce(vcat, [fill(Int(d[k]), nd[k]) for k in eachindex(d)])
        keyed_shuffle!(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.MATING, stage = _STAGE_PAIRING)), dl)
        return hcat(sl, dl)
    elseif mating in (:assortative, :disassortative)
        index === nothing && throw(ArgumentError(":$mating mating needs the index"))
        pair_limit == 0 || throw(ArgumentError("pair_limit applies to :random and :min_coancestry mating"))
        os = sortperm(collect(zip(-index[s], s)))
        od = mating == :assortative ? sortperm(collect(zip(-index[d], d))) : sortperm(collect(zip(index[d], d)))
        sl = reduce(vcat, [fill(Int(s[k]), ns[k]) for k in os])
        dl = reduce(vcat, [fill(Int(d[k]), nd[k]) for k in od])
        return hcat(sl, dl)
    end
    C = if mating == :min_coancestry
        A = nrm(pop.sire, pop.dam, [s; d])
        A[1:length(s), length(s)+1:end] ./ 2
    elseif mating == :random
        reduce(vcat, [begin
            rng = keyed_rng(SimKey(pop.key; generation = g, domain = Domain.MATING, stage = _STAGE_PAIR_COST,
                entity = Int(i)))
            rand(rng, length(d))'
        end for i in s])
    else
        throw(ArgumentError("mating must be :random, :min_coancestry, :assortative or :disassortative"))
    end
    F = transport(ns, nd, C; cap = pair_limit > 0 ? pair_limit : typemax(Int))
    rows = [(Int(s[a]), Int(d[b])) for a in eachindex(s) for b in eachindex(d) for _ = 1:F[a, b]]
    return [first.(rows) last.(rows)]
end

# ---- one generation ----------------------------------------------------

"""
    step!(pop::Population, scheme::Scheme, mm::MeiosisMap) -> NamedTuple

One time step: evaluate, select, mate, and add the offspring (with
genetic values, sexes and phenotypes) as the next cohort. Returns the
summary of the current cohort (see [`generation_summary`](@ref)) with the
accuracy of the evaluation in it (of the index, and with several traits
of each trait as `accuracy_<name>`), the parents used, the generation
interval `L` (mean age of the parents at the birth of their offspring,
over both sexes), the realized selection intensities `intensity_sire`
and `intensity_dam` (the mean index of a sex's parents, weighted by their
numbers of offspring, less the mean index of the eligible candidates of
that sex, in standard deviations of the candidates' index), the mean
reliability ``R^2 = 1 - \\mathrm{PEV}/\\sigma^2_a`` of the candidates'
index (`reliability`, and `reliability_<name>` per trait with several
traits; ``\\sigma^2_a`` of the index ``w'G_0w`` from the base variances
the evaluation used, `NaN` unless the scheme sets `report_reliability`
and the evaluation solves mixed-model equations) and the
variances used by the evaluation, `va` and `ve` (one trait; with
repeated records also `vpe`, the permanent-environment variance, and
`ve` the temporary residual) or `va_<name>` and `ve_<name>`.
"""
function step!(pop::Population, sc::Scheme, mm::MeiosisMap)
    g = current_step(pop)
    _record!(pop, g)
    w = index_weights(pop, sc)
    m = evaluate!(pop, sc)
    cand = generation_ids(pop, g)
    acc = merge(_accuracies(pop, cand, w), _reliabilities(pop, sc, m, cand, w))
    elig = _eligible(pop, sc, g)
    s, ns, d, nd, info = select_parents(pop, sc, m)
    plan = mating_plan(pop, s, ns, d, nd; mating = sc.mating, pair_limit = sc.pair_limit, index = pop.ebv * w)
    L = mean(g + 1 .- pop.generation[plan])
    ints = _intensities(pop.ebv * w, elig, plan)
    _cull_and_survive!(pop, sc, g, [s; d])
    # the parents' haplotypes, gathered from their cohorts; offspring
    # depend on the parents' alleles and the keys only, not on this order
    ups = sort!(unique(plan))
    parents = _haplotypes(pop, ups)
    off = Haplotype(parents.nlc, 2sc.noff)
    local_plan = map(i -> searchsortedfirst(ups, i), plan)
    if sc.record_meiosis
        _, r = make_progeny!(off, parents, local_plan, mm, SimKey(pop.key; generation = g + 1); record = true)
        push!(pop.meioses, MeiosisRecord(ups[r.parent], r.phase, r.ptr, r.pos))
    else
        make_progeny!(off, parents, local_plan, mm, SimKey(pop.key; generation = g + 1))
        push!(pop.meioses, nothing)
    end
    n0 = length(pop)
    ids = (n0+1):(n0+sc.noff)
    k = ntraits(pop)
    append!(pop.sire, Int32.(plan[:, 1]))
    append!(pop.dam, Int32.(plan[:, 2]))
    append!(pop.male, falses(sc.noff))
    append!(pop.generation, fill(Int32(g + 1), sc.noff))
    append!(pop.herd, pop.herd[plan[:, 2]])
    pop.tbv = vcat(pop.tbv, reduce(hcat, [tbv(off, q) for q in pop.qtl]))
    pop.phenotype = vcat(pop.phenotype, fill(NaN, sc.noff, k))
    pop.ebv = vcat(pop.ebv, fill(NaN, sc.noff, k))
    pop.pe = vcat(pop.pe, zeros(sc.noff, k))
    pop.nobs = vcat(pop.nobs, zeros(Int32, sc.noff, k))
    append!(pop.parity, zeros(Int32, sc.noff))
    for i in unique(plan[:, 2])                          # dams with offspring this step
        pop.parity[i] += 1
    end
    append!(pop.exit, fill(typemax(Int32), sc.noff))
    push!(pop.blocks, off)
    push!(pop.first, n0 + 1)
    append!(pop.slot, Int32.(1:sc.noff))
    _sexes!(pop, ids, g + 1)
    _phenotypes!(pop, ids, g + 1)
    segs = sc.record_meiosis ? founder_segments(pop, mm) : nothing
    out = merge(generation_summary(pop, g; weights = w, segs = segs), acc,
        (nsire = length(s), ndam = length(d), L = L), ints, _vc(pop), info)
    sc.retain == :active && _drop_genomes!(pop, sc)
    return out
end

# With retain = :active, after a step: drop the genomes no later step of
# the scheme can use. At step t an animal's genome is used as a parent's
# (alive and at most the maximum age of its sex) or as a genotype of
# ssGBLUP (born in the last genotyped_cohorts cohorts); both stop for
# good once they stop, since ages only grow and exits are final. The
# founders are kept: they give the base allele frequencies of G and of
# Fhet. A cohort's block is compacted to the genomes kept.
function _needs_genome(pop::Population, sc::Scheme, i, t)
    gen = pop.generation[i]
    a = t - gen
    gen == 0 && return true
    pop.exit[i] > t && a <= (pop.male[i] ? sc.max_sire_age : sc.max_dam_age) && return true
    return sc.evaluation == :ssgblup && a < sc.genotyped_cohorts
end

function _drop_genomes!(pop::Population, sc::Scheme)
    t = current_step(pop)
    nlc = pop.blocks[1].nlc
    for c = 1:t-1                                          # cohorts 0 and t are kept
        b = pop.blocks[c+1]
        b === nothing && continue
        ids = generation_ids(pop, c)
        keep = filter(i -> pop.slot[i] > 0 && _needs_genome(pop, sc, i, t), ids)
        length(keep) == count(i -> pop.slot[i] > 0, ids) && continue
        if isempty(keep)
            pop.blocks[c+1] = nothing
        else
            pop.blocks[c+1] = _haplotypes(pop, keep)
        end
        pop.slot[ids] .= 0
        pop.slot[keep] .= Int32.(1:length(keep))
    end
    return pop
end

_no_ocs(sc::Scheme) = sc.selection == :ocs ?
    (ocs_bound = NaN, ocs_kinship = NaN, ocs_L = NaN, ocs_feasible = true) : NamedTuple()

# Between steps g and g + 1: culling of unselected adults, then deaths.
# An animal alive at g of age a survives to g + 1 with probability
# survival[min(a + 1, end)], a keyed draw per animal and step.
function _cull_and_survive!(pop::Population, sc::Scheme, g, parents)
    if sc.cull_unselected
        males, females = _eligible(pop, sc, g)
        sel = Set(parents)
        for i in Iterators.flatten((males, females))
            g - pop.generation[i] >= 1 && !(i in sel) && (pop.exit[i] = g + 1)
        end
    end
    isempty(sc.survival) && return
    for i in alive(pop, g)
        a = g - pop.generation[i]
        q = sc.survival[min(a + 1, end)]
        u = rand(keyed_rng(SimKey(pop.key; generation = g + 1, domain = Domain.SURVIVAL, entity = i)))
        u < q || (pop.exit[i] = g + 1)
    end
end

# The realized selection intensity of each sex: the mean index of the
# parents, each counted once per offspring, less the mean index of the
# eligible candidates of its sex, in standard deviations of the
# candidates' index (n - 1 divisor). On the selection criterion (the EBV
# index), so it is the standardized selection differential the scheme
# achieved, after any contribution constraint (OCS) and death.
function _intensities(index, (males, females), plan)
    i(cand, par) = (mean(index[par]) - mean(index[cand])) / std(index[cand])
    return (intensity_sire = i(males, plan[:, 1]), intensity_dam = i(females, plan[:, 2]))
end

function _accuracies(pop::Population, ids, w)
    k = ntraits(pop)
    acc = (accuracy = cor(pop.ebv[ids, :] * w, pop.tbv[ids, :] * w),)
    k == 1 && return acc
    names = Tuple(Symbol("accuracy_", t.name) for t in pop.traits)
    merge(acc, NamedTuple{names}(Tuple(cor(pop.ebv[ids, j], pop.tbv[ids, j]) for j = 1:k)))
end

# Per candidate (rows `rows` of an evaluation of `n` animals and `k`
# traits), the k × k prediction error covariance of its breeding values,
# from the MME `m`: the diagonal blocks of the inverse of the coefficient
# matrix, solved for 64 candidates at a time so the work array has at most
# 64k columns whatever the size of the cohort.
function _pev_blocks(m::MME, rows, n, k)
    a = first(m.blocks[2]) - 1
    L = m.lhs
    F = cholesky(Symmetric(L); check = false)
    fac = issuccess(F) ? F : lu(L)
    out = Array{Float64}(undef, k, k, length(rows))
    for chunk in Iterators.partition(eachindex(rows), 64)
        idx = [a + (j - 1) * n + rows[c] for c in chunk for j = 1:k]
        E = zeros(size(L, 1), length(idx))
        for (col, r) in enumerate(idx)
            E[r, col] = 1
        end
        S = fac \ E
        for (t, c) in enumerate(chunk)
            cols = (t-1)*k+1:t*k
            out[:, :, c] = S[idx[cols], cols] .* m.scale
        end
    end
    return out
end

# The mean reliability R² = 1 - PEV/σ²_a of the candidates `ids`: of the
# index w'u (σ² = w'G₀w) and, with several traits, of each trait
# (`reliability_<name>`). σ²_a is the base genetic variance the evaluation
# used (`va`), so R² is on the scale of the founders. NaN without an MME
# (`:tbv`, `:phenotype`) or unless the scheme asks for it.
function _reliabilities(pop::Population, sc::Scheme, m, ids, w)
    (sc.report_reliability && m !== nothing) || return _no_reliability(pop)
    k = ntraits(pop)
    v = evaluation_view(pop, sc)
    P = _pev_blocks(m, id_to_row(v, ids), length(v), k)
    G = pop.G0
    σ² = dot(w, G * w)
    r = (reliability = mean(1 - dot(w, view(P, :, :, i) * w) / σ² for i in axes(P, 3)),)
    k == 1 && return r
    names = Tuple(Symbol("reliability_", t.name) for t in pop.traits)
    merge(r, NamedTuple{names}(Tuple(mean(1 - P[j, j, i] / G[j, j] for i in axes(P, 3)) for j = 1:k)))
end

_no_reliability(pop::Population) = (k = ntraits(pop); k == 1 ? (reliability = NaN,) :
    merge((reliability = NaN,), NamedTuple{Tuple(Symbol("reliability_", t.name) for t in pop.traits)}(Tuple(fill(NaN, k)))))

_no_accuracy(pop::Population) = (k = ntraits(pop); k == 1 ? (accuracy = NaN,) :
    merge((accuracy = NaN,), NamedTuple{Tuple(Symbol("accuracy_", t.name) for t in pop.traits)}(Tuple(fill(NaN, k)))))

function _vc(pop::Population)
    ntraits(pop) == 1 && pop.nrec[1] > 1 && return (va = pop.G0[1, 1], vpe = pop.PE0[1, 1], ve = pop.R0[1, 1])
    ntraits(pop) == 1 && return (va = pop.G0[1, 1], ve = pop.R0[1, 1])
    rep = findall(>(1), pop.nrec)
    names = Tuple([[Symbol(p, t.name) for p in ("va_", "ve_") for t in pop.traits];
                   [Symbol("vpe_", pop.traits[j].name) for j in rep]])
    NamedTuple{names}(Tuple([diag(pop.G0); diag(pop.R0); [pop.PE0[j, j] for j in rep]]))
end

# mean coancestry of the founders, self included: they are unrelated and
# non-inbred, so it is 1/(2n₀)
θ0(pop::Population) = 1 / (2 * length(generation_ids(pop, 0)))

# The QTL of the aggregate genotype H = w'g: effects w'a at the union of
# the loci, intercept w'μ.
function _aggregate_qtl(pop::Population, w)
    loci = sort!(unique(reduce(vcat, [q.loci for q in pop.qtl])))
    a = zeros(length(loci))
    for (q, wj) in zip(pop.qtl, w), (l, x) in zip(q.loci, q.additive)
        a[searchsortedfirst(loci, l)] += wj * x
    end
    TraitQTL("H", loci, a; intercept = sum(w .* [q.intercept for q in pop.qtl]))
end

# The QTL of `q` with a non-zero effect, by the state of the favourable
# allele (the allele coded 1 if its effect is positive) in `hap`:
# segregating, fixed or lost. The three add up to the number of such QTL.
function _qtl_status(hap::Haplotype, q::TraitQTL)
    nz = findall(!iszero, q.additive)
    p = vec(sum(dosage(hap, q.loci[nz]); dims = 2)) ./ hap.nhp
    fav = [a > 0 ? x : 1 - x for (a, x) in zip(q.additive[nz], p)]
    return (qtl_seg = count(x -> 0 < x < 1, fav), qtl_fixed = count(==(1), fav), qtl_lost = count(==(0), fav),
        _qtl_architecture(p, q.additive[nz])...)
end

# The architecture of the QTL segregating in a cohort, from their allele
# frequencies `p` and effects `a`: the mean minor allele frequency; α, the
# slope of log a² on log 2p(1-p) (effect variance ∝ [2p(1-p)]^α; negative
# when large effects are rare, as under purifying or stabilizing
# selection; NaN with fewer than three QTL); and the share of the genic
# variance 2p(1-p)a² in the largest tenth of them.
function _qtl_architecture(p, a)
    seg = findall(x -> 0 < x < 1, p)
    isempty(seg) && return (qtl_maf = NaN, qtl_alpha = NaN, qtl_top10 = NaN)
    ps, as = p[seg], a[seg]
    h = 2 .* ps .* (1 .- ps)
    x, y = log.(h), log.(as .^ 2)
    α = length(seg) >= 3 && var(x) > 0 ? cov(x, y) / var(x) : NaN
    v = sort(h .* as .^ 2; rev = true)
    top = sum(v[1:ceil(Int, length(v) / 10)]) / sum(v)
    return (qtl_maf = mean(min.(ps, 1 .- ps)), qtl_alpha = α, qtl_top10 = top)
end

"""
    generation_summary(pop::Population, g; weights = ones(k), segs = nothing) -> NamedTuple

For cohort `g`: its size; mean, variance and genic variance
``\\sum 2p(1-p)a^2`` at the QTL of the genetic values (one trait), or of
the aggregate genotype ``H = w'g`` with `weights` as ``w`` (several
traits, then also `mean_<name>`, `var_<name>` per trait and the genetic
correlation `rg_<a>_<b>` of each pair within the cohort); the QTL with a
non-zero effect on that genetic value (on ``H``, the effects ``w'a``) at
which the favourable allele is segregating (`qtl_seg`), fixed
(`qtl_fixed`) or lost (`qtl_lost`) in the cohort, and the architecture
of the segregating ones: their mean minor allele frequency `qtl_maf`,
`qtl_alpha`, the slope of ``\\log a^2`` on ``\\log 2p(1-p)`` (effect
variance ``\\propto [2p(1-p)]^\\alpha``; negative when large effects are
rare), and `qtl_top10`, the share of the genic variance in the largest
tenth of them; mean pedigree
inbreeding `F` and mean pedigree `coancestry` ``\\bar\\theta_g`` over all
pairs of the cohort, self-coancestry ``(1 + F_i)/2`` included (founders
unrelated and non-inbred); the loss of heterozygosity `Fhet`; the
number of animals `alive` at step `g` (all cohorts); and, given the
founder segments `segs` (`founder_segments(pop, mm)`, which `step!`
passes when the scheme records meioses), the mean realized inbreeding
`Fseg` of the cohort, the proportion of the genome at which an animal's
two haplotypes descend from the same founder haplotype (its expectation
is the pedigree `F`).

`Fgrm` is the mean genomic inbreeding of the cohort, ``G_{ii} - 1`` of
VanRaden's (2008) first ``G = ZZ'/2\\sum p_0(1-p_0)`` with the founders'
allele frequencies ``p_0`` as base, over all loci segregating in the
founders: ``Z_{ij} = M_{ij} - 2p_{0j}`` for the dosage ``M_{ij}``. As
``p_0`` is estimated from the founders, two distinct founder alleles have
mean coancestry ``-1/(2n_0 - 1)`` on this base, so without selection the
expectation of `Fgrm` is about ``F - (1 - F)/(2n_0 - 1)``, ``F`` the
pedigree inbreeding (checked over 200 replicates of
`examples/minimal.json` with random selection, steps 2–6); in cohort 0 it
is the founders' own deviation from Hardy–Weinberg proportions. `Fgrm`
measures the realized homozygosity of each animal, `Fhet` the loss of
heterozygosity of the cohort.

`Fhet` is ``1 - (1 - \\bar\\theta_0)\\,\\bar H_g/\\bar H_0`` over all loci, with
``\\bar H`` the mean of ``2p(1-p)`` from the cohort's own allele
frequencies and ``\\bar\\theta_0 = 1/(2n_0)`` the founders' mean
coancestry. Each ``\\bar H_g`` is a sample estimate whose expectation is
``H(1 - \\bar\\theta_g)``, ``H`` the heterozygosity of the population the
founders were drawn from, so without selection the expectation of `Fhet`
is the mean coancestry ``\\bar\\theta_g``. Under selection,
`Fhet − coancestry` measures the extra loss of heterozygosity at loci
linked to QTL; both carry drift noise that is large when few chromosome
segments segregate.
"""
function generation_summary(pop::Population, g::Integer; weights = ones(ntraits(pop)), segs = nothing)
    ids = generation_ids(pop, g)
    all(i -> pop.slot[i] > 0, ids) || throw(ArgumentError(
        "genomes of cohort $g were dropped (retain = :active): its summary needs all of them"))
    hap = pop.blocks[g+1]
    F = nrm_diag(pop.sire, pop.dam)[ids] .- 1
    θ = mean(nrm(pop.sire, pop.dam, ids)) / 2             # mean coancestry, self included
    H(h) = (p = vec(sum(dosage(h, 1:h.nlc); dims = 2)) ./ h.nhp; mean(2 .* p .* (1 .- p)))
    h0 = pop.blocks[1]
    p0 = vec(sum(dosage(h0, 1:h0.nlc); dims = 2)) ./ h0.nhp
    seg = findall(x -> 0 < x < 1, p0)
    Z = dosage(hap, seg) .- 2 .* p0[seg]
    Fgrm = mean(sum(abs2, Z; dims = 1)) / (2 * sum(p -> p * (1 - p), p0[seg])) - 1
    k = ntraits(pop)
    q = k == 1 ? pop.qtl[1] : _aggregate_qtl(pop, weights)
    s = trait_summary(hap, q)
    out = (generation = Int(g), n = length(ids), mean = s.mean[1], var = s.var[1], genic = s.genic[1],
        _qtl_status(hap, q)...)
    if k > 1
        names = pop.traits
        V = pop.tbv[ids, :]
        C = cor(V)
        per = NamedTuple{Tuple([[Symbol("mean_", t.name) for t in names]; [Symbol("var_", t.name) for t in names]])}(
            Tuple([vec(mean(V; dims = 1)); vec(var(V; dims = 1))]))
        pairs = [(j, l) for j = 1:k for l = j+1:k]
        rg = NamedTuple{Tuple(Symbol("rg_", names[j].name, "_", names[l].name) for (j, l) in pairs)}(
            Tuple(C[j, l] for (j, l) in pairs))
        out = merge(out, per, rg)
    end
    merge(out, (F = mean(F), coancestry = θ, Fhet = 1 - (1 - θ0(pop)) * H(hap) / H(pop.blocks[1]), Fgrm = Fgrm,
        alive = count(i -> pop.generation[i] <= g < pop.exit[i], 1:length(pop))),
        segs === nothing ? NamedTuple() : (Fseg = mean(segment_inbreeding(segs; ids = ids)),))
end

"""
    simulate!(pop::Population, scheme::Scheme, mm::MeiosisMap, nstep; checkpoint = nothing) -> DataFrame

Run `nstep` time steps of `scheme` with [`step!`](@ref); one row per
cohort in which selection took place, plus the summary of the last cohort
(without accuracy, parents and generation interval). With `checkpoint` a
path, the population is saved there with [`save_checkpoint`](@ref) after
every step (atomically: written to a temporary file, then renamed).
"""
function simulate!(pop::Population, sc::Scheme, mm::MeiosisMap, nstep::Integer; checkpoint = nothing)
    rows = map(1:nstep) do _
        r = step!(pop, sc, mm)
        checkpoint === nothing || save_checkpoint(checkpoint, pop)
        r
    end
    return DataFrame([rows; _final_row(pop, sc, mm)])
end

# The summary of the last cohort, which closes the table of `simulate!`
# (no selection took place in it).
function _final_row(pop::Population, sc::Scheme, mm::MeiosisMap)
    segs = sc.record_meiosis ? founder_segments(pop, mm) : nothing
    return merge(generation_summary(pop, current_step(pop); weights = index_weights(pop, sc), segs = segs),
        _no_accuracy(pop), _no_reliability(pop), (nsire = 0, ndam = 0, L = NaN),
        (intensity_sire = NaN, intensity_dam = NaN),
        _vc(pop), _no_ocs(sc))
end

# ---- checkpoints -------------------------------------------------------

const _CHECKPOINT_MAGIC = UInt8['B', 'N', 'G', 'S', 'C', 'K', '1', 0]

"""
    save_checkpoint(path, pop::Population) -> path

Save the whole state of `pop` at a time-step boundary. A run continued
from [`load_checkpoint`](@ref) with the same scheme and `MeiosisMap` is
identical to an uninterrupted one, because every draw is keyed by
coordinate and nothing else carries state between steps.

The file is a magic, the versions of BnGsim and Julia, and the population
in Julia's `Serialization` format. That format is not guaranteed across
Julia or package versions, so a checkpoint is for resuming a run, not
for archiving results; write haplotypes with `BnGStructs.write_haplotypes`
for that.
"""
function save_checkpoint(path::AbstractString, pop::Population)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        write(io, _CHECKPOINT_MAGIC)
        serialize(io, (bngsim = string(pkgversion(@__MODULE__)), julia = string(VERSION)))
        serialize(io, pop)
    end
    mv(tmp, path; force = true)
    return path
end

"""
    load_checkpoint(path; strict = true) -> Population

Restore a population saved by [`save_checkpoint`](@ref). With `strict`, a
checkpoint written by another version of BnGsim or Julia is refused, since
the serialization format may differ; `strict = false` tries anyway.
"""
function load_checkpoint(path::AbstractString; strict::Bool = true)
    open(path) do io
        read(io, 8) == _CHECKPOINT_MAGIC || throw(ArgumentError("$path is not a BnGsim checkpoint"))
        v = deserialize(io)
        now = (bngsim = string(pkgversion(@__MODULE__)), julia = string(VERSION))
        strict && v != now && throw(ArgumentError(
            "$path was written by BnGsim $(v.bngsim) on Julia $(v.julia); this is BnGsim " *
            "$(now.bngsim) on Julia $(now.julia). Use strict = false to try anyway."))
        pop = deserialize(io)
        pop isa Population || throw(ArgumentError("$path does not hold a Population"))
        return pop
    end
end

# ---- realized IBD from founder segments --------------------------------

"""
    meioses(pop::Population) -> RelationshipMatrices.Meioses

The crossover records of all cohorts in the form of
`RelationshipMatrices.Meioses`: gamete `2k − 1` of the `k`-th animal of a
cohort is its haplotype from the sire, `2k` from the dam. Every cohort
must have been recorded (`Scheme(; record_meiosis = true)` from the
first step on).
"""
function meioses(pop::Population)
    any(isnothing, pop.meioses) &&
        throw(ArgumentError("some cohorts were bred without record_meiosis = true; their founder segments are unknown"))
    parent, hap, phases, ptr, pos = Int32[], Int32[], BitMatrix[], [1], UInt32[]
    for (t, r) in enumerate(pop.meioses)
        f = pop.first[t+1]
        for g in eachindex(r.parent)
            push!(parent, r.parent[g])
            push!(hap, 2(f + (g + 1) ÷ 2 - 1) - isodd(g))
            xo = crossovers(r, g)
            append!(pos, xo)
            push!(ptr, ptr[end] + length(xo))
        end
        push!(phases, r.phase)
    end
    phase = isempty(phases) ? falses(0, 0) : reduce(hcat, phases)
    Meioses(parent, hap, phase, ptr, pos)
end

"""
    founder_segments(pop::Population, mm::MeiosisMap) -> RelationshipMatrices.HaplotypeSegments

The founder haplotype of every base of every animal, reconstructed from
the recorded meioses (see [`meioses`](@ref)) and the pedigree, on the
chromosomes of `mm`. Founder haplotype `h` is haplotype `h` of cohort 0.
`RelationshipMatrices.irm(segs)` gives the realized IBD relationship
matrix and `segment_inbreeding(segs)` the realized inbreeding, without
sampling loci.
"""
function founder_segments(pop::Population, mm::MeiosisMap)
    nchr = length(mm.first)
    m = isempty(pop.meioses) ? Meioses(Int32[], Int32[], falses(nchr, 0), [1], UInt32[]) : meioses(pop)
    founder_segments(pop.sire, pop.dam, m, [mm.first[c]:mm.last[c] for c = 1:nchr])
end
