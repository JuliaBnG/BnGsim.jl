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

A trait's `sex` (`0` females, `1` males, `2` both) limits its records to
one sex; its `age`, a whole number of time steps, is the age at which a
record becomes available to evaluations (`0`: at birth). Phenotypes are
drawn for every animal and trait either way, from the animal's own keyed
stream, so these settings change which records are seen, not the values
drawn.
"""
function founder_population(hap::Haplotype, qtl::TraitQTL, trait::Trait, key::SimKey)
    _founders(hap, [qtl], Trait[trait], key, ones(1, 1))
end

function founder_population(hap::Haplotype, qtl::MultiTraitQTL, traits::AbstractVector{<:Trait}, key::SimKey;
    Re = I)
    k = length(traits)
    qtl.trait_names == [t.name for t in traits] ||
        throw(ArgumentError("the traits must be those of the QTL, in order: $(qtl.trait_names)"))
    isempty(qtl.dominance) || all(iszero, qtl.dominance) ||
        throw(ArgumentError("the breeding loop handles additive traits only"))
    Rm = Re isa UniformScaling ? Matrix{Float64}(Re(k)) : Matrix{Float64}(Re)
    size(Rm) == (k, k) && issymmetric(Rm) && all(Rm[j, j] == 1 for j = 1:k) && isposdef(Rm) ||
        throw(ArgumentError("Re must be a $k × $k positive definite correlation matrix"))
    qs = [TraitQTL(qtl.trait_names[j], qtl.loci, qtl.additive[:, j]; intercept = qtl.intercept[j]) for j = 1:k]
    _founders(hap, qs, Vector{Trait}(traits), key, Rm)
end

function _founders(hap, qs, traits, key, Re)
    any(_has_dominance, qs) && throw(ArgumentError("the breeding loop handles additive traits only"))
    all(t -> isinteger(t.age), traits) ||
        throw(ArgumentError("a trait's age must be a whole number of time steps"))
    allunique(t.name for t in traits) || throw(ArgumentError("trait names must be unique"))
    n = hap.nhp ÷ 2
    n >= 2 || throw(ArgumentError("need at least two founders"))
    G = reduce(hcat, [tbv(hap, q) for q in qs])
    σa = [t.σₐ for t in traits]
    σe = [t.σₐ * sqrt((1 - t.h²) / t.h²) for t in traits]
    C = length(traits) == 1 ? ones(1, 1) : cor(G)
    Lr = Diagonal(σe) * Matrix(cholesky(Symmetric(Re)).L)
    G0 = σa .* C .* σa'
    R0 = σe .* Re .* σe'
    pop = Population(zeros(Int32, n), zeros(Int32, n), falses(n), zeros(Int32, n), fill(typemax(Int32), n), G,
        fill(NaN, size(G)), fill(NaN, size(G)), Union{Nothing,Haplotype}[hap], [1], Int32.(1:n), qs, traits, Lr, key, G0, R0,
        Union{Nothing,MeiosisRecord}[])
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
  it lives);
- per trait `<name>`: `tbv_<name>`, the true genetic value;
  `phenotype_<name>`, the phenotype (`missing` in the sex that does not
  express the trait; it is drawn at birth, and is available to an
  evaluation only from the trait's `age` on, see [`observed`](@ref));
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
        for j = 1:k
            e = 0.0
            for l = 1:j
                e += pop.Lr[j, l] * z[l]
            end
            t = pop.traits[j]
            pop.phenotype[i, j] = _expresses(t, pop.male[i]) ? t.μ + pop.tbv[i, j] + e : NaN
        end
    end
end

"""
    Scheme(; noff, nsire, ndam, evaluation = :pblup, selection = :truncation,
           K = NaN, λ = 0.0, κ = 0.0, ub = 1.0, chip = Int32[], blend = 0.05,
           max_sire_age = 0, max_dam_age = 0, reml = false, weights = Float64[],
           survival = Float64[], cull_unselected = false, genotyped_cohorts = 0,
           ocs = :fixed, ΔF = NaN, record_meiosis = false)

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
- `record_meiosis`: keep the crossover records of every cohort in
  `pop.meioses`, for realized IBD from founder segments
  ([`founder_segments`](@ref)); about 8 bytes per crossover;
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
    function Scheme(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, chip, blend,
        max_sire_age, max_dam_age, reml, weights, survival, cull_unselected, genotyped_cohorts, ocs, ΔF,
        record_meiosis, window, ancestor_generations, retain)
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
            genotyped_cohorts, ocs, ΔF, record_meiosis, window, ancestor_generations, retain)
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
    if k == 1
        rec = findall(view(obs, :, 1))
        y = P[rec, 1]
        X = ones(length(rec), 1)
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
        y, X, Zs = _mt_design(P, obs)
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

function _converged(r, pop)
    r.iterations < 100 ||
        throw(ErrorException("REML did not converge in 100 iterates at time step $(current_step(pop))"))
end

# Multi-trait design over the records `obs` (n × k), stacked trait-major.
# Every trait keeps its genetic component, so a trait without records yet
# (e.g. recorded late in life) is predicted through G0; only the
# intercepts of traits with records enter X.
function _mt_design(P, obs)
    n, k = size(obs)
    ot = findall(j -> any(view(obs, :, j)), 1:k)
    y, _ = mt_stack(ifelse.(obs, P, NaN))
    X = reduce(hcat, mt_blocks(obs, [ones(n, 1) for _ = 1:k])[ot])
    Zs = mt_blocks(obs, [incidence(1:n, n) for _ = 1:k])
    return y, X, Zs
end

# ---- selection ---------------------------------------------------------

# Integer numbers of offspring for parents with contributions `w` (summing
# to 1 within the group): largest remainder, ties by ascending id.
function _allocate(ids, w, n)
    raw = w ./ sum(w) .* n
    k = floor.(Int, raw)
    rest = n - sum(k)
    order = sortperm(collect(zip(-(raw .- k), ids)))
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

# living males and females young enough to be parents at step g
function _eligible(pop::Population, sc::Scheme, g)
    young(a) = pop.first[max(0, g - a)+1]:length(pop)
    males = filter(i -> pop.male[i] && pop.exit[i] > g, young(sc.max_sire_age))
    females = filter(i -> !pop.male[i] && pop.exit[i] > g, young(sc.max_dam_age))
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
        "raise the maximum ages or the survival"))
    cand = sort!([males; females])
    if sc.selection in (:truncation, :random)
        crit = sc.selection == :truncation ? index :
               [rand(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.SELECTION,
                    stage = _STAGE_RANDOM_SELECTION, entity = i))) for i = 1:length(pop)]
        s = _best(males, crit[males], sc.nsire)
        d = _best(females, crit[females], sc.ndam)
        s = sort(s)
        d = sort(d)
        return s, _allocate(s, ones(length(s)), sc.noff), d, _allocate(d, ones(length(d)), sc.noff), NamedTuple()
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
    ns, nd = _allocate(s, c[ms], sc.noff), _allocate(d, c[fs], sc.noff)
    keep_s, keep_d = ns .> 0, nd .> 0
    return s[keep_s], ns[keep_s], d[keep_d], nd[keep_d], info
end

"""
    mating_plan(pop, sires, nsire_off, dams, ndam_off) -> Matrix{Int}

Pair sire and dam slots at random (a keyed shuffle of the dam slots,
`Domain.MATING`): row `i` is mating `i`, as ids.
"""
function mating_plan(pop::Population, s, ns, d, nd)
    sum(ns) == sum(nd) || throw(ArgumentError("sire and dam slots differ: $(sum(ns)) and $(sum(nd))"))
    sl = reduce(vcat, [fill(Int(s[k]), ns[k]) for k in eachindex(s)])
    dl = reduce(vcat, [fill(Int(d[k]), nd[k]) for k in eachindex(d)])
    g = current_step(pop)
    keyed_shuffle!(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.MATING, stage = _STAGE_PAIRING)), dl)
    return hcat(sl, dl)
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
over both sexes) and the variances used by the evaluation, `va` and `ve`
(one trait) or `va_<name>` and `ve_<name>`.
"""
function step!(pop::Population, sc::Scheme, mm::MeiosisMap)
    g = current_step(pop)
    w = index_weights(pop, sc)
    m = evaluate!(pop, sc)
    cand = generation_ids(pop, g)
    acc = _accuracies(pop, cand, w)
    s, ns, d, nd, info = select_parents(pop, sc, m)
    plan = mating_plan(pop, s, ns, d, nd)
    L = mean(g + 1 .- pop.generation[plan])
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
    pop.tbv = vcat(pop.tbv, reduce(hcat, [tbv(off, q) for q in pop.qtl]))
    pop.phenotype = vcat(pop.phenotype, fill(NaN, sc.noff, k))
    pop.ebv = vcat(pop.ebv, fill(NaN, sc.noff, k))
    append!(pop.exit, fill(typemax(Int32), sc.noff))
    push!(pop.blocks, off)
    push!(pop.first, n0 + 1)
    append!(pop.slot, Int32.(1:sc.noff))
    _sexes!(pop, ids, g + 1)
    _phenotypes!(pop, ids, g + 1)
    segs = sc.record_meiosis ? founder_segments(pop, mm) : nothing
    out = merge(generation_summary(pop, g; weights = w, segs = segs), acc,
        (nsire = length(s), ndam = length(d), L = L), _vc(pop), info)
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

function _accuracies(pop::Population, ids, w)
    k = ntraits(pop)
    acc = (accuracy = cor(pop.ebv[ids, :] * w, pop.tbv[ids, :] * w),)
    k == 1 && return acc
    names = Tuple(Symbol("accuracy_", t.name) for t in pop.traits)
    merge(acc, NamedTuple{names}(Tuple(cor(pop.ebv[ids, j], pop.tbv[ids, j]) for j = 1:k)))
end

_no_accuracy(pop::Population) = (k = ntraits(pop); k == 1 ? (accuracy = NaN,) :
    merge((accuracy = NaN,), NamedTuple{Tuple(Symbol("accuracy_", t.name) for t in pop.traits)}(Tuple(fill(NaN, k)))))

function _vc(pop::Population)
    ntraits(pop) == 1 && return (va = pop.G0[1, 1], ve = pop.R0[1, 1])
    names = Tuple(Symbol(p, t.name) for p in ("va_", "ve_") for t in pop.traits)
    NamedTuple{names}(Tuple([diag(pop.G0); diag(pop.R0)]))
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

"""
    generation_summary(pop::Population, g; weights = ones(k), segs = nothing) -> NamedTuple

For cohort `g`: its size; mean, variance and genic variance
``\\sum 2p(1-p)a^2`` at the QTL of the genetic values (one trait), or of
the aggregate genotype ``H = w'g`` with `weights` as ``w`` (several
traits, then also `mean_<name>`, `var_<name>` per trait and the genetic
correlation `rg_<a>_<b>` of each pair within the cohort); mean pedigree
inbreeding `F` and mean pedigree `coancestry` ``\\bar\\theta_g`` over all
pairs of the cohort, self-coancestry ``(1 + F_i)/2`` included (founders
unrelated and non-inbred); the loss of heterozygosity `Fhet`; the
number of animals `alive` at step `g` (all cohorts); and, given the
founder segments `segs` (`founder_segments(pop, mm)`, which `step!`
passes when the scheme records meioses), the mean realized inbreeding
`Fseg` of the cohort, the proportion of the genome at which an animal's
two haplotypes descend from the same founder haplotype (its expectation
is the pedigree `F`).

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
    k = ntraits(pop)
    s = trait_summary(hap, k == 1 ? pop.qtl[1] : _aggregate_qtl(pop, weights))
    out = (generation = Int(g), n = length(ids), mean = s.mean[1], var = s.var[1], genic = s.genic[1])
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
    merge(out, (F = mean(F), coancestry = θ, Fhet = 1 - (1 - θ0(pop)) * H(hap) / H(pop.blocks[1]),
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
        _no_accuracy(pop), (nsire = 0, ndam = 0, L = NaN), _vc(pop), _no_ocs(sc))
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
function RelationshipMatrices.founder_segments(pop::Population, mm::MeiosisMap)
    nchr = length(mm.first)
    m = isempty(pop.meioses) ? Meioses(Int32[], Int32[], falses(nchr, 0), [1], UInt32[]) : meioses(pop)
    founder_segments(pop.sire, pop.dam, m, [mm.first[c]:mm.last[c] for c = 1:nchr])
end
