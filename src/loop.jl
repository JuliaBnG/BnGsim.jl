# The breeding loop (roadmap Phase 4, BnGsim 0.4–0.5): one additive trait,
# evaluation through Breeding.jl, selection and keyed mating, in discrete
# or overlapping generations.
#
# Model of one time step t (one birth cohort per step):
# - every animal has one phenotype, y = μ + g + e, drawn at birth; it is
#   recorded only in the sex that expresses the trait (`trait.sex`) and
#   enters evaluations from age `trait.age` (in time steps) on;
# - all records available at t and the full pedigree enter the evaluation,
#   with the true or, optionally, REML-re-estimated variance components;
# - sires and dams are chosen among the males and females of age at most
#   `max_sire_age` and `max_dam_age` (age 0: born at t) and given numbers
#   of offspring; matings pair the sire and dam slots at random;
# - the offspring form cohort t + 1. Nothing is removed: an animal not
#   selected at t can be selected at t + 1 while it is young enough.
# Every random draw is keyed by a `SimKey` (domains MATING, PHENOTYPE,
# SELECTION), so a run is a function of its key and inputs.

# Stage codes within domains; fixed values, never schedule positions.
const _STAGE_SEX = 0x0001
const _STAGE_PAIRING = 0x0002
const _STAGE_RANDOM_SELECTION = 0x0003

"""
    Population

Animals of a breeding program, one row per animal in order of birth.
Animal ids are `1:n` and equal the row; nothing is ever removed, so the
equality holds (invariant I4 holds trivially).

Fields: `sire`, `dam` (`0` = unknown), `male`, `generation` (the birth
cohort, i.e. the time step of birth; founders are 0; with discrete
generations it is the generation number), `tbv` (genetic value),
`phenotype` (`NaN` in the sex that does not express the trait), `ebv`
(the last evaluation), the haplotypes, one `Haplotype` per cohort in
`blocks` with `first[t + 1]` the id of the first animal of cohort `t`,
and `vc = [σ²a, σ²e]`, the variance components of the next evaluation
(the true ones of the founders unless re-estimated). Build one with
[`founder_population`](@ref), extend it with [`step!`](@ref), and save
or restore it with [`save_checkpoint`](@ref) and
[`load_checkpoint`](@ref).
"""
mutable struct Population
    sire::Vector{Int32}
    dam::Vector{Int32}
    male::BitVector
    generation::Vector{Int32}
    tbv::Vector{Float64}
    phenotype::Vector{Float64}
    ebv::Vector{Float64}
    blocks::Vector{Haplotype}
    first::Vector{Int}
    qtl::TraitQTL
    trait::Trait
    σe::Float64
    key::SimKey
    vc::Vector{Float64}
end

Base.length(pop::Population) = length(pop.sire)
ngenerations(pop::Population) = length(pop.blocks)
generation_ids(pop::Population, g::Integer) = pop.first[g+1]:(g + 2 <= length(pop.first) ? pop.first[g+2] - 1 : length(pop))

function Base.show(io::IO, pop::Population)
    print(io, "Population of ", length(pop), " animals in ", ngenerations(pop), " generations (trait \"",
        pop.trait.name, "\", h² = ", pop.trait.h², ")")
end

"""
    founder_population(hap::Haplotype, qtl::TraitQTL, trait::Trait, key::SimKey) -> Population

Cohort 0 of a breeding program: the individuals of `hap`, unrelated
founders, with sexes drawn at random (`Domain.MATING`), genetic values
from `qtl` and one phenotype each (`Domain.PHENOTYPE`). The residual
standard deviation is fixed from the trait as
``\\sigma_e = \\sigma_a\\sqrt{(1 - h^2)/h^2}``, with ``\\sigma_a`` the trait's
`σₐ`, so ``h^2`` is the heritability in the founders. `key` carries the
seed and replicate; its other fields are set per draw.

The trait's `sex` (`0` females, `1` males, `2` both) limits the records to
one sex; its `age`, a whole number of time steps, is the age at which a
record becomes available to evaluations (`0`: at birth). A phenotype is
drawn for every animal either way, from its own keyed stream, so these
settings change which records are seen, not the values drawn.
"""
function founder_population(hap::Haplotype, qtl::TraitQTL, trait::Trait, key::SimKey)
    _has_dominance(qtl) && throw(ArgumentError("the breeding loop handles additive traits only"))
    isinteger(trait.age) || throw(ArgumentError("the trait's age must be a whole number of time steps"))
    n = hap.nhp ÷ 2
    n >= 2 || throw(ArgumentError("need at least two founders"))
    σe = trait.σₐ * sqrt((1 - trait.h²) / trait.h²)
    pop = Population(zeros(Int32, n), zeros(Int32, n), falses(n), zeros(Int32, n),
        tbv(hap, qtl), fill(NaN, n), fill(NaN, n), [hap], [1], qtl, trait, σe, key,
        [trait.σₐ^2, σe^2])
    _sexes!(pop, 1:n, 0)
    _phenotypes!(pop, 1:n, 0)
    return pop
end

_has_dominance(q::TraitQTL) = any(!iszero, q.dominance)

# the current time step: the cohort of the youngest animals
current_step(pop::Population) = ngenerations(pop) - 1

_expresses(pop::Population, i) = pop.trait.sex == 2 || (pop.trait.sex == 1) == pop.male[i]

"""
    recorded(pop::Population, t = current_step(pop)) -> Vector{Int}

Ids of the animals whose phenotype is available to an evaluation at time
step `t`: expressed in their sex and of age at least `trait.age`.
"""
function recorded(pop::Population, t::Integer = current_step(pop))
    a = Int(pop.trait.age)
    filter(i -> !isnan(pop.phenotype[i]) && t - pop.generation[i] >= a, 1:length(pop))
end

function _sexes!(pop::Population, ids, g)
    for i in ids
        pop.male[i] = rand(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.MATING,
            stage = _STAGE_SEX, entity = i)), Bool)
    end
    any(pop.male[ids]) && !all(pop.male[ids]) ||
        throw(ArgumentError("generation $g has one sex only; increase its size"))
end

function _phenotypes!(pop::Population, ids, g)
    for i in ids
        e = randn(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.PHENOTYPE, entity = i)))
        pop.phenotype[i] = _expresses(pop, i) ? pop.trait.μ + pop.tbv[i] + pop.σe * e : NaN
    end
end

"""
    Scheme(; noff, nsire, ndam, evaluation = :pblup, selection = :truncation,
           K = NaN, λ = 0.0, κ = 0.0, ub = 1.0, chip = Int32[], blend = 0.05)

A selection and mating scheme, applied every generation:

- `noff`: offspring per generation (one per mating);
- `nsire`, `ndam`: parents selected by truncation, and the number of each
  sex used by random selection;
- `evaluation`: `:tbv` (true values, an oracle), `:phenotype` (own
  record), `:pblup` (animal model with the pedigree relationship matrix),
  `:gblup` (animal model with ``G^* = (1-w)G + wA``, `w = blend`, `G`
  VanRaden at the loci `chip` with base allele frequencies from the
  founders);
- `selection`: `:truncation` (best `nsire` males and `ndam` females, equal
  numbers of offspring), `:ocs` (`Breeding.TM1997` with ``c'Ac \\le K``,
  i.e. average coancestry ``K/2``), `:rocs` (`Breeding.rocs` with `λ`, `κ`,
  the PEV of the candidates as `Ω`, and an upper bound `ub` on each
  parent's contribution, in `rocs` units where each sex contributes ½;
  `ub = 0.5/nsire` with `λ = κ = 0` is truncation selection), or `:random` (`nsire` and `ndam` chosen at
  random, a neutral reference);
- `max_sire_age`, `max_dam_age`: the oldest age, in time steps, at which an
  animal can be selected as a parent (age 0: born in the current step).
  `0` for both gives discrete generations; larger values give overlapping
  generations, with all eligible ages competing on one criterion. With
  `:ocs` and `:rocs` the contributions are those of the candidates to the
  next cohort, without the age-class weighting of Meuwissen & Sonesson
  (1998) or Bijma & Woolliams (2000); the constraint then does not bound
  the long-term rate of inbreeding exactly;
- `reml`: if `true`, the variance components of `:pblup` and `:gblup` are
  re-estimated by AI-REML (Breeding.jl) from all records before every
  evaluation, starting from the previous estimates, once some recorded
  animal has a known parent (before that, as with the founders' records
  alone, ``\\sigma^2_a`` and ``\\sigma^2_e`` are not separable and the
  true values are used); otherwise the true
  values of the founders, ``\\sigma^2_a`` and ``\\sigma^2_e``, are used.
  With all records and the full pedigree, REML estimates the base-population
  variance even under selection (Sorensen & Kennedy 1984).
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
    function Scheme(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, chip, blend,
        max_sire_age, max_dam_age, reml)
        noff > 1 && nsire > 0 && ndam > 0 || throw(ArgumentError("noff > 1, nsire > 0 and ndam > 0 are required"))
        max_sire_age >= 0 && max_dam_age >= 0 || throw(ArgumentError("maximum ages must be non-negative"))
        reml && !(evaluation in (:pblup, :gblup)) &&
            throw(ArgumentError("reml = true needs evaluation :pblup or :gblup"))
        evaluation in (:tbv, :phenotype, :pblup, :gblup) ||
            throw(ArgumentError("evaluation must be :tbv, :phenotype, :pblup or :gblup"))
        selection in (:truncation, :ocs, :rocs, :random) ||
            throw(ArgumentError("selection must be :truncation, :ocs, :rocs or :random"))
        selection == :ocs && !(K > 0) && throw(ArgumentError(":ocs needs a constraint K > 0"))
        evaluation == :gblup && isempty(chip) && throw(ArgumentError(":gblup needs chip loci"))
        0 <= blend < 1 || throw(ArgumentError("blend must lie in [0, 1)"))
        0 < ub <= 1 || throw(ArgumentError("ub must lie in (0, 1]"))
        new(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, Int32.(chip), blend,
            max_sire_age, max_dam_age, reml)
    end
end

# ---- evaluation --------------------------------------------------------

function _haplotypes(pop::Population, ids)
    nlc = pop.blocks[1].nlc
    out = Haplotype(nlc, 2length(ids))
    for (k, i) in enumerate(ids)
        g = pop.generation[i] + 1
        j = i - pop.first[g] + 1
        for h = 0:1
            copyto!(packed(out), (2k - 2 + h) * size(packed(out), 1) + 1, packed(pop.blocks[g]),
                (2j - 2 + h) * size(packed(pop.blocks[g]), 1) + 1, size(packed(out), 1))
        end
    end
    return out
end

"""
    evaluate!(pop::Population, scheme::Scheme) -> Union{Nothing, MME}

Fill `pop.ebv` for every animal by the scheme's evaluation at the current
time step. For `:pblup` and `:gblup` the mixed-model equations
``y = 1\\mu + Za + e`` over the records available (see [`recorded`](@ref))
are built with Breeding.jl and solved directly, with the variance
components `pop.vc`, first re-estimated by REML if `scheme.reml`; the
equations are returned for prediction error variances. `:phenotype` needs
a trait recorded at birth in both sexes.
"""
function evaluate!(pop::Population, sc::Scheme)
    n = length(pop)
    if sc.evaluation == :tbv
        pop.ebv .= pop.tbv
        return nothing
    elseif sc.evaluation == :phenotype
        pop.trait.sex == 2 && pop.trait.age == 0 ||
            throw(ArgumentError(":phenotype needs a trait recorded at birth in both sexes"))
        pop.ebv .= pop.phenotype
        return nothing
    end
    rec = recorded(pop)
    isempty(rec) && throw(ArgumentError("no records are available at time step $(current_step(pop))"))
    y = pop.phenotype[rec]
    X = ones(length(rec), 1)
    Z = incidence(rec, n)
    Kinv = if sc.evaluation == :pblup
        ainv(pop.sire, pop.dam)
    else
        H = _haplotypes(pop, 1:n)
        base = pop.blocks[1]                                  # the founders
        p = zeros(base.nlc)                                   # grm wants one per locus
        p[sc.chip] = vec(sum(dosage(base, sc.chip); dims = 2)) ./ base.nhp
        G = grm(H; p = p, loci = sc.chip)
        sc.blend > 0 && (G = (1 - sc.blend) .* G .+ sc.blend .* nrm(pop.sire, pop.dam))
        inv(cholesky(Symmetric(G)))
    end
    # σ²a and σ²e are not separable from records of unrelated animals only
    # (the founders at step 0): REML waits for a recorded animal with a
    # known parent
    if sc.reml && any(i -> pop.sire[i] > 0 || pop.dam[i] > 0, rec)
        r = reml(y, X, [RandomEffect(Z, Kinv, pop.vc[1])]; σ²e = pop.vc[2], tol = 1e-6, maxiter = 100)
        r.iterations < 100 ||
            throw(ErrorException("REML did not converge in 100 iterates at time step $(current_step(pop))"))
        pop.vc .= (r.G0[1][1, 1], r.σ²e)
    end
    m = mme(y, X, [RandomEffect(Z, Kinv, pop.vc[1])]; σ²e = pop.vc[2])
    x = solve_mme(m)
    pop.ebv .= x[2:end]
    return m
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

"""
    select_parents(pop::Population, scheme::Scheme, m) -> (sires, nsire_off, dams, ndam_off)

Parents and their numbers of offspring (each list summing to
`scheme.noff`), among the males of age at most `scheme.max_sire_age` and
the females of age at most `scheme.max_dam_age` at the current time step
(the youngest cohort only, with the default ages 0). `m` is the `MME`
returned by [`evaluate!`](@ref) (`nothing` for `:tbv` and `:phenotype`);
`:rocs` takes the candidates' prediction error covariance from it. Every
ranking breaks ties by ascending id (invariant I1a).
"""
function select_parents(pop::Population, sc::Scheme, m)
    g = current_step(pop)
    eligible(a) = pop.first[max(0, g - a)+1]:length(pop)
    males = filter(i -> pop.male[i], eligible(sc.max_sire_age))
    females = filter(i -> !pop.male[i], eligible(sc.max_dam_age))
    cand = sort!([males; females])
    if sc.selection in (:truncation, :random)
        crit = sc.selection == :truncation ? pop.ebv :
               [rand(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.SELECTION,
                    stage = _STAGE_RANDOM_SELECTION, entity = i))) for i = 1:length(pop)]
        s = _best(males, crit[males], sc.nsire)
        d = _best(females, crit[females], sc.ndam)
        s = sort(s)
        d = sort(d)
        return s, _allocate(s, ones(length(s)), sc.noff), d, _allocate(d, ones(length(d)), sc.noff)
    end
    A = nrm(pop.sire, pop.dam, cand)
    u = pop.ebv[cand]
    sex = Int.(pop.male[cand])
    w = if sc.selection == :ocs
        2 .* TM1997(u, sex, A, sc.K)
    else
        Ω = m === nothing ? zeros(length(cand)) : pev_candidates(m, 1, cand)
        2 .* rocs(u, A; Ω = Ω, sex = sex, λ = sc.λ, κ = sc.κ, ub = sc.ub).w
    end
    ms = findall(k -> pop.male[cand[k]] && w[k] > 1e-9, eachindex(cand))
    fs = findall(k -> !pop.male[cand[k]] && w[k] > 1e-9, eachindex(cand))
    s, d = cand[ms], cand[fs]
    ns, nd = _allocate(s, w[ms], sc.noff), _allocate(d, w[fs], sc.noff)
    keep_s, keep_d = ns .> 0, nd .> 0
    return s[keep_s], ns[keep_s], d[keep_d], nd[keep_d]
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
accuracy of the evaluation in it, the parents used, the generation
interval `L` (mean age of the parents at the birth of their offspring,
over both sexes) and the variance components used, `va` and `ve`.
"""
function step!(pop::Population, sc::Scheme, mm::MeiosisMap)
    g = current_step(pop)
    m = evaluate!(pop, sc)
    cand = generation_ids(pop, g)
    acc = cor(pop.ebv[cand], pop.tbv[cand])
    s, ns, d, nd = select_parents(pop, sc, m)
    plan = mating_plan(pop, s, ns, d, nd)
    L = mean(g + 1 .- pop.generation[plan])
    # the parents' haplotypes, gathered from their cohorts; offspring
    # depend on the parents' alleles and the keys only, not on this order
    ups = sort!(unique(plan))
    parents = _haplotypes(pop, ups)
    off = Haplotype(parents.nlc, 2sc.noff)
    local_plan = map(i -> searchsortedfirst(ups, i), plan)
    make_progeny!(off, parents, local_plan, mm, SimKey(pop.key; generation = g + 1))
    n0 = length(pop)
    ids = (n0+1):(n0+sc.noff)
    append!(pop.sire, Int32.(plan[:, 1]))
    append!(pop.dam, Int32.(plan[:, 2]))
    append!(pop.male, falses(sc.noff))
    append!(pop.generation, fill(Int32(g + 1), sc.noff))
    append!(pop.tbv, tbv(off, pop.qtl))
    append!(pop.phenotype, fill(NaN, sc.noff))
    append!(pop.ebv, fill(NaN, sc.noff))
    push!(pop.blocks, off)
    push!(pop.first, n0 + 1)
    _sexes!(pop, ids, g + 1)
    _phenotypes!(pop, ids, g + 1)
    return merge(generation_summary(pop, g), (accuracy = acc, nsire = length(s), ndam = length(d),
        L = L, va = pop.vc[1], ve = pop.vc[2]))
end

# mean coancestry of the founders, self included: they are unrelated and
# non-inbred, so it is 1/(2n₀)
θ0(pop::Population) = 1 / (2 * length(generation_ids(pop, 0)))

"""
    generation_summary(pop::Population, g) -> NamedTuple

For generation `g`: its size, mean and variance of the genetic values,
genic variance ``\\sum 2p(1-p)a^2`` at the QTL, mean pedigree inbreeding
`F` and mean pedigree `coancestry` ``\\bar\\theta_g`` over all pairs of
the generation, self-coancestry ``(1 + F_i)/2`` included (founders
unrelated and non-inbred), and the loss of heterozygosity `Fhet`.

`Fhet` is ``1 - (1 - \\bar\\theta_0)\\,\\bar H_g/\\bar H_0`` over all loci, with
``\\bar H`` the mean of ``2p(1-p)`` from the generation's own allele
frequencies and ``\\bar\\theta_0 = 1/(2n_0)`` the founders' mean
coancestry. Each ``\\bar H_g`` is a sample estimate whose expectation is
``H(1 - \\bar\\theta_g)``, ``H`` the heterozygosity of the population the
founders were drawn from, so without selection the expectation of `Fhet`
is the mean coancestry ``\\bar\\theta_g``. Under selection,
`Fhet − coancestry` measures the extra loss of heterozygosity at loci
linked to QTL; both carry drift noise that is large when few chromosome
segments segregate.
"""
function generation_summary(pop::Population, g::Integer)
    ids = generation_ids(pop, g)
    hap = pop.blocks[g+1]
    F = nrm_diag(pop.sire, pop.dam)[ids] .- 1
    θ = mean(nrm(pop.sire, pop.dam, ids)) / 2             # mean coancestry, self included
    H(h) = (p = vec(sum(dosage(h, 1:h.nlc); dims = 2)) ./ h.nhp; mean(2 .* p .* (1 .- p)))
    s = trait_summary(hap, pop.qtl)
    (generation = Int(g), n = length(ids), mean = s.mean[1], var = s.var[1], genic = s.genic[1],
        F = mean(F), coancestry = θ, Fhet = 1 - (1 - θ0(pop)) * H(hap) / H(pop.blocks[1]))
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
    last = merge(generation_summary(pop, current_step(pop)),
        (accuracy = NaN, nsire = 0, ndam = 0, L = NaN, va = pop.vc[1], ve = pop.vc[2]))
    return DataFrame([rows; last])
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
