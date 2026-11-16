# The breeding loop (roadmap Phase 4, BnGsim 0.4): discrete generations,
# one additive trait, evaluation through Breeding.jl, selection and keyed
# mating.
#
# Model of one generation:
# - every animal of the youngest generation is a candidate and has one
#   phenotype, y = μ + g + e, recorded at birth (both sexes);
# - all records and the full pedigree enter the evaluation;
# - sires and dams are chosen among the candidates and given numbers of
#   offspring; matings pair the sire and dam slots at random;
# - the offspring form the next generation (no overlap, no culling).
# Every random draw is keyed by a `SimKey` (domains MATING, PHENOTYPE,
# SELECTION), so a run is a function of its key and inputs.

# Stage codes within domains; fixed values, never schedule positions.
const _STAGE_SEX = 0x0001
const _STAGE_PAIRING = 0x0002
const _STAGE_RANDOM_SELECTION = 0x0003

"""
    Population

Animals of a discrete-generation breeding program, one row per animal in
order of birth. Animal ids are `1:n` and equal the row; nothing is ever
removed, so the equality holds (invariant I4 holds trivially until culling
is added).

Fields: `sire`, `dam` (`0` = unknown), `male`, `generation` (founders are
0), `tbv` (genetic value), `phenotype` (`NaN` without a record), `ebv`
(the last evaluation), and the haplotypes, one `Haplotype` per generation
in `blocks` with `first[g + 1]` the id of the first animal of generation
`g`. Build one with [`founder_population`](@ref) and extend it with
[`step!`](@ref).
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

Generation 0 of a breeding program: the individuals of `hap`, unrelated
founders, with sexes drawn at random (`Domain.MATING`), genetic values
from `qtl` and one phenotype each (`Domain.PHENOTYPE`). The residual
standard deviation is fixed from the trait as
``\\sigma_e = \\sigma_a\\sqrt{(1 - h^2)/h^2}``, with ``\\sigma_a`` the trait's
`σₐ`, so ``h^2`` is the heritability in the founders. `key` carries the
seed and replicate; its other fields are set per draw.
"""
function founder_population(hap::Haplotype, qtl::TraitQTL, trait::Trait, key::SimKey)
    _has_dominance(qtl) && throw(ArgumentError("the breeding loop handles additive traits only"))
    n = hap.nhp ÷ 2
    n >= 2 || throw(ArgumentError("need at least two founders"))
    σe = trait.σₐ * sqrt((1 - trait.h²) / trait.h²)
    pop = Population(zeros(Int32, n), zeros(Int32, n), falses(n), zeros(Int32, n),
        tbv(hap, qtl), fill(NaN, n), fill(NaN, n), [hap], [1], qtl, trait, σe, key)
    _sexes!(pop, 1:n, 0)
    _phenotypes!(pop, 1:n, 0)
    return pop
end

_has_dominance(q::TraitQTL) = any(!iszero, q.dominance)

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
        pop.phenotype[i] = pop.trait.μ + pop.tbv[i] + pop.σe * e
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
  random, a neutral reference).

Variance components of the evaluation are the true ones of the founders
(``\\sigma^2_a``, ``\\sigma^2_e``); they are not re-estimated.
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
    function Scheme(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, chip, blend)
        noff > 1 && nsire > 0 && ndam > 0 || throw(ArgumentError("noff > 1, nsire > 0 and ndam > 0 are required"))
        evaluation in (:tbv, :phenotype, :pblup, :gblup) ||
            throw(ArgumentError("evaluation must be :tbv, :phenotype, :pblup or :gblup"))
        selection in (:truncation, :ocs, :rocs, :random) ||
            throw(ArgumentError("selection must be :truncation, :ocs, :rocs or :random"))
        selection == :ocs && !(K > 0) && throw(ArgumentError(":ocs needs a constraint K > 0"))
        evaluation == :gblup && isempty(chip) && throw(ArgumentError(":gblup needs chip loci"))
        0 <= blend < 1 || throw(ArgumentError("blend must lie in [0, 1)"))
        0 < ub <= 1 || throw(ArgumentError("ub must lie in (0, 1]"))
        new(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, Int32.(chip), blend)
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

Fill `pop.ebv` for every animal by the scheme's evaluation. For `:pblup`
and `:gblup` the mixed-model equations ``y = 1\\mu + Za + e`` over all
records are built with Breeding.jl and solved directly; they are returned
for prediction error variances.
"""
function evaluate!(pop::Population, sc::Scheme)
    n = length(pop)
    if sc.evaluation == :tbv
        pop.ebv .= pop.tbv
        return nothing
    elseif sc.evaluation == :phenotype
        pop.ebv .= pop.phenotype
        return nothing
    end
    rec = findall(!isnan, pop.phenotype)
    y = pop.phenotype[rec]
    X = ones(length(rec), 1)
    Z = incidence(rec, n)
    σa2, σe2 = pop.trait.σₐ^2, pop.σe^2
    Kinv = if sc.evaluation == :pblup
        ainv(pop.sire, pop.dam)
    else
        H = _haplotypes(pop, 1:n)
        base = pop.blocks[1]                                  # the founders
        p = vec(sum(dosage(base, sc.chip); dims = 2)) ./ base.nhp
        G = grm(H; p = p, loci = sc.chip)
        sc.blend > 0 && (G = (1 - sc.blend) .* G .+ sc.blend .* nrm(pop.sire, pop.dam))
        inv(cholesky(Symmetric(G)))
    end
    m = mme(y, X, [RandomEffect(Z, Kinv, σa2)]; σ²e = σe2)
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

Parents among the youngest generation and their numbers of offspring
(each list summing to `scheme.noff`). `m` is the `MME` returned by
[`evaluate!`](@ref) (`nothing` for `:tbv` and `:phenotype`); `:rocs`
takes the candidates' prediction error covariance from it. Every ranking
breaks ties by ascending id (invariant I1a).
"""
function select_parents(pop::Population, sc::Scheme, m)
    g = ngenerations(pop) - 1
    cand = collect(generation_ids(pop, g))
    males = cand[pop.male[cand]]
    females = cand[.!pop.male[cand]]
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
    g = ngenerations(pop) - 1
    keyed_shuffle!(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.MATING, stage = _STAGE_PAIRING)), dl)
    return hcat(sl, dl)
end

# ---- one generation ----------------------------------------------------

"""
    step!(pop::Population, scheme::Scheme, mm::MeiosisMap) -> NamedTuple

One generation: evaluate, select, mate, and add the offspring (with
genetic values, sexes and phenotypes) as the next generation. Returns the
summary of the selected generation (see [`generation_summary`](@ref)) with
the accuracy and the parents used.
"""
function step!(pop::Population, sc::Scheme, mm::MeiosisMap)
    g = ngenerations(pop) - 1
    m = evaluate!(pop, sc)
    cand = generation_ids(pop, g)
    acc = cor(pop.ebv[cand], pop.tbv[cand])
    s, ns, d, nd = select_parents(pop, sc, m)
    plan = mating_plan(pop, s, ns, d, nd)
    parent_block = pop.blocks[end]
    off = Haplotype(parent_block.nlc, 2sc.noff)
    local_plan = plan .- (pop.first[end] - 1)
    make_progeny!(off, parent_block, local_plan, mm, SimKey(pop.key; generation = g + 1))
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
    return merge(generation_summary(pop, g), (accuracy = acc, nsire = length(s), ndam = length(d)))
end

"""
    generation_summary(pop::Population, g) -> NamedTuple

For generation `g`: its size, mean and variance of the genetic values,
genic variance ``\\sum 2p(1-p)a^2`` at the QTL, mean pedigree inbreeding
``F`` (founders unrelated and non-inbred) and the loss of heterozygosity
``1 - \\bar H_g / \\bar H_0`` over all loci, ``\\bar H`` the mean expected
heterozygosity ``2p(1-p)``; under drift alone its expectation is the
mean ``F`` of the gametes that form the generation.
"""
function generation_summary(pop::Population, g::Integer)
    ids = generation_ids(pop, g)
    hap = pop.blocks[g+1]
    F = nrm_diag(pop.sire, pop.dam)[ids] .- 1
    H(h) = (p = vec(sum(dosage(h, 1:h.nlc); dims = 2)) ./ h.nhp; mean(2 .* p .* (1 .- p)))
    s = trait_summary(hap, pop.qtl)
    (generation = Int(g), n = length(ids), mean = s.mean[1], var = s.var[1], genic = s.genic[1],
        F = mean(F), Fhet = 1 - H(hap) / H(pop.blocks[1]))
end

"""
    simulate!(pop::Population, scheme::Scheme, mm::MeiosisMap, ngen) -> DataFrame

Run `ngen` generations of `scheme` with [`step!`](@ref); one row per
selected generation (0 to `ngen − 1`), plus the summary of the last
generation (without accuracy and parents).
"""
function simulate!(pop::Population, sc::Scheme, mm::MeiosisMap, ngen::Integer)
    rows = [step!(pop, sc, mm) for _ = 1:ngen]
    last = merge(generation_summary(pop, ngenerations(pop) - 1), (accuracy = NaN, nsire = 0, ndam = 0))
    return DataFrame([rows; last])
end
