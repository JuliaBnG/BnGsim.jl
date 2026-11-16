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

Fields:

- `sire`, `dam` (`0` = unknown), `male`, `generation` (the birth cohort,
  i.e. the time step of birth; founders are 0; with discrete generations
  it is the generation number);
- `tbv`, `phenotype`, `ebv`: `n × k` matrices, one column per trait, of
  genetic values, phenotypes (`NaN` in the sex that does not express the
  trait) and the last evaluation;
- `blocks`, `first`: the haplotypes, one `Haplotype` per cohort, with
  `first[t + 1]` the id of the first animal of cohort `t`;
- `qtl` (one `TraitQTL` per trait), `traits`, `Lr` (lower-triangular
  factor of the true residual covariance matrix, ``R = L_rL_r'``), `key`;
- `G0`, `R0`: the genetic and residual covariance matrices of the next
  evaluation (the true ones of the founders unless re-estimated).

Build one with [`founder_population`](@ref), extend it with
[`step!`](@ref), and save or restore it with [`save_checkpoint`](@ref)
and [`load_checkpoint`](@ref).
"""
mutable struct Population
    sire::Vector{Int32}
    dam::Vector{Int32}
    male::BitVector
    generation::Vector{Int32}
    tbv::Matrix{Float64}
    phenotype::Matrix{Float64}
    ebv::Matrix{Float64}
    blocks::Vector{Haplotype}
    first::Vector{Int}
    qtl::Vector{TraitQTL}
    traits::Vector{Trait}
    Lr::Matrix{Float64}
    key::SimKey
    G0::Matrix{Float64}
    R0::Matrix{Float64}
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
    pop = Population(zeros(Int32, n), zeros(Int32, n), falses(n), zeros(Int32, n), G,
        fill(NaN, size(G)), fill(NaN, size(G)), [hap], [1], qs, traits, Lr, key, G0, R0)
    _sexes!(pop, 1:n, 0)
    _phenotypes!(pop, 1:n, 0)
    return pop
end

_has_dominance(q::TraitQTL) = any(!iszero, q.dominance)

# the current time step: the cohort of the youngest animals
current_step(pop::Population) = ngenerations(pop) - 1

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
sex) and is at least the trait's `age` old.
"""
function observed(pop::Population, t::Integer = current_step(pop))
    obs = falses(length(pop), ntraits(pop))
    for j = 1:ntraits(pop), i = 1:length(pop)
        obs[i, j] = !isnan(pop.phenotype[i, j]) && t - pop.generation[i] >= Int(pop.traits[j].age)
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
           max_sire_age = 0, max_dam_age = 0, reml = false, weights = Float64[])

A selection and mating scheme, applied every time step:

- `noff`: offspring per step (one per mating);
- `nsire`, `ndam`: parents selected by truncation, and the number of each
  sex used by random selection;
- `evaluation`: `:tbv` (true values, an oracle), `:phenotype` (own
  records), `:pblup` (animal model with the pedigree relationship matrix),
  `:gblup` (animal model with ``G^* = (1-w)G + wA``, `w = blend`, `G`
  VanRaden at the loci `chip` with base allele frequencies from the
  founders). With several traits the animal model is multi-trait, with
  an intercept per trait and residuals correlated within animals;
- `weights`: the index weights ``w`` (one per trait; empty: all `1`).
  Candidates are ranked on ``I = w'\\hat a``, which with economic values
  as weights and multi-trait BLUP EBVs is the optimal index for the
  aggregate genotype ``H = w'g`` (Henderson 1963). With `:phenotype` it is
  ``w'y``, a fixed-weight index, not the optimal phenotypic index;
- `selection`: `:truncation` (best `nsire` males and `ndam` females on
  the index, equal numbers of offspring), `:ocs` (`Breeding.TM1997` with
  ``c'Ac \\le K``, i.e. average coancestry ``K/2``), `:rocs`
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
  after one cohort).
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
    function Scheme(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, chip, blend,
        max_sire_age, max_dam_age, reml, weights)
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
        all(isfinite, weights) && (isempty(weights) || any(!iszero, weights)) ||
            throw(ArgumentError("weights must be finite and not all zero"))
        new(noff, nsire, ndam, evaluation, selection, K, λ, κ, ub, Int32.(chip), blend,
            max_sire_age, max_dam_age, reml, Float64.(weights))
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

function _kinv(pop::Population, sc::Scheme)
    sc.evaluation == :pblup && return ainv(pop.sire, pop.dam)
    n = length(pop)
    H = _haplotypes(pop, 1:n)
    base = pop.blocks[1]                                  # the founders
    p = zeros(base.nlc)                                   # grm wants one per locus
    p[sc.chip] = vec(sum(dosage(base, sc.chip); dims = 2)) ./ base.nhp
    G = grm(H; p = p, loci = sc.chip)
    sc.blend > 0 && (G = (1 - sc.blend) .* G .+ sc.blend .* nrm(pop.sire, pop.dam))
    inv(cholesky(Symmetric(G)))
end

"""
    evaluate!(pop::Population, scheme::Scheme) -> Union{Nothing, MME}

Fill `pop.ebv` for every animal and trait by the scheme's evaluation at
the current time step. For `:pblup` and `:gblup` the mixed-model
equations over the records available (see [`observed`](@ref)) are built
with Breeding.jl and solved directly: ``y = 1\\mu + Za + e`` with
``\\sigma^2_e`` scaling for one trait; for several traits, records stacked
trait-major with an intercept per recorded trait,
``\\operatorname{var}(a) = G_0 \\otimes K`` and ``R_0`` among the records
of one animal (missing records handled exactly). The covariance
components are `pop.G0` and `pop.R0`, first re-estimated by REML if
`scheme.reml`. The equations are returned for prediction error variances.
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
    obs = observed(pop)
    any(obs) || throw(ArgumentError("no records are available at time step $(current_step(pop))"))
    Kinv = _kinv(pop, sc)
    # Genetic and residual components are not separable from records of
    # unrelated animals only (the founders at step 0): REML waits for every
    # trait to be recorded on an animal with a known parent.
    known = (pop.sire .> 0) .| (pop.dam .> 0)
    do_reml = sc.reml && all(j -> any(view(obs, :, j) .& known), 1:k)
    if k == 1
        rec = findall(view(obs, :, 1))
        y = pop.phenotype[rec, 1]
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
        y, X, Zs = _mt_design(pop, obs)
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
    pop.ebv .= reshape(x[a:a+k*n-1], n, k)
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
function _mt_design(pop::Population, obs)
    n, k = size(obs)
    ot = findall(j -> any(view(obs, :, j)), 1:k)
    y, _ = mt_stack(ifelse.(obs, pop.phenotype, NaN))
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

"""
    select_parents(pop::Population, scheme::Scheme, m) -> (sires, nsire_off, dams, ndam_off)

Parents and their numbers of offspring (each list summing to
`scheme.noff`), among the males of age at most `scheme.max_sire_age` and
the females of age at most `scheme.max_dam_age` at the current time step
(the youngest cohort only, with the default ages 0), ranked on the index
``w'\\hat a`` (see [`Scheme`](@ref)). `m` is the `MME` returned by
[`evaluate!`](@ref) (`nothing` for `:tbv` and `:phenotype`); `:rocs` takes
the prediction error covariance of the candidates' index from it. Every
ranking breaks ties by ascending id (invariant I1a).
"""
function select_parents(pop::Population, sc::Scheme, m)
    g = current_step(pop)
    w = index_weights(pop, sc)
    index = pop.ebv * w
    eligible(a) = pop.first[max(0, g - a)+1]:length(pop)
    males = filter(i -> pop.male[i], eligible(sc.max_sire_age))
    females = filter(i -> !pop.male[i], eligible(sc.max_dam_age))
    cand = sort!([males; females])
    if sc.selection in (:truncation, :random)
        crit = sc.selection == :truncation ? index :
               [rand(keyed_rng(SimKey(pop.key; generation = g, domain = Domain.SELECTION,
                    stage = _STAGE_RANDOM_SELECTION, entity = i))) for i = 1:length(pop)]
        s = _best(males, crit[males], sc.nsire)
        d = _best(females, crit[females], sc.ndam)
        s = sort(s)
        d = sort(d)
        return s, _allocate(s, ones(length(s)), sc.noff), d, _allocate(d, ones(length(d)), sc.noff)
    end
    A = nrm(pop.sire, pop.dam, cand)
    u = index[cand]
    sex = Int.(pop.male[cand])
    c = if sc.selection == :ocs
        2 .* TM1997(u, sex, A, sc.K)
    else
        Ω = m === nothing ? zeros(length(cand), length(cand)) : _pev_index(m, cand, w, length(pop))
        2 .* rocs(u, A; Ω = Ω, sex = sex, λ = sc.λ, κ = sc.κ, ub = sc.ub).w
    end
    ms = findall(k -> pop.male[cand[k]] && c[k] > 1e-9, eachindex(cand))
    fs = findall(k -> !pop.male[cand[k]] && c[k] > 1e-9, eachindex(cand))
    s, d = cand[ms], cand[fs]
    ns, nd = _allocate(s, c[ms], sc.noff), _allocate(d, c[fs], sc.noff)
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
    k = ntraits(pop)
    append!(pop.sire, Int32.(plan[:, 1]))
    append!(pop.dam, Int32.(plan[:, 2]))
    append!(pop.male, falses(sc.noff))
    append!(pop.generation, fill(Int32(g + 1), sc.noff))
    pop.tbv = vcat(pop.tbv, reduce(hcat, [tbv(off, q) for q in pop.qtl]))
    pop.phenotype = vcat(pop.phenotype, fill(NaN, sc.noff, k))
    pop.ebv = vcat(pop.ebv, fill(NaN, sc.noff, k))
    push!(pop.blocks, off)
    push!(pop.first, n0 + 1)
    _sexes!(pop, ids, g + 1)
    _phenotypes!(pop, ids, g + 1)
    return merge(generation_summary(pop, g; weights = w), acc,
        (nsire = length(s), ndam = length(d), L = L), _vc(pop))
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
    generation_summary(pop::Population, g; weights = ones(k)) -> NamedTuple

For cohort `g`: its size; mean, variance and genic variance
``\\sum 2p(1-p)a^2`` at the QTL of the genetic values (one trait), or of
the aggregate genotype ``H = w'g`` with `weights` as ``w`` (several
traits, then also `mean_<name>`, `var_<name>` per trait and the genetic
correlation `rg_<a>_<b>` of each pair within the cohort); mean pedigree
inbreeding `F` and mean pedigree `coancestry` ``\\bar\\theta_g`` over all
pairs of the cohort, self-coancestry ``(1 + F_i)/2`` included (founders
unrelated and non-inbred); and the loss of heterozygosity `Fhet`.

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
function generation_summary(pop::Population, g::Integer; weights = ones(ntraits(pop)))
    ids = generation_ids(pop, g)
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
    merge(out, (F = mean(F), coancestry = θ, Fhet = 1 - (1 - θ0(pop)) * H(hap) / H(pop.blocks[1])))
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
    last = merge(generation_summary(pop, current_step(pop); weights = index_weights(pop, sc)),
        _no_accuracy(pop), (nsire = 0, ndam = 0, L = NaN), _vc(pop))
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
