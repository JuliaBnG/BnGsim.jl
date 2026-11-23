"""
    RandomSign(d)

Effects drawn from the distribution `d` with a random sign: `x` from `d`,
then `-x` with probability ½ (one draw from `d`, then one `Bool`, per
effect). With a Gamma `d`, the effect distribution of Meuwissen, Hayes and
Goddard (2001): symmetric, with the heavy tails of the Gamma. Used as a
`Trait`'s `da` (scenario family `"gamma"` with `"sign": "random"`). Only
`rand`, `mean` (0) and `var` (``E[x^2]`` of `d`) are defined.
"""
struct RandomSign{D<:UnivariateDistribution} <: ContinuousUnivariateDistribution
    d::D
end
Statistics.mean(::RandomSign) = 0.0
Statistics.var(s::RandomSign) = var(s.d) + mean(s.d)^2
function Base.rand(rng::AbstractRNG, s::RandomSign)
    x = rand(rng, s.d)
    return rand(rng, Bool) ? -x : x
end

"""
    SpikeSlab(π, sd; spike_sd = 0.0)

A spike-and-slab effect distribution: with probability `π` an effect from
the slab `Normal(0, sd)`, otherwise from the spike, exactly `0` (or
`Normal(0, spike_sd)` with `spike_sd > 0`). One uniform, then one normal
when the effect is not exactly zero, per effect. With a point-mass spike,
about a fraction `1 - π` of the QTL have no effect (scenario family
`"spike_slab"`). Only `rand`, `mean` (0) and `var` are defined.
"""
struct SpikeSlab <: ContinuousUnivariateDistribution
    π::Float64
    sd::Float64
    spike_sd::Float64
    function SpikeSlab(π::Real, sd::Real; spike_sd::Real = 0.0)
        0 < π <= 1 || throw(ArgumentError("π must lie in (0, 1]"))
        sd > 0 && spike_sd >= 0 || throw(ArgumentError("sd must be positive and spike_sd non-negative"))
        new(π, sd, spike_sd)
    end
end
Statistics.mean(::SpikeSlab) = 0.0
Statistics.var(s::SpikeSlab) = s.π * s.sd^2 + (1 - s.π) * s.spike_sd^2
function Base.rand(rng::AbstractRNG, s::SpikeSlab)
    rand(rng) < s.π && return s.sd * randn(rng)
    return s.spike_sd > 0 ? s.spike_sd * randn(rng) : 0.0
end

# Extract nqtl × nid Int8 dosage matrix at `loci` from a Haplotype
_qtl_dosage_matrix(hap::Haplotype, loci::AbstractVector{<:Integer}) = dosage(hap, loci)

"""
    eqtl(hap::Haplotype, loci, t::AbstractTrait; rng = Random.default_rng(), atol = 1e-8) -> TraitQTL

Sample additive QTL effects of trait `t` at `loci` (a `LocusSet` or locus
indices) and scale them so that the genetic values of the individuals in
`hap`, the base population, have standard deviation `t.σₐ` and mean 0.

- Effects are drawn from `t.da` and multiplied by one factor, so their
  distribution keeps its shape and sign (a Gamma `da` stays positive).
- The mean is set by the intercept of the returned `TraitQTL`
  (`BnGStructs.TraitQTL(...; intercept)`), not by shifting the effects.
- `t.σₐ` is the standard deviation (with the `n − 1` divisor) of the
  genetic values in `hap`. It includes the covariances between QTL from
  linkage disequilibrium in that sample; the genic standard deviation,
  ``\\sqrt{\\sum 2p(1-p)a^2}``, differs from it. [`trait_summary`](@ref)
  reports both.

The realised values are checked after scaling and an `ArgumentError` is
thrown if they miss the target by more than `atol` (relative), or if the
QTL are monomorphic in `hap`.
"""
function eqtl(
    hap::Haplotype,
    loci::Union{LocusSet,AbstractVector{<:Integer}},
    t::AbstractTrait;
    rng::AbstractRNG = Random.default_rng(),
    atol::Real = 1e-8,
)
    loc = loci isa LocusSet ? loci.loci : Int32.(loci)
    isempty(loc) && throw(ArgumentError("no QTL loci"))
    return scale_qtl(hap, loc, Float64.(rand(rng, t.da, length(loc))), t; atol = atol)
end

"""
    scale_qtl(hap::Haplotype, loci, a, t::AbstractTrait; atol = 1e-8) -> TraitQTL

QTL of trait `t` at `loci` with the given effects `a` (of the allele
coded 1), multiplied by one factor so that the genetic values of the
individuals in `hap` have standard deviation `t.σₐ`, and mean 0 through
the intercept: [`eqtl`](@ref) with effects given instead of drawn, e.g.
those of a trait under selection in the history.
"""
function scale_qtl(hap::Haplotype, loc::AbstractVector{<:Integer}, a::AbstractVector{<:Real}, t::AbstractTrait;
    atol::Real = 1e-8)
    loc = Int32.(loc)
    length(a) == length(loc) || throw(DimensionMismatch("one effect per locus"))
    a = Float64.(a)
    Z = Float64.(dosage(hap, loc))                       # nqtl × nid
    g = Z' * a
    s = std(g)
    (isfinite(s) && s > 0) ||
        throw(ArgumentError("the sampled QTL give no genetic variance among the $(size(Z, 2)) individuals"))
    a .*= t.σₐ / s
    qtl = TraitQTL(t.name, loc, a; intercept = -mean(g) * t.σₐ / s)
    v = tbv(hap, qtl)
    (abs(mean(v)) <= atol * t.σₐ && abs(std(v) - t.σₐ) <= atol * t.σₐ) ||
        throw(ArgumentError("realised mean $(mean(v)) and sd $(std(v)) miss the target 0 and $(t.σₐ)"))
    return qtl
end

"""
    eqtl(hap::Haplotype, loci, R, ts::AbstractVector{<:AbstractTrait};
         rng = Random.default_rng(), atol = 1e-8) -> MultiTraitQTL

Sample additive QTL effects of several traits at shared `loci` so that the
genetic values of the individuals in `hap` have mean 0, standard
deviations `ts[i].σₐ` and correlation matrix `R`, exactly in that sample.

Effects of trait `i` are drawn from `ts[i].da`, independently across
traits. With ``\\mathbf{C}`` the sample covariance matrix of the resulting
genetic values and ``\\boldsymbol\\Sigma = \\mathbf{D}\\mathbf{R}\\mathbf{D}``
(``\\mathbf{D}`` the diagonal of the σₐ) the target, the effect matrix
``\\mathbf{A}`` is replaced by ``\\mathbf{A}\\mathbf{L}_C^{-\\top}\\mathbf{L}_\\Sigma^\\top``
(``\\mathbf{L}`` the Cholesky factors). The transform is upper triangular:
the effects of trait 1 are only rescaled and keep the shape of their
distribution; those of trait `i > 1` are linear combinations of the draws
for traits `1…i`, so only their first two moments are controlled. Means
go to the intercepts.

Needs at least as many QTL as traits and more individuals than traits.
The realised covariance is checked and an `ArgumentError` is thrown if it
misses `Σ` by more than `atol` (relative to the largest variance).
"""
function eqtl(
    hap::Haplotype,
    loci::Union{LocusSet,AbstractVector{<:Integer}},
    R::AbstractMatrix{<:Real},
    ts::AbstractVector{<:AbstractTrait};
    rng::AbstractRNG = Random.default_rng(),
    atol::Real = 1e-8,
)
    nt = length(ts)
    issymmetric(R) && size(R) == (nt, nt) ||
        throw(DimensionMismatch("R must be a symmetric $nt × $nt correlation matrix"))
    all(isone, diag(R)) || throw(ArgumentError("R must have a unit diagonal"))
    loc = loci isa LocusSet ? loci.loci : Int32.(loci)
    length(loc) >= nt || throw(ArgumentError("$(length(loc)) QTL cannot carry $nt traits"))
    Z = Float64.(dosage(hap, loc))                       # nqtl × nid
    size(Z, 2) > nt || throw(ArgumentError("need more individuals than traits"))
    A = reduce(hcat, [Float64.(rand(rng, t.da, length(loc))) for t in ts])
    G = Z' * A
    FC = cholesky(Symmetric(cov(G)); check = false)
    issuccess(FC) || throw(ArgumentError(
        "the sampled genetic values of the traits are collinear in the base population"))
    σ = [Float64(t.σₐ) for t in ts]
    Σ = Symmetric(Diagonal(σ) * Matrix{Float64}(R) * Diagonal(σ))
    FΣ = cholesky(Σ; check = false)
    issuccess(FΣ) || throw(ArgumentError("R is not positive definite"))
    T = FC.U \ Matrix(FΣ.U)                             # L_C^{-T} L_Σ^T
    A = A * T
    μ = vec(mean(G; dims = 1) * T)                      # means of the transformed values
    qtl = MultiTraitQTL([t.name for t in ts], loc, A; intercept = -μ)
    V = tbv(hap, qtl)
    (maximum(abs, mean(V; dims = 1)) <= atol * maximum(σ) &&
     maximum(abs, cov(V) .- Σ) <= atol * maximum(σ)^2) ||
        throw(ArgumentError("realised genetic covariance misses the target"))
    return qtl
end

"""
    trait_summary(hap::Haplotype, qtl::Union{TraitQTL, MultiTraitQTL}) -> NamedTuple

Genetic summary of the individuals in `hap` for the traits of `qtl`
(additive architecture), one entry per trait:

- `mean`, `var`: mean and variance (`n − 1` divisor) of the genetic values;
- `genic`: genic variance ``\\sum_j 2p_j(1-p_j)a_j^2``, with ``p_j`` the
  allele frequency at QTL ``j`` in `hap`; `var − genic` is the
  contribution of linkage disequilibrium (and of the `n − 1` divisor);
- `cor`: correlation matrix of the genetic values (`1 × 1` for one trait).
"""
function trait_summary(hap::Haplotype, qtl::Union{TraitQTL,MultiTraitQTL})
    loc = qtl.loci
    Z = dosage(hap, loc)
    p = vec(sum(Z; dims = 2)) ./ (2size(Z, 2))
    A = qtl isa TraitQTL ? reshape(qtl.additive, :, 1) : qtl.additive
    V = tbv(hap, qtl)
    V = V isa AbstractVector ? reshape(V, :, 1) : V
    (mean = vec(mean(V; dims = 1)), var = vec(var(V; dims = 1)),
        genic = vec(sum(2 .* p .* (1 .- p) .* A .^ 2; dims = 1)), cor = cor(V))
end
