"""
    norm_mv!(Q::AbstractMatrix, v::AbstractVector{Float64}; ϵ::Real = 1e-6, μ::Real = 0.0, σ::Real = 1.0) -> Vector{Float64}

Deprecated: [`eqtl`](@ref) no longer uses it, because shifting every
effect to centre the genetic values changes the shape of the effect
distribution. Use the intercept of a `TraitQTL` instead.

Normalize additive QTL effect vector `v` in-place such that the genetic values `Q' * v`
(where `Q` is the `nqtl × nid` dosage matrix) have sample mean `μ` and standard deviation `σ`
within tolerance `ϵ`. Returns the normalized breeding values `Q' * v`.
"""
function norm_mv!(
    Q::AbstractMatrix,
    v::AbstractVector{Float64};
    ϵ::Real = 1e-6,
    μ::Real = 0.0,
    σ::Real = 1.0,
)
    nqtl = size(Q, 1)
    nqtl > 0 || throw(ArgumentError("Q must have at least one QTL row"))
    vals = Q' * v
    m, s = mean(vals), std(vals)
    s > 0 || throw(ArgumentError("QTL dosages have zero variance across individuals"))
    for _ = 1:100
        (abs(m - μ) <= ϵ && abs(s - σ) <= ϵ) && break
        v .-= (m - μ) / nqtl
        v .*= σ / s
        vals = Q' * v
        m, s = mean(vals), std(vals)
    end
    return vals
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
    Z = Float64.(dosage(hap, loc))                       # nqtl × nid
    a = Float64.(rand(rng, t.da, length(loc)))
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

"""
    tbv!(ped::DataFrame, hap::Union{Haplotype, Genotype}, qtl::TraitQTL) -> DataFrame
    tbv!(ped::DataFrame, hap::Union{Haplotype, Genotype}, mqtl::MultiTraitQTL) -> DataFrame

Compute True Breeding Values (TBV) using `BnGStructs.tbv` and store them in `ped`
as column(s) `"tbv_" * trait_name`.
"""
function tbv!(ped::DataFrame, hap::Union{Haplotype,Genotype}, qtl::TraitQTL)
    vals = tbv(hap, qtl)
    nrow(ped) == length(vals) ||
        throw(DimensionMismatch("ped has $(nrow(ped)) rows, hap has $(length(vals)) individuals"))
    ped[!, trait_columns(qtl.name).tbv] = vals
    return ped
end

function tbv!(ped::DataFrame, hap::Union{Haplotype,Genotype}, mqtl::MultiTraitQTL)
    vals = tbv(hap, mqtl)
    nrow(ped) == size(vals, 1) ||
        throw(DimensionMismatch("ped has $(nrow(ped)) rows, hap has $(size(vals, 1)) individuals"))
    for (j, name) in enumerate(mqtl.trait_names)
        ped[!, trait_columns(name).tbv] = vals[:, j]
    end
    return ped
end

"""
    phenotype!(
        ped::DataFrame,
        trts::Union{AbstractTrait, AbstractVector{<:AbstractTrait}};
        R_e::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
        rng::AbstractRNG = Random.default_rng(),
    ) -> DataFrame

Generate phenotypes in `ped` (columns `"ft_" * t.name`) from true breeding values
(`"tbv_" * t.name`) for one or more `BnGStructs.AbstractTrait` objects (`aTrait` or `tTrait`).

- Environmental standard deviation for trait `t` is ``\\sigma_e = \\sigma_a \\sqrt{(1 - h^2) / h^2}``.
- If `R_e` (a residual correlation matrix across `trts`) is given, environmental deviates are drawn jointly with correlation `R_e`.
- Sex-limited traits (`t.sex ∈ (0, 1)`) set records for the non-expressing sex in `ped.sex` to `missing`.
- Threshold traits (`tTrait`) discretize the latent liability ``L = \\text{TBV} + e`` against `t.threshold` into categories `0, 1, ..., length(t.threshold)`.
"""
function phenotype!(
    ped::DataFrame,
    trt::AbstractTrait;
    R_e::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)
    return phenotype!(ped, [trt]; R_e = R_e, rng = rng)
end

function phenotype!(
    ped::DataFrame,
    trts::AbstractVector{<:AbstractTrait};
    R_e::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)
    nid = nrow(ped)
    nt = length(trts)
    σe = [t.σₐ * sqrt((1.0 - t.h²) / t.h²) for t in trts]

    E = if R_e === nothing
        randn(rng, nid, nt) .* reshape(σe, 1, nt)
    else
        size(R_e) == (nt, nt) && issymmetric(R_e) ||
            throw(DimensionMismatch("R_e must be a symmetric $nt × $nt matrix"))
        U = cholesky(Symmetric(Matrix{Float64}(R_e))).U
        (randn(rng, nid, nt) * U) .* reshape(σe, 1, nt)
    end

    has_sex = hasproperty(ped, :sex)
    for (j, t) in enumerate(trts)
        gt = ped[!, trait_columns(t).tbv]
        ft = Vector{Union{Missing,Float64}}(undef, nid)
        @inbounds for i = 1:nid
            ft[i] = gt[i] + E[i, j]
        end
        if t.sex != 2 && has_sex
            for i = 1:nid
                if Int(ped.sex[i]) == 1 - t.sex
                    ft[i] = missing
                end
            end
        end
        if t isa tTrait
            tt = Vector{Union{Missing,Int8}}(undef, nid)
            nt_cat = length(t.threshold)
            for i = 1:nid
                if ismissing(ft[i])
                    tt[i] = missing
                else
                    idx = findfirst(x -> x > ft[i], t.threshold)
                    tt[i] = isnothing(idx) ? Int8(nt_cat) : Int8(idx - 1)
                end
            end
            ped[!, trait_columns(t).phenotype] = tt
        else
            for i = 1:nid
                if !ismissing(ft[i])
                    ft[i] = ft[i] + t.μ
                end
            end
            ped[!, trait_columns(t).phenotype] = ft
        end
    end
    return ped
end

"""
    dummy_traits(n::Int, nlc::Int; rng::AbstractRNG = Random.default_rng())

Generate `n` continuous additive `BnGStructs.Trait` definitions (`1 ≤ n ≤ 5`) plus a
4-category threshold litter-size trait (`"lts"`), a positive-definite `(n + 1) × (n + 1)`
correlation matrix `R`, and a dictionary of index weights and variance ratios.
"""
function dummy_traits(n::Int, nlc::Int; rng::AbstractRNG = Random.default_rng())
    1 <= n <= 5 || error("n must be in [1, 5]")
    trts = AbstractTrait[]
    name = 'a'
    for _ = 1:n
        t = Trait(
            "t$name";
            h² = rand(rng) * 0.5 + 0.05,
            sex = rand(rng, 0:2),
            σₐ = rand(rng) + 0.5,
        )
        push!(trts, t)
        name += 1
    end
    pp = [0.05, 0.2, 0.35, 0.4]
    push!(trts, Trait("lts", pp; h² = 0.2, sex = 0))

    m = n + 1
    R = rand(rng, m, m) .* 2.0 .- 1.0
    R = R * R'
    d = sqrt.(diag(R))
    R ./= (d * d')
    for i = 1:m
        R[i, i] = 1.0
    end

    w = rand(rng, n)
    w ./= sum(w)
    itts = Dict{String,NamedTuple}()
    for i = 1:n
        t = trts[i]
        σᵤ² = t.σₐ^2
        σₐ² = 2σᵤ² / nlc
        σₑ² = σᵤ² / t.h² * (1 - t.h²)
        λ = σₑ² / σₐ²
        itts[t.name] = (w = w[i], λ = λ, h² = t.h²)
    end

    return trts, R, itts
end
