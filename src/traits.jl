"""
    norm_mv!(Q::AbstractMatrix, v::AbstractVector{Float64}; ϵ::Real = 1e-6, μ::Real = 0.0, σ::Real = 1.0) -> Vector{Float64}

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
function _qtl_dosage_matrix(hap::Haplotype, loci::AbstractVector{<:Integer})
    nid = hap.nhp ÷ 2
    nqtl = length(loci)
    Q = Matrix{Int8}(undef, nqtl, nid)
    Threads.@threads for i = 1:nid
        c1, c2 = 2i - 1, 2i
        @inbounds for (k, l) in enumerate(loci)
            Q[k, i] = Int8(hap.gt[l, c1]) + Int8(hap.gt[l, c2])
        end
    end
    return Q
end

"""
    eqtl(
        hap::Haplotype,
        loci::Union{LocusSet, AbstractVector{<:Integer}},
        t::AbstractTrait;
        rng::AbstractRNG = Random.default_rng(),
    ) -> TraitQTL

Simulate additive QTL effects for a single `BnGStructs.AbstractTrait` `t` at `loci`,
normalized so that true breeding values in `hap` have standard deviation `t.σₐ`, and return a
`BnGStructs.TraitQTL`.
"""
function eqtl(
    hap::Haplotype,
    loci::Union{LocusSet,AbstractVector{<:Integer}},
    t::AbstractTrait;
    rng::AbstractRNG = Random.default_rng(),
)
    loc_vec = loci isa LocusSet ? loci.loci : Int32.(loci)
    Q = _qtl_dosage_matrix(hap, loc_vec)
    q = Float64.(rand(rng, t.da, length(loc_vec)))
    norm_mv!(Q, q; σ = t.σₐ)
    return TraitQTL(t.name, loc_vec, q)
end

"""
    eqtl(
        hap::Haplotype,
        loci::Union{LocusSet, AbstractVector{<:Integer}},
        R::AbstractMatrix{<:Real},
        ts::AbstractVector{<:AbstractTrait};
        rng::AbstractRNG = Random.default_rng(),
    ) -> MultiTraitQTL

Simulate correlated additive QTL effects across multiple traits `ts` at shared `loci`
with target genetic correlation matrix `R`, normalizing each trait `i` so its true breeding
values in `hap` have standard deviation `ts[i].σₐ`. Returns a `BnGStructs.MultiTraitQTL`.
"""
function eqtl(
    hap::Haplotype,
    loci::Union{LocusSet,AbstractVector{<:Integer}},
    R::AbstractMatrix{<:Real},
    ts::AbstractVector{<:AbstractTrait};
    rng::AbstractRNG = Random.default_rng(),
)
    nt = length(ts)
    issymmetric(R) && size(R) == (nt, nt) ||
        throw(DimensionMismatch("R must be a symmetric $nt × $nt correlation matrix"))
    loc_vec = loci isa LocusSet ? loci.loci : Int32.(loci)
    nqtl = length(loc_vec)
    Q = _qtl_dosage_matrix(hap, loc_vec)
    U = cholesky(Symmetric(Matrix{Float64}(R))).U
    qs = randn(rng, nqtl, nt) * U
    for (i, t) in enumerate(ts)
        norm_mv!(Q, view(qs, :, i); σ = t.σₐ)
    end
    return MultiTraitQTL([t.name for t in ts], loc_vec, qs)
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
    ped[!, "tbv_" * qtl.name] = vals
    return ped
end

function tbv!(ped::DataFrame, hap::Union{Haplotype,Genotype}, mqtl::MultiTraitQTL)
    vals = tbv(hap, mqtl)
    nrow(ped) == size(vals, 1) ||
        throw(DimensionMismatch("ped has $(nrow(ped)) rows, hap has $(size(vals, 1)) individuals"))
    for (j, name) in enumerate(mqtl.trait_names)
        ped[!, "tbv_" * name] = vals[:, j]
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
        gt = ped[!, "tbv_" * t.name]
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
            ped[!, "ft_" * t.name] = tt
        else
            for i = 1:nid
                if !ismissing(ft[i])
                    ft[i] = ft[i] + t.μ
                end
            end
            ped[!, "ft_" * t.name] = ft
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
