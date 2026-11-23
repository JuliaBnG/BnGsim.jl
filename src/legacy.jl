# The building blocks of BnGsim 0.1–0.3, outside the 1.0 API (package
# roadmap §11, decision a): they draw from an explicit or the global
# generator rather than from the keyed streams (invariant I1), so their
# results depend on the order of calls. The breeding loop, the scenarios
# and the command line do not use them. Load them with
# `using BnGsim.Legacy`.
"""
    BnGsim.Legacy

The building blocks of BnGsim 0.1–0.3: a recombination and gene-drop
toolkit (`sum_map`, `crossovers(sp, vm)`, `gamete!`, `gene_drop!`,
`random_mate`, `rnd_ped`, `expand`, `sample_id`), trait helpers
(`norm_mv!`, `tbv!`, `phenotype!`, `dummy_traits`), founder sampling
(`FounderPopulation`, `sample_founder_population`), historical
populations (`HistoricalSpec`, `historical_population`,
`continue_population`, `read_packed_haplotypes`, `write_xy`),
`vanraden_g` and the text `provenance` header.

They draw from an explicit `rng` or the global generator, not from the
keyed streams of the breeding loop, so their results depend on the order
of calls and are outside the reproducibility invariants and the 1.0
promise of BnGsim. Use them with `using BnGsim.Legacy`.
"""
module Legacy

using BnGStructs
using Dates
using LinearAlgebra
using DataFrames
using Distributions
using Random
using Statistics
using StatsBase: sample
using ..BnGsim
using ..BnGsim: _to_variant_map
using ..FisherWrightAdapter
using ..RelationshipsAdapter
using ..BreedingAdapter
import ..BnGsim: crossovers

export sum_map, gamete!, gene_drop!, random_mate, rnd_ped, expand, sample_id
export norm_mv!, tbv!, phenotype!, dummy_traits
export FounderPopulation, sample_founder_population
export HistoricalSpec, historical_population, continue_population, read_packed_haplotypes, write_xy
export vanraden_g, provenance

# ---- recombination and gene drop (0.1) ----

"""
    sum_map(sp::Species) -> DataFrame
    sum_map(chr_len::AbstractVector{<:Integer}; M::Real = 1e8) -> DataFrame
    sum_map(vm::VariantMap; cM::Real = 1e6) -> DataFrame

Construct a chromosome linkage map summary `DataFrame`.

- For `chr_len` or `sp::Species`, returns columns `chr::Int8`, `len::Float64` (length in Morgans), and `ebp::UInt32` (cumulative ending base-pair coordinate).
- For `vm::VariantMap`, summarizes per-chromosome span in Morgans (`len`), locus count (`nlc::UInt32`), and 1-based starting locus index (`bgn::UInt32`).
"""
function sum_map(chr_len::AbstractVector{<:Integer}; M::Real = 1e8)
    nchr = length(chr_len)
    ebp = UInt32.(cumsum(chr_len))
    return DataFrame(
        chr = Int8.(1:nchr),
        len = Float64.(chr_len) ./ Float64(M),
        ebp = ebp,
    )
end

# Genetic lengths come from the species (`map_length`), which defaults to
# `chromosome ./ M` and may be set per chromosome.
sum_map(sp::Species) = DataFrame(
    chr = Int8.(1:length(sp)),
    len = copy(map_length(sp)),
    ebp = copy(cbp(sp)),
)

function sum_map(vm::VariantMap; cM::Real = 1e6)
    lms = DataFrame(
        chr = Int8[],
        len = Float64[],
        nlc = UInt32[],
        bgn = UInt32[],
    )
    n = length(vm)
    n == 0 && return lms
    i = 1
    while i <= n
        c = vm.chr[i]
        j = i
        while j < n && vm.chr[j + 1] == c
            j += 1
        end
        len = (Float64(vm.pos[j]) - Float64(vm.pos[i])) / (100.0 * Float64(cM))
        nlc = UInt32(j - i + 1)
        push!(lms, (c, len, nlc, UInt32(i)))
        i = j + 1
    end
    return lms
end


"""
    crossovers(lms::DataFrame, tbp::AbstractVector{<:Integer}; rng::AbstractRNG = Random.default_rng()) -> (Vector{UInt32}, Vector{Int32})
    crossovers(sp::Species, vm::VariantMap; rng::AbstractRNG = Random.default_rng()) -> (Vector{UInt32}, Vector{Int32})

Sample Poisson-distributed crossover positions along the genome without interference,
with independent assortment (`p = 0.5`) at chromosome boundaries.

Returns `(gcvs, mcps)`:
- `gcvs::Vector{UInt32}`: Leading element is the starting parental homologue (`1` or `2`), followed by genome-wide crossover base-pair positions.
- `mcps::Vector{Int32}`: Leading element is the starting parental homologue (`1` or `2`), followed by right-inclusive locus indices in `tbp` / `vm.pos`.
"""
function crossovers(
    lms::DataFrame,
    tbp::AbstractVector{<:Integer};
    rng::AbstractRNG = Random.default_rng(),
)
    start_hap = UInt32(rand(rng, 1:2))
    gcvs = UInt32[start_hap]
    mcps = Int32[Int32(start_hap)]
    bgn = UInt32(1)
    for row in eachrow(lms)
        λ = row.len
        ebp = UInt32(row.ebp)
        nc = rand(rng, Poisson(λ))
        if nc > 0
            pts = [rand(rng, bgn:ebp) for _ = 1:nc]
            sort!(pts)
            append!(gcvs, pts)
        end
        rand(rng, Bool) && push!(gcvs, ebp)
        bgn = ebp + UInt32(1)
    end
    !isempty(lms) && gcvs[end] == lms.ebp[end] && pop!(gcvs)
    for bp in view(gcvs, 2:length(gcvs))
        i = searchsortedlast(tbp, bp)
        push!(mcps, Int32(i))
    end
    return gcvs, mcps
end

function crossovers(sp::Species, vm::VariantMap; rng::AbstractRNG = Random.default_rng())
    check_coordinates(vm, sp)
    return crossovers(sum_map(sp), vm.pos; rng = rng)
end

"""
    gamete!(
        hap::AbstractVector,
        prt::AbstractMatrix,
        lms::DataFrame,
        tbp::AbstractVector{<:Integer};
        rng::AbstractRNG = Random.default_rng(),
    ) -> Vector{UInt32}

Generate a recombinant gamete `hap` in-place from a parent's two haplotypes `prt` (`nlc × 2` view),
returning the genome-wide crossover positions `gcs`.
"""
function gamete!(
    hap::AbstractVector,
    prt::AbstractMatrix,
    lms::DataFrame,
    tbp::AbstractVector{<:Integer};
    rng::AbstractRNG = Random.default_rng(),
)
    gcs, mcs = crossovers(lms, tbp; rng = rng)
    nlc = size(prt, 1)
    h = Int(mcs[1])
    l = 1
    @inbounds for cv in view(mcs, 2:length(mcs))
        if cv >= l
            copyto!(view(hap, l:cv), view(prt, l:cv, h))
            l = Int(cv) + 1
        end
        h = 3 - h
    end
    if l <= nlc
        copyto!(view(hap, l:nlc), view(prt, l:nlc, h))
    end
    return gcs
end

"""
    gene_drop!(
        og::Haplotype,
        pg::Haplotype,
        pm::AbstractMatrix{<:Integer},
        sp::Species,
        vm::VariantMap;
        rng::AbstractRNG = Random.default_rng(),
    ) -> DataFrame
    gene_drop!(
        og::AbstractMatrix,
        pg::AbstractMatrix,
        pm::AbstractMatrix{<:Integer},
        lms::DataFrame,
        tbp::AbstractVector{<:Integer};
        rng::AbstractRNG = Random.default_rng(),
    ) -> DataFrame

Drop recombined parental gametes from parental haplotypes `pg` into offspring haplotypes `og`
according to the `nf × 2` parent matrix `pm` (`[sire dam]` per offspring).

Uses per-offspring `Xoshiro` streams seeded deterministically from `rng`, so results are
invariant to `Threads.nthreads()`. Returns a `DataFrame` of crossover records (`generation`, `parent`, `ihap`, `crossovers`).
"""
function gene_drop!(
    og::AbstractMatrix,
    pg::AbstractMatrix,
    pm::AbstractMatrix{<:Integer},
    lms::DataFrame,
    tbp::AbstractVector{<:Integer};
    rng::AbstractRNG = Random.default_rng(),
)
    nf = size(pm, 1)
    size(og, 2) >= 2nf || throw(DimensionMismatch("og must have at least $(2nf) haplotype columns"))
    size(pg, 1) == size(og, 1) == length(tbp) ||
        throw(DimensionMismatch("locus count mismatch between pg, og, and tbp"))

    # Seed one independent RNG stream per gamete (2 * nf) for thread-count invariance
    streams = [Xoshiro(rand(rng, UInt64)) for _ = 1:(2nf)]
    gcs_list = Vector{Vector{UInt32}}(undef, 2nf)

    # Note: BitMatrix columns may share 64-bit chunks when nlc is not a multiple of 64.
    # Thread over individuals when each haplotype is a separate column of a non-BitMatrix
    # or when nlc % 64 == 0; otherwise copy sequentially or into chunk-safe views.
    if !(og isa BitMatrix) || (size(og, 1) % 64 == 0)
        Threads.@threads for id = 1:nf
            ip = pm[id, 1]
            im = pm[id, 2]
            pa = view(pg, :, (2ip - 1):(2ip))
            ma = view(pg, :, (2im - 1):(2im))
            gcs_list[2id - 1] = gamete!(view(og, :, 2id - 1), pa, lms, tbp; rng = streams[2id - 1])
            gcs_list[2id] = gamete!(view(og, :, 2id), ma, lms, tbp; rng = streams[2id])
        end
    else
        for id = 1:nf
            ip = pm[id, 1]
            im = pm[id, 2]
            pa = view(pg, :, (2ip - 1):(2ip))
            ma = view(pg, :, (2im - 1):(2im))
            gcs_list[2id - 1] = gamete!(view(og, :, 2id - 1), pa, lms, tbp; rng = streams[2id - 1])
            gcs_list[2id] = gamete!(view(og, :, 2id), ma, lms, tbp; rng = streams[2id])
        end
    end

    parents = Vector{Int32}(undef, 2nf)
    ihaps = Vector{Int8}(undef, 2nf)
    xovers = Vector{Vector{UInt32}}(undef, 2nf)
    for id = 1:nf
        g_p = gcs_list[2id - 1]
        g_m = gcs_list[2id]
        parents[2id - 1] = Int32(pm[id, 1])
        parents[2id] = Int32(pm[id, 2])
        ihaps[2id - 1] = Int8(g_p[1])
        ihaps[2id] = Int8(g_m[1])
        xovers[2id - 1] = g_p[2:end]
        xovers[2id] = g_m[2:end]
    end

    return DataFrame(
        generation = zeros(Int16, 2nf),
        parent = parents,
        ihap = ihaps,
        crossovers = xovers,
    )
end

function gene_drop!(
    og::Haplotype,
    pg::Haplotype,
    pm::AbstractMatrix{<:Integer},
    sp::Species,
    vm::VariantMap;
    rng::AbstractRNG = Random.default_rng(),
)
    check_coordinates(vm, sp)
    # the bit storage is padded to a multiple of 64 rows; pass the loci only
    nlc = og.nlc
    return gene_drop!(view(og.gt, 1:nlc, :), view(pg.gt, 1:nlc, :), pm, sum_map(sp), vm.pos; rng = rng)
end

"""
    random_mate(n1::Integer, n2::Integer; rng::AbstractRNG = Random.default_rng()) -> (Matrix{Int}, Vector{Bool})
    random_mate(sex::AbstractVector{Bool}, n2::Integer; rng::AbstractRNG = Random.default_rng()) -> Matrix{Int}

Randomly draw sire and dam pairs from `n1` parents (or parents with sexes `sex`) to produce `n2` offspring.
"""
function random_mate(n1::T, n2::T; rng::AbstractRNG = Random.default_rng()) where {T<:Integer}
    (n1 > 1 && n2 > 0) || throw(ArgumentError("n1 > 1 and n2 > 0 required"))
    sex = rand(rng, Bool, n1)
    while all(sex) || !any(sex)
        rand!(rng, sex)
    end
    return random_mate(sex, n2; rng = rng), sex
end

function random_mate(
    sex::AbstractVector{Bool},
    n2::T;
    rng::AbstractRNG = Random.default_rng(),
) where {T<:Integer}
    n2 > 0 || throw(ArgumentError("n2 must be a positive integer"))
    sire_pool = findall(sex)
    dam_pool = findall(.!sex)
    (!isempty(sire_pool) && !isempty(dam_pool)) ||
        throw(ArgumentError("both sexes must be present in parent cohort"))
    sirs = rand(rng, sire_pool, n2)
    dams = rand(rng, dam_pool, n2)
    return sortslices([sirs dams]; dims = 1)
end

"""
    rnd_ped(nn::Integer, ng::Integer; rng::AbstractRNG = Random.default_rng()) -> DataFrame

Generate a random-mating discrete-generation pedigree with `nn` individuals per generation
across a founder cohort (`grt = 0`) and `ng` descendant generations (`1:ng`).
"""
function rnd_ped(nn::Integer, ng::Integer; rng::AbstractRNG = Random.default_rng())
    0 <= ng <= 255 || throw(ArgumentError("ng must be in 0:255"))
    nn > 1 || throw(ArgumentError("nn must exceed 1"))
    sex0 = rand(rng, Bool, nn)
    while all(sex0) || !any(sex0)
        rand!(rng, sex0)
    end
    ped = DataFrame(
        id = zeros(Int32, nn),
        sire = zeros(Int32, nn),
        dam = zeros(Int32, nn),
        sex = sex0,
        grt = fill(UInt8(0), nn),
    )
    for i = 1:ng
        sex_prev = ped.sex[ped.grt .== UInt8(i - 1)]
        pm = random_mate(sex_prev, nn; rng = rng)
        acc = Int32((i - 1) * nn)
        sex_i = rand(rng, Bool, nn)
        while all(sex_i) || !any(sex_i)
            rand!(rng, sex_i)
        end
        off = DataFrame(
            id = zeros(Int32, nn),
            sire = Int32.(pm[:, 1]) .+ acc,
            dam = Int32.(pm[:, 2]) .+ acc,
            sex = sex_i,
            grt = fill(UInt8(i), nn),
        )
        append!(ped, off)
    end
    ped.id = Int32.(1:nrow(ped))
    return ped
end

"""
    expand(
        hap::Haplotype,
        vm::VariantMap,
        sp::Species,
        n_t::Integer,
        t::Integer;
        maf::Real = 0.0,
        rng::AbstractRNG = Random.default_rng(),
    ) -> (Haplotype, VariantMap, DataFrame)

Exponentially expand a founder population `hap` (`n_0 = hap.nhp ÷ 2` individuals) to `n_t`
individuals over `t` discrete generations (`n_i = round(Int, n_0 * (n_t / n_0)^(i / t))`),
updating allele frequencies in the returned `VariantMap` and filtering loci by `maf`.
"""
function expand(
    hap::Haplotype,
    vm::VariantMap,
    sp::Species,
    n_t::Integer,
    t::Integer;
    maf::Real = 0.0,
    rng::AbstractRNG = Random.default_rng(),
)
    nlc = hap.nlc
    nhp = hap.nhp
    n0 = nhp ÷ 2
    nlc == length(vm) || throw(DimensionMismatch("Haplotype ($nlc loci) and VariantMap ($(length(vm)) loci) do not match"))
    t >= 1 || throw(ArgumentError("t must be at least 1"))
    n_t >= 2 || throw(ArgumentError("n_t must be at least 2"))

    check_coordinates(vm, sp)
    lms = sum_map(sp)
    p = (n_t / n0)^(1 / t)
    # Haplotype pads its BitMatrix to 64·cld(nlc, 64) rows; the loci only.
    prt = hap.gt[1:nlc, :]
    np = n0
    cps = DataFrame()

    for i = 1:t
        ni = i == t ? Int(n_t) : max(2, round(Int, n0 * p^i))
        off = BitMatrix(undef, nlc, 2ni)
        sex_p = rand(rng, Bool, np)
        while all(sex_p) || !any(sex_p)
            rand!(rng, sex_p)
        end
        pm = random_mate(sex_p, ni; rng = rng)
        cp = gene_drop!(off, prt, pm, lms, vm.pos; rng = rng)
        cp.generation .= Int16(i)
        append!(cps, cp)
        prt = off
        np = ni
    end

    nht = size(prt, 2)
    frq = Float32.(vec(sum(prt; dims = 2)) ./ nht)
    vld = if 0.0 <= maf < 0.5
        findall(f -> maf < f < 1.0 - maf, frq)
    else
        collect(1:nlc)
    end
    vm_out = VariantMap(vm.chr[vld], vm.pos[vld], vm.ref[vld], vm.alt[vld], frq[vld])
    return Haplotype(prt[vld, :]), vm_out, cps
end

"""
    sample_id(
        hap::Haplotype,
        vm::VariantMap,
        ne::Integer;
        maf::Real = 0.0,
        rng::AbstractRNG = Random.default_rng(),
    ) -> (Haplotype, VariantMap)

Sample `ne` diploid individuals without replacement from `hap`, recompute sample allele
frequencies in the returned `VariantMap`, and retain loci satisfying `maf < frq < 1 - maf`.
"""
function sample_id(
    hap::Haplotype,
    vm::VariantMap,
    ne::Integer;
    maf::Real = 0.0,
    rng::AbstractRNG = Random.default_rng(),
)
    0.0 <= maf < 0.5 || throw(ArgumentError("maf must be in [0.0, 0.5)"))
    nlc = hap.nlc
    nlc == length(vm) || throw(DimensionMismatch("Haplotype and VariantMap locus counts do not match"))
    nid = hap.nhp ÷ 2
    0 < ne <= nid || throw(ArgumentError("ne = $ne must lie in 1:$nid"))

    ids = sort(sample(rng, 1:nid, ne; replace = false))
    cols = Vector{Int}(undef, 2ne)
    for (k, id) in enumerate(ids)
        cols[2k - 1] = 2id - 1
        cols[2k] = 2id
    end
    sub_gt = hap.gt[1:nlc, cols]
    frq = Float32.(vec(mean(sub_gt; dims = 2)))
    vld = findall(f -> maf < f < 1.0 - maf, frq)
    vm_out = VariantMap(vm.chr[vld], vm.pos[vld], vm.ref[vld], vm.alt[vld], frq[vld])
    return Haplotype(sub_gt[vld, :]), vm_out
end


# ---- traits (0.1) ----

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


# ---- founder sampling (0.1) ----


"""
    struct FounderPopulation
Holds the complete state of a sampled founder population for breeding simulations.

# Fields
- `species::Species`: Biological species definition (e.g. `Cattle`, `Pig`).
- `haplotype::Haplotype`: Base genome haplotypes (`nlc × 2nid`).
- `map::VariantMap`: Variant genomic map (`VariantMap` from `BnGStructs`).
- `markers::LocusSet`: Marker chip panel (`LocusSet` from `BnGStructs`).
- `qtl::TraitQTL`: Single-trait QTL architecture.
- `h2::Float64`: Target heritability.
- `tbv::Vector{Float64}`: True Breeding Values for the `n` founder individuals.
- `phenotypes::Vector{Float64}`: Simulated phenotypes (`y = TBV + e`).
- `pedigree::DataFrame`: Pedigree with columns `:id`, `:sire`, `:dam` (`sire`/`dam` = 0 for founders).

Note: `fp.marker_indices` is also supported as a property alias returning `fp.markers.loci`.
"""
struct FounderPopulation
    species::Species
    haplotype::Haplotype
    map::VariantMap
    markers::LocusSet
    qtl::TraitQTL
    h2::Float64
    tbv::Vector{Float64}
    phenotypes::Vector{Float64}
    pedigree::DataFrame
end

function Base.getproperty(fp::FounderPopulation, s::Symbol)
    s === :marker_indices && return getfield(fp, :markers).loci
    return getfield(fp, s)
end


"""
    sample_founder_population(
        species::Species;
        ne::Integer,
        n_id::Int = 100,
        n_marker::Int = 1000,
        n_qtl::Int = 100,
        h2::Float64 = 0.5,
        trait_name::String = "trait1",
        qtl_dist::Distribution = Normal(0.0, 1.0),
        burnin_generations::Int = 500,
        mutation_rate::Float64 = 1.0,
        rng::AbstractRNG = Random.default_rng(),
        verbose::Bool = false,
    ) -> FounderPopulation

Simulate a base founder population with `n_marker` markers, `n_qtl` QTLs, and heritability `h2`
by running a Fisher-Wright burn-in of `burnin_generations` generations of `ne` diploid individuals
on the genome of `species` (a species carries no population size since BnGStructs 0.4, so
`ne` is required).
"""
function sample_founder_population(
    species::Species;
    ne::Integer,
    n_id::Int = 100,
    n_marker::Int = 1000,
    n_qtl::Int = 100,
    h2::Float64 = 0.5,
    trait_name::String = "trait1",
    qtl_dist::Distribution = Normal(0.0, 1.0),
    burnin_generations::Int = 500,
    mutation_rate::Float64 = 1.0,
    rng::AbstractRNG = Random.default_rng(),
    verbose::Bool = false,
)
    0.0 < h2 < 1.0 || throw(ArgumentError("Heritability h2 must be in (0, 1)"))
    n_id > 0 || throw(ArgumentError("n_id must be positive"))

    res = fisher_wright(
        species,
        Int(ne),
        burnin_generations,
        mutation_rate;
        verbose = verbose,
        result = true,
        rng = rng,
    )

    bit_gt, map_df = dense_haplotypes(res.active_haplotypes, res.chromosome_ends; rng=copy(rng))
    vm = _to_variant_map(map_df)
    return sample_founder_population(
        species,
        Haplotype(bit_gt),
        vm;
        n_id = n_id,
        n_marker = n_marker,
        n_qtl = n_qtl,
        h2 = h2,
        trait_name = trait_name,
        qtl_dist = qtl_dist,
        rng = rng,
    )
end

"""
    sample_founder_population(
        species::Species,
        hap::Haplotype,
        vm_or_df::Union{VariantMap, DataFrame};
        n_id::Int = hap.nhp ÷ 2,
        n_marker::Int = 1000,
        n_qtl::Int = 100,
        h2::Float64 = 0.5,
        trait_name::String = "trait1",
        qtl_dist::Distribution = Normal(0.0, 1.0),
        rng::AbstractRNG = Random.default_rng(),
    ) -> FounderPopulation

Construct a founder population directly from an existing pre-computed `Haplotype` and `VariantMap` (or map `DataFrame`).
"""
function sample_founder_population(
    species::Species,
    hap::Haplotype,
    vm_or_df::Union{VariantMap,DataFrame};
    n_id::Int = hap.nhp ÷ 2,
    n_marker::Int = 1000,
    n_qtl::Int = 100,
    h2::Float64 = 0.5,
    trait_name::String = "trait1",
    qtl_dist::Distribution = Normal(0.0, 1.0),
    rng::AbstractRNG = Random.default_rng(),
)
    0.0 < h2 < 1.0 || throw(ArgumentError("Heritability h2 must be in (0, 1)"))
    vm = _to_variant_map(vm_or_df)
    total_avail_id = hap.nhp ÷ 2
    n_id <= total_avail_id ||
        throw(ArgumentError("Requested n_id ($n_id) exceeds available individuals in haplotype ($total_avail_id)"))

    nlc = hap.nlc
    (n_marker + n_qtl) <= nlc ||
        throw(ArgumentError("Requested markers ($n_marker) + QTLs ($n_qtl) exceeds total loci ($nlc)"))

    if n_id < total_avail_id
        id_indices = sample(rng, 1:total_avail_id, n_id; replace = false, ordered = true)
        sampled_cols = Int[]
        for idx in id_indices
            push!(sampled_cols, 2idx - 1)
            push!(sampled_cols, 2idx)
        end
        hap_sub = Haplotype(hap.gt[1:nlc, sampled_cols])
    else
        hap_sub = hap
    end

    all_loci = collect(Int32(1):Int32(nlc))
    sampled_loci = sample(rng, all_loci, n_marker + n_qtl; replace = false)
    qtl_loci = sort(sampled_loci[1:n_qtl])
    marker_loci = sort(sampled_loci[(n_qtl + 1):end])

    markers_set = LocusSet("markers", marker_loci)
    qtl_effects = rand(rng, qtl_dist, n_qtl)
    trait_qtl = TraitQTL(trait_name, qtl_loci, qtl_effects)

    tbvs = tbv(hap_sub, trait_qtl)
    var_g = var(tbvs)
    if var_g <= 0.0 || isnan(var_g)
        var_g = 1e-6
    end
    var_e = var_g * ((1.0 - h2) / h2)
    noise = rand(rng, Normal(0.0, sqrt(var_e)), n_id)
    phenotypes = tbvs .+ noise

    ped = DataFrame(
        id = collect(Int32(1):Int32(n_id)),
        sire = zeros(Int32, n_id),
        dam = zeros(Int32, n_id),
    )

    return FounderPopulation(
        species,
        hap_sub,
        vm,
        markers_set,
        trait_qtl,
        h2,
        tbvs,
        phenotypes,
        ped,
    )
end


# ---- historical populations (0.2) ----

"""
    struct HistoricalSpec

A historical (pre-founder) population run with `FisherWright.fisher_wright`.

# Fields
- `name::String`: arm label, e.g. `"F2"`.
- `ne::Int`: diploid population size, constant over the run.
- `nt::Int`: number of discrete, non-overlapping generations.
- `mr::Float64`: mutations per `1e8` bp per meiosis (μ = `mr × 1e-8`).
- `chromosome::Vector{UInt32}`: autosome lengths in bp.
- `M::Float64`: bp per Morgan (`1e8` is 1 cM/Mb).

Model, as implemented by FisherWright with its neutral defaults: diploid, dioecious Wright-Fisher
population with random union of gametes; infinite-sites mutation; Poisson
crossovers without interference; independent assortment of chromosomes;
the genome starts monomorphic, so `nt` determines how close the run is to
mutation-drift equilibrium (about `10ne` generations for the frequency
spectrum, much less for LD at short distances).
"""
struct HistoricalSpec
    name::String
    ne::Int
    nt::Int
    mr::Float64
    chromosome::Vector{UInt32}
    M::Float64
end

"""
    HistoricalSpec(name, sp::Species, ne, nt, mr)

The run of `ne` diploid individuals for `nt` generations on the genome of
`sp` (its chromosome lengths and bp per Morgan). The size is an argument:
a species carries none since BnGStructs 0.4.
"""
HistoricalSpec(name::AbstractString, sp::Species, ne::Integer, nt::Integer, mr::Real) =
    HistoricalSpec(String(name), Int(ne), Int(nt), Float64(mr), copy(sp.chromosome), Float64(sp.M))

"""
    historical_population(spec::HistoricalSpec; seed = nothing, rng = Xoshiro(seed),
                          verbose = false, progress_interval = 1000) -> FisherWrightResult

Run the historical population described by `spec`. Give either an integer
`seed`, which starts a fresh `Xoshiro(seed)` stream, or an `rng`. The
global RNG is not touched. FisherWright draws every haplotype from its
own stream seeded from `rng`, so the result does not depend on the number
of threads.

A given `seed` gives the same population as in BnGsim 0.1.0, where the
global RNG was seeded with it.
"""
function historical_population(
    spec::HistoricalSpec;
    seed::Union{Nothing,Integer} = nothing,
    rng::Union{Nothing,AbstractRNG} = nothing,
    verbose::Bool = false,
    progress_interval::Integer = 1000,
)
    (seed === nothing) == (rng === nothing) &&
        throw(ArgumentError("give exactly one of seed or rng"))
    return fisher_wright(
        spec.ne,
        spec.nt,
        spec.chromosome,
        spec.mr;
        M = spec.M,
        result = true,
        verbose = verbose,
        progress_interval = progress_interval,
        rng = rng === nothing ? Xoshiro(seed) : rng,
    )
end

"""
    continue_population(res::FisherWrightResult, ns, mr; M = 1e8, mut_base = 1e8,
                        rng = Random.default_rng()) -> FisherWrightResult

Continue a FisherWright population for `length(ns)` further generations with
census size `ns[g]` in generation `g`, e.g. a decline from the historical
size to a breed-like size before founders are sampled. This is
`fisher_wright(...; init = res, popsize = ns)`: the same model (dioecious
Wright–Fisher, Poisson crossovers at `1/M` per bp, infinite-sites mutation
at `mr / mut_base` per bp per meiosis) with a census size that may change
every generation.

The result does not depend on the number of threads. Mutations that fix
during the continuation move to `substitutions`, and the per-generation
records of `res` are extended. An empty `ns` returns `res`.

Realizations for a given `rng` differ from BnGsim 0.1.0, which ran its own
serial loop on FisherWright internals.
"""
function continue_population(
    res::FisherWrightResult,
    ns::AbstractVector{<:Integer},
    mr::Real;
    M::Real = 1e8,
    mut_base::Real = 1e8,
    rng::AbstractRNG = Random.default_rng(),
)
    all(>(1), ns) || throw(ArgumentError("every generation needs at least 2 individuals"))
    isempty(ns) && return res
    chr = diff([0; Int.(res.chromosome_ends)])
    return fisher_wright(
        length(res.active_haplotypes) ÷ 2,
        length(ns),
        chr,
        Float64(mr);
        M = M,
        mut_base = mut_base,
        result = true,
        init = res,
        popsize = collect(Int, ns),
        rng = rng,
    )
end

"""
    read_packed_haplotypes(path, cbp; rng = Random.default_rng())
        -> (xy::BitMatrix, lmp::DataFrame)

Read phased haplotypes written by an external simulator (e.g. the msprime
driver of the founder-design paper) in the `BNGH1` format: 8-byte magic
`"BNGH1\\0\\0\\0"`, `Int64` `nloci` and `nhap`, `UInt32` genome-wide
positions, then each locus as `cld(nhap, 8)` bytes with haplotype `h` in bit
`(h-1) % 8` of byte `(h-1) ÷ 8` (numpy `packbits(bitorder = "little")`).
1 is the derived allele; haplotypes `2i-1`, `2i` form individual `i`.

Returns the same layout as the `FisherWrightResult` method: loci segregating
in the sample, sorted by position, `chr` from the cumulative chromosome ends
`cbp`, `frq` the sample frequency of allele 1, and random `ref`/`alt` bases.
"""
function read_packed_haplotypes(
    path::AbstractString,
    cbp::AbstractVector{<:Integer};
    rng::AbstractRNG = Random.default_rng(),
)
    hap, pos = read_haplotypes(path)        # BnGStructs owns the BNGH1 format
    xy = packed(hap)[1:hap.nlc, :]
    frq = vec(sum(xy; dims = 2)) ./ size(xy, 2)
    keep = findall(f -> 0 < f < 1, frq)
    ord = keep[sortperm(pos[keep])]
    allunique(pos[ord]) || error("$path has duplicate positions")
    bases = ('A', 'C', 'G', 'T')
    ref = [rand(rng, bases) for _ in ord]
    alt = [rand(rng, filter(!=(r), bases)) for r in ref]
    lmp = DataFrame(
        chr = Int8.(searchsortedfirst.(Ref(cbp), pos[ord])),
        pos = pos[ord],
        ref = ref,
        alt = alt,
        frq = Float32.(frq[ord]),
    )
    return xy[ord, :], lmp
end

"""
    write_xy(path, xy::BitMatrix)

Write a locus-major `BitMatrix` in xyBnG's `xy` format (8-byte header
`xy-` with zero type bytes, `Int64` `nrow` and `ncol`, then the bits), the
format read by xyBnG's pipeline and by `founder_haplotypes(fxy, fmp, …)`.
"""
function write_xy(path::AbstractString, xy::BitMatrix)
    open(path, "w") do io
        write(io, Int8['x', 'y', '-', 0, 0, 0, 0, 0])
        write(io, Int[size(xy)...])
        write(io, xy)
    end
    return path
end


# ---- genomic relationships and run headers ----

"""
    vanraden_g(xy, loci) -> Matrix{Float64}

VanRaden (2008) method-1 genomic relationships among the individuals of `xy`
(columns `2i-1`, `2i`) at `loci`, centred on the allele frequencies of the
same individuals: `G = ZZ' / 2Σpq`. Computed by `RelationshipMatrices.grm`.
"""
function vanraden_g(xy::AbstractMatrix{Bool}, loci::AbstractVector{<:Integer})
    nid = size(xy, 2) ÷ 2
    Z = Matrix{Int8}(undef, length(loci), nid)
    Threads.@threads for i = 1:nid
        @inbounds for (k, l) in enumerate(loci)
            Z[k, i] = Int8(xy[l, 2i-1]) + Int8(xy[l, 2i])
        end
    end
    return grm(Z)
end

"""
    provenance(io::IO = stdout; seed = nothing, extra = Dict{String,Any}())

Write a run header recording what is needed to reproduce a simulation: the
date, host, CPU, Julia version, thread count, BLAS threads, the seed, and for
each JuliaBnG package its version, git commit and whether its working tree
was dirty. Lines are prefixed with `# ` so the header can sit on top of a
plain-text log or table.

Results from `FisherWright.fisher_wright` are reproducible for a given seed
only at a fixed thread count, which is why the thread count is recorded.
"""
function provenance(io::IO = stdout; seed = nothing, extra = Dict{String,Any}())
    println(io, "# date:     ", Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))
    println(io, "# host:     ", gethostname())
    println(io, "# cpu:      ", Sys.cpu_info()[1].model, " × ", Sys.CPU_THREADS)
    println(io, "# memory:   ", round(Sys.total_memory() / 2^30; digits = 1), " GiB")
    println(io, "# julia:    ", VERSION)
    println(io, "# threads:  ", Threads.nthreads(), " (BLAS ", BLAS.get_num_threads(), ")")
    seed === nothing || println(io, "# seed:     ", seed)
    for m in _provenance_modules()
        println(io, "# package:  ", rpad(string(m), 21), _pkg_state(m))
    end
    for (k, v) in extra
        println(io, "# ", rpad(string(k) * ":", 10), v)
    end
    flush(io)
    return nothing
end

function _provenance_modules()
    return Module[BnGStructs, FisherWrightAdapter.PACKAGE, RelationshipsAdapter.PACKAGE, BreedingAdapter.PACKAGE,
        BnGsim]
end

function _pkg_state(m::Module)
    dir = pkgdir(m)
    ver = pkgversion(m)
    commit = "no-git"
    if dir !== nothing && isdir(joinpath(dir, ".git"))
        try
            commit = readchomp(`git -C $dir rev-parse --short=10 HEAD`)
            isempty(readchomp(`git -C $dir status --porcelain --untracked-files=no`)) ||
                (commit *= "-dirty")
        catch
        end
    end
    return "v$ver $commit"
end

end # module Legacy
