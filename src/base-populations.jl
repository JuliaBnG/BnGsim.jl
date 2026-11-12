using BnGStructs
using FisherWright
using RelationshipMatrices
using DataFrames
using Distributions
using Random
using Statistics

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

_to_variant_map(vm::VariantMap) = vm
_to_variant_map(df::DataFrame) = VariantMap(
    Int8.(df.chr),
    UInt32.(df.pos),
    Char.(df.ref),
    Char.(df.alt),
    hasproperty(df, :frq) ? Float32.(df.frq) : Float32[],
)

"""
    sample_founder_population(
        species::Species;
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
by running a Fisher-Wright burn-in of `burnin_generations` using the specified `species` parameters.
"""
function sample_founder_population(
    species::Species;
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
        burnin_generations,
        mutation_rate;
        verbose = verbose,
        result = true,
        rng = rng,
    )

    bit_gt, map_df = muts2bitarray(res.active_haplotypes, res.chromosome_ends)
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
