"""
    sample_loci(
        hap::Haplotype,
        vm::VariantMap,
        sets::SNPSet...;
        rng::AbstractRNG,
    ) -> Vector{LocusSet}

Sample one or more `BnGStructs.LocusSet` panels from `vm` according to the
specifications in `sets::SNPSet...` (locus count `nlc`, `maf` threshold, and `exclusive` flag).

Allele frequencies are taken from `vm.frq` if non-empty, or computed directly from `hap`.

`rng` is required (0.38): pass a keyed generator, `keyed_rng(SimKey(…))`, to
draw inside the keyed streams of the run; a scenario does so.
"""
function sample_loci(
    hap::Haplotype,
    vm::VariantMap,
    sets::SNPSet...;
    rng::AbstractRNG,
)
    nlc = hap.nlc
    nlc == length(vm) || throw(DimensionMismatch("Haplotype and VariantMap locus counts do not match"))
    names_seen = Set{String}()
    for s in sets
        s.name in names_seen && throw(ArgumentError("Duplicate SNPSet name: $(s.name)"))
        push!(names_seen, s.name)
    end

    frq = isempty(vm.frq) ? vec(mean(view(hap.gt, 1:nlc, :); dims = 2)) : Float64.(vm.frq)
    all_loci = collect(Int32(1):Int32(nlc))
    used = Int32[]
    out = Vector{LocusSet}(undef, length(sets))

    for (k, s) in enumerate(sets)
        vld = findall(f -> s.maf < f < 1.0 - s.maf, frq)
        pool = all_loci[vld]
        s.exclusive && setdiff!(pool, used)
        length(pool) >= s.nlc ||
            throw(ArgumentError("Not enough qualifying loci ($(length(pool))) to sample $(s.nlc) for SNPSet \"$(s.name)\""))
        chosen = sort(sample(rng, pool, s.nlc; replace = false))
        append!(used, chosen)
        out[k] = LocusSet(s.name, chosen)
    end
    return out
end
