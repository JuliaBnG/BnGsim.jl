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

sum_map(sp::Species) = sum_map(sp.chromosome; M = sp.M)

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

Base.@deprecate_binding sumMap sum_map

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

crossovers(sp::Species, vm::VariantMap; rng::AbstractRNG = Random.default_rng()) =
    crossovers(sum_map(sp), vm.pos; rng = rng)

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
    return gene_drop!(og.gt, pg.gt, pm, sum_map(sp), vm.pos; rng = rng)
end

"""
    gene_drop(pg, og, pm, lms, tbp; rng)

Deprecated: use `gene_drop!(og, pg, pm, lms, tbp; rng)`, which names the
mutation of `og` and takes the offspring first, as the `Haplotype`
method does.
"""
function gene_drop(
    pg::AbstractMatrix,
    og::AbstractMatrix,
    pm::AbstractMatrix{<:Integer},
    lms::DataFrame,
    tbp::AbstractVector{<:Integer};
    rng::AbstractRNG = Random.default_rng(),
)
    Base.depwarn(
        "`gene_drop(pg, og, …)` is deprecated; use `gene_drop!(og, pg, …)`.",
        :gene_drop,
    )
    return gene_drop!(og, pg, pm, lms, tbp; rng = rng)
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

    lms = sum_map(sp)
    p = (n_t / n0)^(1 / t)
    prt = copy(hap.gt)
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
    sub_gt = hap.gt[:, cols]
    frq = Float32.(vec(mean(sub_gt; dims = 2)))
    vld = findall(f -> maf < f < 1.0 - maf, frq)
    vm_out = VariantMap(vm.chr[vld], vm.pos[vld], vm.ref[vld], vm.alt[vld], frq[vld])
    return Haplotype(sub_gt[vld, :]), vm_out
end

Base.@deprecate_binding sampleID sample_id

"""
    sample_loci(
        hap::Haplotype,
        vm::VariantMap,
        sets::SNPSet...;
        rng::AbstractRNG = Random.default_rng(),
    ) -> Vector{LocusSet}

Sample one or more `BnGStructs.LocusSet` panels from `vm` according to the
specifications in `sets::SNPSet...` (locus count `nlc`, `maf` threshold, and `exclusive` flag).

Allele frequencies are taken from `vm.frq` if non-empty, or computed directly from `hap`.
"""
function sample_loci(
    hap::Haplotype,
    vm::VariantMap,
    sets::SNPSet...;
    rng::AbstractRNG = Random.default_rng(),
)
    nlc = hap.nlc
    nlc == length(vm) || throw(DimensionMismatch("Haplotype and VariantMap locus counts do not match"))
    names_seen = Set{String}()
    for s in sets
        s.name in names_seen && throw(ArgumentError("Duplicate SNPSet name: $(s.name)"))
        push!(names_seen, s.name)
    end

    frq = isempty(vm.frq) ? vec(mean(hap.gt; dims = 2)) : Float64.(vm.frq)
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

Base.@deprecate_binding sampleLoci sample_loci
