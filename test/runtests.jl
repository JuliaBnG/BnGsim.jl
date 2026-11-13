using Test
using Aqua
using BnGsim
using BnGStructs
using Distributions
using DataFrames
using Random

@testset "Aqua" begin
    Aqua.test_all(BnGsim)
end

@testset "BnGsim Founder Population Tests" begin
    # Fast test with small species
    sp = Cattle(50) # small Ne for fast test
    
    @testset "Sample from fresh Fisher-Wright simulation" begin
        founder = sample_founder_population(
            sp;
            n_id = 10,
            n_marker = 50,
            n_qtl = 10,
            h2 = 0.4,
            burnin_generations = 10,
            rng = MersenneTwister(42),
        )

        @test founder isa FounderPopulation
        @test founder.haplotype.nhp == 20
        @test founder.map isa VariantMap
        @test founder.markers isa LocusSet
        @test length(founder.markers) == 50
        @test length(founder.marker_indices) == 50
        @test length(founder.qtl) == 10
        @test founder.h2 == 0.4
        @test length(founder.tbv) == 10
        @test length(founder.phenotypes) == 10
        @test size(founder.pedigree, 1) == 10
        @test all(founder.pedigree.sire .== 0)
        @test all(founder.pedigree.dam .== 0)

        # Markers and QTLs should not overlap
        @test isempty(intersect(founder.marker_indices, founder.qtl.loci))
    end

    @testset "Sample from pre-existing Haplotype and VariantMap" begin
        # Synthetic Haplotype (100 loci, 20 haps = 10 individuals)
        gt = BitMatrix(rand(MersenneTwister(11), Bool, 100, 20))
        hap = Haplotype(gt)
        vm = VariantMap(
            fill(Int8(1), 100),
            UInt32.(1000:1000:100000),
            fill('A', 100),
            fill('G', 100),
            fill(0.5f0, 100),
        )

        founder = sample_founder_population(
            sp,
            hap,
            vm;
            n_id = 5,
            n_marker = 40,
            n_qtl = 10,
            h2 = 0.6,
            rng = MersenneTwister(123),
        )

        @test founder.haplotype.nhp == 10
        @test length(founder.tbv) == 5
        @test length(founder.phenotypes) == 5
        @test length(founder.markers) == 40
        @test length(founder.qtl) == 10
    end

    @testset "BnGStructs-native reproduction, QTL, TBV, and phenotypes" begin
        rng = MersenneTwister(99)
        sp_toy = GenericSpecies("toy", 20, [1_000_000, 1_000_000])
        nlc, nid = 128, 20
        gt = BitMatrix(rand(rng, Bool, nlc, 2nid))
        # Ensure all loci are polymorphic
        gt[:, 1] .= true
        gt[:, 2] .= false
        hap0 = Haplotype(gt)
        frq0 = Float32.(vec(mean(gt; dims = 2)))
        vm0 = VariantMap(
            vcat(fill(Int8(1), 64), fill(Int8(2), 64)),
            UInt32.(round.(Int, range(1_000, 1_900_000; length = nlc))),
            fill('A', nlc),
            fill('C', nlc),
            frq0,
        )

        # sum_map
        lms_sp = sum_map(sp_toy)
        @test nrow(lms_sp) == 2
        lms_vm = sum_map(vm0)
        @test nrow(lms_vm) == 2 && sum(lms_vm.nlc) == nlc

        # sample_id and sample_loci
        hap_sub, vm_sub = sample_id(hap0, vm0, 12; maf = 0.01, rng = MersenneTwister(101))
        @test hap_sub.nhp == 24 && hap_sub.nlc == length(vm_sub)
        chip_set = SNPSet("chip", 40; maf = 0.05, exclusive = true)
        qtl_set = SNPSet("qtl", 20; maf = 0.01, exclusive = true)
        panels = sample_loci(hap0, vm0, chip_set, qtl_set; rng = MersenneTwister(102))
        @test length(panels) == 2
        @test length(panels[1]) == 40 && length(panels[2]) == 20
        @test isempty(intersect(panels[1].loci, panels[2].loci))

        # gene_drop! and expand
        pm, sex0 = random_mate(nid, 16; rng = MersenneTwister(103))
        @test size(pm) == (16, 2) && length(sex0) == nid
        og1 = Haplotype(nlc, 32)
        og2 = Haplotype(nlc, 32)
        df_cp1 = gene_drop!(og1, hap0, pm, sp_toy, vm0; rng = MersenneTwister(104))
        df_cp2 = gene_drop!(og2, hap0, pm, sp_toy, vm0; rng = MersenneTwister(104))
        @test og1.gt == og2.gt
        @test nrow(df_cp1) == 32

        hap_exp, vm_exp, cps_exp = expand(hap0, vm0, sp_toy, 30, 2; maf = 0.01, rng = MersenneTwister(105))
        @test hap_exp.nhp == 60 && hap_exp.nlc == length(vm_exp)
        @test nrow(cps_exp) > 0

        # rnd_ped
        ped_r = rnd_ped(10, 3; rng = MersenneTwister(106))
        @test nrow(ped_r) == 40

        # Single-trait and multi-trait eqtl, tbv!, phenotype!, dummy_traits
        t1 = Trait("milk"; h² = 0.3, sex = 0, μ = 10.0, σₐ = 2.0)
        t2 = Trait("growth"; h² = 0.5, sex = 2, μ = 50.0, σₐ = 3.0)
        t_cat = Trait("ce", [0.7, 0.2, 0.1]; h² = 0.2, sex = 2)

        qtl1 = eqtl(hap0, panels[2], t1; rng = MersenneTwister(107))
        @test qtl1 isa TraitQTL
        @test std(tbv(hap0, qtl1)) ≈ 2.0 atol = 1e-4

        R_g = [1.0 0.4 0.1; 0.4 1.0 -0.2; 0.1 -0.2 1.0]
        mqtl = eqtl(hap0, panels[2], R_g, [t1, t2, t_cat]; rng = MersenneTwister(108))
        @test mqtl isa MultiTraitQTL
        ped_t = DataFrame(
            id = 1:nid,
            sire = zeros(Int32, nid),
            dam = zeros(Int32, nid),
            sex = Int8.([i <= 10 ? 0 : 1 for i = 1:nid]),
        )
        tbv!(ped_t, hap0, mqtl)
        @test std(ped_t.tbv_milk) ≈ 2.0 atol = 1e-4
        @test std(ped_t.tbv_growth) ≈ 3.0 atol = 1e-4
        @test std(ped_t.tbv_ce) ≈ 1.0 atol = 1e-4

        phenotype!(ped_t, [t1, t2, t_cat]; rng = MersenneTwister(109))
        # milk is female-only (sex == 0), so males (rows 11:20) must be missing
        @test all(!ismissing, ped_t.ft_milk[1:10])
        @test all(ismissing, ped_t.ft_milk[11:20])
        @test all(!ismissing, ped_t.ft_growth)
        @test all(x -> x in (0, 1, 2), ped_t.ft_ce)

        d_trts, d_R, d_itts = dummy_traits(3, 100; rng = MersenneTwister(110))
        @test length(d_trts) == 4 && size(d_R) == (4, 4) && length(d_itts) == 3
        @test_throws ErrorException dummy_traits(0, 100)
        @test_throws ErrorException dummy_traits(6, 100)
    end
end

@testset "Founder design: historical populations and summaries" begin
    using FisherWright, Statistics
    sp = GenericSpecies("toy", 100, [5_000_000, 5_000_000])
    spec = HistoricalSpec("toy", sp, 300, 5.0)
    res = historical_population(spec; seed = 1)
    @test res isa FisherWrightResult
    @test length(res.active_haplotypes) == 200

    # Decline to 20 then 10 individuals through `fisher_wright(; init, popsize)`.
    rng = MersenneTwister(2)
    dec = continue_population(res, [20, 10], 5.0; M = spec.M, rng = rng)
    @test length(dec.active_haplotypes) == 20
    @test all(issorted, dec.active_haplotypes)
    @test all(h -> allunique(h), dec.active_haplotypes)
    @test dec.chromosome_ends == res.chromosome_ends
    # Same rng state gives the same population.
    dec2 = continue_population(res, [20, 10], 5.0; M = spec.M, rng = MersenneTwister(2))
    @test dec2.active_haplotypes == dec.active_haplotypes
    @test dec.census[end-2:end] == [spec.ne, 20, 10]
    @test continue_population(res, Int[], 5.0) === res
    @test_throws ArgumentError continue_population(res, [20, 1], 5.0)
    # A seed and the stream it starts are the same request, and the stream of
    # the global RNG is not consumed.
    Random.seed!(99)
    next_draw = rand(UInt64)
    Random.seed!(99)
    h1 = historical_population(spec; seed = 11)
    h2 = historical_population(spec; rng = Xoshiro(11))
    @test h1.active_haplotypes == h2.active_haplotypes
    @test rand(UInt64) == next_draw
    @test_throws ArgumentError historical_population(spec)
    @test_throws ArgumentError historical_population(spec; seed = 1, rng = Xoshiro(1))
    # vanraden_g against the definition
    let xyg = rand(MersenneTwister(8), Bool, 120, 30), lc = collect(1:2:120)
        Zg = [xyg[l, 2i-1] + xyg[l, 2i] for i = 1:15, l in lc] .+ 0.0
        pg = vec(sum(Zg; dims = 1)) ./ 30
        Zg .-= 2pg'
        @test vanraden_g(xyg, lc) ≈ (Zg * Zg') ./ (2sum(pg .* (1 .- pg)))
    end

    xy, lmp = founder_haplotypes(res, 50; rng = MersenneTwister(3))
    @test size(xy, 2) == 100
    @test size(xy, 1) == nrow(lmp)
    @test all(0 .< lmp.frq .< 1)

    s = founder_summary(xy, lmp; nchip = 50, nqtl = 100, rng = MersenneTwister(4),
        ld = (npairs = 500, edges = [1e4, 1e5, 1e6]))
    @test s.scalars.nid == 50
    @test isempty(intersect(s.chip, s.qtl))
    @test sum(s.sfs.n) == s.scalars.nseg
    @test 0 < s.scalars.chip_r2_adj < 1
    @test s.scalars.g_diag ≈ 1 atol = 0.2
    @test issorted(s.ld.r2; rev = true)   # LD decays with distance
    # QTL come from what the presample × nchip chip candidates leave
    maf = min.(lmp.frq, 1 .- lmp.frq)
    ncom = count(>(0.2), maf)
    @test count(>(0.2), maf[s.qtl]) ≤ ncom - 3 * 50
    k = ncom ÷ 50
    sk = founder_summary(xy, lmp; nchip = 50, nqtl = 100, presample = k, rng = MersenneTwister(4),
        ld = (npairs = 500, edges = [1e4, 1e5, 1e6]))
    @test count(>(0.2), maf[sk.qtl]) ≤ ncom - k * 50 < 50
    @test_throws ErrorException founder_summary(xy, lmp; nchip = 50, nqtl = 100, presample = k + 1)

    # r² against a direct computation
    W = BnGsim.pack_loci(xy, [1, 2])
    a, b = xy[1, :], xy[2, :]
    @test BnGsim.r2(W, 1, 2, 100) ≈ cor(a, b)^2
end

@testset "Founder design: xyBnG base reader" begin
    using Serialization
    dir = mktempdir()
    nr, nc = 7, 8
    bits = BitMatrix(rand(MersenneTwister(5), Bool, nr, nc))
    bits[1, :] .= true                       # monomorphic, must be dropped
    open(joinpath(dir, "b.xy"), "w") do io
        write(io, Int8['x', 'y', '-', 0, 0, 0, 0, 0])
        write(io, [nr, nc])
        write(io, bits)
    end
    pos = UInt32[10, 30, 20, 150, 40, 120, 160]   # unsorted, as from v0.1.x
    serialize(joinpath(dir, "b.lmp"), DataFrame(chr = Int8[1, 1, 1, 2, 1, 2, 1],
        pos = pos, ref = fill('A', nr), alt = fill('C', nr), frq = zeros(nr)))
    xy, lmp = founder_haplotypes(joinpath(dir, "b.xy"), joinpath(dir, "b.lmp"), 4, [100, 200];
        rng = MersenneTwister(6))
    @test issorted(lmp.pos)
    @test 10 ∉ lmp.pos
    @test lmp.chr == Int8.(ifelse.(lmp.pos .≤ 100, 1, 2))
    for (k, p) in enumerate(lmp.pos)
        @test xy[k, :] == bits[findfirst(==(p), pos), :]
    end
    @test lmp.frq ≈ vec(sum(xy; dims = 2)) ./ 8
end

@testset "Founder design: BNGH1 reader" begin
    nlc, nhp = 5, 10
    bits = BitMatrix(rand(MersenneTwister(7), Bool, nlc, nhp))
    bits[3, :] .= false                      # monomorphic, must be dropped
    pos = UInt32[5, 40, 60, 150, 190]
    path = joinpath(mktempdir(), "h.bin")
    open(path, "w") do io
        write(io, UInt8['B', 'N', 'G', 'H', '1', 0, 0, 0])
        write(io, [nlc, nhp])
        write(io, pos)
        for l = 1:nlc, b = 1:cld(nhp, 8)
            byte = 0x00
            for s = 0:7
                h = 8(b - 1) + s + 1
                h ≤ nhp && bits[l, h] && (byte |= UInt8(1) << s)
            end
            write(io, byte)
        end
    end
    xy, lmp = read_packed_haplotypes(path, [100, 200])
    keep = [k for k = 1:nlc if 0 < sum(bits[k, :]) < nhp]
    @test lmp.pos == pos[keep]
    @test xy == bits[keep, :]
    @test lmp.chr == Int8.(ifelse.(lmp.pos .≤ 100, 1, 2))
    @test all(lmp.ref .!= lmp.alt)
end

@testset "Founder design: xy writer round trip" begin
    dir = mktempdir()
    bits = BitMatrix(rand(MersenneTwister(8), Bool, 9, 6))
    bits[:, 1] .= .!bits[:, 2]               # no monomorphic rows
    write_xy(joinpath(dir, "r.xy"), bits)
    serialize(joinpath(dir, "r.lmp"), DataFrame(chr = Int8.(fill(1, 9)), pos = UInt32.(1:9),
        ref = fill('A', 9), alt = fill('C', 9), frq = zeros(9)))
    xy, lmp = founder_haplotypes(joinpath(dir, "r.xy"), joinpath(dir, "r.lmp"), 3, [100])
    @test xy == bits
end

@testset "provenance header" begin
    io = IOBuffer()
    provenance(io; seed = 42, extra = Dict("arm" => "test"))
    txt = String(take!(io))
    @test occursin("# seed:     42", txt)
    @test occursin("BnGsim", txt) && occursin("FisherWright", txt)
    @test occursin("# arm:", txt)
    @test all(startswith("# "), filter(!isempty, split(txt, '\n')))
end
