# The public reading of a Population (package roadmap §11, decision c) and
# the split of the 0.1–0.3 building blocks into BnGsim.Legacy (decision a).

@testset "Population accessors" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000])
    rng = Xoshiro(41)
    vm = VariantMap(unique(sort(rand(rng, 1:30_000_000, 800))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 400)))
    mm = MeiosisMap(sp, vm)
    ts = [Trait("y"; h² = 0.3), Trait("m"; h² = 0.2, sex = 0, age = 1)]
    mq = eqtl(hap, 1:5:length(vm), [1 0.2; 0.2 1], ts; rng = Xoshiro(2))
    sc = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup, max_sire_age = 2, max_dam_age = 2)

    # one record per trait
    pop = founder_population(hap, mq, ts, SimKey(3))
    simulate!(pop, sc, mm, 3)
    n = length(pop)
    @test tbv(pop) == pop.tbv && isequal(ebv(pop), pop.ebv) && isequal(phenotypes(pop), pop.phenotype)
    @test tbv(pop, [5, 2]) == pop.tbv[[5, 2], :] && size(ebv(pop, 1:3)) == (3, 2)
    @test_throws ArgumentError tbv(pop, [0])
    T = tbv(pop); T[1, 1] += 1
    @test pop.tbv[1, 1] != T[1, 1]                              # a copy
    @test [t.name for t in traits(pop)] == ["y", "m"] && trait_qtl(pop) == pop.qtl
    c = covariances(pop)
    @test c.G == pop.G0 && c.R == pop.R0 && c.PE == pop.PE0
    @test parity(pop) == pop.parity && all(iszero, parity(pop, findall(pop.male)))
    r = records(pop)
    obs = observed(pop)
    @test nrow(r) == count(obs)
    @test all(i -> obs[r.id[i], r.trait[i] == "y" ? 1 : 2], 1:nrow(r))
    @test r.step == [pop.generation[r.id[i]] + (r.trait[i] == "y" ? 0 : 1) for i = 1:nrow(r)]
    @test r.value == [pop.phenotype[r.id[i], r.trait[i] == "y" ? 1 : 2] for i = 1:nrow(r)]
    @test issorted(collect(zip(r.step, r.id)))
    @test nrow(records(pop, 0)) == count(observed(pop, 0))

    # repeated records: records() is the record list up to t
    pr = founder_population(hap, mq, ts, SimKey(3); records = [3, 1], repeatability = [0.5, nothing])
    simulate!(pr, sc, mm, 3)
    rr = records(pr)
    R = pr.records
    @test rr.id == R.id && rr.step == R.step && rr.value == R.y && rr.trait == ["y", "m"][R.trait]
    @test records(pr, 1) == rr[rr.step .<= 1, :]
    @test sum(==("y"), rr.trait) == sum(pr.nobs[:, 1])
end

@testset "Legacy module" begin
    legacy = [:sum_map, :gamete!, :gene_drop!, :random_mate, :rnd_ped, :expand, :sample_id, :norm_mv!, :tbv!,
        :phenotype!, :dummy_traits, :FounderPopulation, :sample_founder_population, :HistoricalSpec,
        :historical_population, :continue_population, :read_packed_haplotypes, :write_xy, :vanraden_g, :provenance]
    @test all(x -> x in names(BnGsim.Legacy), legacy)
    @test !any(x -> x in names(BnGsim), legacy)
    @test !any(x -> isdefined(BnGsim, x), (:sumMap, :sampleID, :sampleLoci, :gene_drop))
    # the old crossover sampler is Legacy's own function (0.38), not a method
    # of the exported crossovers, so it cannot be reached without Legacy
    sp = GenericSpecies("toy", [20_000_000, 10_000_000])
    vm = VariantMap(collect(UInt32, 1_000_000:1_000_000:29_000_000), sp)
    @test BnGsim.Legacy.crossovers !== crossovers
    @test !hasmethod(crossovers, (typeof(sp), typeof(vm)))
    @test hasmethod(BnGsim.Legacy.crossovers, (typeof(sp), typeof(vm)))
end
