# Contemporary groups (Phase 2): herds, group effects in phenotypes, groups
# fitted as fixed effects, the scenario keys and the archive columns.

@testset "Contemporary groups" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(11)
    vm = VariantMap(unique(sort(rand(rng, 1:30_000_000, 1500))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 300)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 2.0)
    q = eqtl(hap, 1:5:length(vm), t; rng = Xoshiro(3))
    key = SimKey(4)

    pop = founder_population(hap, q, t, key; herds = 8, c2 = 0.25)
    @test sort(unique(pop.herd)) == 1:8
    @test pop.σcg[1] ≈ sqrt(0.25 * 4.0 / 0.3)               # c² σ²a/h²
    # the group effect is added to the phenotype drawn without it
    p0 = founder_population(hap, q, t, key; herds = 8, c2 = 0.0)
    @test p0.herd == pop.herd && p0.tbv == pop.tbv
    cg = [BnGsim._cg(pop, i, 1) for i = 1:length(pop)]
    @test pop.phenotype[:, 1] ≈ p0.phenotype[:, 1] .+ cg
    # one effect per group (herd, step): equal within a herd of cohort 0
    @test all(cg[i] == cg[j] for i = 1:150, j = 1:150 if pop.herd[i] == pop.herd[j])
    @test length(unique(cg)) == 8
    @test founder_population(hap, q, t, key).herd == ones(Int32, 150)
    @test_throws ArgumentError founder_population(hap, q, t, key; herds = 0)
    @test_throws ArgumentError founder_population(hap, q, t, key; c2 = 1.0)

    # offspring are born in their dam's herd
    sc = Scheme(; noff = 120, nsire = 6, ndam = 30, evaluation = :pblup)
    simulate!(pop, sc, mm, 2)
    @test all(pop.herd[i] == pop.herd[pop.dam[i]] for i = 151:length(pop))

    # fixed groups: adding a constant to every record of one group leaves
    # the EBVs unchanged; with an overall mean only, it does not
    shifted(p, h, g, δ) = (c = deepcopy(p); m = (c.herd .== h) .& (c.generation .== g); c.phenotype[m, :] .+= δ; c)
    e1 = (evaluate!(pop, sc); copy(pop.ebv))
    ps = shifted(pop, pop.herd[1], 0, 5.0)
    evaluate!(ps, sc)
    @test ps.ebv ≈ e1
    q0 = founder_population(hap, q, t, key; herds = 8)
    simulate!(q0, sc, mm, 2)
    e0 = (evaluate!(q0, sc); copy(q0.ebv))
    qs = shifted(q0, q0.herd[1], 0, 5.0)
    evaluate!(qs, sc)
    @test !(qs.ebv ≈ e0)

    # two traits, the second recorded a step later: group of the record step
    ts = [Trait("a"; h² = 0.3), Trait("b"; h² = 0.4, age = 1)]
    mq = eqtl(hap, 1:5:length(vm), [1 0.3; 0.3 1], ts; rng = Xoshiro(5))
    p2 = founder_population(hap, mq, ts, key; herds = 5, c2 = [0.1, 0.3])
    @test p2.σcg ≈ sqrt.([0.1, 0.3] .* (1 ./ [0.3, 0.4]))
    @test BnGsim._groups(p2, [1, 2], 2) == BnGsim._groups(p2, [1, 2], 1)
    simulate!(p2, sc, mm, 2)
    e2 = (evaluate!(p2, sc); copy(p2.ebv))
    ps2 = shifted(p2, p2.herd[1], 0, 3.0)
    evaluate!(ps2, sc)
    @test ps2.ebv ≈ e2

    # scenario keys and archive columns
    d = JSON.parse(read(MINIMAL, String))
    d["founders"]["herds"] = 4
    d["traits"][1]["c2"] = 0.15
    s = parse_scenario(d)
    @test s.founders.herds == 4 && s.traits[1].c2 == 0.15
    @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    rep = run_replicate(s)
    @test maximum(rep.population.herd) <= 4 && rep.population.σcg[1] > 0
    dir = mktempdir()
    write_archive(dir, rep.population, rep.map, rep.species)
    @test readline(joinpath(dir, "pedigree.csv")) == "id,sire,dam,sex,generation,exit,herd"
    tl = readlines(joinpath(dir, "traits.csv"))
    @test endswith(tl[1], ",c2") && parse(Float64, split(tl[2], ',')[end]) ≈ 0.15
    for (k, v) in (("herds", 0), ("c2", 1.0))
        e = JSON.parse(read(MINIMAL, String))
        k == "herds" ? (e["founders"]["herds"] = v) : (e["traits"][1]["c2"] = v)
        @test_throws ScenarioError parse_scenario(e)
    end
end
