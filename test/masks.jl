# Missing records (Phase 2): records missing at random and herds that do
# not record a trait, keyed so that the other records keep their values.

@testset "Missing records" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(41)
    vm = VariantMap(unique(sort(rand(rng, 1:30_000_000, 1500))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 1200)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3)
    q = eqtl(hap, 1:5:length(vm), t; rng = Xoshiro(3))
    key = SimKey(12)

    # at random: about 30 % missing, the others unchanged
    full = founder_population(hap, q, t, key; herds = 20)
    part = founder_population(hap, q, t, key; herds = 20, missing = 0.3)
    m = isnan.(part.phenotype[:, 1])
    @test abs(mean(m) - 0.3) < 4 * sqrt(0.21 / 600)
    @test part.phenotype[.!m, 1] == full.phenotype[.!m, 1]
    # whole herds: the pattern is constant within a herd
    hh = founder_population(hap, q, t, key; herds = 20, unrecorded_herds = 0.5)
    hm = isnan.(hh.phenotype[:, 1])
    @test all(h -> length(unique(hm[hh.herd .== h])) == 1, unique(hh.herd))
    @test 0 < count(h -> all(hm[hh.herd .== h]), unique(hh.herd)) < 20
    @test hh.phenotype[.!hm, 1] == full.phenotype[.!hm, 1]
    @test_throws ArgumentError founder_population(hap, q, t, key; missing = 1.0)

    # an evaluation uses the records there are; :phenotype refuses missing records
    sc = Scheme(; noff = 200, nsire = 10, ndam = 60, evaluation = :pblup)
    simulate!(part, sc, mm, 2)
    @test all(isfinite, part.ebv[part.generation .< current_step(part), :])      # the newest cohort is not evaluated yet
    @test_throws ArgumentError evaluate!(part, Scheme(; noff = 200, nsire = 10, ndam = 60, evaluation = :phenotype))

    # repeated records: missing records are skipped, the others unchanged,
    # the phenotype the mean of the records made
    rs = Scheme(; noff = 200, nsire = 10, ndam = 60, evaluation = :pblup, max_sire_age = 3, max_dam_age = 3)
    rf = founder_population(hap, q, t, key; records = 3, repeatability = 0.5)
    rp = founder_population(hap, q, t, key; records = 3, repeatability = 0.5, missing = 0.25)
    simulate!(rf, rs, mm, 3)
    simulate!(rp, rs, mm, 3)
    @test abs(length(rp.records.id) / length(rf.records.id) - 0.75) < 4 * sqrt(0.1875 / length(rf.records.id))
    # the founders' records (later cohorts differ: selection saw other records)
    made = Dict((rf.records.id[k], rf.records.step[k]) => rf.records.y[k] for k in eachindex(rf.records.id))
    fk = filter(k -> rp.records.id[k] <= 600, eachindex(rp.records.id))
    @test length(fk) > 1000 && all(made[(rp.records.id[k], rp.records.step[k])] == rp.records.y[k] for k in fk)
    for i in (1, 2, 650, 800)
        mine = findall(==(i), rp.records.id)
        @test isempty(mine) ? isnan(rp.phenotype[i, 1]) : rp.phenotype[i, 1] ≈ mean(rp.records.y[mine])
        @test rp.nobs[i, 1] == length(mine)
    end

    # scenario keys
    d = JSON.parse(read(MINIMAL, String))
    d["traits"][1]["missing"] = 0.2
    d["traits"][1]["unrecorded_herds"] = 0.1
    d["founders"]["herds"] = 5
    d["scheme"]["evaluation"] = "pblup"
    s = parse_scenario(d)
    @test s.traits[1].missing == 0.2 && s.traits[1].unrecorded_herds == 0.1
    @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    @test nrow(run_scenario(s)) == s.steps + 1
    d["traits"][1]["missing"] = 1.0
    @test_throws ScenarioError parse_scenario(d)
end
