# Repeated records (Phase 2): records per step from the trait's age, the
# permanent environment, the repeatability model against a dense solve of
# Henderson's equations, the scenario keys and the archive.
using RelationshipMatrices: ainv

@testset "Repeated records" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(21)
    vm = VariantMap(unique(sort(rand(rng, 1:30_000_000, 1500))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 200)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.2, σₐ = 1.0)
    q = eqtl(hap, 1:5:length(vm), t; rng = Xoshiro(3))
    pop = founder_population(hap, q, t, SimKey(8); records = 3, repeatability = 0.5, herds = 3, c2 = 0.1)
    Vp = 5.0
    @test pop.σpe[1] ≈ sqrt(0.3Vp) && pop.R0[1, 1] ≈ 0.5Vp && pop.PE0[1] ≈ 0.3Vp && pop.σcg[1] ≈ sqrt(0.1Vp)
    @test all(isnan, pop.phenotype) && !all(iszero, pop.pe)
    sc = Scheme(; noff = 80, nsire = 6, ndam = 25, evaluation = :pblup, max_sire_age = 3, max_dam_age = 3,
        survival = [1.0, 0.9])
    steps = 6
    simulate!(pop, sc, mm, steps)
    R = pop.records
    # one record per step alive from the trait's age, at most 3; the
    # phenotype is the mean of the records
    for i = 1:length(pop)
        mine = findall(==(i), R.id)
        first_ = pop.generation[i] + Int(t.age)
        last_ = min(pop.exit[i] - 1, steps - 1, first_ + 2)
        @test R.step[mine] == collect(first_:last_)
        @test isempty(mine) ? isnan(pop.phenotype[i, 1]) : pop.phenotype[i, 1] ≈ mean(R.y[mine])
    end

    # the repeatability model equals Henderson's equations solved densely
    m = evaluate!(pop, sc)
    n = length(pop)
    y = R.y
    groups = BnGsim._groups(pop, Int.(R.id), 1; steps = R.step)
    X = Float64.([g == c for g in groups, c = 1:maximum(groups)])
    Z = Float64.([R.id[k] == i for k in eachindex(y), i = 1:n])
    λa, λp = pop.R0[1, 1] / pop.G0[1, 1], pop.R0[1, 1] / pop.PE0[1]
    Ai = Matrix(ainv(pop.sire, pop.dam))
    C = [X'X X'Z X'Z; Z'X Z'Z+λa*Ai Z'Z; Z'X Z'Z Z'Z+λp*I]
    sol = C \ [X'y; Z'y; Z'y]
    p = size(X, 2)
    @test pop.ebv[:, 1] ≈ sol[p+1:p+n]

    # the same population without repeated records draws no records and no pe
    p1 = founder_population(hap, q, t, SimKey(8))
    @test p1.nrec == [1] && all(iszero, p1.pe) && isempty(p1.records.id)
    @test_throws ArgumentError founder_population(hap, q, t, SimKey(8); repeatability = 0.5)
    @test_throws ArgumentError founder_population(hap, q, t, SimKey(8); records = 2, repeatability = 0.1)
    @test_throws ArgumentError founder_population(hap, q, t, SimKey(8); records = 0)

    # scenario keys and the archive
    d = JSON.parse(read(MINIMAL, String))
    d["traits"][1]["records"] = 3
    d["traits"][1]["repeatability"] = 0.6
    d["scheme"]["max_dam_age"] = 2
    s = parse_scenario(d)
    @test s.traits[1].records == 3 && s.traits[1].repeatability == 0.6
    @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    rep = run_replicate(s)
    @test "vpe" in names(rep.summary)
    dir = mktempdir()
    write_archive(dir, rep.population, rep.map, rep.species)
    rl = readlines(joinpath(dir, "repeated.csv"))
    @test rl[1] == "id,trait,step,value" && length(rl) == length(rep.population.records.id) + 1
    tl = split.(readlines(joinpath(dir, "traits.csv")), ',')
    @test tl[1][end-1:end] == ["records", "repeatability"] && tl[2][end-1] == "3" && parse(Float64, tl[2][end]) ≈ 0.6
    d2 = JSON.parse(read(MINIMAL, String))
    d2["traits"][1]["repeatability"] = 0.6
    @test_throws ScenarioError parse_scenario(d2)               # needs records > 1
end
