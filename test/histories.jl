# Several historical populations per run, stored histories, and the
# between-history components of the aggregation; QTL drawn from and
# correlated with a historical trait.

@testset "Histories: several, stored, between-history variance" begin
    # replicates map to histories 1, 2, 3, 1, ...; history 1 is the single
    # shared history of a scenario without `histories`
    d = JSON.parse(read(MINIMAL, String))
    d["founders"]["histories"] = 3
    d["replicates"] = 6
    s = parse_scenario(d)
    @test [BnGsim._history_index(s, r) for r = 1:6] == [1, 2, 3, 1, 2, 3]
    df = run_scenario(s)
    @test all(df.history .== (df.replicate .- 1) .% 3 .+ 1)
    s1 = parse_scenario(JSON.parse(read(MINIMAL, String)))
    @test BnGsim._history_key(s, 1) == BnGsim._history_key(s1, 1)
    @test isequal(df[df.replicate .<= 1, Not(:history)], run_scenario(s1)[:, Not(:history)])
    @test !isequal(df[df.replicate .== 2, 3:end], df[df.replicate .== 4, 3:end] )
    # without shared histories, history = replicate
    d2 = JSON.parse(read(MINIMAL, String))
    d2["founders"]["shared_history"] = false
    d2["replicates"] = 2
    @test run_scenario(parse_scenario(d2)).history == run_scenario(parse_scenario(d2)).replicate
    d2["founders"]["histories"] = 2
    @test_throws ScenarioError parse_scenario(d2)

    # between-history components, against hand values: 2 histories × 2
    # replicates at one step
    t = DataFrame(replicate=[1, 2, 3, 4], history=[1, 2, 1, 2], generation=zeros(Int, 4), x=[1.0, 5.0, 3.0, 7.0])
    a = aggregate(t)
    @test a.histories == [2] && a.x_mean == [4.0]
    @test a.x_sdh ≈ [std([2.0, 6.0])] && a.x_seh ≈ [std([2.0, 6.0]) / sqrt(2)]
    @test a.x_se ≈ [std([1.0, 5.0, 3.0, 7.0]) / 2]
    @test !("x_sdh" in names(aggregate(t[:, Not(:history)])))
    @test !("x_sdh" in names(aggregate(DataFrame(replicate=1:3, history=1:3, generation=zeros(Int, 3), x=[1.0, 2, 4]))))
    @test isequal(aggregate(t[[4, 2, 3, 1], :]), a)

    # stored histories: written once, reused with another seed, the same
    # simulation as with the histories simulated in the run
    dir = mktempdir()
    hdir = joinpath(dir, "hist")
    sc = joinpath(dir, "s.json")
    write(sc, JSON.json(d))
    @test run_cli("history", sc, "--out", hdir)[1] == 0
    @test sort(readdir(hdir)) == ["history-1.fwr", "history-2.fwr", "history-3.fwr", "manifest.json"]
    m = JSON.parse(read(joinpath(hdir, "manifest.json"), String))
    @test m["schema"] == "bngsim-history" && length(m["histories"]) == 3
    d["founders"]["history_dir"] = hdir
    st = parse_scenario(d)
    @test check_histories(st) === nothing
    @test isequal(run_scenario(st), df)
    d["seed"] = 77
    @test check_histories(parse_scenario(d)) === nothing                 # the seed is not a history parameter
    out = joinpath(dir, "run")
    write(sc, JSON.json(d))
    @test run_cli("run", sc, "--out", out)[1] == 0
    p = JSON.parse(read(joinpath(out, "provenance.json"), String))
    @test p["history"]["manifest_sha256"] == BnGsim._sha256_file(joinpath(hdir, "manifest.json"))
    @test "x_sdh" ∉ names(aggregate(BnGsim._read_summary(joinpath(out, "summary.csv")))) &&
          "mean_seh" in names(aggregate(BnGsim._read_summary(joinpath(out, "summary.csv"))))
    # mismatches are refused: other parameters, another number, a changed file
    bad = deepcopy(d)
    bad["founders"]["nt"] = 499
    @test occursin("other historical parameters", check_histories(parse_scenario(bad)))
    write(sc, JSON.json(bad))
    @test run_cli("validate", sc)[1] == 3
    bad = deepcopy(d)
    bad["founders"]["histories"] = 2
    @test occursin("holds 3 histories", check_histories(parse_scenario(bad)))
    open(joinpath(hdir, "history-2.fwr"), "a") do io
        write(io, 0x00)
    end
    @test occursin("SHA-256", check_histories(parse_scenario(d)))
    clean = deepcopy(d)
    delete!(clean["founders"], "history_dir")
    write(sc, JSON.json(clean))
    @test run_cli("history", sc, "--out", hdir)[1] == 4                  # never overwrites
    rm(dir; recursive=true)
end

@testset "QTL from history: subset and correlated effects" begin
    d = _history_doc()
    d["loci"]["qtl"] = Dict("from" => "history", "maf" => 0.0, "nlc" => 10, "rho" => 0.6)
    s = parse_scenario(d)
    @test s.qtl.nlc == 10 && s.qtl.rho == 0.6
    @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    rep = BnGsim._setup_replicate(s, 1, nothing)
    @test length(rep.qtl) == 10
    h = BnGsim._history_of(s, 1)
    full = BnGsim._setup_replicate(parse_scenario(_history_doc()), 1, nothing)
    @test issubset(rep.qtl, full.qtl)                                    # a subset of the segregating trait loci
    # rho = 1 keeps the historical effects; rho < 1 mixes in new ones
    @test BnGsim._breeding_effects([1.0, -2.0, 3.0], 1.0, Trait("y"), Xoshiro(1)) == [1.0, -2.0, 3.0]
    a = randn(Xoshiro(2), 20_000)
    b = BnGsim._breeding_effects(a, 0.6, Trait("y"), Xoshiro(3))
    @test abs(cor(a, b) - 0.6) < 4 * (1 - 0.36) / sqrt(20_000)
    d["loci"]["qtl"]["nlc"] = 100_000
    @test_throws ArgumentError BnGsim._setup_replicate(parse_scenario(d), 1, nothing)
    d["loci"]["qtl"] = Dict("nlc" => 10, "rho" => 0.5)
    @test_throws ScenarioError parse_scenario(d)                         # rho with neutral QTL
end
