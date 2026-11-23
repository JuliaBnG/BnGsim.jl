# Aggregation over replicates and merging of runs.

@testset "Aggregate" begin
    df = DataFrame(replicate=[1, 1, 1, 2, 2, 2], generation=[0, 1, 2, 0, 1, 2],
        mean=[0.0, 1.0, 3.0, 0.0, 2.0, 2.0], F=[0.0, 0.1, 0.19, 0.0, 0.2, 0.2],
        accuracy=[0.5, 0.6, NaN, 0.7, 0.8, NaN])
    a = aggregate(df)
    @test a.generation == [0, 1, 2] && a.replicates == [2, 2, 2]
    @test a.mean_mean == [0.0, 1.5, 2.5] && a.mean_sd ≈ [0.0, sqrt(0.5), sqrt(0.5)]
    @test a.mean_se ≈ a.mean_sd ./ sqrt(2)
    # changes within each replicate: dG = 1, 2 then 2, 0; dF = 0.1, 0.1 then 0.2, 0
    @test isnan(a.dG_mean[1]) && a.dG_mean[2:3] ≈ [1.5, 1.0]
    @test a.dF_mean[2:3] ≈ [0.15, 0.05]
    # NaN left out: accuracy at the last step has no value
    @test a.accuracy_mean[1:2] ≈ [0.6, 0.7] && isnan(a.accuracy_mean[3]) && isnan(a.accuracy_sd[3])
    # one replicate: no spread
    @test all(isnan, aggregate(df[1:3, :]).mean_sd)
    # the order of the rows does not matter, to the bit
    for seed = 1:5
        p = shuffle(Xoshiro(seed), 1:nrow(df))
        @test isequal(aggregate(df[p, :]), a)
    end
    @test_throws ArgumentError aggregate(DataFrame(x=[1]))
    # variance depletion relative to each replicate's first step, and the
    # rate of loss of heterozygosity, shuffled rows included
    dv = DataFrame(replicate=[2, 1, 1, 2], generation=[1, 0, 1, 0], var=[1.0, 2.0, 1.0, 4.0],
        genic=[3.0, 2.0, 1.5, 4.0], Fhet=[0.5, 0.0, 0.1, 0.2])
    av = aggregate(dv)
    @test av.var_ratio_mean == [1.0, (0.5 + 0.25) / 2] && av.genic_ratio_mean == [1.0, (0.75 + 0.75) / 2]
    @test isnan(av.dFhet_mean[1]) && av.dFhet_mean[2] ≈ (0.1 + 0.3 / 0.8) / 2
    # several traits: dG for the index and for each trait
    t = parse_scenario(scenario_doc())
    s2 = run_scenario(t)
    @test all(c -> c in names(aggregate(s2)), ["dG_mean", "dG_a_mean", "dG_b_mean", "dF_mean", "rg_a_b_mean"])
end

@testset "CLI: aggregate.csv and merge" begin
    three = tempname() * ".json"
    write(three, replace(read(MINIMAL, String), "\"steps\": 3," => "\"steps\": 4, \"replicates\": 3,"))
    dir = mktempdir()
    all3, r12, r3 = joinpath(dir, "all"), joinpath(dir, "r12"), joinpath(dir, "r3")
    @test run_cli("run", three, "--out", all3)[1] == 0
    @test isfile(joinpath(all3, "aggregate.csv"))
    p = JSON.parse(read(joinpath(all3, "provenance.json"), String))
    @test haskey(p["outputs"], "aggregate.csv")
    @test isequal(BnGsim._read_summary(joinpath(all3, "summary.csv")).mean, run_scenario(read_scenario(three)).mean)
    @test run_cli("run", three, "--out", r12, "--replicates", "1:2")[1] == 0
    @test run_cli("run", three, "--out", r3, "--replicates", "3")[1] == 0
    # merged in either order: the bytes of one run of all three
    for (k, order) in enumerate(((r12, r3), (r3, r12)))
        m = joinpath(dir, "m$k")
        code, v, _ = run_cli("merge", order..., "--out", m)
        @test code == 0 && v["replicates"] == [1, 2, 3]
        for f in ("scenario.json", "summary.csv", "aggregate.csv")
            @test read(joinpath(m, f)) == read(joinpath(all3, f))
        end
        q = JSON.parse(read(joinpath(m, "provenance.json"), String))
        @test length(q["sources"]) == 2 && q["scenario_sha256"] == p["scenario_sha256"]
    end
    # refused: shared replicates, another scenario, an unfinished run, a
    # non-empty output, no --out
    @test run_cli("merge", r12, all3, "--out", joinpath(dir, "x1"))[1] == 4
    other = joinpath(dir, "other")
    @test run_cli("run", MINIMAL, "--out", other)[1] == 0
    @test run_cli("merge", r3, other, "--out", joinpath(dir, "x2"))[1] == 4
    mkpath(joinpath(dir, "half", ".bngsim-run"))
    @test run_cli("merge", r12, joinpath(dir, "half"), "--out", joinpath(dir, "x3"))[1] == 4
    @test run_cli("merge", r12, r3, "--out", all3)[1] == 4
    @test run_cli("merge", r12, r3)[1] == 2
    rm(dir; recursive=true)
end
