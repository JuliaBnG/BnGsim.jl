# The resource estimate: the sizes it is computed from, its structure, the
# Phase 1 gate against the measured calibration runs, and the CLI.

@testset "Estimate" begin
    s = read_scenario(MINIMAL)
    z = ScenarioSize(s)
    μL = 1e-8 * 30_000_000
    @test z.alleles ≈ μL * min(500, 400) && z.animals == 50 + 3 * 50 && z.evaluated == 50 + 2 * 50
    @test z.loci ≈ 4 * 100 * μL * log1p(2 * 50 * 500 / 400) && z.genotyped == 0 && z.cohort == 50
    e = resource_estimate(s)
    @test e.peak_rss_bytes == sum(values(e.memory_terms)) && e.memory_terms.dense == 0
    @test e.peak_rss_bound_bytes ≈ 1.3 * e.peak_rss_bytes atol = 1
    @test e.disk_bytes > 0 && e.seconds >= 0
    # GBLUP: every evaluated animal genotyped; ssGBLUP with genotyped_cohorts.
    t = parse_scenario(scenario_doc())
    @test ScenarioSize(t).genotyped == ScenarioSize(t).evaluated && resource_estimate(t).memory_terms.dense > 0
    d = scenario_doc()
    d["scheme"]["evaluation"] = "ssgblup"
    d["scheme"]["genotyped_cohorts"] = 1
    @test ScenarioSize(parse_scenario(d)).genotyped == 40
    # More offspring, more memory.
    big = parse_scenario(edited(["scheme", "noff"], 400))
    @test resource_estimate(big).peak_rss_bytes > resource_estimate(t).peak_rss_bytes

    # The Phase 1 gate on the measured runs (validation/estimate): the mean
    # peak of every scenario within ±25 % of the estimate, and no run above
    # the bound.
    cal = joinpath(@__DIR__, "..", "validation", "estimate")
    runs = Dict{String,Vector{Float64}}()
    for l in eachline(joinpath(cal, "measured.tsv"))
        occursin(".json", l) || continue
        f = split(l, '\t')
        push!(get!(runs, f[1], Float64[]), parse(Float64, f[2]) * 2^20)
    end
    @test length(runs) == 22
    for (f, m) in runs
        est = resource_estimate(read_scenario(joinpath(cal, "scenarios", f)))
        @test abs(mean(m) / est.peak_rss_bytes - 1) <= 0.25
        @test maximum(m) <= est.peak_rss_bound_bytes
    end

    # The terms added after the calibration (0.33): the same gate on their
    # scenarios (several histories, a growing history, repeated records,
    # reliability)
    later = joinpath(cal, "later")
    lr = Dict{String,Vector{Float64}}()
    for l in eachline(joinpath(later, "measured.tsv"))
        occursin(".json", l) || continue
        f = split(l, '\t')
        push!(get!(lr, f[1], Float64[]), parse(Float64, f[2]) * 2^20)
    end
    @test length(lr) == 4
    for (f, m) in lr
        est = resource_estimate(read_scenario(joinpath(later, f)))
        @test abs(mean(m) / est.peak_rss_bytes - 1) <= 0.25
        @test maximum(m) <= est.peak_rss_bound_bytes
    end

    # Command line
    code, v, _ = run_cli("estimate", MINIMAL)
    @test code == 0 && v["status"] == "estimated" && v["peak_rss_bytes"] == e.peak_rss_bytes
    @test v["peak_rss_bound_bytes"] == e.peak_rss_bound_bytes && v["size"]["animals"] == 200
    @test run_cli("estimate")[1] == 2
    dir = mktempdir()
    code, v, _ = run_cli("run", MINIMAL, "--out", joinpath(dir, "a"), "--memory-limit", "100M")
    @test code == 5 && v["status"] == "rejected" && v["memory_limit_bytes"] == 100 * 2^20 && !ispath(joinpath(dir, "a"))
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "b"), "--memory-limit", "64G")[1] == 0
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "c"), "--memory-limit", "lots")[1] == 2
    @test BnGsim._parse_size("512M") == 512 * 2^20 && BnGsim._parse_size("1.5G") == round(Int, 1.5 * 2^30)
    rm(dir; recursive=true)
end
