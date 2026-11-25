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

    # The time model (0.37): the parts add up to the total; the features sum
    # over the steps g = 0, 1 with nid + g·noff animals (GBLUP: all
    # genotyped; ssGBLUP with one genotyped cohort: that cohort only).
    et = resource_estimate(t)
    @test keys(et.time_terms) == (:history, :setup, :meiosis, :summary, :evaluation, :mating, :reliability, :rest,
        :compile)
    @test sum(values(et.time_terms)) ≈ et.seconds atol = 0.1 * length(et.time_terms)
    x = BnGsim._time_features(t)
    @test x.dense == (2.0 * 40)^3 + (2.0 * 80)^3 && x.inverse == 40.0^3 + 80.0^3 && x.grm == (40.0^2 + 80.0^2) * 200
    @test x.equations == 2.0 * (40 + 80) && x.reml_equations == 0 && x.mating == 0
    @test BnGsim._time_features(parse_scenario(d)).inverse == 2 * 40.0^3
    @test x.meiosis == 2 * 2.0 * 40 * ScenarioSize(t).loci
    @test x.dense2 == (2.0 * 40)^2 + (2.0 * 80)^2 && x.equations2 == 4.0 * (40 + 80) && x.animals == 40 + 80
    @test x.inverse_ss == 0 && BnGsim._time_features(parse_scenario(d)).inverse_ss == 2 * 40.0^3
    # Dense memory (0.37): k² - 1 for k traits on top of one trait, k² for REML.
    @test BnGsim._dense_factor(t) == BnGsim._MEM.dense_gblup + 3 * BnGsim._MEM_TRAITS.gblup
    @test BnGsim._dense_factor(s) == 0
    @test BnGsim._time_features(s).dense == 0 && BnGsim._time_features(s).equations == 0   # minimal: phenotype
    # The phase timer of the calibration: every phase of the loop timed while
    # on, nothing collected while off.
    ph = BnGsim._phase_timing!()
    run_replicate(t, 1)
    BnGsim._phase_timing!(false)
    @test issubset((:setup, :evaluate, :kinv, :reliability, :select, :mating, :meiosis, :summary), keys(ph))
    @test all(>=(0), values(ph)) && BnGsim._PHASES[] === nothing

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

    # The time calibration (0.37) on its recorded timings (validation/
    # estimate/phases.tsv: replicate 2 of two processes per scenario):
    # every run over 1 s within 0.75 to 1.5 times the estimate of one warm
    # replicate with its history, except l_records (3 traits with repeated
    # records, 2.2 times: docs/src/estimate.md); and the peak of every
    # process, multi-trait and REML scenarios included, within the bound.
    dirs = [joinpath(cal, d) for d in ("scenarios", "later", "phases", "check")]
    push!(dirs, joinpath(@__DIR__, "..", "examples"))
    warm = Dict{String,Vector{Float64}}(); peak = Dict{String,Vector{Float64}}()
    for l in eachline(joinpath(cal, "phases.tsv"))
        r = split(l, '\t')
        r[3] == "2" || continue
        push!(get!(warm, r[2], Float64[]), parse(Float64, r[4]) + parse(Float64, r[5]) + parse(Float64, r[13]))
        push!(get!(peak, r[2], Float64[]), parse(Float64, r[end]) * 2^20)
    end
    @test length(warm) == 56
    for f in keys(warm)
        sc = read_scenario(only(filter(isfile, joinpath.(dirs, f))))
        tt = BnGsim._time_terms(sc)
        R = sc.replicates
        nh = sc.founders.shared_history ? min(sc.founders.histories, R) : R
        one = (sum(values(tt)) - tt.compile - tt.history) / R + tt.history / nh
        m = mean(warm[f])
        m > 1 && f != "l_records.json" && @test 0.75 <= m / one <= 1.5
        @test maximum(peak[f]) <= resource_estimate(sc).peak_rss_bound_bytes
    end

    # Command line
    code, v, _ = run_cli("estimate", MINIMAL)
    @test code == 0 && v["status"] == "estimated" && v["peak_rss_bytes"] == e.peak_rss_bytes
    @test v["peak_rss_bound_bytes"] == e.peak_rss_bound_bytes && v["size"]["animals"] == 200
    @test v["seconds"] == e.seconds && v["time_terms"]["compile"] == e.time_terms.compile
    @test run_cli("estimate")[1] == 2
    dir = mktempdir()
    code, v, _ = run_cli("run", MINIMAL, "--out", joinpath(dir, "a"), "--memory-limit", "100M")
    @test code == 5 && v["status"] == "rejected" && v["memory_limit_bytes"] == 100 * 2^20 && !ispath(joinpath(dir, "a"))
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "b"), "--memory-limit", "64G")[1] == 0
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "c"), "--memory-limit", "lots")[1] == 2
    @test BnGsim._parse_size("512M") == 512 * 2^20 && BnGsim._parse_size("1.5G") == round(Int, 1.5 * 2^30)
    rm(dir; recursive=true)
end
