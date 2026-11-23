# The example scenarios of examples/: each reads, round-trips, validates
# and has an estimate, and the README lists it.

@testset "Example scenarios" begin
    dir = joinpath(@__DIR__, "..", "examples")
    files = filter(endswith(".json"), readdir(dir))
    @test sort(files) == ["dairy.json", "history-qtl.json", "minimal.json", "sheep.json"]
    for f in files
        s = read_scenario(joinpath(dir, f))
        @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
        @test run_cli("validate", joinpath(dir, f))[1] == 0
        @test resource_estimate(s).peak_rss_bytes > 0
        @test occursin("`$f`", read(joinpath(dir, "README.md"), String))
    end
end
