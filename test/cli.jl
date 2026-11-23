# The command line: arguments, exit codes and output shapes (the CLI
# contract, cli_interface 1).
import JSON

# Run `cli(args)`, returning the exit code, the parsed stdout (one JSON
# object, or the raw text) and stderr.
function run_cli(args...)
    out, err = IOBuffer(), IOBuffer()
    code = cli(collect(String, args); stdout=out, stderr=err)
    o = String(take!(out))
    parsed = try
        JSON.parse(o)
    catch
        o
    end
    return code, parsed, String(take!(err))
end

# An IO whose reader has gone away.
struct ClosedPipe <: IO end
Base.unsafe_write(::ClosedPipe, ::Ptr{UInt8}, ::UInt) = throw(Base.IOError("write: broken pipe", Base.UV_EPIPE))
Base.write(::ClosedPipe, ::UInt8) = throw(Base.IOError("write: broken pipe", Base.UV_EPIPE))

@testset "CLI: exit codes and output shapes" begin
    @test (ExitCode.OK, ExitCode.USAGE, ExitCode.INVALID_SCENARIO, ExitCode.REJECTED, ExitCode.BUDGET_EXCEEDED,
        ExitCode.SIMULATION_FAILED, ExitCode.CANCELLED, ExitCode.INTERNAL_ERROR) == (0, 2, 3, 4, 5, 6, 7, 70)

    code, v, _ = run_cli("version", "--json")
    @test code == 0 && v["cli_interface"] == 1 && v["scenario_schema"] == 1 && v["bngsim"] == string(pkgversion(BnGsim))
    @test all(haskey(v["packages"], p) for p in ("BnGsim", "FisherWright", "BnGStructs", "Breeding"))
    @test haskey(v, "julia") && haskey(v, "manifest_sha256")
    code, text, _ = run_cli("version")
    @test code == 0 && occursin("CLI interface 1", text)
    @test run_cli()[1] == 2 && run_cli("frobnicate")[1] == 2 && run_cli("version", "--yaml")[1] == 2
    @test run_cli("--help")[1] == 0

    # validate
    code, v, _ = run_cli("validate", MINIMAL)
    @test code == 0 && v["status"] == "valid" && v["scenario"] == "minimal"
    @test v["scenario_sha256"] == bytes2hex(BnGsim.sha256(scenario_json(read_scenario(MINIMAL))))
    bad = tempname() * ".json"
    write(bad, replace(read(MINIMAL, String), "\"nsire\"" => "\"nsires\"", "\"seed\": 2026" => "\"seed\": -1"))
    code, v, err = run_cli("validate", bad)
    # an unknown key, the required key it misspells, and the seed
    @test code == 3 && v["status"] == "invalid" && length(v["errors"]) == 3 && occursin("scheme.nsires", err)
    @test run_cli("validate", tempname())[1] == 3
    @test run_cli("validate")[1] == 2 && run_cli("validate", MINIMAL, MINIMAL)[1] == 2

    # run: outputs, provenance, refusal to overwrite
    dir = mktempdir()
    out = joinpath(dir, "out")
    code, v, _ = run_cli("run", MINIMAL, "--out", out)
    @test code == 0 && v["status"] == "completed" && v["replicates"] == [1]
    @test sort(readdir(out)) == ["aggregate.csv", "provenance.json", "scenario.json", "summary.csv"]
    @test read_scenario(joinpath(out, "scenario.json")) == read_scenario(MINIMAL)
    lines = readlines(joinpath(out, "summary.csv"))
    df = run_scenario(read_scenario(MINIMAL))
    @test split(lines[1], ',') == names(df) && length(lines) == nrow(df) + 1
    @test parse.(Float64, split(lines[3], ',')[5:6]) == [df.mean[2], df.var[2]]   # values round-trip
    p = JSON.parse(read(joinpath(out, "provenance.json"), String))
    @test p["scenario_sha256"] == v["scenario_sha256"] && p["seed"] == 2026 && p["threads"] == Threads.nthreads()
    @test p["scenario"]["name"] == "minimal" && haskey(p, "started") && haskey(p, "host")
    for f in ("scenario.json", "summary.csv")
        @test p["outputs"][f] == bytes2hex(open(BnGsim.sha256, joinpath(out, f)))
    end
    code, v, _ = run_cli("run", MINIMAL, "--out", out)          # not empty: never overwritten
    @test code == 4 && v["status"] == "rejected"
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "r"), "--replicates", "2:3")[1] == 4
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "s"), "--replicates", "x")[1] == 2
    @test run_cli("run", MINIMAL)[1] == 2 && run_cli("run", "--out", joinpath(dir, "t"))[1] == 2
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "u"), "--fast")[1] == 2
    @test run_cli("run", bad, "--out", joinpath(dir, "v"))[1] == 3 && !ispath(joinpath(dir, "v"))
    # A scenario that validates but cannot run: more QTL than segregating loci.
    many = tempname() * ".json"
    write(many, replace(read(MINIMAL, String), "\"nlc\": 50" => "\"nlc\": 1000000"))
    code, v, _ = run_cli("run", many, "--out", joinpath(dir, "w"))
    @test code == 6 && v["status"] == "failed" && !ispath(joinpath(dir, "w"))

    # A closed stdout ends the command quietly.
    @test cli(["version", "--json"]; stdout=ClosedPipe(), stderr=devnull) == 0

    # Through a separate process, as bin/bngsim runs it: the exit code
    # reaches the shell.
    jl = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e 'using BnGsim; exit(BnGsim.cli(ARGS))' --`
    p1 = run(pipeline(`$jl version --json`; stdout=devnull, stderr=devnull); wait=true)
    @test p1.exitcode == 0
    p2 = run(pipeline(ignorestatus(`$jl validate $bad`); stdout=devnull, stderr=devnull); wait=true)
    @test p2.exitcode == 3
    rm(dir; recursive=true)
end
