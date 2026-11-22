# Stopping a run at a step boundary and resuming it: the outputs are those
# of an uninterrupted run, byte for byte.

@testset "Resume" begin
    three = tempname() * ".json"
    write(three, replace(read(MINIMAL, String), "\"steps\": 3," => "\"steps\": 4, \"replicates\": 3,"))
    s = read_scenario(three)
    dir = mktempdir()
    ref = joinpath(dir, "ref")
    @test run_cli("run", three, "--out", ref)[1] == 0
    function start(out, every)
        mkpath(joinpath(out, ".bngsim-run"))
        run = (scenario=scenario_json(s), replicates=[1, 2, 3], checkpoint_every=every, sessions=Any[])
        BnGsim._save_state(joinpath(out, ".bngsim-run", "run.jls"), run)
        return run
    end

    # Stopped after step 3 of replicate 2 (as a file DIR/CANCEL would): exit
    # code 7 and the state; `run` refuses the directory; `resume` finishes it.
    a = joinpath(dir, "a")
    r = start(a, 0)
    out, err = IOBuffer(), IOBuffer()
    @test BnGsim._drive(a, s, r, out, err; created=false, should_stop=(r, g) -> r == 2 && g == 3) == 7
    v = JSON.parse(String(take!(out)))
    @test v["status"] == "cancelled" && v["replicate"] == 2 && v["step"] == 3
    @test sort(readdir(joinpath(a, ".bngsim-run"))) == ["current.jls", "done.jls", "history.jls", "run.jls"]
    code, v, _ = run_cli("run", three, "--out", a)
    @test code == 4 && occursin("bngsim resume", v["message"])
    touch(joinpath(a, "CANCEL"))                      # served by the stop; resuming withdraws it
    code, v, _ = run_cli("resume", a)
    @test code == 0 && v["sessions"] == 2
    @test sort(readdir(a)) == ["aggregate.csv", "provenance.json", "scenario.json", "summary.csv"]
    @test read(joinpath(a, "summary.csv")) == read(joinpath(ref, "summary.csv"))
    @test read(joinpath(a, "scenario.json")) == read(joinpath(ref, "scenario.json"))
    p = JSON.parse(read(joinpath(a, "provenance.json"), String))
    @test [x["status"] for x in p["sessions"]] == ["cancelled", "completed"] && p["replicates"] == [1, 2, 3]

    # A crash (the process dies) after the periodic checkpoint of step 2 in
    # replicate 3: resumed from there, the same bytes.
    c = joinpath(dir, "c")
    start(c, 2)
    @test_throws ErrorException BnGsim._execute(c, s, [1, 2, 3]; checkpoint_every=2,
        should_stop=(r, g) -> (r == 3 && g == 3 && error("killed"); false))
    @test isfile(joinpath(c, ".bngsim-run", "current.jls"))
    @test run_cli("resume", c)[1] == 0
    @test read(joinpath(c, "summary.csv")) == read(joinpath(ref, "summary.csv"))

    # The cancel file is what `run` checks.
    @test !BnGsim._cancel_requested(ref)
    touch(joinpath(ref, "CANCEL"))
    @test BnGsim._cancel_requested(ref)

    # Nothing to resume; a state from another version is refused.
    @test run_cli("resume", ref)[1] == 4 && run_cli("resume")[1] == 2
    d = joinpath(dir, "d")
    mkpath(joinpath(d, ".bngsim-run"))
    open(joinpath(d, ".bngsim-run", "run.jls"), "w") do io
        write(io, BnGsim._STATE_MAGIC)
        BnGsim.serialize(io, (bngsim="0.0.1", julia=string(VERSION)))
        BnGsim.serialize(io, nothing)
    end
    code, v, _ = run_cli("resume", d)
    @test code == 4 && occursin("does not carry over", v["message"])

    # --checkpoint-every: a positive number; a whole run with checkpoints
    # gives the same output and removes its state.
    @test run_cli("run", three, "--out", joinpath(dir, "e"), "--checkpoint-every", "0")[1] == 2
    @test run_cli("run", three, "--out", joinpath(dir, "f"), "--checkpoint-every", "1")[1] == 0
    @test read(joinpath(dir, "f", "summary.csv")) == read(joinpath(ref, "summary.csv"))
    @test !isdir(joinpath(dir, "f", ".bngsim-run"))
    rm(dir; recursive=true)
end
