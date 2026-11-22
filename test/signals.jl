# SIGTERM and SIGINT to bin/bngsim stop a run at a step boundary with its
# state saved (exit code 7), through DIR/CANCEL; a resumed run is
# byte-identical.

@testset "Signals" begin
    dir = mktempdir()
    # A run in another process, signalled once its state exists, stops with
    # exit code 7 and resumes to the bytes of an uninterrupted run.
    if Sys.isunix()
        long = tempname() * ".json"
        write(long, replace(read(MINIMAL, String), "\"steps\": 3," => "\"steps\": 120,"))
        launcher = joinpath(@__DIR__, "..", "bin", "bngsim")
        env = copy(ENV)
        env["JULIA"] = joinpath(Sys.BINDIR, "julia")
        env["BNGSIM_PROJECT"] = dirname(Base.active_project())
        env["BNGSIM_THREADS"] = "1"
        # Every run of this test goes through the launcher: the test process
        # runs with --check-bounds=yes, which changes the last bits of
        # floating-point sums, and the launcher's Julia does not.
        cli_run(args...) = run(pipeline(ignorestatus(setenv(`$launcher $args`, env)); stdout=devnull, stderr=devnull)).exitcode
        lref = joinpath(dir, "lref")
        @test cli_run("run", long, "--out", lref) == 0
        for sig in (Base.SIGTERM, Base.SIGINT)
            b = joinpath(dir, "b$sig")
            p = run(pipeline(setenv(`$launcher run $long --out $b`, env); stdout=devnull, stderr=devnull); wait=false)
            t0 = time()
            while !isdir(joinpath(b, ".bngsim-run")) && time() - t0 < 300 && process_running(p)
                sleep(0.05)
            end
            kill(p, sig)
            wait(p)
            @test p.exitcode == 7
            @test isfile(joinpath(b, ".bngsim-run", "current.jls")) && isfile(joinpath(b, "CANCEL"))
            @test cli_run("resume", b) == 0 && !isfile(joinpath(b, "CANCEL"))
            @test read(joinpath(b, "summary.csv")) == read(joinpath(lref, "summary.csv"))
        end
        # Before the run has saved state the launcher stops it at once: no
        # output directory is left.
        c = joinpath(dir, "c")
        p = run(pipeline(setenv(`$launcher run $long --out $c`, env); stdout=devnull, stderr=devnull); wait=false)
        sleep(0.5)
        kill(p, Base.SIGTERM)
        wait(p)
        @test p.exitcode != 0 && !isdir(joinpath(c, ".bngsim-run"))
    end
    rm(dir; recursive=true)
end
