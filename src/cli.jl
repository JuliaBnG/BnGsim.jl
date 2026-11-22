# The command line (docs/src/roadmap.md §4.8): `bngsim validate | run |
# version`. The CLI interface is versioned separately from the package;
# exit codes are enumerated and stable; structured output (one JSON object)
# goes to stdout and diagnostics to stderr.

const CLI_INTERFACE_VERSION = 1

"""
    ExitCode

The exit codes of the `bngsim` command line. They are part of its interface
(`cli_interface` 1) and do not change meaning.

| Code | Name | Meaning |
| :--- | :--- | :--- |
| 0 | `OK` | success |
| 2 | `USAGE` | unknown command, missing or invalid arguments |
| 3 | `INVALID_SCENARIO` | the scenario file is missing, unreadable or does not validate |
| 4 | `REJECTED` | refused before running: the output directory exists and is not empty |
| 5 | `BUDGET_EXCEEDED` | refused before running: the memory bound of `resource_estimate` exceeds `--memory-limit` |
| 6 | `SIMULATION_FAILED` | the simulation stopped with an explicit error (e.g. too few loci pass a `maf` filter) |
| 7 | `CANCELLED` | stopped at a step boundary by a file `DIR/CANCEL` (which `bin/bngsim` creates on SIGINT or SIGTERM); `bngsim resume DIR` continues |
| 70 | `INTERNAL_ERROR` | an unexpected error (a bug) |

"""
module ExitCode
const OK = 0
const USAGE = 2
const INVALID_SCENARIO = 3
const REJECTED = 4
const BUDGET_EXCEEDED = 5
const SIMULATION_FAILED = 6
const CANCELLED = 7
const INTERNAL_ERROR = 70
end

const _CLI_USAGE = """
usage: bngsim validate SCENARIO.json
       bngsim estimate SCENARIO.json
       bngsim run SCENARIO.json --out DIR [--replicates A:B] [--memory-limit SIZE] [--checkpoint-every K]
                  [--archive none|chip|qtl|all]
       bngsim resume DIR
       bngsim merge DIR... --out DIR
       bngsim version [--json]

SIZE is bytes, or a number with K, M, G or T (powers of 1024), e.g. 16G.
A file DIR/CANCEL stops a run at the next step boundary (exit code 7);
`bngsim resume DIR` continues it.

Writes one JSON object to stdout and diagnostics to stderr; see the exit
codes in the BnGsim documentation (ExitCode)."""

_sha256(bytes) = bytes2hex(sha256(bytes))
_sha256_file(path) = open(io -> bytes2hex(sha256(io)), path)

function _emit(io::IO, x)
    JSON.json(io, x; pretty=2)
    println(io)
    return nothing
end

# Versions of BnGsim and the packages whose code a run executes.
function _package_versions()
    O = JSON.Object{String,Any}
    out = O("BnGsim" => string(pkgversion(BnGsim)))
    for m in (BnGStructs, FisherWrightAdapter.PACKAGE, RelationshipsAdapter.PACKAGE, BreedingAdapter.PACKAGE, JSON,
        DataFrames, Distributions)
        out[string(nameof(m))] = string(pkgversion(m))
    end
    return out
end

function _manifest_sha256()
    p = Base.active_project()
    p === nothing && return nothing
    for name in ("Manifest-v$(VERSION.major).$(VERSION.minor).toml", "Manifest.toml",
        "JuliaManifest.toml")
        m = joinpath(dirname(p), name)
        isfile(m) && return _sha256_file(m)
    end
    return nothing
end

function _version_info()
    O = JSON.Object{String,Any}
    return O("bngsim" => string(pkgversion(BnGsim)), "cli_interface" => CLI_INTERFACE_VERSION,
        "scenario_schema" => SCENARIO_SCHEMA_VERSION, "julia" => string(VERSION),
        "packages" => _package_versions(), "manifest_sha256" => _manifest_sha256())
end

# A data frame as CSV: a header, then one line per row; numbers in Julia's
# shortest round-trip form, missing values empty. The same values give the
# same bytes.
function _write_csv(io::IO, df::DataFrame)
    println(io, join(names(df), ','))
    for r in eachrow(df)
        println(io, join((_csv(x) for x in r), ','))
    end
end
_csv(::Missing) = ""
_csv(x::AbstractString) = occursin(r"[,\"\n]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x
_csv(x) = string(x)

"""
    cli(args::Vector{String}; stdout = stdout, stderr = stderr) -> Int

The `bngsim` command line: run the command in `args` and return its exit
code ([`ExitCode`](@ref)). `bin/bngsim` calls it as
`exit(BnGsim.cli(ARGS))`.

- `validate SCENARIO.json`: read and check the scenario. Prints
  `{"status": "valid", "scenario_sha256", ...}` or
  `{"status": "invalid", "errors": [...]}` (exit code 3).
- `estimate SCENARIO.json`: the expected peak memory, a conservative
  bound, the disk written and a rough run time ([`resource_estimate`](@ref)).
- `run SCENARIO.json --out DIR [--replicates A:B] [--memory-limit SIZE]`:
  run the scenario (replicates `A:B`, all by default) and write to `DIR`,
  which must not exist or be empty: `scenario.json` (the canonical
  scenario), `summary.csv` (the per-step summaries of every replicate) and
  `provenance.json`. Prints `{"status": "completed", "outputs": ...}`.
  With `--memory-limit`, a scenario whose memory bound exceeds `SIZE` is
  refused before it runs (exit code 5). While it runs, the state lives in
  `DIR/.bngsim-run/`; a file `DIR/CANCEL` stops it after the current
  time step, with the state saved (exit code 7), and with
  `--checkpoint-every K` the state is also saved every `K` steps, so that
  a killed process can be resumed.
  With `--archive none|chip|qtl|all`, an archive of each replicate
  ([`write_archive`](@ref)) is written to `DIR/archive/r<N>/` when it
  finishes, and listed in `provenance.json`.
- `merge DIR... --out OUT`: combine completed runs of the same scenario with
  disjoint replicates (e.g. split over hosts with `--replicates`) into
  `OUT`: `scenario.json`, `summary.csv` (rows ordered by replicate, the
  bytes of a single run of all of them), `aggregate.csv` and
  `provenance.json` listing the runs merged.
- `resume DIR`: continue a stopped or killed run from its last saved
  state (removing `DIR/CANCEL`). The outputs are those of an
  uninterrupted run, byte for byte; `provenance.json` lists every session.
  The state must come from the same BnGsim and Julia versions (exit code
  4 otherwise).
- `version [--json]`: the BnGsim version, the CLI interface and scenario
  schema versions, the Julia and package versions and the hash of the
  active manifest.

`summary.csv` depends only on the scenario, the replicates and the
package versions: not on the number of threads or the host.
`provenance.json` records the rest of the run (roadmap invariant I2): the
canonical scenario and its SHA-256, the versions and manifest hash, the
host, CPU, threads, start and end times, and the SHA-256 of each output.
"""
function cli(args::AbstractVector{<:AbstractString}; stdout::IO=Base.stdout, stderr::IO=Base.stderr)
    try
        isempty(args) && return _usage(stderr, "no command")
        cmd, rest = args[1], args[2:end]
        cmd == "validate" && return _cli_validate(rest, stdout, stderr)
        cmd == "run" && return _cli_run(rest, stdout, stderr)
        cmd == "estimate" && return _cli_estimate(rest, stdout, stderr)
        cmd == "resume" && return _cli_resume(rest, stdout, stderr)
        cmd == "merge" && return _cli_merge(rest, stdout, stderr)
        cmd == "version" && return _cli_version(rest, stdout, stderr)
        cmd in ("-h", "--help", "help") && (println(stdout, _CLI_USAGE); return ExitCode.OK)
        return _usage(stderr, "unknown command \"$cmd\"")
    catch e
        e isa InterruptException && rethrow()
        # The reader of stdout went away (e.g. `bngsim version --json | head`):
        # stop quietly, as other command-line tools do.
        e isa Base.IOError && e.code == Base.UV_EPIPE && return ExitCode.OK
        try
            println(stderr, "bngsim: internal error: ", sprint(showerror, e, catch_backtrace()))
            _emit(stdout, JSON.Object{String,Any}("status" => "error", "message" => sprint(showerror, e)))
        catch
        end
        return ExitCode.INTERNAL_ERROR
    end
end

function _usage(stderr, msg)
    println(stderr, "bngsim: ", msg)
    println(stderr, _CLI_USAGE)
    return ExitCode.USAGE
end

# The scenario at `path`, or an exit code after reporting why not.
function _load(path, stdout, stderr)
    O = JSON.Object{String,Any}
    if !isfile(path)
        println(stderr, "bngsim: no scenario file ", path)
        _emit(stdout, O("status" => "invalid", "errors" => ["no such file: $path"]))
        return ExitCode.INVALID_SCENARIO
    end
    try
        return parse_scenario(read(path, String))
    catch e
        e isa ScenarioError || rethrow()
        println(stderr, sprint(showerror, e))
        _emit(stdout, O("status" => "invalid", "errors" => e.errors))
        return ExitCode.INVALID_SCENARIO
    end
end

function _cli_validate(args, stdout, stderr)
    length(args) == 1 || return _usage(stderr, "validate takes one scenario file")
    s = _load(args[1], stdout, stderr)
    s isa Integer && return s
    j = scenario_json(s)
    _emit(stdout, JSON.Object{String,Any}("status" => "valid", "scenario" => s.name,
        "schema_version" => SCENARIO_SCHEMA_VERSION, "scenario_sha256" => _sha256(j),
        "replicates" => s.replicates, "steps" => s.steps))
    return ExitCode.OK
end

function _cli_version(args, stdout, stderr)
    isempty(args) || args == ["--json"] || return _usage(stderr, "version takes only --json")
    if args == ["--json"]
        _emit(stdout, _version_info())
    else
        println(stdout, "bngsim ", pkgversion(BnGsim), " (CLI interface ", CLI_INTERFACE_VERSION,
            ", scenario schema ", SCENARIO_SCHEMA_VERSION, ", Julia ", VERSION, ")")
    end
    return ExitCode.OK
end

function _parse_run_args(args)
    scen, out, reps, limit, every, archive = nothing, nothing, nothing, nothing, 0, nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--out", "--replicates", "--memory-limit", "--checkpoint-every", "--archive")
            i < length(args) || return "$a needs a value"
            v = args[i+1]
            if a == "--out"
                out = v
            elseif a == "--archive"
                v in ("none", "chip", "qtl", "all") || return "--archive takes none, chip, qtl or all"
                archive = Symbol(v)
            elseif a == "--checkpoint-every"
                k = tryparse(Int, v)
                (k === nothing || k < 1) && return "--checkpoint-every takes a positive number of steps"
                every = k
            elseif a == "--memory-limit"
                limit = _parse_size(v)
                limit === nothing && return "--memory-limit takes bytes or a number with K, M, G or T, e.g. 16G"
            else
                m = match(r"^(\d+)(?::(\d+))?$", v)
                m === nothing && return "--replicates takes A or A:B, e.g. 1:10"
                lo = parse(Int, m[1])
                hi = m[2] === nothing ? lo : parse(Int, m[2])
                reps = lo:hi
            end
            i += 2
        elseif startswith(a, "-")
            return "unknown option $a"
        else
            scen === nothing || return "run takes one scenario file"
            scen = a
            i += 1
        end
    end
    scen === nothing && return "run needs a scenario file"
    out === nothing && return "run needs --out DIR"
    return (scen, out, reps, limit, every, archive)
end

function _cli_run(args, stdout, stderr)
    O = JSON.Object{String,Any}
    p = _parse_run_args(args)
    p isa String && return _usage(stderr, p)
    path, out, reps, limit, every, archive = p
    s = _load(path, stdout, stderr)
    s isa Integer && return s
    reps = something(reps, 1:s.replicates)
    if isempty(reps) || first(reps) < 1 || last(reps) > s.replicates
        return _reject(stdout, stderr, "--replicates $(first(reps)):$(last(reps)) is outside 1:$(s.replicates)")
    end
    archive === :chip && s.chip === nothing &&
        return _reject(stdout, stderr, "--archive chip needs a chip (loci.chip) in the scenario")
    if isdir(_state_dir(out))
        return _reject(stdout, stderr, "$out holds an unfinished run; continue it with `bngsim resume $out`")
    elseif ispath(out) && !(isdir(out) && isempty(readdir(out)))
        return _reject(stdout, stderr, "output directory $out exists and is not empty; results are never overwritten")
    end
    if limit !== nothing
        e = resource_estimate(s)
        if e.peak_rss_bound_bytes > limit
            msg = "the memory bound of the scenario, $(_human(e.peak_rss_bound_bytes)), exceeds --memory-limit " *
                  "$(_human(limit)) (expected peak $(_human(e.peak_rss_bytes)))"
            println(stderr, "bngsim: ", msg)
            _emit(stdout, O("status" => "rejected", "message" => msg, "peak_rss_bytes" => e.peak_rss_bytes,
                "peak_rss_bound_bytes" => e.peak_rss_bound_bytes, "memory_limit_bytes" => limit))
            return ExitCode.BUDGET_EXCEEDED
        end
    end
    created = !ispath(out)
    mkpath(_state_dir(out))
    run = (scenario=scenario_json(s), replicates=collect(reps), checkpoint_every=every, archive=archive,
        sessions=Any[])
    _save_state(joinpath(_state_dir(out), "run.jls"), run)
    return _drive(out, s, run, stdout, stderr; created)
end

_reject(stdout, stderr, msg) = (println(stderr, "bngsim: ", msg);
    _emit(stdout, JSON.Object{String,Any}("status" => "rejected", "message" => msg)); ExitCode.REJECTED)

function _cli_resume(args, stdout, stderr)
    O = JSON.Object{String,Any}
    length(args) == 1 || return _usage(stderr, "resume takes one run directory")
    out = args[1]
    p = joinpath(_state_dir(out), "run.jls")
    isfile(p) || return _reject(stdout, stderr, "$out holds no unfinished run")
    run = try
        _load_state(p)
    catch e
        e isa ArgumentError || rethrow()
        return _reject(stdout, stderr, e.msg)
    end
    s = parse_scenario(run.scenario)
    # The request to stop has been served; resuming withdraws it.
    rm(joinpath(out, "CANCEL"); force=true)
    return _drive(out, s, run, stdout, stderr; created=false)
end

# A request to stop the run in `out` at the next step boundary.
_cancel_requested(out) = isfile(joinpath(out, "CANCEL"))

# Run (or continue) the run in `out` and finish it: write the outputs, or
# the state and exit code 7 when stopped.
function _drive(out, s, run, stdout, stderr; created, should_stop=(r, g) -> _cancel_requested(out))
    O = JSON.Object{String,Any}
    st = _state_dir(out)
    reps = run.replicates
    started = now(UTC)
    t0 = time()
    result = try
        _execute(out, s, reps; checkpoint_every=run.checkpoint_every, should_stop=should_stop,
            on_replicate=(r, rep) -> _archive(out, _run_archive(run), r, rep))
    catch e
        e isa ArgumentError || e isa DimensionMismatch || e isa ErrorException || rethrow()
        msg = sprint(showerror, e)
        println(stderr, "bngsim: the simulation stopped: ", msg)
        _emit(stdout, O("status" => "failed", "message" => msg))
        # A run that never saved state leaves nothing; one with a checkpoint
        # keeps it, for inspection or another `resume`.
        if isempty(run.sessions) && !isfile(joinpath(st, "current.jls")) && !isfile(joinpath(st, "done.jls"))
            created ? rm(out; recursive=true) : rm(st; recursive=true)
        end
        return ExitCode.SIMULATION_FAILED
    end
    session = O("started" => string(started, "Z"), "finished" => string(now(UTC), "Z"),
        "seconds" => round(time() - t0; digits=3), "host" => gethostname(), "cpu" => Sys.cpu_info()[1].model,
        "threads" => Threads.nthreads(), "status" => result[1] === :cancelled ? "cancelled" : "completed")
    sessions = [run.sessions; session]
    j = run.scenario
    if result[1] === :cancelled
        _, r, g = result
        _save_state(joinpath(st, "run.jls"), merge(run, (sessions=sessions,)))
        println(stderr, "bngsim: stopped after step $g of replicate $r; continue with `bngsim resume $out`")
        _emit(stdout, O("status" => "cancelled", "replicate" => r, "step" => g, "directory" => out))
        return ExitCode.CANCELLED
    end
    df = result[2]
    files = ["scenario.json", "summary.csv", "aggregate.csv"]
    write(joinpath(out, "scenario.json"), j)
    open(io -> _write_csv(io, df), joinpath(out, "summary.csv"), "w")
    open(io -> _write_csv(io, aggregate(df)), joinpath(out, "aggregate.csv"), "w")
    hashes = O(f => _sha256_file(joinpath(out, f)) for f in files)
    seconds = round(sum(x -> x["seconds"], sessions); digits=3)
    prov = O("bngsim" => string(pkgversion(BnGsim)), "cli_interface" => CLI_INTERFACE_VERSION,
        "scenario_schema" => SCENARIO_SCHEMA_VERSION, "scenario_sha256" => _sha256(j),
        "scenario" => scenario_dict(s), "seed" => s.seed, "replicates" => reps,
        "julia" => string(VERSION), "packages" => _package_versions(), "manifest_sha256" => _manifest_sha256(),
        "os" => string(Sys.KERNEL, " ", Sys.MACHINE), "cpu" => Sys.cpu_info()[1].model, "host" => gethostname(),
        "threads" => Threads.nthreads(), "started" => sessions[1]["started"], "finished" => session["finished"],
        "seconds" => seconds, "checkpoint_every" => run.checkpoint_every, "sessions" => sessions,
        "outputs" => hashes)
    if _run_archive(run) !== nothing
        dir = joinpath(out, "archive")
        prov["archive"] = O("panel" => String(_run_archive(run)), "files" => O(relpath(f, out) => _sha256_file(f)
            for f in sort(_files_under(dir))))
    end
    open(io -> _emit(io, prov), joinpath(out, "provenance.json"), "w")
    rm(st; recursive=true)
    rm(joinpath(out, "CANCEL"); force=true)
    _emit(stdout, O("status" => "completed", "scenario_sha256" => _sha256(j), "replicates" => reps,
        "outputs" => [joinpath(out, f) for f in [files; "provenance.json"]], "seconds" => seconds,
        "sessions" => length(sessions)))
    return ExitCode.OK
end

# Bytes from "123", "512M", "16G", "1.5T" (powers of 1024).
function _parse_size(v::AbstractString)
    m = match(r"^(\d+(?:\.\d+)?)([KMGT]?)i?B?$"i, strip(v))
    m === nothing && return nothing
    x = parse(Float64, m[1]) * 1024.0^(findfirst(==(uppercase(m[2])), ["", "K", "M", "G", "T"]) - 1)
    return x >= 1 ? round(Int, x) : nothing
end

_human(b) = b >= 2^30 ? string(round(b / 2^30; digits=2), " GiB") : string(round(b / 2^20; digits=1), " MiB")

function _cli_estimate(args, stdout, stderr)
    length(args) == 1 || return _usage(stderr, "estimate takes one scenario file")
    s = _load(args[1], stdout, stderr)
    s isa Integer && return s
    e = resource_estimate(s)
    z = e.size
    O = JSON.Object{String,Any}
    _emit(stdout, O("status" => "estimated", "scenario" => s.name, "scenario_sha256" => _sha256(scenario_json(s)),
        "peak_rss_bytes" => e.peak_rss_bytes, "peak_rss_bound_bytes" => e.peak_rss_bound_bytes,
        "memory_terms" => O(string(k) => v for (k, v) in pairs(e.memory_terms)),
        "disk_bytes" => e.disk_bytes, "seconds" => e.seconds,
        "size" => O(string(f) => getfield(z, f) for f in fieldnames(ScenarioSize)),
        "calibration" => "AMD Threadripper 3970X, Linux x86_64, Julia 1.13.4, 8 threads; memory within ±15 % " *
                         "(scenario means), time within a factor of about two"))
    return ExitCode.OK
end

# The archive of replicate `r` in `out/archive/r<r>/`.
function _archive(out, panel, r, rep)
    panel === nothing && return nothing
    loci = panel === :chip ? rep.chip : panel === :qtl ? rep.qtl : panel === :all ? (1:length(rep.map)) : Int32[]
    write_archive(joinpath(out, "archive", "r$r"), rep.population, rep.map, rep.species; panel=loci, qtl=rep.qtl)
    return nothing
end

_run_archive(run) = hasproperty(run, :archive) ? run.archive : nothing

_files_under(dir) = isdir(dir) ? [joinpath(root, f) for (root, _, fs) in walkdir(dir) for f in fs] : String[]

function _cli_merge(args, stdout, stderr)
    O = JSON.Object{String,Any}
    i = findfirst(==("--out"), args)
    (i === nothing || i == length(args)) && return _usage(stderr, "merge needs --out DIR")
    out = args[i+1]
    dirs = [args[1:i-1]; args[i+2:end]]
    length(dirs) >= 1 || return _usage(stderr, "merge needs the run directories to merge")
    any(startswith("-"), dirs) && return _usage(stderr, "merge takes run directories and --out DIR")
    ispath(out) && !(isdir(out) && isempty(readdir(out))) &&
        return _reject(stdout, stderr, "output directory $out exists and is not empty; results are never overwritten")
    provs, hdr = [], nothing
    rows = Dict{Int,Vector{String}}()
    for d in dirs
        p = joinpath(d, "provenance.json")
        (isfile(p) && isfile(joinpath(d, "summary.csv")) && !isdir(_state_dir(d))) ||
            return _reject(stdout, stderr, "$d is not a completed run (provenance.json and summary.csv, no unfinished state)")
        prov = JSON.parse(read(p, String))
        push!(provs, (dir=d, prov=prov, sha=_sha256_file(p)))
        lines = readlines(joinpath(d, "summary.csv"))
        hdr === nothing ? (hdr = lines[1]) : hdr == lines[1] ||
            return _reject(stdout, stderr, "$d has other summary columns than $(dirs[1])")
        for l in lines[2:end]
            r = parse(Int, split(l, ',')[1])
            push!(get!(rows, r, String[]), l)
        end
    end
    shas = unique(String(x.prov["scenario_sha256"]) for x in provs)
    length(shas) == 1 || return _reject(stdout, stderr, "the runs are of different scenarios (scenario_sha256 differ)")
    reps = [Int.(x.prov["replicates"]) for x in provs]
    all_reps = reduce(vcat, reps)
    allunique(all_reps) || return _reject(stdout, stderr, "the runs share replicates; merge needs disjoint ones")
    mkpath(out)
    cp(joinpath(dirs[1], "scenario.json"), joinpath(out, "scenario.json"))
    open(joinpath(out, "summary.csv"), "w") do io
        println(io, hdr)
        for r in sort(collect(keys(rows))), l in rows[r]
            println(io, l)
        end
    end
    df = _read_summary(joinpath(out, "summary.csv"))
    open(io -> _write_csv(io, aggregate(df)), joinpath(out, "aggregate.csv"), "w")
    files = ["scenario.json", "summary.csv", "aggregate.csv"]
    hashes = O(f => _sha256_file(joinpath(out, f)) for f in files)
    prov = O("bngsim" => string(pkgversion(BnGsim)), "cli_interface" => CLI_INTERFACE_VERSION,
        "scenario_sha256" => shas[1], "replicates" => sort(all_reps), "merged" => string(now(UTC), "Z"),
        "host" => gethostname(), "julia" => string(VERSION),
        "sources" => [O("directory" => abspath(x.dir), "provenance_sha256" => x.sha, "replicates" => r,
            "bngsim" => x.prov["bngsim"]) for (x, r) in zip(provs, reps)], "outputs" => hashes)
    open(io -> _emit(io, prov), joinpath(out, "provenance.json"), "w")
    _emit(stdout, O("status" => "merged", "scenario_sha256" => shas[1], "replicates" => sort(all_reps),
        "outputs" => [joinpath(out, f) for f in [files; "provenance.json"]]))
    return ExitCode.OK
end
