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
| 5 | `BUDGET_EXCEEDED` | reserved: a resource budget was exceeded |
| 6 | `SIMULATION_FAILED` | the simulation stopped with an explicit error (e.g. too few loci pass a `maf` filter) |
| 7 | `CANCELLED` | reserved: cancelled at a checkpoint |
| 70 | `INTERNAL_ERROR` | an unexpected error (a bug) |

Codes 5 and 7 are reserved for `estimate`/budgets and `resume`, which the
CLI does not have yet.
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
       bngsim run SCENARIO.json --out DIR [--replicates A:B]
       bngsim version [--json]

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
    for m in (BnGStructs, FisherWright, RelationshipMatrices, Breeding, JSON, DataFrames, Distributions)
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
- `run SCENARIO.json --out DIR [--replicates A:B]`: run the scenario
  (replicates `A:B`, all by default) and write to `DIR`, which must not
  exist or be empty: `scenario.json` (the canonical scenario),
  `summary.csv` (the per-step summaries of every replicate) and
  `provenance.json`. Prints `{"status": "completed", "outputs": ...}`.
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
    scen, out, reps = nothing, nothing, nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out" || a == "--replicates"
            i < length(args) || return "$a needs a value"
            v = args[i+1]
            if a == "--out"
                out = v
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
    return (scen, out, reps)
end

function _cli_run(args, stdout, stderr)
    O = JSON.Object{String,Any}
    p = _parse_run_args(args)
    p isa String && return _usage(stderr, p)
    path, out, reps = p
    s = _load(path, stdout, stderr)
    s isa Integer && return s
    reps = something(reps, 1:s.replicates)
    if isempty(reps) || first(reps) < 1 || last(reps) > s.replicates
        msg = "--replicates $(first(reps)):$(last(reps)) is outside 1:$(s.replicates)"
        println(stderr, "bngsim: ", msg)
        _emit(stdout, O("status" => "rejected", "message" => msg))
        return ExitCode.REJECTED
    end
    if ispath(out) && !(isdir(out) && isempty(readdir(out)))
        msg = "output directory $out exists and is not empty; results are never overwritten"
        println(stderr, "bngsim: ", msg)
        _emit(stdout, O("status" => "rejected", "message" => msg))
        return ExitCode.REJECTED
    end
    j = scenario_json(s)
    started = now(UTC)
    t = @elapsed df = try
        run_scenario(s; replicates=reps)
    catch e
        e isa ArgumentError || e isa DimensionMismatch || e isa ErrorException || rethrow()
        msg = sprint(showerror, e)
        println(stderr, "bngsim: the simulation stopped: ", msg)
        _emit(stdout, O("status" => "failed", "message" => msg))
        return ExitCode.SIMULATION_FAILED
    end
    finished = now(UTC)
    mkpath(out)
    files = ["scenario.json", "summary.csv"]
    write(joinpath(out, "scenario.json"), j)
    open(io -> _write_csv(io, df), joinpath(out, "summary.csv"), "w")
    hashes = O(f => _sha256_file(joinpath(out, f)) for f in files)
    prov = O("bngsim" => string(pkgversion(BnGsim)), "cli_interface" => CLI_INTERFACE_VERSION,
        "scenario_schema" => SCENARIO_SCHEMA_VERSION, "scenario_sha256" => _sha256(j),
        "scenario" => scenario_dict(s), "seed" => s.seed, "replicates" => collect(reps),
        "julia" => string(VERSION), "packages" => _package_versions(), "manifest_sha256" => _manifest_sha256(),
        "os" => string(Sys.KERNEL, " ", Sys.MACHINE), "cpu" => Sys.cpu_info()[1].model, "host" => gethostname(),
        "threads" => Threads.nthreads(), "started" => string(started, "Z"), "finished" => string(finished, "Z"),
        "seconds" => round(t; digits=3), "outputs" => hashes)
    open(io -> _emit(io, prov), joinpath(out, "provenance.json"), "w")
    _emit(stdout, O("status" => "completed", "scenario_sha256" => _sha256(j), "replicates" => collect(reps),
        "outputs" => [joinpath(out, f) for f in [files; "provenance.json"]], "seconds" => round(t; digits=3)))
    return ExitCode.OK
end
