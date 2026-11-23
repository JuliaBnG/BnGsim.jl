# Historical populations stored once and reused by several runs
# (docs/src/roadmap.md; decision of 2026-11-23): `write_histories` writes
# the K histories of a scenario as FisherWright result files with a
# manifest; a scenario with `founders.history_dir` reads them instead of
# simulating. A run refuses a manifest whose historical parameters differ
# from its own, so results never mix histories silently.

const HISTORY_SCHEMA = "bngsim-history"
const HISTORY_SCHEMA_VERSION = 1

# What determines a history: the species and the historical part of
# `founders` (not the founders sampled, nor the seed, nor the directory).
function _history_params(s::Scenario)
    d = scenario_dict(s)
    f = d["founders"]
    O = JSON.Object{String,Any}
    return O("species" => d["species"], "founders" => O(k => f[k] for k in
        ("ne", "nt", "mr", "decline", "popsize", "sexratio", "selection")))
end

"""
    write_histories(dir, s::Scenario) -> String

Simulate the `founders.histories` historical populations of `s` and write
them to `dir` (created; it must not exist or be empty): `history-<k>.fwr`
in FisherWright's result format and `manifest.json` with the historical
parameters, the SHA-256 of every file and the loci of the trait under
selection (if any). Returns the path of the manifest. A scenario with
`founders.history_dir = dir` and the same historical parameters (and
`histories`) then uses them, whatever its seed. Needs FisherWright 0.10.
"""
function write_histories(dir::AbstractString, s::Scenario)
    s.founders.shared_history || throw(ArgumentError("stored histories are shared: set shared_history = true"))
    isempty(s.founders.history_dir) || throw(ArgumentError("the scenario reads its histories from history_dir"))
    ispath(dir) && !(isdir(dir) && isempty(readdir(dir))) &&
        throw(ArgumentError("$dir exists and is not empty"))
    mkpath(dir)
    O = JSON.Object{String,Any}
    entries = Any[]
    for k = 1:s.founders.histories
        h = _history_of(s, k)
        f = "history-$k.fwr"
        write_history_file(joinpath(dir, f), h.res)
        push!(entries, O("file" => f, "sha256" => _sha256_file(joinpath(dir, f)),
            "trait_loci" => h.loci === nothing ? nothing : O("pos" => Int.(h.loci.pos), "a" => h.loci.a[:, 1],
                "mu" => h.loci.μ)))
    end
    m = O("schema" => HISTORY_SCHEMA, "schema_version" => HISTORY_SCHEMA_VERSION,
        "bngsim" => string(pkgversion(@__MODULE__)), "seed" => s.seed, "parameters" => _history_params(s),
        "histories" => entries)
    p = joinpath(dir, "manifest.json")
    open(io -> (JSON.json(io, m; pretty=true); println(io)), p, "w")
    return p
end

_manifest(s::Scenario) = JSON.parse(read(joinpath(s.founders.history_dir, "manifest.json"), String))

"""
    check_histories(s::Scenario) -> Union{Nothing,String}

For a scenario with `founders.history_dir`: `nothing` if the manifest
there is one of BnGsim, its historical parameters equal those of `s`, it
holds `founders.histories` histories and every file has the SHA-256 the
manifest records; a message saying what does not fit otherwise.
"""
function check_histories(s::Scenario)
    dir = s.founders.history_dir
    isempty(dir) && return nothing
    p = joinpath(dir, "manifest.json")
    isfile(p) || return "no manifest.json in history_dir $dir"
    m = try
        JSON.parse(read(p, String))
    catch
        return "$p is not JSON"
    end
    get(m, "schema", nothing) == HISTORY_SCHEMA && get(m, "schema_version", nothing) == HISTORY_SCHEMA_VERSION ||
        return "$p is not a BnGsim history manifest (schema $HISTORY_SCHEMA, version $HISTORY_SCHEMA_VERSION)"
    JSON.json(m["parameters"]) == JSON.json(JSON.parse(JSON.json(_history_params(s)))) ||
        return "the histories in $dir were made with other historical parameters (species or founders) than this scenario"
    n = length(m["histories"])
    n == s.founders.histories || return "$dir holds $n histories; the scenario asks for founders.histories = $(s.founders.histories)"
    for e in m["histories"]
        f = joinpath(dir, e["file"])
        isfile(f) || return "missing history file $f"
        _sha256_file(f) == e["sha256"] || return "$f does not match the SHA-256 in its manifest"
    end
    return nothing
end

# History `k` read from `founders.history_dir` (checked by `check_histories`).
function _load_history(s::Scenario, k::Integer)
    e = _manifest(s)["histories"][k]
    res = read_history_file(joinpath(s.founders.history_dir, e["file"]))
    tl = e["trait_loci"]
    loci = tl === nothing ? nothing : TraitLoci(UInt32.(tl["pos"]), Float64.(tl["a"]); h=0.5, μ=Float64(tl["mu"]))
    return (res=res, loci=loci)
end
