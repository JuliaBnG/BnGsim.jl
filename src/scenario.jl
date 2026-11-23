# Scenario schema (docs/src/roadmap.md, Phase 1): a versioned JSON
# description of a whole simulation (species, historical population,
# founders, loci, traits, breeding scheme, steps, seed and replicates), its
# validation, its canonical JSON form, and its execution through the API.

const SCENARIO_SCHEMA = "bngsim-scenario"
const SCENARIO_SCHEMA_VERSION = 1

"""
    ScenarioError(errors)

A scenario that cannot be read or does not validate. `errors` lists every
problem found, each prefixed by its path in the JSON document (for
example `scheme.nsire`).
"""
struct ScenarioError <: Exception
    errors::Vector{String}
end

function Base.showerror(io::IO, e::ScenarioError)
    n = length(e.errors)
    print(io, "ScenarioError: ", n, " problem", n == 1 ? "" : "s", " in the scenario")
    for m in e.errors
        print(io, "\n  ", m)
    end
end

"""
    Scenario

A validated simulation scenario, read from JSON with [`read_scenario`](@ref)
or [`parse_scenario`](@ref), written with [`write_scenario`](@ref), and run
with [`run_scenario`](@ref). Every optional value is filled with its
default, so the canonical JSON ([`scenario_json`](@ref)) states the whole
simulation. Two scenarios are equal when all their values are.

The JSON document (schema `"bngsim-scenario"`, version 1) has these keys;
those without a default are required.

**Top level**
- `schema`: `"bngsim-scenario"`; `schema_version`: `1`;
- `name`: a label;
- `seed`: the master seed, an integer in `0:2^64-1` (every random stream
  derives from it through [`SimKey`](@ref));
- `replicates` (`1`): replicates `1:replicates`;
- `steps`: time steps of the breeding scheme.

**`species`** either `{"preset": "Cattle"}` (one of `Cattle`, `Sheep`,
`Goat`, `Pig`, `Horse`, `Dog`, `Cat`, `Chicken`, `Rabbit`), or
`{"name", "chromosome"}` with autosome lengths in bp; and either `M`, bp per
Morgan (`1e8`), or `morgan`, the genetic length of each chromosome. The
canonical form lists `name`, `chromosome` and `morgan`.

**`founders`**: the historical population, run with FisherWright from no
variation: `ne` diploids for `nt` generations at `mr` mutations per 1e8 bp
per meiosis (`1.0`); then `decline` (`[]`), the census sizes of further
generations (e.g. a decline to a breed); and `nid`, the founders sampled
from the last generation. With `shared_history` (`true`) one historical
population serves every replicate, each sampling its own founders;
otherwise each replicate runs its own. `herds` (`1`): the number of herds
the founders are spread over at random; offspring are born in their
dam's herd (contemporary groups, see `c2` of the traits).

**`loci`**: `qtl`, `{"nlc", "maf"}` (`maf` `0.0`: loci with minor allele
frequency above `maf` among the founders); `chip` (`null`), the same for
the genotyping panel, needed by `:gblup` and `:ssgblup`; `overlap`
(`"disjoint"`): whether chip loci may also be QTL (`"allowed"`).

**`traits`**: a non-empty array of `{"name", "h2" (0.25), "sigma_a" (1.0),
"mu" (0.0), "sex" (2), "age" (0), "da"}`, as `BnGStructs.Trait`: `sex` `0`
females, `1` males, `2` both; `age` the whole time step of the record; `da`
the distribution of the QTL effects before scaling,
`{"family": "normal", "mean": 0.0, "sd": 1.0}` (default),
`{"family": "gamma", "shape", "scale": 1.0, "sign": "positive"}` (`"sign":
"random"` for effects of either sign, [`RandomSign`](@ref)),
`{"family": "laplace", "scale": 1.0}` (location 0) or
`{"family": "spike_slab", "pi", "sd": 1.0, "spike_sd": 0.0}`
([`SpikeSlab`](@ref): a fraction `pi` of the effects from `Normal(0, sd)`,
the others `0`, or from `Normal(0, spike_sd)`). Effects are scaled to
`sigma_a` afterwards, so only the shape of the distribution matters; with
several traits, only the first keeps it (see [`eqtl`](@ref)). `c2` (`0.0`, in `[0, 1)`): the variance of the contemporary-group
effects (herd × time step of the record) as a fraction of the founders'
`sigma_a² + sigma_e²`; `h2` is then the heritability within groups. Any
`c2 > 0` makes evaluations fit the groups as fixed effects.

**`genetic_correlation`**, **`residual_correlation`** (identity): `k × k`
correlation matrices of `k` traits, as arrays of rows.

**`scheme`**: the keywords of [`Scheme`](@ref), in ASCII: `noff`, `nsire`,
`ndam`, `evaluation`, `selection`, `K`, `lambda`, `kappa`, `ub`, `blend`,
`max_sire_age`, `max_dam_age`, `reml`, `weights`, `survival`,
`cull_unselected`, `genotyped_cohorts`, `ocs`, `delta_F`, `record_meiosis`,
`window`, `ancestor_generations`, `retain`, `report_reliability`, with
the same defaults;
`null` for `K` and `delta_F` (unset) and for `ancestor_generations` (all).
The chip loci come from `loci.chip`.

Unknown keys, duplicate keys and values of the wrong type are errors.
"""
struct Scenario
    name::String
    seed::UInt64
    replicates::Int
    steps::Int
    species::@NamedTuple{name::String, chromosome::Vector{Int}, morgan::Vector{Float64}}
    founders::@NamedTuple{ne::Int, nt::Int, mr::Float64, decline::Vector{Int}, nid::Int, shared_history::Bool,
        herds::Int}
    qtl::@NamedTuple{nlc::Int, maf::Float64}
    chip::Union{Nothing,@NamedTuple{nlc::Int, maf::Float64}}
    overlap::Symbol
    traits::Vector{@NamedTuple{name::String, h2::Float64, sigma_a::Float64, mu::Float64, sex::Int, age::Int,
        da::@NamedTuple{family::Symbol, a::Float64, b::Float64, c::Float64, sign::Symbol}, c2::Float64}}
    genetic_correlation::Matrix{Float64}
    residual_correlation::Matrix{Float64}
    scheme::NamedTuple
end

Base.:(==)(a::Scenario, b::Scenario) = all(f -> isequal(getfield(a, f), getfield(b, f)), fieldnames(Scenario))
Base.hash(s::Scenario, h::UInt) = hash(map(f -> getfield(s, f), fieldnames(Scenario)), h)
ntraits(s::Scenario) = length(s.traits)

Base.show(io::IO, s::Scenario) = print(io, "Scenario(\"", s.name, "\", ", ntraits(s), " trait",
    ntraits(s) == 1 ? "" : "s", ", ", s.steps, " steps, ", s.replicates, " replicate",
    s.replicates == 1 ? "" : "s", ")")

# ---- reading ---------------------------------------------------------------

# A JSON object that refuses a repeated key (JSON.jl keeps the last one).
struct _StrictDict <: AbstractDict{String,Any}
    d::Dict{String,Any}
    order::Vector{String}
end
_StrictDict() = _StrictDict(Dict{String,Any}(), String[])
function Base.setindex!(s::_StrictDict, v, k::AbstractString)
    haskey(s.d, k) && throw(ArgumentError("duplicate key \"$k\""))
    s.d[k] = v
    push!(s.order, k)
    return s
end
Base.length(s::_StrictDict) = length(s.d)
Base.getindex(s::_StrictDict, k) = s.d[k]
Base.get(s::_StrictDict, k, d) = get(s.d, k, d)
Base.haskey(s::_StrictDict, k) = haskey(s.d, k)
Base.keys(s::_StrictDict) = s.order
Base.iterate(s::_StrictDict, i=1) = i > length(s.order) ? nothing : (s.order[i] => s.d[s.order[i]], i + 1)

const _PRESETS = Dict("Cattle" => Cattle, "Sheep" => Sheep, "Goat" => Goat, "Pig" => Pig, "Horse" => Horse,
    "Dog" => Dog, "Cat" => Cat, "Chicken" => Chicken, "Rabbit" => Rabbit)

const _NONE = Val(:required)    # marks a required key

# Collects the problems of one document.
struct _Ctx
    errors::Vector{String}
end
_err!(c::_Ctx, path, msg) = (push!(c.errors, isempty(path) ? msg : "$path: $msg"); nothing)
_join(path, key) = isempty(path) ? key : "$path.$key"

_isobject(x) = x isa AbstractDict
function _known!(c::_Ctx, d, path, keys)
    for k in Base.keys(d)
        k in keys || _err!(c, _join(path, k), "unknown key (expected one of: $(join(keys, ", ")))")
    end
end

# The value of `key` in object `d`, converted to `kind`; `nothing` (with an
# error recorded) if it is missing or invalid. `null` is accepted where
# `nullable`, as `null_as`.
function _get(c::_Ctx, d, path, key, kind; default=_NONE, nullable=false, null_as=nothing,
    choices=nothing, check=nothing)
    p = _join(path, key)
    if !haskey(d, key)
        default === _NONE && (_err!(c, p, "required"); return nothing)
        return default
    end
    v = d[key]
    if v === nothing
        nullable && return null_as
        _err!(c, p, "must not be null")
        return nothing
    end
    x = _convert(c, p, v, kind)
    x === nothing && return nothing
    if choices !== nothing && !(x in choices)
        _err!(c, p, "must be one of: $(join(repr.(choices), ", ")); got $(repr(x))")
        return nothing
    end
    if check !== nothing
        msg = check(x)
        msg === nothing || (_err!(c, p, msg); return nothing)
    end
    return x
end

_typename(v) = v isa AbstractString ? "a string" : v isa Bool ? "a boolean" : v isa Number ? "a number" :
               v isa AbstractVector ? "an array" : v isa AbstractDict ? "an object" : "null"

function _convert(c::_Ctx, p, v, kind)
    if kind === :int
        v isa Number && !(v isa Bool) && isfinite(v) && isinteger(v) && abs(v) < 2.0^62 && return Int(v)
        _err!(c, p, "must be an integer; got $(repr(v))")
    elseif kind === :uint64
        v isa Integer && !(v isa Bool) && 0 <= v <= typemax(UInt64) && return UInt64(v)
        v isa AbstractFloat && isinteger(v) && 0 <= v < 2.0^53 && return UInt64(v)
        _err!(c, p, "must be an integer in 0:2^64-1; got $(repr(v))")
    elseif kind === :float
        v isa Number && !(v isa Bool) && isfinite(v) && return Float64(v)
        _err!(c, p, "must be a finite number; got $(repr(v))")
    elseif kind === :bool
        v isa Bool && return v
        _err!(c, p, "must be true or false; got $(_typename(v))")
    elseif kind === :string
        v isa AbstractString && return String(v)
        _err!(c, p, "must be a string; got $(_typename(v))")
    elseif kind === :symbol
        v isa AbstractString && return Symbol(v)
        _err!(c, p, "must be a string; got $(_typename(v))")
    elseif kind === :ints || kind === :floats
        v isa AbstractVector || (_err!(c, p, "must be an array; got $(_typename(v))"); return nothing)
        out = kind === :ints ? Int[] : Float64[]
        for (i, x) in pairs(v)
            y = _convert(c, "$p[$i]", x, kind === :ints ? :int : :float)
            y === nothing && return nothing
            push!(out, y)
        end
        return out
    elseif kind === :matrix
        (v isa AbstractVector && all(r -> r isa AbstractVector, v)) ||
            (_err!(c, p, "must be an array of rows"); return nothing)
        k = length(v)
        all(r -> length(r) == k, v) || (_err!(c, p, "must be square ($k rows)"); return nothing)
        M = zeros(k, k)
        for i = 1:k, j = 1:k
            y = _convert(c, "$p[$i][$j]", v[i][j], :float)
            y === nothing && return nothing
            M[i, j] = y
        end
        return M
    else
        error("unknown kind $kind")
    end
    return nothing
end

_positive(x) = x > 0 ? nothing : "must be positive"
_atleast(n) = x -> x >= n ? nothing : "must be at least $n"
_nonneg(x) = x >= 0 ? nothing : "must be non-negative"
_maf(x) = 0 <= x < 0.5 ? nothing : "must lie in [0, 0.5)"

function _object(c::_Ctx, d, path, key; required=true)
    p = _join(path, key)
    haskey(d, key) || (required && _err!(c, p, "required"); return nothing)
    v = d[key]
    v === nothing && !required && return nothing
    _isobject(v) || (_err!(c, p, "must be an object; got $(_typename(v))"); return nothing)
    return v
end

# The scheme: JSON key, keyword of Scheme, kind, default, extra options.
const _SCHEME_FIELDS = [
    ("noff", :noff, :int, _NONE, (;)),
    ("nsire", :nsire, :int, _NONE, (;)),
    ("ndam", :ndam, :int, _NONE, (;)),
    ("evaluation", :evaluation, :symbol, :pblup, (choices=[:tbv, :phenotype, :pblup, :gblup, :ssgblup],)),
    ("selection", :selection, :symbol, :truncation, (choices=[:truncation, :ocs, :rocs, :random],)),
    ("K", :K, :float, NaN, (nullable=true, null_as=NaN)),
    ("lambda", :λ, :float, 0.0, (;)),
    ("kappa", :κ, :float, 0.0, (;)),
    ("ub", :ub, :float, 1.0, (;)),
    ("blend", :blend, :float, 0.05, (;)),
    ("max_sire_age", :max_sire_age, :int, 0, (;)),
    ("max_dam_age", :max_dam_age, :int, 0, (;)),
    ("reml", :reml, :bool, false, (;)),
    ("weights", :weights, :floats, Float64[], (;)),
    ("survival", :survival, :floats, Float64[], (;)),
    ("cull_unselected", :cull_unselected, :bool, false, (;)),
    ("genotyped_cohorts", :genotyped_cohorts, :int, 0, (;)),
    ("ocs", :ocs, :symbol, :fixed, (choices=[:fixed, :cohort, :age_year, :age_generation],)),
    ("delta_F", :ΔF, :float, NaN, (nullable=true, null_as=NaN)),
    ("record_meiosis", :record_meiosis, :bool, false, (;)),
    ("window", :window, :int, 0, (;)),
    ("ancestor_generations", :ancestor_generations, :int, typemax(Int), (nullable=true, null_as=typemax(Int))),
    ("retain", :retain, :symbol, :all, (choices=[:all, :active],)),
    ("report_reliability", :report_reliability, :bool, false, (;)),
]

const _TOP_KEYS = ["schema", "schema_version", "name", "seed", "replicates", "steps", "species", "founders", "loci",
    "traits", "genetic_correlation", "residual_correlation", "scheme"]

"""
    parse_scenario(json::AbstractString) -> Scenario
    parse_scenario(d::AbstractDict) -> Scenario

Read a scenario from JSON text (or a parsed object) and validate it. Throws
a [`ScenarioError`](@ref) that lists every problem found.
"""
function parse_scenario(text::AbstractString)
    d = try
        JSON.parse(text; dicttype=_StrictDict)
    catch e
        msg = e isa ArgumentError ? e.msg : sprint(showerror, e)
        throw(ScenarioError(["not valid JSON: " * first(split(msg, '\n'))]))
    end
    return parse_scenario(d)
end

function parse_scenario(d::AbstractDict)
    c = _Ctx(String[])
    s = _scenario(c, d)
    isempty(c.errors) || throw(ScenarioError(c.errors))
    return s
end
parse_scenario(x) = throw(ScenarioError(["the scenario must be a JSON object; got $(_typename(x))"]))

"""
    validate_scenario(json) -> Vector{String}

The problems of a scenario (JSON text or a parsed object), each prefixed
by its path; empty if it is valid.
"""
function validate_scenario(x)
    try
        x isa AbstractString ? parse_scenario(x) : parse_scenario(x)
        return String[]
    catch e
        e isa ScenarioError && return e.errors
        rethrow()
    end
end

"""
    read_scenario(path) -> Scenario

Read and validate the scenario in the JSON file `path`.
"""
read_scenario(path::AbstractString) = parse_scenario(read(path, String))

function _scenario(c::_Ctx, d)
    _isobject(d) || (_err!(c, "", "the scenario must be a JSON object"); return nothing)
    _known!(c, d, "", _TOP_KEYS)
    schema = _get(c, d, "", "schema", :string)
    schema === nothing || schema == SCENARIO_SCHEMA ||
        _err!(c, "schema", "must be \"$SCENARIO_SCHEMA\"; got \"$schema\"")
    version = _get(c, d, "", "schema_version", :int)
    version === nothing || version == SCENARIO_SCHEMA_VERSION || _err!(c, "schema_version",
        "this BnGsim reads version $SCENARIO_SCHEMA_VERSION; got $version")
    name = _get(c, d, "", "name", :string)
    seed = _get(c, d, "", "seed", :uint64)
    reps = _get(c, d, "", "replicates", :int; default=1, check=_positive)
    steps = _get(c, d, "", "steps", :int; check=_positive)

    sp = _species(c, _object(c, d, "", "species"))
    fo = _founders(c, _object(c, d, "", "founders"))
    traits = _traits(c, d)
    k = traits === nothing ? 0 : length(traits)
    lo = _loci(c, _object(c, d, "", "loci"), k)
    G = _correlation(c, d, "genetic_correlation", k)
    R = _correlation(c, d, "residual_correlation", k)
    sc = _scheme(c, _object(c, d, "", "scheme"), k, lo)

    if fo !== nothing && traits !== nothing
        fo.nid > k || _err!(c, "founders.nid", "must exceed the number of traits ($k)")
    end
    isempty(c.errors) || return nothing
    return Scenario(name, seed, reps, steps, sp, fo, lo.qtl, lo.chip, lo.overlap, traits, G, R, sc)
end

function _species(c, d)
    d === nothing && return nothing
    path = "species"
    if haskey(d, "preset")
        _known!(c, d, path, ["preset", "M", "morgan"])
        preset = _get(c, d, path, "preset", :string; choices=sort!(collect(keys(_PRESETS))))
        preset === nothing && return nothing
        base = _PRESETS[preset]()
        name, chr = base.name, Int.(base.chromosome)
    else
        _known!(c, d, path, ["name", "chromosome", "M", "morgan"])
        name = _get(c, d, path, "name", :string)
        chr = _get(c, d, path, "chromosome", :ints;
            check=v -> isempty(v) ? "must not be empty" : all(>(0), v) ? nothing : "lengths must be positive")
        (name === nothing || chr === nothing) && return nothing
        length(chr) <= typemax(Int8) || (_err!(c, "$path.chromosome", "at most 127 chromosomes"); return nothing)
        sum(chr) < 2^32 || (_err!(c, "$path.chromosome", "the genome must be shorter than 2^32 bp"); return nothing)
    end
    haskey(d, "M") && haskey(d, "morgan") && (_err!(c, path, "give M or morgan, not both"); return nothing)
    morgan = if haskey(d, "morgan")
        m = _get(c, d, path, "morgan", :floats;
            check=v -> all(>(0), v) ? nothing : "genetic lengths must be positive")
        m === nothing && return nothing
        length(m) == length(chr) ||
            (_err!(c, "$path.morgan", "needs one length per chromosome ($(length(chr)))"); return nothing)
        m
    else
        M = _get(c, d, path, "M", :float; default=1e8, check=_positive)
        M === nothing && return nothing
        chr ./ M
    end
    return (name=name, chromosome=chr, morgan=morgan)
end

function _founders(c, d)
    d === nothing && return nothing
    path = "founders"
    _known!(c, d, path, ["ne", "nt", "mr", "decline", "nid", "shared_history", "herds"])
    ne = _get(c, d, path, "ne", :int; check=_atleast(2))
    nt = _get(c, d, path, "nt", :int; check=_positive)
    mr = _get(c, d, path, "mr", :float; default=1.0,
        check=x -> 0.01 < x < 20 ? nothing : "must lie in (0.01, 20), FisherWright's bounds")
    decline = _get(c, d, path, "decline", :ints; default=Int[],
        check=v -> all(>=(2), v) ? nothing : "every size must be at least 2")
    nid = _get(c, d, path, "nid", :int; check=_atleast(2))
    shared = _get(c, d, path, "shared_history", :bool; default=true)
    herds = _get(c, d, path, "herds", :int; default=1, check=_positive)
    any(isnothing, (ne, nt, mr, decline, nid, shared, herds)) && return nothing
    last = isempty(decline) ? ne : decline[end]
    nid <= last || (_err!(c, "$path.nid", "exceeds the $last individuals of the last historical generation");
        return nothing)
    return (ne=ne, nt=nt, mr=mr, decline=decline, nid=nid, shared_history=shared, herds=herds)
end

function _panel(c, d, path)
    _known!(c, d, path, ["nlc", "maf"])
    nlc = _get(c, d, path, "nlc", :int; check=_positive)
    maf = _get(c, d, path, "maf", :float; default=0.0, check=_maf)
    (nlc === nothing || maf === nothing) && return nothing
    return (nlc=nlc, maf=maf)
end

function _loci(c, d, k)
    bad = (qtl=(nlc=0, maf=0.0), chip=nothing, overlap=:disjoint)
    d === nothing && return bad
    path = "loci"
    _known!(c, d, path, ["qtl", "chip", "overlap"])
    q = _object(c, d, path, "qtl")
    qtl = q === nothing ? nothing : _panel(c, q, "$path.qtl")
    ch = _object(c, d, path, "chip"; required=false)
    chip = ch === nothing ? nothing : _panel(c, ch, "$path.chip")
    overlap = _get(c, d, path, "overlap", :symbol; default=:disjoint, choices=[:disjoint, :allowed])
    qtl !== nothing && k > 0 && qtl.nlc < k && _err!(c, "$path.qtl.nlc", "needs at least as many QTL as traits ($k)")
    return (qtl=something(qtl, bad.qtl), chip=chip, overlap=something(overlap, :disjoint))
end

function _traits(c, d)
    haskey(d, "traits") || (_err!(c, "traits", "required"); return nothing)
    v = d["traits"]
    v isa AbstractVector || (_err!(c, "traits", "must be an array"); return nothing)
    isempty(v) && (_err!(c, "traits", "needs at least one trait"); return nothing)
    out = []
    for (i, t) in pairs(v)
        path = "traits[$i]"
        _isobject(t) || (_err!(c, path, "must be an object"); continue)
        _known!(c, t, path, ["name", "h2", "sigma_a", "mu", "sex", "age", "da", "c2"])
        name = _get(c, t, path, "name", :string; check=x -> occursin(r"^[[:alpha:]_\p{L}][\p{L}\p{N}_]*$", x) ?
            nothing : "must be a valid identifier (letters, digits, _; not starting with a digit)")
        h2 = _get(c, t, path, "h2", :float; default=0.25, check=x -> 0 < x <= 1 ? nothing : "must lie in (0, 1]")
        sa = _get(c, t, path, "sigma_a", :float; default=1.0, check=_positive)
        mu = _get(c, t, path, "mu", :float; default=0.0)
        sex = _get(c, t, path, "sex", :int; default=2, check=x -> x in 0:2 ? nothing : "must be 0 (females), 1 (males) or 2 (both)")
        age = _get(c, t, path, "age", :int; default=0, check=_nonneg)
        da = _da(c, t, path)
        c2 = _get(c, t, path, "c2", :float; default=0.0, check=x -> 0 <= x < 1 ? nothing : "must lie in [0, 1)")
        any(isnothing, (name, h2, sa, mu, sex, age, da, c2)) && continue
        push!(out, (name=name, h2=h2, sigma_a=sa, mu=mu, sex=sex, age=age, da=da, c2=c2))
    end
    names = [t.name for t in out]
    allunique(names) || _err!(c, "traits", "trait names must be unique")
    length(out) == length(v) || return nothing
    return [t for t in out]
end

function _da(c, t, path)
    haskey(t, "da") || return (family=:normal, a=0.0, b=1.0, c=0.0, sign=:positive)
    p = "$path.da"
    d = t["da"]
    _isobject(d) || (_err!(c, p, "must be an object"); return nothing)
    fam = _get(c, d, p, "family", :symbol; choices=[:normal, :gamma, :laplace, :spike_slab])
    fam === nothing && return nothing
    x, sign = 0.0, :positive
    if fam === :normal
        _known!(c, d, p, ["family", "mean", "sd"])
        a = _get(c, d, p, "mean", :float; default=0.0)
        b = _get(c, d, p, "sd", :float; default=1.0, check=_positive)
    elseif fam === :gamma
        _known!(c, d, p, ["family", "shape", "scale", "sign"])
        a = _get(c, d, p, "shape", :float; check=_positive)
        b = _get(c, d, p, "scale", :float; default=1.0, check=_positive)
        sign = _get(c, d, p, "sign", :symbol; default=:positive, choices=[:positive, :random])
        sign === nothing && return nothing
    elseif fam === :laplace
        _known!(c, d, p, ["family", "scale"])
        a = 0.0
        b = _get(c, d, p, "scale", :float; default=1.0, check=_positive)
    else
        _known!(c, d, p, ["family", "pi", "sd", "spike_sd"])
        a = _get(c, d, p, "pi", :float; check=x -> 0 < x <= 1 ? nothing : "must lie in (0, 1]")
        b = _get(c, d, p, "sd", :float; default=1.0, check=_positive)
        x = _get(c, d, p, "spike_sd", :float; default=0.0, check=x -> x >= 0 ? nothing : "must not be negative")
        x === nothing && return nothing
    end
    (a === nothing || b === nothing) && return nothing
    return (family=fam, a=a, b=b, c=x, sign=sign)
end

function _correlation(c, d, key, k)
    k == 0 && return zeros(0, 0)
    haskey(d, key) || return Matrix{Float64}(I, k, k)
    M = _get(c, d, "", key, :matrix)
    M === nothing && return nothing
    size(M, 1) == k || (_err!(c, key, "must be $k × $k (one row per trait)"); return nothing)
    issymmetric(M) || (_err!(c, key, "must be symmetric"); return nothing)
    all(M[i, i] == 1 for i = 1:k) || (_err!(c, key, "must have a unit diagonal"); return nothing)
    isposdef(M) || (_err!(c, key, "must be positive definite"); return nothing)
    return M
end

function _scheme(c, d, k, lo)
    d === nothing && return nothing
    path = "scheme"
    _known!(c, d, path, [f[1] for f in _SCHEME_FIELDS])
    vals = Pair{Symbol,Any}[]
    ok = true
    for (key, kw, kind, default, opts) in _SCHEME_FIELDS
        v = _get(c, d, path, key, kind; default=default, opts...)
        v === nothing && (ok = false; continue)
        push!(vals, kw => v)
    end
    ok || return nothing
    nt = (; vals...)
    isempty(nt.weights) || length(nt.weights) == k ||
        (_err!(c, "$path.weights", "needs one weight per trait ($k)"); ok = false)
    nt.evaluation in (:gblup, :ssgblup) && lo.chip === nothing &&
        (_err!(c, "$path.evaluation", ":$(nt.evaluation) needs loci.chip"); ok = false)
    nt.max_sire_age >= 0 && nt.max_dam_age >= 0 || (_err!(c, path, "maximum ages must be non-negative"); ok = false)
    ok || return nothing
    # The remaining combinations are checked by Scheme itself.
    try
        Scheme(; nt..., chip=lo.chip === nothing ? Int32[] : Int32[1])
    catch e
        e isa ArgumentError || rethrow()
        _err!(c, path, e.msg)
        return nothing
    end
    return nt
end

# ---- writing ---------------------------------------------------------------

_jnum(x::Float64) = isnan(x) ? nothing : x
_rows(M) = [collect(M[i, :]) for i = 1:size(M, 1)]

"""
    scenario_dict(s::Scenario) -> OrderedDict-like JSON object

The canonical JSON object of `s`: every key, with defaults filled.
"""
function scenario_dict(s::Scenario)
    O = JSON.Object{String,Any}
    sp = O("name" => s.species.name, "chromosome" => s.species.chromosome, "morgan" => s.species.morgan)
    f = s.founders
    fo = O("ne" => f.ne, "nt" => f.nt, "mr" => f.mr, "decline" => f.decline, "nid" => f.nid,
        "shared_history" => f.shared_history, "herds" => f.herds)
    panel(p) = p === nothing ? nothing : O("nlc" => p.nlc, "maf" => p.maf)
    lo = O("qtl" => panel(s.qtl), "chip" => panel(s.chip), "overlap" => String(s.overlap))
    da(d) = d.family === :normal ? O("family" => "normal", "mean" => d.a, "sd" => d.b) :
            d.family === :gamma ? O("family" => "gamma", "shape" => d.a, "scale" => d.b, "sign" => String(d.sign)) :
            d.family === :laplace ? O("family" => "laplace", "scale" => d.b) :
            O("family" => "spike_slab", "pi" => d.a, "sd" => d.b, "spike_sd" => d.c)
    tr = [O("name" => t.name, "h2" => t.h2, "sigma_a" => t.sigma_a, "mu" => t.mu, "sex" => t.sex, "age" => t.age,
        "da" => da(t.da), "c2" => t.c2) for t in s.traits]
    sc = O()
    for (key, kw, kind, _, _) in _SCHEME_FIELDS
        v = getfield(s.scheme, kw)
        sc[key] = kind === :symbol ? String(v) : kind === :float ? _jnum(v) :
                  kw === :ancestor_generations && v == typemax(Int) ? nothing : v
    end
    return O("schema" => SCENARIO_SCHEMA, "schema_version" => SCENARIO_SCHEMA_VERSION, "name" => s.name,
        "seed" => s.seed, "replicates" => s.replicates, "steps" => s.steps, "species" => sp, "founders" => fo,
        "loci" => lo, "traits" => tr, "genetic_correlation" => _rows(s.genetic_correlation),
        "residual_correlation" => _rows(s.residual_correlation), "scheme" => sc)
end

"""
    scenario_json(s::Scenario) -> String

The canonical JSON text of `s` (two-space indentation, every key present).
Reading it back gives a scenario equal to `s`, and writing that gives the
same text.
"""
scenario_json(s::Scenario) = JSON.json(scenario_dict(s); pretty=2) * "\n"

"""
    write_scenario(path, s::Scenario)

Write the canonical JSON of `s` to `path`.
"""
write_scenario(path::AbstractString, s::Scenario) = write(path, scenario_json(s))

# ---- running ---------------------------------------------------------------

_species_of(s::Scenario) = GenericSpecies(s.species.name, s.species.chromosome; morgan=s.species.morgan)

# Stages of the founder domain.
const _STAGE_HISTORY, _STAGE_DECLINE, _STAGE_SAMPLE, _STAGE_NUCLEOTIDES = 0x0001, 0x0002, 0x0003, 0x0004

# Placeholder nucleotides (`ref` and a different `alt`) from `rng`. The map
# of `founder_haplotypes` draws them from the global generator, so they
# would differ between runs of the same scenario; nothing else depends on
# them.
function _nucleotides!(vm::VariantMap, rng::AbstractRNG)
    bases = ('A', 'C', 'G', 'T')
    for i in eachindex(vm.ref)
        r = rand(rng, 1:4)
        vm.ref[i] = bases[r]
        vm.alt[i] = bases[mod1(r + rand(rng, 1:3), 4)]
    end
    return vm
end

# The historical population of replicate `r` (0: shared by all).
function _history(s::Scenario, sp::Species, r::Integer)
    f = s.founders
    key = SimKey(s.seed; replicate=r, domain=Domain.FOUNDERS)
    res = fisher_wright(sp, f.ne, f.nt, f.mr; result=true, rng=keyed_rng(SimKey(key; stage=_STAGE_HISTORY)))
    isempty(f.decline) && return res
    return fisher_wright(f.ne, length(f.decline), s.species.chromosome, f.mr; morgan=s.species.morgan,
        init=res, popsize=f.decline, result=true, rng=keyed_rng(SimKey(key; stage=_STAGE_DECLINE)))
end

_dist(d) = d.family === :normal ? Normal(d.a, d.b) :
           d.family === :gamma ? (d.sign === :random ? RandomSign(Gamma(d.a, d.b)) : Gamma(d.a, d.b)) :
           d.family === :laplace ? Laplace(0.0, d.b) : SpikeSlab(d.a, d.b; spike_sd=d.c)

"""
    run_replicate(s::Scenario, r = 1; history = nothing)
        -> (; population, summary, species, map, qtl, chip)

Run replicate `r` of scenario `s`: the historical population (or
`history`, a `FisherWrightResult` to start from), the founders, the QTL and
chip loci, the traits, then [`simulate!`](@ref) for `s.steps` steps with
the [`Scheme`](@ref) of the scenario. Every random draw comes from a key of
`SimKey(s.seed; replicate = r)` (the historical population from replicate
`0` when it is shared), so a replicate does not depend on the others or on
the number of threads.
"""
function run_replicate(s::Scenario, r::Integer=1; history=nothing)
    r >= 1 || throw(ArgumentError("replicates are numbered from 1"))
    rep = _setup_replicate(s, r, history)
    summary = simulate!(rep.population, rep.scheme, rep.meiosis, s.steps)
    return (; population=rep.population, summary, species=rep.species, map=rep.map, qtl=rep.qtl, chip=rep.chip)
end

# Replicate `r` of `s` before its first step: the founders, loci and traits,
# the scheme and the meiosis map.
function _setup_replicate(s::Scenario, r::Integer, history)
    sp = _species_of(s)
    h = history !== nothing ? history : _history(s, sp, s.founders.shared_history ? 0 : r)
    key = SimKey(s.seed; replicate=r)
    xy, lmp = founder_haplotypes(h, s.founders.nid;
        rng=keyed_rng(SimKey(key; domain=Domain.FOUNDERS, stage=_STAGE_SAMPLE)))
    hap = Haplotype(xy)
    vm = _to_variant_map(lmp)
    _nucleotides!(vm, keyed_rng(SimKey(key; domain=Domain.FOUNDERS, stage=_STAGE_NUCLEOTIDES)))
    sets = [SNPSet("qtl", s.qtl.nlc; maf=s.qtl.maf, exclusive=true)]
    s.chip === nothing || push!(sets, SNPSet("chip", s.chip.nlc; maf=s.chip.maf, exclusive=s.overlap === :disjoint))
    panels = sample_loci(hap, vm, sets...; rng=keyed_rng(SimKey(key; domain=Domain.GENOTYPING)))
    ts = [Trait(t.name; sex=t.sex, age=t.age, h²=t.h2, μ=t.mu, σₐ=t.sigma_a, da=_dist(t.da)) for t in s.traits]
    qrng = keyed_rng(SimKey(key; domain=Domain.QTL))
    pop = if length(ts) == 1
        founder_population(hap, eqtl(hap, panels[1], ts[1]; rng=qrng), ts[1], key;
            herds=s.founders.herds, c2=s.traits[1].c2)
    else
        founder_population(hap, eqtl(hap, panels[1], s.genetic_correlation, ts; rng=qrng), ts, key;
            Re=s.residual_correlation, herds=s.founders.herds, c2=[t.c2 for t in s.traits])
    end
    chip = s.chip === nothing ? Int32[] : Int32.(panels[2].loci)
    sc = Scheme(; s.scheme..., chip=chip)
    return (; population=pop, scheme=sc, meiosis=MeiosisMap(sp, vm), species=sp, map=vm,
        qtl=Int32.(panels[1].loci), chip)
end

"""
    run_scenario(s::Scenario; replicates = 1:s.replicates) -> DataFrame
    run_scenario(path::AbstractString; kw...) -> DataFrame

Run the replicates of `s` and return their per-step summaries
([`simulate!`](@ref)) in one table, with a first column `replicate`. A
shared historical population is simulated once.
"""
function run_scenario(s::Scenario; replicates=1:s.replicates)
    all(r -> 1 <= r <= s.replicates, replicates) ||
        throw(ArgumentError("replicates must lie in 1:$(s.replicates)"))
    h = s.founders.shared_history ? _history(s, _species_of(s), 0) : nothing
    parts = map(collect(replicates)) do r
        df = run_replicate(s, r; history=h).summary
        insertcols!(df, 1, :replicate => fill(r, nrow(df)))
        df
    end
    return reduce(vcat, parts)
end
run_scenario(path::AbstractString; kw...) = run_scenario(read_scenario(path); kw...)
