# Scenario schema: reading, validation, canonical JSON and runs.

const MINIMAL = joinpath(@__DIR__, "..", "examples", "minimal.json")

# A valid two-trait scenario with a chip, small enough to run in seconds.
function scenario_doc()
    return Dict{String,Any}(
        "schema" => "bngsim-scenario", "schema_version" => 1, "name" => "two-traits", "seed" => 7,
        "replicates" => 2, "steps" => 2,
        "species" => Dict{String,Any}("name" => "toy", "chromosome" => [20_000_000, 10_000_000],
            "morgan" => [0.3, 0.1]),
        "founders" => Dict{String,Any}("ne" => 100, "nt" => 400, "decline" => [80, 60], "nid" => 40),
        "loci" => Dict{String,Any}("qtl" => Dict{String,Any}("nlc" => 40, "maf" => 0.02),
            "chip" => Dict{String,Any}("nlc" => 200, "maf" => 0.05)),
        "traits" => [Dict{String,Any}("name" => "a", "h2" => 0.3),
            Dict{String,Any}("name" => "b", "h2" => 0.5, "sex" => 0, "da" => Dict{String,Any}("family" => "gamma",
                "shape" => 0.5))],
        "genetic_correlation" => [[1.0, -0.3], [-0.3, 1.0]],
        "scheme" => Dict{String,Any}("noff" => 40, "nsire" => 4, "ndam" => 12, "evaluation" => "gblup",
            "weights" => [1.0, 0.5]))
end

# A copy of the document with the value at `path` (a vector of keys and
# indices) replaced, or deleted with `delete`.
function edited(path, value; delete=false)
    d = deepcopy(scenario_doc())
    x = d
    for k in path[1:end-1]
        x = x[k]
    end
    delete ? Base.delete!(x, path[end]) : (x[path[end]] = value)
    return d
end

@testset "Scenario: reading, defaults and canonical JSON" begin
    s = read_scenario(MINIMAL)
    @test s isa Scenario && s.name == "minimal" && s.seed == 2026 && s.replicates == 1 && s.steps == 3
    @test s.species.chromosome == [20_000_000, 10_000_000] && s.species.morgan == [0.2, 0.1]
    @test isequal(s.founders, (ne=100, nt=500, mr=1.0, decline=Int[], nid=50, shared_history=true, herds=1,
        popsize=nothing, sexratio=NaN, selection=nothing))
    @test s.qtl == (nlc=50, maf=0.01, from=:neutral) && s.chip === nothing && s.overlap === :disjoint
    @test s.traits[1] == (name="y", h2=0.3, sigma_a=1.0, mu=0.0, sex=2, age=0, da=(family=:normal, a=0.0, b=1.0, c=0.0, sign=:positive), c2=0.0, records=1,
        repeatability=0.3, missing=0.0, unrecorded_herds=0.0)
    @test s.genetic_correlation == ones(1, 1) && s.residual_correlation == ones(1, 1)
    # Scheme defaults are those of Scheme itself.
    d = Scheme(; noff=50, nsire=5, ndam=20)
    for f in keys(s.scheme)
        f in (:evaluation,) || @test isequal(getfield(s.scheme, f), getfield(d, f))
    end
    @test s.scheme.evaluation === :phenotype && isnan(s.scheme.K) && s.scheme.ancestor_generations == typemax(Int)

    # The canonical JSON states every key; it reads back equal, and writes
    # the same text again.
    j = scenario_json(s)
    @test occursin("\"ancestor_generations\": null", j) && occursin("\"K\": null", j) && occursin("\"chip\": null", j)
    s2 = parse_scenario(j)
    @test s2 == s && scenario_json(s2) == j && hash(s2) == hash(s)
    path = tempname()
    write_scenario(path, s)
    @test read_scenario(path) == s
    rm(path)
    t = parse_scenario(scenario_doc())
    @test parse_scenario(scenario_json(t)) == t && ntraits(t) == 2
    @test t.traits[2].da == (family=:gamma, a=0.5, b=1.0, c=0.0, sign=:positive) && t.genetic_correlation == [1.0 -0.3; -0.3 1.0]
    # Presets resolve to the species' genome.
    c = parse_scenario(edited(["species"], Dict{String,Any}("preset" => "Cattle")))
    @test c.species.chromosome == Int.(Cattle().chromosome) && c.species.morgan == Cattle().morgan
    @test parse_scenario(scenario_json(c)) == c
    @test validate_scenario(scenario_json(t)) == String[]
end

@testset "Scenario: malformed manifests" begin
    # Each case: the document (or text) and a fragment of the error expected.
    cases = Any[
        ("{\"schema\": ", "not valid JSON"),
        ("[1, 2]", "must be a JSON object"),
        ("{\"schema\": \"bngsim-scenario\", \"schema\": \"x\"}", "duplicate key \"schema\""),
        (edited(["schema"], nothing; delete=true), "schema: required"),
        (edited(["schema"], "other"), "schema: must be \"bngsim-scenario\""),
        (edited(["schema_version"], 2), "schema_version: this BnGsim reads version 1"),
        (edited(["seeds"], 1), "seeds: unknown key"),
        (edited(["seed"], -1), "seed: must be an integer in 0:2^64-1"),
        (edited(["seed"], 1.5), "seed: must be an integer"),
        (edited(["replicates"], 0), "replicates: must be positive"),
        (edited(["steps"], nothing; delete=true), "steps: required"),
        (edited(["steps"], "3"), "steps: must be an integer"),
        (edited(["species", "chromosome"], [100, -5]), "species.chromosome: lengths must be positive"),
        (edited(["species", "chromosome"], fill(3_000_000_000, 2)), "shorter than 2^32"),
        (edited(["species", "morgan"], [0.3]), "species.morgan: needs one length per chromosome"),
        (edited(["species", "M"], 1e8), "give M or morgan, not both"),
        (edited(["species"], Dict{String,Any}("preset" => "Yak")), "species.preset: must be one of"),
        (edited(["founders", "ne"], 1), "founders.ne: must be at least 2"),
        (edited(["founders", "mr"], 50.0), "founders.mr: must lie in (0.01, 20)"),
        (edited(["founders", "decline"], [80, 1]), "founders.decline: every size must be at least 2"),
        (edited(["founders", "nid"], 70), "founders.nid: exceeds the 60 individuals"),
        (edited(["loci", "qtl", "maf"], 0.5), "loci.qtl.maf: must lie in [0, 0.5)"),
        (edited(["loci", "qtl", "nlc"], 1), "loci.qtl.nlc: needs at least as many QTL as traits"),
        (edited(["loci", "overlap"], "some"), "loci.overlap: must be one of"),
        (edited(["loci", "chip"], nothing; delete=true), "scheme.evaluation: :gblup needs loci.chip"),
        (edited(["traits"], []), "traits: needs at least one trait"),
        (edited(["traits", 2, "name"], "a"), "traits: trait names must be unique"),
        (edited(["traits", 1, "name"], "1y"), "traits[1].name: must be a valid identifier"),
        (edited(["traits", 1, "h2"], 0.0), "traits[1].h2: must lie in (0, 1]"),
        (edited(["traits", 1, "sex"], 3), "traits[1].sex: must be 0 (females), 1 (males) or 2"),
        (edited(["traits", 1, "age"], 0.5), "traits[1].age: must be an integer"),
        (edited(["traits", 2, "da", "family"], "cauchy"), "traits[2].da.family: must be one of"),
        (edited(["traits", 2, "da", "shape"], nothing; delete=true), "traits[2].da.shape: required"),
        (edited(["genetic_correlation"], [[1.0, 0.3]]), "genetic_correlation: must be square"),
        (edited(["genetic_correlation"], [[1.0, 0.3], [0.2, 1.0]]), "genetic_correlation: must be symmetric"),
        (edited(["genetic_correlation"], [[1.0, 1.2], [1.2, 1.0]]), "genetic_correlation: must be positive definite"),
        (edited(["residual_correlation"], [[2.0, 0.0], [0.0, 1.0]]), "residual_correlation: must have a unit diagonal"),
        (edited(["scheme", "nsires"], 4), "scheme.nsires: unknown key"),
        (edited(["scheme", "evaluation"], "blup"), "scheme.evaluation: must be one of"),
        (edited(["scheme", "noff"], 40.5), "scheme.noff: must be an integer"),
        (edited(["scheme", "reml"], 1), "scheme.reml: must be true or false"),
        (edited(["scheme", "weights"], [1.0]), "scheme.weights: needs one weight per trait (2)"),
        (edited(["scheme", "selection"], "ocs"), "scheme: :ocs needs a constraint K > 0"),
        (edited(["scheme", "survival"], [1.5]), "scheme: survival probabilities must lie in [0, 1]"),
        (edited(["scheme", "window"], 1) |> d -> (d["scheme"]["max_sire_age"] = 2; d), "scheme: window = 1"),
        (edited(["scheme", "noff"], nothing), "scheme.noff: must not be null"),
    ]
    for (doc, msg) in cases
        errs = validate_scenario(doc)
        ok = any(e -> occursin(msg, e), errs)
        ok || @info "expected \"$msg\"" errs
        @test ok
        @test_throws ScenarioError (doc isa AbstractString ? parse_scenario(doc) : parse_scenario(doc))
    end
    @test length(cases) >= 30
    # Every problem of a document is reported at once.
    bad = edited(["steps"], 0)
    bad["seed"] = -1
    bad["scheme"]["nsires"] = 4
    errs = validate_scenario(bad)
    @test length(errs) == 3
    e = try
        parse_scenario(bad)
    catch err
        err
    end
    @test e isa ScenarioError && occursin("3 problems", sprint(showerror, e))
end

@testset "Scenario: runs" begin
    s = read_scenario(MINIMAL)
    df = run_scenario(s)
    @test nrow(df) == s.steps + 1 && all(==(1), df.replicate) && df.n == fill(50, 4)
    @test isequal(run_scenario(MINIMAL), df)     # from the file, the same run
    # Replicates are independent of each other and of the batch they run in.
    t = parse_scenario(scenario_doc())
    all2 = run_scenario(t)
    @test sort(unique(all2.replicate)) == [1, 2]
    one = run_replicate(t, 2)
    @test isequal(one.summary, all2[all2.replicate .== 2, 2:end])
    @test !isequal(all2[all2.replicate .== 1, 2:end], all2[all2.replicate .== 2, 2:end])
    @test length(one.chip) == 200 && isempty(intersect(one.chip, one.qtl)) && length(one.qtl) == 40
    @test ntraits(one.population) == 2 && length(one.population) == 40 + 2 * 40
    # Without a shared history each replicate runs its own.
    u = parse_scenario(edited(["founders", "shared_history"], false))
    @test !isequal(run_replicate(u, 1).summary, run_replicate(t, 1).summary)
    @test_throws ArgumentError run_scenario(t; replicates=[3])
end
