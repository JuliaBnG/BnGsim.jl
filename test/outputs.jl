# The contract of the results (package roadmap §11, decisions b and d):
# the results hash, which ignores keys at their defaults, and the version
# and columns of the output files (docs/src/outputs.md).

# Golden: changes only if the results-relevant content of minimal.json, or
# the way it is hashed, changes.
const GOLDEN_RESULTS_SHA = "1c06ae566ffef0f7e684ed402c2fdd3e68dd75fc3a99d8e37956ada4fcd2d287"

@testset "Results hash" begin
    s = read_scenario(MINIMAL)
    h = results_sha256(s)
    @test h == GOLDEN_RESULTS_SHA
    d = results_dict(s)
    @test !any(haskey(d, k) for k in ("schema", "schema_version", "name", "replicates"))
    @test d["scheme"] == Dict("noff" => 50, "nsire" => 5, "ndam" => 20, "evaluation" => "phenotype")
    txt = read(MINIMAL, String)
    variant(pairs...) = results_sha256(parse_scenario(replace(txt, pairs...)))
    # not results-relevant: the name, the number of replicates
    @test variant("\"name\": \"minimal\"" => "\"name\": \"other\"", "\"steps\": 3" => "\"steps\": 3, \"replicates\": 7") == h
    # keys written at their defaults (the canonical form of every key)
    @test results_sha256(parse_scenario(scenario_json(s))) == h
    @test variant("\"h2\": 0.3}" => "\"h2\": 0.3, \"c2\": 0.0, \"records\": 1, \"missing\": 0.0}",
        "\"M\": 100000000}" => "\"M\": 100000000, \"interference\": 1}") == h
    # results-relevant: any value off its default
    @test variant("\"seed\": 2026" => "\"seed\": 2027") != h
    @test variant("\"h2\": 0.3}" => "\"h2\": 0.3, \"c2\": 0.1}") != h
    @test variant("\"steps\": 3" => "\"steps\": 4") != h
    # a dependent default: a bottleneck's recovery is ne unless given
    bn(x) = variant("\"nid\": 50}" => "\"nid\": 50, \"popsize\": {\"bottleneck\": {\"size\": 50, \"start\": 400, \"duration\": 10$x}}}")
    @test bn("") == bn(", \"recovery\": 100") != bn(", \"recovery\": 200")
    # every example: the pruned scenario reads back as the scenario
    for f in readdir(joinpath(@__DIR__, "..", "examples"); join=true)
        endswith(f, ".json") || continue
        e = read_scenario(f)
        r = results_dict(e)
        for k in ("schema", "schema_version", "name", "replicates")
            r[k] = scenario_dict(e)[k]
        end
        @test parse_scenario(r) == e
    end
end

@testset "Output files: version and columns" begin
    @test BnGsim.OUTPUTS_VERSION == 1
    code, v, _ = run_cli("version", "--json")
    @test code == 0 && v["outputs_version"] == 1
    code, v, _ = run_cli("validate", MINIMAL)
    @test v["results_sha256"] == GOLDEN_RESULTS_SHA

    # merge: runs of scenarios that differ only in name and replicates
    dir = mktempdir()
    a, b = joinpath(dir, "a.json"), joinpath(dir, "b.json")
    write(a, replace(read(MINIMAL, String), "\"steps\": 3," => "\"steps\": 2, \"replicates\": 2,"))
    write(b, replace(read(a, String), "\"name\": \"minimal\"" => "\"name\": \"renamed\"", "\"replicates\": 2" => "\"replicates\": 3"))
    @test run_cli("run", a, "--out", joinpath(dir, "ra"), "--replicates", "1")[1] == 0
    @test run_cli("run", b, "--out", joinpath(dir, "rb"), "--replicates", "2")[1] == 0
    code, v, _ = run_cli("merge", joinpath(dir, "ra"), joinpath(dir, "rb"), "--out", joinpath(dir, "m"))
    @test code == 0 && v["replicates"] == [1, 2] && v["results_sha256"] == results_sha256(read_scenario(a))
    p = JSON.parse(read(joinpath(dir, "m", "provenance.json"), String))
    @test p["outputs_version"] == 1 && length(p["scenario_sha256"]) == 2
    pa = JSON.parse(read(joinpath(dir, "ra", "provenance.json"), String))
    @test pa["outputs_version"] == 1 && pa["results_sha256"] == v["results_sha256"]
    # the merged summary is that of one run of both replicates
    @test run_cli("run", a, "--out", joinpath(dir, "all"))[1] == 0
    @test read(joinpath(dir, "m", "summary.csv")) == read(joinpath(dir, "all", "summary.csv"))
    rm(dir; recursive=true)

    # the columns of summary.csv (and so of aggregate.csv) for outputs
    # version 1: a change here is a change of the output contract
    base(tr, sc) = """
    {"schema": "bngsim-scenario", "schema_version": 1, "name": "c", "seed": 1, "steps": 1, "replicates": 1,
     "species": {"name": "toy", "chromosome": [2000000, 1000000]},
     "founders": {"ne": 60, "nt": 60, "mr": 5.0, "nid": 40},
     "loci": {"qtl": {"nlc": 10}},
     "traits": $tr, "scheme": {"noff": 30, "nsire": 3, "ndam": 10, "evaluation": "pblup" $sc}}"""
    common = "replicate,history,generation,n,mean,var,genic,qtl_seg,qtl_fixed,qtl_lost,qtl_maf,qtl_alpha,qtl_top10"
    one = run_scenario(parse_scenario(base("""[{"name": "y", "h2": 0.3}]""", "")))
    @test join(names(one), ",") == common * ",F,coancestry,Fhet,Fgrm,alive,accuracy,reliability,nsire,ndam,L," *
                                   "intensity_sire,intensity_dam,va,ve"
    rep1 = run_scenario(parse_scenario(base("""[{"name": "y", "h2": 0.3, "records": 2}]""", "")))
    @test names(rep1)[end-2:end] == ["va", "vpe", "ve"]
    two = run_scenario(parse_scenario(base("""[{"name": "a", "h2": 0.3, "records": 2}, {"name": "b", "h2": 0.4}]""",
        ", \"selection\": \"ocs\", \"K\": 0.2, \"record_meiosis\": true")))
    @test join(names(two), ",") == common * ",mean_a,mean_b,var_a,var_b,rg_a_b,F,coancestry,Fhet,Fgrm,alive,Fseg," *
          "accuracy,accuracy_a,accuracy_b,reliability,reliability_a,reliability_b,nsire,ndam,L,intensity_sire," *
          "intensity_dam,va_a,va_b,ve_a,ve_b,vpe_a,ocs_bound,ocs_kinship,ocs_L,ocs_feasible"
    ag = names(aggregate(one))
    @test ag[1:5] == ["generation", "replicates", "n_mean", "n_sd", "n_se"]
    @test ag[end-14:end] == [c * s for c in ("dG", "dF", "dFhet", "var_ratio", "genic_ratio") for s in ("_mean", "_sd", "_se")]
end
