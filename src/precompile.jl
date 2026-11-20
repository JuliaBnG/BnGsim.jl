# Compile the scenario and command-line paths into the package image, so
# that `bngsim run` does not spend its first half-minute compiling. The
# workload runs two tiny scenarios (evaluation by phenotype and by PBLUP)
# through `cli`, writing to a temporary directory.
@setup_workload begin
    doc(evaluation) = """
    {"schema": "bngsim-scenario", "schema_version": 1, "name": "precompile", "seed": 1, "steps": 2,
     "species": {"name": "toy", "chromosome": [2000000, 1000000]},
     "founders": {"ne": 30, "nt": 60, "mr": 5.0, "nid": 20},
     "loci": {"qtl": {"nlc": 10}},
     "traits": [{"name": "y", "h2": 0.3}],
     "scheme": {"noff": 20, "nsire": 3, "ndam": 8, "evaluation": "$evaluation"}}
    """
    @compile_workload begin
        dir = mktempdir()
        for ev in ("phenotype", "pblup")
            path = joinpath(dir, "$ev.json")
            write(path, doc(ev))
            cli(["validate", path]; stdout=devnull, stderr=devnull)
            cli(["run", path, "--out", joinpath(dir, ev)]; stdout=devnull, stderr=devnull)
        end
        cli(["version", "--json"]; stdout=devnull, stderr=devnull)
        validate_scenario("{\"schema\": 1}")
        rm(dir; recursive=true)
    end
end
