# Compile the scenario and command-line paths into the package image, so
# that `bngsim run` does not spend its first half-minute compiling. The
# workload runs tiny scenarios through `cli`, writing to a temporary
# directory: evaluation by phenotype and by PBLUP; two traits with repeated
# records, contemporary groups, optimum contributions and minimum-
# coancestry mating; and GBLUP on a chip.
@setup_workload begin
    doc(evaluation) = """
    {"schema": "bngsim-scenario", "schema_version": 1, "name": "precompile", "seed": 1, "steps": 2,
     "species": {"name": "toy", "chromosome": [2000000, 1000000]},
     "founders": {"ne": 30, "nt": 60, "mr": 5.0, "nid": 20},
     "loci": {"qtl": {"nlc": 10}},
     "traits": [{"name": "y", "h2": 0.3}],
     "scheme": {"noff": 20, "nsire": 3, "ndam": 8, "evaluation": "$evaluation"}}
    """
    features = """
    {"schema": "bngsim-scenario", "schema_version": 1, "name": "features", "seed": 1, "steps": 2,
     "species": {"name": "toy", "chromosome": [2000000, 1000000]},
     "founders": {"ne": 30, "nt": 60, "mr": 5.0, "nid": 20, "herds": 3},
     "loci": {"qtl": {"nlc": 10}},
     "traits": [{"name": "a", "h2": 0.3, "records": 3, "repeatability": 0.5, "c2": 0.1}, {"name": "b", "h2": 0.4}],
     "scheme": {"noff": 20, "nsire": 3, "ndam": 8, "evaluation": "pblup", "selection": "ocs", "K": 0.2,
                "mating": "min_coancestry", "max_dam_age": 1}}
    """
    genomic = """
    {"schema": "bngsim-scenario", "schema_version": 1, "name": "gblup", "seed": 1, "steps": 2,
     "species": {"name": "toy", "chromosome": [2000000, 1000000]},
     "founders": {"ne": 30, "nt": 60, "mr": 5.0, "nid": 20},
     "loci": {"qtl": {"nlc": 10}, "chip": {"nlc": 30}},
     "traits": [{"name": "y", "h2": 0.3}],
     "scheme": {"noff": 20, "nsire": 3, "ndam": 8, "evaluation": "gblup"}}
    """
    @compile_workload begin
        dir = mktempdir()
        for ev in ("phenotype", "pblup")
            path = joinpath(dir, "$ev.json")
            write(path, doc(ev))
            cli(["validate", path]; stdout=devnull, stderr=devnull)
            cli(["estimate", path]; stdout=devnull, stderr=devnull)
            cli(["run", path, "--out", joinpath(dir, ev)]; stdout=devnull, stderr=devnull)
        end
        for (name, extra) in (("features", features), ("gblup", genomic))
            path = joinpath(dir, "$name.json")
            write(path, extra)
            cli(["run", path, "--out", joinpath(dir, name)]; stdout=devnull, stderr=devnull)
        end
        cli(["version", "--json"]; stdout=devnull, stderr=devnull)
        validate_scenario("{\"schema\": 1}")
        rm(dir; recursive=true)
    end
end
