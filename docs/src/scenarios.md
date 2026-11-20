# Scenarios

A scenario is a JSON file that describes a whole simulation: the species,
the historical population and the founders sampled from it, the QTL and
chip loci, the traits, the breeding scheme, the number of time steps, the
master seed and the number of replicates. It is the input of the planned
`bngsim` command line and web service, and runs in Julia with
[`run_scenario`](@ref).

```julia
using BnGsim
s = read_scenario("examples/minimal.json")   # read and validate
df = run_scenario(s)                         # per-step summaries of every replicate
```

The smallest example, `examples/minimal.json`:

```@eval
using Markdown, BnGsim
Markdown.MD(Markdown.Code("json", read(joinpath(pkgdir(BnGsim), "examples", "minimal.json"), String)))
```

## Schema

The keys, their defaults and their meaning are listed in the docstring of
[`Scenario`](@ref) (schema `"bngsim-scenario"`, version 1). They mirror the
Julia API: the `scheme` object takes the keywords of [`Scheme`](@ref),
with ASCII names for the Unicode ones (`lambda`, `kappa`, `delta_F`), and a
trait takes those of `BnGStructs.Trait` (`h2`, `sigma_a`, `mu`, `sex`,
`age`, and `da` for the distribution of the QTL effects).

## Validation

[`read_scenario`](@ref) and [`parse_scenario`](@ref) check the whole
document before anything runs and throw a [`ScenarioError`](@ref) that
lists every problem with its path, rather than the first one;
[`validate_scenario`](@ref) returns the list. Unknown keys, repeated keys,
values of the wrong type and out-of-range values are errors, never
ignored or repaired. So are combinations: a `:gblup` scheme without chip
loci, weights that do not match the traits, a correlation matrix that is
not a positive definite correlation matrix, more founders than the last
historical generation holds, or a scheme that [`Scheme`](@ref) itself
refuses. For example, a scheme with a misspelt key and a negative seed:

```@example scenarios
using BnGsim
text = replace(read(joinpath(pkgdir(BnGsim), "examples", "minimal.json"), String),
    "\"seed\": 2026" => "\"seed\": -1", "\"nsire\"" => "\"nsires\"")
validate_scenario(text)
```

The number of loci that pass the `maf` filters is only known once the
founders are sampled; a panel larger than that fails then, with an
`ArgumentError`.

## Canonical form

A [`Scenario`](@ref) holds every value, with the defaults filled.
[`scenario_json`](@ref) and [`write_scenario`](@ref) write it with every
key present, so the file states the whole simulation; reading it back
gives an equal scenario, and writing that gives the same text.

```@example scenarios
s = read_scenario(joinpath(pkgdir(BnGsim), "examples", "minimal.json"))
print(scenario_json(s))
```

## Runs and reproducibility

[`run_replicate`](@ref) runs one replicate and returns its population and
summary; [`run_scenario`](@ref) runs replicates `1:replicates` and stacks
their summaries with a column `replicate`. Every random draw comes from a
[`SimKey`](@ref) of the scenario's `seed`: the historical population from
replicate `0` when it is shared (`founders.shared_history`, the default),
then for each replicate `r` the founders (`Domain.FOUNDERS`), the chip and
QTL loci (`Domain.GENOTYPING`), the QTL effects (`Domain.QTL`) and the
breeding loop. A replicate is therefore the same whether it runs alone or
with others, and on any number of threads.

## Limits of version 1

- Founders come from a FisherWright historical population; a scenario
  takes no file paths (the service accepts no user paths, roadmap I9).
- Traits are additive (as the breeding loop).
- One scheme for the whole run.
- The scenario records no provenance yet; the command line will write it
  with the outputs (roadmap I2).
