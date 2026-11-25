# API Reference

## What 1.0 promises

From 1.0 on, these do not change incompatibly within a major version:

- the names exported by `BnGsim` (below), with their documented behaviour;
- the command line: interface 1 and its exit codes ("Command line");
- the scenario schema, version 1, under the policy of "Scenarios";
- the output files, `outputs_version` 1 ("Output files");
- the reproducibility invariants: results independent of the number of
  threads and of stopping and resuming, and `results_sha256`.

Outside the promise: names that are documented here but not exported
(reach them as `BnGsim.name`; for example `BnGsim.philox4x64`,
`BnGsim.keyed_int` and `BnGsim.transport`, unexported in 0.38), the
`BnGsim.Legacy` module (its crossover sampler is
`BnGsim.Legacy.crossovers`, separate from the exported `crossovers`), the coefficients of [`resource_estimate`](@ref)
(recalibrated as the code changes), and numerical results of a scenario
across BnGsim versions (each run records its versions and hashes).

Types of BnGStructs.jl appear in this API (`Haplotype`, `VariantMap`,
`LocusSet`, for example in [`founder_haplotypes`](@ref) and
[`genomes`](@ref)). They are promised as far as BnGsim's compat bound on
BnGStructs allows: a breaking release of BnGStructs can change them, and
BnGsim then states the change in its release notes.

```@autodocs
Modules = [BnGsim, BnGsim.Domain, BnGsim.ExitCode, BnGsim.FisherWrightAdapter]
```

## Legacy building blocks (0.1–0.3)

`BnGsim.Legacy` holds the building blocks of BnGsim 0.1–0.3. They draw
from an explicit `rng` or the global generator, not from the keyed
streams of the breeding loop, so they are outside the reproducibility
invariants and the 1.0 promise. Load them with `using BnGsim.Legacy`.

```@autodocs
Modules = [BnGsim.Legacy]
```
