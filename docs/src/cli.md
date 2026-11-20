# Command line

`bin/bngsim` runs a [scenario](scenarios.md) from the shell. It is the
interface the planned web service will call, so its arguments, output and
exit codes are versioned (CLI interface 1) separately from the package.

```text
bngsim validate SCENARIO.json
bngsim estimate SCENARIO.json
bngsim run SCENARIO.json --out DIR [--replicates A:B] [--memory-limit SIZE]
bngsim version [--json]
```

The script starts Julia in the environment of the checkout (or
`BNGSIM_PROJECT`) with `BNGSIM_THREADS` threads (default: all) and calls
[`cli`](@ref). From Julia, `BnGsim.cli(["run", "s.json", "--out", "out"])`
does the same and returns the exit code.

## Output

Every command writes one JSON object to stdout and diagnostics to
stderr, so a caller reads the result from stdout and the reason for a
failure from both.

- `validate` prints `{"status": "valid", "scenario_sha256", ...}`, or
  `{"status": "invalid", "errors": [...]}` with every problem of the
  scenario and its path.
- `estimate` prints the expected peak memory, a conservative bound, the
  disk written and a rough run time ([Resource estimate](estimate.md)).
- `run` writes three files to `DIR`, which must not exist or must be
  empty (results are never overwritten):
  - `scenario.json`, the canonical scenario ([`scenario_json`](@ref));
  - `summary.csv`, the per-step summaries of every replicate
    ([`run_scenario`](@ref)), numbers in Julia's shortest round-trip form;
  - `provenance.json`: the canonical scenario and its SHA-256, the seed and
    replicates, the BnGsim, CLI and schema versions, the Julia and package
    versions, the SHA-256 of the active manifest, the operating system,
    CPU, host and threads, the start and end times, and the SHA-256 of the
    other two files.

  It then prints `{"status": "completed", "outputs": [...], ...}`. With
  `--memory-limit SIZE` (bytes, or a number with `K`, `M`, `G` or `T`),
  a scenario whose memory bound exceeds `SIZE` is refused before it
  runs.
- `version --json` prints the versions and the manifest hash.

`scenario.json` and `summary.csv` depend only on the scenario, the
replicates and the package versions. They are byte-identical whatever
the number of threads: for `examples/minimal.json`, checked at 1 and 8
threads. Everything about the run itself is in `provenance.json`.

## Exit codes

| Code | Meaning |
| :--- | :--- |
| 0 | success |
| 2 | usage: unknown command, missing or invalid arguments |
| 3 | the scenario is missing, unreadable or does not validate |
| 4 | refused before running: the output directory is not empty, or the replicates are outside the scenario's |
| 5 | refused before running: the memory bound exceeds `--memory-limit` |
| 6 | the simulation stopped with an explicit error, e.g. fewer segregating loci than the panels ask for |
| 7 | reserved: cancelled at a checkpoint |
| 70 | internal error (a bug) |

Code 7 belongs to `resume`, which comes later ([`ExitCode`](@ref)). A run
that fails or is refused writes no output directory.

## Start-up time

The scenario and command-line paths are compiled into the package image
when BnGsim is precompiled (a workload of two small scenarios), which
takes about half a minute once per installation. After that,
`bngsim run examples/minimal.json` takes about 3 s from the shell,
start-up included, instead of about 35 s.
