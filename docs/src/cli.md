# Command line

`bin/bngsim` runs a [scenario](scenarios.md) from the shell. It is the
interface the planned web service will call, so its arguments, output and
exit codes are versioned (CLI interface 1) separately from the package.

```text
bngsim validate SCENARIO.json
bngsim estimate SCENARIO.json
bngsim run SCENARIO.json --out DIR [--replicates A:B] [--memory-limit SIZE] [--checkpoint-every K]
           [--archive none|chip|qtl|all]
bngsim resume DIR
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
- With `--archive PANEL`, `run` also writes an [archive](archive.md) of each
  replicate (pedigree, records, true and estimated breeding values, QTL,
  and PLINK genotypes at the loci of `PANEL`) to `DIR/archive/r<N>/`.
- `resume DIR` continues a run that was stopped or killed (below).
- `version --json` prints the versions and the manifest hash.

`scenario.json` and `summary.csv` depend only on the scenario, the
replicates and the package versions. They are byte-identical whatever
the number of threads: for `examples/minimal.json`, checked at 1 and 8
threads. Everything about the run itself is in `provenance.json`.

## Stopping and resuming

While a run is unfinished its state lives in `DIR/.bngsim-run/`: the
scenario, the shared historical population, the summaries of the
replicates done, and a checkpoint of the one in progress. A file
`DIR/CANCEL` (created by a user or the service) stops the run after the
current time step: the state is saved and `bngsim` exits with code 7.
`bngsim resume DIR` removes `DIR/CANCEL` and continues. With
`--checkpoint-every K` the state is also saved every `K` steps, so that a
process that is killed, or a host that goes down, loses at most `K` steps.

A resumed run gives the same `scenario.json` and `summary.csv`, byte for
byte, as one that was never interrupted, since every random draw is keyed
by its coordinate (replicate, step, animal) and nothing else carries over
between steps. `provenance.json` lists each session (start, end, host,
CPU, threads, and whether it was stopped or completed). `run` refuses a
directory that holds an unfinished run. The state uses Julia's
serialization, so it resumes only under the same BnGsim and Julia
versions; it is for continuing a run, not for archiving.

## Exit codes

| Code | Meaning |
| :--- | :--- |
| 0 | success |
| 2 | usage: unknown command, missing or invalid arguments |
| 3 | the scenario is missing, unreadable or does not validate |
| 4 | refused: the output directory is not empty or holds an unfinished run, the replicates are outside the scenario's, or there is nothing to resume (or its state is from another version) |
| 5 | refused before running: the memory bound exceeds `--memory-limit` |
| 6 | the simulation stopped with an explicit error, e.g. fewer segregating loci than the panels ask for |
| 7 | stopped at a step boundary by `DIR/CANCEL`; `bngsim resume DIR` continues |
| 70 | internal error (a bug) |

A run that fails before it saved any state, or is refused, writes no
output directory ([`ExitCode`](@ref)).

## Start-up time

The scenario and command-line paths are compiled into the package image
when BnGsim is precompiled (a workload of two small scenarios), which
takes about half a minute once per installation. After that,
`bngsim run examples/minimal.json` takes about 3 s from the shell,
start-up included, instead of about 35 s.
