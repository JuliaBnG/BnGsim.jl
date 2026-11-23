# BnGsim.jl release notes

## v0.29.0 (2026-11-23)

### New

- Founder histories with demography and selection (Phase 2), scenario
  keys of `founders`: `popsize` (`{"bottleneck": {size, start, duration,
  recovery}}` or `{"growth": {rate, start, limit}}`, the census of the
  `nt` historical generations), `sexratio` (a fixed proportion of males),
  and `selection`, a quantitative trait under stabilizing or directional
  selection in the history (`nsites`, effects `da`, mutation rate `mu`,
  `sigma_e`, `fitness`; FisherWright's `TraitFitness`). With
  `loci.qtl.from = "history"` the loci of that trait segregating among
  the founders are the QTL of the (one) breeding trait, with their
  historical effects times one factor (new `scale_qtl`, which `eqtl` now
  uses); the chip comes from the other loci. The trait's alleles are
  merged into the founder haplotypes (`founder_haplotypes(res, nid,
  loci)`); a neutral mutation at a trait locus is dropped.
- Checks: the census follows the trajectory and the number of male
  parents is fixed; the merged haplotypes carry exactly the alleles of
  `selected_haplotypes` at trait loci and the neutral ones elsewhere;
  the QTL effects are the historical ones times one factor. Over 40
  independent histories each (300 diploids, 1500 generations, 2000
  trait loci): with near-neutral fitness cor(|a|, MAF) = 0.000 ± 0.011;
  stabilizing selection (ω² = 1) gives −0.110 ± 0.011 and fewer
  segregating loci (232 against 304); directional selection (β = 0.3)
  gives +0.108 ± 0.010, the favourable allele at 0.516 ± 0.006 and lower
  neutral diversity (MAF 0.037 against 0.093).

### Changed

- The canonical JSON of every scenario now lists `founders.popsize`,
  `founders.sexratio`, `founders.selection` and `loci.qtl.from`, so its
  SHA-256 changes; results with the defaults do not.

## v0.28.0 (2026-11-23)

### New

- Missing records (Phase 2). `founder_population(...; missing,
  unrecorded_herds)` and the scenario keys `traits[].missing` and
  `traits[].unrecorded_herds` (both `0.0`, in `[0, 1)`): a record is
  missing at random with probability `missing` (a keyed draw per animal,
  trait and time step of the record), and a herd never records a trait
  with probability `unrecorded_herds` (a keyed draw per herd and trait).
  A missing record is not made; with repeated records the opportunity
  still counts towards `records`, and `pop.nobs` counts the records made.
  The masks have their own keyed streams, so the other records keep
  their values. Evaluations use the records there are; `:phenotype`
  refuses missing records.
- Checks: missing fractions at the expected rate (single and repeated
  records); records not masked equal those of the unmasked population;
  herd masks constant within a herd.

### Changed

- The canonical JSON of every scenario now lists each trait's `missing`
  and `unrecorded_herds`, so its SHA-256 changes; with the defaults the
  simulation does not (golden loop checksum unchanged).

## v0.27.0 (2026-11-23)

### New

- Repeated records of several traits (e.g. methane emission with growth
  traits). `founder_population(...; records, repeatability, Rpe)` take one
  value per trait, and the scenario key `pe_correlation` (identity) gives
  the correlation of the permanent environments among the traits with
  repeated records; traits with single and repeated records can be mixed.
  The temporary residuals of the records an animal gets at one step are
  drawn jointly, with the residual correlation; records at different
  steps are independent. PBLUP, GBLUP and ssGBLUP fit the multi-trait
  repeatability model: one row per record unit (animal and step), the
  animal effect G₀ ⊗ K, the permanent environment P₀ ⊗ I over the
  repeated traits, R₀ among the records of a unit (missing records
  exact), contemporary groups per record; REML estimates G₀, P₀ and R₀.
  The summary gains `vpe_<trait>` for the repeated traits.
- Checks: the EBVs equal Henderson's equations solved densely (three
  traits: 3, 2 and 1 records, permanent environment on two, groups); the
  simulated permanent environments and same-step residuals have the
  specified correlations; REML under random selection (two traits, 3
  and 2 records, 100 replicates) recovers all nine components within
  1.2 SE, e.g. P₀ 1.410 ± 0.015, 0.662 ± 0.013, 0.386 ± 0.009 (true
  1.400, 0.667, 0.386) and R₀[1,2] 0.484 ± 0.005 (0.490).

### Changed

- `pop.PE0` is a matrix (traits × traits).
- The canonical JSON of every scenario now lists `pe_correlation`, so
  its SHA-256 changes; results do not (single-trait repeated runs
  checked byte for byte against 0.26.0, the golden loop unchanged).

## v0.26.0 (2026-11-23)

### New

- Repeated records (Phase 2), one trait. `founder_population(...;
  records, repeatability)` and the scenario keys `traits[].records` (`1`)
  and `traits[].repeatability` (`null`: `h2`): with `records = m > 1` an
  animal gets a record at every time step it is alive from the trait's
  `age` on, at most `m`, y = μ + a + pe + cg + e, with a permanent
  environment drawn at birth, σ²pe = (r − h²)Vp, and a temporary residual
  per record, σ²e = (1 − r)Vp (Vp = σ²a/h²). Records are made at the
  start of each step (`pop.records`); `pop.phenotype` is the mean of an
  animal's records. PBLUP, GBLUP and ssGBLUP fit the repeatability model
  (permanent environment as a random effect with I, when r > h²; groups
  per record when contemporary groups are fitted), and REML estimates
  σ²a, σ²pe and σ²e; the summary gains `vpe`.
- Checks: the EBVs equal Henderson's equations of the repeatability
  model solved densely; every animal has one record per step alive from
  its age up to the maximum; REML under random selection (h² = 0.2,
  r = 0.5, 4 records, 60 replicates) gives σ²a 1.020 ± 0.023 (true 1),
  σ²pe 1.517 ± 0.020 (1.5) and σ²e 2.491 ± 0.009 (2.5); under truncation
  σ²a 1.014 ± 0.032.
- Archive: `repeated.csv` (`id, trait, step, value`) for traits with
  repeated records; `traits.csv` gains `records` and `repeatability`.

### Changed

- With contemporary groups, `σ²cg = c2 · σ²a/h²` (unchanged for single
  records, where σ²a + σ²e = σ²a/h²).
- The canonical JSON of every scenario now lists each trait's `records`
  and `repeatability`, so its SHA-256 changes; with the defaults the
  simulation does not (golden loop checksum unchanged).
- New fields of `Population`; run state and checkpoints of earlier
  versions are refused (as for every version).

## v0.25.0 (2026-11-23)

### New

- Mating policies (Phase 3). `Scheme(; mating, pair_limit)`, scenario
  keys `scheme.mating` and `scheme.pair_limit`: `:random` (default, as
  before), `:min_coancestry` (the pairing of the selected sire and dam
  slots with least total coancestry), `:assortative` and
  `:disassortative` (by the index); `pair_limit` caps the offspring per
  sire–dam pair (`1`: no full sibs) with `:random` or `:min_coancestry`.
  Every parent keeps the number of offspring its selection gave it, so
  the policies combine with any selection, OCS included (contributions,
  then mate allocation).
- `transport(supply, demand, cost; cap)`: the exact minimum-cost
  transportation problem (successive shortest paths with potentials;
  deterministic ties), matching brute force on 300 random small cases.
- Checks over 60 replicates (100 offspring, 10 sires, 50 dams, PBLUP
  truncation, 8 steps): minimum-coancestry mating gives F 0.000, 0.004,
  0.021 at steps 2–4 against 0.013, 0.029, 0.046 at random, and
  0.087 ± 0.001 against 0.109 ± 0.001 at step 8, with the same
  asymptotic ΔF (0.017 per step) and gain (3.09 ± 0.08 against
  3.07 ± 0.09). With OCS (per cohort, ΔF = 0.01) it gives 0.046 against
  0.069 at step 8, for the same constrained kinship of the contributions.
- New dependency: DataStructures (a binary heap; already in the manifest
  through DataFrames).

### Changed

- The canonical JSON of every scenario now lists `mating` and
  `pair_limit`, so its SHA-256 changes; with the defaults the simulation
  does not.

## v0.24.0 (2026-11-23)

### New

- Contemporary groups (Phase 2). `founder_population(...; herds, c2)` and
  the scenario keys `founders.herds` (`1`) and `traits[].c2` (`0.0`):
  founders are placed in `herds` herds at random (keyed), offspring are
  born in their dam's herd, and each record of trait `j` gets the effect
  of its group, herd × time step of the record (birth plus the trait's
  `age`), drawn keyed per group from ``N(0, c_j^2 V_{p,j})`` with
  ``V_p = σ²_a + σ²_e`` of the founders; `h2` stays the heritability
  within groups. Effects of different traits are independent. When some
  `c2 > 0`, PBLUP, GBLUP and ssGBLUP fit the groups as fixed effects (a
  sparse incidence matrix) in place of the overall mean; `:phenotype`
  selection ranks on phenotypes that include the group effects.
- Checks: adding a constant to every record of one group leaves the EBVs
  unchanged (one and two traits), and changes them without groups;
  REML with groups fitted (c2 = 0.2, 20 herds, random selection, 80
  replicates) gives σ²e = 2.338 ± 0.017 (true 2.333; 3.000 if the group
  variance were absorbed) and σ²a = 0.965 ± 0.023, against 0.982 ± 0.024
  for the same design without groups (true 1).
- `pedigree` and the archive's `pedigree.csv` have a column `herd`;
  `traits.csv` a column `c2`.

### Changed

- The canonical JSON of every scenario now lists `founders.herds` and
  each trait's `c2`, so its SHA-256 in the provenance changes; with the
  defaults the simulation does not (the golden loop checksum is
  unchanged).
- Columns added to `pedigree` (a DataFrame) and to the archive files:
  readers by position must be updated.

## v0.23.0 (2026-11-23)

### New

- QTL effect distributions (Phase 2). Scenario families for `traits[].da`:
  `"laplace"` (`scale`), `"gamma"` with `"sign": "random"` (effects of
  either sign, as in Meuwissen, Hayes and Goddard 2001) and
  `"spike_slab"` (`pi`, `sd`, `spike_sd`: a fraction `pi` of the effects
  from `Normal(0, sd)`, the others exactly `0`, or from
  `Normal(0, spike_sd)`). The samplers `RandomSign(d)` and
  `SpikeSlab(π, sd; spike_sd)` are exported for the Julia API (`rand`,
  `mean`, `var` only). `eqtl` only rescales the draws of a single trait,
  so zeros and signs are kept (tested); with several traits only the
  first keeps the shape, as before.

### Changed

- The canonical form of a `"gamma"` effect distribution now includes
  `"sign": "positive"`, so the canonical JSON (and its SHA-256 in the
  provenance) of such scenarios changes; the simulation does not.
- `:random` selection draws its keyed uniforms for the eligible
  candidates only, not for every animal ever born; the values, keyed by
  id, and so the runs are unchanged (checked byte for byte against
  0.22.0).
- Removal of animals' records from the store is no longer planned (see
  the package roadmap).

## v0.22.0 (2026-11-23)

### New

- `Fgrm` in the summary: the mean genomic inbreeding of the cohort,
  ``G_ii - 1`` of VanRaden's first G with the founders' allele
  frequencies as base, over the loci segregating in the founders. As
  those frequencies are estimated from the founders, its expectation
  without selection is about ``F - (1 - F)/(2n_0 - 1)``: over 200
  replicates of `examples/minimal.json` with random selection
  (n₀ = 50, predicted offset −0.010), steps 2–4 gave −0.0106, −0.0102
  and −0.0103 (SE 0.003–0.004).
- `Scheme(; report_reliability = true)` (scenario key
  `report_reliability`): `reliability`, the mean ``1 - PEV/σ²_a`` of the
  candidates of each step for the index (σ² = w'G₀w) and, with several
  traits, `reliability_<name>` per trait, from the prediction error
  covariances of the PBLUP, GBLUP or ssGBLUP equations. Solved for 64
  candidates at a time, so memory stays bounded; checked against a dense
  inverse of the coefficient matrix. Off by default (`NaN`): it costs a
  solve per candidate and trait, and the resource estimate does not
  include it.

### Changed

- `summary.csv` and `aggregate.csv` have the new columns (`Fgrm`,
  `reliability`): readers by position must be updated; not breaking for
  readers by name.

## v0.21.0 (2026-11-23)

### New

- Report set (Phase 5). New per-cohort columns of the summary: `qtl_seg`,
  `qtl_fixed` and `qtl_lost`, the QTL with a non-zero effect on the
  genetic value (the index, with several traits) at which the favourable
  allele is segregating, fixed or lost; `intensity_sire` and
  `intensity_dam`, the realized selection intensity of each sex (the mean
  index of its parents, each counted once per offspring, less the mean
  index of the eligible candidates, in their standard deviations; `NaN`
  in the last row). Under random selection the mean is 0 (0.013 ± 0.027
  and −0.004 ± 0.008 over 200 steps of `examples/minimal.json`); under
  its truncation of 5 of about 25 males, 1.36.
- `aggregate` adds, derived within each replicate, `dFhet` (the rate of
  loss of heterozygosity, as `dF`) and the depletion of the genetic
  variance, `var_ratio` and `genic_ratio` (relative to the first step).

### Changed

- `summary.csv` and `aggregate.csv` have the new columns: readers that
  index columns by position must be updated (not breaking for readers by
  name).

## v0.20.0 (2026-11-22)

### New

- SIGINT (Ctrl-C) and SIGTERM to `bin/bngsim run` or `resume` stop the
  run at the next time step with the state saved, exit code 7: the
  launcher, which now runs Julia as a child process, creates
  `DIR/CANCEL`. Before the run has saved any state, or on a second
  signal, it stops Julia at once. The launcher handles the signals rather
  than Julia, whose delivery of interrupts to threaded tasks differs
  between versions. A resumed run is byte-identical to an uninterrupted
  one (checked after each signal).

## v0.19.0 (2026-11-20)

### New

- `aggregate(summary)` and `aggregate.csv`, written by every `bngsim run`:
  per time step, over the replicates, the mean, standard deviation and
  standard error of every summary column, and of the changes per step
  within each replicate (`dG` for the genetic means, `dF` for the rate of
  inbreeding). Rows are sorted by replicate first, so the result does not
  depend on their order.
- `bngsim merge DIR... --out DIR`: combine completed runs of the same
  scenario with disjoint replicates (e.g. `--replicates` on several
  hosts). The merged `summary.csv` and `aggregate.csv` are byte-identical
  to a single run of all the replicates, in whatever order the runs are
  given; `provenance.json` lists the runs merged.

## v0.18.0 (2026-11-20)

### New

- Adapter layer (roadmap §4.1): FisherWright, RelationshipMatrices and
  Breeding are imported only in `src/adapters/`, one module each, which
  re-exports the names BnGsim uses and absorbs differences between their
  versions (`census(res)` for the census record before and after
  FisherWright 0.10; `dense_haplotypes` passing an `rng` where
  FisherWright accepts one). A layering test fails if any other source
  file names these packages, and contract tests pin the conventions
  BnGsim relies on (haplotype layout and coordinates, `ainv = inv(nrm)`,
  `grm` and `hinv` against their formulas, `mme` against a dense solve
  of Henderson's equations, `TM1997`). Docs page "Upstream contract".

### Changed

- `founder_haplotypes` and `sample_founder_population` draw the
  placeholder nucleotides of their map from a copy of their `rng` (with
  FisherWright 0.10 or later) instead of the global generator, so the map
  is reproducible; the caller's stream does not advance, so the founders,
  QTL and everything else drawn from it are unchanged.

## v0.17.0 (2026-11-20)

### New

- Archives: `bngsim run ... --archive none|chip|qtl|all` writes, for each
  replicate, `pedigree.csv`, `records.csv` (phenotypes as recorded),
  `truth.csv` (true and estimated breeding values), `traits.csv`,
  `qtl.csv` (QTL and their effects) and PLINK 1 binary genotypes
  (`.bed/.bim/.fam`) at the loci of the panel, to `DIR/archive/r<N>/`,
  listed with their SHA-256 in `provenance.json`; a resumed run writes
  the same files. `write_archive` writes one from Julia. The genotypes
  were checked with PLINK 1.9 and plink2. Docs page "Archives".
- `infinitesimal(sire, dam, G0, R0, key)`: breeding values exactly
  `N(0, G0 ⊗ A)` and phenotypes on a pedigree, by Mendelian sampling, each
  animal from its own keyed stream (new `Domain.INFINITESIMAL`).

### Fixed

- The placeholder nucleotides (`ref`, `alt`) of the map of a scenario run
  came from Julia's global generator and differed between runs; they are
  now keyed from the scenario's seed. Nothing else depended on them.

## v0.16.0 (2026-11-20)

### New

- Stopping and resuming runs. While `bngsim run` is unfinished its state
  lives in `DIR/.bngsim-run/`. A file `DIR/CANCEL` stops it after the
  current time step with the state saved (exit code 7, `CANCELLED`,
  reserved until now); `--checkpoint-every K` also saves it every `K`
  steps, so a killed process loses at most `K` steps. `bngsim resume DIR`
  continues the run: its `scenario.json` and `summary.csv` are
  byte-identical to an uninterrupted run, and `provenance.json` lists
  every session. `run` refuses a directory holding an unfinished run.
  The state resumes only under the same BnGsim and Julia versions.

## v0.15.0 (2026-11-20)

### New

- `resource_estimate(scenario)` and `bngsim estimate SCENARIO.json`: the
  expected peak resident memory of a run, a conservative bound (1.3 × the
  expectation), the disk written and a rough run time, from the scenario
  alone. Memory is a component model (the process, the historical
  population, the packed genomes, the dense relationship matrices of
  GBLUP and ssGBLUP) whose sizes follow from neutral theory
  (`ScenarioSize`) and whose multiples were fitted on 16 scenarios and
  checked on 6 held out: every scenario's mean peak within 15 %, no run
  above the bound. Calibrated on one host (AMD Threadripper 3970X, Julia
  1.13.4, 8 threads); the measurements and the fit are in
  `validation/estimate/`. Docs page "Resource estimate".
- `bngsim run --memory-limit SIZE` refuses a scenario whose memory bound
  exceeds `SIZE` before it runs, with exit code 5 (`BUDGET_EXCEEDED`,
  reserved until now).

## v0.14.0 (2026-11-20)

### New

- Command line: `bin/bngsim validate | run | version` (CLI interface 1),
  calling `BnGsim.cli(args)`. One JSON object goes to stdout and
  diagnostics to stderr; exit codes are enumerated (`ExitCode`: 0 success,
  2 usage, 3 invalid scenario, 4 refused, 6 simulation failed, 70 internal
  error; 5 and 7 reserved). `run SCENARIO.json --out DIR` writes the
  canonical `scenario.json`, `summary.csv` (the summaries of every
  replicate, byte-identical at any number of threads) and
  `provenance.json` (the scenario and its SHA-256, versions, manifest
  hash, host, CPU, threads, times and the SHA-256 of the outputs), and
  never overwrites a non-empty directory. Docs page "Command line".
- A precompile workload compiles the scenario and command-line paths into
  the package image: a cold `bngsim run examples/minimal.json` takes about
  3 s instead of about 35 s; precompiling BnGsim takes about half a
  minute longer.
- New dependencies: SHA (standard library) and PrecompileTools.

## v0.13.0 (2026-11-20)

### New

- Scenarios: a simulation described in JSON (schema `"bngsim-scenario"`,
  version 1): species, historical population and founders, QTL and chip
  loci, traits, the breeding scheme, steps, seed and replicates; keys
  mirror the Julia API (`scheme` takes the keywords of `Scheme`).
  `read_scenario`, `parse_scenario` and `validate_scenario` check the
  whole document and report every problem with its path (unknown,
  repeated and mistyped keys, out-of-range values, inconsistent
  combinations); `scenario_json` and `write_scenario` write the
  canonical form, every key present, which reads back equal.
  `run_scenario` and `run_replicate` run it, every draw keyed from the
  scenario's seed, so that a replicate is the same alone or with others
  and on any number of threads. Example: `examples/minimal.json`; docs
  page "Scenarios".
- New dependency: JSON.jl (1).

### Breaking

- Requires FisherWright 0.8.2 or later (the per-chromosome genetic
  lengths of the species in the historical population); 0.8.1 is
  dropped.

## v0.12.5 (2026-11-19)

- Allows FisherWright 0.11 as well as 0.8.1, 0.9 and 0.10. FisherWright
  0.11 adds assortative mating, a `mate_correlation` record and result
  file format 2; BnGsim uses none of these, and its seeded output is the
  same with every allowed version. No breaking changes.

## v0.12.4 (2026-11-19)

- Fixed: `expand` threw `DimensionMismatch` in its first generation
  whenever the number of loci was not a multiple of 64. `Haplotype`
  pads its `BitMatrix` to `64·cld(nlc, 64)` rows, and `expand` copied
  the padded matrix into a gene drop that expects `nlc` rows. It now
  copies the loci only; output for a multiple of 64 is unchanged.
- `sample_id` and `sample_loci` (when `vm.frq` is empty) read the
  padding rows too; with the zero padding that `Haplotype` writes,
  frequency 0 kept them out, but non-zero padding gave a `BoundsError`.
  Both now read the loci only; output is unchanged.

## v0.12.3 (2026-11-19)

- Allows FisherWright 0.10 as well as 0.8.1 and 0.9. FisherWright 0.10
  moved the per-generation record of a result into `result.records`,
  made `to_haplotype(result)` return a `VariantMap` and removed
  `quickGT`/`quickHap`; BnGsim uses none of these, and its seeded output
  is the same with every allowed version. No breaking changes.

## v0.12.2 (2026-11-19)

- Allows FisherWright 0.9 as well as 0.8.1. FisherWright 0.9 changed
  only its trait records (matrices for several traits) and added
  sex-specific fitness; BnGsim runs neutral founder populations and reads
  none of those records, so its output is the same with either version.
  No breaking changes.

## v0.12.1 (2026-11-18)

- Requires Breeding 0.4, whose `TM1997` returns its contributions with a
  `feasible` flag and the coancestry reached, and no longer returns
  negative contributions when `K` cannot be met.
- `selection = :ocs` with `ocs = :fixed` reports the same columns as the
  dynamic constraints: `ocs_bound` (`K/2`), `ocs_kinship` (``c'Ac/2`` of
  the contributions, the same kinship scale), `ocs_L` (`NaN`) and
  `ocs_feasible`. When `K` cannot be met, the minimum-coancestry
  contributions are used and `ocs_feasible` is `false`, without a
  warning.

## v0.12.0 (2026-11-18)

### New

- `Scheme(; retain = :active)`: after each step, drop the genomes no
  later step can use, those of animals that can no longer be parents
  (dead, culled, or past the maximum age of their sex) and are outside
  the genotyped cohorts of `:ssgblup`. Founders are kept (base allele
  frequencies). Pedigree and records stay for every animal, so memory
  follows the animals still in use, not all animals born. A run with
  `:active` is identical to one with `:all` (the default). Not with
  `:gblup` or `:ssgblup` with every animal genotyped, which need every
  genome.
- `genomes(pop[, ids])`: the phased genomes of animals, in the order of
  `ids` (an `ArgumentError` for a dropped one); `retained(pop)`: the ids
  whose genomes are kept. Use these instead of reading `pop.blocks`.

### Breaking

- `Population` has a new field `slot` (each animal's position in its
  cohort's block, `0` once dropped), and `blocks` is a
  `Vector{Union{Nothing,Haplotype}}`. Checkpoints written by earlier
  versions do not load (they were tied to the BnGsim version already).
- `generation_summary` of a cohort whose genomes were dropped throws an
  `ArgumentError`; the summaries of a run are taken during it.

## v0.11.1 (2026-11-18)

### New

- `pedigree(pop[, ids])`: the animals as a `DataFrame` with the id
  explicit (`id`, `sire`, `dam`, `male`, `generation`, `exit`, and per
  trait `tbv_<name>`, `phenotype_<name>`, `ebv_<name>`; `missing` for an
  animal still alive, a trait not expressed, or no EBV). Inside a
  `Population` ids equal rows; a table that leaves it keeps the identity
  of every animal through a subset, a reordering or a concatenation, and
  `pedigree_view(pop, df.id)` rebuilds the relationships of a subset.

## v0.11.0 (2026-11-18)

### New

- `PedigreeView`: a pedigree restricted to a set of ids, with dense rows
  in topological order and `0` for a parent that is unknown or outside
  the view, the form RelationshipMatrices expects; `row_to_id` and
  `id_to_row` translate (an `ArgumentError` for an id not in the view).
  `pedigree_view(pop, ids; generations)` adds the ancestors of `ids`
  traced for `generations` generations (default: all, which keeps the
  relationships among `ids` exact).
- `Scheme(; window, ancestor_generations)`: `:pblup`, `:gblup` and
  `:ssgblup` on the records of the last `window` cohorts, on the
  pedigree of these animals and their ancestors traced for
  `ancestor_generations` generations. Animals outside the view get no
  EBV (`NaN`); `window` must exceed the maximum ages. The default
  `window = 0` is the whole population, as before (seeded runs are
  unchanged). Discarding earlier records and parents biases BLUP under
  selection; the window trades accuracy for the size of the equations.
- `evaluation_view(pop, scheme)`: the view an evaluation uses; the
  animal effects of the equations `evaluate!` returns are in its rows.

### Tests

- Property tests on random pedigrees: round trips `row_to_id ∘
  id_to_row` on sparse, unsorted and repeated ids; parents as rows;
  topological order; relationships among the ids equal to those of the
  full pedigree when every ancestor is traced.
- The window's PBLUP against an independent dense solve; a window over
  every cohort equal to the full evaluation; ROCS, ssGBLUP and REML with
  overlapping generations on a window.

## v0.10.0 (2026-11-17)

BnGStructs 0.4 (M1, second tranche): a species carries no population
size.

### Breaking

- `HistoricalSpec(name, sp, ne, nt, mr)`: the population size is an
  argument; it used to come from `sp.nid`.
- `sample_founder_population(species; ne, ...)`: `ne`, the size of the
  Fisher–Wright burn-in, is a required keyword.
- Requires BnGStructs 0.4 and FisherWright 0.8.1.

### Documentation

- New page "Relationships: A, G and IBD": A, realized IBD from founder
  segments (exact and at the markers), G, and four measures of
  inbreeding (pedigree, segments, ROH, G) on one simulated population,
  computed during the docs build.
- `validation/ibd-selection.jl`: realized IBD against the pedigree
  expectation over 20 replicates. Without selection mean(I_seg − A) is
  −0.004 (s.e. 0.004); under truncation on PBLUP EBV it is +0.033 (s.e.
  0.008), and F_seg − F is +0.015 (s.e. 0.004): selected relatives share
  the segments around favourable QTL alleles more often than Mendelian
  sampling implies.

## v0.9.1 (2026-11-17)

### New

- With `Scheme(; record_meiosis = true)` the per-generation summaries
  have a column `Fseg`: the cohort's mean realized inbreeding from
  founder segments, next to the pedigree `F` and `Fhet`.
  `generation_summary(pop, g; segs)` takes the segments.

## v0.9.0 (2026-11-16)

Realized IBD (M5, BnGsim side): founder segments from the loop's
crossover records.

### Breaking

- `Population` has a new field `meioses` (per cohort, the crossover
  records or `nothing`).
- Requires RelationshipMatrices 0.7.1 (segment-level IBD).

### New

- `Scheme(; record_meiosis = true)`: keep each cohort's `MeiosisRecord`,
  parents as population ids (about 8 bytes per crossover). The run is
  unchanged by recording.
- `meioses(pop)`: the records as `RelationshipMatrices.Meioses`.
- `founder_segments(pop, mm)`: the founder haplotype of every base of
  every animal (`RelationshipMatrices.HaplotypeSegments`); with it
  `irm(segs)` is the realized IBD relationship matrix and
  `segment_inbreeding(segs)` the realized inbreeding, without sampling
  loci.

### Validation

- Tests: at every locus every haplotype carries the allele of the
  founder haplotype its segment names (60 founders, 4 steps, overlapping
  generations); founders unrelated and parent–offspring relationships of
  founders' offspring exactly ½.
- `validation/ibd.jl`: realized relationships of paternal half sibs and
  full sibs from BnGsim's meiosis against their exact mean and variance
  under Haldane's map, on a genome with unequal cM/Mb (50/100/150 Mb at
  0.5/1.5/1.0 M), 400 replicates
  (`validation/results/ibd-2026-11-16.txt`): half sibs (n = 32 734) mean
  0.2495 (theory 0.25, s.e. 0.0005), variance 0.008812 (theory 0.008811,
  s.e. 0.000059); full sibs (n = 7248) mean 0.4989 (0.5, 0.0016), variance
  0.01740 (0.01762, 0.00027).

## v0.8.0 (2026-11-16)

Sixth part of M3: optimum contribution selection with a target rate of
inbreeding, for overlapping generations, and a comparison of three
constraints.

### Breaking

- `select_parents` returns a fifth element, `info` (empty unless a
  dynamic OCS constraint is used); code destructuring four values is
  unaffected.

### New

- `Scheme(; selection = :ocs, ocs, ΔF)`: OCS holding the rate of
  inbreeding at `ΔF` per generation. The bound on the constrained kinship
  rises from 1/(2n₀) by 1 − (1 − ΔF)^(1/L) per step.
  - `ocs = :cohort`: mean kinship of the newborn cohort, c'Θc; L estimated
    from the parental ages used so far.
  - `ocs = :age_year` (Meuwissen & Sonesson 1998; Wellmann 2019, optiSel):
    mean kinship of the population at the next step, w = W₀c + v, age
    classes weighted by their remaining expected contributions (estimated
    from realised usage), equal weights within a class; L as for
    `:cohort`.
  - `ocs = :age_generation`: the same with L the generation interval of
    the contributions being chosen, iterated to a fixed point (a
    per-generation constraint, Sonesson et al. 2000).
  - `ocs = :fixed` (default): `TM1997` with a fixed `K`, as before.
  - With discrete generations the three dynamic constraints coincide.
  - Summaries report `ocs_bound`, `ocs_kinship`, `ocs_L`, `ocs_feasible`.
- `ocs_contributions(u, Θ, b, K, male)`: OCS with a linear term in the
  constraint, solved exactly by homotopy in the weight of the gain from
  the best male and female down to the bound (or to the minimum of the
  kinship, flagged infeasible). The free set stays as small as the number
  of parents, so a step costs well under a second with hundreds of
  candidates.

### Validation

- Tests: KKT conditions of the solver at several bounds, with and without
  the linear term, and at the infeasible minimum; agreement with
  `TM1997` where its solution has no negative contributions; age weights
  summing to 1; the reported kinship equal to w'Θw recomputed; the
  per-generation fixed point; the three constraints identical with
  discrete generations, with the bound 1 − (1 − 1/2n₀)(1 − ΔF)^t.
- `validation/ocs-overlap.jl`: the three constraints and truncation,
  100 paired replicates per scenario, 20 steps, parents of ages 0–3,
  PBLUP, target ΔF = 0.01 per generation
  (`validation/results/ocs-overlap-2026-11-16.txt`, per-replicate CSVs):

  | scenario | method | ΔF/gen | gain/year | L | parents |
  | :--- | :--- | :--- | :--- | :--- | :--- |
  | records at birth | cohort | 0.0103 | 0.280 | 1.27 | 63 |
  | | age_year | 0.0104 | 0.271 | 1.28 | 59 |
  | | age_generation | 0.0100 | 0.282 | 1.29 | 60 |
  | | truncation | 0.0499 | 0.206 | 1.28 | 33 |
  | records at age 1 | cohort | 0.0107 | 0.219 | 2.36 | 61 |
  | | age_year | 0.0109 | 0.220 | 2.19 | 43 |
  | | age_generation | 0.0104 | 0.229 | 2.19 | 44 |
  | | truncation | 0.0585 | 0.205 | 1.66 | 33 |

  Standard errors are about 0.0001 for ΔF and 0.004–0.005 for gain.
  All three constraints hold ΔF within 10 % of the target; truncation
  without control is 5–6 times above it. `:age_generation` comes closest
  (0.0100 and 0.0104) and gives the most gain: +0.012 (s.e. 0.005) and
  +0.010 (0.004) per year over `:age_year` (paired). `:cohort` and
  `:age_year` give the same gain, but with records at age 1 `:cohort`
  uses older parents (L 2.36 against 2.19) and more of them. The
  differences between the constraints are small next to their common
  advantage over truncation, in gain per year as well as in ΔF over this
  horizon.

## v0.7.0 (2026-11-16)

Fifth part of M3: survival, culling and single-step GBLUP.

### Breaking

- `Population` has a new field `exit` (the first step at which an animal
  is dead or culled); the summaries have a new column `alive`.
- A record taken at an age the animal did not reach is missing (only
  relevant with `survival` or `cull_unselected`).

### New

- `Scheme(; survival)`: per-age probabilities of surviving to the next
  step, keyed draws (`Domain.SURVIVAL`), independent of genotype and
  records. Dead animals are not candidates and make no later records.
- `Scheme(; cull_unselected)`: eligible animals of age 1 or more that
  are not selected are culled; young candidates are kept for the next
  step.
- `evaluation = :ssgblup` with `Scheme(; genotyped_cohorts)`: single-step
  GBLUP (`RelationshipMatrices.hinv`) with the last cohorts genotyped,
  for one or several traits and with REML; with every animal genotyped
  it equals `:gblup`.
- `alive(pop, t)`, `genotyped(pop, scheme, t)`.

### Validation

- Tests: survival frequencies by age class against the probabilities
  (4 s.e.); dead and culled animals never parents; records missing for
  animals that did not reach the recording age; culling keeps exactly the
  selected adults; checkpoint–resume repeats the deaths; `:ssgblup` with
  everyone genotyped equal to `:gblup`; single-step with genotyped young
  animals, REML and two traits.
- `validation/accuracy.jl` (30 replicates, Fisher–Wright founders) now
  includes `:ssgblup` with the last 2 cohorts genotyped: accuracy 0.703
  (s.e. 0.006), above PBLUP by 0.166 (s.e. 0.009) and below GBLUP with
  every animal genotyped by 0.028 (s.e. 0.004). The other three methods
  repeat the v0.4.0 values exactly (`accuracy-2026-11-16-ssgblup.txt`).
- `validation/reml.jl ... cull` (300 replicates): records at age 1,
  parents up to age 2, survival 0.9, unselected adults culled; σ̂²a =
  0.983 (s.e. 0.007), σ̂²e = 2.332 (0.007), both within 3 s.e. of 1 and
  2.333 (`reml-cull-2026-11-16.txt`).
- `validation/reml-overlap.jl`, a diagnostic: 2–4 replicates of the
  culling scenario had suggested σ̂²a ≈ 0.75–0.85. Over 50 replicates for
  each of five configurations (selection rule, record age, parental ages,
  300 or 2900 QTL), no estimate is off by more than about 2 s.e. (largest
  deficit: records at birth, parents up to age 2, 0.966, s.e. 0.016).
  With records at age 1 and parents up to age 2, σ̂²a is 0.980 (0.019)
  with 300 QTL and 1.003 (0.017) with 2900, consistent with a small
  finite-locus effect (`reml-overlap-2026-11-16.txt`).

## v0.6.0 (2026-11-16)

Fourth part of M3: several traits, selection on an index.

### Breaking

- `Population` holds `tbv`, `phenotype` and `ebv` as `n × k` matrices, one
  column per trait (single-trait code: `pop.tbv[:, 1]`). `trait`, `σe` and
  `vc` are replaced by `traits`, `Lr` (factor of the true residual
  covariance), `G0` and `R0` (the covariance components of the next
  evaluation); `qtl` is a vector of one `TraitQTL` per trait.
- With several traits, `step!`/`simulate!` report `mean`, `var`, `genic`
  and `accuracy` of the index or aggregate genotype, and per-trait columns
  (`mean_<name>`, `var_<name>`, `accuracy_<name>`, `va_<name>`,
  `ve_<name>`, `rg_<a>_<b>`). Single-trait output is unchanged.

### New

- `founder_population(hap, ::MultiTraitQTL, traits, key; Re)`: several
  traits with residual correlation `Re`; the evaluation starts from the
  founders' genetic covariance (exactly the `R` of a multi-trait `eqtl`).
  The first trait's phenotypes are those of a single-trait run.
- Multi-trait PBLUP and GBLUP (Breeding.jl): records stacked trait-major,
  an intercept per trait, residuals correlated within animals, missing
  records (sex- or age-limited) handled exactly. A trait without records
  yet is predicted through `G0`.
- Multi-trait AI-REML with `reml = true`, from the first step with every
  trait recorded on an animal of known parentage; refused when two traits
  are never recorded on the same animal. Steps are taken on the Cholesky
  factors (`param = :cholesky`): after one cohort the maximum is often on
  the boundary (singular `G0`), where natural-scale AI steps were rejected
  and the EM fall-back did not converge in 100 iterates.
- `Scheme(; weights)`: selection on the index `w'â` (truncation, OCS,
  ROCS); ROCS takes the prediction error covariance of the index over all
  trait components. `generation_summary(pop, g; weights)` reports the
  aggregate genotype `H = w'g`.
- `observed`, `ntraits`, `index_weights`; `recorded(pop, t; trait)`.

### Validation

- Tests: founder covariances; multi-trait PBLUP with female-only records
  of one trait against a dense solve; truncation on the index; the
  prediction error covariance of the index against `pev_candidates` and a
  dense inverse; a trait without records predicted through `G0`; REML
  refusal for traits never recorded together; checkpoint–resume with two
  traits and REML; single-trait golden run unchanged.
- `validation/mtreml.jl`: multi-trait REML under index selection, one
  trait recorded in females only, 200 replicates (amd, 2026-11-16,
  `validation/results/mtreml-2026-11-16.txt`). All six elements of G0 and
  R0 are within 3 s.e. of the founder values under truncation on the
  index and under random selection (largest |z| 2.40 and 1.82). Under
  selection σ²a of the trait recorded in both sexes is 1.6 % low
  (z = −1.85), the same direction and size as in the single-trait gate.
  The mean of the ratio r̂g = Ĝab/√(ĜaaĜbb) is 0.52–0.53 against 0.5: a
  ratio of estimates is not unbiased, and it is not part of the gate.

## v0.5.0 (2026-11-16)

Third part of M3: overlapping generations, sex- and age-limited records,
re-estimated variance components and checkpoints.

### Breaking

- A trait expressed in one sex (`Trait(...; sex = 0)` or `1`) now leaves
  the phenotypes of the other sex `NaN` in `founder_population` and
  `step!`; they used to be recorded in both sexes. The values drawn for
  the expressing sex are unchanged (each animal has its own keyed
  stream).
- A trait's `age` must be a whole number of time steps.
- `Population` has a new field `vc`; `step!` and `simulate!` return the
  extra columns `L`, `va`, `ve`.

### New

- `Scheme(; max_sire_age, max_dam_age)`: parents are chosen among all
  animals of age at most these (in time steps; `0`, the default, gives
  discrete generations), ranked on one criterion. `:ocs` and `:rocs`
  take candidates of all eligible ages, without age-class weighting of
  contributions. `step!` reports the generation interval `L`.
- Records follow the trait: `trait.sex` limits them to one sex,
  `trait.age` delays their use in evaluations; `recorded(pop, t)` lists
  those available at step `t`, `current_step(pop)` the step.
- `Scheme(; reml = true)`: σ²a and σ²e re-estimated by AI-REML
  (Breeding.jl) before every PBLUP or GBLUP evaluation, starting from the
  previous estimates, once a recorded animal has a known parent; a run
  stops if REML does not converge. The components used are reported as
  `va` and `ve`.
- `save_checkpoint`, `load_checkpoint`, `simulate!(...; checkpoint)`: a
  run resumed from a checkpoint is identical to an uninterrupted one. The
  file uses Julia's `Serialization` and is refused by another BnGsim or
  Julia version unless `strict = false`: it is for resuming runs, not for
  archiving.

### Validation

- Tests: selection by age limits against a brute-force ranking; the
  generation interval; PBLUP with female records only against a dense
  solve; records available from an age; REML in the loop against an
  independent call; checkpoint–resume equal to an uninterrupted run
  (with REML and overlapping generations); discrete generations and the
  golden run unchanged.
- `validation/reml.jl`: REML estimates under selection against the
  founder values, 300 replicates (amd, 2026-11-16,
  `validation/results/reml-2026-11-16.txt`). Under truncation on EBV,
  σ̂²a = 0.986 (s.e. 0.007) and σ̂²e = 2.321 (0.006), against 1 and 2.333;
  under random selection, 1.009 (0.009) and 2.325 (0.007). The gate
  (|bias| < 3 s.e.) passes for both. Under selection both estimates sit
  about 2 s.e. low (selection − random for σ²a: −0.023, s.e. 0.012),
  which suggests a small bias of 1–2 % not shown at this precision. A bias
  of that size is expected with 300 QTL: the unbiasedness of REML under
  selection holds for the infinitesimal model, and selection here also
  changes allele frequencies.

## v0.4.0 (2026-11-16)

The breeding loop (second part of M3): discrete generations of one
additive trait, evaluated through Breeding.jl.

### New

- `Population` and `founder_population`: animals in birth order with
  pedigree, sex, genetic value, phenotype and EBV, and one haplotype block
  per generation; founder sexes and phenotypes are keyed.
- `Scheme`: evaluation by true values (`:tbv`), own phenotype,
  `:pblup` or `:gblup` (mixed-model equations of Breeding.jl with the true
  variance components; `G` at chip loci with founder allele frequencies,
  blended with `A`); selection by truncation, `TM1997` (`:ocs`), `rocs`
  (`:rocs`, with the candidates' prediction error covariance) or at random.
- `step!`, `simulate!`, `evaluate!`, `select_parents`, `mating_plan`,
  `generation_summary`: one generation is evaluate → select → pair (keyed
  shuffle) → keyed meiosis → genetic values, sexes and phenotypes. Every
  ranking breaks ties by ascending id. Per-generation summaries: mean,
  variance and genic variance of genetic values, pedigree F, the loss of
  heterozygosity, accuracy, parents used.
- `keyed_int`, `keyed_shuffle!`: integer draws from raw 64-bit outputs.
  Julia's `shuffle!` changed between 1.10 and 1.11, which made seeded runs
  differ across Julia versions; with these, the golden runs are identical
  on 1.10 and 1.13 and at 1 and 4 threads.
- Breeding.jl 0.3 is a dependency.

### Validation

- Tests: PBLUP against a dense solve of the mixed-model equations;
  truncation ranking and tie-breaking; the OCS coancestry constraint;
  `rocs` with λ = κ = 0 equal to truncation; offspring genetic means equal
  to parental means up to Mendelian sampling (4 s.e.); drift under random
  selection; golden checksums.
- `validation/accuracy.jl`: the accuracy ordering phenotype < PBLUP <
  GBLUP on a stated scenario with Fisher–Wright founders, paired over
  replicates.

## v0.3.0 (not registered; released as part of v0.4.0)

First part of the breeding-loop work (M3): keyed random streams, a
faster meiosis kernel, and a trait model that keeps effect distributions.

### Breaking

- `eqtl` samples different effects for a given `rng`. Effects are now
  only multiplied by one factor, so a positive distribution such as
  `Gamma` stays positive; the base-population mean is set through the
  intercept of the returned `TraitQTL`. Before, every effect was shifted
  to centre the genetic values, and the multi-trait method ignored each
  trait's `da` and imposed `R` on the effects rather than on the genetic
  values.
- Multi-trait `eqtl` gives genetic values whose sample correlation matrix
  is `R` exactly in the base population; it used to match `R` only in
  expectation over effects.

### New

- `SimKey` and `keyed_rng`: every random stream is addressed by a
  coordinate (seed, replicate, generation, domain, stage, entity,
  subcounter), mapped by Philox4x64-10 to the state of a `Xoshiro`.
  Draws do not depend on thread count, schedule or call order. A golden
  test fixes the layout; Philox is checked against Random123.jl.
- `MeiosisMap` and `make_progeny!`: keyed meiosis on packed haplotypes,
  one offspring per row of a mating plan, with per-chromosome genetic
  lengths from the species and optional crossover records
  (`MeiosisRecord`). With 50 000 loci and 10 000 offspring on 8 threads
  it is about 28 times faster than `gene_drop!`
  (`bench/results/meiosis-2026-11-16-ryzen3900x.txt`).
- `trait_summary`: mean, variance, genic variance and correlations of the
  genetic values of a population.
- Statistical tests of meiosis against exact theory: Poisson crossover
  counts per chromosome, Haldane recombination fractions (½ between
  chromosomes), unbiased segregation, and the variance of a gamete's
  genetic value under linkage. Seeds and 4-standard-error tolerances are
  declared in `test/meiosis.jl`.

### Fixed

- `gene_drop!(::Haplotype, ...)` threw `DimensionMismatch` unless the
  number of loci was a multiple of 64.

### Deprecated

- `norm_mv!` (no longer used by `eqtl`).

## v0.2.0 (2026-11-15)

Released as 0.2.0 rather than 0.1.2 because it breaks code written for
0.1.0: it requires BnGStructs 0.3, refuses maps that are not on the
genome of the species, and changes the realizations of
`continue_population`.

Requires BnGStructs 0.3 and FisherWright 0.7.1 or 0.8 (0.7.1 is the first
FisherWright release that allows BnGStructs 0.3); works with
RelationshipMatrices 0.6.1 and 0.7.

- The roadmap link in the README points to the deployed page
  (`…/BnGsim/dev/roadmap/`).

- `gene_drop!(og, pg, pm, species, map)`, `crossovers(species, map)` and
  `expand` check that the `VariantMap` is on the genome of the species
  (`BnGStructs.check_coordinates`). A map with within-chromosome
  positions used to give wrong recombinants without an error.
- Meiosis takes the genetic length of each chromosome from the species
  (`BnGStructs.map_length`), so `Species(...; morgan = ...)` is honoured.
  With the default map the crossovers are unchanged.
- `read_packed_haplotypes` reads through `BnGStructs.read_haplotypes`,
  which now owns the `BNGH1` format; its result is unchanged.
- Record columns are named through `BnGStructs.trait_columns`; the names
  (`tbv_<trait>`, `ft_<trait>`) are unchanged.

## v0.1.1 (not registered; released as part of v0.2.0)

- Works with FisherWright 0.6 and 0.7.
- `historical_population` takes `rng` as an alternative to `seed` and no
  longer seeds the global RNG. A given `seed` gives the same population
  as in v0.1.0.
- `continue_population` runs through `fisher_wright(; init, popsize)`
  instead of a serial loop on FisherWright internals. It is threaded,
  its result does not depend on the thread count, mutations that fix go
  to `substitutions`, and the per-generation records are extended.
  **Realizations for a given `rng` differ from v0.1.0.** `mr` is now
  subject to FisherWright's bounds (0.01 < `mr` < 20).
- `vanraden_g` is computed by `RelationshipMatrices.grm`; values agree
  with v0.1.0 to rounding.
- `gene_drop!(og, pg, pm, lms, tbp)` is the matrix method of gene drop.
  `gene_drop(pg, og, …)`, which mutated `og` without saying so, is
  deprecated, as are the camel-case aliases `sumMap`, `sampleID` and
  `sampleLoci`.
- The tests run Aqua.jl.

## v0.1.0

First registered release. It supports Julia 1.10 (LTS) through 1.13 and
requires BnGStructs 0.2.7, FisherWright 0.6 and RelationshipMatrices 0.6.1.

- Forward breeding simulation moved here from Breeding.jl (≤ 0.1.8) and was
  rewritten on the BnGStructs types: `sum_map`, `crossovers`, `gamete!`,
  `gene_drop`/`gene_drop!`, `random_mate`, `rnd_ped`, `expand`,
  `sample_id`, `sample_loci`, `eqtl`, `norm_mv!`, `tbv!`, `phenotype!` and
  `dummy_traits`. The camel-case aliases `sumMap`, `sampleID` and
  `sampleLoci` are kept. The old `BitMatrix`-based files are removed.
- `phenotype!` adds correlated residuals (`R_e`), sex-limited expression
  and threshold (`tTrait`) traits.
- Founders: `HistoricalSpec`, `historical_population`,
  `continue_population`, `founder_haplotypes`, `read_packed_haplotypes`,
  `write_xy`, and `FounderPopulation`/`sample_founder_population`
  (`VariantMap` and `LocusSet` based).
- Founder diagnostics: `founder_summary`, `ld_decay`, `vanraden_g`.
- `provenance` writes a reproducibility header with host, threads, seed,
  and the version and git commit of each package. Breeding.jl is not a
  dependency and is listed only when loaded.
- Docs: a roadmap page replaces the earlier working plans.
