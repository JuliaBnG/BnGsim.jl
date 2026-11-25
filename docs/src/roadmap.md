# Roadmap

This page is the single development plan for BnGsim.jl. It supersedes the
working drafts `plan.md`, `plan-gemini.md`, `plan-claude.md` and
`plan-unified.md` (August 2026). It also takes in the November 2026
reorganization of the JuliaBnG ecosystem, which moved forward simulation
out of Breeding.jl into BnGsim and left genetic evaluation and
optimum-contribution selection in Breeding.jl.

Upstream names cited below were checked against BnGStructs 0.3.0,
FisherWright 0.8.0, RelationshipMatrices 0.7.0 and Breeding 0.3.0. If the
installed packages disagree with this page, the installed packages win and
this page is amended.

## 1. Purpose and package boundaries

BnGsim.jl simulates breeding programs. It models founder populations,
meiosis, traits and records, lifecycle events, selection and mating, and
reproducible execution. It does **not** solve mixed-model equations, build
relationship matrices, or run forward sequence-level burn-in. It calls the
sibling packages for those.

```text
                  ┌──────────────────────────────────────┐
                  │             BnGStructs               │
                  │ (Species, Haplotype, VariantMap,     │
                  │  LocusSet, Trait, TraitQTL, tbv)     │
                  └──────┬───────────┬────────────┬──────┘
                         │           │            │
         ┌───────────────┘           │            └───────────────┐
         ▼                           ▼                            ▼
┌─────────────────┐      ┌────────────────────────┐     ┌──────────────────┐
│  FisherWright   │      │  RelationshipMatrices  │     │   Breeding.jl    │
│ (historical WF  │      │ (A, A⁻¹, G, H⁻¹, IRM,  │────►│ (MME, solvers,   │
│  burn-in, muts) │      │  DRM, APY, kinship)    │     │  REML/Gibbs,     │
└────────┬────────┘      └───────────┬────────────┘     │  GBLUP/Bayes,    │
         │                           │                  │  OCS/ROCS)       │
         └───────────────┬───────────┴──────────────────┴─────────┬────────┘
                         ▼                                        │
                  ┌───────────────────────────────────────────────▼──┐
                  │                     BnGsim                       │
                  │ (founder sampling, meiosis/gene-drop, QTL and    │
                  │  phenotype simulation, selection/mating loops,   │
                  │  historical/founder summaries, provenance)       │
                  └──────────────────────────────────────────────────┘
```

| Package | Responsibility as seen from BnGsim |
| :--- | :--- |
| **FisherWright.jl** | The only founder-population engine. Diploid, dioecious Wright–Fisher with infinite-sites mutation, recombination maps, demography (`bottleneck`, `growth`) and optional fitness or trait selection. BnGsim calls `fisher_wright` and `to_haplotype`; it does not use `quickGT` or `quickHap`. |
| **BnGStructs.jl** | The standard data types: packed phased `Haplotype`, `Genotype`, `VariantMap`, `LocusSet`/`SNPSet`, `Species`, trait definitions (`aTrait`, `tTrait`), `TraitQTL`/`MultiTraitQTL` and `tbv`. |
| **RelationshipMatrices.jl** | Pedigree validation, `A`, `A⁻¹`, `G`, `H⁻¹`, kinship and realized IBD. |
| **Breeding.jl** | (Phase 4 dependency; see §4.6.) Mixed-model equations and solvers, variance components (REML, Gibbs), GBLUP/Bayesian regression, PEV and reliability, optimum-contribution selection (`TM1997`, `OC_fixsex2`), and robust truncation and contribution selection (`rocs`: TS, RTS, OCS, ROCS; Fogg et al. 2026). |
| **BnGsim.jl** | Scenario definitions, founder sampling, meiosis, QTL and phenotype simulation, lifecycle, choosing which selection and mating policy runs when, results, checkpoints and the CLI. |
| **bngsim-service** | A separate application repository for the web UI/API, the PostgreSQL metadata store, the job protocol and the private worker. It runs the BnGsim CLI as a subprocess; it does not `import BnGsim`. |

AlphaSimR and MoBPS show how to build flexible breeding programs from
scenarios, and QMSim and SLiM serve as research references. None of them is
a dependency, and BnGsim does not aim for API compatibility with them.
Sequence-level forward simulation belongs in FisherWright (or SLiM), not
here.

## 2. Current status (v0.38.0)

| Area | State |
| :--- | :--- |
| Founders | Wright–Fisher histories through FisherWright (demography, selection, shared and stored histories, `write_histories`); `founder_haplotypes` samples founders from a `FisherWrightResult`, a packed BNGH1 file, or xyBnG `BitMatrix` output. The 0.2 wrappers `historical_population`/`HistoricalSpec`, `continue_population`, `read_packed_haplotypes`, `write_xy` and `sample_founder_population`/`FounderPopulation` are in `BnGsim.Legacy` (0.34.0). |
| Founder diagnostics | `founder_summary` (SFS, LD decay, chip/QTL MAF, VanRaden `G` statistics, implied `Mₑ`), `ld_decay`; `vanraden_g` in `BnGsim.Legacy`. |
| Keyed streams | `SimKey`, `keyed_rng` (Philox4x64-10 → Xoshiro), golden test (§4.2). |
| Meiosis | `MeiosisMap`, `make_progeny!` (keyed, word-wise, per-chromosome and sex-specific genetic lengths, crossover interference, `MeiosisRecord`); tested against exact theory. The older `crossovers(sp, vm)`, `gamete!`, `gene_drop!` with an explicit `rng`, `random_mate`, `rnd_ped`, `expand` are in `BnGsim.Legacy`. |
| Sampling | `sample_loci` (MAF-filtered chip and QTL `LocusSet`s); `sample_id` in `BnGsim.Legacy`. |
| Traits | `eqtl` (effects scaled only, mean in the intercept, each trait's `da` honoured, realised covariance exact in the base population), `scale_qtl`, `trait_summary`. The 0.1 helpers `tbv!`, `phenotype!`, `norm_mv!`, `dummy_traits` are in `BnGsim.Legacy`. |
| Provenance | `provenance.json` of every command-line run (I2), with `scenario_sha256`, `results_sha256` and `outputs_version` (0.35.0). The text header `provenance` is in `BnGsim.Legacy`. |
| Breeding loop | `Population` (ids equal rows; nothing is removed), `Scheme`, `step!`/`simulate!`: discrete or overlapping generations (parents up to a maximum age), one or several additive traits with sex- and age-limited records, selection on an index of multi-trait EBVs, optional (multi-trait) REML re-estimation of the covariance components, survival and culling, single-step GBLUP with a genotyped window, OCS with a target ΔF per generation for overlapping generations (newborn-cohort, age-weighted per year, age-weighted per generation), crossover records kept for realized IBD from founder segments (`founder_segments`, RelationshipMatrices 0.7.1), checkpoints that resume a run identically, evaluation through Breeding.jl (`:tbv`, `:phenotype`, `:pblup`, `:gblup`, `:ssgblup`), optionally on a window of recent cohorts and their traced ancestors through a `PedigreeView` (id ↔ row, topological, `0` for parents outside), selection by truncation, `TM1997`, `rocs` or at random, `(criterion, id)` ranking, keyed pairing and meiosis, per-generation summaries. |
| Scenarios | A versioned JSON scenario (`Scenario`, schema `"bngsim-scenario"` version 1, 0.13.0): species, historical population and founders, QTL and chip panels, traits, scheme, steps, seed, replicates. Reading and validation that report every problem with its path (unknown, repeated and mistyped keys, ranges, combinations), a canonical form that survives a JSON round trip, and `run_scenario`/`run_replicate`, keyed so that a replicate does not depend on the others or on threads. |
| Command line | `bin/bngsim validate | run | version` (0.14.0, CLI interface 1): one JSON object on stdout, diagnostics on stderr, enumerated exit codes; `run` writes the canonical scenario, `summary.csv` (byte-identical at 1 and 8 threads) and `provenance.json` (scenario and its SHA-256, versions, manifest hash, host, CPU, threads, times, output hashes; I2), and never overwrites. A precompile workload brings a cold `bngsim run examples/minimal.json` to about 3 s (Phase 0 gate: under 10 s). |
| Resource estimate | `resource_estimate` / `bngsim estimate` (0.15.0): expected peak memory from a component model (process, historical population, genomes, dense evaluation; sizes from neutral theory), fitted on 16 scenarios and checked on 6 held out (every mean within 15 %; `examples/minimal.json` within 10 %), a bound of 1.3 × the expectation that no run exceeded (0.37.0: with terms for several traits and REML, checked on 51 scenarios), disk, and a run time by phase (0.37.0: 51 scenarios, mostly within ±30 % per replicate; 0.37.1: 56, up to 20 000 genotyped) with `time_terms`; `run --memory-limit` refuses a scenario whose bound exceeds it (exit code 5). Memory calibrated on a Threadripper 3970X, time on a Ryzen 9 5950X; scripts in `validation/estimate/`. |
| Resume | `bngsim resume` (0.16.0): a file `DIR/CANCEL` stops a run after the current step with its state saved (exit code 7), `--checkpoint-every K` saves it periodically, and `resume DIR` continues; outputs byte-identical to an uninterrupted run, every session in the provenance. |
| Archives | `bngsim run --archive none\|chip\|qtl\|all` and `write_archive` (0.17.0): pedigree, records, true and estimated breeding values, traits, QTL, and PLINK 1 binary genotypes (checked with PLINK 1.9 and plink2), per replicate, in the provenance; `infinitesimal` draws exact `N(0, G0 ⊗ A)` breeding values and phenotypes on a pedigree (Tier 2 outputs of Phase 5, for the benchmark of roadmap §7a). |
| Adapter layer | `src/adapters/` (0.18.0): FisherWright, RelationshipMatrices and Breeding are imported only through one module each, which re-exports the names BnGsim uses and absorbs version differences (`census`, `dense_haplotypes`); a layering test forbids any other reference, and contract tests pin the conventions BnGsim relies on (docs page "Upstream contract"). This closes Phases 0–1. |
| Aggregation | `aggregate.csv` from every run and `bngsim merge` (0.19.0): per step, mean, SD and SE over replicates of every summary and of the per-step gain (`dG`) and rate of inbreeding (`dF`); runs of disjoint replicates (on several hosts) merge into the bytes of one run, whatever the order. Tier 1 summaries stay CSV (no Parquet dependency). |
| Reports | (0.21.0) per cohort, the QTL at which the favourable allele is segregating, fixed or lost (`qtl_seg`, `qtl_fixed`, `qtl_lost`), and the realized selection intensity of each sex (`intensity_sire`, `intensity_dam`: the parents' mean index weighted by offspring, less the eligible candidates' mean, in their SDs); in `aggregate.csv`, variance depletion (`var_ratio`, `genic_ratio`, relative to the first step) and `dFhet`. Already reported: ΔG, accuracy, generation interval, pedigree F and coancestry, `Fhet`, realized IBD (`Fseg`). 0.22.0 completes the list: `Fgrm` (VanRaden G with founder frequencies; E ≈ F − (1 − F)/(2n₀ − 1) without selection, checked over 200 replicates) and, with `report_reliability`, the candidates' mean reliability 1 − PEV/σ²a from the MME (solved in chunks of 64 candidates; checked against a dense inverse). |
| Signals | SIGINT and SIGTERM to `bin/bngsim` (0.20.0) create `DIR/CANCEL`, so the run stops at the next step boundary with its state saved (exit code 7); before any state exists, or on a second signal, Julia is stopped at once. Handled by the launcher, not by Julia, whose delivery of interrupts to threaded tasks differs between versions. Resumed runs are byte-identical. |
| Not planned | removal of animals' records from the store (decided 2026-11-23): a record is about 20 + 24k bytes (k traits), 44 MB per million animals with one trait, while genomes, the large part, are already dropped once unusable (`retain = :active`, 0.12.0); pedigree `F` and coancestry need every ancestor anyway. Ids stay rows of the store. Per-step work over the whole store is kept to cheap passes (0.23.0: `:random` selection draws only for eligible candidates, same values). |

The breeding loop is keyed by coordinate (§4.2): meiosis, sexes,
phenotypes, random selection and pairing each draw from their own
`SimKey`, and integers are sampled from raw 64-bit outputs
(`keyed_int`), so a run is identical at any thread count and on Julia
1.10 and 1.13. The older building blocks still take an explicit `rng`.

## 3. Design invariants

1. **I1 — Bitwise reproducibility.** With a pinned Julia version, manifest,
   CPU architecture and platform, a manifest and seed reproduce Tier 1
   results bit for bit, whatever the thread count, process count, or
   checkpoint/resume history.
   - **I1a — Deterministic decisions.** Every discrete decision derived from
     a floating-point quantity uses a total order with an explicit tie-break
     on ascending `id`. Selection, mating assignment and culling sort by
     `(criterion, id)`.
   - **I1b — Scoped cross-platform claim.** Across platforms, BnGsim claims
     identical discrete decisions only up to the first decision that uses a
     solver result. After that it claims statistical equivalence over
     replicates, never a tolerance on a single run. A 1-ulp change in an EBV
     can select a different animal, and from then on it is a different
     simulation.
2. **I2 — Immutable provenance.** Every output records the scenario JSON and
   its schema version, the BnGsim version, the manifest hash, the Julia
   version, OS/CPU, the worker identity, and the seed.
3. **I3 — Explicit failure.** No event silently drops animals, truncates
   candidate pools, repairs invalid values, or marks a partial result as
   successful.
4. **I4 — Identity stability.** Animal IDs are `Int32`, unique per
   replicate, and never reused or treated as row positions.
5. **I5 — One coordinate system.** Markers, QTL, chip loci and map positions
   all resolve through one ordered `VariantMap`.
6. **I6 — Pedigree monotonicity.** Parents exist before their offspring. Any
   dense row indexing that RelationshipMatrices or Breeding needs is built
   in an adapter and never exposed as an identity.
7. **I7 — Bounded footprint.** Memory, disk, CPU and wall time are estimated
   at admission and enforced at run time.
8. **I8 — Single execution path.** Tests, local CLI runs, web validation and
   worker runs all execute the same code path.
9. **I9 — Zero inbound attack surface.** The compute machine has no public
   listener and accepts no user code, shell commands, paths or worker flags.

## 4. Architectural decisions

### 4.1 Upstream contract and adapter layer

Before Phase 0 code depends on upstream internals, write
`docs/upstream-contract.md`. It records the exported API, the packing and
ploidy layout, the missing-parent sentinel (`0`), allele coding and index
base for each dependency. All upstream calls go through `src/adapters/`,
one thin module per dependency, each with a contract test. Serialized
manifests and result tables use types that BnGsim owns.

As of v0.2.0 BnGsim calls no unexported FisherWright function:
`continue_population` is `fisher_wright(...; init = res, popsize = ns)`.

### 4.2 Counter-based keyed PRNG (`SimKey`)

```julia
struct SimKey
    seed::UInt64        # scenario master seed
    replicate::UInt32
    generation::UInt32
    domain::UInt16      # FOUNDERS, MEIOSIS, PHENOTYPE, MATING, SURVIVAL, …
    stage::UInt16       # stable named event code, never a schedule index
    entity::UInt32      # 0 when population-scoped
    subcounter::UInt64
end
```

Each draw is a pure function of its coordinate. The key is mapped onto a
counter-based generator (Philox/Threefry) with a documented bit-packing:
`key = (seed, replicate)`, `counter = (generation, domain, stage, entity,
subcounter)`. `Base.hash` and `objectid` are prohibited on this path. A
golden test hard-codes the outputs for fixed keys.

Progeny draws are keyed by `(pair_index, within_pair_index)` from the mating
plan, **not** by the offspring `id`. IDs are then assigned in that same
order, so they do not depend on the thread schedule. Draws for animals that
already exist (phenotype, survival, genotyping) are keyed by `id`.

### 4.3 State storage and the ID/row boundary

Hot-loop state lives in a typed structure of arrays (`IndividualStore`:
`Vector{Int32}`, `Vector{Float64}`, `BitVector`). `DataFrame`s are
projections for output only. `PedigreeView` holds `id_to_row`/`row_to_id`
and the 1-based topological order with `0` for unknown parents. A property
test checks the round trip `row_to_id(id_to_row(x)) == x` on pruned, sparse
and reordered populations.

Haplotypes are stored in generation-sized blocks. Two retention policies
exist (0.12.0): *full* (`retain = :all`), and *active-only* (`retain =
:active`), which drops the genomes no later step can use (dead, culled or
too old to breed, and not genotyped) but keeps every animal's pedigree
and records; a cohort's block is compacted to the genomes kept.

### 4.4 Meiosis ownership

BnGsim decides *who and when*: the crossover model, map function,
interference, sex-specific maps and parent choice. The packed-genome layer
decides *how bits move*. In the current code, `gamete!` copies bit segments
on packed `Haplotype`s. The target interface is

```julia
offspring = make_progeny(parents, mating_plan, meiosis_model, key_base)
```

where `key_base` is a `SimKey` prefix, not a shared stream. The result is
then the same however the plan is split across threads.

### 4.5 Reproduction abstraction

`ReproductionModel` is declared in Phase 1 with outcrossing implemented.
Selfing, doubled haploids, clonal propagation and full-sib lines are
declared as variants now, so that plant schemes, where families or lines
rather than individuals are selected, do not force a rewrite.

### 4.6 Evaluation and contribution selection are delegated

BnGsim builds the evaluation inputs (records, pedigree projection, genotype
matrix) and calls Breeding.jl to solve. It does not have its own MME
solvers. Optimum-contribution and robust selection policies call
`Breeding.TM1997`, `Breeding.OC_fixsex2` or `Breeding.rocs`. `rocs` takes
the PEV/PEC matrix Ω from Breeding's `lhs_inverse`/`pev` and the
relationship matrix Σ (`A`, `G`, `H` or `I`), so one call covers TS, RTS,
OCS and ROCS.

*Julia support:* BnGsim supports Julia 1.10 (LTS) through the current
release. In v0.1 Breeding is **not** a dependency, and `provenance` reports
it only if the session has loaded it. Breeding lowered its floor to Julia
1.10 in v0.2.1, so Phase 4 can add `Breeding = "0.2.1"` without raising
BnGsim's Julia floor.

### 4.7 CLI-first boundary and layering

BnGsim.jl is the registered package. It holds the core, results,
checkpoints and the CLI. Heavy optional dependencies (Parquet, plotting)
go in `ext/` as package extensions. The web service is a separate
repository that calls the CLI.

The submodules form a layered DAG. Back edges between layers are not
allowed, and CI enforces this:

| Layer | May import from |
| :--- | :--- |
| L0 `adapters/`, `config/` | nothing in BnGsim (and not each other) |
| L1 `state/` | L0 |
| L2 `genetics/` | L0–L1 |
| L3 `breeding/`, `evaluate/` | L0–L2 |
| L4 `results/` | L0–L3 |
| L5 `cli/` | L0–L4 |

### 4.8 The CLI is versioned public API

```text
bngsim validate manifest.json   # schema + policy, no execution
bngsim estimate manifest.json   # predicted peak RAM, disk, runtime
bngsim run      manifest.json
bngsim resume   checkpoint.dir
bngsim version  --json          # package version, CLI interface version, manifest hash
```

Exit codes are enumerated and stable. They distinguish success, invalid
manifest, policy/resource rejection, budget exceeded, controlled
simulation failure, cancellation at a checkpoint, and internal error.
Structured output goes to stdout and diagnostics to stderr. The CLI
interface version is separate from the package version. A Tier 1 contract
test covers the arguments, exit codes and output shape.

## 5. Capacity model

For diploid packed genomes (1 bit per allele):

```text
haplotype bytes   = n_retained × n_loci / 4
dense G bytes     = 8 × n_genotyped²
sparse A⁻¹ nnz    ≈ 10–16 × n_pedigree
```

| Scenario | Loci | Retained | Haplotypes | Dense `G` | Note |
| :--- | ---: | ---: | ---: | ---: | :--- |
| Unit test | 2 000 | 500 | 0.25 MB | 2 MB | CI, < 10 s |
| Smoke / tutorial | 20 000 | 10 000 | 50 MB | 800 MB | laptop |
| Realistic cattle | 50 000 | 100 000 | 1.25 GB | 20 GB (50 000 genotyped) | dense `G` is the wall |
| Large study | 100 000 | 500 000 | 12.5 GB | infeasible | reference-set pruning, APY, iterative solver in Breeding |

The scaling bottleneck is `n_genotyped`, not haplotype storage.
`bngsim estimate` does this arithmetic at admission, and the runtime fails
cleanly before the OS kills the process.

## 6. Phases

Each phase ends with an exit gate. Phases 2 and 3 are partly done (see §2).

### Phase 0 — Executable vertical slice

1. CLI skeleton, `SimKey`, and a JSON provenance block that extends the
   current `provenance` header.
2. `docs/upstream-contract.md` and the first `src/adapters/` modules. (The
   FisherWright-internal calls were removed in v0.2.0.)
3. A minimal scenario: small genome, one additive trait, random mating,
   phenotypic truncation, 3 generations.
4. A benchmark harness on the smoke scenario, with a recorded baseline.

**Gate:** `bngsim run examples/minimal.json` finishes in under 10 s and
produces byte-identical Tier 1 output at 1 and 8 threads. The `SimKey`
golden test passes.

### Phase 1 — Scenario contracts and state

1. A versioned `SimulationConfig` JSON schema. It validates parameter
   bounds, locus ordering, QTL/chip overlap policy, trait dimensions and
   event schedules.
2. `IndividualStore`, `PedigreeView`, and the ID/row property tests.
3. The `ReproductionModel` and `MeiosisModel` extension points.
4. `bngsim validate`, `estimate`, `version --json`, and the exit-code table.

**Gate:** the config survives a JSON round trip unchanged. About 30 malformed
manifests are rejected with explicit errors. `estimate` is within ±25% of the
measured peak RAM on the smoke scenario. The CLI contract test passes.
(0.15.0: the ±25% applies to the expected peak, which is also checked on 21
other scenarios; the bound used for admission is never below a measured
run. The garbage collector alone moves the peak of one scenario by up to
25% between runs, so the gate is checked on the mean of two.)

### Phase 2 — Founders, genomic architecture, traits *(partly done)*

1. Founders come only from `fisher_wright`. Record the burn-in config, base
   frequencies, LD decay and genetic variance; `founder_summary` already
   computes most of these.
2. Chip `LocusSet`s and QTL architectures in one `VariantMap`.
   *Done:* `sample_loci`, `eqtl`.
3. QTL effect distributions: Gaussian, Gamma, Laplace, Gamma with random
   sign, spike-and-slab *(done, 0.23.0: scenario families and the
   samplers `RandomSign`, `SpikeSlab`)*.
4. Phenotypes: TBV plus fixed/contemporary-group effects plus correlated
   residuals. *Done:* correlated residuals; contemporary groups (0.24.0:
   herds, offspring in the dam's herd, herd × record-step effects with
   variance `c2` of the founders' `Vp`, fitted as fixed effects by the
   evaluations).
5. Threshold traits *(done)*, sex-limited traits *(done)*, age-limited
   records *(done: a trait's `age`)*, repeated records *(done, 0.26.0: one
   trait, records per step alive from its age, permanent environment from
   the repeatability, fitted as a random effect)*, missing-record masks
   *(done, 0.28.0: at random per record, and herds not recording a
   trait)*. Multi-trait repeated records *(done, 0.27.0: correlated
   permanent environments, residuals correlated among the records of one
   animal and step, single-record and repeated traits mixed)*.
6. Founder histories with selection or demography, using the FisherWright
   0.6 `FitnessLoci`/`TraitLoci`, `bottleneck` and `growth` *(done, 0.29.0:
   bottleneck, growth and sex ratio in the history; a quantitative trait
   under stabilizing or directional selection whose loci become the QTL,
   `loci.qtl.from = "history"`)*; several histories per run, stored once
   and reused (`bngsim history`, `founders.history_dir`), the
   between-history standard error in the aggregation, QTL drawn from and
   correlated with the historical trait (`nlc`, `rho`), and QTL
   architecture diagnostics per cohort *(done, 0.30.0)*.

**Gate:** realized `Var(TBV)` and the trait correlations recover their targets
within pre-declared sampling tolerances. LD decays monotonically. The
neutral SFS matches its expectation.

### Phase 3 — Meiosis, lifecycle, mating *(partly done)*

1. `MeiosisModel`: Poisson/Haldane *(done)*, Gamma-count interference,
   sex-specific maps *(done, 0.31.0: `MeiosisMap(...; male, female,
   interference)`, scenario `species.morgan_female` and
   `species.interference`; crossover counts and r(d) checked against the
   gamma model)*.
2. Demography: age classes, parity, survival, reproductive limits, culling,
   and discrete or overlapping generations. *(done: overlapping generations,
   maximum ages (0.5.0), survival and culling (0.7.0); 0.32.0: minimum breeding
   ages, offspring per parent and step (litter size, natural service),
   parity limits of dams)*.

3. Mating: random, assortative/disassortative, minimum-coancestry and a
   limit per sire–dam pair (factorial-type with limit 1) *(done, 0.25.0:
   `mating`, `pair_limit`; minimum-cost pairings solved exactly by
   `transport`)*. Fixed progeny per dam is the equal allocation of
   truncation selection. Joint mate selection (contributions and pairs
   optimised together) *(prototyped and not adopted, 2026-11-23)*: a
   search maximising the parents' index minus μ times the coancestry of
   the matings, under the contribution-kinship bound of `:ocs`, gave over
   40 replicates no more gain than OCS followed by minimum-coancestry
   mating and slightly higher inbreeding (8 steps, K = 0.09: G 1.354 ±
   0.058 and 1.339 ± 0.049 for μ = 1 and 10 against 1.370 ± 0.053; F
   0.0464 and 0.0447 against 0.0438). The two-stage plan already maximises
   gain under the bound and then minimises the coancestry of the matings.
   Worth revisiting only for constraints on the pairs themselves or a
   bound on the offspring's inbreeding.
4. Strict validation that parents are appended before offspring.

**Gate:** crossover counts match the map length in Morgans. Parent–offspring
realized IBD ≈ 0.5. Mendelian sampling variance
``\mathrm{Var}(\mathrm{MS}) \approx \tfrac{1}{2}(1-\bar F)\sigma_A^2``.
Pedigree `A` and realized IBD agree in expectation. Phase forms contiguous
parental segments. Allele dosage is conserved without selection.

### Phase 4 — Evaluation and selection (through Breeding.jl)

1. Evaluations, in this order: TBV oracle → phenotypic truncation → PBLUP
   (`A⁻¹`) → GBLUP (`G` from packed haplotypes: base frequencies, MAF
   filter, blending ``G_w = (1-w)G + wA_{22}``) → ssGBLUP (`H⁻¹`). Solving
   happens in Breeding.jl.
2. Multi-trait economic indices, sex-specific selection intensity,
   multi-stage selection.
3. Contribution policies: truncation, `TM1997` (target coancestry),
   `rocs` (TS/RTS/OCS/ROCS with Ω from PEV).
4. `(criterion, id)` ranking everywhere (I1a). Per-animal PEV and
   reliability ``R^2 = 1 - \mathrm{PEV}/\sigma_a^2``.

**Gate:**
- The TBV oracle gives ρ = 1 by construction.
- PBLUP agrees with an independent MME solve on a small data set.
- GBLUP and ssGBLUP agree when everyone is genotyped and `w = 0`.
- ρ(phenotype) < ρ(PBLUP) < ρ(GBLUP) ≤ ρ(ssGBLUP) on a named scenario
  built so that this should hold (dense panel, large reference). The
  replicate count is declared, and the margin must exceed the Monte-Carlo
  SE.
- `rocs` with κ = λ = 0 reproduces truncation selection.
- Selection does not change with thread count on a scenario with deliberate
  ties.

### Phase 5 — Persistence, checkpoints, replicates

1. Tier 1 compact summaries (JSON/Parquet) and Tier 2 archives (pedigree,
   EBVs, haplotypes).
2. Reports: ΔG, variance depletion, accuracy, PEV/reliability, intensity,
   generation interval, ΔF (pedigree and genomic), realized IBD, fixation and
   loss.
3. Replicate parallelism with aggregation that does not depend on order.
4. Generation-boundary checkpoints, `SIGINT`/`SIGTERM` handling, and
   `bngsim resume` as the only way to resume.

**Gate:** killing a run and then running `bngsim resume` gives output
byte-identical to an uninterrupted run. Aggregates do not depend on the order
in which replicates finish.

### Phase 6 — Web service (`bngsim-service`)

1. Scenario builder and submission UI that runs `validate`/`estimate`
   preflight.
2. A private worker that polls outbound over HTTPS/mTLS, takes atomic
   leases, runs under cgroup limits as a non-root user, and calls
   `bngsim run`.
3. Job state machine:
   `queued → claimed → running → completed | failed | cancelled`. A lease
   that expires goes back to `queued` or ends as `lost`. Terminal state
   comes only from the CLI exit code, and a job is never marked `completed`
   for a worker that stopped reporting.
4. Artifacts are uploaded idempotently under content-addressed keys, with
   retention classes, quota GC, schema migrations and user-safe error
   messages.

**Gate:** all acceptance scenarios (§8) pass with no inbound access to the
worker, against a released and pinned BnGsim. BnGsim.jl itself gains no
HTTP, database or auth dependency.

## 7. Verification tiers

| Tier | Trigger | Content |
| :--- | :--- | :--- |
| 1 — invariants and contracts | every commit, seconds, zero flakiness | schema round trip, ID/row bijection, pedigree monotonicity, thread invariance (1 vs 8), checkpoint/resume identity, adapter and CLI contracts, `SimKey` golden streams, layering, tie determinism |
| 2 — statistical | nightly, minutes, pre-registered tolerances | breeder's equation ``\Delta G = i\, r_{IH}\, \sigma_A``, Bulmer effect, neutral ``\Delta F \approx 1/(2N_e)``, Mendelian sampling variance, crossover rate vs map length, pedigree vs realized IBD |
| 3 — regression and cross-tool | releases, hours | golden scenario hashes, migration fixtures, qualitative trends vs AlphaSimR/QMSim, statistical (not run-level) cross-platform equivalence |

A change to a golden hash is never accepted as a plain test update. It
needs a reviewed diff that explains which decision changed and why.

## 8. Acceptance scenarios

1. On a pinned platform, two runs with the same seed give byte-identical
   Tier 1 output at 1, 2 and 8 threads, including a scenario with
   deliberate EBV ties.
2. Invalid parameters, coordinate mismatches and oversized resource
   requests are rejected before a job is created.
3. The worker claims a job through outbound HTTPS, runs a 20-generation
   ssGBLUP scenario, and uploads the results.
4. The UI shows the state transitions and generation progress, and offers
   the Tier 1 metrics and the manifest for download.
5. Resubmitting the downloaded manifest reproduces the run.
6. When a worker crashes, its lease expires and the job is requeued or
   marked `lost`, never `completed`.
7. A cancelled job stops at a generation boundary, and `bngsim resume`
   continues it to a byte-identical completion.
8. All Tier 2 tests pass nightly within their declared tolerances.

## 9. Deferred scope

| Feature | Re-entry trigger |
| :--- | :--- |
| Iterative (PCG) / APY evaluation | `n_genotyped` pushes a direct solve past the worker RAM budget (§5). The work goes in Breeding.jl. |
| Dominance, heterosis, epistasis | Additive validation is complete and a study needs non-additive variance. |
| Plant reproduction (selfing, DH, clones) | A concrete line/family workflow is scheduled. |
| Live UI (SSE/WebSocket) | Polling measurably degrades the user experience. |
| Slurm / HPC worker | Single-worker throughput is a measured bottleneck. |
| Zarr / multidimensional storage | Parquet Tier 2 output is a measured I/O bottleneck. |

Optimum-contribution selection is no longer deferred. Breeding.jl provides
the solvers (§4.6), so BnGsim only has to wire them in as a selection policy
in Phase 4.

## 10. Risk register

| Risk | Impact | Mitigation |
| :--- | :--- | :--- |
| API churn across four sibling packages | High | §4.1 contract document and adapters; tight `[compat]` bounds; remove calls to unexported internals. |
| Non-reproducibility under threads or resume | High | `SimKey` (§4.2); thread-invariance and resume-identity tests. |
| Floating-point tie divergence | High | I1a `(criterion, id)` ordering; I1b scoped claim; deliberate-ties gate. |
| Unspecified CLI contract leading to wrong job status | High | §4.8 enumerated exit codes; Tier 1 CLI test. |
| Meiosis keys that depend on creation order | High | Key progeny by `(pair_index, within_pair_index)` (§4.2). |
| ID/row conflation corrupting `A`/`G` | High | `PedigreeView`; round-trip property tests. |
| Dense `G` memory wall | Medium | §5 capacity model; `bngsim estimate` admission. |
| Flaky statistical tests getting muted | Medium | Tier 2 kept separate, with pre-registered tolerances and fixed seeds. |
| Duplicated evaluation code in BnGsim and Breeding | Medium | §4.6: BnGsim never solves MMEs. |
| Julia floor mismatch with Breeding | Resolved | Breeding 0.2.1 supports Julia 1.10 (§4.6); CI covers 1.10 and the current release. |
| Scope inflation | Medium | Phase 0 vertical slice first; §9 re-entry triggers. |
| Web dependencies leaking into the core | Medium | CLI subprocess boundary; Phase 6 gate. |
| Disk exhaustion on the worker | Medium | Tiered artifacts, retention classes, quotas, GC. |

## 11. Toward 1.0 (review of 2026-11-23, 0.33.0)

Phases 0–5 and the listed items of Phase 3 are done; Phase 6 (the web
service) is a separate repository. What a 1.0 would promise, and what
stands in the way:

**The contract for other programs** (the service, BnGBench): the scenario
schema, the command line (interface 1, exit codes), the output files and
the reproducibility invariants (I1–I3). These are stable in substance:
every key added since 0.13 is optional with a default that keeps earlier
scenarios' results, and keyed draws make results independent of threads
and of resume history. Three things are not yet stable:

1. **The canonical form grows.** Each release that adds a key changes the
   canonical JSON, and so the scenario SHA-256 in the provenance, of
   every scenario, although results do not change. A schema version
   still reads `1`.
2. **The summary columns change.** Columns were added in 0.21, 0.22 and
   0.30 (and their set depends on the traits and the selection method);
   readers by position broke each time. There is no version of the
   output files.
3. **The time estimate** was calibrated on single-trait scenarios with
   small genomes; for large genomes and multi-trait scenarios with
   repeated records it was about five times too low, and after the
   performance fixes of 0.36.0 two to three times. 0.37.0 models it by
   phase, calibrated on 51 scenarios (56 in 0.37.1, up to 20 000
   genotyped animals): mostly within about ±30 % per replicate on the
   calibration host, with one scenario 2.2 times low. 0.37.0 also
   fixed the memory estimate of several traits and REML, whose peaks had
   exceeded the bound (see "Resource estimate").

**The Julia API** has three layers:

- the breeding loop and its pieces (`Population`, `Scheme`, `step!`,
  `simulate!`, `evaluate!`, `select_parents`, `mating_plan`,
  `founder_population`, `MeiosisMap`, `make_progeny!`, the keys, the
  scenario functions, `aggregate`, `write_archive`, `write_histories`,
  `resource_estimate`): documented and tested, but `Population` is a
  mutable struct whose fields (about 35, one or more added in most
  releases) users read directly;
- the reports and relationship helpers (`generation_summary`,
  `founder_segments`, `irm`, `pedigree_view`, …): stable;
- the building blocks of 0.1–0.3 (`sum_map`, `crossovers`, `gamete!`,
  `gene_drop!`, `random_mate`, `rnd_ped`, `expand`, `sample_id`,
  `norm_mv!`, `tbv!`, `phenotype!`, `dummy_traits`, `FounderPopulation`,
  `sample_founder_population`, `HistoricalSpec`, `historical_population`,
  `continue_population`, `read_packed_haplotypes`, `write_xy`,
  `vanraden_g`, the text `provenance`): they draw from an explicit or the
  global generator, outside the keyed streams (invariant I1), and
  `vanraden_g` duplicates RelationshipMatrices' `grm`. Three camelCase
  aliases (`sumMap`, `sampleID`, `sampleLoci`) and `gene_drop(pg, …)` are
  deprecated since 0.2.

**Decisions before 1.0** (for the maintainer):

- a. The building blocks: remove them, keep them in a `BnGsim.Legacy`
  module outside the 1.0 promise, or keep them as they are. Recommended:
  `Legacy`, and remove the deprecated aliases.
- b. The schema: keep version 1 with the rule "keys may be added, with
  defaults that keep results; the canonical form, and its hash, may gain
  them in minor releases", or freeze the key set as version 2 at 1.0.
  Recommended: the rule, with a separate hash of the scenario's
  *results-relevant* content, so that adding a key with its default does
  not change it.
- c. `Population`: document its fields as API, or give accessors and
  declare the fields internal. Recommended: accessors for what users read
  (records, EBVs, true values, phenotypes, pedigree, genomes), fields
  internal.
- d. Output files: a version in the provenance (`outputs_version`) and a
  documented column table, columns only appended after 1.0.

**Decided (2026-11-23): all four as recommended.** Done in 0.34.0: a
(the building blocks in `BnGsim.Legacy`, the deprecated aliases and
`gene_drop(pg, …)` removed) and c (`tbv`, `ebv`, `phenotypes`, `records`,
`traits`, `trait_qtl`, `covariances`, `parity` beside `pedigree`,
`genomes` and the others; the fields of `Population` internal). Done in
0.35.0: b (`results_sha256`, which `bngsim merge` compares; the schema
policy in "Scenarios") and d (`outputs_version` 1 and the page "Output
files"; the columns of `summary.csv` pinned by tests). What remains open
for 1.0 was item 3 of the contract: the time estimate of large genomes and
multi-trait scenarios with repeated records. 0.36.0 removed two costs
that grew with the genome (the examples then ran two to three times their
estimate, not five); 0.37.0 models the run time by phase (the examples
1.04 and 1.26 times their estimate per replicate on the calibration
host); 0.37.1 measured the cubic cost of large dense evaluations, up to
20 000 genotyped animals, which 0.37.0 had not separated from the
quadratic ones. Nothing in the contract remains open. 0.38.0 narrows
the API to what 1.0 promises (decided 2026-11-25): `sample_loci`
requires `rng`; `philox4x64`, `keyed_int` and `transport` are
unexported; BnGStructs types are promised as far as BnGsim's compat
bound allows; the promise is stated on the API page.
