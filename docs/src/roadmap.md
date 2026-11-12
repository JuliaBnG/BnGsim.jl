# Roadmap

This page is the single development plan for BnGsim.jl. It supersedes the
working drafts `plan.md`, `plan-gemini.md`, `plan-claude.md` and
`plan-unified.md` (August 2026). It also takes in the November 2026
reorganization of the JuliaBnG ecosystem, which moved forward simulation
out of Breeding.jl into BnGsim and left genetic evaluation and
optimum-contribution selection in Breeding.jl.

Upstream names cited below were checked against BnGStructs 0.2.7,
FisherWright 0.6.0, RelationshipMatrices 0.6.1 and Breeding 0.2.0. If the
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

## 2. Current status (v0.1.0)

| Area | State |
| :--- | :--- |
| Founders | `historical_population`/`HistoricalSpec` wrap `fisher_wright`. `continue_population` extends a run. `founder_haplotypes` samples founders from a `FisherWrightResult`, a packed BNGH1 file, or xyBnG `BitMatrix` output. `sample_founder_population` returns a `FounderPopulation` built on `VariantMap` and `LocusSet`. |
| Founder diagnostics | `founder_summary` (SFS, LD decay, chip/QTL MAF, VanRaden `G` statistics, implied `Mₑ`), `ld_decay`, `vanraden_g`. |
| Meiosis | `crossovers` (Poisson, no interference, independent assortment), `gamete!` and `gene_drop`/`gene_drop!` on packed BnGStructs haplotypes. `random_mate`, `rnd_ped`, `expand`. |
| Sampling | `sample_id`, `sample_loci` (MAF-filtered chip and QTL `LocusSet`s). |
| Traits | `eqtl` (single trait, and multi-trait with genetic correlation matrix), `norm_mv!`, `tbv!`, `phenotype!` (correlated residuals, sex-limited expression, threshold traits from a liability scale), `dummy_traits`. |
| Provenance | `provenance` (Breeding is listed when loaded) writes a text header with the date, host, CPU, threads, seed and package versions/commits. This partly covers I2: there is no JSON manifest yet. |
| Not started | Scenario schema, `SimKey`, `IndividualStore`/`PedigreeView`, selection and evaluation loop, lifecycle, checkpoints, CLI, adapter layer. |

All randomness takes an explicit `rng`. Runs reproduce exactly for a given
seed only at a fixed thread count, because FisherWright and the threaded
kernels here use task-local streams. Section 3.2 removes that limit.

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

*Known debt:* `continue_population` calls the unexported FisherWright
internals `random_mate!`, `cobp!`, `recombine` and `merge_sorted`. In
FisherWright 0.6, `fisher_wright(...; init = res)` continues a run through
the public API, and `continue_population` should be rewritten to use it.

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

Haplotypes are stored in generation-sized blocks (`Vector{Haplotype}`).
Two retention policies are planned: *full*, and *active-only*, which drops
the genomes of culled, non-breeding animals but keeps their pedigree and
records.

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
2. `docs/upstream-contract.md` and the first `src/adapters/` modules. Remove
   the FisherWright-internal calls (§4.1).
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

### Phase 2 — Founders, genomic architecture, traits *(partly done)*

1. Founders come only from `fisher_wright`. Record the burn-in config, base
   frequencies, LD decay and genetic variance; `founder_summary` already
   computes most of these.
2. Chip `LocusSet`s and QTL architectures in one `VariantMap`.
   *Done:* `sample_loci`, `eqtl`.
3. QTL effect distributions: Gaussian *(done)*; Gamma, Laplace and
   spike-and-slab *(to do)*.
4. Phenotypes: TBV plus fixed/contemporary-group effects plus correlated
   residuals. *Done:* correlated residuals. *To do:* contemporary groups.
5. Threshold traits *(done)*, sex-limited traits *(done)*, age-limited and
   repeated records *(to do)*, missing-record masks.
6. Founder histories with selection or demography, using the FisherWright
   0.6 `FitnessLoci`/`TraitLoci`, `bottleneck` and `growth`.

**Gate:** realized `Var(TBV)` and the trait correlations recover their targets
within pre-declared sampling tolerances. LD decays monotonically. The
neutral SFS matches its expectation.

### Phase 3 — Meiosis, lifecycle, mating *(partly done)*

1. `MeiosisModel`: Poisson/Haldane *(done)*, Gamma-count interference,
   sex-specific maps.
2. Demography: age classes, parity, survival, reproductive limits, culling,
   and discrete or overlapping generations.
3. Mating: random *(done)*, assortative/disassortative, minimum-coancestry,
   factorial, fixed-progeny.
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
