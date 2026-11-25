# Resource estimate

[`resource_estimate`](@ref) and `bngsim estimate` predict, from the scenario
alone, the peak resident memory of a run, a conservative bound on it, the
disk the command line writes, and a rough run time. `bngsim run
--memory-limit SIZE` refuses a scenario whose bound exceeds `SIZE` before
it starts (exit code 5), so that a run fails cleanly instead of being
killed by the operating system (roadmap invariant I7).

## The model

Peak memory is the sum of four terms, each a fitted multiple of a size
computed from the scenario ([`ScenarioSize`](@ref)):

| Term | Size | Fitted bytes per byte |
| :--- | :--- | ---: |
| process | BnGsim loaded and compiled | 577 MB |
| historical population | one copy of its haplotypes, ``8N\,\bar m`` with ``\bar m = \mu L \min(t, 4N)`` derived alleles per haplotype | 5.42 |
| genomes | the packed genomes of all animals, ``2 \cdot n_{animals} \cdot S / 8``, with ``S = 4N\mu L \ln(1 + 2 n_{id} t / 4N)`` segregating loci among the founders | 4.15 |
| dense evaluation | one ``n_g \times n_g`` matrix, ``8 n_g^2``, for the genotyped animals | 13.3 (`:gblup`), 19.0 (`:ssgblup`), one trait; more for several traits and REML (below) |

``\bar m`` and ``S`` follow from neutral theory for a population that
starts without variation: ``\mu L`` new alleles per haplotype per
generation until fixation balances them near ``4N\mu L``, and the total
branch length of the founders' genealogy truncated at the start of the
run. On the calibration scenarios they matched the realized values within
4 %. The multiples count the copies the code holds at its peak (buffers
of the historical run; the relationship matrix, its blend with A, its
factor and inverse and the mixed-model equations of GBLUP), plus memory
not yet collected. Terms for the pedigree and for the dense cohort matrix
of the per-step summary were tried and fitted to zero.

The bound is the expected peak times 1.3: no calibration or held-out run
exceeded the expectation by more than 15.3 %.

The run time has its own model, by phase of the run (below, "Run time").

## Calibration of the memory model

The memory coefficients were fitted on one host: AMD Threadripper 3970X,
Linux x86_64, Julia 1.13.4, 8 threads, BnGsim 0.14. Sixteen calibration
scenarios cover the historical population (``N`` 500 to 3,000, ``t`` 500 to
4,000), PBLUP with 6,000 to 41,000 animals, and GBLUP and ssGBLUP with
2,000 to 8,000 genotyped animals; six held-out scenarios, not used in the
fit, vary the genome, the number of traits, the genotyped cohorts and the
historical size, and include `examples/minimal.json`. Each ran twice in a
fresh process. The garbage collector makes the peak of one scenario vary
between runs by up to 25 % (GBLUP with 8,000 animals: 4,512 and 5,761
MB), so the expectation is fitted to the mean of the two runs.

```@eval
using Markdown, BnGsim
f = joinpath(pkgdir(BnGsim), "validation", "estimate", "fit.txt")
Markdown.MD(Markdown.Code("text", read(f, String)))
```

The mean of every scenario lies within 15 % of its estimate (calibration
within 11.4 %, held out within 14.3 %), and `examples/minimal.json`, the
roadmap's smoke scenario, within 10 %. The test suite checks the roadmap's
gate on these measurements: every mean within 25 % of the estimate and no
run above the bound.

On another host, Julia version or thread count the coefficients change:
on a 12-core desktop (AMD Ryzen 9 3900X) at 4 threads the same scenarios
peaked 6 % to 28 % lower than on the calibration host at 8 (GBLUP with
6,000 genotyped animals: 2,536 MB against about 3,500 MB), so there the
estimate errs high, on the safe side. The scripts to measure and refit are in
`validation/estimate/` (`run.sh`, `measure.jl`, `fit.jl`).

## Several traits and REML (0.37)

The calibration of 0.15 had one trait and no REML. With ``k`` traits the
mixed-model equations of GBLUP and ssGBLUP are ``k n_g`` square, and REML
holds further matrices of that size, so the dense term is ``8 n_g^2`` times

```math
c + t\,(k^2 - 1) + r\,k^2\,[\text{REML}],
```

with ``c`` the one-trait multiple above, ``t = 9.7`` (`:gblup`) and
``15.1`` (`:ssgblup`), and ``r = 9.95``. Until 0.36 the term ignored
``k`` and REML, and the peaks of such runs exceeded the bound: up to 1.85
times with two traits and 3.06 times with two traits under REML (amd, 8
threads). ``t`` and ``r`` are fitted to the peaks of the time calibration
(below), two scenarios per evaluation, so exactly; REML with ssGBLUP uses
the GBLUP value of ``r``. Five held-out scenarios in
`validation/estimate/check/` (three traits; REML with one trait, and with
two traits under ssGBLUP; repeated records with GBLUP; three traits with
repeated records) had mean peaks 0.94 to 1.13 times their expectation.
Over all 51 scenarios of the time calibration the largest peak is 0.96
times its bound (the table at the end of the fit below). These peaks were
measured on the time calibration host, the one-trait multiples on the
memory calibration host; on the one-trait scenarios the two agree within
the 15 % of the memory calibration.

## Later terms (0.33)

Features added after the calibration contribute terms computed from their
data structures, checked against measured runs rather than fitted to many
(scenarios and runs in `validation/estimate/later/`):

- several shared histories (`founders.histories`): the history being
  simulated as before, and each further one held as a result at 2.4 bytes
  per byte of its haplotypes (fitted to the one scenario with four);
- a census trajectory (`founders.popsize`): the history term uses the
  largest census, since FisherWright sizes its buffers by it, and `N` in
  the allele and locus counts the harmonic mean of the census; the run
  time sums the work of every generation at its census;
- repeated records: 40 bytes per record;
- `report_reliability`: its work array, 64 columns per trait of the
  mixed-model equations, and one solve per candidate in the run time.

Measured on svr (8 threads, Julia 1.13.1, two runs each):

| scenario | measured (MiB) | expected (MiB) | bound (MiB) |
| :--- | :--- | ---: | ---: |
| four shared histories | 1339, 1453 | 1399 | 1818 |
| history growing to 6000 | 1663, 1534 | 1661 | 2159 |
| three traits, repeated records | 1550, 1586 | 1367 | 1777 |
| reliability, 3000 candidates | 1333, 1310 | 1217 | 1582 |

These runs first exceeded their bounds by up to 38 %: the genomic
inbreeding of the summaries (`Fgrm`, 0.22) formed a dense floating-point
matrix of loci × animals of the cohort every step (570 MiB for 25 000 loci
and 3000 animals). 0.33 computes it from per-locus sums of the dosages;
the peaks above are after that change.

## Run time

The run time (`seconds`, and by part `time_terms`) is the sum of the
phases of a run, each a fitted multiple of a size computed from the
scenario (`_time_features`), summed over the time steps ``g = 0, \dots,
T-1`` at which the loop evaluates and selects, when the pedigree holds
``n_{id} + g\,n_{off}`` animals:

| Part | Size, per replicate |
| :--- | :--- |
| history | alleles carried through meiosis in each historical run, about ``2N\mu L \sum_g \min(g, 4N)`` (per history simulated), and the same times ``\log_2`` of the haplotypes of one generation over 16 MB (zero below) |
| setup | founder haplotypes × segregating loci, ``2 n_{id} S`` |
| meiosis | offspring haplotypes × loci, ``2 T n_{off} S`` |
| summary | cohort and founder haplotypes × loci (allele counts), and animals × cohort (pedigree relationships) |
| evaluation | equations and equations × ``k`` (sparse PBLUP), ``(k\,n_g)^2`` and ``(k\,n_g)^3`` (dense GBLUP, ssGBLUP), equations and ``(k\,n_g)^3`` again under REML, repeated records; and the G or H inverse, ``n_g^2`` × chip and ``n_g^3``, separately for ssGBLUP |
| mating | animals × parents (minimum-coancestry mating) |
| reliability | candidates × equations (`report_reliability`) |
| rest | per step, and per animal and step |
| compile | the first replicate: a cost per feature used (BLUP, dense G, ssGBLUP, several traits, repeated records, REML, minimum-coancestry mating) |

with ``S`` the segregating loci, ``k`` the traits and ``n_g`` the
genotyped animals at the step (all evaluated animals with `:gblup`, the
last `genotyped_cohorts` cohorts with `:ssgblup`).

### Calibration of the time model (0.37)

The time model is calibrated on a different host from the memory model:
an AMD Ryzen 9 5950X (16 cores, 93 GB), Linux x86_64, Julia 1.13.4, 8
threads, BnGsim 0.37, otherwise idle. `validation/estimate/run_phases.sh`
runs every scenario of `scenarios/`, `later/`, `phases/` and `check/` and
the examples, 51 in all, in two processes; each process runs replicate 1
(cold, with compilation) and replicate 2 (warm), each with its own
history, timing every phase (`measure_phases.jl`, `phases.tsv`).
`fit_phases.jl` fits each phase to the mean warm time of its scenarios by
non-negative least squares on relative error (runs under 0.2 s counted at
0.2 s), and the compilation to the cold-minus-warm difference. The
examples, the held-out scenarios `v1`–`v6` of the memory calibration and
`check/` are held out.

On every scenario over 1 s but one, measured warm time per replicate was
0.92 to 1.22 times the model on the calibration scenarios, and 0.89 to
1.45 on the held-out ones (`dairy.json` 1.04, `sheep.json` 1.26). The
largest held-out misses are REML (1.45 with one trait under GBLUP, 1.24
with two under ssGBLUP): REML's cost is its number of iterations, which
depends on the data, not on the size. The exception is `l_records`
(three traits, two with repeated records, 3000 offspring per step for
ten steps): 2.1 times the model, all in the evaluation, whose cost per
equation is 5.7 µs against 0.6 to 2.1 µs in every other PBLUP scenario;
a smaller scenario of the same kind (`check/c_rep3_pblup.json`) is at
0.99, so this is the solver's convergence on that scenario rather than a
missing size term. Compilation is modelled to within about 1 s, worse
under REML, whose cold-minus-warm difference also carries the difference
in iterations between the two replicates.

Limits. The cubic terms of the dense equations and of the G and H
inverses fitted to zero: up to the largest calibration scenario (13 000
genotyped animals under GBLUP, 9800 under ssGBLUP, ``k n_g`` up to
13 000) the quadratic terms account for them, and beyond it the estimate
is likely low. The run time depends on the host and on the number of
threads; on others it is a relative guide only. The history term
depends on the processor caches.

The fit, as `fit_phases.jl` prints it:

```@eval
using Markdown, BnGsim
f = joinpath(pkgdir(BnGsim), "validation", "estimate", "fit_phases.txt")
Markdown.MD(Markdown.Code("text", read(f, String)))
```

Until 0.35 the time model had terms for the history, the evaluated
pedigree and the cube of the genotyped set only, fitted to scenarios with
small genomes and one trait; the examples `dairy.json` and `sheep.json`
ran about five times longer than their estimate. 0.36 removed two costs
that grew with the genome (G of GBLUP and ssGBLUP built from the chip
loci only, the allele counts of the summaries taken from the packed
haplotypes); 0.37 models the rest by phase.
