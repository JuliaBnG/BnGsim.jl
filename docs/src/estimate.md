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
| dense evaluation | one ``n_g \times n_g`` matrix, ``8 n_g^2``, for the genotyped animals | 13.3 (`:gblup`), 19.0 (`:ssgblup`) |

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

Run time is the alleles carried through meiosis in the historical run
(about ``2N\mu L \sum_g \min(g, 4N)``), plus, per step, terms in the
evaluated pedigree and the cube of the genotyped set, plus a fixed cost of
the dense evaluations. It is a rough guide: within a factor of about two
on the calibration host.

## Calibration

The coefficients were fitted on one host: AMD Threadripper 3970X, Linux
x86_64, Julia 1.13.4, 8 threads, BnGsim 0.14. Sixteen calibration
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

The run time of large genomes and of multi-trait scenarios with repeated
records is not modelled: the calibration scenarios had small genomes, and
the model has no term for sampling the founders from the history, for the
per-step work over all loci (meiosis and the summaries), or for
multi-trait equations. Until 0.35 the examples `dairy.json` and
`sheep.json` ran about five times longer than their estimate; 0.36
removed two costs that grew with the genome (G of GBLUP and ssGBLUP built
from the chip loci only, the allele counts of the summaries taken from
the packed haplotypes), and at 12 threads on a 12-core desktop they ran
in 53 s and 55 s, against estimates of 24.8 s and 17.4 s (2.1 and 3.2
times). A recalibration with these terms is planned.
