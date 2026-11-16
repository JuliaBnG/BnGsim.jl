# BnGsim.jl release notes

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
