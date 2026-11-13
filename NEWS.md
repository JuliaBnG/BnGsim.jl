# BnGsim.jl release notes

## v0.1.2 (unreleased)

Requires BnGStructs 0.3.

- `gene_drop!(og, pg, pm, species, map)`, `crossovers(species, map)` and
  `expand` check that the `VariantMap` is on the genome of the species
  (`BnGStructs.check_coordinates`). A map with within-chromosome
  positions used to give wrong recombinants without an error.
- `read_packed_haplotypes` reads through `BnGStructs.read_haplotypes`,
  which now owns the `BNGH1` format; its result is unchanged.
- Record columns are named through `BnGStructs.trait_columns`; the names
  (`tbv_<trait>`, `ft_<trait>`) are unchanged.

## v0.1.1 (unreleased)

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
