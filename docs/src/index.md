# BnGsim.jl

BnGsim.jl simulates breeding programs in the JuliaBnG ecosystem. Version
0.1 provides the building blocks:

- **Founders:** Wright–Fisher burn-in through FisherWright.jl
  (`historical_population`, `continue_population`), founder sampling
  (`founder_haplotypes`, `sample_founder_population`), and readers and
  writers for packed haplotype files.
- **Meiosis:** Poisson crossovers without interference, `gamete!` and
  `gene_drop` on packed BnGStructs `Haplotype`s, `random_mate`, `rnd_ped`.
- **Traits:** QTL sampling (`sample_loci`, `eqtl`), true breeding values
  (`tbv!`), and phenotypes with correlated residuals, sex-limited expression
  and threshold traits (`phenotype!`).
- **Diagnostics and provenance:** `founder_summary`, `ld_decay`,
  `vanraden_g`, `provenance`.

Every stochastic function takes an explicit `rng`. See the [Roadmap](@ref)
for the planned scenario, evaluation and selection layers.

## Installation

BnGsim requires Julia 1.10 (LTS) or later.

```julia
using Pkg
Pkg.add("BnGsim")
```

## Quick start

```julia
using BnGStructs, BnGsim, Random

rng = Xoshiro(2026)
fp  = sample_founder_population(Cattle(100); n_id = 20, n_marker = 500, n_qtl = 50,
                                h2 = 0.3, burnin_generations = 200, rng = rng)
fp.tbv           # true breeding values of the founders
fp.markers.loci  # marker chip loci in fp.map
```
