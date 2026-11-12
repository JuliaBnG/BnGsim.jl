# BnGsim.jl

Forward simulation of breeding programs in the
[JuliaBnG](https://juliabng.github.io/) ecosystem: Wright–Fisher founder
populations (FisherWright.jl), meiosis and gene dropping on packed phased
haplotypes, QTL and phenotype simulation, and founder diagnostics, all
built on the BnGStructs.jl data types.

[![CI](https://github.com/JuliaBnG/BnGsim.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/BnGsim.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-juliabng.github.io-blue.svg)](https://juliabng.github.io/BnGsim/)

Requires Julia 1.10 (LTS) or later.

```julia
using Pkg
Pkg.add("BnGsim")
```

```julia
using BnGStructs, BnGsim, Random

rng = Xoshiro(2026)
fp  = sample_founder_population(Cattle(100); n_id = 20, n_marker = 500, n_qtl = 50,
                                h2 = 0.3, burnin_generations = 200, rng = rng)
```

The development plan is in the [roadmap](https://juliabng.github.io/BnGsim/roadmap/).
Release notes are in [`NEWS.md`](NEWS.md).
