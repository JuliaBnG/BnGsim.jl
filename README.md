# BnGsim.jl

Forward simulation of breeding programs in the
[JuliaBnG](https://juliabng.github.io/) ecosystem: Wright–Fisher founder
histories (FisherWright.jl), keyed meiosis on packed phased haplotypes,
several traits with repeated records and contemporary groups, genetic
evaluation (Breeding.jl), selection and mating, run from a versioned
JSON scenario or from Julia, with results independent of thread count and
resume history. Built on the BnGStructs.jl data types.

[![CI](https://github.com/JuliaBnG/BnGsim.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/BnGsim.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-juliabng.github.io-blue.svg)](https://juliabng.github.io/BnGsim/)

Requires Julia 1.10 (LTS) or later.

```julia
using Pkg
Pkg.add("BnGsim")
```

```julia
using BnGsim

s  = read_scenario(joinpath(pkgdir(BnGsim), "examples", "minimal.json"))
df = run_scenario(s)            # one row per replicate and time step
```

or `bngsim run minimal.json --out out` from a shell; scenarios are described in
the [documentation](https://juliabng.github.io/BnGsim/dev/scenarios/).

The development plan is in the [roadmap](https://juliabng.github.io/BnGsim/dev/roadmap/).
Release notes are in [`NEWS.md`](NEWS.md).
