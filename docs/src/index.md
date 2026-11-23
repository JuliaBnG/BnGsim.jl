# BnGsim.jl

BnGsim.jl simulates breeding programs in the JuliaBnG ecosystem, from a
scenario file or from Julia:

- **Scenarios and the command line:** a versioned JSON scenario
  ([Scenarios](@ref)) run by `bngsim run` or `run_scenario`, with
  provenance, a resource estimate, stop and resume, replicate
  aggregation and archives for other software ([Command line](@ref)).
- **Founders:** Wright–Fisher histories through FisherWright.jl, with
  demography and selection, shared and stored histories
  (`write_histories`), and founder sampling (`founder_haplotypes`).
- **Keyed streams:** `SimKey` addresses every random stream by its
  coordinate in the simulation, so results do not depend on thread count,
  call order or resume history (`keyed_rng`).
- **Meiosis:** `MeiosisMap` and `make_progeny!` (keyed, word-wise,
  sex-specific maps, crossover interference, optional crossover records).
- **Traits:** QTL sampling (`sample_loci`, `eqtl`, `scale_qtl`) with
  normal, Gamma, Laplace and spike-and-slab effects, several correlated
  traits, sex- and age-limited, repeated and missing records,
  contemporary groups.
- **Breeding loop:** `founder_population`, `Scheme` and `simulate!`:
  evaluation through Breeding.jl (true values, phenotype, PBLUP, GBLUP,
  single-step GBLUP, optional REML), selection by truncation, optimum
  contributions or at random, mating policies, reproductive limits,
  overlapping generations, survival and culling, per-step summaries, and
  accessors for the population (`pedigree`, `tbv`, `ebv`, `phenotypes`,
  `records`, `genomes`).
- **Diagnostics:** `generation_summary`, `founder_summary`, `ld_decay`,
  `founder_segments`, `irm`, `pedigree_view`.

The building blocks of 0.1–0.3 (gene drop with an explicit `rng`,
`historical_population`, `sample_founder_population`, `vanraden_g`, the
text `provenance` header, …) are in `BnGsim.Legacy`, outside the keyed
streams and the 1.0 promise (see the [API Reference](@ref)).

## A breeding program

```julia
using BnGStructs, BnGsim, Random

sp   = GenericSpecies("toy", [20_000_000, 10_000_000])
rng  = Xoshiro(1)
vm   = VariantMap(sort(unique(rand(rng, 1:30_000_000, 2000))), sp)
hap  = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 400)))      # 200 founders
t    = Trait("y"; h² = 0.3, σₐ = 1.0)
qtl  = eqtl(hap, 1:5:length(vm), t; rng = rng)
pop  = founder_population(hap, qtl, t, SimKey(2026))
sc   = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup)
simulate!(pop, sc, MeiosisMap(sp, vm), 10)    # one row per generation
```

## Installation

BnGsim requires Julia 1.10 (LTS) or later.

```julia
using Pkg
Pkg.add("BnGsim")
```

## Quick start

```julia
using BnGsim

s  = read_scenario(joinpath(pkgdir(BnGsim), "examples", "minimal.json"))
df = run_scenario(s)            # one row per replicate and time step
```

or, from a shell, `bngsim run minimal.json --out out` (see [Command line](@ref)).
