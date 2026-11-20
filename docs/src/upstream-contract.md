# Upstream contract

BnGsim builds on three JuliaBnG packages besides the shared types of
BnGStructs: FisherWright (historical populations), RelationshipMatrices
(relationship matrices and realized IBD) and Breeding (evaluation and
selection). It imports each of them in one place only, an adapter module in
`src/adapters/` that re-exports the names BnGsim uses and absorbs
differences between their versions. A test (`test/adapters.jl`,
"Upstream contract: layering") fails if any other source file imports or
names them, and the contract tests below fail if a version of an upstream
package breaks a convention BnGsim relies on. BnGStructs is the shared
vocabulary of all packages (`Haplotype`, `VariantMap`, `Species`,
`TraitQTL`, `SNPSet`, `LocusSet`) and is used directly. Result tables and
files that BnGsim writes use types and formats BnGsim owns (CSV, JSON,
PLINK).

## FisherWright (`FisherWrightAdapter`)

Names: `fisher_wright`, `FisherWrightResult`, `muts2bitarray`, and
`census`, `dense_haplotypes` (version-independent).

- A population of `N` diploids is `2N` haplotypes; individual `i` is
  haplotypes `2i - 1` and `2i`.
- A haplotype is a sorted, unique `Vector{UInt32}` of the positions of the
  derived alleles it carries: genome-wide, 1-based coordinates on the
  chromosomes laid end to end (the `VariantMap` convention of BnGStructs).
- `chromosome_ends` is `cumsum` of the chromosome lengths, as `UInt32`;
  `substitutions` are the sorted positions fixed and removed.
- A run with the same `rng` gives the same population; `init` with
  `popsize` continues a run and appends to its census record.
- `census(res)`: the census sizes `N₀, N₁, …` (`res.records.census` from
  FisherWright 0.10, `res.census` before).
- `dense_haplotypes(muts, ends; rng)`: the loci segregating among `muts`
  as a `BitMatrix` (loci × haplotypes) and a map with `chr` (1-based),
  genome-wide `pos` and the allele frequency `frq`; the placeholder
  nucleotides come from `rng` from FisherWright 0.10 on (BnGsim passes a
  copy, so the caller's stream does not advance), from the global
  generator before.
- The species method reads the genetic length of each chromosome from
  `Species.morgan` (FisherWright 0.8.2 and later).

## RelationshipMatrices (`RelationshipsAdapter`)

Names: `nrm`, `nrm_diag`, `ainv`, `grm`, `hinv`, `irm`, `Meioses`,
`HaplotypeSegments`, `founder_segments` (extended by BnGsim for a
`Population`), `segment_inbreeding`.

- A pedigree is two vectors `sire`, `dam` of row numbers, `0` for an
  unknown parent, parents before their offspring.
- `nrm(sire, dam[, ids])` is A (or its rows and columns `ids`),
  `nrm_diag` its diagonal `1 + F`, and `ainv` its inverse, sparse.
- `grm(hap; p, loci)` is VanRaden's (2008) first G, `ZZ' / 2Σpq` with
  `Z = M - 2p` for the dosages `M` of the haplotype pairs at `loci`.
- `hinv(sire, dam, G, genotyped)` is
  `A⁻¹ + [0 0; 0 G⁻¹ - A₂₂⁻¹]` for the genotyped rows.

## Breeding (`BreedingAdapter`)

Names: `MME`, `mme`, `RandomEffect`, `solve_mme`, `incidence`, `TM1997`,
`feasible`, `rocs`, `pev_candidates`, `reml`, `MultiTrait`, `mt_stack`,
`mt_blocks`, `mt_rinv`.

- `RandomEffect(Z, Kinv, σ²)` takes the inverse of the covariance
  structure of the effect; `mme(y, X, effects; σ²e)` builds Henderson's
  equations and `solve_mme` returns the fixed effects, then the random
  ones, in one vector, equal to a dense solve of those equations.
- `TM1997(u, sex, A, K)` returns contributions `c` summing to 1, `cAc`,
  and a feasibility flag read by `feasible`.

## Contract tests

`test/adapters.jl` checks each convention above on small cases with exact
answers: the haplotype layout and coordinates, the census across versions,
the map of `dense_haplotypes`, `ainv = inv(nrm)`, `grm` against the formula
of VanRaden, `hinv` against its definition, the solution of `mme` against
Henderson's equations solved densely, and `TM1997`'s contributions. A new
version of an upstream package that changes one of them fails these tests
before it changes a simulation.
