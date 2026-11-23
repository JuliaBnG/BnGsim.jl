# Archives

An archive holds what other software needs to analyse a simulated
population: the pedigree, the records, the true and estimated breeding
values, the QTL with their effects, and genotypes. `bngsim run ...
--archive PANEL` writes one per replicate to `DIR/archive/r<N>/` when the
replicate finishes, and [`write_archive`](@ref) writes one from Julia.

| File | Content |
| :--- | :--- |
| `pedigree.csv` | `id, sire, dam, sex, generation, exit`; `0` for an unknown parent, `sex` `M`/`F`, `exit` the time step the animal died or was culled (empty while alive) |
| `records.csv` | `id` and one column per trait: the phenotype where it was recorded (expressed in the animal's sex, at or after the trait's age, alive then), empty otherwise |
| `truth.csv` | `id, tbv_<trait>..., ebv_<trait>...`: the true genetic values and the last evaluation |
| `traits.csv` | per trait: `name, h2, sigma_a, mu, sex, age, intercept` |
| `qtl.csv` | `locus, chr, bp` (within the chromosome) and the effect of one copy of the allele coded 1 for each trait |
| `genotypes.bed/.bim/.fam` | PLINK 1 binary genotypes of every animal whose genome is kept, at the loci of `PANEL` |

`PANEL` is `chip` (the scenario's chip), `qtl`, `all` (every segregating
locus of the founders) or `none` (no genotypes). The archives are listed
with their SHA-256 in `provenance.json`, and a stopped and resumed run
writes the same files, byte for byte.

In `genotypes.bed` the allele coded 1 is A1, the fifth column of
`genotypes.bim`. PLINK 1.9 counts it with `--recode A --keep-allele-order`;
without that flag it recodes A1 as the minor allele of the sample. plink2
counts it with `--export A --export-allele` on that column, and by default
counts the sixth. Both tools were checked to read back the simulated
genotypes exactly (200 animals, 300 loci).

## Infinitesimal breeding values

The breeding loop simulates finite loci. To test software on data whose
genetic model is exactly the animal model, [`infinitesimal`](@ref) draws
``u \sim N(0, G_0 \otimes A)`` and phenotypes ``y = \mu + u + e``,
``e \sim N(0, R_0)``, on any pedigree, by Mendelian sampling: each animal
gets the mean of its parents' values plus a term with covariance
``(1 - \tfrac14(1 + F_s) - \tfrac14(1 + F_d)) G_0``. Every animal draws from
its own [`SimKey`](@ref) stream, so its values do not depend on the order
of computation. The test suite checks the covariances against
``G_0 \otimes A`` on a pedigree with inbred animals.

```julia
using BnGsim
r = run_replicate(read_scenario("examples/minimal.json"))
ped = pedigree(r.population)
u, y = infinitesimal(ped.sire, ped.dam, [1.0;;], [3.0;;], SimKey(1))
```
