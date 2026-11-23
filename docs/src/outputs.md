# Output files

The files that `bngsim` writes are versioned: `provenance.json` and
`bngsim version --json` carry `outputs_version`, now `1`. Within a
version, a file keeps its columns, in this order, for every scenario;
new columns are only appended at the end of a table, so a reader that
selects columns by name keeps working, and one that reads them by
position keeps working for the columns it knows. A change that renames,
removes, moves or redefines a column raises `outputs_version`.

Which optional columns a table has depends on the scenario (the number
of traits, repeated records, the selection method, recorded meioses), as
listed below; for one scenario the set is fixed.

## Run directory

| File | Content |
| :--- | :--- |
| `scenario.json` | the canonical scenario ([`scenario_json`](@ref)) |
| `summary.csv` | one row per replicate and time step (below) |
| `aggregate.csv` | one row per time step, over the replicates (below) |
| `provenance.json` | the run: versions, hashes, host, sessions ([Command line](cli.md)) |
| `archive/r<N>/` | with `--archive`: the archive of replicate `N` (below) |

`provenance.json` holds `scenario_sha256`, the SHA-256 of the canonical
text, and `results_sha256`, that of the content that determines the
results ([`results_sha256`](@ref)): the canonical scenario without its
name, number of replicates, schema keys and history directory, and
without every key at its default. A key added in a later release, with a
default that keeps earlier results, changes the canonical text of every
scenario but not its `results_sha256`; `bngsim merge` accepts runs whose
scenarios agree on it.

## `summary.csv`

Numbers are in Julia's shortest round-trip form; `NaN` where a value is
undefined. Columns in order ([`generation_summary`](@ref), [`step!`](@ref)):

| Column | Present | Content |
| :--- | :--- | :--- |
| `replicate` | always | the replicate |
| `history` | always | the historical population of the replicate (its index among `founders.histories`; the replicate without shared histories) |
| `generation` | always | the time step (the cohort born at it) |
| `n` | always | animals in the cohort |
| `mean`, `var`, `genic` | always | mean, variance and genic variance ``\sum 2p(1-p)a^2`` of the genetic value (several traits: of the index ``w'g``) |
| `qtl_seg`, `qtl_fixed`, `qtl_lost` | always | QTL at which the favourable allele segregates, is fixed or lost |
| `qtl_maf`, `qtl_alpha`, `qtl_top10` | always | mean MAF of the segregating QTL, slope of ``\log a^2`` on ``\log 2p(1-p)``, share of the genic variance in the largest tenth |
| `mean_<t>`..., `var_<t>`..., `rg_<t>_<u>`... | several traits | per trait mean and variance, then the genetic correlation of each pair |
| `F`, `coancestry` | always | mean pedigree inbreeding and coancestry |
| `Fhet`, `Fgrm` | always | loss of heterozygosity; mean genomic inbreeding (VanRaden, founder base) |
| `alive` | always | animals alive at the step |
| `Fseg` | `record_meiosis` | mean realized inbreeding from founder segments |
| `accuracy` | always | accuracy of the index of the candidates |
| `accuracy_<t>`... | several traits | per trait |
| `reliability` | always | mean reliability of the index (`NaN` unless `report_reliability`) |
| `reliability_<t>`... | several traits | per trait |
| `nsire`, `ndam` | always | parents used |
| `L` | always | generation interval |
| `intensity_sire`, `intensity_dam` | always | realized selection intensities |
| `va`, `ve` | one trait | variances used by the evaluation; with repeated records `va`, `vpe`, `ve` |
| `va_<t>`..., `ve_<t>`..., `vpe_<t>`... | several traits | per trait; `vpe_` for the traits with repeated records |
| `ocs_bound`, `ocs_kinship`, `ocs_L`, `ocs_feasible` | `selection = "ocs"` | the bound, the kinship reached, the interval used, whether the bound was met |

The last row of a replicate (the newest cohort, in which no selection
took place) has `NaN` for the selection columns and `0` parents.

## `aggregate.csv`

One row per time step ([`aggregate`](@ref)): `generation`, `replicates`,
then for each column `c` of `summary.csv` after `generation` (in that
order, `replicate` and `history` excepted) `c_mean`, `c_sd`, `c_se`;
then the same three for the changes within a replicate: `dG` and
`dG_<t>`... (several traits), `dF`, `dFhet`, `var_ratio`,
`genic_ratio`. When the replicates aggregated include at least two
histories with several replicates each (so this depends on the
replicates, not only on the scenario), also `histories` (after
`replicates`) and `c_sdh`, `c_seh` after each `c_se`.

## Archive (`archive/r<N>/`)

| File | Columns |
| :--- | :--- |
| `pedigree.csv` | `id, sire, dam, sex, generation, exit, herd` |
| `records.csv` | `id`, one column per trait (the phenotype, or the mean of the records, recorded by the last step; empty otherwise) |
| `truth.csv` | `id`, then per trait `tbv_<t>`, `ebv_<t>` |
| `traits.csv` | `name, h2, sigma_a, mu, sex, age, intercept, c2, records, repeatability` |
| `repeated.csv` | with repeated records: `id, trait, step, value` |
| `qtl.csv` | `locus, chr, bp`, then the effect of each trait |
| `genotypes.bed/.bim/.fam` | PLINK 1 binary, at the loci of the panel |

Details in [`write_archive`](@ref) and [Archives](archive.md).

## History directory

`bngsim history` writes one file per history in FisherWright's result
format and `manifest.json` (schema `"bngsim-history"`, its own version
`1`): the parameters that determine the histories, the seed, and the
file names with their SHA-256 ([`write_histories`](@ref)).
