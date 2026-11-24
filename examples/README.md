# Example scenarios

Scenario files for `bngsim run` (schema `bngsim-scenario` version 1). Each
validates with `bngsim validate` and runs with
`bngsim run FILE --out DIR`; run times are on a 12-core desktop at 12
threads (BnGsim 0.36).

| File | What it shows | Run time |
| :--- | :--- | ---: |
| `minimal.json` | the smallest useful scenario: one trait, phenotypic truncation, three steps (the roadmap's smoke scenario) | 3 s |
| `dairy.json` | cattle preset; milk recorded on cows from age 2, up to five lactations with a permanent environment and herd × year groups; fertility with 20 % missing records; ssGBLUP with three genotyped cohorts; one calf per cow, at most five parities, survival by age, minimum-coancestry mating; genomes dropped once unusable | 55 s |
| `sheep.json` | sheep preset with crossover interference (ν = 4.3); a bottleneck in the history; a random-sign Gamma effect distribution; GBLUP on two traits; litters of at most two, at most 40 offspring per ram, at most two full sibs per pair | 55 s |
| `history-qtl.json` | three shared histories with a trait under stabilizing selection whose loci become the QTL (100 of them, breeding effects correlated 0.8 with the historical ones); OCS with ΔF = 1 % per cohort and minimum-coancestry mating; the aggregate reports the between-history standard error (`_seh`) | 11 s |

The resource estimate (`bngsim estimate FILE`) is calibrated on
single-trait scenarios; for the multi-trait examples with repeated records
its run time is two to three times too low (see the documentation
page "Resource estimate").
