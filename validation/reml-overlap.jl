# Diagnostic (not a gate): REML σ²a with overlapping generations and
# records delayed to age 1, by selection rule and number of QTL.
#
#   julia -t 16 --project=<env with BnGsim, BnGStructs> validation/reml-overlap.jl
#
# Written when 2–4 replicates of the culling scenario of validation/reml.jl
# gave σ̂²a ≈ 0.75–0.85; with 50 replicates per configuration no
# configuration is biased by more than ~2 s.e., and the largest deficit
# vanishes with 2900 QTL (a finite-locus effect). Scenario otherwise as in
# validation/reml.jl (LE founders, h² = 0.3, 400 offspring, 20 sires,
# 80 dams, PBLUP with REML, 5 evaluated steps; step 0 random when records
# are delayed).
using BnGsim, BnGStructs, Random, Statistics, Printf
sp = GenericSpecies("val", [50_000_000, 50_000_000, 50_000_000])
function one(r; age, maxage, sel, nqtl)
    rng = Xoshiro(30_000 + r)
    vm = VariantMap(sort(unique(rand(rng, 1:150_000_000, 3000))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 800)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 1.0, age = age)
    qtl = eqtl(hap, sort(randperm(rng, length(vm))[1:nqtl]), t; rng = rng)
    pop = founder_population(hap, qtl, t, SimKey(2026; replicate = r))
    age > 0 && step!(pop, Scheme(; noff = 400, nsire = 20, ndam = 80, evaluation = :tbv, selection = :random,
        max_sire_age = maxage, max_dam_age = maxage), mm)
    sc = Scheme(; noff = 400, nsire = 20, ndam = 80, evaluation = :pblup, selection = sel, reml = true,
        max_sire_age = maxage, max_dam_age = maxage)
    simulate!(pop, sc, mm, 5).va[5]
end
nrep = 50
println("# BnGsim ", pkgversion(BnGsim), ", Julia ", VERSION, ", ", Threads.nthreads(), " threads, ", nrep, " replicates per configuration")
for (lab, kw) in (("age1 max2 trunc 300", (age = 1, maxage = 2, sel = :truncation, nqtl = 300)),
                  ("age1 max2 random 300", (age = 1, maxage = 2, sel = :random, nqtl = 300)),
                  ("age1 max2 trunc 2900", (age = 1, maxage = 2, sel = :truncation, nqtl = 2900)),
                  ("age1 max0 trunc 300", (age = 1, maxage = 0, sel = :truncation, nqtl = 300)),
                  ("age0 max2 trunc 300", (age = 0, maxage = 2, sel = :truncation, nqtl = 300)))
    v = [one(r; kw...) for r = 1:nrep]
    @printf("%-22s σ²a mean %.4f s.e. %.4f\n", lab, mean(v), std(v) / sqrt(nrep))
    flush(stdout)
end
