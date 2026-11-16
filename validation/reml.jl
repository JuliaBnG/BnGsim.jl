# Gate for re-estimated variance components: REML under selection.
#
#   julia -t 16 --project=<env with BnGsim, BnGStructs> validation/reml.jl [replicates] [truncation|random|cull]
#
# Claim tested: with all records and the full pedigree back to the
# founders, REML estimates of σ²a and σ²e are unbiased for the founder
# (base) values although every generation is selected on the EBV that the
# estimates feed (Sorensen & Kennedy 1984; Henderson's selection
# argument). The genetic variance within the last cohort is reported next
# to it: it falls under selection (Bulmer effect and drift), and REML is
# not estimating it.
#
# Scenario (stated in full, as the reference for the gate):
# - genome: 3 chromosomes of 50 Mb, 1 M each; 3000 loci at random
#   positions; founders: 400 individuals with independent random
#   haplotypes, allele frequency ½ (linkage equilibrium in expectation);
# - trait: additive, 300 QTL among the loci, normal effects scaled so
#   that σ²a = 1 exactly in the founders, h² = 0.3, so σ²e = 7/3;
# - scheme: 400 offspring per step, 20 sires and 80 dams by truncation on
#   PBLUP EBV (or, as the control, at random), discrete generations,
#   5 steps;
# - `cull`: the same trait recorded at age 1 (so an animal that dies or is
#   culled before then has no record), parents up to age 2, survival 0.9
#   per step, unselected adults culled; selection and culling on EBV make
#   records missing depending on earlier records, which REML with all
#   data and the pedigree ignores correctly (missing at random given the
#   data; Im, Fernando & Gianola 1989); step 0 is random mating, as there
#   are no records yet, then 5 evaluated steps; REML (AI) re-estimates the
#   components before every evaluation from step 1 on (at step 0 only the
#   unrelated founders are recorded and the components are not separable);
# - estimate: the components used at the last step (2000 records:
#   the founders and cohorts 1 to 4).
# Pass: |mean − true| < 3 s.e. for both, s.e. over replicates.
using BnGsim, BnGStructs, Random, Statistics, Printf

nrep = parse(Int, get(ARGS, 1, "100"))
sel = Symbol(get(ARGS, 2, "truncation"))
sp = GenericSpecies("val", [50_000_000, 50_000_000, 50_000_000])
va, ve, vlast, h2 = zeros(nrep), zeros(nrep), zeros(nrep), zeros(nrep)
t0 = time()
for r = 1:nrep
    rng = Xoshiro(30_000 + r)
    vm = VariantMap(sort(unique(rand(rng, 1:150_000_000, 3000))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 800)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 1.0, age = sel == :cull ? 1 : 0)
    qtl = eqtl(hap, sort(randperm(rng, length(vm))[1:300]), t; rng = rng)
    pop = founder_population(hap, qtl, t, SimKey(2026; replicate = r))
    sc = sel == :cull ?
         Scheme(; noff = 400, nsire = 20, ndam = 80, evaluation = :pblup, reml = true, max_sire_age = 2,
             max_dam_age = 2, survival = [0.9], cull_unselected = true) :
         Scheme(; noff = 400, nsire = 20, ndam = 80, evaluation = :pblup, selection = sel, reml = true)
    # records at age 1: none at step 0, so the founders are mated at random
    sel == :cull && step!(pop, Scheme(; noff = 400, nsire = 20, ndam = 80, evaluation = :tbv,
        selection = :random, max_sire_age = 2, max_dam_age = 2, survival = [0.9]), mm)
    df = simulate!(pop, sc, mm, 5)
    va[r], ve[r] = df.va[5], df.ve[5]
    h2[r] = va[r] / (va[r] + ve[r])
    vlast[r] = df.var[end]
    @printf(stderr, "replicate %d of %d: σ²a %.3f σ²e %.3f, var(last cohort) %.3f (%.0f s)\n",
        r, nrep, va[r], ve[r], vlast[r], time() - t0)
end

println("# BnGsim ", pkgversion(BnGsim), ", Julia ", VERSION, ", ", Threads.nthreads(), " threads, ",
    nrep, " replicates, selection ", sel, ", ", round(time() - t0; digits = 0), " s")
se(x) = std(x) / sqrt(length(x))
ok = true
for (name, x, truth) in (("σ²a", va, 1.0), ("σ²e", ve, 7 / 3), ("h²", h2, 0.3))
    pass = abs(mean(x) - truth) < 3se(x)
    name == "h²" || (global ok &= pass)
    @printf("%-4s REML mean %.4f (s.e. %.4f), true %.4f, %s\n", name, mean(x), se(x), truth,
        pass ? "within 3 s.e." : "BIASED")
end
@printf("genetic variance within the last cohort %.4f (s.e. %.4f): not a REML target\n",
    mean(vlast), se(vlast))
println(ok ? "gate passed" : "gate NOT passed")
