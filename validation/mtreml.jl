# Gate for several traits: multi-trait REML under index selection.
#
#   julia -t 16 --project=<env with BnGsim, BnGStructs> validation/mtreml.jl [replicates] [truncation|random]
#
# Claim tested: with all records and the full pedigree, multi-trait REML
# estimates of the genetic (G0) and residual (R0) covariance matrices are
# unbiased for the founder values although every generation is selected on
# an index of the EBVs they feed, and although one trait is recorded in
# females only (missing at random given sex). Unbiasedness under
# selection is an infinitesimal-model result (Sorensen & Kennedy 1984);
# the trait here has 300 QTL, so a small bias of the size seen in
# validation/reml.jl (1–2 % in σ²a) is possible. Random selection is the
# control.
#
# Scenario (stated in full, as the reference for the gate):
# - genome: 3 chromosomes of 50 Mb, 1 M each; 3000 loci at random
#   positions; founders: 400 individuals with independent random
#   haplotypes, allele frequency ½ (linkage equilibrium in expectation);
# - traits: `a` (h² = 0.3, both sexes) and `b` (h² = 0.2, females only),
#   σₐ = 1 each, additive, 300 shared QTL with genetic correlation 0.5
#   exact in the founders (`eqtl`), residual correlation 0.2; so
#   G0 = [1 0.5; 0.5 1] and R0 = [7/3 0.2√(7/3·4); · 4];
# - scheme: 400 offspring per step, 20 sires and 80 dams by truncation on
#   the index â_a + â_b of multi-trait PBLUP EBVs (or at random), discrete
#   generations, 5 steps; REML from step 1 on;
# - estimate: the components used at the last step.
# Pass: every element within 3 s.e. of its true value (six comparisons per
# selection rule, so a single excursion near 3 s.e. is not by itself
# evidence of bias).
using BnGsim, BnGStructs, Random, Statistics, Printf

nrep = parse(Int, get(ARGS, 1, "200"))
sel = Symbol(get(ARGS, 2, "truncation"))
sp = GenericSpecies("val", [50_000_000, 50_000_000, 50_000_000])
Rg = [1 0.5; 0.5 1]
Re = [1 0.2; 0.2 1]
ts = [Trait("a"; h² = 0.3), Trait("b"; h² = 0.2, sex = 0)]
σe = sqrt.([0.7 / 0.3, 0.8 / 0.2])
G0, R0 = Rg, σe .* Re .* σe'
est = zeros(nrep, 6)                 # G11 G12 G22 R11 R12 R22
t0 = time()
for r = 1:nrep
    rng = Xoshiro(40_000 + r)
    vm = VariantMap(sort(unique(rand(rng, 1:150_000_000, 3000))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 800)))
    mm = MeiosisMap(sp, vm)
    qtl = eqtl(hap, sort(randperm(rng, length(vm))[1:300]), Rg, ts; rng = rng)
    pop = founder_population(hap, qtl, ts, SimKey(2026; replicate = r); Re = Re)
    sc = Scheme(; noff = 400, nsire = 20, ndam = 80, evaluation = :pblup, selection = sel,
        weights = [1.0, 1.0], reml = true)
    simulate!(pop, sc, mm, 5)
    G, R = pop.G0, pop.R0
    est[r, :] = [G[1, 1], G[1, 2], G[2, 2], R[1, 1], R[1, 2], R[2, 2]]
    @printf(stderr, "replicate %d of %d: %s (%.0f s)\n", r, nrep, round.(est[r, :]; digits = 3), time() - t0)
end

println("# BnGsim ", pkgversion(BnGsim), ", Julia ", VERSION, ", ", Threads.nthreads(), " threads, ",
    nrep, " replicates, selection ", sel, ", ", round(time() - t0; digits = 0), " s")
truth = [G0[1, 1], G0[1, 2], G0[2, 2], R0[1, 1], R0[1, 2], R0[2, 2]]
names = ["G_aa", "G_ab", "G_bb", "R_aa", "R_ab", "R_bb"]
ok = true
for k = 1:6
    x = est[:, k]
    se = std(x) / sqrt(nrep)
    z = (mean(x) - truth[k]) / se
    pass = abs(z) < 3
    global ok &= pass
    @printf("%-5s REML mean %.4f (s.e. %.4f), true %.4f, z = %+.2f  %s\n", names[k], mean(x), se, truth[k], z,
        pass ? "" : "BIASED")
end
rg = est[:, 2] ./ sqrt.(est[:, 1] .* est[:, 3])
@printf("r_g   mean %.4f (s.e. %.4f), true 0.5000\n", mean(rg), std(rg) / sqrt(nrep))
println(ok ? "gate passed" : "gate NOT passed")
