# Phase 4 gate: accuracy of selection candidates by evaluation method.
#
#   julia -t 16 --project=<env with BnGsim, BnGStructs> validation/accuracy.jl [replicates]
#
# Expected ordering on this scenario: ρ(phenotype) < ρ(PBLUP) < ρ(GBLUP),
# with each margin larger than its Monte-Carlo standard error.
#
# Scenario (stated in full, as the reference for the gate):
# - genome: 3 chromosomes of 50 Mb, 1 M each (1 cM/Mb);
# - founders: 300 individuals from a Fisher–Wright population (Ne = 300,
#   1000 generations, μ = 1e-8 per bp per generation via `mr = 1`), so
#   founders carry LD; loci with MAF ≥ 0.05 in the founders;
# - trait: additive, 200 QTL drawn from the loci, normal effects,
#   h² = 0.25, σₐ = 1 in the founders;
# - chip: all non-QTL loci (dense panel); GBLUP blends G with 5 % of A;
# - scheme: 300 offspring per generation, 15 sires and 60 dams selected by
#   truncation on the method's EBV, 5 generations;
# - accuracy: cor(EBV, TBV) among the candidates of generations 1 to 4,
#   each candidate with its own record, averaged per replicate.
# Every method runs on the same founders, QTL and keys within a replicate,
# so the comparison is paired; the standard error is that of the paired
# differences over replicates.
using BnGsim, BnGStructs, FisherWright, Random, Statistics, Printf

nrep = parse(Int, get(ARGS, 1, "20"))
sp = GenericSpecies("val", [50_000_000, 50_000_000, 50_000_000])
methods = (:phenotype, :pblup, :gblup)
acc = zeros(nrep, length(methods))
t0 = time()
for r = 1:nrep
    res = fisher_wright(300, 1000, Int.(sp.chromosome), 1.0; M = 1e8, result = true, rng = Xoshiro(r))
    xy, lmp = muts2bitarray(res.active_haplotypes, res.chromosome_ends)
    frq = vec(sum(xy; dims = 2)) ./ size(xy, 2)
    keep = findall(f -> 0.05 <= min(f, 1 - f), frq)
    hap = Haplotype(BitMatrix(xy[keep, :]))
    vm = VariantMap(UInt32.(lmp.pos[keep]), sp)
    mm = MeiosisMap(sp, vm)
    rng = Xoshiro(10_000 + r)
    qloci = sort(randperm(rng, length(keep))[1:200])
    chip = setdiff(1:length(keep), qloci)
    t = Trait("y"; h² = 0.25, σₐ = 1.0)
    qtl = eqtl(hap, qloci, t; rng = rng)
    for (k, ev) in enumerate(methods)
        pop = founder_population(hap, qtl, t, SimKey(2026; replicate = r))
        sc = Scheme(; noff = 300, nsire = 15, ndam = 60, evaluation = ev, chip = chip)
        df = simulate!(pop, sc, mm, 5)
        acc[r, k] = mean(df.accuracy[2:5])
    end
    @printf(stderr, "replicate %d of %d: %s (%.0f s)\n", r, nrep, round.(acc[r, :]; digits = 3), time() - t0)
end

println("# BnGsim ", pkgversion(BnGsim), ", Julia ", VERSION, ", ", Threads.nthreads(), " threads, ",
    nrep, " replicates, ", round(time() - t0; digits = 0), " s")
for (k, ev) in enumerate(methods)
    @printf("%-10s mean accuracy %.4f  (s.e. %.4f)\n", ev, mean(acc[:, k]), std(acc[:, k]) / sqrt(nrep))
end
ok = true
for (a, b) in ((1, 2), (2, 3))
    d = acc[:, b] .- acc[:, a]
    se = std(d) / sqrt(nrep)
    pass = mean(d) > 2se
    global ok &= pass
    @printf("%-10s − %-10s %.4f  (s.e. %.4f, %s)\n", methods[b], methods[a], mean(d), se,
        pass ? "margin > 2 s.e." : "NOT SHOWN")
end
println(ok ? "gate passed" : "gate NOT passed")
