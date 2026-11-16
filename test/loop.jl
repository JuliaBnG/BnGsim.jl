# The breeding loop: reproducibility, the evaluation against an
# independent solve, the selection rules, and the response to selection
# against exact theory.
#
# Scenario: two chromosomes (20 Mb, 1 M; 10 Mb, 0.5 M), ~2000 loci, 200
# founders with random haplotypes (linkage equilibrium in expectation),
# one additive trait with h² = 0.3 and σₐ = 1 in the founders, QTL every
# fifth locus. Seeds are fixed; statistical checks use 4 standard errors.

@testset "Breeding loop" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(20261116)
    pos = unique(sort(rand(rng, 1:30_000_000, 2000)))
    vm = VariantMap(pos, sp)
    nlc = length(vm)
    hap = Haplotype(BitMatrix(rand(rng, Bool, nlc, 400)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 1.0, μ = 10.0)
    qtl = eqtl(hap, 1:5:nlc, t; rng = Xoshiro(2))
    key = SimKey(20261116)
    run(sc, ngen; k = key) = (p = founder_population(hap, qtl, t, k); (p, simulate!(p, sc, mm, ngen)))

    # ---- founders ---------------------------------------------------------
    pop = founder_population(hap, qtl, t, key)
    @test length(pop) == 200 && ngenerations(pop) == 1
    @test pop.tbv ≈ tbv(hap, qtl)
    @test std(pop.tbv) ≈ 1.0 atol = 1e-8
    @test any(pop.male) && !all(pop.male)
    @test pop.σe ≈ sqrt(0.7 / 0.3)
    @test founder_population(hap, qtl, t, key).phenotype == pop.phenotype     # keyed
    @test founder_population(hap, qtl, t, SimKey(key; replicate = 1)).phenotype != pop.phenotype
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, evaluation = :blup)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, selection = :ocs)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, evaluation = :gblup)

    # ---- reproducibility and bookkeeping ---------------------------------
    sc = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :phenotype)
    p1, d1 = run(sc, 4)
    p2, d2 = run(sc, 4)
    @test p1.sire == p2.sire && p1.dam == p2.dam && p1.tbv == p2.tbv && p1.phenotype == p2.phenotype
    @test isequal(d1, d2)                     # the last accuracy is NaN
    @test nrow(d1) == 5 && d1.generation == 0:4
    @test length(p1) == 1000 && collect(generation_ids(p1, 2)) == 401:600
    # parents come from the previous generation, as many as selected
    for g = 1:4
        ids = generation_ids(p1, g)
        @test all(p1.generation[p1.sire[ids]] .== g - 1) && all(p1.generation[p1.dam[ids]] .== g - 1)
        @test all(p1.male[p1.sire[ids]]) && !any(p1.male[p1.dam[ids]])
        @test length(unique(p1.sire[ids])) == 10 && length(unique(p1.dam[ids])) == 50
    end
    # golden checksum of a run without a solver (keys, meiosis, phenotypes,
    # truncation, pairing together)
    # (FNV-style fold, not Base.hash, which may change between Julia versions)
    loopsum(p) = foldl((a, x) -> a * 0x100000001b3 ⊻ (x % UInt64),
        vcat(Int64.(p.sire), Int64.(p.dam), round.(Int64, p.tbv .* 1e6)); init = 0xcbf29ce484222325)
    @test loopsum(p1) == GOLDEN_LOOP

    # ---- truncation picks the best by (criterion, id) ----------------------
    pop = founder_population(hap, qtl, t, key)
    evaluate!(pop, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv))
    s, ns, d, nd = select_parents(pop, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv), nothing)
    males = findall(pop.male)
    @test s == sort(males[sortperm(pop.tbv[males]; rev = true)[1:10]])
    @test sum(ns) == 200 && sum(nd) == 200 && all(ns .== 20) && all(nd .== 4)
    # ties are broken by ascending id
    pop.ebv .= 0.0
    s0, = select_parents(pop, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv), nothing)
    @test s0 == males[1:10]

    # ---- PBLUP against an independent dense solve -------------------------
    pb, _ = run(Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup), 2)
    m = evaluate!(pb, Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup))
    n = length(pb)
    A = nrm(pb.sire, pb.dam)
    λ = pb.σe^2 / t.σₐ^2
    Z = Matrix{Float64}(I, n, n)
    X = ones(n, 1)
    C = [X'X X'Z; Z'X Z'Z+λ*inv(A)]
    sol = C \ [X'pb.phenotype; Z'pb.phenotype]
    @test pb.ebv ≈ sol[2:end] rtol = 1e-8
    # GBLUP runs on all loci and predicts the genetic values
    gb = Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :gblup, chip = 1:nlc, blend = 0.0)
    evaluate!(pb, gb)
    @test cor(pb.ebv, pb.tbv) > 0.3
    # a chip that is a subset of the loci, with founder allele frequencies
    gs = Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :gblup, chip = 2:2:nlc)
    evaluate!(pb, gs)
    @test cor(pb.ebv, pb.tbv) > 0.3

    # ---- OCS and ROCS -------------------------------------------------------
    po, _ = run(Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv), 2)
    cand = collect(generation_ids(po, 2))
    A = nrm(po.sire, po.dam, cand)
    evaluate!(po, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv))
    K = 0.25
    ocs = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv, selection = :ocs, K = K)
    s, ns, d, nd = select_parents(po, ocs, nothing)
    @test sum(ns) == 200 && sum(nd) == 200
    c = zeros(length(cand))
    for (ids, k) in ((s, ns), (d, nd)), (i, x) in zip(ids, k)
        c[findfirst(==(i), cand)] = x / 400
    end
    @test c' * A * c <= K + 0.01                  # after rounding to whole offspring
    # rocs with λ = κ = 0 and a cap of one truncation share is truncation
    tr = select_parents(po, Scheme(; noff = 200, nsire = 10, ndam = 10, evaluation = :tbv), nothing)
    ro = select_parents(po, Scheme(; noff = 200, nsire = 10, ndam = 10, evaluation = :tbv,
        selection = :rocs, λ = 0.0, κ = 0.0, ub = 0.05), nothing)
    @test ro[1] == tr[1] && ro[3] == tr[3] && ro[2] == tr[2] && ro[4] == tr[4]

    # ---- response to selection against theory ------------------------------
    # Under an additive model the expected genetic value of an offspring is
    # the mean of its parents', exactly; the deviation of the generation mean
    # is Mendelian sampling, with variance Var(MS)/n.
    for ev in (:tbv, :phenotype, :pblup)
        p, df = run(Scheme(; noff = 400, nsire = 10, ndam = 50, evaluation = ev), 3)
        for g = 1:3
            ids = generation_ids(p, g)
            ms = p.tbv[ids] .- (p.tbv[p.sire[ids]] .+ p.tbv[p.dam[ids]]) ./ 2
            @test abs(mean(ms)) < 4 * std(ms) / sqrt(length(ids))
        end
        @test df.mean[end] > df.mean[1]                        # selection responds
        ev == :tbv && @test all(df.accuracy[1:end-1] .≈ 1)
    end
    # random selection: no expected response; drift raises F and Fhet
    pr, dr = run(Scheme(; noff = 200, nsire = 5, ndam = 20, evaluation = :tbv, selection = :random), 5)
    @test abs(dr.mean[end] - dr.mean[1]) < 4 * sqrt(5 * 1.0 / 10)   # loose: drift variance ≈ t σ²/Ne
    @test dr.F[2] == 0                        # offspring of unrelated founders
    @test dr.F[end] > dr.F[3] > 0
    @test dr.Fhet[end] > 0
    # Without selection the loss of heterozygosity has the pedigree mean
    # coancestry as its expectation: paired over 12 replicates of 6
    # generations on a 10-chromosome genome, within 4 s.e.
    sp10 = GenericSpecies("ten", fill(10_000_000, 10))
    r10 = Xoshiro(7)
    vm10 = VariantMap(sort(unique(rand(r10, 1:100_000_000, 3000))), sp10)
    h10 = Haplotype(BitMatrix(rand(r10, Bool, length(vm10), 300)))
    q10 = eqtl(h10, 1:5:length(vm10), t; rng = r10)
    mm10 = MeiosisMap(sp10, vm10)
    dif = map(1:12) do r
        p = founder_population(h10, q10, t, SimKey(key; replicate = r))
        df = simulate!(p, Scheme(; noff = 150, nsire = 8, ndam = 30, evaluation = :tbv, selection = :random), mm10, 6)
        df.Fhet[end] - df.coancestry[end]
    end
    @test abs(mean(dif)) < 4 * std(dif) / sqrt(12)
end
