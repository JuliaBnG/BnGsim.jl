# The id/row boundary (PedigreeView) and evaluation on a window of
# cohorts. Property tests on random pedigrees (parents born before
# offspring, unknown parents, pruned, sparse and reordered subsets); the
# window's PBLUP against an independent dense solve.

@testset "PedigreeView" begin
    rng = Xoshiro(20261118)
    for rep = 1:50
        n = rand(rng, 2:120)
        sire = [i == 1 || rand(rng) < 0.2 ? 0 : rand(rng, 1:i-1) for i = 1:n]
        dam = [i == 1 || rand(rng) < 0.2 ? 0 : rand(rng, 1:i-1) for i = 1:n]
        ids = rand(rng, 1:n, rand(rng, 1:n))              # sparse, unsorted, repeats
        v = PedigreeView(sire, dam, ids)
        @test v.ids == sort(unique(ids))
        @test length(v) == length(v.ids)
        # round trips, both ways
        @test row_to_id(v, id_to_row(v, ids)) == ids
        @test id_to_row(v, row_to_id(v, 1:length(v))) == 1:length(v)
        # parents as rows: the same animal, or 0 when unknown or outside
        for r = 1:length(v), (p, q) in ((sire, v.sire), (dam, v.dam))
            i = v.ids[r]
            @test q[r] == (p[i] > 0 && p[i] in v.ids ? id_to_row(v, p[i]) : 0)
            @test q[r] < r                                 # topological
        end
        out = setdiff(1:n, v.ids)
        isempty(out) || @test_throws ArgumentError id_to_row(v, first(out))
        @test_throws ArgumentError id_to_row(v, n + 1)
        # a view with every ancestor traced keeps the relationships
        w = pedigree_view(sire, dam, ids)
        @test issubset(v.ids, w.ids)
        r = id_to_row(w, v.ids)
        @test nrm(w.sire, w.dam)[r, r] ≈ nrm(sire, dam)[v.ids, v.ids]
        # depth 0 is the ids alone; depth 1 adds their parents
        @test pedigree_view(sire, dam, ids; generations = 0).ids == v.ids
        @test pedigree_view(sire, dam, ids; generations = 1).ids ==
              sort(unique([v.ids; filter(>(0), [sire[v.ids]; dam[v.ids]])]))
    end
    @test_throws ArgumentError PedigreeView([0, 2], [0, 0], 1:2)       # parent not born before
    @test_throws ArgumentError PedigreeView([0, 1], [0, 0], [3])
    @test_throws ArgumentError pedigree_view([0, 1], [0, 0], [1]; generations = -1)
end

@testset "Evaluation on a window of cohorts" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(7)
    pos = unique(sort(rand(rng, 1:30_000_000, 1000)))
    vm = VariantMap(pos, sp)
    nlc = length(vm)
    hap = Haplotype(BitMatrix(rand(rng, Bool, nlc, 200)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 1.0, μ = 10.0)
    qtl = eqtl(hap, 1:5:nlc, t; rng = Xoshiro(2))
    sc = Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup)
    pop = founder_population(hap, qtl, t, SimKey(11))
    simulate!(pop, sc, mm, 4)                              # cohorts 0–4
    n = length(pop)

    @test_throws ArgumentError Scheme(; noff = 100, nsire = 5, ndam = 20, window = -1)
    @test_throws ArgumentError Scheme(; noff = 100, nsire = 5, ndam = 20, ancestor_generations = -1)
    @test_throws ArgumentError Scheme(; noff = 100, nsire = 5, ndam = 20, window = 2, max_dam_age = 2)

    # window 0 is the identity view; a window over every cohort is the same evaluation
    @test evaluation_view(pop, sc).ids == 1:n
    evaluate!(pop, sc)
    full = copy(pop.ebv)
    evaluate!(pop, Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup, window = 5))
    @test pop.ebv ≈ full rtol = 1e-10

    # the last 2 cohorts with every ancestor: records of cohorts 3–4 only,
    # relationships of the whole pedigree
    sw = Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup, window = 2)
    v = evaluation_view(pop, sw)
    win = collect(generation_ids(pop, 3).start:n)
    @test issubset(win, v.ids) && length(v) < n
    m = evaluate!(pop, sw)
    @test all(isnan, pop.ebv[setdiff(1:n, v.ids), :])
    @test !any(isnan, pop.ebv[v.ids, :])
    A = nrm(pop.sire, pop.dam)[v.ids, v.ids]
    rec = id_to_row(v, win)
    nv = length(v)
    Z = zeros(length(rec), nv)
    for (k, r) in enumerate(rec)
        Z[k, r] = 1
    end
    X = ones(length(rec), 1)
    yw = pop.phenotype[win, 1]
    λ = pop.R0[1, 1] / pop.G0[1, 1]
    C = [X'X X'Z; Z'X Z'Z+λ*inv(A)]
    sol = C \ [X'yw; Z'yw]
    @test pop.ebv[v.ids, 1] ≈ sol[2:end] rtol = 1e-8
    @test size(m.lhs, 1) == nv + 1                        # the equations are in view rows
    # windowing changes the EBV: earlier records carried information
    @test !(pop.ebv[win, 1] ≈ full[win, 1])

    # no ancestors: the window's animals alone, parents outside unknown
    s0 = Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup, window = 2, ancestor_generations = 0)
    v0 = evaluation_view(pop, s0)
    @test v0.ids == win
    @test all(==(0), v0.sire[1:100]) && all(>(0), v0.sire[101:200])
    evaluate!(pop, s0)
    @test all(isnan, pop.ebv[1:first(win)-1, :]) && !any(isnan, pop.ebv[win, :])

    # the loop runs on a window: ROCS (prediction error covariance from the
    # windowed equations), ssGBLUP, and overlapping generations with REML
    for s in (
        Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup, selection = :rocs, ub = 0.1,
            window = 2, ancestor_generations = 2),
        Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :ssgblup, chip = 1:2:nlc,
            genotyped_cohorts = 1, window = 2),
        Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup, max_dam_age = 1, reml = true,
            window = 3, ancestor_generations = 1),
    )
        p = founder_population(hap, qtl, t, SimKey(12))
        d = simulate!(p, s, mm, 4)
        @test nrow(d) == 5 && length(p) == 500
        @test all(isfinite, d.accuracy[1:end-1])
    end
end

@testset "Pedigree table with explicit ids" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000])
    rng = Xoshiro(9)
    vm = VariantMap(unique(sort(rand(rng, 1:30_000_000, 500))), sp)
    nlc = length(vm)
    hap = Haplotype(BitMatrix(rand(rng, Bool, nlc, 100)))
    traits = [Trait("milk"; h² = 0.3, sex = 0, age = 1), Trait("growth"; h² = 0.4)]
    qtl = eqtl(hap, 1:5:nlc, [1.0 0.2; 0.2 1.0], traits; rng = Xoshiro(3))
    pop = founder_population(hap, qtl, traits, SimKey(5))
    simulate!(pop, Scheme(; noff = 60, nsire = 3, ndam = 10, evaluation = :pblup, max_dam_age = 1,
        survival = [0.9], window = 2), MeiosisMap(sp, vm), 3)
    n = length(pop)
    df = pedigree(pop)
    @test names(df) == ["id", "sire", "dam", "male", "generation", "exit", "herd",
        "tbv_milk", "phenotype_milk", "ebv_milk", "tbv_growth", "phenotype_growth", "ebv_growth"]
    @test df.id == 1:n && df.sire == pop.sire && df.dam == pop.dam && df.male == pop.male
    @test df.tbv_growth == pop.tbv[:, 2]
    @test all(ismissing, df.phenotype_milk[df.male]) && !any(ismissing, df.phenotype_milk[.!df.male])
    @test isequal(df.exit, [x == typemax(Int32) ? missing : x for x in pop.exit])
    @test any(ismissing, df.exit) && !all(ismissing, df.exit)
    # outside the last evaluation's window: no EBV
    v = evaluation_view(pop, Scheme(; noff = 60, nsire = 3, ndam = 10, window = 2, max_dam_age = 1), current_step(pop) - 1)
    @test all(ismissing, df.ebv_milk[setdiff(1:n, v.ids)])
    # a reordered subset keeps the ids, and rebuilds the same relationships
    ids = shuffle(Xoshiro(1), collect(n-59:n))[1:30]
    sub = pedigree(pop, ids)
    @test sub.id == ids && sub.tbv_milk == pop.tbv[ids, 1]
    w = pedigree_view(pop, sub.id)
    r = id_to_row(w, sub.id)
    @test nrm(w.sire, w.dam)[r, r] ≈ nrm(pop.sire, pop.dam)[ids, ids]
    @test_throws ArgumentError pedigree(pop, [n + 1])
end
