# The breeding loop with several traits: founders, multi-trait BLUP
# against a dense solve, the index, the prediction error covariance of the
# index, traits without records yet, REML and checkpoints.
#
# Scenario as in loop.jl (two chromosomes, ~2000 loci, 200 founders with
# random haplotypes); two traits with σₐ = 1, genetic correlation 0.5 in
# the founders (exact, from `eqtl`), h² 0.3 and 0.2, residual correlation
# 0.2, QTL every fifth locus shared by both traits.
using Breeding: pev_candidates, mme, RandomEffect, incidence
using RelationshipMatrices: ainv

@testset "Breeding loop: several traits" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(20261116)
    pos = unique(sort(rand(rng, 1:30_000_000, 2000)))
    vm = VariantMap(pos, sp)
    nlc = length(vm)
    hap = Haplotype(BitMatrix(rand(rng, Bool, nlc, 400)))
    mm = MeiosisMap(sp, vm)
    Rg = [1 0.5; 0.5 1]
    Re = [1 0.2; 0.2 1]
    ts = [Trait("a"; h² = 0.3, μ = 10.0), Trait("b"; h² = 0.2, μ = 5.0)]
    mq = eqtl(hap, 1:5:nlc, Rg, ts; rng = Xoshiro(3))
    key = SimKey(7)
    found(tr = ts) = founder_population(hap, mq, tr, key; Re = Re)

    # ---- founders ------------------------------------------------------------
    pop = found()
    @test ntraits(pop) == 2 && size(pop.tbv) == (200, 2)
    @test pop.tbv ≈ tbv(hap, mq)
    @test cor(pop.tbv) ≈ Rg atol = 1e-8
    @test pop.G0 ≈ Rg atol = 1e-8
    σe = sqrt.([0.7 / 0.3, 0.8 / 0.2])
    @test pop.R0 ≈ σe .* Re .* σe'
    @test pop.Lr * pop.Lr' ≈ pop.R0
    # the first trait's phenotypes are those of a single-trait run
    p1 = founder_population(hap, pop.qtl[1], ts[1], key)
    @test p1.phenotype[:, 1] == pop.phenotype[:, 1]
    @test_throws ArgumentError founder_population(hap, mq, reverse(ts), key)
    @test_throws ArgumentError founder_population(hap, mq, ts, key; Re = [1 0.2; 0.3 1])
    @test_throws ArgumentError founder_population(hap, mq, ts, key; Re = [2 0; 0 1])
    @test_throws ArgumentError step!(pop, Scheme(; noff = 200, nsire = 10, ndam = 50, weights = [1.0]), mm)

    # residuals are correlated as Re (4 s.e. of a correlation)
    pr = found()
    simulate!(pr, Scheme(; noff = 400, nsire = 20, ndam = 80, evaluation = :tbv, selection = :random), mm, 4)
    E = pr.phenotype .- [10.0 5.0] .- pr.tbv
    @test abs(cor(E[:, 1], E[:, 2]) - 0.2) < 4 * (1 - 0.2^2) / sqrt(length(pr))
    @test std(E[:, 1]) ≈ σe[1] rtol = 4 / sqrt(2length(pr))

    # ---- multi-trait PBLUP against a dense solve ---------------------------
    # records missing at random by sex: trait b only in females
    tf = [ts[1], Trait("b"; h² = 0.2, μ = 5.0, sex = 0)]
    pm = found(tf)
    sc = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup, weights = [1.0, 2.0])
    simulate!(pm, sc, mm, 2)
    @test all(isnan, pm.phenotype[pm.male, 2]) && !any(isnan, pm.phenotype[.!pm.male, 2])
    m = evaluate!(pm, sc)
    n = length(pm)
    obs = observed(pm)
    Ainv = Matrix(ainv(pm.sire, pm.dam))
    # stacked trait-major, records of one animal correlated through R0
    y = [pm.phenotype[obs[:, 1], 1]; pm.phenotype[obs[:, 2], 2]]
    Z1 = Matrix(incidence(findall(obs[:, 1]), n))
    Z2 = Matrix(incidence(findall(obs[:, 2]), n))
    n1, n2 = size(Z1, 1), size(Z2, 1)
    X = [ones(n1) zeros(n1); zeros(n2) ones(n2)]
    Z = [Z1 zeros(n1, n); zeros(n2, n) Z2]
    ids = [findall(obs[:, 1]); findall(obs[:, 2])]
    tr = [fill(1, n1); fill(2, n2)]
    R = [ids[i] == ids[j] ? pm.R0[tr[i], tr[j]] : 0.0 for i in eachindex(ids), j in eachindex(ids)]
    Ri = inv(R)
    C = [X'Ri*X X'Ri*Z; Z'Ri*X Z'Ri*Z+kron(inv(pm.G0), Ainv)]
    sol = C \ [X'Ri*y; Z'Ri*y]
    @test vec(pm.ebv) ≈ sol[3:end] rtol = 1e-7
    # selection is on the index w'â, ties by id
    s, = select_parents(pm, sc, m)
    cand = collect(generation_ids(pm, 2))
    males = cand[pm.male[cand]]
    I2 = pm.ebv * [1.0, 2.0]
    @test s == sort(males[sortperm(I2[males]; rev = true)[1:10]])

    # ---- prediction error covariance of the index ---------------------------
    c = cand[1:30]
    Γ1 = pev_candidates(m, 1, c; component = 1)
    Γ2 = pev_candidates(m, 1, c; component = 2)
    @test BnGsim._pev_index(m, c, [1.0, 0.0], n) ≈ Γ1
    @test BnGsim._pev_index(m, c, [0.0, 1.0], n) ≈ Γ2
    Ω = BnGsim._pev_index(m, c, [1.0, 2.0], n)
    @test issymmetric(Ω) && isposdef(Ω)
    # Var(w'(a - â)) of one candidate, from the dense inverse
    Ci = inv(C)
    i = c[5]
    Pi = Ci[[2 + i, 2 + n + i], [2 + i, 2 + n + i]]
    @test Ω[5, 5] ≈ [1.0, 2.0]' * Pi * [1.0, 2.0] rtol = 1e-7
    # rocs and OCS run on the index
    for sel in (:rocs, :ocs)
        px = found()
        dx = simulate!(px, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup, selection = sel,
            K = 0.3, ub = 0.1, weights = [1.0, 2.0]), mm, 2)
        @test nrow(dx) == 3 && all(dx.nsire[1:2] .> 0)
    end

    # ---- a trait without records yet is predicted through G0 --------------
    ta = [ts[1], Trait("b"; h² = 0.2, μ = 5.0, age = 1)]
    pa = found(ta)
    @test isempty(recorded(pa; trait = 2)) && length(recorded(pa; trait = 1)) == 200
    evaluate!(pa, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup))
    @test pa.ebv[:, 2] ≈ pa.G0[2, 1] / pa.G0[1, 1] .* pa.ebv[:, 1]   # unrelated founders, no records of b
    @test cor(pa.ebv[:, 2], pa.tbv[:, 2]) > 0

    # ---- summaries ---------------------------------------------------------------
    pd, df = (p = found(); (p, simulate!(p, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv,
        weights = [1.0, 2.0]), mm, 3)))
    @test df.mean[1] ≈ mean(pd.tbv[1:200, :] * [1.0, 2.0]) atol = 1e-10
    @test df.mean[2] ≈ mean(pd.tbv[generation_ids(pd, 1), :] * [1.0, 2.0])
    @test df.mean_a[2] ≈ mean(pd.tbv[generation_ids(pd, 1), 1])
    @test df.rg_a_b[1] ≈ 0.5 atol = 1e-8
    @test all(df.accuracy[1:3] .≈ 1) && all(df.accuracy_b[1:3] .≈ 1)
    @test df.mean[4] > df.mean[1]
    @test hasproperty(df, :va_a) && hasproperty(df, :ve_b) && !hasproperty(df, :va)

    # ---- REML ------------------------------------------------------------------
    sr = Scheme(; noff = 300, nsire = 15, ndam = 60, evaluation = :pblup, weights = [1.0, 1.0], reml = true)
    prm = found()
    dr = simulate!(prm, sr, mm, 3)
    @test dr.va_a[1] == pop.G0[1, 1]          # founders only: true values kept
    @test dr.va_a[2] != dr.va_a[1] && isposdef(prm.G0) && isposdef(prm.R0)
    # traits never recorded together: residual covariance not estimable
    tx = [Trait("a"; h² = 0.3, sex = 0), Trait("b"; h² = 0.2, sex = 1)]
    px = founder_population(hap, mq, tx, key)
    step!(px, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv), mm)
    @test_throws ArgumentError evaluate!(px, sr)

    # ---- checkpoints -------------------------------------------------------------
    sc = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup, weights = [1.0, 2.0], reml = true,
        max_dam_age = 1)
    full = found()
    dfull = simulate!(full, sc, mm, 3)
    path = joinpath(mktempdir(), "pop.ck")
    half = found()
    d1 = simulate!(half, sc, mm, 1; checkpoint = path)
    back = load_checkpoint(path)
    d2 = simulate!(back, sc, mm, 2)
    @test isequal(vcat(d1[1:1, :], d2), dfull)
    @test isequal(back.ebv, full.ebv) && back.G0 == full.G0 && back.R0 == full.R0
end
