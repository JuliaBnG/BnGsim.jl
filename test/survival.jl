# Survival, culling and single-step GBLUP in the breeding loop.
#
# Scenario as in loop.jl: two chromosomes, ~2000 loci, 200 founders with
# random haplotypes, one additive trait (h² = 0.3, σₐ = 1), QTL every
# fifth locus. Statistical checks use 4 standard errors.

@testset "Breeding loop: survival, culling, single-step" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(20261116)
    pos = unique(sort(rand(rng, 1:30_000_000, 2000)))
    vm = VariantMap(pos, sp)
    nlc = length(vm)
    hap = Haplotype(BitMatrix(rand(rng, Bool, nlc, 400)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 1.0, μ = 10.0)
    qtl = eqtl(hap, 1:5:nlc, t; rng = Xoshiro(2))
    key = SimKey(99)
    run(sc, n; tr = t) = (p = founder_population(hap, qtl, tr, key); (p, simulate!(p, sc, mm, n)))

    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, survival = [1.2])
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, evaluation = :ssgblup)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, genotyped_cohorts = -1)

    # ---- survival ----------------------------------------------------------------
    q = [0.8, 0.5]
    sc = Scheme(; noff = 400, nsire = 20, ndam = 60, evaluation = :tbv, selection = :random,
        max_sire_age = 3, max_dam_age = 3, survival = q)
    pop, df = run(sc, 6)
    @test all(pop.exit[generation_ids(pop, 6)] .== typemax(Int32))
    # transitions g → g + 1 of animals alive at g, by age class (0, ≥ 1)
    at_risk, survived = zeros(Int, 2), zeros(Int, 2)
    for i = 1:length(pop), g = pop.generation[i]:5
        pop.exit[i] > g || continue
        c = min(g - pop.generation[i], 1) + 1
        at_risk[c] += 1
        survived[c] += pop.exit[i] > g + 1
    end
    for c = 1:2
        f = survived[c] / at_risk[c]
        @test abs(f - q[c]) < 4 * sqrt(q[c] * (1 - q[c]) / at_risk[c])
    end
    # dead animals are never parents
    for g = 1:6
        ids = generation_ids(pop, g)
        @test all(pop.exit[pop.sire[ids]] .> g - 1) && all(pop.exit[pop.dam[ids]] .> g - 1)
    end
    @test df.alive[1] == 200 && df.alive[2] == count(pop.exit[1:200] .> 1) + 400
    @test alive(pop, 6) == findall(i -> pop.generation[i] <= 6 < pop.exit[i], 1:length(pop))
    # a record at age 1 exists only if the animal reached age 1
    ta = Trait("y"; h² = 0.3, σₐ = 1.0, μ = 10.0, age = 1)
    pa = founder_population(hap, qtl, ta, key)
    simulate!(pa, sc, mm, 3)
    obs = observed(pa)[:, 1]
    @test obs == [3 - pa.generation[i] >= 1 && pa.generation[i] + 1 < pa.exit[i] for i = 1:length(pa)]
    @test count(obs) < count(i -> 3 - pa.generation[i] >= 1, 1:length(pa))
    # nobody survives a step: only the newborn cohort is left to select from
    p0, _ = run(Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv, max_sire_age = 2,
        survival = [0.0]), 3)
    @test all(g -> all(p0.generation[p0.sire[generation_ids(p0, g)]] .== g - 1), 1:3)
    # keyed: resuming from a checkpoint repeats the deaths
    path = joinpath(mktempdir(), "pop.ck")
    full, dfull = run(sc, 4)
    half = founder_population(hap, qtl, t, key)
    d1 = simulate!(half, sc, mm, 2; checkpoint = path)
    back = load_checkpoint(path)
    d2 = simulate!(back, sc, mm, 2)
    @test back.exit == full.exit && isequal(vcat(d1[1:2, :], d2), dfull)

    # ---- culling of unselected adults --------------------------------------------
    scc = Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :tbv, max_sire_age = 2, max_dam_age = 2,
        cull_unselected = true)
    pc, _ = run(scc, 5)
    for g = 0:4
        parents = Set([pc.sire[generation_ids(pc, g + 1)]; pc.dam[generation_ids(pc, g + 1)]])
        for i = 1:length(pc)
            a = g - pc.generation[i]
            (a >= 1 && pc.exit[i] > g && a <= 2) || continue
            # alive adults of eligible age at g: kept if and only if selected
            @test (pc.exit[i] == g + 1) == !(i in parents)
        end
    end
    # culled animals are not selected later
    for g = 1:4, i in [pc.sire[generation_ids(pc, g + 1)]; pc.dam[generation_ids(pc, g + 1)]]
        @test pc.exit[i] > g
    end

    # ---- single-step GBLUP -------------------------------------------------------
    pb, _ = run(Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup), 2)
    chip = collect(2:2:nlc)
    sg = Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :gblup, chip = chip)
    ss = Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :ssgblup, chip = chip)
    evaluate!(pb, sg)
    e_g = copy(pb.ebv)
    evaluate!(pb, ss)
    @test genotyped(pb, ss) == 1:length(pb)
    @test pb.ebv ≈ e_g rtol = 1e-6               # everyone genotyped: H⁻¹ = G*⁻¹
    # genotypes of the last cohort only
    s1 = Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :ssgblup, chip = chip, genotyped_cohorts = 1)
    @test genotyped(pb, s1) == collect(generation_ids(pb, 2))
    evaluate!(pb, s1)
    e_s = copy(pb.ebv)
    evaluate!(pb, Scheme(; noff = 100, nsire = 5, ndam = 20, evaluation = :pblup))
    @test all(isfinite, e_s) && e_s != pb.ebv
    young = generation_ids(pb, 2)
    @test cor(e_s[young], pb.tbv[young]) > 0.3
    # ungenotyped animals change through their genotyped relatives only
    @test e_s[1:5] != pb.ebv[1:5]
    # the loop runs with single-step evaluations, REML and two traits
    _, dss = run(Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :ssgblup, chip = chip,
        genotyped_cohorts = 2, reml = true), 3)
    @test nrow(dss) == 4 && all(dss.va .> 0)
    ts = [Trait("a"; h² = 0.3), Trait("b"; h² = 0.2, sex = 0)]
    mq = eqtl(hap, 1:5:nlc, [1 0.5; 0.5 1], ts; rng = Xoshiro(3))
    pm = founder_population(hap, mq, ts, key)
    dm = simulate!(pm, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :ssgblup, chip = chip,
        genotyped_cohorts = 1, weights = [1.0, 1.0], survival = [0.9]), mm, 2)
    @test nrow(dm) == 3 && all(isfinite, dm.accuracy_b[1:2])
end
