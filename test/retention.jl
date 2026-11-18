# Genome retention: dropping the genomes no later step can use leaves the
# run unchanged, keeps exactly the genomes still usable, and refuses to
# read a dropped one.

@testset "Genome retention" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(20261118)
    vm = VariantMap(unique(sort(rand(rng, 1:30_000_000, 800))), sp)
    nlc = length(vm)
    hap = Haplotype(BitMatrix(rand(rng, Bool, nlc, 200)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 1.0, μ = 10.0)
    qtl = eqtl(hap, 1:5:nlc, t; rng = Xoshiro(2))
    run(sc, n) = (p = founder_population(hap, qtl, t, SimKey(3)); (p, simulate!(p, sc, mm, n)))

    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, retain = :some)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, evaluation = :gblup, chip = 1:10,
        retain = :active)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, evaluation = :ssgblup, chip = 1:10,
        retain = :active)

    # overlapping generations, deaths, culling, ssGBLUP on the last 2 cohorts
    base = (noff = 100, nsire = 5, ndam = 25, evaluation = :ssgblup, chip = 1:2:nlc, genotyped_cohorts = 2,
        max_sire_age = 3, max_dam_age = 1, survival = [0.8], cull_unselected = true, record_meiosis = true)
    pa, da = run(Scheme(; base...), 6)
    pr, dr = run(Scheme(; base..., retain = :active), 6)
    # the same run
    @test isequal(da, dr)
    @test pa.sire == pr.sire && pa.dam == pr.dam && pa.exit == pr.exit
    @test pa.tbv == pr.tbv && isequal(pa.phenotype, pr.phenotype) && isequal(pa.ebv, pr.ebv)
    # the genomes kept are the usable ones, and are the right ones
    tt = current_step(pr)
    usable(i) = (a = tt - pr.generation[i];
        pr.generation[i] == 0 || pr.exit[i] > tt && a <= (pr.male[i] ? 3 : 1) || a < 2)
    @test retained(pa) == 1:length(pa)
    @test retained(pr) == filter(usable, 1:length(pr))
    @test length(retained(pr)) < length(pr)
    keep = retained(pr)
    @test genomes(pr, keep).gt == genomes(pa, keep).gt
    dropped = setdiff(1:length(pr), keep)
    @test_throws ArgumentError genomes(pr, dropped[1:1])
    @test_throws ArgumentError generation_summary(pr, pr.generation[dropped[1]])
    @test generation_summary(pr, tt).F == generation_summary(pa, tt).F       # current cohort kept

    # discrete generations: the founders and the newest cohort only
    pd, _ = run(Scheme(; noff = 100, nsire = 5, ndam = 25, evaluation = :pblup, retain = :active), 4)
    @test retained(pd) == [1:100; collect(generation_ids(pd, 4))]
    @test count(isnothing, pd.blocks) == 3
end
