# The breeding loop beyond discrete generations: parents of several ages,
# sex- and age-limited records, re-estimated variance components, and
# checkpoints.
#
# Scenario as in loop.jl: two chromosomes (20 Mb, 1 M; 10 Mb, 0.5 M),
# ~2000 loci, 200 founders with random haplotypes, one additive trait with
# h² = 0.3 and σₐ = 1, QTL every fifth locus. Seeds are fixed.
using Breeding: reml, RandomEffect, incidence
using RelationshipMatrices: ainv

@testset "Breeding loop: overlapping generations, records, REML, checkpoints" begin
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
    run(sc, n; tr = t) = (p = founder_population(hap, qtl, tr, key); (p, simulate!(p, sc, mm, n)))

    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, max_sire_age = -1)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, evaluation = :tbv, reml = true)
    @test_throws ArgumentError founder_population(hap, qtl, Trait("y"; age = 0.5), key)

    # ---- discrete generations are the special case of age 0 ---------------
    sc = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup)
    p0, d0 = run(sc, 3)
    @test all(d0.L[1:end-1] .== 1)
    @test all(d0.va .== 1) && all(d0.ve .≈ 0.7 / 0.3)
    @test current_step(p0) == 3

    # ---- overlapping generations ------------------------------------------
    # random selection among males of ages 0–2 and females of ages 0–1
    so = Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :tbv, selection = :random,
        max_sire_age = 2, max_dam_age = 1)
    po, dfo = run(so, 5)
    sire_age = Int[]
    dam_age = Int[]
    for g = 1:5
        ids = generation_ids(po, g)
        append!(sire_age, g - 1 .- po.generation[po.sire[ids]])
        append!(dam_age, g - 1 .- po.generation[po.dam[ids]])
        @test all(po.male[po.sire[ids]]) && !any(po.male[po.dam[ids]])
        @test length(unique(po.sire[ids])) == 10 && length(unique(po.dam[ids])) == 40
        # the generation interval reported is the mean parental age at birth
        @test dfo.L[g] ≈ mean(g .- po.generation[[po.sire[ids]; po.dam[ids]]])
    end
    @test extrema(sire_age) == (0, 2) && extrema(dam_age) == (0, 1)
    @test dfo.L[1] == 1                       # only founders are available at step 0
    @test all(dfo.L[3:end-1] .> 1)
    # truncation ranks all eligible ages on one criterion, ties by id
    st = Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :tbv, max_sire_age = 2, max_dam_age = 1)
    pt, _ = run(st, 3)
    evaluate!(pt, st)
    s, ns, d, nd = select_parents(pt, st, nothing)
    elig_m = filter(i -> pt.male[i] && pt.generation[i] >= 1, 1:length(pt))
    elig_f = filter(i -> !pt.male[i] && pt.generation[i] >= 2, 1:length(pt))
    @test s == sort(elig_m[sortperm(pt.tbv[elig_m]; rev = true)[1:10]])
    @test d == sort(elig_f[sortperm(pt.tbv[elig_f]; rev = true)[1:40]])
    @test sum(ns) == sum(nd) == 200
    # ROCS and OCS accept candidates of several ages
    for sel in (:ocs, :rocs)
        sx = Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :pblup, selection = sel, K = 0.3,
            ub = 0.1, max_sire_age = 1, max_dam_age = 1)
        px, dx = run(sx, 3)
        @test nrow(dx) == 4 && all(dx.nsire[1:3] .> 0)
    end

    # ---- sex-limited records ----------------------------------------------
    tf = Trait("y"; h² = 0.3, σₐ = 1.0, μ = 10.0, sex = 0)
    pf, _ = run(Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup), 2; tr = tf)
    @test all(isnan, pf.phenotype[pf.male]) && !any(isnan, pf.phenotype[.!pf.male])
    @test recorded(pf) == findall(.!pf.male)
    # the phenotypes drawn do not depend on which sex is recorded
    pb, _ = run(Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup), 0)
    pf0, _ = run(Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup), 0; tr = tf)
    @test pf0.phenotype[.!pf0.male] == pb.phenotype[.!pb.male]
    # PBLUP with female records against an independent dense solve
    sc = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :pblup)
    evaluate!(pf, sc)
    n, rec = length(pf), recorded(pf)
    Z = Matrix(incidence(rec, n))
    X = ones(length(rec), 1)
    λ = pf.σe^2
    C = [X'X X'Z; Z'X Z'Z+λ*inv(nrm(pf.sire, pf.dam))]
    y = pf.phenotype[rec]
    @test pf.ebv ≈ (C\[X'y; Z'y])[2:end] rtol = 1e-8
    @test_throws ArgumentError evaluate!(pf, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :phenotype))

    # ---- records available from an age ------------------------------------
    ta = Trait("y"; h² = 0.3, σₐ = 1.0, μ = 10.0, age = 1)
    pa = founder_population(hap, qtl, ta, key)
    @test isempty(recorded(pa))
    @test_throws ArgumentError evaluate!(pa, sc)            # nothing to evaluate at step 0
    step!(pa, Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv), mm)
    @test recorded(pa) == collect(generation_ids(pa, 0))
    r = step!(pa, sc, mm)                                    # cohort 1 is evaluated on its parents
    @test isfinite(r.accuracy) && r.accuracy > 0
    @test recorded(pa, 3) == 1:last(generation_ids(pa, 2))

    # ---- REML ---------------------------------------------------------------
    sr = Scheme(; noff = 300, nsire = 15, ndam = 60, evaluation = :pblup, reml = true)
    pr, dr = run(sr, 3)
    @test all(dr.va .> 0) && all(dr.ve .> 0)
    @test dr.va[1] == 1                       # founders only: not separable, true values kept
    @test dr.va[2] != 1 && dr.va[3] != dr.va[2]   # estimated from step 1 on
    # the estimates used at the last step are those of an independent REML
    pr2 = founder_population(hap, qtl, t, key)
    simulate!(pr2, sr, mm, 2)
    vc2 = copy(pr2.vc)
    rec = recorded(pr2)
    rr = reml(pr2.phenotype[rec], ones(length(rec), 1),
        [RandomEffect(incidence(rec, length(pr2)), ainv(pr2.sire, pr2.dam), vc2[1])]; σ²e = vc2[2], tol = 1e-6)
    evaluate!(pr2, sr)
    @test pr2.vc ≈ [rr.G0[1][1, 1], rr.σ²e] rtol = 1e-5
    @test dr.va[3] == pr2.vc[1]

    # ---- checkpoints ----------------------------------------------------------
    sc = Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :pblup, max_sire_age = 1, reml = true)
    full, dfull = run(sc, 4)
    path = joinpath(mktempdir(), "pop.ck")
    half = founder_population(hap, qtl, t, key)
    d1 = simulate!(half, sc, mm, 2; checkpoint = path)
    @test isfile(path) && !isfile(path * ".tmp")
    back = load_checkpoint(path)
    d2 = simulate!(back, sc, mm, 2)
    for f in fieldnames(Population)
        f in (:blocks, :qtl, :trait, :key) && continue
        @test isequal(getfield(back, f), getfield(full, f))
    end
    @test all(packed(back.blocks[k]) == packed(full.blocks[k]) for k in eachindex(full.blocks))
    @test isequal(vcat(d1[1:2, :], d2), dfull)
    # a checkpoint from another version is refused unless asked
    bad = joinpath(dirname(path), "old.ck")
    open(bad, "w") do io
        write(io, read(path, 8))
        BnGsim.serialize(io, (bngsim = "0.0.1", julia = string(VERSION)))
        BnGsim.serialize(io, back)
    end
    @test_throws ArgumentError load_checkpoint(bad)
    @test length(load_checkpoint(bad; strict = false)) == length(back)
    @test_throws ArgumentError load_checkpoint(@__FILE__)
end
