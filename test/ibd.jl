# Realized IBD from the loop's crossover records: founder segments agree
# with the haplotypes at every locus, and realized relationships follow
# the pedigree.
using RelationshipMatrices: RelationshipMatrices, nrm

@testset "Realized IBD from founder segments" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(20261116)
    pos = unique(sort(rand(rng, 1:30_000_000, 2000)))
    vm = VariantMap(pos, sp)
    nlc = length(vm)
    hap = Haplotype(BitMatrix(rand(rng, Bool, nlc, 120)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3)
    qtl = eqtl(hap, 1:5:nlc, t; rng = Xoshiro(2))
    key = SimKey(31)

    sc = Scheme(; noff = 60, nsire = 6, ndam = 15, evaluation = :tbv, selection = :random,
        max_sire_age = 1, max_dam_age = 1, record_meiosis = true)
    pop = founder_population(hap, qtl, t, key)
    @test founder_segments(pop, mm).ptr == collect(1:2:241) .|> Int   # founders: one segment per chromosome
    simulate!(pop, sc, mm, 4)
    @test length(pop.meioses) == 4 && all(!isnothing, pop.meioses)
    segs = founder_segments(pop, mm)
    # every allele is the allele of the founder haplotype its segment names
    H = BnGsim._haplotypes(pop, 1:length(pop))
    F0 = pop.blocks[1]
    ok = true
    for h = 1:2length(pop)
        st, spp, lb = RelationshipMatrices.segments(segs, h)
        k = 1
        for l = 1:nlc
            while mm.pos[l] > spp[k]
                k += 1
            end
            ok &= H.gt[l, h] == F0.gt[l, lb[k]]
        end
    end
    @test ok
    # realized relationships follow the pedigree
    R = irm(segs)
    A = nrm(pop.sire, pop.dam)
    # exact where the genome cannot vary: founders are unrelated and
    # non-inbred, and an offspring of two founders shares exactly one
    # haplotype with each
    @test R[1:60, 1:60] == I
    for i in generation_ids(pop, 1)
        @test R[i, pop.sire[i]] == 0.5 && R[i, pop.dam[i]] == 0.5
    end
    # elsewhere they scatter around A (the variance is checked exactly in
    # RelationshipMatrices and in validation/ibd.jl)
    @test abs(mean(R .- A)) < 0.01
    @test diag(R) .- 1 ≈ segment_inbreeding(segs)
    # unrecorded cohorts: the segments are unknown
    p2 = founder_population(hap, qtl, t, key)
    simulate!(p2, Scheme(; noff = 60, nsire = 6, ndam = 15, evaluation = :tbv, selection = :random), mm, 1)
    @test_throws ArgumentError founder_segments(p2, mm)
    # recording does not change the run
    @test p2.sire == pop.sire[1:120] && p2.tbv == pop.tbv[1:120, :]
end
