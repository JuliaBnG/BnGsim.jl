# Report set (Phase 5): QTL status and realized selection intensity.

@testset "Report: QTL status and selection intensity" begin
    # four QTL, one of them without effect; 4 haplotypes
    #   locus 1 (a > 0): 1 1 1 1 -> favourable fixed
    #   locus 2 (a < 0): 1 1 1 1 -> favourable lost
    #   locus 3 (a > 0): 1 0 1 0 -> segregating
    #   locus 4 (a = 0): 0 0 0 0 -> not counted
    #   locus 5 (a < 0): 0 0 0 0 -> favourable fixed
    hap = Haplotype(BitMatrix([1 1 1 1; 1 1 1 1; 1 0 1 0; 0 0 0 0; 0 0 0 0]))
    q = TraitQTL("y", [1, 2, 3, 4, 5], [0.5, -1.0, 2.0, 0.0, -0.3])
    st = BnGsim._qtl_status(hap, q)
    @test (st.qtl_seg, st.qtl_fixed, st.qtl_lost) == (1, 2, 1) && st.qtl_maf == 0.5 && isnan(st.qtl_alpha)
    # architecture: α is the slope of log a² on log 2p(1 - p); the largest
    # tenth (here 1 of 4) of the genic variance
    p, a = [0.1, 0.2, 0.5, 0.0, 0.3], [2.0, -1.0, 0.5, 9.0, 1.0]       # the monomorphic locus is left out
    ar = BnGsim._qtl_architecture(p, a)
    h = 2 .* p[[1, 2, 3, 5]] .* (1 .- p[[1, 2, 3, 5]])
    x, y = log.(h), log.(a[[1, 2, 3, 5]] .^ 2)
    @test ar.qtl_alpha ≈ cov(x, y) / var(x) && ar.qtl_maf ≈ mean([0.1, 0.2, 0.5, 0.3])
    v = h .* a[[1, 2, 3, 5]] .^ 2
    @test ar.qtl_top10 ≈ maximum(v) / sum(v)
    # intensities: candidates' index, parents counted once per offspring
    index = [1.0, 2.0, 3.0, 4.0, 10.0, 20.0, 30.0]
    males, females = [1, 2, 3, 4], [5, 6, 7]
    plan = [4 7; 4 7; 3 6]
    r = BnGsim._intensities(index, (males, females), plan)
    @test r.intensity_sire ≈ ((4 + 4 + 3) / 3 - 2.5) / std([1.0, 2, 3, 4])
    @test r.intensity_dam ≈ ((30 + 30 + 20) / 3 - 20) / 10
    # in a run: the three counts add up to the QTL, truncation selects above the mean
    sm = run_scenario(read_scenario(MINIMAL))
    @test all(sm.qtl_seg .+ sm.qtl_fixed .+ sm.qtl_lost .== 50)
    @test all(>(0), sm.intensity_sire[1:end-1]) && isnan(sm.intensity_sire[end])
end


@testset "Report: reliability and genomic F" begin
    for (doc, ev) in ((JSON.parse(read(MINIMAL, String)), "pblup"), (scenario_doc(), "gblup"))
        doc["scheme"]["evaluation"] = ev
        doc["scheme"]["report_reliability"] = true
        s = parse_scenario(doc)
        rep = run_replicate(s)
        pop, sm = rep.population, rep.summary
        # the PEV blocks, solved 64 candidates at a time, against the dense
        # inverse of the coefficient matrix (one trait and two)
        sc = Scheme(; s.scheme..., chip = rep.chip)
        m = evaluate!(pop, sc)
        n, k = length(evaluation_view(pop, sc)), ntraits(pop)
        @test n > 64
        P = BnGsim._pev_blocks(m, collect(1:n), n, k)
        C = inv(Matrix(m.lhs)) .* m.scale
        a = first(m.blocks[2]) - 1
        @test all(i -> (r = [a + (j - 1) * n + i for j = 1:k]; P[:, :, i] ≈ C[r, r]), 1:n)
        @test all(x -> 0 < x < 1, sm.reliability[1:end-1]) && isnan(sm.reliability[end])
        k > 1 && @test all(x -> 0 < x < 1, sm.reliability_a[1:end-1])
        # Fgrm: VanRaden's G of the last cohort with the founders' frequencies
        h0 = pop.blocks[1]
        p0 = vec(sum(dosage(h0, 1:h0.nlc); dims = 2)) ./ h0.nhp
        seg = findall(x -> 0 < x < 1, p0)
        G = BnGsim.grm(pop.blocks[end]; p = p0, loci = seg)
        @test sm.Fgrm[end] ≈ mean(diag(G)) - 1
    end
    # off by default: NaN
    @test all(isnan, run_scenario(read_scenario(MINIMAL)).reliability)
end
