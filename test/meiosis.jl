# Keyed meiosis: reproducibility, and its statistics against exact theory.
#
# Model under test (src/meiosis.jl): independent assortment of
# chromosomes, Poisson crossovers with mean L_c (Morgans) placed uniformly
# on each chromosome, no interference. Statistical checks use one parent
# and 2 × 6000 gametes drawn from keyed streams with seed 20261116; a
# check passes within 4 standard errors (normal approximation), so with
# fixed seeds the outcome is deterministic and a failure is a real change.

@testset "Keyed meiosis" begin
    sp = GenericSpecies("toy", [2_000_000, 1_000_000, 500_000]; morgan = [1.5, 0.5, 0.25])
    nlc = 1200
    rng = Xoshiro(20261116)
    pos = sort(sample(rng, 1:3_500_000, nlc; replace = false))
    vm = VariantMap(pos, sp)
    mm = MeiosisMap(sp, vm)
    @test length(mm) == nlc
    @test mm.morgan == [1.5, 0.5, 0.25]
    @test sum(length, mm.loci) == nlc
    @test occursin("2.25 M", sprint(show, mm))
    @test_throws ArgumentError MeiosisMap(sp, VariantMap(fill(Int8(1), nlc), pos, vm.ref, vm.alt))

    key = SimKey(20261116; generation = 1, stage = 1)

    # ---- reproducibility -------------------------------------------------
    nid = 40
    par = Haplotype(BitMatrix(rand(rng, Bool, nlc, 2nid)))
    plan = hcat(rand(rng, 1:nid, 100), rand(rng, 1:nid, 100))
    off1 = make_progeny!(Haplotype(nlc, 200), par, plan, mm, key)
    off2 = make_progeny!(Haplotype(nlc, 200), par, plan, mm, key)
    @test off1 == off2
    # rows are independent: any subset of the plan gives the same offspring
    # when keyed alike, so the result cannot depend on how rows are split
    # among tasks or threads
    sub = [17, 3, 88]
    for i in sub
        one = Haplotype(nlc, 2 * i)
        p = zeros(Int, i, 2)
        p[i, :] = plan[i, :]
        p[1:i-1, :] .= 1
        make_progeny!(one, par, p, mm, key)
        @test one[:, (2i-1):(2i)] == off1[:, (2i-1):(2i)]
    end
    # every gamete is a mosaic of its parent's two haplotypes
    for i = 1:100, k = 1:2
        g, p = 2i - 2 + k, plan[i, k]
        @test all((off1[l, g] == par[l, 2p-1]) | (off1[l, g] == par[l, 2p]) for l = 1:nlc)
    end
    # a different stage or generation gives different offspring
    @test make_progeny!(Haplotype(nlc, 200), par, plan, mm, SimKey(key; stage = 2)) != off1
    # golden checksum of the offspring of this key (fixes the stream layout,
    # the crossover model and the draw order together)
    checksum(h) = foldl((a, w) -> a * 0x100000001b3 ⊻ w, packed(h).chunks; init = 0xcbf29ce484222325)
    @test checksum(off1) == GOLDEN_MEIOSIS
    # records
    off3, rec = make_progeny!(Haplotype(nlc, 200), par, plan, mm, key; record = true)
    @test off3 == off1
    @test rec.parent[1:2] == plan[1, :] && length(rec.parent) == 200
    @test size(rec.phase) == (3, 200)
    @test sum(ncrossovers(rec, g) for g = 1:200) == length(rec.pos)
    # the record reproduces the gamete: phase and crossovers determine it
    for g in (1, 2, 77, 200)
        p = rec.parent[g]
        xo = collect(crossovers(rec, g))
        for (c, r) in enumerate(mm.loci), l in r
            ph = rec.phase[c, g] ⊻ isodd(count(x -> mm.first[c] <= x <= mm.last[c] && x < pos[l], xo))
            @test off3[l, g] == par[l, 2p-1+ph]
        end
    end
    # v0.2 gene_drop! on Haplotypes whose loci are not a multiple of 64 (1200
    # here) used to throw DimensionMismatch
    og = Haplotype(nlc, 200)
    gene_drop!(og, par, plan, sp, vm; rng = Xoshiro(3))
    @test all((og[l, 1] == par[l, 2plan[1, 1]-1]) | (og[l, 1] == par[l, 2plan[1, 1]]) for l = 1:nlc)
    @test !any(packed(og)[(nlc+1):end, :])                     # padding untouched
    @test_throws DimensionMismatch make_progeny!(Haplotype(nlc, 10), par, plan, mm, key)
    @test_throws ArgumentError make_progeny!(Haplotype(nlc, 200), par, plan .+ nid, mm, key)
    @test_throws DimensionMismatch make_progeny!(Haplotype(nlc, 200), par, plan[:, 1:1], mm, key)

    # ---- statistics against theory ----------------------------------------
    N = 6000                                   # offspring, 2N gametes from one parent
    marked = Haplotype(nlc, 2)
    for l = 1:nlc
        marked[l, 1] = true                    # haplotype 1 all 1, haplotype 2 all 0
    end
    self = fill(1, N, 2)
    gam, rec = make_progeny!(Haplotype(nlc, 2N), marked, self, mm, SimKey(key; stage = 9); record = true)
    G = Matrix(gam)                            # G[l, g] = 1 if from haplotype 1
    ng = 2N

    # crossovers per chromosome ~ Poisson(L_c): mean and variance
    for c = 1:3
        n = [count(x -> mm.first[c] <= x <= mm.last[c], crossovers(rec, g)) for g = 1:ng]
        L = mm.morgan[c]
        @test abs(mean(n) - L) < 4 * sqrt(L / ng)
        @test abs(var(n) - L) < 4 * sqrt((L + 2L^2) / ng)       # Var of s² for Poisson
    end

    # unbiased segregation at every locus: P(haplotype 1) = ½ (largest |z|
    # over 1200 linked loci)
    f = vec(mean(G; dims = 2))
    z = (f .- 0.5) ./ sqrt(0.25 / ng)
    @test maximum(abs, z) < 5                                   # 1200 loci, linked

    # recombination fraction against Haldane, and ½ between chromosomes
    ℓ(c) = Int(mm.last[c] - mm.first[c] + 1)
    function rfrac(a, b)
        r̂ = count(g -> G[a, g] != G[b, g], 1:ng) / ng
        ca, cb = vm.chr[a], vm.chr[b]
        r = ca == cb ? (1 - exp(-2 * mm.morgan[ca] * (pos[b] - pos[a]) / ℓ(ca))) / 2 : 0.5
        r̂, r
    end
    r1 = mm.loci[1]
    for (a, b) in ((first(r1), first(r1) + 5), (first(r1), first(r1) + 60),
            (first(r1), first(r1) + 200), (first(r1), last(r1)),
            (last(r1), first(mm.loci[2])), (first(r1), last(mm.loci[3])))
        r̂, r = rfrac(a, b)
        @test abs(r̂ - r) < 4 * sqrt(max(r * (1 - r), 1e-4) / ng)
    end

    # variance of a gamete's genetic value, with linkage, against the exact
    # value Σ_jk α_j α_k d_j d_k Cov(I_j, I_k), Cov = (1 − 2 r_jk)/4 on one
    # chromosome and 0 between chromosomes (I_j = 1 if from haplotype 1)
    hp = Haplotype(BitMatrix(rand(rng, Bool, nlc, 2)))
    α = randn(rng, nlc) .* (rand(rng, nlc) .< 0.05)            # 5 % of loci are QTL
    q = findall(!iszero, α)
    d = [Float64(hp[l, 1]) - Float64(hp[l, 2]) for l in q]
    V = 0.0
    for (s, j) in enumerate(q), (t, k) in enumerate(q)
        vm.chr[j] == vm.chr[k] || continue
        c = vm.chr[j]
        r = (1 - exp(-2 * mm.morgan[c] * abs(Int(pos[k]) - Int(pos[j])) / ℓ(c))) / 2
        V += α[j] * α[k] * d[s] * d[t] * (1 - 2r) / 4
    end
    gv = make_progeny!(Haplotype(nlc, 2N), hp, self, mm, SimKey(key; stage = 10))
    val = [sum(α[l] * gv[l, g] for l in q) for g = 1:ng]
    @test abs(var(val) - V) < 4 * V * sqrt(2 / (ng - 1))
end
