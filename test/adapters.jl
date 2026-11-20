# The upstream contract (docs/src/upstream-contract.md): the behaviour of
# FisherWright, RelationshipMatrices and Breeding that BnGsim relies on,
# through its adapters, and the rule that nothing else imports them.

@testset "Upstream contract: layering" begin
    # Outside src/adapters/, no file imports or qualifies the upstream
    # packages (docstrings and comments may name them).
    src = joinpath(@__DIR__, "..", "src")
    pkgs = r"\b(FisherWright|RelationshipMatrices|Breeding)\b"
    for f in readdir(src; join=true)
        isfile(f) && endswith(f, ".jl") || continue
        indoc = false
        for (n, l) in enumerate(readlines(f))
            if count("\"\"\"", l) == 1
                indoc = !indoc
                continue
            end
            (indoc || startswith(strip(l), "#")) && continue
            code = replace(l, r"#.*$" => "", r"\"[^\"]*\"" => "\"\"", r"`[^`]*`" => "``")
            bad = occursin(pkgs, code)
            bad && @info "upstream reference outside src/adapters/" file = basename(f) line = n text = l
            @test !bad
        end
    end
end

@testset "Upstream contract: FisherWright" begin
    chr = [100_000, 50_000]
    sp = GenericSpecies("c", chr; morgan=[0.4, 0.2])
    res = BnGsim.fisher_wright(sp, 30, 40, 2.0; mut_base=1e5, result=true, rng=Xoshiro(1))
    # haplotypes 2i - 1 and 2i of individual i; sorted, unique, genome-wide,
    # 1-based UInt32 positions; cumulative chromosome ends
    @test res isa BnGsim.FisherWrightResult && length(res.active_haplotypes) == 2 * 30
    @test all(h -> eltype(h) == UInt32 && issorted(h) && allunique(h) && all(p -> 1 <= p <= sum(chr), h),
        res.active_haplotypes)
    @test res.chromosome_ends == UInt32.(cumsum(chr)) && issorted(res.substitutions)
    @test BnGsim.census(res) == fill(30, 41)
    # continuation appends to the record; the same rng gives the same run
    more = BnGsim.fisher_wright(30, 3, chr, 2.0; init=res, popsize=[20, 20, 10], result=true, rng=Xoshiro(2))
    @test BnGsim.census(more)[end-3:end] == [30, 20, 20, 10] && length(more.active_haplotypes) == 20
    again = BnGsim.fisher_wright(sp, 30, 40, 2.0; mut_base=1e5, result=true, rng=Xoshiro(1))
    @test again.active_haplotypes == res.active_haplotypes
    # dense export: loci segregating, rows by haplotype columns, map with
    # 1-based chromosome, genome-wide position and frequency
    xy, map = BnGsim.dense_haplotypes(res.active_haplotypes, res.chromosome_ends; rng=Xoshiro(3))
    @test size(xy) == (nrow(map), 60) && all(0 .< vec(sum(xy; dims=2)) .< 60)
    @test map.pos == sort(unique(reduce(vcat, res.active_haplotypes)))
    @test map.chr == Int8.(ifelse.(map.pos .<= chr[1], 1, 2)) && map.frq ≈ vec(sum(xy; dims=2)) ./ 60
    # a copied rng gives the same map; the caller's stream is untouched
    rng = Xoshiro(9)
    a, _ = founder_haplotypes(res, 10; rng=rng)
    @test rand(rng) == rand(begin r = Xoshiro(9); founder_haplotypes(res, 10; rng=r); r end)
end

@testset "Upstream contract: RelationshipMatrices" begin
    # pedigree: 0 for unknown, parents before offspring
    sire = Int32[0, 0, 1, 1, 3]
    dam = Int32[0, 0, 2, 2, 4]
    A = BnGsim.nrm(sire, dam)
    @test A[5, 5] ≈ 1.25 && A[3, 4] ≈ 0.5 && BnGsim.nrm_diag(sire, dam) ≈ diag(A)
    @test Matrix(BnGsim.ainv(sire, dam)) ≈ inv(A) && BnGsim.nrm(sire, dam, [3, 5]) ≈ A[[3, 5], [3, 5]]
    # grm: VanRaden (2008) G = ZZ'/(2Σpq) from the haplotypes, Z = M - 2p
    rng = Xoshiro(4)
    hap = Haplotype(BitMatrix(rand(rng, Bool, 120, 10)))
    M = Float64[hap[l, 2i-1] + hap[l, 2i] for i = 1:5, l = 1:120]
    p = vec(sum(M; dims=1)) ./ 10
    Z = M .- 2p'
    @test BnGsim.grm(hap; p=p) ≈ Z * Z' ./ (2sum(p .* (1 .- p)))
    @test BnGsim.grm(hap; p=p, loci=1:60) ≈ Z[:, 1:60] * Z[:, 1:60]' ./ (2sum(p[1:60] .* (1 .- p[1:60])))
    # hinv: A⁻¹ + [0 0; 0 G⁻¹ - A₂₂⁻¹] for the genotyped animals 3:5
    G = BnGsim.grm(hap; p=p)[3:5, 3:5] + 0.1I
    H = Matrix(BnGsim.hinv(sire, dam, G, [3, 4, 5]))
    E = inv(A)
    E[3:5, 3:5] .+= inv(G) .- inv(A[3:5, 3:5])
    @test H ≈ E
end

@testset "Upstream contract: Breeding" begin
    # y = 1μ + Za + e with a ~ N(0, Aσ²a): the solution of mme is that of
    # Henderson's equations, solved densely.
    sire = Int32[0, 0, 1, 1, 3, 3]
    dam = Int32[0, 0, 2, 2, 4, 4]
    Ainv = BnGsim.ainv(sire, dam)
    y = [1.2, -0.3, 0.8, 2.1]
    rec = [3, 4, 5, 6]
    X = ones(4, 1)
    Z = Matrix(BnGsim.incidence(rec, 6))
    σa, σe = 0.5, 1.5
    m = BnGsim.mme(y, X, [BnGsim.RandomEffect(BnGsim.incidence(rec, 6), Ainv, σa)]; σ²e=σe)
    x = BnGsim.solve_mme(m)
    C = [X'X X'Z; Z'X Z'Z+Matrix(Ainv)*σe/σa]
    @test x ≈ C \ [X'y; Z'y]          # [μ; a], one fixed effect then the animals
    # TM1997: contributions and feasibility of a coancestry bound
    A = BnGsim.nrm(sire, dam)
    r = BnGsim.TM1997([0.1, 0.2, 0.5, 0.4, 0.9, 0.3], Bool[1, 0, 1, 0, 1, 0], A, 1.0)
    @test BnGsim.feasible(r) && sum(r.c) ≈ 1
end
