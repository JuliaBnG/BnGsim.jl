# OCS with a target rate of inbreeding: the solver against TM1997, the
# age-class weights, the constrained quantity, and the coincidence of the
# three constraints with discrete generations.
using Breeding: TM1997

@testset "OCS with a target rate of inbreeding" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(20261116)
    pos = unique(sort(rand(rng, 1:30_000_000, 2000)))
    vm = VariantMap(pos, sp)
    nlc = length(vm)
    hap = Haplotype(BitMatrix(rand(rng, Bool, nlc, 400)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 1.0)
    qtl = eqtl(hap, 1:5:nlc, t; rng = Xoshiro(2))
    key = SimKey(5)

    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, selection = :ocs, ocs = :cohort)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, selection = :ocs, ocs = :age, ΔF = 0.01)

    # ---- the solver -------------------------------------------------------------
    p = founder_population(hap, qtl, t, key)
    simulate!(p, Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :tbv, selection = :random), mm, 2)
    cand = collect(generation_ids(p, 2))
    A = nrm(p.sire, p.dam, cand)
    u = p.tbv[cand, 1]
    male = Vector{Bool}(p.male[cand])
    # KKT conditions at the solution: with ν the gain weight, the gradient
    # g = 2Θc + 2b − 2νu plus a multiplier per sex vanishes on the positive
    # contributions and is ≥ 0 on the zero ones
    function kkt(c, Θ, b, u)
        P = c .> 1e-12
        G = 2Θ * c .+ 2b
        # least squares for (ν, λm, λf) on the positive set
        X = [-2u[P] Float64.(male[P]) Float64.(.!male[P])]
        θ = X \ -G[P]
        r = G .+ [-2u Float64.(male) Float64.(.!male)] * θ
        maximum(abs, r[P]) / maximum(abs, G), minimum(r[.!P]; init = 0.0) / maximum(abs, G), θ[1]
    end
    # the smallest mean kinship these candidates can give (bound 0 is infeasible)
    _, ok, fmin = ocs_contributions(u, A ./ 2, zeros(length(u)), 0.0, male)
    @test !ok && fmin > 0
    for K in fmin .+ (0.001, 0.005, 0.05)
        c, ok, val = ocs_contributions(u, A ./ 2, zeros(length(u)), K, male)
        @test ok && all(c .>= 0)
        @test sum(c[male]) ≈ 0.5 && sum(c[.!male]) ≈ 0.5
        @test val ≈ K rtol = 1e-8                  # the constraint is active
        st, comp, ν = kkt(c, A ./ 2, zeros(length(u)), u)
        @test st < 1e-8 && comp > -1e-8 && ν > 0
        # TM1997 drops negative contributions without re-admitting them;
        # where its solution has none, it is the same optimum
        ref = TM1997(u, Int.(male), A, 2K)
        all(ref .>= 0) && @test c ≈ ref atol = 1e-6
    end
    # gain rises with the bound
    g(K) = (c = first(ocs_contributions(u, A ./ 2, zeros(length(u)), K, male)); dot(u, c))
    @test g(fmin + 0.001) < g(fmin + 0.005) < g(fmin + 0.05)
    # an unreachable bound gives the minimum-kinship contributions, flagged
    c0, ok0, v0 = ocs_contributions(u, A ./ 2, zeros(length(u)), fmin / 2, male)
    @test !ok0 && all(c0 .>= 0) && v0 ≈ fmin
    st0, comp0, _ = kkt(c0, A ./ 2, zeros(length(u)), zeros(length(u)))
    @test st0 < 1e-8 && comp0 > -1e-8          # the minimum satisfies its own KKT conditions
    # the linear term: constraint on (c + z)'Θ(c + z) with fixed z
    z = rand(Xoshiro(9), length(u)) ./ 400
    Θ = A ./ 2
    _, _, fz = ocs_contributions(u, Θ, Θ * z, 0.0, male)
    c, ok, val = ocs_contributions(u, Θ, Θ * z, fz + 0.01, male)
    @test ok
    @test dot(c .+ z, Θ * (c .+ z)) ≈ fz + 0.01 + dot(z, Θ * z) rtol = 1e-8
    st, comp, _ = kkt(c, Θ, Θ * z, u)
    @test st < 1e-8 && comp > -1e-8

    # ---- age-class weights --------------------------------------------------------
    so = Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :tbv, selection = :ocs, ocs = :age_year,
        ΔF = 0.01, max_sire_age = 2, max_dam_age = 1)
    po = founder_population(hap, qtl, t, key)
    d = simulate!(po, so, mm, 4)
    W0, v, L = BnGsim._age_weights(po, so, 4)
    @test W0 + sum(values(v)) ≈ 1                # the population at t + 1 is a mixture
    @test all(i -> 4 - po.generation[i] + 1 <= (po.male[i] ? 2 : 1), keys(v))
    @test L >= 1
    # usage by age: at step 0 only age 0 exists
    @test BnGsim._age_usage(po, 0, true, 2) == [1.0, 0.0, 0.0]
    # the constrained kinship reported is w'Θw of the population at t + 1
    cand = sort!(vcat(BnGsim._eligible(po, so, 4)...))
    c, info = BnGsim._dynamic_ocs(po, so, 4, cand, po.tbv[cand, 1])
    ids = sort!(union(cand, keys(v)))
    w = zeros(length(ids))
    for (k, i) in enumerate(ids)
        j = findfirst(==(i), cand)
        w[k] = (j === nothing ? 0.0 : W0 * c[j]) + get(v, i, 0.0)
    end
    Θ = nrm(po.sire, po.dam, ids) ./ 2
    @test dot(w, Θ * w) ≈ info.ocs_kinship rtol = 1e-8
    @test !info.ocs_feasible || info.ocs_kinship <= info.ocs_bound + 1e-12
    @test all(d.ocs_bound[1:4] .> 0) && issorted(d.ocs_bound[1:4])
    # the per-generation constraint iterates L to the generation interval of c
    sg = Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :tbv, selection = :ocs, ocs = :age_generation,
        ΔF = 0.01, max_sire_age = 2, max_dam_age = 1)
    cg, ig = BnGsim._dynamic_ocs(po, sg, 4, cand, po.tbv[cand, 1])
    @test ig.ocs_L ≈ sum(cg .* (4 .- po.generation[cand] .+ 1)) rtol = 1e-4

    # ---- discrete generations: the three constraints coincide ---------------------
    runs = map((:cohort, :age_year, :age_generation)) do m
        pm = founder_population(hap, qtl, t, key)
        dm = simulate!(pm, Scheme(; noff = 200, nsire = 10, ndam = 40, evaluation = :pblup, selection = :ocs,
            ocs = m, ΔF = 0.02), mm, 3)
        (pm.sire, pm.dam, dm.ocs_bound)
    end
    @test runs[1][1] == runs[2][1] == runs[3][1] && runs[1][2] == runs[2][2] == runs[3][2]
    @test runs[1][3][1:3] ≈ runs[2][3][1:3] ≈ runs[3][3][1:3]
    # the bound rises by ΔF per generation from 1/(2n₀)
    @test runs[1][3][1] ≈ 1 - (1 - 1 / 400) * (1 - 0.02)
    @test runs[1][3][3] ≈ 1 - (1 - 1 / 400) * (1 - 0.02)^3
end
