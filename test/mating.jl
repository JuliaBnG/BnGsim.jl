# Mating policies (Phase 3): the transportation solver against brute
# force, and the pairings of mating_plan.

# every assignment of the dam slots to the sire slots (small cases)
function _brute_transport(a, b, C, cap)
    sl = reduce(vcat, [fill(s, a[s]) for s in eachindex(a)])
    best = Inf
    rem = copy(b)
    F = zeros(Int, length(a), length(b))
    function go(k, c)
        if k > length(sl)
            best = min(best, c)
            return
        end
        s = sl[k]
        for d in eachindex(b)
            rem[d] > 0 && F[s, d] < cap || continue
            rem[d] -= 1; F[s, d] += 1
            go(k + 1, c + C[s, d])
            rem[d] += 1; F[s, d] -= 1
        end
    end
    go(1, 0.0)
    return best
end

@testset "Mating: transportation problem" begin
    rng = Xoshiro(1)
    for trial = 1:300
        p, q = rand(rng, 1:3), rand(rng, 1:4)
        n = rand(rng, max(p, q):7)
        a = ones(Int, p); for _ = 1:n-p; a[rand(rng, 1:p)] += 1; end
        b = ones(Int, q); for _ = 1:n-q; b[rand(rng, 1:q)] += 1; end
        C = round.(rand(rng, p, q); digits = 1)        # ties on purpose
        cap = rand(rng, (1, 2, typemax(Int)))
        best = _brute_transport(a, b, C, cap)
        if isinf(best)
            @test_throws ArgumentError transport(a, b, C; cap = cap)
        else
            F = transport(a, b, C; cap = cap)
            @test vec(sum(F; dims = 2)) == a && vec(sum(F; dims = 1)) == b && all(0 .<= F .<= cap)
            @test sum(C .* F) ≈ best atol = 1e-9
            @test transport(a, b, C; cap = cap) == F          # deterministic
        end
    end
    @test_throws ArgumentError transport([1], [2], ones(1, 1))
    @test_throws ArgumentError transport([1], [1], -ones(1, 1))
    @test_throws DimensionMismatch transport([1, 1], [2], ones(1, 1))
end

@testset "Mating: plans" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(2)
    vm = VariantMap(unique(sort(rand(rng, 1:30_000_000, 1500))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 300)))
    t = Trait("y"; h² = 0.3)
    q = eqtl(hap, 1:5:length(vm), t; rng = Xoshiro(3))
    pop = founder_population(hap, q, t, SimKey(5))
    simulate!(pop, Scheme(; noff = 120, nsire = 8, ndam = 30, evaluation = :pblup), MeiosisMap(sp, vm), 3)
    sc = Scheme(; noff = 60, nsire = 6, ndam = 20, evaluation = :pblup)
    evaluate!(pop, sc)
    s, ns, d, nd, _ = select_parents(pop, sc, nothing)
    index = vec(pop.ebv)
    θ(plan) = (A = nrm(pop.sire, pop.dam, sort(unique(plan))); r = Dict(x => i for (i, x) in enumerate(sort(unique(plan))));
               sum(A[r[a], r[b]] / 2 for (a, b) in eachrow(plan)))
    counts(plan, col, ids) = [count(==(i), plan[:, col]) for i in ids]
    plans = Dict()
    for (m, pl) in ((:random, 0), (:random, 1), (:min_coancestry, 0), (:min_coancestry, 1), (:assortative, 0),
                    (:disassortative, 0))
        plan = mating_plan(pop, s, ns, d, nd; mating = m, pair_limit = pl, index = index)
        @test counts(plan, 1, s) == ns && counts(plan, 2, d) == nd
        pl > 0 && @test maximum(count(==(r), eachrow(plan)) for r in eachrow(plan)) <= pl
        @test plan == mating_plan(pop, s, ns, d, nd; mating = m, pair_limit = pl, index = index)
        plans[(m, pl)] = plan
    end
    # least total coancestry, exactly, against the transportation solution
    # and below the random pairings
    @test θ(plans[(:min_coancestry, 0)]) <= θ(plans[(:random, 0)]) + 1e-12
    @test θ(plans[(:min_coancestry, 1)]) <= θ(plans[(:random, 1)]) + 1e-12
    @test θ(plans[(:min_coancestry, 0)]) <= θ(plans[(:min_coancestry, 1)]) + 1e-12
    # best with best; best with worst
    pa, pd = plans[(:assortative, 0)], plans[(:disassortative, 0)]
    @test issorted(index[pa[:, 1]]; rev = true) && issorted(index[pa[:, 2]]; rev = true)
    @test issorted(index[pd[:, 1]]; rev = true) && issorted(index[pd[:, 2]])
    @test_throws ArgumentError mating_plan(pop, s, ns, d, nd; mating = :assortative)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, mating = :assortative, pair_limit = 1)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, mating = :nearest)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, pair_limit = -1)

    # scenario keys
    doc = JSON.parse(read(MINIMAL, String))
    doc["scheme"]["mating"] = "min_coancestry"
    doc["scheme"]["pair_limit"] = 1
    sm = parse_scenario(doc)
    @test sm.scheme.mating == :min_coancestry && sm.scheme.pair_limit == 1
    @test scenario_json(parse_scenario(scenario_json(sm))) == scenario_json(sm)
    @test nrow(run_scenario(sm)) == sm.steps + 1
    doc["scheme"]["mating"] = "nearest"
    @test_throws ScenarioError parse_scenario(doc)
end
