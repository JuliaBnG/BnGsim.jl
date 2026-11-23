# Repeated records of several traits: single-record and repeated traits
# mixed, permanent environments correlated across traits, residuals
# correlated among the records of one animal and step; the evaluation
# against Henderson's equations solved densely.

@testset "Repeated records: several traits" begin
    sp = GenericSpecies("toy", [20_000_000, 10_000_000]; morgan = [1.0, 0.5])
    rng = Xoshiro(31)
    vm = VariantMap(unique(sort(rand(rng, 1:30_000_000, 1500))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 1200)))
    mm = MeiosisMap(sp, vm)
    ts = [Trait("ch4"; h² = 0.25, age = 1), Trait("growth"; h² = 0.3), Trait("carcass"; h² = 0.4, age = 2)]
    Rg = [1 0.3 0.2; 0.3 1 0.5; 0.2 0.5 1]
    Re = [1 0.3 0; 0.3 1 0.1; 0 0.1 1]
    Rpe = [1 0.4 0; 0.4 1 0; 0 0 1]
    mq = eqtl(hap, 1:5:length(vm), Rg, ts; rng = Xoshiro(4))
    pop = founder_population(hap, mq, ts, SimKey(6); Re = Re, herds = 4, c2 = 0.1,
        records = [3, 2, 1], repeatability = [0.6, 0.4, nothing], Rpe = Rpe)
    Vp = 1 ./ [0.25, 0.3, 0.4]
    @test pop.nrec == [3, 2, 1]
    @test pop.σpe ≈ sqrt.(Vp .* [0.35, 0.1, 0.0])
    @test pop.R0 ≈ sqrt.(Vp .* [0.4, 0.6, 0.6]) .* Re .* sqrt.(Vp .* [0.4, 0.6, 0.6])'
    @test pop.PE0[1:2, 1:2] ≈ (pop.σpe[1:2] .* Rpe[1:2, 1:2] .* pop.σpe[1:2]') && all(iszero, pop.PE0[3, :])
    # permanent environments with the specified correlation (600 founders)
    @test abs(cor(pop.pe[:, 1], pop.pe[:, 2]) - 0.4) < 4 * (1 - 0.16) / sqrt(600)
    @test all(iszero, pop.pe[:, 3])

    sc = Scheme(; noff = 600, nsire = 30, ndam = 200, evaluation = :pblup, max_sire_age = 3, max_dam_age = 3)
    steps = 4
    simulate!(pop, sc, mm, steps)
    R = pop.records
    # records per trait: one per step alive from the trait's age, at most `records`
    for i in (1, 2, 601, 1300, length(pop) - 1), j = 1:3
        mine = findall(q -> R.id[q] == i && R.trait[q] == j, eachindex(R.id))
        a = pop.generation[i] + Int(ts[j].age)
        @test R.step[mine] == collect(a:min(pop.exit[i] - 1, steps - 1, a + pop.nrec[j] - 1))
    end
    # temporary residuals of records made at the same step: covariance R0
    e = [R.y[q] - ts[R.trait[q]].μ - pop.tbv[R.id[q], R.trait[q]] - pop.pe[R.id[q], R.trait[q]] -
         BnGsim._cg(pop, R.id[q], R.trait[q], R.step[q]) for q in eachindex(R.id)]
    key(q) = (R.id[q], R.step[q])
    byunit = Dict{Tuple{Int32,Int32},Vector{Int}}()
    for q in eachindex(R.id)
        push!(get!(byunit, key(q), Int[]), q)
    end
    pairs12 = [(e[first(filter(q -> R.trait[q] == 1, v))], e[first(filter(q -> R.trait[q] == 2, v))])
               for v in values(byunit) if count(q -> R.trait[q] in (1, 2), v) == 2]
    x1, x2 = first.(pairs12), last.(pairs12)
    np = length(pairs12)
    @test np > 500
    @test abs(cor(x1, x2) - 0.3) < 4 * (1 - 0.09) / sqrt(np)
    @test abs(var(x1) / pop.R0[1, 1] - 1) < 4 * sqrt(2 / np) && abs(var(x2) / pop.R0[2, 2] - 1) < 4 * sqrt(2 / np)

    # the evaluation equals Henderson's equations solved densely
    small = founder_population(Haplotype(BitMatrix(rand(Xoshiro(7), Bool, length(vm), 80))), mq, ts, SimKey(9);
        Re = Re, herds = 3, c2 = 0.1, records = [3, 2, 1], repeatability = [0.6, 0.4, nothing], Rpe = Rpe)
    ssc = Scheme(; noff = 40, nsire = 4, ndam = 12, evaluation = :pblup, max_sire_age = 3, max_dam_age = 3)
    simulate!(small, ssc, mm, 4)
    evaluate!(small, ssc)
    Rs = small.records
    n, k = length(small), 3
    recs = collect(eachindex(Rs.id))
    y = Rs.y[recs]
    gcols = Dict{Tuple{Int,Int32,Int32},Int}()
    for q in recs
        get!(gcols, (Int(Rs.trait[q]), small.herd[Rs.id[q]], Rs.step[q]), length(gcols) + 1)
    end
    X = zeros(length(recs), length(gcols))
    Za = zeros(length(recs), k * n)
    pej = [1, 2]
    Zp = zeros(length(recs), length(pej) * n)
    for (r, q) in enumerate(recs)
        j, i = Int(Rs.trait[q]), Int(Rs.id[q])
        X[r, gcols[(j, small.herd[i], Rs.step[q])]] = 1
        Za[r, (j-1)*n+i] = 1
        j in pej && (Zp[r, (findfirst(==(j), pej)-1)*n+i] = 1)
    end
    Rm = [Rs.id[a] == Rs.id[b] && Rs.step[a] == Rs.step[b] ? small.R0[Rs.trait[a], Rs.trait[b]] : 0.0
          for a in recs, b in recs]
    Ri = inv(Rm)
    W = [X Za Zp]
    Ainv = Matrix(ainv(small.sire, small.dam))
    pa, pp = size(X, 2), k * n
    C = W' * Ri * W
    C[pa+1:pa+pp, pa+1:pa+pp] += kron(inv(small.G0), Ainv)
    C[pa+pp+1:end, pa+pp+1:end] += kron(inv(small.PE0[pej, pej]), Matrix(1.0I, n, n))
    sol = C \ (W' * Ri * y)
    @test small.ebv ≈ reshape(sol[pa+1:pa+pp], n, k)

    # scenario keys
    d = scenario_doc()
    d["traits"][1]["records"] = 3
    d["traits"][1]["repeatability"] = 0.5
    d["traits"][2]["records"] = 2
    d["traits"][2]["repeatability"] = 0.7
    d["pe_correlation"] = [[1.0, 0.2], [0.2, 1.0]]
    d["scheme"]["evaluation"] = "pblup"
    s = parse_scenario(d)
    @test s.pe_correlation == [1.0 0.2; 0.2 1.0] && [t.records for t in s.traits] == [3, 2]
    @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    rep = run_replicate(s)
    @test all(c -> c in names(rep.summary), ["vpe_a", "vpe_b"])
    d["pe_correlation"] = [[1.0, 1.2], [1.2, 1.0]]
    @test_throws ScenarioError parse_scenario(d)
    @test_throws ArgumentError founder_population(hap, mq, ts, SimKey(6); records = [3, 2, 1],
        repeatability = [0.6, 0.4, 0.5])                  # carcass has one record
end
