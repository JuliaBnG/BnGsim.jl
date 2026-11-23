# Reproductive limits (Phase 3): minimum breeding ages, offspring per parent
# and time step, and the parity of dams.

@testset "Reproductive limits" begin
    # capped allocation: shares above the cap are set to it and the rest
    # spread in proportion; rounding never passes the cap; ties by id
    @test BnGsim._allocate([1, 2, 3, 4], [1.0, 1, 1, 7], 10; cap = 3) == [3, 2, 2, 3]
    @test BnGsim._allocate([1, 2, 3, 4], [1.0, 1, 1, 7], 10) == [1, 1, 1, 7]
    for seed = 1:50
        w = rand(Xoshiro(seed), 8)
        k = BnGsim._allocate(collect(1:8), w, 20; cap = 3)
        @test sum(k) == 20 && all(<=(3), k)
    end
    @test_throws ArgumentError BnGsim._allocate([1, 2, 3], ones(3), 10; cap = 3)

    d = JSON.parse(read(MINIMAL, String))
    d["founders"] = Dict("ne" => 200, "nt" => 300, "nid" => 200)
    d["steps"] = 8
    d["scheme"] = Dict{String,Any}("noff" => 50, "nsire" => 5, "ndam" => 30, "evaluation" => "pblup",
        "max_sire_age" => 4, "max_dam_age" => 5, "min_sire_age" => 1, "min_dam_age" => 2,
        "max_sire_offspring" => 12, "max_dam_offspring" => 2, "max_parity" => 2)
    s = parse_scenario(d)
    @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    p = run_replicate(s).population
    kids = findall(>(0), p.sire)
    for g in unique(p.generation[kids])
        born = filter(i -> p.generation[i] == g, kids)
        @test maximum(count(==(x), p.sire[born]) for x in unique(p.sire[born])) <= 12
        @test maximum(count(==(x), p.dam[born]) for x in unique(p.dam[born])) <= 2
    end
    # parity: the number of time steps with offspring, at most 2
    for i in unique(p.dam[kids])
        @test p.parity[i] == length(unique(p.generation[kids[p.dam[kids] .== i]])) <= 2
    end
    # minimum ages: a parent of offspring born at g + 1 was at least the
    # minimum age at g (founders excepted)
    for i in kids, (par, a) in ((p.sire[i], 1), (p.dam[i], 2))
        p.generation[par] == 0 || @test p.generation[i] - 1 - p.generation[par] >= a
    end
    # too few dams for the litter limit: refused, not truncated
    d["scheme"]["ndam"] = 10
    @test_throws ArgumentError run_replicate(parse_scenario(d))
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, max_dam_age = 1, min_dam_age = 2)
    @test_throws ArgumentError Scheme(; noff = 10, nsire = 2, ndam = 2, max_parity = -1)
end
