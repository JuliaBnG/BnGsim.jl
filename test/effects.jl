# QTL effect distributions (Phase 2): the samplers, eqtl keeping their
# shape, and the scenario families.

@testset "Effect distributions" begin
    n = 10^6
    rng = Xoshiro(2026)
    # random-sign Gamma: symmetric, |x| ~ Gamma(0.4, 1)
    rs = RandomSign(Gamma(0.4, 1.0))
    x = rand(rng, rs, n)
    @test abs(count(<(0), x) / n - 0.5) < 5 * sqrt(0.25 / n)
    @test abs(mean(abs.(x)) - 0.4) < 5 * sqrt(0.4 / n)
    @test abs(mean(x)) < 5 * sqrt(var(rs) / n) && var(rs) ≈ 0.4 + 0.16
    # spike and slab: a fraction 1 - π exactly zero, the rest N(0, sd)
    ss = SpikeSlab(0.2, 2.0)
    x = rand(rng, ss, n)
    @test abs(count(iszero, x) / n - 0.8) < 5 * sqrt(0.16 / n)
    @test abs(var(x) - var(ss)) < 5 * var(ss) * sqrt(2 / (0.2n)) && var(ss) ≈ 0.8
    x = rand(rng, SpikeSlab(0.5, 1.0; spike_sd = 0.1), n)
    @test count(iszero, x) == 0 && var(SpikeSlab(0.5, 1.0; spike_sd = 0.1)) ≈ 0.505
    @test_throws ArgumentError SpikeSlab(0.0, 1.0)
    @test_throws ArgumentError SpikeSlab(0.5, 1.0; spike_sd = -1)

    # eqtl only rescales the draws: zeros stay zero, signs are kept
    hap = Haplotype(BitMatrix(rand(Xoshiro(1), Bool, 400, 300)))
    for da in (rs, ss, Laplace(0.0, 1.0))
        t = Trait("y"; h² = 0.3, da = da)
        q = eqtl(hap, 1:400, t; rng = Xoshiro(5))
        raw = rand(Xoshiro(5), da, 400)
        c = q.additive[findfirst(!iszero, raw)] / raw[findfirst(!iszero, raw)]
        @test c > 0 && q.additive ≈ c .* raw
    end

    # scenario families: canonical round trip and errors
    base = JSON.parse(read(MINIMAL, String))
    for da in (Dict("family" => "laplace", "scale" => 2.0),
               Dict("family" => "gamma", "shape" => 0.4, "sign" => "random"),
               Dict("family" => "spike_slab", "pi" => 0.1, "spike_sd" => 0.01))
        d = deepcopy(base)
        d["traits"][1]["da"] = da
        s = parse_scenario(d)
        @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    end
    d = deepcopy(base)
    d["traits"][1]["da"] = Dict("family" => "gamma", "shape" => 0.4, "sign" => "random")
    @test BnGsim._dist(parse_scenario(d).traits[1].da) isa RandomSign
    for bad in (Dict("family" => "spike_slab", "pi" => 0.0), Dict("family" => "spike_slab", "pi" => 1.5),
                Dict("family" => "spike_slab", "pi" => 0.5, "spike_sd" => -1.0),
                Dict("family" => "laplace", "location" => 1.0), Dict("family" => "gamma", "shape" => 1.0, "sign" => "both"),
                Dict("family" => "cauchy"))
        d = deepcopy(base)
        d["traits"][1]["da"] = bad
        @test_throws ScenarioError parse_scenario(d)
    end
end
