# Meiosis models (Phase 3): sex-specific genetic maps and the gamma model
# of crossover interference, against crossover counts and the map
# function r(d) = ½(1 − P₀(d)) of a stationary gamma renewal process.

@testset "Meiosis: sex-specific maps and interference" begin
    sp = GenericSpecies("toy", [100_000_000]; morgan = [1.25])
    pos = collect(500_000:1_000_000:99_500_000)
    vm = VariantMap(pos, sp)
    nl = length(pos)
    par = Haplotype(BitMatrix([falses(nl) trues(nl)]))       # haplotype 1 all 0, haplotype 2 all 1
    rtheory(d, ν) = (θ = 1 / (2ν); 0.5 * (1 - 2 * (ν * θ * ccdf(Gamma(ν + 1, θ), d) - d * ccdf(Gamma(ν, θ), d))))
    @test rtheory(0.3, 1.0) ≈ 0.5 * (1 - exp(-0.6))                # Haldane
    n = 20_000
    for ν in (1.0, 4.3)
        mm = MeiosisMap(sp, vm; male = [1.0], female = [1.5], interference = ν)
        off = Haplotype(falses(nl, 2n))
        _, rec = make_progeny!(off, par, ones(Int, n, 2), mm, SimKey(11); record = true)
        xs = [ncrossovers(rec, g) for g = 1:2n]
        m, f = xs[1:2:end], xs[2:2:end]
        @test abs(mean(m) - 1.0) < 4 * std(m) / sqrt(n) && abs(mean(f) - 1.5) < 4 * std(f) / sqrt(n)
        ν == 1 ? (@test abs(var(m) - 1.0) < 4 * sqrt(2 / n) * 1.5) : (@test var(m) < 0.8)   # underdispersed
        H = off.gt[1:nl, 1:2:end]                                   # sire gametes
        for k in (20, 50)
            d = k / 100
            r = mean(H[1, :] .!= H[1+k, :])
            @test abs(r - rtheory(d, ν)) < 4 * sqrt(r * (1 - r) / n)
        end
    end
    # ν = 1 with equal maps is the model of earlier versions, draw for draw
    a, b = Haplotype(falses(nl, 200)), Haplotype(falses(nl, 200))
    make_progeny!(a, par, ones(Int, 100, 2), MeiosisMap(sp, vm), SimKey(3))
    make_progeny!(b, par, ones(Int, 100, 2), MeiosisMap(sp, vm; male = [1.25], female = [1.25], interference = 1), SimKey(3))
    @test a.gt == b.gt
    @test_throws ArgumentError MeiosisMap(sp, vm; female = [1.0, 2.0])
    @test_throws ArgumentError MeiosisMap(sp, vm; interference = 0)

    # scenario keys: the history uses the average map, the steps the sexes' maps
    d = JSON.parse(read(MINIMAL, String))
    d["species"]["morgan_female"] = [0.3, 0.2]
    d["species"]["morgan"] = [0.2, 0.1]
    delete!(d["species"], "M")
    d["species"]["interference"] = 2.5
    s = parse_scenario(d)
    @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    @test BnGsim.map_length(BnGsim._species_of(s)) ≈ [0.25, 0.15]
    rep = BnGsim._setup_replicate(s, 1, nothing)
    @test rep.meiosis.morgan == [0.2, 0.1] && rep.meiosis.fmorgan == [0.3, 0.2] && rep.meiosis.ν == 2.5
    @test nrow(run_replicate(s).summary) == s.steps + 1
    for (k, v) in (("morgan_female", [0.3]), ("interference", 0.0))
        e = JSON.parse(read(MINIMAL, String))
        e["species"][k] = v
        @test_throws ScenarioError parse_scenario(e)
    end
end
