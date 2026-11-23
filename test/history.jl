# Founder histories (Phase 2): demography and a quantitative trait under
# selection in the historical population, whose loci become the QTL.

function _history_doc(; selection = true)
    d = JSON.parse(read(MINIMAL, String))
    d["founders"] = Dict{String,Any}("ne" => 80, "nt" => 150, "nid" => 40, "sexratio" => 0.3,
        "popsize" => Dict("bottleneck" => Dict("size" => 30, "start" => 100, "duration" => 5)))
    if selection
        d["founders"]["selection"] = Dict("nsites" => 400, "da" => Dict("family" => "normal", "sd" => 0.3),
            "mu" => 2e-4, "sigma_e" => 1.0, "fitness" => Dict("model" => "stabilizing", "theta" => 0.0, "omega2" => 2.0))
        d["loci"]["qtl"] = Dict("from" => "history", "maf" => 0.0)
    end
    d["loci"]["chip"] = Dict("nlc" => 40, "maf" => 0.05)
    d["scheme"]["evaluation"] = "gblup"
    return d
end

@testset "Founder histories" begin
    s = parse_scenario(_history_doc())
    @test s.founders.popsize.model === :bottleneck && s.founders.sexratio == 0.3 && s.qtl.from === :history
    @test scenario_json(parse_scenario(scenario_json(s))) == scenario_json(s)
    # demography: the census follows the trajectory, the sex ratio is fixed
    sp = BnGsim._species_of(s)
    h = BnGsim._history(s, sp, 1)
    traj = BnGsim._trajectory(s.founders.popsize, 80)
    N = BnGsim.census(h)
    @test N == [80; traj.(1:150)] && N[101:105] == fill(30, 5) && N[end] == 80
    @test h.records.males == round.(Int, 0.3 .* N[1:end-1])
    g = _history_doc(; selection = false)
    g["founders"]["popsize"] = Dict("growth" => Dict("rate" => 0.02, "limit" => 200))
    sg = parse_scenario(g)
    Ng = BnGsim.census(BnGsim._history(sg, sp, 1))
    @test issorted(Ng) && Ng[end] == BnGsim._trajectory(sg.founders.popsize, 80)(150) && Ng[end] <= 200

    # the trait's alleles merged into the founder haplotypes: rows of trait
    # loci carry the derived alleles of `selected_haplotypes`, neutral rows
    # the neutral ones, and the effects are those of the loci
    loci = BnGsim._trait_loci(s, 1)
    xy, lmp, rows, eff = BnGsim.founder_haplotypes(h, 40, loci; rng = Xoshiro(3))
    xy0, lmp0 = BnGsim.founder_haplotypes(h, 40; rng = Xoshiro(3))           # the same founders, neutral only
    @test !isempty(rows) && all(r -> insorted(lmp.pos[r], loci.pos), rows)
    @test eff == [loci.a[searchsortedfirst(loci.pos, lmp.pos[r]), 1] for r in rows]
    neutral = setdiff(1:size(xy, 1), rows)
    keep0 = findall(p -> !insorted(p, loci.pos), lmp0.pos)
    @test lmp.pos[neutral] == lmp0.pos[keep0] && xy[neutral, :] == xy0[keep0, :]
    ids = sort(BnGsim.sample(Xoshiro(3), 1:80, 40; replace = false))
    hs = vec([2ids' .- 1; 2ids'])
    for r in rows[1:min(end, 20)]
        @test xy[r, :] == [insorted(UInt32(lmp.pos[r]), sort(h.selected_haplotypes[k])) for k in hs]
    end

    # a replicate: the QTL are the segregating trait loci, their effects the
    # historical ones times one factor, the chip apart from them
    rep = BnGsim._setup_replicate(s, 1, h)
    q = rep.population.qtl[1]
    @test q.loci == rep.qtl && isempty(intersect(rep.qtl, rep.chip))
    a = [loci.a[searchsortedfirst(loci.pos, rep.map.pos[l]), 1] for l in q.loci]
    c = q.additive[1] / a[1]
    @test c > 0 && q.additive ≈ c .* a
    @test std(tbv(rep.population.blocks[1], q)) ≈ 1.0
    @test nrow(run_replicate(s).summary) == s.steps + 1

    # errors
    for (path, v) in ((["founders", "popsize"], Dict("bottleneck" => Dict("size" => 30, "start" => 1, "duration" => 1),
                          "growth" => Dict("rate" => 0.1))),
                      (["loci", "qtl"], Dict("from" => "history", "rho" => 1.5)),
                      (["founders", "sexratio"], 1.0),
                      (["founders", "selection", "mu"], 0.0))
        d = _history_doc()
        x = d
        for k in path[1:end-1]
            x = x[k]
        end
        x[path[end]] = v
        @test_throws ScenarioError parse_scenario(d)
    end
    d = _history_doc(; selection = false)
    d["loci"]["qtl"] = Dict("from" => "history")
    @test_throws ScenarioError parse_scenario(d)                        # needs founders.selection
    d = _history_doc()
    d["traits"] = [Dict("name" => "a"), Dict("name" => "b")]
    @test_throws ScenarioError parse_scenario(d)                        # one trait
end
