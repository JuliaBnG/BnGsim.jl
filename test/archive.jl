# Tier 2 archives and infinitesimal breeding values.

# Dosages from a PLINK 1 .bed file (SNP-major), read independently of the
# writer: two bits per animal from the low bits, 00 = two copies of A1,
# 10 = one, 11 = none, 01 = missing.
function read_bed(path, nid, nsnp)
    b = read(path)
    @test b[1:3] == UInt8[0x6c, 0x1b, 0x01]
    nb = cld(nid, 4)
    D = zeros(Int, nid, nsnp)
    for l = 1:nsnp, k = 1:nid
        byte = b[3+(l-1)*nb+(k-1)÷4+1]
        c = (byte >> (2 * ((k - 1) % 4))) & 0x03
        D[k, l] = c == 0x00 ? 2 : c == 0x02 ? 1 : c == 0x03 ? 0 : -1
    end
    return D
end

@testset "Archive" begin
    t = parse_scenario(scenario_doc())
    r = run_replicate(t, 1)
    pop = r.population
    dir = mktempdir()
    files = write_archive(dir, pop, r.map, r.species; panel=r.chip, qtl=r.qtl)
    @test basename.(files) == ["pedigree.csv", "records.csv", "truth.csv", "traits.csv", "qtl.csv",
        "genotypes.fam", "genotypes.bim", "genotypes.bed"]
    ped = readlines(joinpath(dir, "pedigree.csv"))
    @test ped[1] == "id,sire,dam,sex,generation,exit" && length(ped) == length(pop) + 1
    f = split(ped[end], ',')
    n = length(pop)
    @test parse.(Int, f[1:3]) == [n, pop.sire[n], pop.dam[n]] && f[4] == (pop.male[n] ? "M" : "F")
    # records: present exactly where observed
    rec = [split(l, ',') for l in readlines(joinpath(dir, "records.csv"))]
    @test rec[1] == ["id", "a", "b"]
    obs = BnGsim.observed(pop)
    @test all((rec[i+1][j+1] != "") == obs[i, j] for i = 1:n, j = 1:2)
    @test all(parse(Float64, rec[i+1][2]) == pop.phenotype[i, 1] for i = 1:n if obs[i, 1])
    tr = [split(l, ',') for l in readlines(joinpath(dir, "truth.csv"))]
    @test tr[1] == ["id", "tbv_a", "tbv_b", "ebv_a", "ebv_b"] && parse(Float64, tr[2][2]) == pop.tbv[1, 1]
    q = [split(l, ',') for l in readlines(joinpath(dir, "qtl.csv"))]
    @test length(q) == length(r.qtl) + 1 && parse(Float64, q[2][4]) == pop.qtl[1].additive[1]
    # genotypes: the .bed decoded independently equals the genomes
    ids = retained(pop)
    h = genomes(pop, ids)
    D = read_bed(joinpath(dir, "genotypes.bed"), length(ids), length(r.chip))
    @test D == [Int(h[l, 2i-1]) + Int(h[l, 2i]) for i in eachindex(ids), l in r.chip]
    bim = [split(l, '\t') for l in readlines(joinpath(dir, "genotypes.bim"))]
    @test length(bim) == length(r.chip) && bim[1][2] == "l$(r.chip[1])" && bim[1][5] == string(r.map.alt[r.chip[1]])
    @test parse(Int, bim[1][4]) == chrpos(r.map, r.species)[r.chip[1]]
    @test length(readlines(joinpath(dir, "genotypes.fam"))) == length(ids)
    rm(dir; recursive=true)
end

@testset "Infinitesimal breeding values" begin
    # Three generations with full and half sibs; 7 is from a full-sib mating
    # (F = 1/4) and 9 from parents related by 1/2 (F = 1/4).
    sire = Int32[0, 0, 0, 1, 1, 1, 4, 4, 7]
    dam = Int32[0, 0, 0, 2, 2, 3, 5, 6, 8]
    A = nrm(sire, dam)
    G0 = [1.0 0.5; 0.5 2.0]
    R0 = [1.0 -0.3; -0.3 0.5]
    reps = 4000
    U = zeros(reps, 9, 2)
    E = zeros(reps, 9, 2)
    for k = 1:reps
        u, y = infinitesimal(sire, dam, G0, R0, SimKey(11; replicate=k); mean=[10.0, 0.0])
        U[k, :, :] = u
        E[k, :, :] = y .- u .- [10.0 0.0]
    end
    # Cov(u_i^s, u_j^t) = A_ij G0_st, within 4 standard errors; the same for
    # the residuals (R0, independent between animals).
    z(x, e) = (mean(x) - e) / (std(x) / sqrt(length(x)))
    for (i, j) in ((1, 1), (4, 4), (9, 9), (1, 4), (4, 5), (4, 6), (7, 8), (1, 9)), (s, tt) in ((1, 1), (1, 2), (2, 2))
        @test abs(z(U[:, i, s] .* U[:, j, tt], A[i, j] * G0[s, tt])) < 4
    end
    @test A[7, 7] ≈ 1.25 && A[9, 9] ≈ 1.25
    for (s, tt) in ((1, 1), (1, 2), (2, 2))
        @test abs(z(E[:, 3, s] .* E[:, 3, tt], R0[s, tt])) < 4
    end
    @test abs(z(E[:, 3, 1] .* E[:, 4, 1], 0.0)) < 4
    # Keyed by animal: the values of an animal do not depend on the others.
    u1, _ = infinitesimal(sire, dam, G0, R0, SimKey(5))
    u2, _ = infinitesimal(sire[1:3], dam[1:3], G0, R0, SimKey(5))
    @test u1[1:3, :] == u2
    @test_throws ArgumentError infinitesimal(Int32[0, 3, 0], Int32[0, 0, 0], G0, R0, SimKey(1))
    @test_throws ArgumentError infinitesimal(sire, dam, G0, ones(3, 3), SimKey(1))
end

@testset "CLI: --archive" begin
    dir = mktempdir()
    chipped = tempname() * ".json"
    write(chipped, replace(read(MINIMAL, String),
        "\"loci\": {\"qtl\": {\"nlc\": 50, \"maf\": 0.01}}" =>
            "\"loci\": {\"qtl\": {\"nlc\": 50, \"maf\": 0.01}, \"chip\": {\"nlc\": 200, \"maf\": 0.05}}",
        "\"steps\": 3," => "\"steps\": 3, \"replicates\": 2,"))
    a = joinpath(dir, "a")
    code, v, _ = run_cli("run", chipped, "--out", a, "--archive", "chip")
    @test code == 0
    @test sort(readdir(joinpath(a, "archive"))) == ["r1", "r2"] && isfile(joinpath(a, "archive", "r2", "genotypes.bed"))
    p = JSON.parse(read(joinpath(a, "provenance.json"), String))
    @test p["archive"]["panel"] == "chip" && haskey(p["archive"]["files"], joinpath("archive", "r1", "pedigree.csv"))
    @test p["archive"]["files"][joinpath("archive", "r1", "genotypes.bed")] ==
          bytes2hex(open(BnGsim.sha256, joinpath(a, "archive", "r1", "genotypes.bed")))
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "b"), "--archive", "chip")[1] == 4   # no chip
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "c"), "--archive", "everything")[1] == 2
    @test run_cli("run", MINIMAL, "--out", joinpath(dir, "d"), "--archive", "none")[1] == 0
    @test sort(readdir(joinpath(dir, "d", "archive", "r1"))) == ["pedigree.csv", "qtl.csv", "records.csv", "traits.csv", "truth.csv"]
    # Stopped in replicate 2 and resumed: the same archives, byte for byte.
    s = read_scenario(chipped)
    e = joinpath(dir, "e")
    mkpath(joinpath(e, ".bngsim-run"))
    run = (scenario=scenario_json(s), replicates=[1, 2], checkpoint_every=0, archive=:chip, sessions=Any[])
    BnGsim._save_state(joinpath(e, ".bngsim-run", "run.jls"), run)
    @test BnGsim._drive(e, s, run, devnull, devnull; created=false, should_stop=(r, g) -> r == 2 && g == 2) == 7
    @test run_cli("resume", e)[1] == 0
    for f in BnGsim._files_under(joinpath(a, "archive"))
        @test read(f) == read(joinpath(e, relpath(f, a)))
    end
    rm(dir; recursive=true)
end
