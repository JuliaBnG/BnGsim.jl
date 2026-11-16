# Realized relationships of half and full sibs from BnGsim's meiosis
# against their exact distribution moments under Haldane's map.
#
#   julia -t 4 --project=<env with BnGsim, BnGStructs> validation/ibd.jl [replicates]
#
# Theory. Two gametes of one parent carry the same parental haplotype at
# a position t with probability ½; along a chromosome this indicator is a
# two-state Markov chain switching at rate 2 per Morgan (each gamete
# switches at rate 1, independently; independent assortment gives the
# stationary start). The proportion p_c of chromosome c (L_c Morgans)
# shared has mean ½ and
#     Var(p_c) = (L_c/2 − (1 − e^{−4L_c})/8) / (4 L_c²).
# Realized relationships weight bases: p_G = Σ B_c p_c / B, B_c bases of
# chromosome c, so Var(p_G) = Σ B_c² Var(p_c) / B². Paternal half sibs
# with unrelated dams: R = p_G/2 (mean ¼, variance Var(p_G)/4); full sibs
# of unrelated parents: R = (p_G + p_G′)/2 (mean ½, variance Var(p_G)/2).
#
# Scenario: 3 chromosomes of 50, 100 and 150 Mb with 0.5, 1.5 and 1.0 M
# (unequal cM/Mb), 2000 loci; 200 unrelated founders; one step of random
# mating, 10 sires and 50 dams, 200 offspring. Per replicate the pairs
# taken share no offspring, so they are independent given the founders:
# two full-sib pairs per dam family, and half-sib pairs of one sire with
# different dams. Pass: |mean − theory| < 4 s.e. and
# |variance − theory| < 4 s.e. for both kinds.
using BnGsim, BnGStructs, Random, Statistics, Printf

nrep = parse(Int, get(ARGS, 1, "50"))
B = [50_000_000, 100_000_000, 150_000_000]
L = [0.5, 1.5, 1.0]
sp = GenericSpecies("val", B; morgan = L)
vp(l) = (l / 2 - (1 - exp(-4l)) / 8) / (4l^2)
VpG = sum(B .^ 2 .* vp.(L)) / sum(B)^2
theory = Dict(:half => (0.25, VpG / 4), :full => (0.5, VpG / 2))

rng = Xoshiro(1)
vm = VariantMap(sort(unique(rand(rng, 1:sum(B), 2000))), sp)
mm = MeiosisMap(sp, vm)
hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 400)))
t = Trait("y"; h² = 0.3)
qtl = eqtl(hap, 1:10:length(vm), t; rng = Xoshiro(2))
sc = Scheme(; noff = 200, nsire = 10, ndam = 50, evaluation = :tbv, selection = :random, record_meiosis = true)
r = Dict(:half => Float64[], :full => Float64[])
t0 = time()
for rep = 1:nrep
    pop = founder_population(hap, qtl, t, SimKey(77; replicate = rep))
    step!(pop, sc, mm)
    ids = collect(generation_ids(pop, 1))
    R = irm(founder_segments(pop, mm); ids = ids)
    pos = Dict(i => k for (k, i) in enumerate(ids))
    used = Set{Int}()
    # full sibs: disjoint pairs within each dam family
    for d in unique(pop.dam[ids])
        fam = filter(i -> pop.dam[i] == d && !(i in used), ids)
        for k = 1:2:length(fam)-1
            a, b = fam[k], fam[k+1]
            pop.sire[a] == pop.sire[b] || continue
            push!(r[:full], R[pos[a], pos[b]]); push!(used, a, b)
        end
    end
    # half sibs: disjoint pairs of one sire with different dams, among the rest
    for s in unique(pop.sire[ids])
        fam = filter(i -> pop.sire[i] == s && !(i in used), ids)
        while length(fam) >= 2
            a = popfirst!(fam)
            k = findfirst(i -> pop.dam[i] != pop.dam[a], fam)
            k === nothing && break
            b = fam[k]; deleteat!(fam, k)
            push!(r[:half], R[pos[a], pos[b]]); push!(used, a, b)
        end
    end
end

println("# BnGsim ", pkgversion(BnGsim), ", Julia ", VERSION, ", ", nrep, " replicates, ",
    round(time() - t0; digits = 0), " s")
@printf("genome: B = %s bases, L = %s M; Var(p_G) = %.6f\n", B, L, VpG)
ok = true
for kind in (:half, :full)
    x = r[kind]
    n = length(x)
    μ, V = theory[kind]
    m, v = mean(x), var(x)
    se_m = sqrt(V / n)
    # s.e. of a sample variance: sqrt((μ₄ − σ⁴(n−3)/(n−1))/n), with the sample fourth moment
    μ4 = mean((x .- m) .^ 4)
    se_v = sqrt(max(μ4 - v^2 * (n - 3) / (n - 1), 0) / n)
    pm, pv = abs(m - μ) < 4se_m, abs(v - V) < 4se_v
    global ok &= pm & pv
    @printf("%-5s n = %5d  mean %.5f (theory %.5f, s.e. %.5f) %s   variance %.6f (theory %.6f, s.e. %.6f) %s\n",
        kind, n, m, μ, se_m, pm ? "ok" : "OFF", v, V, se_v, pv ? "ok" : "OFF")
end
println(ok ? "gate passed" : "gate NOT passed")
