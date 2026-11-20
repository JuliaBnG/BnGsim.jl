# Fits the coefficients of `resource_estimate` (src/estimate.jl) to measured.tsv and
# prints them with the error on every calibration and held-out scenario:
#
#   julia --project=../.. fit.jl .  > fit.txt
#
# Memory: relative least squares, non-negative, on the mean of the two runs
# of each calibration scenario (h*, p*, gblup*, ssgblup*); the held-out
# scenarios (v*) are not used in the fit. Time: relative least squares on
# the faster run.
using BnGsim, Statistics, Printf, LinearAlgebra
dir = ARGS[1]
rows = [split(l, '\t') for l in eachline(joinpath(dir, "measured.tsv")) if occursin(".json", l)]
scen(f) = joinpath(dir, "scenarios", f)
data = Dict{String,Vector{Float64}}(); sec = Dict{String,Vector{Float64}}()
for r in rows
    push!(get!(data, r[1], Float64[]), parse(Float64, r[2]) * 2^20)
    push!(get!(sec, r[1], Float64[]), parse(Float64, r[11]))
end
feat(s) = begin
    z = ScenarioSize(s); ev = s.scheme.evaluation; d = 8.0 * z.genotyped^2
    [1.0, 8.0 * s.founders.ne * z.alleles, z.animals * 2 * z.loci / 8, z.evaluated, 8.0 * z.cohort^2,
     ev === :gblup ? d : 0.0, ev === :ssgblup ? d : 0.0]
end
names_ = ["baseline", "history", "genomes", "pedigree", "cohort", "dense_gblup", "dense_ssgblup"]
cal = sort([k for k in keys(data) if !startswith(k, "v")]); hold = sort([k for k in keys(data) if startswith(k, "v")])
X = reduce(vcat, [feat(read_scenario(scen(f)))' for f in cal]); y = [mean(data[f]) for f in cal]
# relative least squares, non-negative by active set
W = Diagonal(1 ./ y)
active = collect(1:size(X, 2)); β = zeros(size(X, 2))
while true
    b = (W * X[:, active]) \ (W * y)
    if all(>=(0), b); β .= 0; β[active] = b; break; end
    deleteat!(active, findfirst(<(0), b))
end
for (n, b) in zip(names_, β); @printf("  %-14s %.6g\n", n, b); end
pred(f) = dot(feat(read_scenario(scen(f))), β)
worst = 0.0
println("\nscenario           runs (MB)              expected (MB)   rel. error of mean   max ratio")
for f in [cal; hold]
    p = pred(f); m = data[f]
    global worst = max(worst, maximum(m) / p)
    @printf("  %-16s %-22s %8.0f   %+6.1f%%   %.3f %s\n", f, join(round.(Int, m ./ 2^20), ", "), p / 2^20,
        100(mean(m) - p) / p, maximum(m) / p, f in hold ? "(held out)" : "")
end
@printf("\nlargest measured/expected over all runs: %.3f\n", worst)
# time: seconds ≈ a·carried + steps·(b·evaluated + c·genotyped³)
tf(s) = begin
    z = ScenarioSize(s); f = s.founders; t = f.nt; μL = f.mr * 1e-8 * sum(s.species.chromosome)
    [2f.ne * μL * sum(min(g, 4f.ne) for g = 1:t), s.steps * z.evaluated, s.steps * float(z.genotyped)^3,
     s.scheme.evaluation in (:gblup, :ssgblup) ? 1.0 : 0.0]
end
T = reduce(vcat, [tf(read_scenario(scen(f)))' for f in cal]); ty = [minimum(sec[f]) for f in cal]
γ = (Diagonal(1 ./ ty) * T) \ (Diagonal(1 ./ ty) * ty); println("time coefficients ", γ)
rat = Float64[]
for f in [cal; hold]; m = minimum(sec[f]); p = dot(tf(read_scenario(scen(f))), γ); m > 1 && push!(rat, m / p); end
@printf("ratio measured/model for runs over 1 s: %.2f to %.2f\n", extrema(rat)...)
for f in [cal; hold]; @printf("  %-16s measured %5.1f s  model %5.1f s\n", f, minimum(sec[f]), dot(tf(read_scenario(scen(f))), γ)); end
