# Fits the time model of `resource_estimate` (0.37, src/estimate.jl) to the
# phase timings of run_phases.sh and prints the coefficients and the error
# on every scenario:
#
#   julia --project=../.. fit_phases.jl phases.tsv > fit_phases.txt
#
# Each phase is fitted on its own timings (the warm replicate, replicate 2,
# mean of the two processes) by non-negative least squares on relative
# error, with errors of runs under 0.2 s counted at 0.2 s so that timer
# noise does not drive the fit; compilation is fitted to the difference
# between the cold and the warm replicate. Held out: the examples, the
# held-out scenarios v1-v6 of the memory calibration, and check/.
using BnGsim, Statistics, Printf, LinearAlgebra
const COLS = (:history, :setup, :evaluate, :kinv, :reliability, :select, :mating, :meiosis, :summary, :loop)
const DIRS = ["scenarios", "later", "phases", "check", joinpath("..", "..", "examples")]
here = @__DIR__
findscen(f) = only(p for p in (joinpath(here, d, f) for d in DIRS) if isfile(p))
heldout(f) = occursin(r"^v\d\.json$", f) || dirname(findscen(f)) in (joinpath(here, "check"), joinpath(here, "..", "..", "examples"))

# rows: run, file, replicate, then the columns of COLS
T = Dict{Tuple{String,Int},Vector{NamedTuple}}()
rss = Dict{String,Vector{Float64}}()                  # peak RSS of each process (MiB), after replicate 2
for l in eachline(ARGS[1])
    r = split(l, '\t')
    length(r) >= 13 && endswith(r[2], ".json") || continue
    v = NamedTuple{COLS}(Tuple(parse.(Float64, r[4:13])))
    push!(get!(T, (r[2], parse(Int, r[3])), NamedTuple[]), v)
    r[3] == "2" && push!(get!(rss, r[2], Float64[]), parse(Float64, r[end]))
end
files = sort(unique(first.(collect(keys(T)))))
S = Dict(f => read_scenario(findscen(f)) for f in files)
F = Dict(f => BnGsim._time_features(S[f]) for f in files)
avg(f, r, c) = mean(getfield(x, c) for x in T[(f, r)])
cal = filter(!heldout, files)

function nnls_rel(X, y; floor = 0.2)
    W = Diagonal(1 ./ max.(y, floor))
    active = findall(j -> any(!iszero, view(X, :, j)), 1:size(X, 2)); β = zeros(size(X, 2))
    while !isempty(active)
        b = (W * X[:, active]) \ (W * y)
        if all(>=(0), b); β[active] = b; break; end
        deleteat!(active, argmin(b))
    end
    return β
end

# the rest of the loop (records, phenotypes, selection, culling)
rest(f) = avg(f, 2, :loop) - sum(avg(f, 2, c) for c in (:evaluate, :reliability, :select, :mating, :meiosis, :summary))
# phase => (features, measured warm seconds)
phases = [
    :history => (f -> (c = BnGsim._carried(S[f]); [c, c * BnGsim._cache_factor(BnGsim.ScenarioSize(S[f]))]),
        f -> avg(f, 2, :history)),
    :setup => (f -> [F[f].setup, 1.0], f -> avg(f, 2, :setup)),
    :meiosis => (f -> [F[f].meiosis], f -> avg(f, 2, :meiosis)),
    :summary => (f -> [F[f].summary_genome, F[f].summary_pedigree], f -> avg(f, 2, :summary)),
    :evaluate => (f -> [F[f].equations, F[f].equations2, F[f].dense, F[f].dense2, F[f].reml_equations,
                        F[f].reml_dense, F[f].records],
        f -> avg(f, 2, :evaluate) - avg(f, 2, :kinv)),
    :kinv => (f -> [F[f].grm, F[f].inverse, F[f].grm_ss, F[f].inverse_ss], f -> avg(f, 2, :kinv)),
    :mating => (f -> [F[f].mating], f -> avg(f, 2, :mating)),
    :reliability => (f -> [F[f].reliability], f -> avg(f, 2, :reliability)),
    :rest => (f -> [float(S[f].steps), F[f].animals], f -> rest(f)),
]
coef = Dict{Symbol,Vector{Float64}}()
println("coefficients (seconds per unit, 8 threads)")
for (p, (x, y)) in phases
    use = filter(f -> any(!iszero, x(f)) || y(f) > 0.05, cal)
    X = reduce(vcat, [x(f)' for f in use]); Y = [y(f) for f in use]
    coef[p] = nnls_rel(X, Y)
    @printf("  %-9s %s\n", p, join((@sprintf("%.4g", b) for b in coef[p]), "  "))
end
# compilation: cold minus warm (setup and loop), by what the scenario uses
cflags(f) = (s = S[f]; ev = s.scheme.evaluation;
    [1.0, ev in (:pblup, :gblup, :ssgblup), ev in (:gblup, :ssgblup), ev === :ssgblup, length(s.traits) > 1,
     any(t -> t.records > 1, s.traits), s.scheme.reml, s.scheme.mating === :min_coancestry])
compile(f) = (avg(f, 1, :setup) + avg(f, 1, :loop)) - (avg(f, 2, :setup) + avg(f, 2, :loop))
coef[:compile] = nnls_rel(reduce(vcat, [cflags(f)' for f in cal]), [compile(f) for f in cal]; floor = 0.5)
@printf("  %-9s %s\n", :compile, join((@sprintf("%.3g", b) for b in coef[:compile]), "  "))

model(f, p) = dot(first(Dict(phases)[p])(f), coef[p])
function total(f)                                     # one replicate, warm, with its history
    sum(model(f, p) for p in first.(phases))
end
meas(f) = avg(f, 2, :history) + avg(f, 2, :setup) + avg(f, 2, :loop)
println("\nper replicate, warm (s): measured, model, ratio; then compile (cold - warm) measured, model")
ratios = Float64[]
for f in [cal; filter(heldout, files)]
    m, p = meas(f), total(f)
    m > 1 && push!(ratios, m / p)
    @printf("  %-24s %7.2f %7.2f %5.2f   %6.2f %6.2f %s\n", f, m, p, m / p, compile(f), dot(cflags(f), coef[:compile]),
        heldout(f) ? "(held out)" : "")
end
@printf("\nmeasured/model for replicates over 1 s: %.2f to %.2f\n", extrema(ratios)...)

# The dense matrices of several traits and of REML (0.37): what the peak
# exceeds the rest of the memory model by, per byte of the dense
# genotyped × genotyped matrix, against k² - 1 and REML × k², by
# evaluation (calibration scenarios with several traits or REML).
println("\nmemory of several traits and REML, bytes per byte of the dense matrix: k² - 1, REML × k²")
for ev in (:gblup, :ssgblup)
    use = filter(f -> S[f].scheme.evaluation === ev && (length(S[f].traits) > 1 || S[f].scheme.reml), cal)
    isempty(use) && continue
    k(f) = length(S[f].traits)
    X = reduce(vcat, [[k(f)^2 - 1.0, S[f].scheme.reml * k(f)^2]' for f in use])
    nd(f) = 8.0 * BnGsim.ScenarioSize(S[f]).genotyped^2
    base(f) = (t = BnGsim._memory_terms(S[f], BnGsim.ScenarioSize(S[f])); sum(values(t)) - t.dense +
               (ev === :gblup ? BnGsim._MEM.dense_gblup : BnGsim._MEM.dense_ssgblup) * nd(f))
    Y = [(mean(rss[f]) * 2^20 - base(f)) / nd(f) for f in use]
    b = nnls_rel(X, Y; floor = 1.0)
    @printf("  %-8s %.3g  %.3g   (%s)\n", ev, b..., join(use, ", "))
end

# The memory model (calibrated on svr, 0.15) on this host: the peak of each
# process against the expectation and the bound.
println("\npeak resident memory (MiB): processes, expected, mean/expected, max/bound")
mr = Float64[]
for f in files
    e = resource_estimate(S[f]); m = rss[f]
    push!(mr, maximum(m) * 2^20 / e.peak_rss_bound_bytes)
    @printf("  %-24s %-16s %7.0f %6.2f %6.2f\n", f, join(round.(Int, m), ", "), e.peak_rss_bytes / 2^20,
        mean(m) * 2^20 / e.peak_rss_bytes, maximum(m) * 2^20 / e.peak_rss_bound_bytes)
end
@printf("largest max/bound: %.2f\n", maximum(mr))
