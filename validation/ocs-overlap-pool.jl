# Pool the per-replicate metrics written by validation/ocs-overlap.jl on
# several hosts and summarise them (paired over replicates).
#
#   julia --project=<env with DataFrames> validation/ocs-overlap-pool.jl file1.csv file2.csv ...
using Statistics, Printf

rows = [split(l, ",") for f in ARGS for l in eachline(f)]
methods = ["cohort", "age_year", "age_generation", "truncation"]
get_(m, k) = Dict(parse(Int, r[1]) => parse(Float64, r[k]) for r in rows if r[2] == m)
reps = sort(collect(intersect((keys(get_(m, 3)) for m in methods)...)))
se(x) = std(x) / sqrt(length(x))
println("# pooled: ", length(reps), " replicates (", minimum(reps), " to ", maximum(reps), ")")
@printf("%-15s %18s %18s %18s %14s %10s %10s\n", "method", "ΔF/generation", "Δθ/generation", "gain/year",
    "L", "parents", "infeasible")
for m in methods
    v(k) = (d = get_(m, k); [d[r] for r in reps])
    @printf("%-15s %9.4f (%.4f) %9.4f (%.4f) %9.4f (%.4f) %6.2f (%.2f) %10.1f %10.2f\n", m,
        mean(v(3)), se(v(3)), mean(v(4)), se(v(4)), mean(v(5)), se(v(5)), mean(v(6)), se(v(6)), mean(v(7)), mean(v(8)))
end
for m in ("cohort", "age_generation")
    d(k) = (a = get_(m, k); b = get_("age_year", k); [a[r] - b[r] for r in reps])
    @printf("%-15s − age_year: ΔF/gen %+.4f (%.4f), Δθ/gen %+.4f (%.4f), gain/year %+.4f (%.4f), L %+.3f (%.3f)\n",
        m, mean(d(3)), se(d(3)), mean(d(4)), se(d(4)), mean(d(5)), se(d(5)), mean(d(6)), se(d(6)))
end
