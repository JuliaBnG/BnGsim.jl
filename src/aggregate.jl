# Aggregation of replicates (docs/src/roadmap.md, Phase 5): per time step,
# the mean, standard deviation and standard error of every summary over
# replicates, with the change per step of the genetic means, the rates of
# inbreeding and the depletion of the genetic variance. The result does not depend on the order in which replicates
# finished or were merged: the rows are sorted by replicate first.

# The changes per step derived within each replicate: `dG` for `mean` (the
# index) and `dG_<trait>` for `mean_<trait>`; `dF`, the rate of inbreeding
# (F_t - F_{t-1}) / (1 - F_{t-1}), and `dFhet` the same of `Fhet`; NaN at
# a replicate's first step. `var_ratio` and `genic_ratio`: `var` and
# `genic` relative to the replicate's first step.
function _with_changes(df::DataFrame)
    d = sort(df, [:replicate, :generation])
    means = [c for c in names(d) if c == "mean" || startswith(c, "mean_")]
    starts = [i == 1 || d.replicate[i] != d.replicate[i-1] for i = 1:nrow(d)]
    for c in means
        x = Float64.(d[!, c])
        d[!, replace(c, "mean" => "dG"; count=1)] = [starts[i] ? NaN : x[i] - x[i-1] for i = 1:nrow(d)]
    end
    for (c, r) in (("F", :dF), ("Fhet", :dFhet))
        c in names(d) || continue
        F = Float64.(d[!, c])
        d[!, r] = [starts[i] ? NaN : (F[i] - F[i-1]) / (1 - F[i-1]) for i = 1:nrow(d)]
    end
    first_row = accumulate((j, i) -> starts[i] ? i : j, 1:nrow(d); init = 1)
    for c in ("var", "genic")
        c in names(d) || continue
        x = Float64.(d[!, c])
        d[!, c*"_ratio"] = x ./ x[first_row]
    end
    return d
end

"""
    aggregate(summary::DataFrame) -> DataFrame

Aggregate the per-step summaries of several replicates (the table of
[`run_scenario`](@ref), with a column `replicate`). One row per time step
(`generation`), with `replicates`, the number of replicates at that step,
and for every other column `c` of the summary `c_mean`, `c_sd` (with the
`n - 1` divisor) and `c_se` (`c_sd / √n`) over the replicates, `NaN`
values left out (`n` then counts the others). Also aggregated are the
changes per step within each replicate: `dG` (of `mean`, the index of
the traits) and `dG_<trait>` (of `mean_<trait>`) for the genetic means,
`dF`, the rate of inbreeding ``(F_t - F_{t-1})/(1 - F_{t-1})``, and
`dFhet`, the same rate of the loss of heterozygosity `Fhet`; and the
depletion of the genetic variance within each replicate, `var_ratio` and
`genic_ratio`, the variance and the genic variance relative to the
replicate's first step (with selection, `genic_ratio` falls with allele
frequencies changing at the QTL, and `var_ratio` below it with the
negative disequilibrium of the Bulmer effect).

With a column `history` (the historical population of each replicate)
in which at least two histories have several replicates each, also
`histories`, the number of histories at that step, and for every column
`c` `c_sdh`, the standard deviation of the history means of `c`, and
`c_seh` (`c_sdh / √K` over the `K` histories): with replicates nested in
histories, `c_seh` is the standard error of the mean over histories,
which includes the variation between histories that `c_se` leaves out
(replicates of one history are not independent draws of the history).

The rows are sorted by replicate before summing, so the result is the
same, to the bit, whatever the order of the replicates in `summary`.
"""
function aggregate(df::DataFrame)
    "replicate" in names(df) && "generation" in names(df) ||
        throw(ArgumentError("the summary needs columns replicate and generation"))
    d = _with_changes(df)
    metrics = [c for c in names(d) if !(c in ("replicate", "generation", "history"))]
    gens = sort(unique(d.generation))
    out = DataFrame(generation=gens, replicates=[count(==(g), d.generation) for g in gens])
    nested = "history" in names(d) && length(unique(d.history)) >= 2 &&
             any(h -> length(unique(d.replicate[d.history .== h])) >= 2, unique(d.history))
    nested && (out[!, :histories] = [length(unique(d.history[d.generation .== g])) for g in gens])
    for c in metrics
        μ, σ, se = Float64[], Float64[], Float64[]
        for g in gens
            x = [Float64(v) for (v, gg) in zip(d[!, c], d.generation) if gg == g && !isnan(Float64(v))]
            n = length(x)
            m = n == 0 ? NaN : sum(x) / n
            s = n < 2 ? NaN : sqrt(sum(abs2, x .- m) / (n - 1))
            push!(μ, m)
            push!(σ, s)
            push!(se, n < 2 ? NaN : s / sqrt(n))
        end
        out[!, c*"_mean"] = μ
        out[!, c*"_sd"] = σ
        out[!, c*"_se"] = se
        if nested
            sdh, seh = Float64[], Float64[]
            for g in gens
                hs = sort(unique(d.history[d.generation .== g]))
                hm = Float64[]
                for h in hs
                    x = [Float64(v) for (v, gg, hh) in zip(d[!, c], d.generation, d.history)
                         if gg == g && hh == h && !isnan(Float64(v))]
                    isempty(x) || push!(hm, sum(x) / length(x))
                end
                K = length(hm)
                m = K == 0 ? NaN : sum(hm) / K
                s = K < 2 ? NaN : sqrt(sum(abs2, hm .- m) / (K - 1))
                push!(sdh, s)
                push!(seh, K < 2 ? NaN : s / sqrt(K))
            end
            out[!, c*"_sdh"] = sdh
            out[!, c*"_seh"] = seh
        end
    end
    return out
end

# A summary written by `_write_csv`, read back: every value as a Float64
# except `replicate`, `history` and `generation` (Julia's shortest
# round-trip form makes this exact).
function _read_summary(path)
    lines = readlines(path)
    hdr = split(lines[1], ',')
    cols = [Float64[] for _ in hdr]
    for l in lines[2:end]
        for (j, v) in enumerate(split(l, ','))
            push!(cols[j], parse(Float64, v))
        end
    end
    df = DataFrame([h => c for (h, c) in zip(hdr, cols)])
    df.replicate = Int.(df.replicate)
    df.generation = Int.(df.generation)
    "history" in names(df) && (df.history = Int.(df.history))
    return df
end
