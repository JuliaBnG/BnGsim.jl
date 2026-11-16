# Comparison of three OCS constraints with overlapping generations.
#
#   julia -t 4 --project=<env with BnGsim, BnGStructs, DataFrames> validation/ocs-overlap.jl \
#         [replicates] [steps] [first replicate] [csv] [birth|age1]
#
# Replicates first … first + replicates − 1; with `csv`, the per-replicate
# metrics are also written there, so runs on several hosts can be pooled
# (`validation/ocs-overlap-pool.jl`).
#
# Question: which way of turning a target ΔF per generation into a
# per-step constraint holds the realised rate of inbreeding at the target
# when parents of several ages compete, and at what genetic gain per year?
#
# Methods (BnGsim `Scheme(; selection = :ocs, ocs, ΔF)`):
# - `:cohort`: mean kinship of the newborn cohort;
# - `:age_year`: mean kinship of the population at the next step, age
#   classes weighted by their remaining expected contributions
#   (Meuwissen & Sonesson 1998; Wellmann 2019, optiSel), per-step increase
#   from the estimated generation interval;
# - `:age_generation`: the same, per-step increase from the generation
#   interval of the chosen contributions (Sonesson et al. 2000).
# Truncation selection (`nsire`, `ndam`, same ages) is the reference
# without inbreeding control.
#
# Scenario (stated in full):
# - genome: 3 chromosomes of 50 Mb, 1 M each; 3000 loci at random
#   positions; founders: 200 individuals with independent random
#   haplotypes (linkage equilibrium in expectation);
# - trait: additive, 300 QTL, normal effects, h² = 0.3, σₐ = 1 in the
#   founders, recorded in both sexes at birth (`birth`) or at age 1
#   (`age1`: young candidates are evaluated through relatives only, so
#   older parents compete; step 0 is random mating, as nothing is
#   recorded yet);
# - scheme: 100 offspring per step, parents of age 0 to 3 (both sexes),
#   PBLUP with the true variance components, ΔF = 0.01 per generation
#   (Nₑ = 50); truncation: 8 sires and 25 dams;
# - metrics over the second half of the run (steps ≥ nstep/2), per
#   replicate: ΔF per year from the regression of ln(1 − F̄) of each birth
#   cohort on the cohort, converted to ΔF per generation with the
#   realised generation interval L, (1 − ΔF_g) = (1 − ΔF_y)^L; the same from
#   the mean coancestry of each cohort (Δθ); genetic gain per year (slope
#   of the cohort means of the genetic values, in founder σₐ); L; parents
#   used; steps where the bound could not be met.
# Every method runs on the same founders, QTL and keys within a replicate
# (paired comparison).
using BnGsim, BnGStructs, Random, Statistics, Printf, DataFrames

nrep = parse(Int, get(ARGS, 1, "50"))
nstep = parse(Int, get(ARGS, 2, "20"))
first_rep = parse(Int, get(ARGS, 3, "1"))
csv = get(ARGS, 4, "")
scenario = get(ARGS, 5, "birth")
scenario in ("birth", "age1") || error("scenario must be birth or age1")
sp = GenericSpecies("val", [50_000_000, 50_000_000, 50_000_000])
ΔF = 0.01
methods = (:cohort, :age_year, :age_generation, :truncation)

slope(x, y) = (xm = mean(x); sum((x .- xm) .* (y .- mean(y))) / sum((x .- xm) .^ 2))

function metrics(df, h)
    rows = df[h:end-1, :]                       # selection steps of the second half
    coh = df[h:end, :]                          # cohorts born in the second half
    L = mean(rows.L)
    by = slope(coh.generation, log.(1 .- coh.F))
    bθ = slope(coh.generation, log.(1 .- coh.coancestry))
    dFg = 1 - exp(by)^L
    dθg = 1 - exp(bθ)^L
    gain = slope(coh.generation, coh.mean)
    infeas = hasproperty(rows, :ocs_feasible) ? count(!, rows.ocs_feasible) : 0
    (ΔF_gen = dFg, Δθ_gen = dθg, gain_year = gain, L = L, parents = mean(rows.nsire .+ rows.ndam),
        infeasible = infeas)
end

res = Dict(m => NamedTuple[] for m in methods)
t0 = time()
for r = first_rep:first_rep+nrep-1
    rng = Xoshiro(50_000 + r)
    vm = VariantMap(sort(unique(rand(rng, 1:150_000_000, 3000))), sp)
    hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 400)))
    mm = MeiosisMap(sp, vm)
    t = Trait("y"; h² = 0.3, σₐ = 1.0, age = scenario == "age1" ? 1 : 0)
    qtl = eqtl(hap, sort(randperm(rng, length(vm))[1:300]), t; rng = rng)
    for m in methods
        pop = founder_population(hap, qtl, t, SimKey(2026; replicate = r))
        sc = m == :truncation ?
             Scheme(; noff = 100, nsire = 8, ndam = 25, evaluation = :pblup, max_sire_age = 3, max_dam_age = 3) :
             Scheme(; noff = 100, nsire = 8, ndam = 25, evaluation = :pblup, selection = :ocs, ocs = m, ΔF = ΔF,
                 max_sire_age = 3, max_dam_age = 3)
        scenario == "age1" && step!(pop, Scheme(; noff = 100, nsire = 8, ndam = 25, evaluation = :tbv,
            selection = :random, max_sire_age = 3, max_dam_age = 3), mm)
        df = simulate!(pop, sc, mm, nstep)
        push!(res[m], metrics(df, nstep ÷ 2 + 1))
    end
    if !isempty(csv)
        open(csv, "a") do io
            for m in methods
                x = res[m][end]
                println(io, join((r, m, x.ΔF_gen, x.Δθ_gen, x.gain_year, x.L, x.parents, x.infeasible), ","))
            end
        end
    end
    @printf(stderr, "replicate %d (%d to %d) (%.0f s): %s\n", r, first_rep, first_rep + nrep - 1, time() - t0,
        join(("$m ΔF $(round(res[m][end].ΔF_gen; digits = 4)) gain $(round(res[m][end].gain_year; digits = 3))"
              for m in methods), "; "))
end

println("# BnGsim ", pkgversion(BnGsim), ", Julia ", VERSION, ", ", Threads.nthreads(), " threads, ",
    nrep, " replicates (from ", first_rep, "), ", nstep, " steps, scenario ", scenario, ", target ΔF = ", ΔF,
    " per generation, ",
    round(time() - t0; digits = 0), " s")
se(x) = std(x) / sqrt(length(x))
@printf("%-15s %18s %18s %18s %14s %10s %10s\n", "method", "ΔF/generation", "Δθ/generation", "gain/year",
    "L", "parents", "infeasible")
for m in methods
    x = res[m]
    f(k) = [getproperty(v, k) for v in x]
    @printf("%-15s %9.4f (%.4f) %9.4f (%.4f) %9.4f (%.4f) %6.2f (%.2f) %10.1f %10.2f\n", m,
        mean(f(:ΔF_gen)), se(f(:ΔF_gen)), mean(f(:Δθ_gen)), se(f(:Δθ_gen)), mean(f(:gain_year)),
        se(f(:gain_year)), mean(f(:L)), se(f(:L)), mean(f(:parents)), mean(f(:infeasible)))
end
# paired differences of the OCS methods from :age_year
for m in (:cohort, :age_generation)
    d(k) = [getproperty(a, k) - getproperty(b, k) for (a, b) in zip(res[m], res[:age_year])]
    @printf("%-15s − age_year: ΔF/gen %+.4f (%.4f), gain/year %+.4f (%.4f), L %+.2f (%.2f)\n", m,
        mean(d(:ΔF_gen)), se(d(:ΔF_gen)), mean(d(:gain_year)), se(d(:gain_year)), mean(d(:L)), se(d(:L)))
end
