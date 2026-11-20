# Resource estimate of a scenario (docs/src/roadmap.md §5, invariant I7):
# the expected peak resident memory and a conservative bound, the disk
# written, and a rough run time, from the scenario alone.

"""
    ScenarioSize

The quantities of a scenario that its memory and time grow with, computed
before it runs:

- `alleles`: the expected number of derived alleles a haplotype of the
  historical population carries, ``μL \\min(t, 4N)`` for `N = ne` over
  `t` generations (`nt` plus the decline): ``μL`` new ones per
  generation, until fixation balances them at the neutral equilibrium
  ``4NμL``;
- `loci`: the expected number of loci segregating among the `2nid` founder
  haplotypes, ``4NμL \\ln(1 + 2n_{id} t / 4N)``, the total branch length of
  their genealogy truncated at the start of the run (the population starts
  without variation);
- `animals`: founders and offspring of the run, `nid + steps × noff`;
- `evaluated`: the largest pedigree evaluated, `nid + (steps - 1) × noff`;
- `genotyped`: the animals in the dense relationship matrices of `:gblup`
  (all evaluated) and `:ssgblup` (the last `genotyped_cohorts` cohorts, or
  all); 0 otherwise;
- `cohort`: the largest cohort (the per-step summary forms its dense
  relationship matrix).
"""
struct ScenarioSize
    alleles::Float64
    loci::Float64
    animals::Int
    evaluated::Int
    genotyped::Int
    cohort::Int
end

function ScenarioSize(s::Scenario)
    f, sc = s.founders, s.scheme
    μL = f.mr * 1e-8 * sum(s.species.chromosome)
    t = f.nt + length(f.decline)
    N = f.ne
    alleles = μL * min(t, 4N)
    loci = 4N * μL * log1p(2f.nid * t / 4N)
    animals = f.nid + s.steps * sc.noff
    evaluated = f.nid + (s.steps - 1) * sc.noff
    genotyped = if sc.evaluation === :gblup || (sc.evaluation === :ssgblup && sc.genotyped_cohorts == 0)
        evaluated
    elseif sc.evaluation === :ssgblup
        k = sc.genotyped_cohorts
        min(k, s.steps) * sc.noff + (k >= s.steps ? f.nid : 0)
    else
        0
    end
    return ScenarioSize(alleles, loci, animals, evaluated, genotyped, max(f.nid, sc.noff))
end

# Coefficients of the memory model, fitted by relative least squares (non-
# negative) to the mean peak resident set size of 16 calibration scenarios,
# two runs each, on the calibration host (docs/src/estimate.md: AMD
# Threadripper 3970X, Linux x86_64, Julia 1.13.4, 8 threads, BnGsim 0.14).
# Bytes per unit of each term. Terms for the pedigree (per animal evaluated)
# and for the dense cohort matrix of the summaries were tried and fitted to
# zero: the genome term absorbs them.
const _MEM = (
    baseline = 5.767e8,    # bytes, the process with BnGsim loaded and compiled
    history = 5.416,       # per byte of one copy of the historical haplotypes, 8N × alleles
    genomes = 4.150,       # per byte of the packed genomes of all animals, animals × 2 × loci / 8
    dense_gblup = 13.28,   # per byte of one dense genotyped × genotyped matrix, :gblup
    dense_ssgblup = 19.01, # the same, :ssgblup
)
# The bound is the expected peak times this factor: the largest ratio of
# measured to expected peak over every calibration and held-out run (1.153,
# 44 runs), with 10 % added and rounded up.
const _MEM_BOUND_FACTOR = 1.3
# Run time, fitted to the same runs: seconds per allele carried through one
# meiosis in the historical run; per evaluated animal and step; per cube of
# the genotyped set and step; and a fixed cost of the dense evaluations
# (mostly compilation, as they are not in the precompile workload).
const _TIME = (history = 7.63e-10, pedigree = 1.139e-5, dense = 8.446e-12, dense_fixed = 4.48)

# The memory terms of a scenario, in bytes.
function _memory_terms(s::Scenario, z::ScenarioSize)
    ev = s.scheme.evaluation
    dense = 8.0 * z.genotyped^2
    return (
        baseline = _MEM.baseline,
        history = _MEM.history * 8.0 * s.founders.ne * z.alleles,
        genomes = _MEM.genomes * z.animals * 2 * z.loci / 8,
        dense = (ev === :gblup ? _MEM.dense_gblup : ev === :ssgblup ? _MEM.dense_ssgblup : 0.0) * dense,
    )
end

"""
    resource_estimate(s::Scenario) -> NamedTuple

The resources scenario `s` needs, from the scenario alone:

- `peak_rss_bytes`: the expected peak resident memory of a run of one
  replicate at a time, the sum of `memory_terms`: the process itself, the
  historical population, the packed genomes of all animals (and with them
  the pedigree and records), and the dense relationship matrices of
  `:gblup` and `:ssgblup`;
- `peak_rss_bound_bytes`: a conservative bound for admission: the
  expectation times a factor that no calibration or held-out run
  exceeded, with a margin;
- `disk_bytes`: the files `bngsim run` writes;
- `seconds`: a rough run time of all replicates (historical population and
  breeding loop), on the calibration host within a factor of about two;
- `size`: the [`ScenarioSize`](@ref) the terms are computed from.

The coefficients are calibrated on one host and thread count (see the
documentation page "Resource estimate"); memory and time differ on
others, and the garbage collector makes the peak of a given run vary by
about 10 %. Replicates run one after the other, so memory does not grow
with their number.
"""
function resource_estimate(s::Scenario)
    z = ScenarioSize(s)
    terms = _memory_terms(s, z)
    peak = sum(values(terms))
    f = s.founders
    t = f.nt + length(f.decline)
    μL = f.mr * 1e-8 * sum(s.species.chromosome)
    # alleles carried through meiosis over the historical run
    carried = 2f.ne * μL * sum(min(g, 4f.ne) for g = 1:t)
    nh = f.shared_history ? 1 : s.replicates
    loop = s.replicates * s.steps * (_TIME.pedigree * z.evaluated + _TIME.dense * float(z.genotyped)^3) +
           (s.scheme.evaluation in (:gblup, :ssgblup) ? _TIME.dense_fixed : 0.0)
    seconds = nh * _TIME.history * carried + loop
    rows = s.replicates * (s.steps + 1)
    disk = 4096 + 2length(scenario_json(s)) + 400 * rows
    return (; peak_rss_bytes=round(Int, peak), peak_rss_bound_bytes=round(Int, _MEM_BOUND_FACTOR * peak),
        memory_terms=map(x -> round(Int, x), terms), disk_bytes=disk, seconds=round(seconds; digits=1), size=z)
end
