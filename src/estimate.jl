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
  relationship matrix);
- `histories`: the historical populations held in memory at once (`K` of
  `founders.histories` when shared; 1 otherwise), and `census`, the
  largest census of the history (FisherWright sizes its buffers by it);
- `records`: the repeated records made over the run (record mode: about
  one per animal, trait and step of its recorded life);
- `equations`: the mixed-model equations of the evaluations (animal and,
  with repeated records, permanent-environment effects).

With a census trajectory (`founders.popsize`), `N` in `alleles` and
`loci` is the harmonic mean of the census over the historical
generations.
"""
struct ScenarioSize
    alleles::Float64
    loci::Float64
    animals::Int
    evaluated::Int
    genotyped::Int
    cohort::Int
    histories::Int
    census::Int
    records::Float64
    equations::Int
end

function ScenarioSize(s::Scenario)
    f, sc = s.founders, s.scheme
    μL = f.mr * 1e-8 * sum(s.species.chromosome)
    t = f.nt + length(f.decline)
    sizes = [f.popsize === nothing ? fill(f.ne, f.nt) : _trajectory(f.popsize, f.ne).(1:f.nt); f.decline]
    N = isempty(sizes) ? f.ne : round(Int, length(sizes) / sum(1 ./ sizes))       # harmonic mean
    census = maximum([f.ne; sizes])
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
    k = length(s.traits)
    life = max(sc.max_sire_age, sc.max_dam_age) + 1
    recmode = any(t -> t.records > 1, s.traits)
    records = recmode ? animals * sum(min(t.records, life, s.steps) for t in s.traits) : 0.0
    pe = any(t -> t.records > 1 && t.repeatability > t.h2, s.traits)
    equations = evaluated * k * (pe ? 2 : 1)
    histories = f.shared_history ? f.histories : 1
    return ScenarioSize(alleles, loci, animals, evaluated, genotyped, max(f.nid, sc.noff), histories, census,
        records, equations)
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
const _TIME = (history = 7.63e-10, pedigree = 1.139e-5, dense = 8.446e-12, dense_fixed = 4.48,
    reliability = 1.2e-8)   # per candidate × equation and step (report_reliability; 0.22, PBLUP)
# Terms added after the calibration, from the data structures rather than a
# fit (docs/src/estimate.md, "Later terms"): bytes per repeated record
# (four vectors of 4, 4, 4 and 8 bytes, with up to twice that from their
# growth) and per element of the work array of `report_reliability`.
const _MEM_RECORD = 40.0
# Bytes per byte of each further history held in memory (the results only,
# without FisherWright's working buffers that `_MEM.history` includes);
# fitted to one scenario with four histories (validation/estimate/later).
const _MEM_HISTORY_STORED = 2.4
const _MEM_RELIABILITY = 8.0 * 64

# The memory terms of a scenario, in bytes.
function _memory_terms(s::Scenario, z::ScenarioSize)
    ev = s.scheme.evaluation
    dense = 8.0 * z.genotyped^2
    k = length(s.traits)
    return (
        baseline = _MEM.baseline,
        history = 8.0 * z.census * z.alleles * (_MEM.history + _MEM_HISTORY_STORED * (z.histories - 1)),
        genomes = _MEM.genomes * z.animals * 2 * z.loci / 8,
        dense = (ev === :gblup ? _MEM.dense_gblup : ev === :ssgblup ? _MEM.dense_ssgblup : 0.0) * dense,
        records = _MEM_RECORD * z.records,
        reliability = s.scheme.report_reliability && ev in (:pblup, :gblup, :ssgblup) ?
            _MEM_RELIABILITY * k * z.equations : 0.0,
    )
end

"""
    resource_estimate(s::Scenario) -> NamedTuple

The resources scenario `s` needs, from the scenario alone:

- `peak_rss_bytes`: the expected peak resident memory of a run of one
  replicate at a time, the sum of `memory_terms`: the process itself, the
  historical populations held at once, the packed genomes of all animals
  (and with them the pedigree and single records), the dense relationship
  matrices of `:gblup` and `:ssgblup`, the repeated records, and the work
  array of `report_reliability`;
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
    # alleles carried through meiosis over the historical run, generation by
    # generation with the census trajectory
    sizes = [f.popsize === nothing ? fill(f.ne, f.nt) : _trajectory(f.popsize, f.ne).(1:f.nt); f.decline]
    carried = sum(2 * sizes[g] * μL * min(g, 4z.census) for g in eachindex(sizes); init=0.0)
    nh = f.shared_history ? min(f.histories, s.replicates) : s.replicates
    isempty(f.history_dir) || (nh = 0)                     # read, not simulated
    rel = s.scheme.report_reliability ? _TIME.reliability * s.scheme.noff * length(s.traits) * z.equations : 0.0
    loop = s.replicates * s.steps * (_TIME.pedigree * z.evaluated + _TIME.dense * float(z.genotyped)^3 + rel) +
           (s.scheme.evaluation in (:gblup, :ssgblup) ? _TIME.dense_fixed : 0.0)
    seconds = nh * _TIME.history * carried + loop
    rows = s.replicates * (s.steps + 1)
    disk = 4096 + 2length(scenario_json(s)) + 400 * rows
    return (; peak_rss_bytes=round(Int, peak), peak_rss_bound_bytes=round(Int, _MEM_BOUND_FACTOR * peak),
        memory_terms=map(x -> round(Int, x), terms), disk_bytes=disk, seconds=round(seconds; digits=1), size=z)
end
