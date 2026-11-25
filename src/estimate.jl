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
# Run time (0.37): seconds per unit of each feature of `_time_features`,
# each phase fitted to its own timings (validation/estimate/fit_phases.jl
# on the phase timings of run_phases.sh, phases.tsv: 56 scenarios, two
# processes each, on an AMD Ryzen 9 5950X, Linux x86_64, Julia 1.13.4, 8
# threads); `history` per allele carried through one meiosis of the
# historical run and `history_cache` the same times `_cache_factor`,
# `rest` per time step and `rest_animal` per animal and step, `compile`
# the cold-minus-warm cost of the first replicate by what the scenario
# uses (`_compile_flags`). Terms fitted to zero are kept so that the
# fit can be rerun. The cubic terms are identified by scenarios up to
# 20 000 genotyped animals with chips of 1000 to 20 000 loci (0.37.1).
const _TIME = (
    history = 2.223e-10, history_cache = 9.748e-11,
    setup = 2.53e-9, setup_fixed = 0.013,
    meiosis = 1.741e-11,
    summary_genome = 7.183e-11, summary_pedigree = 2.268e-9,
    equations = 1.514e-6, equations2 = 0.0, dense = 2.213e-13, dense2 = 4.061e-8,
    reml_equations = 1.742e-5, reml_dense = 8.44e-11, records = 0.0,
    grm = 2.923e-12, inverse = 1.655e-12, grm_ss = 8.514e-12, inverse_ss = 0.0,
    mating = 2.364e-8,
    reliability = 8.329e-9,
    rest = 0.003428, rest_animal = 4.635e-7,
    compile = (0.414, 0.0, 0.0, 0.0, 0.0, 0.683, 0.0, 0.369),
)
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
# The dense matrices of several traits and of REML (0.37), which the 0.15
# calibration (one trait, no REML) did not cover: bytes per byte of the
# genotyped × genotyped matrix per k² - 1 for k traits, and per k² with
# REML; fitted to the peaks of the multi-trait and REML scenarios of the
# time calibration (validation/estimate/fit_phases.jl, Ryzen 9 5950X).
# Fitted on the multi-trait and REML scenarios of each evaluation (GBLUP:
# mt2_gblup, mt2_gblup_n10k, reml2_gblup; ssGBLUP: mt2_ssgblup and
# rep2_ssgblup, so exactly); REML with ssGBLUP uses the GBLUP value. Held out in check/: three traits, REML with one and two
# traits, ssGBLUP with REML, repeated records (mean/expected 0.94 to 1.13).
const _MEM_TRAITS = (gblup = 9.83, ssgblup = 15.1)
const _MEM_REML = 9.85

# Bytes per byte of one dense genotyped × genotyped matrix.
function _dense_factor(s::Scenario)
    ev, k = s.scheme.evaluation, length(s.traits)
    ev in (:gblup, :ssgblup) || return 0.0
    c, t = ev === :gblup ? (_MEM.dense_gblup, _MEM_TRAITS.gblup) : (_MEM.dense_ssgblup, _MEM_TRAITS.ssgblup)
    return c + t * (k^2 - 1) + (s.scheme.reml ? _MEM_REML * k^2 : 0.0)
end

# The memory terms of a scenario, in bytes.
function _memory_terms(s::Scenario, z::ScenarioSize)
    ev = s.scheme.evaluation
    dense = 8.0 * z.genotyped^2
    k = length(s.traits)
    return (
        baseline = _MEM.baseline,
        history = 8.0 * z.census * z.alleles * (_MEM.history + _MEM_HISTORY_STORED * (z.histories - 1)),
        genomes = _MEM.genomes * z.animals * 2 * z.loci / 8,
        dense = _dense_factor(s) * dense,
        records = _MEM_RECORD * z.records,
        reliability = s.scheme.report_reliability && ev in (:pblup, :gblup, :ssgblup) ?
            _MEM_RELIABILITY * k * z.equations : 0.0,
    )
end

# The alleles carried through meiosis over one historical run, generation
# by generation with the census trajectory.
function _carried(s::Scenario, z::ScenarioSize = ScenarioSize(s))
    f = s.founders
    μL = f.mr * 1e-8 * sum(s.species.chromosome)
    sizes = [f.popsize === nothing ? fill(f.ne, f.nt) : _trajectory(f.popsize, f.ne).(1:f.nt); f.decline]
    return sum(2 * sizes[g] * μL * min(g, 4z.census) for g in eachindex(sizes); init=0.0)
end

# The cost per carried allele grows as the haplotypes of one generation,
# about 2 × census × alleles × 8 bytes, outgrow the processor caches: the
# factor log2 of that size over 16 MB, zero below (validation/estimate).
_cache_factor(z::ScenarioSize) = max(0.0, log2(16.0 * z.census * z.alleles / 2^24))

# The sizes the time model of the breeding loop grows with (0.37), per
# replicate, summed over the time steps g = 0, …, steps - 1 at which the
# loop evaluates and selects (the pedigree then holds nid + g·noff
# animals). Shared by `resource_estimate` and validation/estimate/fit_phases.jl.
function _time_features(s::Scenario, z::ScenarioSize = ScenarioSize(s))
    f, sc = s.founders, s.scheme
    nid, noff, T, k = f.nid, sc.noff, s.steps, length(s.traits)
    ev = sc.evaluation
    dense = ev in (:gblup, :ssgblup)
    blup = ev in (:pblup, :gblup, :ssgblup)
    pe = any(t -> t.records > 1 && t.repeatability > t.h2, s.traits) ? 2 : 1
    size(c) = c == 0 ? nid : noff                          # cohort c
    function ngen(g)                                       # genotyped at step g
        dense || return 0
        (ev === :gblup || sc.genotyped_cohorts == 0) && return nid + g * noff
        return sum(size(c) for c = max(0, g - sc.genotyped_cohorts + 1):g)
    end
    L = z.loci
    steps = 0:T-1
    animals(g) = nid + g * noff
    return (
        setup = 2.0 * nid * L,                             # founder haplotypes × loci
        meiosis = T * 2.0 * noff * L,                      # offspring haplotypes × loci
        summary_genome = sum(2.0 * (size(g) + nid) * L for g = 0:T),   # counts over cohort and founders
        summary_pedigree = sum(float(animals(g)) * size(g) for g = 0:T), # relationships of the cohort
        equations = blup ? sum(float(animals(g)) * k * pe for g in steps) : 0.0,
        equations2 = blup ? sum(float(animals(g)) * (k * pe)^2 for g in steps) : 0.0,  # k × k blocks
        dense = sum((float(k) * ngen(g))^3 for g in steps),
        dense2 = sum((float(k) * ngen(g))^2 for g in steps),
        reml_equations = blup && sc.reml ? sum(float(animals(g)) * k * pe for g in steps) : 0.0,
        reml_dense = sc.reml ? sum((float(k) * ngen(g))^3 for g in steps) : 0.0,
        records = z.records * (T + 1) / 2,                 # records made by each step, summed
        grm = s.chip === nothing ? 0.0 : sum(float(ngen(g))^2 * s.chip.nlc for g in steps),
        inverse = sum(float(ngen(g))^3 for g in steps),
        grm_ss = ev === :ssgblup && s.chip !== nothing ? sum(float(ngen(g))^2 * s.chip.nlc for g in steps) : 0.0,
        inverse_ss = ev === :ssgblup ? sum(float(ngen(g))^3 for g in steps) : 0.0,   # H: A22 and its inverse
        mating = sc.mating === :min_coancestry ? sum(float(animals(g)) * (sc.nsire + sc.ndam) for g in steps) : 0.0,
        reliability = sc.report_reliability && blup ? T * float(noff) * k * z.equations : 0.0,
        animals = sum(float(animals(g)) for g in steps),  # records, phenotypes, selection, culling
    )
end

# What the first replicate compiles: always; BLUP; dense G; ssGBLUP;
# several traits; repeated records; REML; minimum-coancestry mating.
function _compile_flags(s::Scenario)
    sc = s.scheme
    ev = sc.evaluation
    return (1.0, float(ev in (:pblup, :gblup, :ssgblup)), float(ev in (:gblup, :ssgblup)), float(ev === :ssgblup),
        float(length(s.traits) > 1), float(any(t -> t.records > 1, s.traits)), float(sc.reml),
        float(sc.mating === :min_coancestry))
end

# The run time by part, in seconds: the historical runs, and over all
# replicates the founder setup, meiosis, summaries, evaluation (with the
# G or H inverse), mating, reliabilities, the rest of the loop, and the
# compilation of the first replicate.
function _time_terms(s::Scenario, z::ScenarioSize = ScenarioSize(s))
    f, c = s.founders, _TIME
    x = _time_features(s, z)
    nh = f.shared_history ? min(f.histories, s.replicates) : s.replicates
    isempty(f.history_dir) || (nh = 0)                     # read, not simulated
    R = s.replicates
    return (
        history = nh * _carried(s, z) * (c.history + c.history_cache * _cache_factor(z)),
        setup = R * (c.setup * x.setup + c.setup_fixed),
        meiosis = R * c.meiosis * x.meiosis,
        summary = R * (c.summary_genome * x.summary_genome + c.summary_pedigree * x.summary_pedigree),
        evaluation = R * (c.equations * x.equations + c.equations2 * x.equations2 + c.dense * x.dense +
                          c.dense2 * x.dense2 + c.reml_equations * x.reml_equations + c.reml_dense * x.reml_dense +
                          c.records * x.records + c.grm * x.grm + c.inverse * x.inverse + c.grm_ss * x.grm_ss +
                          c.inverse_ss * x.inverse_ss),
        mating = R * c.mating * x.mating,
        reliability = R * c.reliability * x.reliability,
        rest = R * (c.rest * s.steps + c.rest_animal * x.animals),
        compile = sum(c.compile .* _compile_flags(s)),
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
  breeding loop), on the calibration host of the time model mostly within
  about ±30 % (up to 2.2 times low on one scenario; calibrated up to
  20 000 genotyped animals);
- `time_terms`: `seconds` by part: the historical runs, and over all
  replicates the founder setup, meiosis, summaries, evaluation (with the
  G or H inverse), mating, reliabilities, the rest of the loop, and the
  compilation of the first replicate;
- `size`: the [`ScenarioSize`](@ref) the terms are computed from.

The coefficients are calibrated at one thread count, memory and time on
two hosts (see the documentation page "Resource estimate"); memory and
time differ on others, and the garbage collector makes the peak of a given run vary by
about 10 %. Replicates run one after the other, so memory does not grow
with their number.
"""
function resource_estimate(s::Scenario)
    z = ScenarioSize(s)
    terms = _memory_terms(s, z)
    peak = sum(values(terms))
    tt = _time_terms(s, z)
    seconds = sum(values(tt))
    rows = s.replicates * (s.steps + 1)
    disk = 4096 + 2length(scenario_json(s)) + 400 * rows
    return (; peak_rss_bytes=round(Int, peak), peak_rss_bound_bytes=round(Int, _MEM_BOUND_FACTOR * peak),
        memory_terms=map(x -> round(Int, x), terms), disk_bytes=disk, seconds=round(seconds; digits=1),
        time_terms=map(x -> round(x; digits=1), tt), size=z)
end
