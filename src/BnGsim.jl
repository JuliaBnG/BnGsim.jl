module BnGsim

using BnGStructs
using Dates
using LinearAlgebra
using Mmap
using Serialization
using FisherWright
using JSON
using SHA: sha256
using PrecompileTools: @setup_workload, @compile_workload
using RelationshipMatrices
using Breeding: Breeding, MME, mme, RandomEffect, solve_mme, incidence, TM1997, feasible, rocs, pev_candidates, reml,
    MultiTrait, mt_stack, mt_blocks, mt_rinv
using DataFrames
using Distributions
using Random
using Statistics
using StatsBase: sample

include("keys.jl")
include("reproduction.jl")
include("meiosis.jl")
include("pedigree-view.jl")
include("loop.jl")
include("ocs.jl")
include("traits.jl")
include("base-populations.jl")
include("provenance.jl")
include("historical.jl")
include("founder-summary.jl")
include("scenario.jl")
include("estimate.jl")
include("cli.jl")

export sum_map, sumMap, crossovers, gamete!, gene_drop, gene_drop!
export random_mate, rnd_ped, expand, sample_id, sampleID, sample_loci, sampleLoci
export norm_mv!, eqtl, trait_summary, tbv!, phenotype!, dummy_traits
export FounderPopulation, sample_founder_population
export provenance
export SimKey, Domain, keyed_rng, philox4x64, keyed_int, keyed_shuffle!
export MeiosisMap, MeiosisRecord, make_progeny!, ncrossovers
export Population, founder_population, Scheme, evaluate!, select_parents, mating_plan, step!, simulate!
export generation_summary, ngenerations, generation_ids, current_step, recorded, observed, ntraits, index_weights,
    alive, genotyped, ocs_contributions,
    meioses, founder_segments, segment_inbreeding, irm
export PedigreeView, pedigree_view, row_to_id, id_to_row, evaluation_view, pedigree, genomes, retained
export save_checkpoint, load_checkpoint
export HistoricalSpec, historical_population, continue_population, founder_haplotypes
export read_packed_haplotypes, write_xy
export ld_decay, vanraden_g, founder_summary
export Scenario, ScenarioError, read_scenario, parse_scenario, validate_scenario, write_scenario, scenario_json,
    scenario_dict, run_scenario, run_replicate
export cli, ExitCode
export resource_estimate, ScenarioSize

include("precompile.jl")

end # module BnGsim
