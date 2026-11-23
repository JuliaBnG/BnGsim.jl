module BnGsim

using BnGStructs
using Dates
using LinearAlgebra
using Mmap
using Serialization
using JSON
using SHA: sha256
using PrecompileTools: @setup_workload, @compile_workload
using DataFrames
using DataStructures: BinaryMinHeap
using Distributions
using Random
using Statistics
using StatsBase: sample

# Upstream packages are imported only through these adapters
# (docs/src/upstream-contract.md).
include("adapters/fisherwright.jl")
include("adapters/relationships.jl")
include("adapters/breeding.jl")
using .FisherWrightAdapter
using .RelationshipsAdapter
import .RelationshipsAdapter: founder_segments
using .BreedingAdapter

include("keys.jl")
include("loci.jl")
include("meiosis.jl")
include("pedigree-view.jl")
include("loop.jl")
include("mating.jl")
include("ocs.jl")
include("traits.jl")
include("historical.jl")
include("founder-summary.jl")
include("scenario.jl")
include("estimate.jl")
include("runner.jl")
include("histories.jl")
include("archive.jl")
include("aggregate.jl")
include("cli.jl")
include("legacy.jl")
using .Legacy: vanraden_g   # for founder_summary

export crossovers, sample_loci
export eqtl, scale_qtl, trait_summary, RandomSign, SpikeSlab
export SimKey, Domain, keyed_rng, philox4x64, keyed_int, keyed_shuffle!
export MeiosisMap, MeiosisRecord, make_progeny!, ncrossovers
export Population, founder_population, Scheme, evaluate!, select_parents, mating_plan, transport, step!, simulate!
export generation_summary, ngenerations, generation_ids, current_step, recorded, observed, ntraits, index_weights,
    alive, genotyped, ocs_contributions,
    meioses, founder_segments, segment_inbreeding, irm
export PedigreeView, pedigree_view, row_to_id, id_to_row, evaluation_view, pedigree, genomes, retained
export tbv, ebv, phenotypes, records, traits, trait_qtl, covariances, parity
export save_checkpoint, load_checkpoint
export founder_haplotypes
export ld_decay, founder_summary
export Scenario, ScenarioError, read_scenario, parse_scenario, validate_scenario, write_scenario, scenario_json,
    scenario_dict, results_dict, results_sha256, run_scenario, run_replicate
export cli, ExitCode
export resource_estimate, ScenarioSize
export write_archive, infinitesimal
export aggregate, write_histories, check_histories

include("precompile.jl")

end # module BnGsim
