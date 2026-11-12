module BnGsim

using BnGStructs
using Dates
using LinearAlgebra
using Mmap
using Serialization
using FisherWright
using RelationshipMatrices
using DataFrames
using Distributions
using Random
using Statistics
using StatsBase: sample

include("reproduction.jl")
include("traits.jl")
include("base-populations.jl")
include("provenance.jl")
include("historical.jl")
include("founder-summary.jl")

export sum_map, sumMap, crossovers, gamete!, gene_drop, gene_drop!
export random_mate, rnd_ped, expand, sample_id, sampleID, sample_loci, sampleLoci
export norm_mv!, eqtl, tbv!, phenotype!, dummy_traits
export FounderPopulation, sample_founder_population
export provenance
export HistoricalSpec, historical_population, continue_population, founder_haplotypes
export read_packed_haplotypes, write_xy
export ld_decay, vanraden_g, founder_summary

end # module BnGsim
