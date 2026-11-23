# The FisherWright side of BnGsim (docs/src/upstream-contract.md). BnGsim
# imports FisherWright only here; the rest of the package uses these names.
module FisherWrightAdapter

using FisherWright: FisherWright, fisher_wright, muts2bitarray, FisherWrightResult, bottleneck, growth,
    RandomMating, trait_loci, TraitLoci, TraitFitness, Stabilizing, Directional
using Random: AbstractRNG

# The upstream package (for its version).
const PACKAGE = FisherWright

export fisher_wright, muts2bitarray, FisherWrightResult, census, dense_haplotypes, bottleneck, growth,
    RandomMating, trait_loci, TraitLoci, TraitFitness, Stabilizing, Directional, write_history_file,
    read_history_file

"""
    census(res::FisherWrightResult) -> Vector{Int}

The census sizes of the recorded generations, `N₀, N₁, …`: `res.records.census`
from FisherWright 0.10, `res.census` before.
"""
census(res::FisherWrightResult) =
    hasfield(FisherWrightResult, :records) ? getfield(res, :records).census : getfield(res, :census)

"""
    dense_haplotypes(muts, ends; rng) -> (BitMatrix, DataFrame)

`muts2bitarray(muts, ends)`: the loci segregating among the haplotypes
`muts` (rows, genome-wide positions) by haplotype (columns), and the map.
The placeholder nucleotides of the map come from `rng` where FisherWright
accepts one (0.10 and later), from the global generator before.
"""
dense_haplotypes(muts, ends; rng::AbstractRNG) =
    applicable_rng() ? muts2bitarray(muts, ends; rng=rng) : muts2bitarray(muts, ends)

"""
    write_history_file(path, res::FisherWrightResult)
    read_history_file(path) -> FisherWrightResult

FisherWright's binary result format (`write_result`, `read_result`, from
FisherWright 0.10): every haplotype, the selected alleles, the records and
the pedigree. An `ArgumentError` with an earlier FisherWright.
"""
function write_history_file(path, res::FisherWrightResult)
    isdefined(FisherWright, :write_result) ||
        throw(ArgumentError("storing a history needs FisherWright 0.10 or later"))
    return FisherWright.write_result(path, res)
end
function read_history_file(path)
    isdefined(FisherWright, :read_result) ||
        throw(ArgumentError("reading a stored history needs FisherWright 0.10 or later"))
    return FisherWright.read_result(path)
end

# Whether `muts2bitarray` takes an `rng` (FisherWright 0.10 and later).
applicable_rng() = pkgversion(FisherWright) >= v"0.10.0"

end
