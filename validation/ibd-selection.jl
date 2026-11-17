# Realized IBD against the pedigree expectation, with and without
# selection: the population of docs/src/relationships.md (100 founders,
# 3 chromosomes of 1 M in all, 6 generations of 100 offspring from 10
# sires and 25 dams), 20 replicates per selection rule (keys 1..20).
#
#   julia -t 4 --project=docs validation/ibd-selection.jl
#
# Reports, for the last generation, the mean over pairs of I_seg − A and
# the mean of F_seg − F, with the s.e. over replicates. Without selection
# E[I_seg | pedigree] = A; selection on EBV makes selected relatives share
# the segments around favourable QTL alleles more often than Mendelian
# sampling implies, which A ignores.
using BnGsim, BnGStructs, RelationshipMatrices, Random, Statistics, LinearAlgebra, Printf
sp  = GenericSpecies("tutorial", [40_000_000, 30_000_000, 30_000_000])
rng = Xoshiro(2026)
vm  = VariantMap(sort(unique(rand(rng, 1:100_000_000, 6000))), sp)
hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 200)))
mm  = MeiosisMap(sp, vm)
t   = Trait("y"; h² = 0.3)
qtl = eqtl(hap, 1:20:length(vm), t; rng = rng)
for sel in (:truncation, :random)
    dI, dF = Float64[], Float64[]
    for r = 1:20
        pop = founder_population(hap, qtl, t, SimKey(1; replicate = r))
        simulate!(pop, Scheme(; noff = 100, nsire = 10, ndam = 25, evaluation = :pblup, selection = sel,
            record_meiosis = true), mm, 6)
        last = collect(generation_ids(pop, 6))
        A = RelationshipMatrices.nrm(pop.sire, pop.dam, last)
        segs = founder_segments(pop, mm)
        I = irm(segs; ids = last)
        off = [(a, b) for a in eachindex(last) for b in a+1:length(last)]
        push!(dI, mean(I[a, b] - A[a, b] for (a, b) in off))
        push!(dF, mean(segment_inbreeding(segs; ids = last)) - mean(diag(A) .- 1))
    end
    @printf("%-10s mean(I − A) off-diagonal %.4f (s.e. %.4f, sd over replicates %.4f); mean(Fseg − F) %.4f (s.e. %.4f)\n",
        sel, mean(dI), std(dI) / sqrt(20), std(dI), mean(dF), std(dF) / sqrt(20))
end
