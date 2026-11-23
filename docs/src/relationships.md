# Relationships: A, G and realized IBD

This page compares four relationship matrices on one simulated
population, and four measures of inbreeding:

- ``A``, the numerator relationship matrix from the pedigree: the
  *expected* proportion of the genome shared identical by descent (IBD),
  twice the coancestry;
- ``I_{\mathrm{seg}}``, the *realized* IBD relationships: the proportion
  of the genome at which two animals carry the same founder haplotype,
  from the founder segments that the loop's crossover records determine
  (`founder_segments`, `irm`). It is exact: no loci are sampled;
- ``I_{\mathrm{locus}}``, the same at the marker loci only;
- ``G``, the genomic relationship matrix of VanRaden (2008), from the
  genotypes at the markers, centred on the founders' allele frequencies.

``E[I_{\mathrm{seg}} \mid \text{pedigree}] = A``: the pedigree gives the
expectation, Mendelian sampling of whole chromosome segments the
deviation. ``G`` estimates ``I_{\mathrm{seg}}`` from identity by state,
with a base given by the allele frequencies used for centring.

## The population

Model, stated in full: diploid; 3 chromosomes of 40, 30 and 30 Mb with
1 cM/Mb (1 Morgan in all); meiosis with independent assortment and
Poisson crossovers (Haldane, no interference); 100 founders with
independent random haplotypes at 6000 loci (allele frequency ½, linkage
equilibrium in expectation), unrelated and non-inbred by definition; an
additive trait (h² = 0.3, QTL every 20th locus); 6 discrete generations
of 100 offspring from 10 sires and 25 dams selected by truncation on
PBLUP EBV. The scheme records the meioses.

```@example rel
using BnGsim, BnGStructs, RelationshipMatrices, DataFrames, Random, Statistics, LinearAlgebra

sp  = GenericSpecies("tutorial", [40_000_000, 30_000_000, 30_000_000])
rng = Xoshiro(2026)
vm  = VariantMap(sort(unique(rand(rng, 1:100_000_000, 6000))), sp)
hap = Haplotype(BitMatrix(rand(rng, Bool, length(vm), 200)))
mm  = MeiosisMap(sp, vm)
t   = Trait("y"; h² = 0.3)
qtl = eqtl(hap, 1:20:length(vm), t; rng = rng)
pop = founder_population(hap, qtl, t, SimKey(1))
sc  = Scheme(; noff = 100, nsire = 10, ndam = 25, evaluation = :pblup, record_meiosis = true)
summary = simulate!(pop, sc, mm, 6)
summary[:, [:generation, :mean, :F, :Fseg, :Fhet]]
```

`F` is the mean pedigree inbreeding of each generation and `Fseg` the
mean realized inbreeding from the founder segments; `Fhet` is the loss of
heterozygosity at the loci (see [`generation_summary`](@ref)).

## Four relationship matrices

```@example rel
n    = length(pop)
ped  = pedigree(pop)
A    = RelationshipMatrices.nrm(ped.sire, ped.dam)
segs = founder_segments(pop, mm)
Iseg = irm(segs)

# founder labels at the loci, from the segments
labels = zeros(UInt32, length(vm), 2n)
for h = 1:2n, k = segs.ptr[h]:segs.ptr[h+1]-1
    r = searchsortedfirst(vm.pos, segs.start[k]):searchsortedlast(vm.pos, segs.stop[k])
    labels[r, h] .= segs.label[k]
end
Ilocus = irm(labels)

# genotypes of all animals (`gt` is padded to whole words: take the loci);
# G centred on the founders' allele frequencies
nlc = length(vm)
H = genomes(pop)
p = vec(sum(hap.gt[1:nlc, :]; dims = 2)) ./ hap.nhp
G = RelationshipMatrices.grm(H; p = p)
nothing # hide
```

The pairs of the last generation, by pedigree class:

```@example rel
last = collect(generation_ids(pop, 6))
s, d = ped.sire, ped.dam
cls(i, j) = s[i] == s[j] && d[i] == d[j] ? "full sibs" : s[i] == s[j] || d[i] == d[j] ? "half sibs" : "other"
pairs = [(i, j) for (a, i) in enumerate(last) for j in last[a+1:end]]
df = DataFrame(class = [cls(i, j) for (i, j) in pairs],
    A = [A[i, j] for (i, j) in pairs], Iseg = [Iseg[i, j] for (i, j) in pairs],
    Ilocus = [Ilocus[i, j] for (i, j) in pairs], G = [G[i, j] for (i, j) in pairs])
combine(groupby(df, :class), nrow => :pairs,
    [c => (x -> round(mean(x); digits = 3)) => Symbol(c, "_mean") for c in (:A, :Iseg, :G)]...,
    [c => (x -> round(std(x); digits = 3)) => Symbol(c, "_sd") for c in (:A, :Iseg, :G)]...)
```

Within a class the pedigree assigns nearly the same ``A`` to every pair;
the realized relationships spread around it, several times more. Under
Haldane's map the spread of half-sib and full-sib relationships has an
exact value, which the loop reproduces (`validation/ibd.jl`): it depends
on the genome's length in Morgans, here only 1 Morgan, so the spread is
large. ``G`` spreads only slightly more than ``I_{\mathrm{seg}}``: with
6000 loci, the sampling of loci and identity by state add little.

The realized relationships are also *higher* than ``A`` on average. That
is selection, not chance: over 20 replicates of this scheme the mean of
``I_{\mathrm{seg}} - A`` over the pairs of the last generation is
+0.033 (s.e. 0.008) under truncation on PBLUP EBV and −0.004 (s.e.
0.004) under random selection (`validation/ibd-selection.jl`). Selected
parents carry favourable QTL alleles, so their descendants share the
segments around them more often than Mendelian sampling implies; ``A``
assumes Mendelian sampling and ignores this. This replicate lies within
the spread of the replicates under selection (s.d. 0.034).

How the matrices relate over all pairs of the last generation:

```@example rel
off = [(i, j) for (a, i) in enumerate(last) for j in last[a+1:end]]
v(M) = [M[i, j] for (i, j) in off]
slope(x, y) = cov(x, y) / var(x)
DataFrame(comparison = ["Iseg on A", "Ilocus on Iseg", "G on Iseg"],
    correlation = round.([cor(v(A), v(Iseg)), cor(v(Iseg), v(Ilocus)), cor(v(Iseg), v(G))]; digits = 3),
    slope = round.([slope(v(A), v(Iseg)), slope(v(Iseg), v(Ilocus)), slope(v(Iseg), v(G))]; digits = 3),
    max_abs_difference = round.([maximum(abs, v(Iseg) .- v(A)), maximum(abs, v(Ilocus) .- v(Iseg)),
        maximum(abs, v(G) .- v(Iseg))]; digits = 3))
```

``I_{\mathrm{locus}}`` differs from ``I_{\mathrm{seg}}`` only by
sampling 6000 positions of the genome. ``G`` follows
``I_{\mathrm{seg}}`` with a regression slope near one because it is
centred on the founders' frequencies, the base of the pedigree; centred
on the frequencies of the current generation it would be shifted down by
the coancestry accumulated since the founders.

## Inbreeding

```@example rel
Fped = diag(A)[last] .- 1
Fseg = segment_inbreeding(segs; ids = last)
chrs = [mm.first[c]:mm.last[c] for c in eachindex(mm.first)]
Froh = froh(H.gt[1:nlc, vcat([[2i - 1, 2i] for i in last]...)], vm.pos, chrs; min_len = 4_000_000)
Fseg4 = segment_inbreeding(segs; ids = last, min_len = 4_000_000)
FG = diag(G)[last] .- 1
DataFrame(measure = ["pedigree F", "realized F (segments)", "realized F, runs ≥ 4 Mb",
        "F_ROH, runs ≥ 4 Mb", "diag(G) − 1"],
    mean = round.(mean.([Fped, Fseg, Fseg4, Froh, FG]); digits = 4),
    sd = round.(std.([Fped, Fseg, Fseg4, Froh, FG]); digits = 4),
    cor_with_Fseg = round.([cor(Fped, Fseg), 1.0, cor(Fseg4, Fseg), cor(Froh, Fseg), cor(FG, Fseg)]; digits = 3))
```

Realized inbreeding varies among animals with the same pedigree
inbreeding, for the same reason as the relationships, so the pedigree
``F`` explains little of it. It is higher than the pedigree ``F`` on
average, again because of selection (+0.015, s.e. 0.004, over the 20
replicates; +0.001, s.e. 0.003, without selection). Runs of homozygosity
of at least 4 Mb (4 cM) estimate the autozygosity in IBD runs of that
length, the recent inbreeding; they miss shorter, older autozygous
segments. With a marker every 17 kb on average, a run differs from its
autozygous segment by a few marker intervals, so here ``F_{\mathrm{ROH}}``
and the realized ``F`` in runs of at least 4 Mb agree to the precision
shown. The diagonal of ``G``, centred on the founders' frequencies,
follows the realized inbreeding closely with dense markers.

## Reproducing

Every random draw is keyed by `SimKey(1)`, so the page gives the same
numbers on any machine and number of threads (for a given set of package
versions). The variance check of realized relationships against Haldane's
map is `validation/ibd.jl`; segment IBD is validated against `A` and
locus-level IBD in RelationshipMatrices.
