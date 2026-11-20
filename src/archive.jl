# Tier 2 archives (docs/src/roadmap.md, Phase 5): the pedigree, records,
# true and estimated breeding values, QTL and genotypes of a population,
# as plain tables and PLINK 1 binary files, for analysis by other
# software; and infinitesimal (multivariate normal) breeding values on a
# pedigree, for data whose genetic model is exactly the animal model.

const _ARCHIVE_PANELS = (:none, :chip, :qtl, :all)

_csvnum(x::Real) = isnan(x) ? "" : string(x)

"""
    write_archive(dir, pop::Population, vm::VariantMap, sp::Species;
                  panel = Int32[], qtl = Int32[], t = current_step(pop)) -> Vector{String}

Write the population `pop` to directory `dir` (created if needed) and
return the files written:

- `pedigree.csv`: `id, sire, dam, sex, generation, exit` (`0` for an
  unknown parent; `sex` `M` or `F`; `generation` the birth cohort; `exit`
  the first time step at which the animal was dead or culled, empty while
  it lives);
- `records.csv`: `id` and one column per trait, the phenotype where it
  was recorded by time step `t` (expressed in the animal's sex, at or
  after the trait's age, alive at that age; [`observed`](@ref)), empty
  otherwise;
- `truth.csv`: `id`, then `tbv_<trait>` (the true genetic value) and
  `ebv_<trait>` (the last evaluation, empty if none) for each trait;
- `traits.csv`: one row per trait: `name, h2, sigma_a, mu, sex, age,
  intercept` (the constant of its genetic value);
- `qtl.csv` (if `qtl` is given, the locus indices of the QTL in `vm`):
  `locus, chr, bp`, then the additive effect of one copy of the allele
  coded 1 for each trait;
- with `panel` (locus indices in `vm`, e.g. the chip), the genotypes of
  every animal whose genome is kept at those loci, in PLINK 1 binary
  format: `genotypes.bed` (SNP-major), `genotypes.bim` (`chr`, the id
  `l<locus>`, the position in centimorgans within the chromosome from the
  map or, without one, from the species' genetic length, the position in
  bp within the chromosome, the allele coded 1 as A1 and the other as A2)
  and `genotypes.fam` (family `1`, the animal id, sire, dam, sex `1`/`2`,
  phenotype `-9`).

Numbers are written in Julia's shortest round-trip form. Animals keep the
ids of the population. In the `.bed` file the allele coded 1 is A1, the
fifth column of the `.bim`: PLINK 1.9 counts it with `--recode A
--keep-allele-order` (without that flag it makes A1 the minor allele of
the sample), plink2 with `--export A --export-allele` on that column (by
default it counts the sixth). Both were checked to read back the
genotypes exactly.
"""
function write_archive(dir::AbstractString, pop::Population, vm::VariantMap, sp::Species;
    panel::AbstractVector{<:Integer}=Int32[], qtl::AbstractVector{<:Integer}=Int32[],
    t::Integer=current_step(pop))
    mkpath(dir)
    n, k = length(pop), ntraits(pop)
    names_ = [tr.name for tr in pop.traits]
    files = String[]
    out(f) = (p = joinpath(dir, f); push!(files, p); p)
    open(out("pedigree.csv"), "w") do io
        println(io, "id,sire,dam,sex,generation,exit")
        for i = 1:n
            e = pop.exit[i] == typemax(Int32) ? "" : string(pop.exit[i])
            println(io, i, ',', pop.sire[i], ',', pop.dam[i], ',', pop.male[i] ? 'M' : 'F', ',',
                pop.generation[i], ',', e)
        end
    end
    obs = observed(pop, t)
    open(out("records.csv"), "w") do io
        println(io, join(["id"; names_], ','))
        for i = 1:n
            println(io, i, ',', join((obs[i, j] ? _csvnum(pop.phenotype[i, j]) : "" for j = 1:k), ','))
        end
    end
    open(out("truth.csv"), "w") do io
        println(io, join(["id"; ["tbv_" * m for m in names_]; ["ebv_" * m for m in names_]], ','))
        for i = 1:n
            println(io, i, ',', join(_csvnum.(pop.tbv[i, :]), ','), ',', join(_csvnum.(pop.ebv[i, :]), ','))
        end
    end
    open(out("traits.csv"), "w") do io
        println(io, "name,h2,sigma_a,mu,sex,age,intercept")
        for (tr, q) in zip(pop.traits, pop.qtl)
            println(io, tr.name, ',', tr.h², ',', tr.σₐ, ',', tr.μ, ',', tr.sex, ',', tr.age, ',', q.intercept)
        end
    end
    bp = chrpos(vm, sp)
    if !isempty(qtl)
        open(out("qtl.csv"), "w") do io
            println(io, join(["locus", "chr", "bp", names_...], ','))
            for (r, l) in enumerate(qtl)
                eff = [begin
                    j = findfirst(==(Int32(l)), q.loci)
                    j === nothing ? 0.0 : q.additive[j]
                end for q in pop.qtl]
                println(io, l, ',', vm.chr[l], ',', bp[l], ',', join(string.(eff), ','))
            end
        end
    end
    if !isempty(panel)
        ids = retained(pop)
        _write_plink(out, pop, ids, genomes(pop, ids), panel, vm, sp, bp)
    end
    return files
end

# PLINK 1 binary genotypes of animals `ids` (haplotypes `hap`, columns
# 2i - 1 and 2i) at loci `panel`. The allele coded 1 is A1: two copies are
# 00, one is 10, none is 11 (01 would be missing), four animals per byte
# from the low bits (PLINK 1.9 file format).
function _write_plink(out, pop, ids, hap::Haplotype, panel, vm, sp, bp)
    nid = length(ids)
    open(out("genotypes.fam"), "w") do io
        for i in ids
            println(io, "1 ", i, ' ', pop.sire[i], ' ', pop.dam[i], ' ', pop.male[i] ? 1 : 2, " -9")
        end
    end
    ends = cbp(sp)
    open(out("genotypes.bim"), "w") do io
        for l in panel
            c = vm.chr[l]
            len = c == 1 ? ends[1] : ends[c] - ends[c-1]
            cm = isempty(vm.morgan) ? 100 * sp.morgan[c] * bp[l] / len : 100 * vm.morgan[l]
            println(io, c, "\tl", l, '\t', cm, '\t', bp[l], '\t', vm.alt[l], '\t', vm.ref[l])
        end
    end
    code = (0x03, 0x02, 0x00)          # dosage 0, 1, 2 of the allele coded 1
    nb = cld(nid, 4)
    row = zeros(UInt8, nb)
    open(out("genotypes.bed"), "w") do io
        write(io, UInt8[0x6c, 0x1b, 0x01])
        gt = hap.gt
        for l in panel
            fill!(row, 0x00)
            @inbounds for k = 1:nid
                d = gt[l, 2k-1] + gt[l, 2k]
                row[(k-1)>>2+1] |= code[d+1] << (2 * ((k - 1) & 3))
            end
            write(io, row)
        end
    end
    return nothing
end

"""
    infinitesimal(sire, dam, G0, R0, key::SimKey; mean = zeros(size(G0, 1)))
        -> (u::Matrix{Float64}, y::Matrix{Float64})

Breeding values and phenotypes of `k = size(G0, 1)` traits under the
infinitesimal model, exactly ``u \\sim N(0, G_0 \\otimes A)``, for the
animals of a pedigree (`sire`, `dam`: `0` for unknown; parents before
offspring). Each animal gets ``u_i = \\tfrac12 (u_s + u_d) + m_i`` with the
Mendelian sampling term ``m_i \\sim N(0, d_i G_0)``,
``d_i = 1 - \\tfrac14 (1 + F_s) - \\tfrac14 (1 + F_d)`` (a term only for each
known parent; `F` the inbreeding coefficients), and
``y_i = \\mu + u_i + e_i`` with ``e_i \\sim N(0, R_0)``. `u` and `y` are
`n × k`.

The draws of animal `i` come from `SimKey(key; domain = Domain.INFINITESIMAL,
entity = i)`, so the values of an animal do not depend on the others or
on the order of computation. These are data whose genetic model is the
animal model itself (an "Architecture A" for benchmarks), unlike the
finite-locus genetic values of the breeding loop.
"""
function infinitesimal(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer}, G0::AbstractMatrix,
    R0::AbstractMatrix, key::SimKey; mean::AbstractVector=zeros(size(G0, 1)))
    n, k = length(sire), size(G0, 1)
    length(dam) == n || throw(ArgumentError("sire and dam must have the same length"))
    size(G0) == (k, k) && size(R0) == (k, k) && length(mean) == k ||
        throw(ArgumentError("G0 and R0 must be $k × $k and mean of length $k"))
    for i = 1:n
        (0 <= sire[i] < i && 0 <= dam[i] < i) ||
            throw(ArgumentError("animal $i: parents must be known animals listed before it, or 0"))
    end
    Lg = cholesky(Symmetric(Matrix{Float64}(G0))).L
    Le = cholesky(Symmetric(Matrix{Float64}(R0))).L
    F = nrm_diag(sire, dam) .- 1
    u = zeros(n, k)
    y = zeros(n, k)
    for i = 1:n
        rng = keyed_rng(SimKey(key; domain=Domain.INFINITESIMAL, entity=i))
        s, d = sire[i], dam[i]
        δ = 1.0 - (s > 0 ? (1 + F[s]) / 4 : 0.0) - (d > 0 ? (1 + F[d]) / 4 : 0.0)
        z = randn(rng, k)
        ε = randn(rng, k)
        for j = 1:k
            pa = (s > 0 ? u[s, j] / 2 : 0.0) + (d > 0 ? u[d, j] / 2 : 0.0)
            u[i, j] = pa + sqrt(δ) * dot(view(Lg, j, :), z)
            y[i, j] = mean[j] + u[i, j] + dot(view(Le, j, :), ε)
        end
    end
    return u, y
end
