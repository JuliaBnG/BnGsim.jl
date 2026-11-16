# Meiosis benchmark: the keyed kernel `make_progeny!` against the v0.2
# `gene_drop!` on the same genome, parents and mating plan.
#
#   julia -t 8 --project=<env with BnGsim, BnGStructs> bench/meiosis.jl [nloci] [noffspring]
#
# Genome: cattle autosomes (29, 2.49 Gb, 1 cM/Mb). Loci: `nloci` positions
# drawn uniformly. Parents: 1000 individuals with random haplotypes.
# Times are the minimum of 3 runs after a warm-up; one machine, no
# replication across machines.
using BnGsim, BnGStructs, Random, Printf

nlc = parse(Int, get(ARGS, 1, "50000"))
noff = parse(Int, get(ARGS, 2, "10000"))
sp = Cattle()
rng = Xoshiro(1)
pos = sort(unique(rand(rng, 1:Int(total_bp(sp)), nlc)))
vm = VariantMap(pos, sp)
nlc = length(vm)
par = Haplotype(BitMatrix(rand(rng, Bool, nlc, 2000)))
plan = hcat(rand(rng, 1:1000, noff), rand(rng, 1:1000, noff))
mm = MeiosisMap(sp, vm)
key = SimKey(1; generation = 1)
off = Haplotype(nlc, 2noff)
best(f) = (f(); minimum(@elapsed(f()) for _ = 1:3))
told = best(() -> gene_drop!(off, par, plan, sp, vm; rng = Xoshiro(2)))
tnew = best(() -> make_progeny!(off, par, plan, mm, key))
trec = best(() -> make_progeny!(off, par, plan, mm, key; record = true))
tmap = @elapsed MeiosisMap(sp, vm)
@printf("# %d loci, %d offspring, %d threads, Julia %s, %s\n", nlc, noff, Threads.nthreads(),
    VERSION, Sys.cpu_info()[1].model)
@printf("gene_drop! (v0.2)          %8.3f s\n", told)
@printf("make_progeny!              %8.3f s  (%.1f× faster)\n", tnew, told / tnew)
@printf("make_progeny!, record      %8.3f s\n", trec)
@printf("MeiosisMap, once           %8.3f s\n", tmap)
