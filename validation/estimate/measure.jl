# One run of a scenario in this process: prints the file, the peak resident
# set size (MiB), the mean number of derived alleles per historical
# haplotype, the segregating loci among the founders, the animals, ne, nt,
# nid, the evaluation, the chip size and the run time (s), tab-separated.
using BnGsim, Statistics
f = ARGS[1]
s = read_scenario(f)
sp = BnGsim._species_of(s)
t = @elapsed begin
    h = BnGsim._history(s, sp, 0)
    r = run_replicate(s, 1; history=h)
end
mbar = mean(length, h.active_haplotypes)
println(join([f, round(Sys.maxrss() / 2^20, digits=1), round(mbar, digits=1), length(r.map), length(r.population),
    s.founders.ne, s.founders.nt, s.founders.nid, s.scheme.evaluation, s.chip === nothing ? 0 : s.chip.nlc,
    round(t, digits=1)], "\t"))
