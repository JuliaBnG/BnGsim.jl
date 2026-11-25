# The run time of a scenario by phase (0.37 time model). In this process,
# replicates 1 (cold: includes compilation) and 2 (warm) of the scenario,
# each with its own history; per replicate one tab-separated row:
#   file, replicate, history, setup, evaluate, kinv, reliability, select,
#   mating, meiosis, summary, loop (all of simulate!), in seconds,
# then the sizes the terms use: steps, nid, noff, traits, chip, reml,
# mating, evaluation, the fields of ScenarioSize, and the peak resident
# set size of the process so far (MiB).
# `kinv` (G or H inverse) is part of `evaluate`.
using BnGsim
const PH = (:setup, :evaluate, :kinv, :reliability, :select, :mating, :meiosis, :summary)
f = ARGS[1]
s = read_scenario(f)
sp = BnGsim._species_of(s)
z = BnGsim.ScenarioSize(s)
sc = s.scheme
for r = 1:2
    th = @elapsed h = BnGsim._history(s, sp, r)
    d = BnGsim._phase_timing!()
    tr = @elapsed run_replicate(s, r; history = h)
    BnGsim._phase_timing!(false)
    loop = tr - d[:setup]
    sizes = Any[s.steps, s.founders.nid, sc.noff, length(s.traits), s.chip === nothing ? 0 : s.chip.nlc,
        sc.reml, sc.mating, sc.evaluation, (getfield(z, k) for k in fieldnames(typeof(z)))...]
    println(join(Any[basename(f), r, round(th, digits = 3), (round(get(d, p, 0.0), digits = 3) for p in PH)...,
        round(loop, digits = 3), sizes..., round(Sys.maxrss() / 2^20, digits = 1)], "\t"))
    flush(stdout)
end
