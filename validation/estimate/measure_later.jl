# One run of a scenario (all its replicates, as `bngsim run` does) in this
# process: prints the file, the peak resident set size (MiB), the expected
# peak and bound of `resource_estimate` (MiB), and the run time (s) against
# its estimate, tab-separated. For the terms added after the calibration.
using BnGsim
f = ARGS[1]
s = read_scenario(f)
e = resource_estimate(s)
t = @elapsed run_scenario(s)
println(join([basename(f), round(Sys.maxrss() / 2^20, digits=1), round(e.peak_rss_bytes / 2^20, digits=1),
    round(e.peak_rss_bound_bytes / 2^20, digits=1), round(t, digits=1), e.seconds], "\t"))
