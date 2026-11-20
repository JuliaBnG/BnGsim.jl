"""
    provenance(io::IO = stdout; seed = nothing, extra = Dict{String,Any}())

Write a run header recording what is needed to reproduce a simulation: the
date, host, CPU, Julia version, thread count, BLAS threads, the seed, and for
each JuliaBnG package its version, git commit and whether its working tree
was dirty. Lines are prefixed with `# ` so the header can sit on top of a
plain-text log or table.

Results from `FisherWright.fisher_wright` are reproducible for a given seed
only at a fixed thread count, which is why the thread count is recorded.
"""
function provenance(io::IO = stdout; seed = nothing, extra = Dict{String,Any}())
    println(io, "# date:     ", Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))
    println(io, "# host:     ", gethostname())
    println(io, "# cpu:      ", Sys.cpu_info()[1].model, " × ", Sys.CPU_THREADS)
    println(io, "# memory:   ", round(Sys.total_memory() / 2^30; digits = 1), " GiB")
    println(io, "# julia:    ", VERSION)
    println(io, "# threads:  ", Threads.nthreads(), " (BLAS ", BLAS.get_num_threads(), ")")
    seed === nothing || println(io, "# seed:     ", seed)
    for m in _provenance_modules()
        println(io, "# package:  ", rpad(string(m), 21), _pkg_state(m))
    end
    for (k, v) in extra
        println(io, "# ", rpad(string(k) * ":", 10), v)
    end
    flush(io)
    return nothing
end

function _provenance_modules()
    return Module[BnGStructs, FisherWrightAdapter.PACKAGE, RelationshipsAdapter.PACKAGE, BreedingAdapter.PACKAGE,
        BnGsim]
end

function _pkg_state(m::Module)
    dir = pkgdir(m)
    ver = pkgversion(m)
    commit = "no-git"
    if dir !== nothing && isdir(joinpath(dir, ".git"))
        try
            commit = readchomp(`git -C $dir rev-parse --short=10 HEAD`)
            isempty(readchomp(`git -C $dir status --porcelain --untracked-files=no`)) ||
                (commit *= "-dirty")
        catch
        end
    end
    return "v$ver $commit"
end
