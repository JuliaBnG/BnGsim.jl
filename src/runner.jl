# A run of a scenario that can stop at a time-step boundary and resume
# (docs/src/roadmap.md §4.8, `bngsim resume`). The state lives in
# `DIR/.bngsim-run/` while the run is unfinished; a run continued from it is
# identical to an uninterrupted one, because every draw is keyed by
# coordinate (SimKey) and nothing else carries state between steps.

const _STATE_DIR = ".bngsim-run"
const _STATE_MAGIC = UInt8['B', 'N', 'G', 'S', 'R', 'U', 'N', '1']

_state_dir(dir) = joinpath(dir, _STATE_DIR)
_versions() = (bngsim = string(pkgversion(@__MODULE__)), julia = string(VERSION))

# Write `x` with a magic and the versions, atomically.
function _save_state(path, x)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        write(io, _STATE_MAGIC)
        serialize(io, _versions())
        serialize(io, x)
    end
    mv(tmp, path; force=true)
    return path
end

# Read what `_save_state` wrote; refuses another BnGsim or Julia version,
# whose serialization format may differ.
function _load_state(path)
    open(path) do io
        read(io, 8) == _STATE_MAGIC || throw(ArgumentError("$path is not BnGsim run state"))
        v = deserialize(io)
        v == _versions() || throw(ArgumentError(
            "the run was started by BnGsim $(v.bngsim) on Julia $(v.julia); this is BnGsim " *
            "$(_versions().bngsim) on Julia $(_versions().julia), and run state does not carry over"))
        return deserialize(io)
    end
end

"""
    _execute(dir, s, reps; checkpoint_every = 0, should_stop = (r, g) -> false)

Run replicates `reps` of `s` with its state in `dir/.bngsim-run`, or
continue a run whose state is there. After step `g` of replicate `r`, if
`should_stop(r, g)` the state is saved and `(:cancelled, r, g)` returned;
with `checkpoint_every = k > 0` it is also saved every `k` steps, so that a
killed process can be resumed from there. Returns `(:completed, summary)`
with the summaries of all replicates otherwise.

State files: `run.jls` (scenario JSON, replicates, options, the sessions
of the run), `history.jls` (the shared historical population, once saved),
`done.jls` (the summaries of finished replicates) and `current.jls` (the
replicate in progress: its population, map, loci and the rows so far).
"""
function _execute(dir, s::Scenario, reps; checkpoint_every::Integer=0, should_stop=(r, g) -> false)
    st = _state_dir(dir)
    p_hist, p_done, p_cur = joinpath(st, "history.jls"), joinpath(st, "done.jls"), joinpath(st, "current.jls")
    sp = _species_of(s)
    hist = if !s.founders.shared_history
        nothing
    elseif isfile(p_hist)
        _load_state(p_hist)
    else
        _history(s, sp, 0)
    end
    done = isfile(p_done) ? _load_state(p_done) : Pair{Int,DataFrame}[]
    cur = isfile(p_cur) ? _load_state(p_cur) : nothing
    function save!(r, g, rep, rows)
        hist === nothing || isfile(p_hist) || _save_state(p_hist, hist)
        _save_state(p_done, done)
        _save_state(p_cur, (replicate=r, step=g, population=rep.population, map=rep.map, qtl=rep.qtl,
            chip=rep.chip, rows=rows))
    end
    for r in reps
        any(d -> d.first == r, done) && continue
        if cur !== nothing && cur.replicate == r
            sc = Scheme(; s.scheme..., chip=cur.chip)
            rep = (; population=cur.population, scheme=sc, meiosis=MeiosisMap(sp, cur.map), map=cur.map,
                qtl=cur.qtl, chip=cur.chip)
            g0, rows = cur.step, cur.rows
            cur = nothing
        else
            rep = _setup_replicate(s, r, hist)
            g0, rows = 0, Any[]
        end
        for g = g0+1:s.steps
            push!(rows, step!(rep.population, rep.scheme, rep.meiosis))
            if should_stop(r, g)
                save!(r, g, rep, rows)
                return (:cancelled, r, g)
            end
            checkpoint_every > 0 && g % checkpoint_every == 0 && g < s.steps && save!(r, g, rep, rows)
        end
        df = DataFrame([rows; _final_row(rep.population, rep.scheme, rep.meiosis)])
        insertcols!(df, 1, :replicate => fill(r, nrow(df)))
        push!(done, r => df)
        if checkpoint_every > 0 || isfile(p_cur)
            _save_state(p_done, done)
            isfile(p_cur) && rm(p_cur)
        end
    end
    order = Dict(r => i for (i, r) in enumerate(reps))
    sort!(done; by=d -> order[d.first])
    return (:completed, reduce(vcat, last.(done)))
end
