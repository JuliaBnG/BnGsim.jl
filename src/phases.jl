# Phase timing of the breeding loop, for calibrating the time model of the
# resource estimate (validation/estimate). Off by default: `@_phase` then
# costs one load of a Ref. Not part of the API.

const _PHASES = Ref{Union{Nothing,Dict{Symbol,Float64}}}(nothing)

# Start (or stop, `on = false`) accumulating the wall time of each phase;
# returns the dictionary of seconds by phase.
function _phase_timing!(on::Bool = true)
    _PHASES[] = on ? Dict{Symbol,Float64}() : nothing
    return _PHASES[]
end

macro _phase(name, ex)
    quote
        local d = _PHASES[]
        if d === nothing
            $(esc(ex))
        else
            local t0 = time_ns()
            local v = $(esc(ex))
            d[$(QuoteNode(name))] = get(d, $(QuoteNode(name)), 0.0) + (time_ns() - t0) / 1e9
            v
        end
    end
end
