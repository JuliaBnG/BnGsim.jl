# Optimum contribution selection with a target rate of inbreeding per
# generation, for discrete or overlapping generations (BnGsim 0.8).
#
# Three ways of turning a target ΔF per generation into a constraint at
# each time step t (selection of the parents of cohort t + 1):
#
# - `:cohort`: the mean kinship of the newborn cohort, c'Θc, may rise by
#   ΔF per generation, i.e. by 1 − (1 − ΔF)^(1/L) per step, L the
#   generation interval estimated from the parental ages used so far.
# - `:age_year` (Meuwissen & Sonesson 1998; Wellmann 2019, optiSel): the
#   mean kinship of the whole population at t + 1, w'Θw with
#   w = W₀c + v a mixture of age classes weighted by their remaining
#   expected contributions, may rise by the same per-step amount.
# - `:age_generation`: as `:age_year`, but the per-step amount uses the
#   generation interval of the contributions being chosen, L(c), instead
#   of the estimate (a per-generation constraint; Sonesson et al. 2000).
#
# Θ = A/2 is the pedigree kinship (coancestry) matrix.

"""
    ocs_contributions(u, Θ, b, K, male) -> (c, feasible, value)

Contributions `c` (each sex summing to ½, all ``\\ge 0``) maximising
``u'c`` subject to ``f(c) = c'\\Theta c + 2b'c \\le K``: optimum
contribution selection (Meuwissen 1997) with the linear term ``b`` that
the fixed contributions of older age classes add.

The problem is convex. It is solved exactly by homotopy in the weight
``\\nu`` of the gain: ``c(\\nu)`` minimises ``f(c) - 2\\nu u'c`` under the
sex constraints and ``c \\ge 0``. For a fixed set of zero contributions
``c(\\nu) = y + \\nu x`` is linear. The path starts at ``\\nu = \\infty``,
the best male and the best female (ties by position), and ``\\nu`` is
lowered from one breakpoint to the next (a candidate entering when its
multiplier reaches zero, or a contribution reaching zero) until
``f(c(\\nu)) = K``. The free set stays as small as the solution, so a
breakpoint costs a factorisation of the order of the number of parents.
If the path reaches ``\\nu = 0`` first, ``c(0)`` minimises ``f``, the bound
cannot be met, and it is returned with `feasible = false`. `value` is
``f(c)`` at the solution.
"""
function ocs_contributions(u::AbstractVector, Θ::AbstractMatrix, b::AbstractVector, K::Real,
    male::AbstractVector{Bool})
    n = length(u)
    length(b) == length(male) == n == size(Θ, 1) == size(Θ, 2) ||
        throw(DimensionMismatch("u, Θ, b and male must agree"))
    any(male) && !all(male) || throw(ArgumentError("candidates of both sexes are needed"))
    T = Matrix{Float64}(Θ)
    uv, bv = Vector{Float64}(u), Vector{Float64}(b)
    Q = Float64.(hcat(male, .!male))
    s = [0.5, 0.5]
    f(c) = dot(c, T * c) + 2dot(bv, c)
    zero_set = trues(n)
    for sex in (true, false)              # the best of each sex, ties by position
        ids = findall(==(sex), male)
        zero_set[ids[argmax(uv[ids])]] = false
    end
    c = zeros(n)
    ν = Inf
    for _ = 1:10n
        F = findall(.!zero_set)
        y, x, e0, e1, TF = _kkt_path(T, bv, uv, Q, s, F, zero_set)
        α = dot(x, TF * x)
        β = 2dot(x, TF * y) + 2dot(bv[F], x)
        γ = dot(y, TF * y) + 2dot(bv[F], y)
        fν(t) = isinf(t) ? (α > 0 ? Inf : γ + (β == 0 ? 0.0 : sign(β) * Inf)) : α * t^2 + β * t + γ
        if fν(ν) <= K                      # the bound holds already: c(ν) is optimal
            c[F] = isinf(ν) ? y : y .+ ν .* x
            return max.(c, 0), true, f(max.(c, 0))
        end
        # next breakpoint below ν
        νnext, kind, who = 0.0, :zero, 0
        for (k, i) in enumerate(F)         # a contribution reaches zero
            if x[k] > 0
                νi = -y[k] / x[k]
                νi < ν && νi > νnext && ((νnext, kind, who) = (νi, :leave, i))
            end
        end
        Z = findall(zero_set)
        for (k, j) in enumerate(Z)         # a candidate's multiplier reaches zero
            if e1[k] > 0
                νj = -e0[k] / e1[k]
                νj < ν && νj > νnext && ((νnext, kind, who) = (νj, :enter, j))
            end
        end
        # the bound on this piece: the largest root of f(t) = K in (νnext, ν)
        if fν(νnext) <= K
            t = if α > 0
                d = β^2 - 4α * (γ - K)
                (-β + sqrt(max(d, 0.0))) / 2α
            else
                (K - γ) / β
            end
            t = clamp(t, νnext, isinf(ν) ? t : ν)
            fill!(c, 0)
            c[F] = y .+ t .* x
            c .= max.(c, 0)
            return c, true, f(c)
        end
        if kind == :zero                   # ν = 0 reached: the minimum of f, above the bound
            fill!(c, 0)
            c[F] = y
            c .= max.(c, 0)
            return c, false, f(c)
        end
        ν = νnext
        zero_set[who] = kind == :leave
    end
    error("OCS homotopy did not terminate")
end

# Solution of the KKT system on the free set F with the contributions of
# the zero set fixed at 0, as a linear function of ν:
# c_F = y + ν x, and the multipliers of the zero set η = e0 + ν e1.
# Stationarity: 2Θc + 2b − 2νu + Qλ − η = 0.
function _kkt_path(T, b, u, Q, s, F, zero_set)
    nf = length(F)
    TF = T[F, F]
    QF = Q[F, :]
    M = [2TF QF; QF' zeros(2, 2)]
    Fm = lu(M)
    z0 = Fm \ [-2b[F]; s]
    z1 = Fm \ [2u[F]; zeros(2)]
    y, x = z0[1:nf], z1[1:nf]
    λ0, λ1 = z0[nf+1:end], z1[nf+1:end]
    Z = findall(zero_set)
    TZF = T[Z, F]
    # η = 2Θc + 2b − 2νu + Qλ on the zero set
    e0 = 2TZF * y .+ 2b[Z] .+ Q[Z, :] * λ0
    e1 = 2TZF * x .- 2u[Z] .+ Q[Z, :] * λ1
    return y, x, e0, e1, TF
end

# Proportions of each sex's parental contributions by age at selection,
# from the parents of cohorts 1..t (age 0 only when there is no history),
# for ages 0..amax.
function _age_usage(pop::Population, t::Integer, male::Bool, amax::Integer)
    u = zeros(amax + 1)
    for g = 1:t, i in generation_ids(pop, g)
        p = male ? pop.sire[i] : pop.dam[i]
        u[clamp(g - 1 - pop.generation[p], 0, amax)+1] += 1
    end
    sum(u) == 0 && (u[1] = 1)
    return u ./ sum(u)
end

# Weights of the population at t + 1 for the age-weighted constraints:
# the newborn cohort's share W₀ (its genes come from the contributions c)
# and the fixed shares v of the animals alive at t that are still of
# breeding age at t + 1. Age class k of a sex has weight
# ½ Σ_{b ≥ k} u_b / L, u the usage by age and L = Σ (b + 1) u_b, i.e. its
# share of the sex's remaining expected contributions; within a class the
# weight is shared equally (optiSel). Also returns the generation
# interval estimate ½(L_m + L_f).
function _age_weights(pop::Population, sc::Scheme, t::Integer)
    W0 = 0.0
    v = Dict{Int,Float64}()
    Ls = Float64[]
    for (male, amax) in ((true, sc.max_sire_age), (false, sc.max_dam_age))
        u = _age_usage(pop, t, male, amax)
        L = sum(((0:amax) .+ 1) .* u)
        push!(Ls, L)
        W = [0.5 * sum(u[k+1:end]) / L for k = 0:amax]
        W0 += W[1]
        ids = filter(i -> pop.male[i] == male && pop.generation[i] <= t < pop.exit[i] &&
                          t - pop.generation[i] + 1 <= amax, 1:length(pop))
        for a = 0:amax-1
            cls = filter(i -> t - pop.generation[i] == a, ids)
            for i in cls
                v[i] = W[a+2] / length(cls)
            end
        end
    end
    return W0, v, sum(Ls) / 2
end

# Realised generation interval of step s: mean parental age at the birth
# of cohort s + 1, over both sexes.
_realised_L(pop::Population, s::Integer) =
    mean(i -> (s + 1 - pop.generation[pop.sire[i]] + s + 1 - pop.generation[pop.dam[i]]) / 2,
        generation_ids(pop, s + 1))

# The bound of the constrained kinship at step t + 1: from the founders'
# mean kinship 1/(2n₀), it rises by 1 − (1 − ΔF)^(1/L_s) at every step s,
# L_s the estimated (`:cohort`, `:age_year`) or realised
# (`:age_generation`) generation interval; `Lt` is the interval used at
# the current step t.
function _ocs_bound(pop::Population, sc::Scheme, t::Integer, Lt::Real)
    B = θ0(pop)
    for s = 0:t
        L = s == t ? Lt :
            sc.ocs == :age_generation ? _realised_L(pop, s) : last(_age_weights(pop, sc, s))
        B = 1 - (1 - B) * (1 - sc.ΔF)^(1 / L)
    end
    return B
end

# Contributions of the candidates `cand` (index values `u`) under the
# scheme's dynamic OCS constraint. Returns (c, info) with info the bound,
# the value of the constrained kinship reached, the generation interval
# used and feasibility.
function _dynamic_ocs(pop::Population, sc::Scheme, t::Integer, cand, u)
    W0, v, Lest = sc.ocs == :cohort ? (1.0, Dict{Int,Float64}(), last(_age_weights(pop, sc, t))) :
                  _age_weights(pop, sc, t)
    others = sort!(collect(setdiff(keys(v), cand)))
    ids = [cand; others]
    Θ = nrm(pop.sire, pop.dam, ids) ./ 2
    nc = length(cand)
    vv = [get(v, i, 0.0) for i in ids]
    Θcc = Θ[1:nc, 1:nc]
    b = (Θ[1:nc, :] * vv) ./ W0
    vΘv = dot(vv, Θ * vv)
    male = Vector{Bool}(pop.male[cand])
    age = t .- pop.generation[cand]
    L = Lest
    local c, feasible, val, B
    for _ = 1:(sc.ocs == :age_generation ? 30 : 1)
        B = _ocs_bound(pop, sc, t, L)
        K = (B - vΘv) / W0^2
        c, feasible, val = ocs_contributions(u, Θcc, b, K, male)
        sc.ocs == :age_generation || break
        Lc = sum(c .* (age .+ 1))
        abs(Lc - L) < 1e-6 && break
        L = (L + Lc) / 2                  # damped fixed-point iteration on L(c)
    end
    θ = W0^2 * val + vΘv                  # the constrained kinship reached
    return c, (ocs_bound = B, ocs_kinship = θ, ocs_L = L, ocs_feasible = feasible)
end
