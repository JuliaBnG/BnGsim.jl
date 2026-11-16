# Counter-based keyed random streams (roadmap §4.2).
#
# Every random draw of a simulation is addressed by a coordinate, a
# `SimKey`. The coordinate is mapped by the Philox4x64-10 bijection
# (Salmon et al., 2011) to 256 bits, which seed a Xoshiro256++ stream used
# for the draws of that coordinate. A draw therefore depends only on its
# coordinate, not on the order in which coordinates are visited, the
# number of threads, or how work is split among them.

"""
    Domain

Codes of the `domain` field of a [`SimKey`](@ref BnGsim.SimKey): which part of a
simulation a stream belongs to. The values are part of the stream layout
and never change; new domains get new codes.
"""
module Domain
const FOUNDERS = 0x0001
const MEIOSIS = 0x0002
const PHENOTYPE = 0x0003
const MATING = 0x0004
const SURVIVAL = 0x0005
const QTL = 0x0006
const SELECTION = 0x0007
const GENOTYPING = 0x0008
end

# Version of the key → stream layout. It is mixed into the Philox key, so
# a change of layout gives new streams instead of silently reusing old ones.
const KEY_LAYOUT = 0x0001

"""
    SimKey(seed; replicate = 0, generation = 0, domain = 0, stage = 0, entity = 0, subcounter = 0)
    SimKey(base::SimKey; field = value, ...)

The coordinate of one random stream:

| Field | Type | Meaning |
| :--- | :--- | :--- |
| `seed` | `UInt64` | master seed of the scenario |
| `replicate` | `UInt32` | replicate number |
| `generation` | `UInt32` | generation or time step |
| `domain` | `UInt16` | [`Domain`](@ref) code: meiosis, phenotype, mating, … |
| `stage` | `UInt16` | a stable code of the event within the domain, never a position in a schedule |
| `entity` | `UInt32` | the animal, mating or locus the draws are for; `0` when population-scoped |
| `subcounter` | `UInt64` | a further index, e.g. sire (0) or dam (1) gamete |

The second form copies `base` and replaces the named fields, so callers
derive keys from a scenario key.

[`keyed_rng`](@ref) turns a key into a random number generator. Two keys
that differ in any field give unrelated streams.

# Layout
The Philox4x64-10 key is `(seed, replicate | KEY_LAYOUT << 32)` and the
counter is `(generation << 32 | domain << 16 | stage, entity, subcounter, 0)`.
The four 64-bit outputs are the state of a `Random.Xoshiro`. This layout
is fixed; the golden test in `test/keys.jl` guards it.
"""
struct SimKey
    seed::UInt64
    replicate::UInt32
    generation::UInt32
    domain::UInt16
    stage::UInt16
    entity::UInt32
    subcounter::UInt64
end

SimKey(seed::Integer; replicate::Integer = 0, generation::Integer = 0, domain::Integer = 0,
    stage::Integer = 0, entity::Integer = 0, subcounter::Integer = 0) =
    SimKey(UInt64(seed), UInt32(replicate), UInt32(generation), UInt16(domain),
        UInt16(stage), UInt32(entity), UInt64(subcounter))

SimKey(k::SimKey; seed::Integer = k.seed, replicate::Integer = k.replicate,
    generation::Integer = k.generation, domain::Integer = k.domain, stage::Integer = k.stage,
    entity::Integer = k.entity, subcounter::Integer = k.subcounter) =
    SimKey(seed; replicate, generation, domain, stage, entity, subcounter)

function Base.show(io::IO, k::SimKey)
    print(io, "SimKey(seed = ", k.seed, ", replicate = ", k.replicate, ", generation = ",
        k.generation, ", domain = ", k.domain, ", stage = ", k.stage, ", entity = ",
        k.entity, ", subcounter = ", k.subcounter, ")")
end

# ---- Philox4x64-10 ---------------------------------------------------------
# Constants and round function of Salmon et al. (2011), as in the Random123
# library; `test/keys.jl` checks the output against Random123.jl.

const _PHILOX_M0 = 0xD2E7470EE14C6C93
const _PHILOX_M1 = 0xCA5A826395121157
const _PHILOX_W0 = 0x9E3779B97F4A7C15
const _PHILOX_W1 = 0xBB67AE8584CAA73B

@inline function _mulhilo(a::UInt64, b::UInt64)
    p = UInt128(a) * UInt128(b)
    return (p >> 64) % UInt64, p % UInt64
end

"""
    philox4x64(key::NTuple{2,UInt64}, ctr::NTuple{4,UInt64}) -> NTuple{4,UInt64}

The Philox4x64-10 block function: a keyed bijection of 256-bit counters.
"""
@inline function philox4x64(key::NTuple{2,UInt64}, ctr::NTuple{4,UInt64})
    k0, k1 = key
    c0, c1, c2, c3 = ctr
    for r = 1:10
        if r > 1
            k0 += _PHILOX_W0
            k1 += _PHILOX_W1
        end
        hi0, lo0 = _mulhilo(_PHILOX_M0, c0)
        hi1, lo1 = _mulhilo(_PHILOX_M1, c2)
        c0, c1, c2, c3 = hi1 ⊻ c1 ⊻ k0, lo1, hi0 ⊻ c3 ⊻ k1, lo0
    end
    return (c0, c1, c2, c3)
end

_philox_key(k::SimKey) = (k.seed, UInt64(k.replicate) | (UInt64(KEY_LAYOUT) << 32))
_philox_ctr(k::SimKey) = (
    (UInt64(k.generation) << 32) | (UInt64(k.domain) << 16) | UInt64(k.stage),
    UInt64(k.entity),
    k.subcounter,
    UInt64(0),
)

"""
    keyed_rng(k::SimKey) -> Random.Xoshiro

The random number generator of coordinate `k`: a `Xoshiro` whose state is
the Philox4x64-10 image of the coordinate. It is a pure function of `k`;
use it for all the draws of that coordinate and for nothing else.
"""
function keyed_rng(k::SimKey)
    s0, s1, s2, s3 = philox4x64(_philox_key(k), _philox_ctr(k))
    # an all-zero state is invalid for Xoshiro; Philox maps exactly one
    # counter to it for a given key, so this branch is never expected
    (s0 | s1 | s2 | s3) == 0 && (s0 = 0x0000000000000001)
    return Xoshiro(s0, s1, s2, s3)
end

# ---- version-independent integer sampling ----------------------------------
# `rand(rng, a:b)` and `shuffle!` use algorithms that Julia may change
# between versions (`shuffle!` did between 1.10 and 1.11). Keyed results
# must not depend on them, so integers are drawn here from raw `UInt64`
# outputs of the stream, whose algorithm (Xoshiro256++) is fixed.

# Uniform integer in 0:(n-1), n ≥ 1, by Lemire's multiply-and-reject
# method (exact, unbiased).
@inline function _uniform_below(rng::AbstractRNG, n::UInt64)
    m = UInt128(rand(rng, UInt64)) * n
    lo = m % UInt64
    if lo < n
        t = (-n) % n                       # 2^64 mod n
        while lo < t
            m = UInt128(rand(rng, UInt64)) * n
            lo = m % UInt64
        end
    end
    return (m >> 64) % UInt64
end

"""
    keyed_int(rng, a, b)

A uniform integer in `a:b` drawn from raw 64-bit outputs of `rng`, by an
algorithm that does not change with the Julia version.
"""
@inline keyed_int(rng::AbstractRNG, a::T, b::T) where {T<:Integer} =
    a + T(_uniform_below(rng, UInt64(b - a) + one(UInt64)))

"""
    keyed_shuffle!(rng, v)

Fisher–Yates shuffle of `v` with [`keyed_int`](@ref); the same `rng` state
gives the same order on every Julia version.
"""
function keyed_shuffle!(rng::AbstractRNG, v::AbstractVector)
    for i = length(v):-1:2
        j = keyed_int(rng, 1, i)
        v[i], v[j] = v[j], v[i]
    end
    return v
end
