# Keyed streams: the Philox block function, the key layout and the
# independence of streams. The golden values below fix the layout; a
# change to any of them is a change of every seeded result and must be
# reviewed as such, never updated as a plain test fix.
using Random123: Random123

@testset "SimKey and keyed streams" begin
    z, f = UInt64(0), typemax(UInt64)
    # Known-answer vectors of Philox4x64-10 (Random123 kat_vectors)
    @test philox4x64((z, z), (z, z, z, z)) ==
          (0x16554d9eca36314c, 0xdb20fe9d672d0fdc, 0xd7e772cee186176b, 0x7e68b68aec7ba23b)
    @test philox4x64((f, f), (f, f, f, f)) ==
          (0x87b092c3013fe90b, 0x438c3c67be8d0224, 0x9cc7d7c69cd777b6, 0xa09caebf594f0ba0)
    @test philox4x64((0xa4093822299f31d0, 0x082efa98ec4e6c89),
        (0x243f6a8885a308d3, 0x13198a2e03707344, z, z)) ==
          (0x4876440167f74106, 0xf268c71745f14876, 0x2aeeaf1df0e3941c, 0xa0b2ec3a9d392db5)
    # and agrees with Random123.jl on arbitrary inputs
    rng = Xoshiro(1)
    for _ = 1:500
        key = (rand(rng, UInt64), rand(rng, UInt64))
        ctr = (rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64))
        @test philox4x64(key, ctr) == Random123.philox(key, ctr, Val(10))
    end

    # Golden stream: the layout of a key and the Xoshiro draws it starts
    k = SimKey(2026; generation = 3, domain = Domain.MEIOSIS, stage = 1, entity = 7, subcounter = 1)
    r = keyed_rng(k)
    @test [rand(r, UInt64) for _ = 1:3] == [0x5b5a5fd02115167e, 0x68e4218105a35751, 0x52cbe1a2dc87bcb3]
    @test [rand(keyed_rng(k), UInt64) for _ = 1:2] == fill(0x5b5a5fd02115167e, 2)   # pure function

    # Every field addresses a different stream
    first_draw(k) = rand(keyed_rng(k), UInt64)
    variants = [SimKey(k; seed = 2027), SimKey(k; replicate = 1), SimKey(k; generation = 4),
        SimKey(k; domain = Domain.PHENOTYPE), SimKey(k; stage = 2), SimKey(k; entity = 8),
        SimKey(k; subcounter = 0)]
    draws = first_draw.([k; variants])
    @test allunique(draws)
    @test SimKey(k; entity = 8).entity == 8 && SimKey(k; entity = 8).seed == k.seed
    @test occursin("entity = 7", sprint(show, k))

    # Streams of neighbouring coordinates are not correlated: the first
    # uniform of 20 000 consecutive entities, mean ½ and lag-1 correlation 0
    # (4 standard errors)
    u = [rand(keyed_rng(SimKey(k; entity = e))) for e = 1:20_000]
    @test abs(mean(u) - 0.5) < 4 * sqrt(1 / 12 / 20_000)
    @test abs(cor(u[1:end-1], u[2:end])) < 4 / sqrt(20_000)
end
