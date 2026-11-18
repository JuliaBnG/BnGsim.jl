"""
    PedigreeView(sire, dam, ids)

A pedigree restricted to the animals `ids`, in the form
RelationshipMatrices expects: dense rows `1:length(ids)`, parents before
offspring, and `0` for an unknown parent. Row `r` is the animal
`v.ids[r]`; its parents are `v.sire[r]` and `v.dam[r]`, as rows, `0` when
the parent is unknown or not in the view.

`sire` and `dam` are the full pedigree by id (`0` unknown), with every
parent's id below its offspring's, as in a [`Population`](@ref); `ids`
need not be sorted or unique and are stored in ascending order, which is
then a topological order. The view is the boundary between the domain
ids, which are never reused, and the row indices of the relationship
matrices and mixed-model equations built on it: a matrix built from
`v.sire` and `v.dam` is indexed by rows, and [`row_to_id`](@ref) and
[`id_to_row`](@ref) translate.

See also [`pedigree_view`](@ref).
"""
struct PedigreeView
    ids::Vector{Int}
    sire::Vector{Int32}
    dam::Vector{Int32}
    row::Vector{Int}          # id → row, 0 when the id is not in the view
end

function PedigreeView(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer}, ids)
    n = length(sire)
    length(dam) == n || throw(ArgumentError("sire and dam differ in length"))
    v = sort!(unique(Int.(collect(ids))))
    isempty(v) || 1 <= first(v) && last(v) <= n ||
        throw(ArgumentError("ids must lie in 1:$n"))
    row = zeros(Int, n)
    row[v] = eachindex(v)
    rs, rd = zeros(Int32, length(v)), zeros(Int32, length(v))
    for (r, i) in enumerate(v)
        s, d = sire[i], dam[i]
        s < i && d < i || throw(ArgumentError("animal $i is not born after its parents"))
        rs[r] = s > 0 ? row[s] : 0
        rd[r] = d > 0 ? row[d] : 0
    end
    return PedigreeView(v, rs, rd, row)
end

Base.length(v::PedigreeView) = length(v.ids)

"""
    row_to_id(v::PedigreeView, r) -> Int

The id of row `r` of the view (elementwise for a vector of rows).
"""
row_to_id(v::PedigreeView, r::Integer) = v.ids[r]
row_to_id(v::PedigreeView, r::AbstractVector{<:Integer}) = v.ids[r]

"""
    id_to_row(v::PedigreeView, id) -> Int

The row of animal `id` in the view; an `ArgumentError` if it is not in it
(elementwise for a vector of ids).
"""
function id_to_row(v::PedigreeView, i::Integer)
    r = 1 <= i <= length(v.row) ? v.row[i] : 0
    r > 0 || throw(ArgumentError("animal $i is not in the pedigree view"))
    return r
end
id_to_row(v::PedigreeView, i::AbstractVector{<:Integer}) = [id_to_row(v, x) for x in i]

"""
    pedigree_view(sire, dam, ids; generations = typemax(Int)) -> PedigreeView
    pedigree_view(pop::Population, ids; generations = typemax(Int)) -> PedigreeView

The view of the animals `ids` and their ancestors, traced through the
pedigree for up to `generations` generations (parents are generation 1,
grandparents 2; `0` gives `ids` alone). With all generations traced, the
relationships among `ids` computed on the view equal those of the full
pedigree; parents beyond the traced depth become unknown, and their
offspring are treated as base animals.
"""
function pedigree_view(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer}, ids;
    generations::Integer = typemax(Int))
    generations >= 0 || throw(ArgumentError("generations must be non-negative"))
    n = length(sire)
    keep = falses(n)
    front = unique(Int.(collect(ids)))
    all(i -> 1 <= i <= n, front) || throw(ArgumentError("ids must lie in 1:$n"))
    keep[front] .= true
    depth = 0
    while !isempty(front) && depth < generations
        next = Int[]
        for i in front, p in (sire[i], dam[i])
            p > 0 && !keep[p] && (keep[p] = true; push!(next, p))
        end
        front = next
        depth += 1
    end
    return PedigreeView(sire, dam, findall(keep))
end
