# Ordinary roots of the default continuation on the spin grid a = 1 - u^2, u uniform in
# [0.01, 1] (60 nodes, u descending), for s = -2..2, |s| <= l <= 15 and n <= nmax(l) (10 for
# l <= 5, decreasing to 2 at l = 15); roots above a = 0.9 polished in 84-bit arithmetic, and
# stagnated near-extremal damped nodes filled from the expansion in 1 - a^2. Read from
# data/spin_root_table.bin: Int32 mode count, then per mode Int8 (s, l, m, n), 60 complex
# Float64 roots and 60 Int8 Leaver inversion indices (NaN and -1 where the default path
# did not accept the root), little endian. Used only as continuation seeds by
# `_tabulated_spin_seed`.
const _SPIN_TABLE_U = collect(range(1.0, 0.01; length=60))
const _SPIN_TABLE_MAXIMUM_SPIN = 0.9995
const _SPIN_ROOT_TABLE_FILE = joinpath(@__DIR__, "data", "spin_root_table.bin")
include_dependency(_SPIN_ROOT_TABLE_FILE)

function _read_spin_root_table(path)
    table = Dict{NTuple{4,Int},Tuple{Vector{ComplexF64},Vector{Int}}}()
    open(path) do io
        count = ltoh(read(io, Int32))
        for _ in 1:count
            key = Tuple(Int.(read!(io, Vector{Int8}(undef, 4))))
            parts = ltoh.(read!(io, Vector{Float64}(undef, 120)))
            roots = [complex(parts[2k-1], parts[2k]) for k in 1:60]
            inversions = Int.(read!(io, Vector{Int8}(undef, 60)))
            table[key] = (roots, inversions)
        end
    end
    return table
end

const _SPIN_ROOT_TABLE = _read_spin_root_table(_SPIN_ROOT_TABLE_FILE)
