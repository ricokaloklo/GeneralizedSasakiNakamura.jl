# Helpers shared by the group files. Each group includes this file into its own module.

if !@isdefined(GSN_TEST_LEVEL)
    const GSN_TEST_LEVEL = Symbol(get(ENV, "GSN_TEST_LEVEL", "full"))
end

const _LEVEL_RANK = Dict(:quick => 1, :full => 2, :extended => 3)

"True when the current level includes cases of level `level` (:quick < :full < :extended)."
at_level(level::Symbol) = _LEVEL_RANK[GSN_TEST_LEVEL] >= _LEVEL_RANK[level]

"Every `stride`-th element at level :quick, everything otherwise."
level_subset(xs, stride::Int) = GSN_TEST_LEVEL === :quick ? xs[1:stride:end] : xs

const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
data_file(name) = joinpath(DATA_DIR, name)

relerr(value, reference) = abs(value - reference) / abs(reference)

"Error in units of `scale` (e.g. the larger of the two amplitudes of a branch)."
scaled_err(value, reference, scale) = abs(value - reference) / scale

"Teukolsky amplitude ratios (incidence/transmission, reflection/transmission) of a radial function."
teukolsky_ratios(R) = (R.incidence_amplitude / R.transmission_amplitude,
    R.reflection_amplitude / R.transmission_amplitude)

"Kerr quantities in the conventions of Homogeneous/Kerr.jl."
kerr_rplus(a) = 1 + sqrt(1 - a^2)
kerr_kappa(a) = sqrt(1 - a^2)
kerr_delta(a, r) = r^2 - 2r + a^2
horizon_k(a, m, omega) = omega - m * a / (2 * kerr_rplus(a))

"""
Wronskian identity between the IN and UP Teukolsky amplitudes (unit transmissions), derived from the
package's asymptotic forms:
    C_up (4 i k r_+ + 2 s kappa) = 2 i omega B_inc .
Returns (lhs, rhs).
"""
function amplitude_wronskian_sides(s, m, a, omega, binc, cup)
    k = horizon_k(a, m, omega)
    lhs = cup * (4im * k * kerr_rplus(a) + 2s * kerr_kappa(a))
    rhs = 2im * omega * binc
    return lhs, rhs
end

"Radial Wronskian Delta^(s+1) (Rin Rup' - Rup Rin') from the Teukolsky solution states at real r."
function radial_wronskian(Rin, Rup, s, a, r)
    vin = Rin.Teukolsky_solution(r)
    vup = Rup.Teukolsky_solution(r)
    return kerr_delta(a, r)^(s + 1) * (vin[1] * vup[2] - vup[1] * vin[2])
end

"Exact frequency of a point-particle harmonic from the orbit's Mino frequencies."
function harmonic_frequency(a, p, e, x, m, n, k)
    frequencies = GSN.KerrGeodesics.kerr_geo_orbit(a, p, e, x)["Frequencies"]
    return (m * frequencies["ϒϕ"] + n * frequencies["ϒr"] + k * frequencies["ϒθ"]) / frequencies["ϒt"]
end

"Quietly evaluate f(); warnings and info logged inside are discarded (they are tested separately)."
quietly(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
