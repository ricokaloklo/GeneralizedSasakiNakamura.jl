# Complex-frequency candidate solver for ISEM. The Teukolsky equation is written for
# P(x) = R(r) / T(x) in the compact coordinate x = (r_+ - r)/(r_+ - r_-), and local series
# solutions at the horizon, at infinity and at ordinary points are matched on a rotated contour
# in the complex x plane. ISEM uses the result as an independent cross-check and as a fallback
# for the GSN route at complex frequencies.
module PEquation

using Logging
using ...GeneralizedSasakiNakamura: Mode, TeukolskyRadialFunction, BoundaryCondition,
    IN, UP, UNIT_TEUKOLSKY_TRANS
using ..TeukolskyTransformation

include("Matching.jl")
using .Matching

export build_teukolsky_function

# Default matching tolerance, and the matching mismatch above which a build is rescued (or,
# failing that, reported with a warning).
const _DEFAULT_TOL = 1e-13
const _DEFAULT_RHOM = 1.0
const _MATCHING_WARNING_TOL = 1e-8

# Truncation order N (at least 15), growing like -20 log10(1 - |a|) near extremality.
function _default_N(a, l)
    return max(Int64(floor(20 - 20 * log10(1 - abs(a)) - 2 * log(l))), 15)
end

_metadata_property(metadata, key::Symbol, default) =
    metadata !== missing && hasproperty(metadata, key) ? getproperty(metadata, key) : default

function _metadata_mismatch(metadata)
    value = _metadata_property(metadata, :split_mismatch, Inf)
    value === missing && return Inf
    mismatch = Float64(value)
    return isfinite(mismatch) ? mismatch : Inf
end

# Matching controls for one build: xm (nothing = automatic), contour radius rhom, truncation
# order N and tolerance.
_controls(; xm, rhom, N, tol) = (xm = xm, rhom = rhom, N = N, tol = tol)
_with_N(c, N::Int) = _controls(; xm = c.xm, rhom = c.rhom, N = N, tol = c.tol)

# N values tried by the local-N rescue: n0, then n0 +- d for increasing d, kept within [10, 120].
function _local_n_rescue_values(n0::Int)
    vals = Int[n0]
    for d in (1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50)
        push!(vals, n0 + d)
        push!(vals, max(10, n0 - d))
    end
    return (n for n in unique(vals) if 10 <= n <= 120)
end

# Compact-coordinate parameters (kappa, epsilon, tau, lambda) and the contour rotation z.
function _rho_parameters(s, l, m, a, omega)
    params = Matching.Parameters.rho_parameters(Matching.Parameters.TeukolskyParameters(s, l, m, a, omega))
    return params.kappa, params.epsilon, params.tau, params.lambda, params.z
end

# R(r) = T(x) P(x) with x = (r_+ - r)/(r_+ - r_-); returns (R, dR/dr, d2R/dr2, residual estimate).
function _teukolsky_solution(P, s, a, omega, m, lambda, kappa, epsilon, tau)
    T, Tp, _ = TeukolskyTransformation.Tx(s, epsilon, kappa, tau)
    function Rsoln(r)
        x = TeukolskyTransformation.r_to_x(r, kappa)
        r = TeukolskyTransformation.x_to_r(x, kappa)
        Pval, Px, _, error = P(x)
        R = T(x) * Pval
        Rp = - (Tp(x) * Pval + T(x) * Px) / (2 * kappa)
        Rpp = TeukolskyTransformation.d2R(r, R, Rp, s, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    return Rsoln
end

# One IN or UP build with fixed controls; returns the Teukolsky function and the matching metadata.
function _build(s, l, m, a, omega, boundary_condition, c)
    if boundary_condition == IN
        P, (incidence_amplitude, reflection_amplitude), metadata = _Pin(s, l, m, a, omega; xm=c.xm, rhom=c.rhom, N=c.N, tol=c.tol, components=true)
    elseif boundary_condition == UP
        P, (incidence_amplitude, reflection_amplitude), metadata = _Pup(s, l, m, a, omega; xm=c.xm, rhom=c.rhom, N=c.N, tol=c.tol, components=true)
    else
        error("the P-equation solver provides the IN and UP basis solutions only.")
    end
    kappa, epsilon, tau, lambda, _ = _rho_parameters(s, l, m, a, omega)
    Rsoln = _teukolsky_solution(P, s, a, omega, m, lambda, kappa, epsilon, tau)
    mode = Mode(s, l, m, a, omega, lambda)
    teuk_func = TeukolskyRadialFunction(mode, boundary_condition, one(incidence_amplitude), incidence_amplitude, reflection_amplitude, P, missing, Rsoln, UNIT_TEUKOLSKY_TRANS)
    return teuk_func, metadata
end

# Builds with controls `c` and warnings silenced; returns the function, metadata, mismatch,
# controls and any caught error.
function _try_build(s, l, m, a, omega, boundary_condition, c; rescue_source = :default)
    try
        teuk_func, metadata = Logging.with_logger(Logging.NullLogger()) do
            _build(s, l, m, a, omega, boundary_condition, c)
        end
        metadata = _attach_control_metadata(metadata, c, rescue_source)
        return (teuk_func = teuk_func, metadata = metadata, mismatch = _metadata_mismatch(metadata), controls = c, error = nothing)
    catch err
        # An invalid requested option (e.g. xm outside the matching range) is
        # reported to the caller, not retried with other controls.
        err isa ArgumentError && rethrow()
        return (teuk_func = nothing, metadata = missing, mismatch = Inf, controls = c, error = err)
    end
end

function _attach_control_metadata(metadata, c, rescue_source)
    metadata === missing && return metadata
    return merge(metadata, (
        control_rescue_source = rescue_source,
        control_xm = c.xm,
        control_rhom = c.rhom,
        control_N = c.N,
        control_tol = c.tol,
        control_small_frequency_expansion = false,
        control_large_frequency_expansion = false,
    ))
end

# Local-N rescue: scan N around base.N and return the first build within tolerance, else the
# one with the smallest mismatch.
function _rescue_by_local_N(s, l, m, a, omega, boundary_condition, base)
    best_any = nothing
    for n in _local_n_rescue_values(Int(base.N))
        result = _try_build(s, l, m, a, omega, boundary_condition, _with_N(base, n); rescue_source = :local_N)
        if result.error === nothing
            if best_any === nothing || result.mismatch < best_any.mismatch
                best_any = result
            end
            result.mismatch <= _MATCHING_WARNING_TOL && return result
        elseif best_any === nothing
            best_any = result
        end
    end
    return best_any
end

# R(r_+ + 1/4) divided by the incidence amplitude, used to compare the normalization of builds.
function _normalized_radial_probe(result, a)
    result.error !== nothing && return nothing
    result.teuk_func.incidence_amplitude === missing && return nothing
    r = 1 + sqrt(1 - a^2) + 0.25
    value = result.teuk_func.Teukolsky_solution(r)[1] / result.teuk_func.incidence_amplitude
    return isfinite(abs(value)) ? value : nothing
end

function _relative_probe_difference(a, b)
    (a === nothing || b === nothing) && return Inf
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
end

function _downward_n_values(n0::Int)
    values = Int[]
    for n in n0:-4:10
        push!(values, n)
    end
    10 in values || push!(values, 10)
    return unique(values)
end

# For N >= 45: scan N upward from 10 in steps of 4 for two consecutive builds whose
# R(r_+ + 1/4)/incidence_amplitude agree to 1e-6, and return that build if `result` disagrees with it.
function _guard_normalization_stability(s, l, m, a, omega, boundary_condition, base, result)
    (result.error === nothing && Int(base.N) >= 45) || return result
    reference_probe = _normalized_radial_probe(result, a)
    reference_probe === nothing && return result

    stable_result = nothing
    previous_probe = nothing
    for n in reverse(_downward_n_values(Int(base.N)))
        candidate = _try_build(s, l, m, a, omega, boundary_condition, _with_N(base, n); rescue_source = :normalization_guard)
        candidate.error === nothing || continue
        candidate.mismatch <= _MATCHING_WARNING_TOL || continue
        probe = _normalized_radial_probe(candidate, a)
        probe === nothing && continue
        if previous_probe !== nothing && _relative_probe_difference(probe, previous_probe) <= 1e-6
            stable_result = candidate
            break
        end
        previous_probe = probe
    end

    stable_result === nothing && return result
    current_probe = _normalized_radial_probe(result, a)
    stable_probe = _normalized_radial_probe(stable_result, a)
    _relative_probe_difference(current_probe, stable_probe) > 1e-6 ? stable_result : result
end

"""
    build_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, return_metadata=false)

Teukolsky `IN` or `UP` solution for complex `omega` from the P equation, as a
`TeukolskyRadialFunction` (with the matching metadata when `return_metadata = true`).
When `N` is not given it starts from a spin-dependent default, rescues a failed or
mismatched build by scanning nearby `N`, and for large `N` checks that the normalization
is stable under a change of `N`. Warns if the final matching mismatch exceeds `1e-8`.
"""
function build_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, return_metadata=false)
    isreal(omega) && throw(ArgumentError("the P-equation solver is used for complex frequencies only."))
    base = _controls(;
        xm = xm,
        rhom = rhom === nothing ? _DEFAULT_RHOM : rhom,
        N = N === nothing ? _default_N(a, l) : N,
        tol = tol === nothing ? _DEFAULT_TOL : tol)
    result = _try_build(s, l, m, a, omega, boundary_condition, base)
    if N === nothing && (result.error !== nothing || result.mismatch > _MATCHING_WARNING_TOL)
        rescued = _rescue_by_local_N(s, l, m, a, omega, boundary_condition, base)
        if rescued !== nothing && rescued.error === nothing && rescued.mismatch <= _MATCHING_WARNING_TOL
            result = rescued
        elseif result.error !== nothing
            result = rescued
        end
    end
    if N === nothing && result !== nothing && result.error === nothing && result.mismatch <= _MATCHING_WARNING_TOL
        result = _guard_normalization_stability(s, l, m, a, omega, boundary_condition, result.controls, result)
    end
    if result.error !== nothing
        error("P-equation radial construction failed after the local-N rescue. Original error: $(sprint(showerror, result.error))")
    end
    if result.mismatch > _MATCHING_WARNING_TOL
        @warn "P-equation radial matching mismatch above tolerance after the local-N rescue" s=s l=l m=m a=a omega=omega boundary_condition=boundary_condition mismatch=result.mismatch controls=result.controls
    end
    return return_metadata ? (result.teuk_func, result.metadata) : result.teuk_func
end

end
