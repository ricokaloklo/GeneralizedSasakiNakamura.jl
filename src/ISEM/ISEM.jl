module ISEM

using ..GeneralizedSasakiNakamura
using ..Coordinates
using ..ConversionFactors
using ..Solutions
using ..InitialConditions
using ..AsymptoticExpansionCoefficients
using SpinWeightedSpheroidalHarmonics
using StaticArrays
using Logging
using LoggingExtras
using Serialization

include("Matching.jl")
using .Matching

export Matching, _Pin, _Pup, _P
export BoundaryCondition, IN, UP, OUT, DOWN
export NormalizationConvention, UNIT_GSN_TRANS, UNIT_TEUKOLSKY_TRANS
export Mode, GSNRadialFunction, TeukolskyRadialFunction, GSN_radial, Teukolsky_radial
export YRadialFunction, Y_radial
export matching_controls

const _DEFAULT_XM = -1.0
const _DEFAULT_RHOM = 1.0
const _DEFAULT_N = 40
const _DEFAULT_TOL = 1e-13
const _GSN_HORIZON_ASYMPTOTIC_ORDER = 12
const _GSN_INFINITY_ASYMPTOTIC_ORDER = 20
const _GSN_HORIZON_DELTA_R_MAX = 0.03
const _GSN_INFINITY_PHASE_MIN = 30.0
const _GSN_ASYMPTOTIC_PATCH_RELATIVE_MISMATCH_MAX = 1e-4
const _RADIAL_MATCHING_WARNING_TOL = 1e-8
const _STATIC_OMEGA_TOL = 1e-12
const _SELECTOR_MODEL_FILE = joinpath(@__DIR__, "selector_mlp_model.jls")
const _SELECTOR_PAYLOAD = Ref{Any}(nothing)
const _CONTROL_SOURCE_SELECTOR = :selector
const _CONTROL_SOURCE_SECOND_LAYER = :second_layer
const _CONTROL_SOURCE_LOCAL_N = :local_N
const _CONTROL_SOURCE_NORMALIZATION_GUARD = :normalization_guard

_is_static_frequency(omega) = abs(omega) < _STATIC_OMEGA_TOL
_horizon_frequency(a, m) = m * a / (2 * (1 + sqrt(1 - a^2)))
_is_horizon_superradiance_frequency(a, m, omega) = isreal(omega) && !_is_static_frequency(omega) && abs(omega - _horizon_frequency(a, m)) < _STATIC_OMEGA_TOL
_use_large_frequency_expansion(a, omega) = isreal(omega) && abs(sqrt(1 - a^2) * omega) > 3

function _is_omega_complex(omega)
    return !iszero(imag(complex(omega)))
end

function _omega_phase_ratio(omega)
    ω = complex(omega)
    denom = max(abs(real(ω)), abs(imag(ω)))
    iszero(denom) && return 0.0
    return abs(imag(ω)) / denom
end

function _selector_available(s, l, m, a, omega, boundary_condition)
    return isreal(omega) &&
        s in (-2, 2) &&
        2 <= l <= 30 &&
        abs(m) <= l &&
        abs(m) <= 30 &&
        abs(a) <= 0.99 &&
        real(omega) > 0 &&
        boundary_condition in (IN, UP) &&
        isfile(_SELECTOR_MODEL_FILE)
end

function _selector_payload()
    payload = _SELECTOR_PAYLOAD[]
    payload !== nothing && return payload
    payload = deserialize(_SELECTOR_MODEL_FILE)
    _SELECTOR_PAYLOAD[] = payload
    return payload
end

function _selector_features(s::Int, l::Int, m::Int, a, omega, boundary_condition)
    lf = Float32(max(l, 1))
    logw = Float32((log10(max(Float64(real(omega)), 1e-12)) + 4.0) / 6.0)
    bup = boundary_condition == UP ? 1.0f0 : 0.0f0
    x = Float32[
        Float32(s) / 2.0f0,
        Float32(l) / 30.0f0,
        Float32(m) / 30.0f0,
        Float32(m) / lf,
        Float32(abs(m)) / lf,
        Float32(a),
        Float32(abs(a)),
        Float32(sign(a)),
        logw,
        bup,
        Float32(l) / 30.0f0 * logw,
        Float32(abs(a)) * logw,
    ]
    append!(x, s == -2 ? Float32[1, 0] : Float32[0, 1])
    for li in 2:30
        push!(x, l == li ? 1.0f0 : 0.0f0)
    end
    for mi in -30:30
        push!(x, m == mi ? 1.0f0 : 0.0f0)
    end
    append!(x, boundary_condition == IN ? Float32[1, 0] : Float32[0, 1])
    return reshape(x, :, 1)
end

function _selector_predict_controls(s::Int, l::Int, m::Int, a, omega, boundary_condition)
    payload = _selector_payload()
    model = payload.model
    X = _selector_features(s, l, m, a, omega, boundary_condition)
    H = max.(model.W1 * X .+ model.b1, 0.0f0)
    ZN = model.WN * H .+ model.bN
    ZB = model.WB * H .+ model.bB
    nidx = argmax(@view ZN[:, 1])
    flags = ntuple(i -> ZB[i, 1] >= 0.0f0, 7)
    return (
        N = Int(payload.N_classes[nidx]),
        sfe = flags[1],
        lfe = flags[2],
        TSinInf = flags[3],
        TSoutInf = flags[4],
        TSinHor = flags[5],
        TSoutHor = flags[6],
        fallback = flags[7],
    )
end

function _is_lfe_stall_risk(s::Int, l::Int, m::Int, omega, boundary_condition)
    return boundary_condition == UP &&
        isreal(omega) &&
        abs(s) == 2 &&
        4.0 <= real(omega) <= 5.0 &&
        8 <= l <= 10 &&
        abs(m) >= l - 1
end

function _selector_lfe_override(s::Int, l::Int, m::Int, omega, boundary_condition)
    return boundary_condition == UP &&
        isreal(omega) &&
        abs(s) == 2 &&
        4.0 <= real(omega) <= 9.0 &&
        2 <= l <= 10 &&
        abs(m) >= l - 1 &&
        !_is_lfe_stall_risk(s, l, m, omega, boundary_condition)
end

function matching_controls(s, l, m, a, omega; boundary_condition=IN)
    if _selector_available(s, l, m, a, omega, boundary_condition)
        pred = _selector_predict_controls(s, l, m, a, omega, boundary_condition)
        lfe = _selector_lfe_override(s, l, m, omega, boundary_condition) ? true : pred.lfe
        return (
            xm = _DEFAULT_XM,
            rhom = _DEFAULT_RHOM,
            N = pred.N,
            tol = _DEFAULT_TOL,
            sfe = pred.sfe,
            lfe = lfe,
            TSinInf = pred.TSinInf,
            TSoutInf = pred.TSoutInf,
            TSinHor = pred.TSinHor,
            TSoutHor = pred.TSoutHor,
            boundary_condition = boundary_condition,
            selector_fallback = pred.fallback,
        )
    end
    absomega = abs(omega)
    small_real_omega = isreal(omega) && absomega < 0.01
    sfe = small_real_omega && abs(s) == 2
    lfe = false
    N = _legacy_matching_N(a, l)
    xm = _DEFAULT_XM
    if boundary_condition == IN && s == 1 && isreal(omega)
        N = max(N, 50)
        xm = -0.6
        TSinInf = false
    elseif boundary_condition == IN && s == 2 && isreal(omega) && iszero(a) && l == 2
        if m == 2
            N = 55
            xm = -0.6
            sfe = absomega < 0.1 ? true : sfe
            TSinInf = false
        elseif m == 1
            N = 60
            xm = -0.6
            TSinInf = false
        elseif m == 0
            N = 80
            xm = -0.8
            TSinInf = false
        elseif m == -1
            N = 25
            xm = -0.4
            TSinInf = false
        else
            N = 40
            xm = -0.5
            TSinInf = false
        end
    elseif boundary_condition == IN && s == 2 && isreal(omega) && l == 2 && m != 2
        N = 10
        xm = -0.9
        TSinInf = false
    elseif boundary_condition == IN && s == 2 && isreal(omega)
        # Direct +2 infinity TSI is inaccurate for small/moderate real frequencies.
        N = max(N, 50)
        xm = -0.6
        TSinInf = false
    else
        TSinInf = s == 2 && isreal(omega)
    end
    TSoutInf = s == -2 && isreal(omega)
    TSinHor = false
    TSoutHor = false
    return (
        xm = xm,
        rhom = _DEFAULT_RHOM,
        N = N,
        tol = _DEFAULT_TOL,
        sfe = sfe,
        lfe = lfe,
        TSinInf = TSinInf,
        TSoutInf = TSoutInf,
        TSinHor = TSinHor,
        TSoutHor = TSoutHor,
        boundary_condition = boundary_condition,
    )
end

function _resolve_option(user_value, default_value)
    return user_value === nothing ? default_value : user_value
end

function _legacy_matching_N(a, l)
    return max(Int64(floor(20 - 20 * log10(1 - abs(a)) - 2 * log(l))), 15)
end

function _y_default_N(a, l, omega)
    return _legacy_matching_N(a, l)
end

function _radial_controls(defaults; xm, rhom, N, tol, sfe, lfe, TSinInf, TSoutInf, TSinHor, TSoutHor)
    return (
        xm = _resolve_option(xm, defaults.xm),
        rhom = _resolve_option(rhom, defaults.rhom),
        N = _resolve_option(N, defaults.N),
        tol = _resolve_option(tol, defaults.tol),
        sfe = _resolve_option(sfe, defaults.sfe),
        lfe = _resolve_option(lfe, defaults.lfe),
        TSinInf = _resolve_option(TSinInf, defaults.TSinInf),
        TSoutInf = _resolve_option(TSoutInf, defaults.TSoutInf),
        TSinHor = _resolve_option(TSinHor, defaults.TSinHor),
        TSoutHor = _resolve_option(TSoutHor, defaults.TSoutHor),
    )
end

function _metadata_mismatch(metadata)
    value = _metadata_property(metadata, :split_mismatch, Inf)
    value === missing && return Inf
    mismatch = Float64(value)
    return isfinite(mismatch) ? mismatch : Inf
end

function _local_n_rescue_values(n0::Int)
    vals = Int[n0]
    for d in (1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50)
        push!(vals, n0 + d)
        push!(vals, max(10, n0 - d))
    end
    return (n for n in unique(vals) if 10 <= n <= 120)
end

function _with_N(c, N::Int)
    return (
        xm = c.xm,
        rhom = c.rhom,
        N = N,
        tol = c.tol,
        sfe = c.sfe,
        lfe = c.lfe,
        TSinInf = c.TSinInf,
        TSoutInf = c.TSoutInf,
        TSinHor = c.TSinHor,
        TSoutHor = c.TSoutHor,
    )
end

function _with_lfe(c, lfe::Bool)
    return (
        xm = c.xm,
        rhom = c.rhom,
        N = c.N,
        tol = c.tol,
        sfe = c.sfe,
        lfe = lfe,
        TSinInf = c.TSinInf,
        TSoutInf = c.TSoutInf,
        TSinHor = c.TSinHor,
        TSoutHor = c.TSoutHor,
    )
end

function _unique_control_candidates(candidates)
    seen = Set{Tuple{Int,Bool,Bool,Bool,Bool,Bool,Bool}}()
    out = Any[]
    for c in candidates
        key = (Int(c.N), Bool(c.sfe), Bool(c.lfe), Bool(c.TSinInf), Bool(c.TSoutInf), Bool(c.TSinHor), Bool(c.TSoutHor))
        key in seen && continue
        push!(seen, key)
        push!(out, c)
    end
    return out
end

function _second_layer_rule_A(s::Int, l::Int, m::Int, omega, boundary_condition)
    return boundary_condition == UP &&
        isreal(omega) &&
        s == 2 &&
        l == 2 &&
        abs(m) >= 1 &&
        0.04 <= real(omega) <= 0.80 &&
        !_is_lfe_stall_risk(s, l, m, omega, boundary_condition)
end

function _second_layer_rule_B(l::Int, omega, boundary_condition, base)
    return boundary_condition == UP &&
        isreal(omega) &&
        l >= 23 &&
        !Bool(base.lfe) &&
        real(omega) <= 1.60
end

function _second_layer_repair_candidates(s, l, m, omega, boundary_condition, base)
    candidates = Any[]
    if _second_layer_rule_A(Int(s), Int(l), Int(m), omega, boundary_condition)
        lfe_base = _with_lfe(base, true)
        push!(candidates, lfe_base)
        for n in (20, 27, 28, 50)
            push!(candidates, _with_N(lfe_base, n))
        end
    elseif _second_layer_rule_B(Int(l), omega, boundary_condition, base)
        push!(candidates, _with_N(base, 20))
        push!(candidates, _with_lfe(_with_N(base, 20), true))
        push!(candidates, base)
        push!(candidates, _with_lfe(base, true))
    end
    return _unique_control_candidates(candidates)
end

function _rescue_teukolsky_by_second_layer_candidates(s, l, m, a, omega, boundary_condition, base)
    best_any = nothing
    for c in _second_layer_repair_candidates(s, l, m, omega, boundary_condition, base)
        result = _try_build_teukolsky_function_with_controls(s, l, m, a, omega, boundary_condition, c; rescue_source = _CONTROL_SOURCE_SECOND_LAYER)
        if result.error === nothing
            if best_any === nothing || result.mismatch < best_any.mismatch
                best_any = result
            end
            result.mismatch <= _RADIAL_MATCHING_WARNING_TOL && return result
        elseif best_any === nothing
            best_any = result
        end
    end
    return best_any
end

function _try_build_teukolsky_function_with_controls(s, l, m, a, omega, boundary_condition, c; rescue_source = _CONTROL_SOURCE_SELECTOR)
    try
        teuk_func, metadata = Logging.with_logger(Logging.NullLogger()) do
            _build_teukolsky_function(
                s,
                l,
                m,
                a,
                omega,
                boundary_condition;
                xm = c.xm,
                rhom = c.rhom,
                N = c.N,
                tol = c.tol,
                sfe = c.sfe,
                lfe = c.lfe,
                TSinInf = c.TSinInf,
                TSoutInf = c.TSoutInf,
                TSinHor = c.TSinHor,
                TSoutHor = c.TSoutHor,
                return_metadata = true,
            )
        end
        metadata = _attach_control_metadata(metadata, c, rescue_source)
        return (teuk_func = teuk_func, metadata = metadata, mismatch = _metadata_mismatch(metadata), controls = c, error = nothing)
    catch err
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
        control_sfe = Bool(c.sfe),
        control_lfe = Bool(c.lfe),
        control_TSinInf = Bool(c.TSinInf),
        control_TSoutInf = Bool(c.TSoutInf),
        control_TSinHor = Bool(c.TSinHor),
        control_TSoutHor = Bool(c.TSoutHor),
    ))
end

function _normalization_probe_radius(a)
    rp = 1 + sqrt(1 - a^2)
    return rp + 0.25
end

function _normalized_radial_probe(result, s, m, a, omega)
    result.error !== nothing && return nothing
    result.teuk_func.incidence_amplitude === missing && return nothing
    r = _normalization_probe_radius(a)
    r <= 1 + sqrt(1 - a^2) && return nothing
    value = result.teuk_func.Teukolsky_solution(r)[1] / result.teuk_func.incidence_amplitude
    return isfinite(abs(value)) ? value : nothing
end

function _relative_probe_difference(a, b)
    (a === nothing || b === nothing) && return Inf
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
end

function _needs_normalization_stability_guard(s, l, m, omega, boundary_condition, base, result)
    result.error !== nothing && return false
    a = result.teuk_func.mode.a
    !_selector_available(s, l, m, a, omega, boundary_condition) && Int(base.N) >= 45 && return true
    return false
end

function _downward_n_values(n0::Int)
    values = Int[]
    for n in n0:-4:10
        push!(values, n)
    end
    10 in values || push!(values, 10)
    return unique(values)
end

function _guard_teukolsky_normalization_stability(s, l, m, a, omega, boundary_condition, base, result)
    _needs_normalization_stability_guard(s, l, m, omega, boundary_condition, base, result) || return result
    reference_probe = _normalized_radial_probe(result, s, m, a, omega)
    reference_probe === nothing && return result

    stable_result = nothing
    previous_result = nothing
    previous_probe = nothing
    for n in reverse(_downward_n_values(Int(base.N)))
        c = _with_N(base, n)
        candidate = _try_build_teukolsky_function_with_controls(s, l, m, a, omega, boundary_condition, c; rescue_source = _CONTROL_SOURCE_NORMALIZATION_GUARD)
        candidate.error === nothing || continue
        candidate.mismatch <= _RADIAL_MATCHING_WARNING_TOL || continue
        probe = _normalized_radial_probe(candidate, s, m, a, omega)
        probe === nothing && continue
        if previous_probe !== nothing && _relative_probe_difference(probe, previous_probe) <= 1e-6
            stable_result = candidate
            break
        end
        previous_result = candidate
        previous_probe = probe
    end

    stable_result === nothing && return result
    current_probe = _normalized_radial_probe(result, s, m, a, omega)
    stable_probe = _normalized_radial_probe(stable_result, s, m, a, omega)
    _relative_probe_difference(current_probe, stable_probe) > 1e-6 ? stable_result : result
end

function _rescue_teukolsky_by_local_N(s, l, m, a, omega, boundary_condition, base)
    best_valid = nothing
    best_any = nothing
    high_lfe_scale = _use_large_frequency_expansion(a, omega)
    candidates = if high_lfe_scale
        Bool(base.lfe) ? (base, _with_lfe(base, false)) : (base, _with_lfe(base, true))
    elseif isreal(omega) && !Bool(base.lfe)
        (base, _with_lfe(base, true))
    else
        (base,)
    end
    for candidate in candidates
        for n in _local_n_rescue_values(Int(candidate.N))
            c = _with_N(candidate, n)
            result = _try_build_teukolsky_function_with_controls(s, l, m, a, omega, boundary_condition, c; rescue_source = _CONTROL_SOURCE_LOCAL_N)
            if result.error === nothing
                if best_any === nothing || result.mismatch < best_any.mismatch
                    best_any = result
                end
                if result.mismatch <= _RADIAL_MATCHING_WARNING_TOL
                    best_valid = result
                    break
                end
            elseif best_any === nothing
                best_any = result
            end
        end
        best_valid !== nothing && break
    end
    return best_valid === nothing ? best_any : best_valid
end

function _validate_y_branch(s, boundary_condition)
    if !(s == -2 && boundary_condition == IN || s == 2 && boundary_condition == UP)
        error("Y_radial currently supports only s = -2 with IN and s = +2 with UP.")
    end
end

function _y_phase_argument(r, a, m, omega)
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    rs = rstar_from_r(a, r)
    return omega * rs + a * m * log((r - rm) / (r - rp)) / (rp - rm)
end

function _y_phase_derivative(r, a, m, omega)
    return ((r^2 + a^2) * omega - a * m) / (r^2 - 2 * r + a^2)
end

function _y_branch_transmission_amplitude(s, m, a, omega, lambda)
    if s == -2
        return ConversionFactors.Btrans(s, m, a, omega, lambda)
    elseif s == 2
        return ConversionFactors.Ctrans(s, m, a, omega, lambda)
    else
        error("Unsupported Y branch.")
    end
end

struct YRadialFunction
    mode::Mode
    boundary_condition::BoundaryCondition
    transmission_amplitude
    incidence_amplitude
    reflection_amplitude
    P_solution::Union{Function, Missing}
    Teukolsky_solution::Union{Function, Missing}
    X_solution::Union{Function, Missing}
    Y_scalar_solution::Function
    Y_solution::Function
    normalization_convention::NormalizationConvention
end

function Base.show(io::IO, ::MIME"text/plain", y_func::YRadialFunction)
    println(io, "YRadialFunction(")
    print(io, "    mode = "); show(io, "text/plain", y_func.mode); println(io, ",")
    println(io, "    boundary_condition = $(y_func.boundary_condition),")
    println(io, "    transmission_amplitude = $(y_func.transmission_amplitude),")
    println(io, "    incidence_amplitude = $(y_func.incidence_amplitude),")
    println(io, "    reflection_amplitude = $(y_func.reflection_amplitude),")
    println(io, "    normalization_convention = $(y_func.normalization_convention)")
    print(io, ")")
end

function Base.show(io::IO, y_func::YRadialFunction)
    print(io, "YRadialFunction(mode = Mode(s = $(y_func.mode.s), l = $(y_func.mode.l), m = $(y_func.mode.m), a = $(y_func.mode.a), omega = $(y_func.mode.omega), lambda = $(y_func.mode.lambda)), boundary_condition = $(y_func.boundary_condition))")
end

function _isem_parameters(s, l, m, a, omega)
    params = Matching.Parameters.isem_parameters(Matching.Parameters.TeukolskyParameters(s, l, m, a, omega))
    return params.kappa, params.epsilon, params.tau, params.lambda
end

function _rho_parameters(s, l, m, a, omega)
    params = Matching.Parameters.rho_parameters(Matching.Parameters.TeukolskyParameters(s, l, m, a, omega))
    return params.kappa, params.epsilon, params.tau, params.lambda, params.z
end

_is_negative_real_frequency(omega) = isreal(omega) && real(omega) < -_STATIC_OMEGA_TOL
_conjugate_or_missing(x) = x === missing ? missing : conj(x)
_conjugate_solution_value(x) = x isa Tuple ? map(_conjugate_solution_value, x) : _conjugate_or_missing(x)
_conjugate_solution_function(f) = f === missing ? missing : (args...; kwargs...) -> _conjugate_solution_value(f(args...; kwargs...))

function _conjugate_teukolsky_radial_function(f::TeukolskyRadialFunction, m::Int, omega)
    mode = Mode(f.mode.s, f.mode.l, m, f.mode.a, omega, f.mode.lambda)
    return TeukolskyRadialFunction(
        mode,
        f.boundary_condition,
        _conjugate_or_missing(f.transmission_amplitude),
        _conjugate_or_missing(f.incidence_amplitude),
        _conjugate_or_missing(f.reflection_amplitude),
        _conjugate_solution_function(f.P_solution),
        f.GSN_solution === missing ? missing : _conjugate_gsn_radial_function(f.GSN_solution, m, omega),
        _conjugate_solution_function(f.Teukolsky_solution),
        f.normalization_convention,
    )
end

function _conjugate_gsn_radial_function(f::GSNRadialFunction, m::Int, omega)
    mode = Mode(f.mode.s, f.mode.l, m, f.mode.a, omega, f.mode.lambda)
    return GSNRadialFunction(
        mode,
        f.boundary_condition,
        f.rsin,
        f.rsout,
        f.rsmp,
        f.horizon_expansion_order,
        f.infinity_expansion_order,
        _conjugate_or_missing(f.transmission_amplitude),
        _conjugate_or_missing(f.incidence_amplitude),
        _conjugate_or_missing(f.reflection_amplitude),
        f.numerical_GSN_solution,
        _conjugate_solution_function(f.numerical_Riccati_solution),
        _conjugate_solution_function(f.GSN_solution),
        f.normalization_convention,
        f.method,
    )
end

function _conjugate_y_radial_function(f::YRadialFunction, m::Int, omega)
    mode = Mode(f.mode.s, f.mode.l, m, f.mode.a, omega, f.mode.lambda)
    return YRadialFunction(
        mode,
        f.boundary_condition,
        _conjugate_or_missing(f.transmission_amplitude),
        _conjugate_or_missing(f.incidence_amplitude),
        _conjugate_or_missing(f.reflection_amplitude),
        _conjugate_solution_function(f.P_solution),
        _conjugate_solution_function(f.Teukolsky_solution),
        _conjugate_solution_function(f.X_solution),
        _conjugate_solution_function(f.Y_scalar_solution),
        _conjugate_solution_function(f.Y_solution),
        f.normalization_convention,
    )
end

_kerr_delta(r, a) = r^2 + a^2 - 2 * r
_kerr_rplus(a) = 1 + sqrt(1 - a^2)
_kerr_omega_horizon(a) = a / (2 * _kerr_rplus(a))

function _eval_inf_ansatz(coeffs, omega, r)
    zinv = inv(omega * r)
    zpow = one(zinv)
    ans = zero(eltype(coeffs))
    @inbounds for c in coeffs
        ans += c * zpow
        zpow *= zinv
    end
    return ans
end

function _eval_dinf_ansatz_dr(coeffs, omega, r)
    ans = zero(eltype(coeffs))
    omega_inv = inv(omega)
    omega_pow = omega_inv
    r_pow_inv = inv(r^2)
    r_inv = inv(r)
    @inbounds for i in 1:(length(coeffs) - 1)
        ans -= i * coeffs[i + 1] * omega_pow * r_pow_inv
        omega_pow *= omega_inv
        r_pow_inv *= r_inv
    end
    return ans
end

function _eval_hor_ansatz(coeffs, dr)
    drpow = one(dr)
    ans = zero(eltype(coeffs))
    @inbounds for c in coeffs
        ans += c * drpow
        drpow *= dr
    end
    return ans
end

function _eval_dhor_ansatz_dr(coeffs, dr)
    ans = zero(eltype(coeffs))
    drpow = one(dr)
    @inbounds for i in 1:(length(coeffs) - 1)
        ans += i * coeffs[i + 1] * drpow
        drpow *= dr
    end
    return ans
end

function _lazy_coefficients!(coeffs, builder)
    if coeffs[] === nothing
        coeffs[] = builder()
    end
    return coeffs[]
end

function _gsn_horizon_ansatz(mode::Mode, boundary_condition::BoundaryCondition, transmission_amplitude, incidence_amplitude, reflection_amplitude, rs; order = _GSN_HORIZON_ASYMPTOTIC_ORDER)
    s = mode.s
    m = mode.m
    a = mode.a
    omega = mode.omega
    lambda = mode.lambda
    r = r_from_rstar(a, rs)
    p = omega - m * _kerr_omega_horizon(a)

    gin(r0) = InitialConditions.gansatz(
        ord -> AsymptoticExpansionCoefficients.ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r0;
        order = order,
    )
    dgin_dr(r0) = InitialConditions.dgansatz_dr(
        ord -> AsymptoticExpansionCoefficients.ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r0;
        order = order,
    )
    gout(r0) = InitialConditions.gansatz(
        ord -> AsymptoticExpansionCoefficients.outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r0;
        order = order,
    )
    dgout_dr(r0) = InitialConditions.dgansatz_dr(
        ord -> AsymptoticExpansionCoefficients.outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r0;
        order = order,
    )

    Δfac = _kerr_delta(r, a) / (r^2 + a^2)
    if boundary_condition == IN
        phase = exp(-1im * p * rs)
        X = transmission_amplitude * gin(r) * phase
        Xp = transmission_amplitude * (Δfac * dgin_dr(r) - 1im * p * gin(r)) * phase
        return X, Xp
    elseif boundary_condition == UP
        phase_out = exp(1im * p * rs)
        phase_in = exp(-1im * p * rs)
        X = incidence_amplitude * gout(r) * phase_out + reflection_amplitude * gin(r) * phase_in
        Xp = incidence_amplitude * (Δfac * dgout_dr(r) + 1im * p * gout(r)) * phase_out +
             reflection_amplitude * (Δfac * dgin_dr(r) - 1im * p * gin(r)) * phase_in
        return X, Xp
    end
    return nothing
end

function _gsn_infinity_ansatz(mode::Mode, boundary_condition::BoundaryCondition, transmission_amplitude, incidence_amplitude, reflection_amplitude, rs; order = _GSN_INFINITY_ASYMPTOTIC_ORDER)
    s = mode.s
    m = mode.m
    a = mode.a
    omega = mode.omega
    lambda = mode.lambda
    r = r_from_rstar(a, rs)

    fin(r0) = InitialConditions.fansatz(
        ord -> AsymptoticExpansionCoefficients.ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r0;
        order = order,
    )
    dfin_dr(r0) = InitialConditions.dfansatz_dr(
        ord -> AsymptoticExpansionCoefficients.ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r0;
        order = order,
    )
    fout(r0) = InitialConditions.fansatz(
        ord -> AsymptoticExpansionCoefficients.outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r0;
        order = order,
    )
    dfout_dr(r0) = InitialConditions.dfansatz_dr(
        ord -> AsymptoticExpansionCoefficients.outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r0;
        order = order,
    )

    Δfac = _kerr_delta(r, a) / (r^2 + a^2)
    if boundary_condition == UP
        phase = exp(1im * omega * rs)
        X = transmission_amplitude * fout(r) * phase
        Xp = transmission_amplitude * (Δfac * dfout_dr(r) + 1im * omega * fout(r)) * phase
        return X, Xp
    elseif boundary_condition == IN
        phase_out = exp(1im * omega * rs)
        phase_in = exp(-1im * omega * rs)
        X = reflection_amplitude * fout(r) * phase_out + incidence_amplitude * fin(r) * phase_in
        Xp = reflection_amplitude * (Δfac * dfout_dr(r) + 1im * omega * fout(r)) * phase_out +
             incidence_amplitude * (Δfac * dfin_dr(r) - 1im * omega * fin(r)) * phase_in
        return X, Xp
    end
    return nothing
end

function _relative_mismatch(a, b)
    scale = max(abs(a), abs(b), eps(float(real(abs(a) + abs(b)))))
    return abs(a - b) / scale
end

function _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
    raw = gsn_solution(rs)
    mismatch = max(
        _relative_mismatch(X, raw[1]),
        _relative_mismatch(Xp, raw[2]),
    )
    mismatch > _GSN_ASYMPTOTIC_PATCH_RELATIVE_MISMATCH_MAX && return raw
    return (X, Xp, 0.0)
end

function _with_gsn_asymptotic_patches(gsn_solution, mode::Mode, boundary_condition::BoundaryCondition, transmission_amplitude, incidence_amplitude, reflection_amplitude; horizon_delta_r_max = _GSN_HORIZON_DELTA_R_MAX, infinity_phase_min = _GSN_INFINITY_PHASE_MIN)
    _is_omega_complex(mode.omega) && return gsn_solution
    s = mode.s
    m = mode.m
    a = mode.a
    omega = mode.omega
    lambda = mode.lambda
    rp = _kerr_rplus(a)
    p = omega - m * _kerr_omega_horizon(a)

    hor_order = _GSN_HORIZON_ASYMPTOTIC_ORDER
    inf_order = _GSN_INFINITY_ASYMPTOTIC_ORDER
    real_type = promote_type(typeof(float(real(a))), typeof(float(real(omega))), typeof(float(real(lambda))))
    coeff_type = Complex{real_type}
    hor_in_coeffs = Ref{Union{Nothing, Vector{coeff_type}}}(nothing)
    hor_out_coeffs = Ref{Union{Nothing, Vector{coeff_type}}}(nothing)
    inf_in_coeffs = Ref{Union{Nothing, Vector{coeff_type}}}(nothing)
    inf_out_coeffs = Ref{Union{Nothing, Vector{coeff_type}}}(nothing)
    get_hor_in_coeffs() = _lazy_coefficients!(hor_in_coeffs, () -> coeff_type[AsymptoticExpansionCoefficients.ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord) for ord in 0:hor_order])
    get_hor_out_coeffs() = _lazy_coefficients!(hor_out_coeffs, () -> coeff_type[AsymptoticExpansionCoefficients.outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord) for ord in 0:hor_order])
    get_inf_in_coeffs() = _lazy_coefficients!(inf_in_coeffs, () -> coeff_type[AsymptoticExpansionCoefficients.ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord) for ord in 0:inf_order])
    get_inf_out_coeffs() = _lazy_coefficients!(inf_out_coeffs, () -> coeff_type[AsymptoticExpansionCoefficients.outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord) for ord in 0:inf_order])

    return rs -> begin
        r = r_from_rstar(a, rs)
        Δfac = _kerr_delta(r, a) / (r^2 + a^2)
        if r - rp <= horizon_delta_r_max
            dr = r - rp
            coeffs_in = get_hor_in_coeffs()
            gin = _eval_hor_ansatz(coeffs_in, dr)
            dgin = _eval_dhor_ansatz_dr(coeffs_in, dr)
            if boundary_condition == IN
                phase = exp(-1im * p * rs)
                X = transmission_amplitude * gin * phase
                Xp = transmission_amplitude * (Δfac * dgin - 1im * p * gin) * phase
                return _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
            elseif boundary_condition == UP
                coeffs_out = get_hor_out_coeffs()
                gout = _eval_hor_ansatz(coeffs_out, dr)
                dgout = _eval_dhor_ansatz_dr(coeffs_out, dr)
                phase_out = exp(1im * p * rs)
                phase_in = exp(-1im * p * rs)
                X = incidence_amplitude * gout * phase_out + reflection_amplitude * gin * phase_in
                Xp = incidence_amplitude * (Δfac * dgout + 1im * p * gout) * phase_out +
                     reflection_amplitude * (Δfac * dgin - 1im * p * gin) * phase_in
                return _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
            end
        elseif abs(omega * r) >= infinity_phase_min
            coeffs_out = get_inf_out_coeffs()
            fout = _eval_inf_ansatz(coeffs_out, omega, r)
            dfout = _eval_dinf_ansatz_dr(coeffs_out, omega, r)
            if boundary_condition == UP
                phase = exp(1im * omega * rs)
                X = transmission_amplitude * fout * phase
                Xp = transmission_amplitude * (Δfac * dfout + 1im * omega * fout) * phase
                return _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
            elseif boundary_condition == IN
                coeffs_in = get_inf_in_coeffs()
                fin = _eval_inf_ansatz(coeffs_in, omega, r)
                dfin = _eval_dinf_ansatz_dr(coeffs_in, omega, r)
                phase_out = exp(1im * omega * rs)
                phase_in = exp(-1im * omega * rs)
                X = reflection_amplitude * fout * phase_out + incidence_amplitude * fin * phase_in
                Xp = reflection_amplitude * (Δfac * dfout + 1im * omega * fout) * phase_out +
                     incidence_amplitude * (Δfac * dfin - 1im * omega * fin) * phase_in
                return _continuous_gsn_patch_or_raw(gsn_solution, rs, X, Xp)
            end
        end
        return gsn_solution(rs)
    end
end

function _teukolsky_solution_real(P, s, a, omega, m, lambda, kappa, epsilon, tau)
    T, Tp, _ = Matching.TeukolskyTransformation.Tx(s, epsilon, kappa, tau)
    function Rsoln(r)
        x = Matching.TeukolskyTransformation.r_to_x(r, kappa)
        Pval, Px, _, error = P(x)
        R = T(x) * Pval
        Rp = - (Tp(x) * Pval + T(x) * Px) / (2 * kappa)
        Rpp = Matching.TeukolskyTransformation.d2R(r, R, Rp, s, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    return Rsoln
end

function _teukolsky_solution_contour(P, s, a, omega, m, lambda, kappa, epsilon, tau, z)
    T, Tp, _ = Matching.TeukolskyTransformation.Tx(s, epsilon, kappa, tau)
    function Rsoln(r)
        x = Matching.TeukolskyTransformation.r_to_x(r, kappa)
        r = Matching.TeukolskyTransformation.x_to_r(x, kappa)
        Pval, Px, _, error = P(x)
        R = T(x) * Pval
        Rp = - (Tp(x) * Pval + T(x) * Px) / (2 * kappa)
        Rpp = Matching.TeukolskyTransformation.d2R(r, R, Rp, s, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    return Rsoln
end

function _build_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, return_metadata=false)
    if boundary_condition == IN
        P, (incidence_amplitude, reflection_amplitude), metadata = _Pin(s, l, m, a, omega; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, components=true)
    elseif boundary_condition == UP
        P, (incidence_amplitude, reflection_amplitude), metadata = _Pup(s, l, m, a, omega; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, components=true)
    else
        error("ISEM currently provides the IN and UP basis solutions only.")
    end

    if _is_omega_complex(omega)
        kappa, epsilon, tau, lambda, z = _rho_parameters(s, l, m, a, omega)
        Rsoln = _teukolsky_solution_contour(P, s, a, omega, m, lambda, kappa, epsilon, tau, z)
        mode = Mode(s, l, m, a, omega, lambda)
        teuk_func = TeukolskyRadialFunction(mode, boundary_condition, one(incidence_amplitude), incidence_amplitude, reflection_amplitude, P, missing, Rsoln, UNIT_TEUKOLSKY_TRANS)
        return return_metadata ? (teuk_func, metadata) : teuk_func
    end

    kappa, epsilon, tau, lambda = _isem_parameters(s, l, m, a, omega)
    Rsoln = _teukolsky_solution_real(P, s, a, omega, m, lambda, kappa, epsilon, tau)
    mode = Mode(s, l, m, a, omega, lambda)
    teuk_func = TeukolskyRadialFunction(mode, boundary_condition, one(incidence_amplitude), incidence_amplitude, reflection_amplitude, P, missing, Rsoln, UNIT_TEUKOLSKY_TRANS)
    return return_metadata ? (teuk_func, metadata) : teuk_func
end

function _build_superradiance_threshold_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, return_metadata=false)
    boundary_condition == UP || error("Superradiance-threshold ISEM special branch currently supports only UP.")
    defaults = matching_controls(s, l, m, a, omega; boundary_condition = boundary_condition)
    controls = _radial_controls(defaults; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    P, _, metadata = _Pup(s, l, m, a, omega; xm=controls.xm, rhom=controls.rhom, N=controls.N, tol=controls.tol, sfe=controls.sfe, lfe=controls.lfe, TSoutInf=controls.TSoutInf, TSinHor=controls.TSinHor, TSoutHor=controls.TSoutHor, components=true)
    kappa, epsilon, tau, lambda = _isem_parameters(s, l, m, a, omega)
    Rsoln = _teukolsky_solution_real(P, s, a, omega, m, lambda, kappa, epsilon, tau)
    mode = Mode(s, l, m, a, omega, lambda)
    teuk_func = TeukolskyRadialFunction(mode, boundary_condition, 1.0 + 0.0im, missing, missing, P, missing, Rsoln, UNIT_TEUKOLSKY_TRANS)
    return return_metadata ? (teuk_func, metadata) : teuk_func
end

function _build_teukolsky_function_fixed(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, return_metadata=false)
    defaults = matching_controls(s, l, m, a, omega; boundary_condition = boundary_condition)
    base = _radial_controls(defaults; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    result = _try_build_teukolsky_function_with_controls(s, l, m, a, omega, boundary_condition, base)
    if N === nothing && (result.error !== nothing || result.mismatch > _RADIAL_MATCHING_WARNING_TOL)
        second_layer = _rescue_teukolsky_by_second_layer_candidates(s, l, m, a, omega, boundary_condition, base)
        if second_layer !== nothing && second_layer.error === nothing && second_layer.mismatch <= _RADIAL_MATCHING_WARNING_TOL
            result = second_layer
        end
    end
    if N === nothing && (result.error !== nothing || result.mismatch > _RADIAL_MATCHING_WARNING_TOL)
        rescued = _rescue_teukolsky_by_local_N(s, l, m, a, omega, boundary_condition, base)
        if rescued !== nothing && rescued.error === nothing && rescued.mismatch <= _RADIAL_MATCHING_WARNING_TOL
            result = rescued
        elseif result.error !== nothing
            result = rescued
        end
    end
    if N === nothing && result !== nothing && result.error === nothing && result.mismatch <= _RADIAL_MATCHING_WARNING_TOL
        result = _guard_teukolsky_normalization_stability(s, l, m, a, omega, boundary_condition, result.controls, result)
    end
    if result.error !== nothing
        error("ISEM radial construction failed after selector/second-layer/local-N controls; consider method = \"linear\" or method = \"Riccati\". Original error: $(sprint(showerror, result.error))")
    end
    if result.mismatch > _RADIAL_MATCHING_WARNING_TOL
        @warn "ISEM radial matching mismatch above tolerance after selector/second-layer/local-N controls; method = \"auto\" will switch to legacy auto (Riccati, then linear)" s=s l=l m=m a=a omega=omega boundary_condition=boundary_condition mismatch=result.mismatch controls=result.controls
    end
    return return_metadata ? (result.teuk_func, result.metadata) : result.teuk_func
end

function _combine_down(Rin::TeukolskyRadialFunction, Rup::TeukolskyRadialFunction)
    Binc = Rin.incidence_amplitude
    Bref = Rin.reflection_amplitude
    Cinc = Rup.incidence_amplitude
    Cref = Rup.reflection_amplitude

    function teuk_func(coord)
        Rin_vals = Rin.Teukolsky_solution(coord)
        Rup_vals = Rup.Teukolsky_solution(coord)
        R = Binc^-1 * (Rin_vals[1] - Bref * Rup_vals[1])
        dR = Binc^-1 * (Rin_vals[2] - Bref * Rup_vals[2])
        d2R = Binc^-1 * (Rin_vals[3] - Bref * Rup_vals[3])
        err = max(Rin_vals[4], Rup_vals[4])
        return (R, dR, d2R, err)
    end

    return TeukolskyRadialFunction(
        Rin.mode,
        DOWN,
        one(Binc),
        Binc^-1 - Bref * Cref / Binc,
        -Bref * Cinc / Binc,
        missing,
        missing,
        teuk_func,
        UNIT_TEUKOLSKY_TRANS,
    )
end

function _combine_out(Rin::TeukolskyRadialFunction, Rup::TeukolskyRadialFunction)
    Binc = Rin.incidence_amplitude
    Cinc = Rup.incidence_amplitude
    Cref = Rup.reflection_amplitude

    function teuk_func(coord)
        Rin_vals = Rin.Teukolsky_solution(coord)
        Rup_vals = Rup.Teukolsky_solution(coord)
        R = Cinc^-1 * (Rup_vals[1] - Cref * Rin_vals[1])
        dR = Cinc^-1 * (Rup_vals[2] - Cref * Rin_vals[2])
        d2R = Cinc^-1 * (Rup_vals[3] - Cref * Rin_vals[3])
        err = max(Rin_vals[4], Rup_vals[4])
        return (R, dR, d2R, err)
    end

    return TeukolskyRadialFunction(
        Rin.mode,
        OUT,
        one(Cinc),
        Cinc^-1 - Rin.reflection_amplitude * Cref / Cinc,
        -Binc * Cref / Cinc,
        missing,
        missing,
        teuk_func,
        UNIT_TEUKOLSKY_TRANS,
    )
end

function _gsn_amplitudes_from_teukolsky(teuk_func::TeukolskyRadialFunction)
    s = teuk_func.mode.s
    m = teuk_func.mode.m
    a = teuk_func.mode.a
    omega = teuk_func.mode.omega
    lambda = teuk_func.mode.lambda

    if teuk_func.boundary_condition == IN
        trans_factor = ConversionFactors.Btrans(s, m, a, omega, lambda)
        incidence_factor = ConversionFactors.Binc(s, m, a, omega, lambda)
        reflection_factor = ConversionFactors.Bref(s, m, a, omega, lambda)
    elseif teuk_func.boundary_condition == UP
        trans_factor = ConversionFactors.Ctrans(s, m, a, omega, lambda)
        incidence_factor = ConversionFactors.Cinc(s, m, a, omega, lambda)
        reflection_factor = ConversionFactors.Cref(s, m, a, omega, lambda)
    elseif teuk_func.boundary_condition == DOWN
        trans_factor = ConversionFactors.Binc(s, m, a, omega, lambda)
        incidence_factor = ConversionFactors.Btrans(s, m, a, omega, lambda)
        reflection_factor = ConversionFactors.Cinc(s, m, a, omega, lambda)
    elseif teuk_func.boundary_condition == OUT
        trans_factor = ConversionFactors.Cinc(s, m, a, omega, lambda)
        incidence_factor = ConversionFactors.Bref(s, m, a, omega, lambda)
        reflection_factor = ConversionFactors.Binc(s, m, a, omega, lambda)
    else
        error("Does not understand the boundary condition applied to the solution")
    end

    transmission_amplitude = teuk_func.transmission_amplitude / trans_factor
    incidence_amplitude = teuk_func.incidence_amplitude === missing ? missing : teuk_func.incidence_amplitude / incidence_factor
    reflection_amplitude = teuk_func.reflection_amplitude === missing ? missing : teuk_func.reflection_amplitude / reflection_factor
    return transmission_amplitude, incidence_amplitude, reflection_amplitude
end

_metadata_property(metadata, key::Symbol, default) = metadata !== missing && hasproperty(metadata, key) ? getproperty(metadata, key) : default

function _isem_parameter_summary(; xm, rhom, N, tol, sfe, lfe, TSinInf, TSoutInf, TSinHor, TSoutHor, metadata=missing)
    return (
        xm = _metadata_property(metadata, :xsplit, _metadata_property(metadata, :xm_match, xm)),
        rhom = _metadata_property(metadata, :control_rhom, _metadata_property(metadata, :rhom, rhom)),
        N = _metadata_property(metadata, :control_N, _metadata_property(metadata, :N, N)),
        tol = _metadata_property(metadata, :control_tol, tol),
        sfe = Bool(_metadata_property(metadata, :control_sfe, sfe)),
        lfe = Bool(_metadata_property(metadata, :control_lfe, lfe)),
        TSinInf = Bool(_metadata_property(metadata, :control_TSinInf, TSinInf)),
        TSoutInf = Bool(_metadata_property(metadata, :control_TSoutInf, TSoutInf)),
        TSinHor = Bool(_metadata_property(metadata, :control_TSinHor, TSinHor)),
        TSoutHor = Bool(_metadata_property(metadata, :control_TSoutHor, TSoutHor)),
        rescue_source = _metadata_property(metadata, :control_rescue_source, _CONTROL_SOURCE_SELECTOR),
        xm_requested = _metadata_property(metadata, :xm, xm),
        xm_min = _metadata_property(metadata, :xm_min, missing),
        xm_max = _metadata_property(metadata, :xm_max, missing),
        xm_match = _metadata_property(metadata, :xm_match, missing),
        split_mismatch = _metadata_property(metadata, :split_mismatch, missing),
    )
end

function _build_gsn_function_from_teukolsky(teuk_func::TeukolskyRadialFunction; isem_parameters = missing, use_asymptotic_patches = true, gsn_horizon_delta_r_max = _GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min = _GSN_INFINITY_PHASE_MIN)
    mode = teuk_func.mode
    teukolsky_from_gsn_matrix = r -> begin
        return Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(mode.s, mode.m, mode.a, mode.omega, mode.lambda, r)
    end
    if teuk_func.P_solution !== missing && !_is_omega_complex(mode.omega)
        gsn_solution = Matching.TeukolskyTransformation.P_to_GSN_solution_from_matrix(rs -> r_from_rstar(mode.a, rs), teukolsky_from_gsn_matrix, teuk_func.P_solution, mode.s, mode.m, mode.a, mode.omega, mode.lambda)
    else
        teuk_to_gsn_matrix = r -> inv(teukolsky_from_gsn_matrix(r))
        gsn_solution = Matching.TeukolskyTransformation.Teukolsky_to_GSN_solution(rs -> r_from_rstar(mode.a, rs), teuk_to_gsn_matrix, teuk_func.Teukolsky_solution)
    end
    transmission_amplitude, incidence_amplitude, reflection_amplitude = _gsn_amplitudes_from_teukolsky(teuk_func)
    gsn_transmission_amplitude = transmission_amplitude
    raw_gsn_solution = gsn_solution
    gsn_solution_unit_trans = rs -> begin
        vals = raw_gsn_solution(rs)
        return vals ./ gsn_transmission_amplitude
    end
    incidence_amplitude = incidence_amplitude === missing ? missing : incidence_amplitude / gsn_transmission_amplitude
    reflection_amplitude = reflection_amplitude === missing ? missing : reflection_amplitude / gsn_transmission_amplitude
    transmission_amplitude = one(gsn_transmission_amplitude)
    can_use_asymptotic_patches = use_asymptotic_patches && incidence_amplitude !== missing && reflection_amplitude !== missing
    gsn_solution = can_use_asymptotic_patches ?
        _with_gsn_asymptotic_patches(
            gsn_solution_unit_trans,
            mode,
            teuk_func.boundary_condition,
            transmission_amplitude,
            incidence_amplitude,
            reflection_amplitude,
            horizon_delta_r_max = gsn_horizon_delta_r_max,
            infinity_phase_min = gsn_infinity_phase_min,
        ) :
        gsn_solution_unit_trans

    return GSNRadialFunction(
        mode,
        teuk_func.boundary_condition,
        missing,
        missing,
        missing,
        missing,
        missing,
        transmission_amplitude,
        incidence_amplitude,
        reflection_amplitude,
        isem_parameters,
        missing,
        gsn_solution,
        UNIT_GSN_TRANS,
        "ISEM",
    )
end

function Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition::BoundaryCondition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    if _is_static_frequency(omega)
        lambda = spin_weighted_spherical_eigenvalue(s, l, m)
        mode = Mode(s, l, m, a, zero(omega), lambda)
        if boundary_condition == IN
            teuk_func = Solutions.solve_static_Rin(s, l, m, a)
        elseif boundary_condition == UP
            teuk_func = Solutions.solve_static_Rup(s, l, m, a)
        else
            error("Static Teukolsky modes currently support only IN and UP boundary conditions.")
        end

        return TeukolskyRadialFunction(
            mode,
            boundary_condition,
            1,
            missing,
            missing,
            missing,
            missing,
            teuk_func,
            UNIT_TEUKOLSKY_TRANS,
        )
    end

    if _is_negative_real_frequency(omega)
        positive_func = Teukolsky_radial(s, l, -m, a, -omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        return _conjugate_teukolsky_radial_function(positive_func, m, omega)
    end

    if _is_horizon_superradiance_frequency(a, m, omega) && boundary_condition == UP
        return _build_superradiance_threshold_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    end

    if boundary_condition == DOWN || boundary_condition == OUT
        Rin = Teukolsky_radial(s, l, m, a, omega, IN; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        Rup = Teukolsky_radial(s, l, m, a, omega, UP; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        return boundary_condition == DOWN ? _combine_down(Rin, Rup) : _combine_out(Rin, Rup)
    end

    return _build_teukolsky_function_fixed(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
end

function GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition::BoundaryCondition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, use_gsn_asymptotic_patches=false, gsn_horizon_delta_r_max=_GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min=_GSN_INFINITY_PHASE_MIN)
    if _is_negative_real_frequency(omega)
        positive_func = GSN_radial(s, l, -m, a, -omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
        return _conjugate_gsn_radial_function(positive_func, m, omega)
    end

    if _is_horizon_superradiance_frequency(a, m, omega) && boundary_condition == UP
        teuk_func, metadata = _build_superradiance_threshold_teukolsky_function(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, return_metadata=true)
        defaults = matching_controls(s, l, m, a, omega; boundary_condition = boundary_condition)
        controls = _radial_controls(defaults; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        isem_parameters = _isem_parameter_summary(; xm=controls.xm, rhom=controls.rhom, N=controls.N, tol=controls.tol, sfe=controls.sfe, lfe=controls.lfe, TSinInf=controls.TSinInf, TSoutInf=controls.TSoutInf, TSinHor=controls.TSinHor, TSoutHor=controls.TSoutHor, metadata=metadata)
        return _build_gsn_function_from_teukolsky(teuk_func; isem_parameters=isem_parameters, use_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    end

    if boundary_condition == DOWN || boundary_condition == OUT
        teuk_func = Teukolsky_radial(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        defaults = matching_controls(s, l, m, a, omega; boundary_condition = boundary_condition)
        controls = _radial_controls(defaults; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        isem_parameters = _isem_parameter_summary(; xm=controls.xm, rhom=controls.rhom, N=controls.N, tol=controls.tol, sfe=controls.sfe, lfe=controls.lfe, TSinInf=controls.TSinInf, TSoutInf=controls.TSoutInf, TSinHor=controls.TSinHor, TSoutHor=controls.TSoutHor)
        return _build_gsn_function_from_teukolsky(teuk_func; isem_parameters=isem_parameters, use_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    end
    teuk_func, metadata = _build_teukolsky_function_fixed(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, return_metadata=true)
    defaults = matching_controls(s, l, m, a, omega; boundary_condition = boundary_condition)
    controls = _radial_controls(defaults; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    isem_parameters = _isem_parameter_summary(; xm=controls.xm, rhom=controls.rhom, N=controls.N, tol=controls.tol, sfe=controls.sfe, lfe=controls.lfe, TSinInf=controls.TSinInf, TSoutInf=controls.TSoutInf, TSinHor=controls.TSinHor, TSoutHor=controls.TSoutHor, metadata=metadata)
    return _build_gsn_function_from_teukolsky(teuk_func; isem_parameters=isem_parameters, use_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
end

function GSN_radial(s::Int, l::Int, m::Int, a, omega; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, use_gsn_asymptotic_patches=false, gsn_horizon_delta_r_max=_GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min=_GSN_INFINITY_PHASE_MIN)
    Xin = GSN_radial(s, l, m, a, omega, IN; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    Xup = GSN_radial(s, l, m, a, omega, UP; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    return Xin, Xup
end

function Teukolsky_radial(s::Int, l::Int, m::Int, a, omega; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    Rin = Teukolsky_radial(s, l, m, a, omega, IN; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    Rup = Teukolsky_radial(s, l, m, a, omega, UP; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    return Rin, Rup
end

function _build_y_solution(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    _validate_y_branch(s, boundary_condition)
    iszero(omega) && error("Y_radial currently requires nonzero omega.")
    isreal(omega) || error("Y_radial currently supports only real omega.")

    teuk_func = Teukolsky_radial(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    lambda = teuk_func.mode.lambda
    kappa = sqrt(1 - a^2)
    branch_transmission_amplitude = _y_branch_transmission_amplitude(s, m, a, omega, lambda)
    if s == -2 && boundary_condition == IN
        func1, func2, func3, func4, func5 = Matching.TeukolskyTransformation.P_to_Y_coeffs_minus2_in(m, a, omega, lambda, branch_transmission_amplitude)
    elseif s == 2 && boundary_condition == UP
        func1, func2, func3, func4, func5 = Matching.TeukolskyTransformation.P_to_Y_coeffs_plus2_up(m, a, omega, lambda, branch_transmission_amplitude)
    else
        error("Unsupported Y branch.")
    end

    x_of_r = r -> Matching.TeukolskyTransformation.r_to_x(r, kappa)
    pullback_P = r -> teuk_func.P_solution(x_of_r(r))
    Xsoln = r -> begin
        P, Pp, _, _ = pullback_P(r)
        return func4(r) * P + func5(r) * Pp
    end

    function Yscalar(r)
        P, _, _, _ = pullback_P(r)
        return func1(r) * P
    end

    function Ysoln(r)
        P, Pp, _, error = pullback_P(r)
        Y = func1(r) * P
        Yp = func2(r) * P + func3(r) * Pp
        X = func4(r) * P + func5(r) * Pp
        return (Y, Yp, X, error)
    end

    mode = Mode(s, l, m, a, omega, lambda)
    return YRadialFunction(mode, boundary_condition, branch_transmission_amplitude, branch_transmission_amplitude * teuk_func.incidence_amplitude, branch_transmission_amplitude * teuk_func.reflection_amplitude, teuk_func.P_solution, teuk_func.Teukolsky_solution, Xsoln, Yscalar, Ysoln, UNIT_TEUKOLSKY_TRANS)
end

function _y_fallback_horizon_expansion_order(a)
    denom = max(1 - abs(float(a)), eps(Float64))
    return max(20, ceil(Int, 5 * max(log10(inv(denom)), 0)))
end

function _y_fallback_infinity_expansion_order(omega)
    ω = abs(float(real(omega)))
    if ω <= 0
        return 30
    end
    return max(20, ceil(Int, 10 * max(-log10(ω), 0) + 10))
end

function _build_y_solution_from_gsn(s, l, m, a, omega, boundary_condition; method="Riccati", tol=nothing)
    _validate_y_branch(s, boundary_condition)
    iszero(omega) && error("Y_radial currently requires nonzero omega.")
    isreal(omega) || error("Y_radial currently supports only real omega.")

    parent = parentmodule(@__MODULE__)
    X_func = getfield(parent, :GSN_radial)(
        s,
        l,
        m,
        a,
        omega,
        boundary_condition,
        getfield(parent, :_DEFAULT_rsin),
        getfield(parent, :_DEFAULT_rsout);
        method=method,
        tolerance = tol === nothing ? Solutions._DEFAULTTOLERANCE : tol,
        horizon_expansion_order = _y_fallback_horizon_expansion_order(a),
        infinity_expansion_order = _y_fallback_infinity_expansion_order(omega),
    )
    lambda = X_func.mode.lambda
    branch_transmission_amplitude = _y_branch_transmission_amplitude(s, m, a, omega, lambda)
    teukolsky_from_gsn_matrix = r -> Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(s, m, a, omega, lambda, r)
    Ysoln = Matching.TeukolskyTransformation.GSN_to_Y_solution_from_matrix(
        r -> rstar_from_r(a, r),
        teukolsky_from_gsn_matrix,
        X_func.GSN_solution,
        s,
        m,
        a,
        omega,
        lambda,
        branch_transmission_amplitude,
    )

    conv_incidence = s == -2 ?
        ConversionFactors.Binc(s, m, a, omega, lambda) :
        ConversionFactors.Cinc(s, m, a, omega, lambda)
    conv_reflection = s == -2 ?
        ConversionFactors.Bref(s, m, a, omega, lambda) :
        ConversionFactors.Cref(s, m, a, omega, lambda)
    incidence_amplitude = X_func.incidence_amplitude === missing ? missing : conv_incidence * X_func.incidence_amplitude
    reflection_amplitude = X_func.reflection_amplitude === missing ? missing : conv_reflection * X_func.reflection_amplitude
    Xsoln = r -> X_func.GSN_solution(rstar_from_r(a, r))[1]
    Yscalar = r -> Ysoln(r)[1]

    mode = Mode(s, l, m, a, omega, lambda)
    return YRadialFunction(mode, boundary_condition, branch_transmission_amplitude, incidence_amplitude, reflection_amplitude, missing, missing, Xsoln, Yscalar, Ysoln, UNIT_TEUKOLSKY_TRANS)
end

function _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
    saw_warn = Ref(false)
    logger = EarlyFilteredLogger(
        log -> begin
            if log.level == Logging.Warn
                saw_warn[] = true
            end
            return log.level >= Logging.Error
        end,
        current_logger(),
    )

    try
        sol = with_logger(logger) do
            riccati_build()
        end
        if saw_warn[]
            @info "Y_radial method = \"auto\" tried Riccati, received a warning, and switched to method = \"linear\"." context=context
            return linear_build()
        end
        return sol
    catch err
        @info "Y_radial method = \"auto\" tried Riccati, received an error, and switched to method = \"linear\"." context=context error=sprint(showerror, err)
        return linear_build()
    end
end

function _try_y_isem_then_legacy(context, isem_build, riccati_build, linear_build)
    saw_warn = Ref(false)
    logger = EarlyFilteredLogger(
        log -> begin
            if log.level == Logging.Warn
                saw_warn[] = true
            end
            return log.level >= Logging.Error
        end,
        current_logger(),
    )

    try
        sol = with_logger(logger) do
            isem_build()
        end
        if saw_warn[]
            @info "Y_radial method = \"auto\" tried ISEM, received an ISEM warning, and switched to legacy auto (Riccati, then linear). Use method = \"ISEM\" to force ISEM." context=context
            return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
        end
        return sol
    catch err
        @info "Y_radial method = \"auto\" tried ISEM, received an ISEM error, and switched to legacy auto (Riccati, then linear). Use method = \"ISEM\" to force ISEM." context=context error=sprint(showerror, err)
        return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
    end
end

function _y_compare_values!(vals::Vector, y::YRadialFunction)
    for value in (y.transmission_amplitude, y.incidence_amplitude, y.reflection_amplitude)
        value === missing || push!(vals, value)
    end
    return vals
end

function _y_relative_difference(a::YRadialFunction, b::YRadialFunction)
    av = _y_compare_values!(Any[], a)
    bv = _y_compare_values!(Any[], b)
    n = min(length(av), length(bv))
    n == 0 && return Inf
    diffs = Float64[]
    for i in 1:n
        ai = av[i]
        bi = bv[i]
        (ai === missing || bi === missing) && continue
        denom = max(abs(ai), abs(bi), 1e-300)
        push!(diffs, abs(ai - bi) / denom)
    end
    return isempty(diffs) ? Inf : maximum(diffs)
end

function _try_y_high_spin_isem_sanity_then_legacy(
    context,
    selector_build,
    low_n_build,
    mid_n_build,
    riccati_build,
    linear_build;
    sanity_tol = 1e-8,
)
    function try_build_with_warning_flag(build)
        saw_warn = Ref(false)
        logger = EarlyFilteredLogger(
            log -> begin
                if log.level == Logging.Warn
                    saw_warn[] = true
                end
                return log.level >= Logging.Error
            end,
            current_logger(),
        )
        sol = with_logger(logger) do
            build()
        end
        return sol, saw_warn[]
    end

    selector_sol, selector_warn = try
        try_build_with_warning_flag(selector_build)
    catch err
        @info "Y_radial method = \"auto\" high-spin ISEM selector build failed and switched to legacy auto (Riccati, then linear)." context=context error=sprint(showerror, err)
        return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
    end
    if selector_warn
        @info "Y_radial method = \"auto\" high-spin ISEM selector build warned and switched to legacy auto (Riccati, then linear)." context=context
        return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
    end

    low_sol, low_warn = try
        try_build_with_warning_flag(low_n_build)
    catch err
        @info "Y_radial method = \"auto\" high-spin low-N ISEM sanity build failed and switched to legacy auto (Riccati, then linear)." context=context error=sprint(showerror, err)
        return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
    end
    if low_warn
        @debug "Y_radial method = \"auto\" high-spin low-N ISEM sanity build warned and switched to legacy auto (Riccati, then linear)." context=context low_N=12
        return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
    end

    selector_low_diff = _y_relative_difference(selector_sol, low_sol)
    if selector_low_diff <= sanity_tol
        return selector_sol
    end

    mid_sol, mid_warn = try
        try_build_with_warning_flag(mid_n_build)
    catch err
        @info "Y_radial method = \"auto\" high-spin mid-N ISEM sanity build failed and switched to legacy auto (Riccati, then linear)." context=context selector_low_diff=selector_low_diff error=sprint(showerror, err)
        return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
    end
    if mid_warn
        @debug "Y_radial method = \"auto\" high-spin mid-N ISEM sanity build warned and switched to legacy auto (Riccati, then linear)." context=context selector_low_diff=selector_low_diff mid_N=20
        return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
    end

    low_mid_diff = _y_relative_difference(low_sol, mid_sol)
    if low_mid_diff <= sanity_tol
        @debug "Y_radial method = \"auto\" high-spin ISEM selector/low-N sanity check selected low-N result." context=context selector_low_diff=selector_low_diff low_mid_diff=low_mid_diff low_N=12 mid_N=20
        return low_sol
    end

    @debug "Y_radial method = \"auto\" high-spin ISEM sanity check failed and switched to legacy auto (Riccati, then linear)." context=context selector_low_diff=selector_low_diff low_mid_diff=low_mid_diff low_N=12 mid_N=20
    return _try_y_legacy_riccati_then_linear(context, riccati_build, linear_build)
end

function Y_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition::BoundaryCondition; method="auto", xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    if _is_negative_real_frequency(omega)
        positive_func = Y_radial(s, l, -m, a, -omega, boundary_condition; method=method, xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        return _conjugate_y_radial_function(positive_func, m, omega)
    end

    if method == "ISEM"
        return _build_y_solution(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    elseif method == "Riccati" || method == "linear"
        return _build_y_solution_from_gsn(s, l, m, a, omega, boundary_condition; method=method, tol=tol)
    elseif method == "auto"
        context = (function_name = "Y_radial", s = s, l = l, m = m, a = a, omega = omega, boundary_condition = boundary_condition)
        if N === nothing && abs(a) > 0.95
            return _try_y_high_spin_isem_sanity_then_legacy(
                context,
                () -> _build_y_solution(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor),
                () -> _build_y_solution(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=12, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor),
                () -> _build_y_solution(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=20, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor),
                () -> _build_y_solution_from_gsn(s, l, m, a, omega, boundary_condition; method="Riccati", tol=tol),
                () -> _build_y_solution_from_gsn(s, l, m, a, omega, boundary_condition; method="linear", tol=tol),
            )
        end
        return _try_y_isem_then_legacy(
            context,
            () -> _build_y_solution(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor),
            () -> _build_y_solution_from_gsn(s, l, m, a, omega, boundary_condition; method="Riccati", tol=tol),
            () -> _build_y_solution_from_gsn(s, l, m, a, omega, boundary_condition; method="linear", tol=tol),
        )
    else
        error("Method must be 'auto', 'ISEM', 'Riccati', or 'linear'")
    end
end

(y_func::YRadialFunction)(coord) = y_func.Y_scalar_solution(coord)

end
