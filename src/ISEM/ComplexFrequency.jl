module ComplexFrequency

using Logging: NullLogger, with_logger
using ...Coordinates: r_from_rstar, rstar_from_r
using ...ConversionFactors
using ...Solutions: Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix
using ..Parameters:
    _require_binary64_input,
    gsn_parameters,
    gsn_controls
using ..CoefficientTables: gsn_coefficients
using ..ComplexRational:
    ComplexRationalEvaluator,
    complex_rational_build,
    complex_rational_state_r,
    complex_up_horizon_in_candidate,
    complex_up_initial_match_candidates,
    complex_offpole_match_plateau_candidate
import ..Matching:
    ConjugatedRoute,
    amplitude_errors,
    radial_route,
    match_route,
    route_state,
    route_patch_count
import ..MSTInfinity:
    mst_physical_amplitudes,
    mst_physical_plan,
    mst_physical_state,
    mst_nia_jump,
    nia_strength_fast
import ..RadialFunctions:
    teukolsky_horizon_pair,
    TeukolskyRadialEvaluator,
    TeukolskySolution,
    _teukolsky_converter,
    dx_drstar,
    gsn_radial_function,
    teukolsky_radial_function,
    route_teukolsky_solution,
    y_branch_supported,
    y_radial_function,
    route_y_solution

export ComplexParameters, ComplexControls
export ComplexContour, ComplexRoutePlan, ComplexRoute
export ComplexMSTEvaluator
export complex_parameters, complex_controls
export complex_contour_parameters, complex_x_from_rho, complex_r_from_rho
export complex_endpoint_states
export complex_route, complex_route_plan_of, complex_route_patch_count
export complex_amplitudes
export complex_nia_jump
export complex_route_state_r, complex_route_state_rstar
export complex_route_evaluate_r, complex_route_evaluate_rstar
export evaluate_complex_route_on_real_axis
export complex_gsn_solution_rstar, complex_gsn_radial_function
export complex_radial_route

struct ComplexParameters
    s::Int
    l::Int
    m::Int
    a::Float64
    omega::ComplexF64
    lambda::ComplexF64
    kappa::Float64
    rplus::Float64
    rminus::Float64
    omega_horizon::Float64
    epsilon::ComplexF64
    tau::ComplexF64
    p::ComplexF64
    z::ComplexF64
end

struct ComplexControls
    xm::Union{Nothing,Float64}
    rhom::Union{Nothing,Float64}
    N::Union{Nothing,Int}
    tolerance::Union{Nothing,Float64}
    small_frequency_expansion::Union{Nothing,Bool}
    large_frequency_expansion::Union{Nothing,Bool}
end

struct ComplexContour
    omega::ComplexF64
    p::ComplexF64
    z::ComplexF64
    dxdrho::ComplexF64
    infinity_direction::ComplexF64
    horizon_direction::ComplexF64
    rho_match::Float64
    ray_angle::Float64
end

struct ComplexRoutePlan{C,M,K,A}
    contour::C
    matching::M
    controls::K
    amplitudes::A
    patch_count::Int
end

struct ComplexStateEvaluator{P,F,C,T}
    params::P
    p_solution::F
    coefficients::C
    scale::T
end

struct ComplexSFEvaluator{R,P}
    route::R
    params::P
end

struct ComplexMSTEvaluator{P,M}
    params::P
    plan::M
end

# A complex-frequency solution: mode parameters, requested controls, state
# evaluator on the real axis, backend metadata, route plan, and the GSN and
# Teukolsky amplitudes in unit-GSN-transmission normalization.
struct ComplexRoute{P,C,F,E,M,Q,T}
    branch::Symbol
    params::P
    controls::C
    p_solution::F
    state_evaluator::E
    metadata::M
    plan::Q
    transmission::T
    incidence::T
    reflection::T
    teukolsky_transmission::T
    teukolsky_incidence::T
    teukolsky_reflection::T
end

_finite_complex(value) = isfinite(real(value)) && isfinite(imag(value))
_optional_float(value) = value === nothing ? nothing : Float64(value)
_optional_int(value) = value === nothing ? nothing : Int(value)
_optional_bool(value) = value === nothing ? nothing : Bool(value)

const _RATIONAL_STEP_LIMIT_ERROR =
    "two-ray coordinate path exceeded step limit."
const _RATIONAL_ORDINARY_Q0_ERROR =
    "noninvertible q0 in ISEM GSN ordinary P/Q recurrence"
const _RATIONAL_PATCH_LIMIT_ERROR =
    "complex rational patch limit reached"
const _RATIONAL_INFINITY_STEP_ERROR =
    "infinity path propagation failed:"
const _RATIONAL_INFINITY_CERTIFICATE_ERROR =
    "no certified two-ray infinity path;"

# Symbol naming the rational build failure that allows the P-equation
# fallback, or `nothing` for any other error.
@inline function _rational_build_fallback_trigger(error)
    error isa ErrorException || return nothing
    error.msg == _RATIONAL_STEP_LIMIT_ERROR &&
        return :two_ray_coordinate_step_limit
    occursin(_RATIONAL_ORDINARY_Q0_ERROR, error.msg) &&
        return :ordinary_q0_degeneracy
    occursin(_RATIONAL_PATCH_LIMIT_ERROR, error.msg) &&
        return :complex_rational_patch_limit
    startswith(error.msg, _RATIONAL_INFINITY_STEP_ERROR) &&
        return :infinity_path_step_failure
    startswith(error.msg, _RATIONAL_INFINITY_CERTIFICATE_ERROR) &&
        return :infinity_path_certificate_failure
    return nothing
end

"""
    complex_controls(; xm, rhom, N, tol, small_frequency_expansion, large_frequency_expansion)

Optional controls for `complex_route`; `nothing` leaves a control to
the backend. Requires `N >= 10` and positive `rhom` and `tol` when given.
"""
function complex_controls(;
    xm=nothing,
    rhom=nothing,
    N=nothing,
    tol=nothing,
    small_frequency_expansion=nothing,
    large_frequency_expansion=nothing,
)
    controls = ComplexControls(
        _optional_float(xm),
        _optional_float(rhom),
        _optional_int(N),
        _optional_float(tol),
        _optional_bool(small_frequency_expansion),
        _optional_bool(large_frequency_expansion),
    )
    controls.N !== nothing && controls.N < 10 &&
        throw(ArgumentError("complex-frequency expansion order N must be at least 10."))
    controls.rhom !== nothing && controls.rhom <= 0 &&
        throw(ArgumentError("complex-frequency contour matching radius rhom must be positive."))
    controls.tolerance !== nothing && controls.tolerance <= 0 &&
        throw(ArgumentError("complex-frequency tolerance must be positive."))
    return controls
end

function _validate_mode(s::Int, l::Int, m::Int, a, omega)
    s in -2:2 || throw(DomainError(s, "ISEM complex GSN supports s = -2,-1,0,1,2."))
    l >= max(abs(s), abs(m)) ||
        throw(DomainError((s, l, m), "Require l >= max(abs(s), abs(m))."))
    _require_binary64_input("a", a)
    isreal(a) || throw(ArgumentError("Kerr spin a must be real."))
    af = Float64(real(a))
    isfinite(af) && abs(af) < 1 ||
        throw(DomainError(a, "The ISEM complex backend requires finite abs(a) < 1."))
    _require_binary64_input("omega", omega)
    omegac = ComplexF64(omega)
    _finite_complex(omegac) || throw(ArgumentError("omega must be finite."))
    iszero(omegac) && throw(DomainError(omega, "omega = 0 is excluded from the complex-frequency route."))
    iszero(imag(omegac)) &&
        throw(ArgumentError("complex_route requires nonzero imag(omega); use the real-frequency route otherwise."))
    return af, omegac
end

"""
    complex_parameters(s, l, m, a, omega, lambda)

Validated mode parameters and derived Kerr quantities (horizons, horizon
frequency, `p = omega - m*Omega_H`) for a complex-frequency route.
"""
function complex_parameters(s::Int, l::Int, m::Int, a, omega, lambda)
    af, omegac = _validate_mode(s, l, m, a, omega)
    _require_binary64_input("lambda", lambda)
    lambdac = ComplexF64(lambda)
    _finite_complex(lambdac) || throw(ArgumentError("lambda must be finite."))
    kappa = sqrt((1 - af) * (1 + af))
    rplus = 1 + kappa
    rminus = 1 - kappa
    omega_horizon = af / (2 * rplus)
    epsilon = 2 * omegac
    tau = (epsilon - m * af) / kappa
    p = omegac - m * omega_horizon
    z = omegac / abs(omegac)
    return ComplexParameters(
        s,
        l,
        m,
        af,
        omegac,
        lambdac,
        kappa,
        rplus,
        rminus,
        omega_horizon,
        epsilon,
        tau,
        p,
        z,
    )
end

"""
    complex_contour_parameters(params; rho_match=1.0)

Straight complex contour x = dxdrho*rho used by the P-equation backend, with
`dxdrho = -conj(omega)/|omega|` and matching radius `rho_match`.
"""
function complex_contour_parameters(params::ComplexParameters; rho_match=1.0)
    rho = Float64(rho_match)
    isfinite(rho) && rho > 0 ||
        throw(ArgumentError("complex contour matching radius must be finite and positive."))
    dxdrho = -inv(params.z)
    return ComplexContour(
        params.omega,
        params.p,
        params.z,
        dxdrho,
        dxdrho,
        dxdrho,
        rho,
        angle(dxdrho),
    )
end

complex_x_from_rho(contour::ComplexContour, rho::Real) = contour.dxdrho * rho

function complex_r_from_rho(params::ComplexParameters, contour::ComplexContour, rho::Real)
    x = complex_x_from_rho(contour, rho)
    return params.rplus - 2 * params.kappa * x
end

function _isem_module()
    return parentmodule(@__MODULE__)
end

_root_module() = parentmodule(_isem_module())

function _matching_transformation_module()
    return getfield(_isem_module(), :TeukolskyTransformation)
end

function _normalize_branch(branch)
    branch === :IN && return :IN
    branch === :UP && return :UP
    root = _root_module()
    branch == getfield(root, :IN) && return :IN
    branch == getfield(root, :UP) && return :UP
    throw(ArgumentError("ISEM complex route supports only IN and UP."))
end

function _boundary_condition(branch::Symbol)
    return getfield(_root_module(), branch)
end

_metadata_value(metadata, name::Symbol, default) =
    metadata !== missing && hasproperty(metadata, name) ? getproperty(metadata, name) : default

function _resolved_controls(metadata, requested::ComplexControls)
    return (
        xm = _metadata_value(metadata, :control_xm, requested.xm),
        rhom = _metadata_value(metadata, :rhom,
            _metadata_value(metadata, :control_rhom, requested.rhom)),
        N = _metadata_value(metadata, :control_N, requested.N),
        tolerance = _metadata_value(metadata, :control_tol, requested.tolerance),
        small_frequency_expansion = _metadata_value(metadata, :control_small_frequency_expansion, requested.small_frequency_expansion),
        large_frequency_expansion = _metadata_value(metadata, :control_large_frequency_expansion, requested.large_frequency_expansion),
        rescue_source = _metadata_value(metadata, :control_rescue_source, missing),
    )
end

function _patch_count(metadata)
    names = (
        :horizon_in_patch_count,
        :horizon_out_patch_count,
        :infinity_in_patch_count,
        :infinity_out_patch_count,
    )
    return sum(Int(_metadata_value(metadata, name, 0)) for name in names)
end

function _amplitude_factors(params::ComplexParameters, branch::Symbol)
    p = params
    if branch == :IN
        return (
            transmission = ConversionFactors.Btrans(p.s, p.m, p.a, p.omega, p.lambda),
            incidence = ConversionFactors.Binc(p.s, p.m, p.a, p.omega, p.lambda),
            reflection = ConversionFactors.Bref(p.s, p.m, p.a, p.omega, p.lambda),
        )
    end
    return (
        transmission = ConversionFactors.Ctrans(p.s, p.m, p.a, p.omega, p.lambda),
        incidence = ConversionFactors.Cinc(p.s, p.m, p.a, p.omega, p.lambda),
        reflection = ConversionFactors.Cref(p.s, p.m, p.a, p.omega, p.lambda),
    )
end

function _p_to_gsn_coefficients(params::ComplexParameters)
    p = params
    matrix_function = r -> Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(
        p.s,
        p.m,
        p.a,
        p.omega,
        p.lambda,
        r,
    )
    tt = _matching_transformation_module()
    return getfield(tt, :P_to_GSN_coefficients_from_matrix)(
        matrix_function,
        p.s,
        p.m,
        p.a,
        p.omega,
        p.lambda,
    )
end

function _state_at_r(evaluator::ComplexStateEvaluator, r::Real)
    p = evaluator.params
    rf = Float64(r)
    isfinite(rf) && rf > p.rplus ||
        throw(DomainError(r, "real-axis ISEM complex evaluation requires r > r_plus."))
    x = (p.rplus - rf) / (2 * p.kappa)
    P, Px, _, error = evaluator.p_solution(x)
    A0, A1, B0, B1 = evaluator.coefficients(rf)
    X = evaluator.scale * (A0 * P + A1 * Px)
    dXdrstar = evaluator.scale * (B0 * P + B1 * Px)
    return (X=ComplexF64(X), dXdrstar=ComplexF64(dXdrstar), error=Float64(error))
end

function _state_at_r(evaluator::ComplexRationalEvaluator, r::Real)
    return complex_rational_state_r(evaluator, r)
end

function _state_at_r(evaluator::ComplexSFEvaluator, r::Real)
    p = evaluator.params
    rf = Float64(r)
    isfinite(rf) && rf > p.rplus ||
        throw(DomainError(
            r, "real-axis ISEM complex evaluation requires r > r_plus."))
    x = (rf - p.rplus) / (rf - p.rminus)
    state = route_state(evaluator.route, x)
    dxdrstar = p.kappa * (x - 1.0)^2 * x /
        (p.kappa + (x - 1.0)^2 + 2 * p.kappa^2 * x -
         p.kappa * x^2)
    return (
        X=ComplexF64(state.X),
        dXdrstar=ComplexF64(state.dXdx * dxdrstar),
        error=0.0,
    )
end

function _state_at_r(evaluator::ComplexMSTEvaluator, r::Real)
    p = evaluator.params
    rf = Float64(r)
    isfinite(rf) && rf > p.rplus ||
        throw(DomainError(
            r, "real-axis ISEM complex evaluation requires r > r_plus."))
    x = (rf - p.rplus) / (rf - p.rminus)
    state = mst_physical_state(evaluator.plan, x)
    dxdrstar = p.kappa * (x - 1.0)^2 * x /
        (p.kappa + (x - 1.0)^2 + 2 * p.kappa^2 * x -
         p.kappa * x^2)
    return (
        X=ComplexF64(state.X),
        dXdrstar=ComplexF64(state.dXdx * dxdrstar),
        error=Float64(state.estimated_relerr),
    )
end

function (evaluator::ComplexStateEvaluator)(rstar::Real)
    r = r_from_rstar(evaluator.params.a, Float64(rstar))
    state = _state_at_r(evaluator, r)
    return state.X, state.dXdrstar, state.error
end

function (evaluator::ComplexSFEvaluator)(rstar::Real)
    r = r_from_rstar(evaluator.params.a, Float64(rstar))
    state = _state_at_r(evaluator, r)
    return state.X, state.dXdrstar, state.error
end

function (evaluator::ComplexMSTEvaluator)(rstar::Real)
    r = r_from_rstar(evaluator.params.a, Float64(rstar))
    state = _state_at_r(evaluator, r)
    return state.X, state.dXdrstar, state.error
end

function _state_vectors(evaluator, rstars::AbstractVector{<:Real})
    states = map(evaluator, rstars)
    return (
        getindex.(states, 1),
        getindex.(states, 2),
        getindex.(states, 3),
    )
end

(evaluator::ComplexStateEvaluator)(
    rstars::AbstractVector{<:Real},
) = _state_vectors(evaluator, rstars)

(evaluator::ComplexSFEvaluator)(
    rstars::AbstractVector{<:Real},
) = _state_vectors(evaluator, rstars)

(evaluator::ComplexMSTEvaluator)(
    rstars::AbstractVector{<:Real},
) = _state_vectors(evaluator, rstars)

# Route from the P-equation builder (backend=:p_equation); the state is
# normalized to unit GSN transmission.
function _complex_p_route(
    s::Int,
    l::Int,
    m::Int,
    a,
    omega,
    branch;
    controls::Union{Nothing,ComplexControls}=nothing,
    xm=nothing,
    rhom=nothing,
    N=nothing,
    tol=nothing,
    small_frequency_expansion=nothing,
    large_frequency_expansion=nothing,
)
    af, omegac = _validate_mode(s, l, m, a, omega)
    branch_symbol = _normalize_branch(branch)
    requested = controls === nothing ? complex_controls(
        xm=xm,
        rhom=rhom,
        N=N,
        tol=tol,
        small_frequency_expansion=small_frequency_expansion,
        large_frequency_expansion=large_frequency_expansion,
    ) : controls
    controls !== nothing && any(value !== nothing for value in
        (xm, rhom, N, tol, small_frequency_expansion, large_frequency_expansion)) &&
        throw(ArgumentError("pass either controls or individual complex-route control keywords, not both."))

    requested.small_frequency_expansion === true && throw(ArgumentError(
        "small_frequency_expansion=true is not available for backend=:p_equation."))
    requested.large_frequency_expansion === true && throw(ArgumentError(
        "large_frequency_expansion=true is not available for backend=:p_equation."))

    # PEquation is included after this module, so it is looked up at run time.
    builder = getfield(getfield(_isem_module(), :PEquation), :build_teukolsky_function)
    teukolsky, metadata = builder(
        s,
        l,
        m,
        af,
        omegac,
        _boundary_condition(branch_symbol);
        xm=requested.xm,
        rhom=requested.rhom,
        N=requested.N,
        tol=requested.tolerance,
        return_metadata=true,
    )
    teukolsky.P_solution === missing &&
        error("complex-frequency P-equation construction did not return a real-axis P solution.")

    params = complex_parameters(s, l, m, af, omegac, teukolsky.mode.lambda)
    factors = _amplitude_factors(params, branch_symbol)
    all(_finite_complex, (factors.transmission, factors.incidence, factors.reflection)) ||
        error("nonfinite Teukolsky-to-GSN asymptotic conversion factor.")
    raw_teukolsky_transmission =
        ComplexF64(teukolsky.transmission_amplitude)
    raw_teukolsky_incidence = ComplexF64(teukolsky.incidence_amplitude)
    raw_teukolsky_reflection = ComplexF64(teukolsky.reflection_amplitude)

    raw_gsn_transmission =
        raw_teukolsky_transmission / factors.transmission
    state_scale = inv(raw_gsn_transmission)
    incidence = ComplexF64(
        (raw_teukolsky_incidence / factors.incidence) * state_scale)
    reflection = ComplexF64(
        (raw_teukolsky_reflection / factors.reflection) * state_scale)
    transmission = ComplexF64(1)
    teukolsky_transmission =
        ComplexF64(raw_teukolsky_transmission * state_scale)
    teukolsky_incidence =
        ComplexF64(raw_teukolsky_incidence * state_scale)
    teukolsky_reflection =
        ComplexF64(raw_teukolsky_reflection * state_scale)
    all(_finite_complex, (state_scale, incidence, reflection,
        teukolsky_transmission, teukolsky_incidence,
        teukolsky_reflection)) ||
        error("nonfinite GSN unit-transmission normalization in ISEM complex construction.")

    coefficient_evaluator = _p_to_gsn_coefficients(params)
    state_evaluator = ComplexStateEvaluator(
        params,
        teukolsky.P_solution,
        coefficient_evaluator,
        ComplexF64(state_scale),
    )
    resolved = _resolved_controls(metadata, requested)
    rho_match = resolved.rhom === nothing ? 1.0 : resolved.rhom
    contour = complex_contour_parameters(params; rho_match=rho_match)
    matching = (
        requested_x = _metadata_value(metadata, :xm, missing),
        allowed_min_x = _metadata_value(metadata, :xm_min,
            _metadata_value(metadata, :xmin, missing)),
        allowed_max_x = _metadata_value(metadata, :xm_max,
            _metadata_value(metadata, :xmax, missing)),
        candidate_x = _metadata_value(metadata, :xm_match, missing),
        selected_x = _metadata_value(metadata, :xsplit, missing),
        split_mismatch = _metadata_value(metadata, :split_mismatch, missing),
    )
    amplitude_summary = (
        teukolsky = (
            transmission=teukolsky_transmission,
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        ),
        conversion_factors=factors,
        raw_gsn_transmission=ComplexF64(raw_gsn_transmission),
        gsn=(transmission=transmission, incidence=incidence, reflection=reflection),
    )
    plan = ComplexRoutePlan(
        contour,
        matching,
        resolved,
        amplitude_summary,
        _patch_count(metadata),
    )
    return ComplexRoute(
        branch_symbol,
        params,
        requested,
        teukolsky.P_solution,
        state_evaluator,
        metadata,
        plan,
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
end

function _tag_p_build_fallback(
    route::ComplexRoute,
    fallback_trigger::Symbol,
    rational_build_error::String,
)
    metadata = merge(route.metadata, (;
        backend=:direct_gsn_p_build_fallback,
        fallback_from=:direct_gsn_two_ray_rational,
        fallback_trigger,
        rational_build_error,
    ))
    return ComplexRoute(
        route.branch,
        route.params,
        route.controls,
        route.p_solution,
        route.state_evaluator,
        metadata,
        route.plan,
        route.transmission,
        route.incidence,
        route.reflection,
        route.teukolsky_transmission,
        route.teukolsky_incidence,
        route.teukolsky_reflection,
    )
end

@inline function _physical_ray_direction(q)
    iszero(q) && throw(DomainError(q, "complex physical ray frequency cannot be zero."))
    sign = real(q) >= 0 ? 1.0 : -1.0
    return ComplexF64(sign * cis(-angle(q)))
end

@inline _principal_small_frequency_expansion_axis(omega) =
    iszero(real(omega)) && imag(omega) < 0 && abs(omega) <= 0.1

function _principal_small_frequency_expansion_match(params, branch, requested)
    order = requested.N === nothing ? 40 : requested.N
    match_x = requested.xm === nothing ? 0.5 : requested.xm
    tolerance = requested.tolerance === nothing ? 1.0e-14 :
        requested.tolerance
    direct_params = gsn_parameters(
        params.s,
        params.l,
        params.m,
        params.a,
        params.omega;
        lambda=params.lambda,
    )
    controls = gsn_controls(
        direct_params;
        horizon_order=order,
        ordinary_order=order,
        infinity_order=order,
        xm=match_x,
        tol=tolerance,
        small_frequency_expansion=true,
        large_frequency_expansion=false,
    )
    coefficients = gsn_coefficients(direct_params; controls)
    return match_route(coefficients, branch; controls)
end

# Replace the incidence and reflection amplitudes of a route on the principal
# small-frequency-expansion (principal MST) axis by the matched GSN amplitudes, unless small_frequency_expansion=false was requested.
function _principal_small_frequency_expansion_overlay(route, requested)
    _principal_small_frequency_expansion_axis(route.params.omega) &&
        requested.small_frequency_expansion !== false || return route
    matched = _principal_small_frequency_expansion_match(route.params, route.branch, requested)
    transmission = route.transmission
    incidence = ComplexF64(matched.incidence)
    reflection = ComplexF64(matched.reflection)
    factors = route.plan.amplitudes.conversion_factors
    teukolsky_transmission = ComplexF64(
        factors.transmission * transmission)
    teukolsky_incidence = ComplexF64(factors.incidence * incidence)
    teukolsky_reflection = ComplexF64(factors.reflection * reflection)
    metadata = merge(route.metadata, (;
        control_small_frequency_expansion=true,
        amplitude_backend=:principal_mst_logscaled,
        amplitude_patch_count=route_patch_count(matched),
    ))
    controls = merge(route.plan.controls, (;
        small_frequency_expansion=true,
        rescue_source=:principal_mst_logscaled,
    ))
    amplitudes = (;
        teukolsky=(
            transmission=teukolsky_transmission,
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        ),
        conversion_factors=factors,
        raw_gsn_transmission=transmission,
        gsn=(; transmission, incidence, reflection),
    )
    plan = ComplexRoutePlan(
        route.plan.contour,
        route.plan.matching,
        controls,
        amplitudes,
        route.plan.patch_count + route_patch_count(matched),
    )
    return ComplexRoute(
        route.branch,
        route.params,
        route.controls,
        route.p_solution,
        route.state_evaluator,
        metadata,
        plan,
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
end

# GSN match with small-frequency-expansion (principal MST) endpoint solutions, used on the
# negative imaginary axis for |omega| <= 0.1.
function _complex_principal_small_frequency_expansion_route(
    s,
    l,
    m,
    a,
    omega,
    branch,
    requested;
    lambda=nothing,
    nu=nothing,
)
    af, omegac = _validate_mode(s, l, m, a, omega)
    branch_symbol = _normalize_branch(branch)
    order = requested.N === nothing ? 40 : requested.N
    match_x = requested.xm === nothing ? 0.5 : requested.xm
    tolerance = requested.tolerance === nothing ? 1.0e-14 :
        requested.tolerance
    params = gsn_parameters(
        s, l, m, af, omegac; lambda, nu)
    controls = gsn_controls(
        params;
        horizon_order=order,
        ordinary_order=order,
        infinity_order=order,
        xm=match_x,
        tol=tolerance,
        small_frequency_expansion=true,
        large_frequency_expansion=false,
    )
    coefficients = gsn_coefficients(params; controls)
    matched = match_route(coefficients, branch_symbol; controls)
    complex_params = complex_parameters(
        s, l, m, af, omegac, params.lambda)
    factors = _amplitude_factors(complex_params, branch_symbol)
    transmission = ComplexF64(1)
    incidence = ComplexF64(matched.incidence)
    reflection = ComplexF64(matched.reflection)
    teukolsky_transmission = ComplexF64(factors.transmission)
    teukolsky_incidence = ComplexF64(factors.incidence * incidence)
    teukolsky_reflection = ComplexF64(factors.reflection * reflection)
    all(_finite_complex, (
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )) || error("nonfinite principal-MST complex amplitude.")
    patch_count = route_patch_count(matched)
    contour = (;
        model=:principal_mst_logscaled_real_axis,
        omega=complex_params.omega,
        p=complex_params.p,
        horizon_direction=ComplexF64(1),
        infinity_direction=ComplexF64(1),
        horizon_rho=NaN,
        infinity_rho=NaN,
        rho_match=NaN,
        ray_angle=0.0,
    )
    resolved = (;
        xm=match_x,
        rhom=nothing,
        N=order,
        tolerance,
        small_frequency_expansion=true,
        large_frequency_expansion=false,
        rescue_source=:principal_mst_logscaled,
    )
    matching = (;
        requested_x=requested.xm,
        allowed_min_x=0.0,
        allowed_max_x=1.0,
        candidate_x=match_x,
        selected_x=match_x,
        split_mismatch=0.0,
        condition=NaN,
        policy=:principal_mst_logscaled,
    )
    amplitudes = (;
        teukolsky=(
            transmission=teukolsky_transmission,
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        ),
        conversion_factors=factors,
        raw_gsn_transmission=transmission,
        gsn=(; transmission, incidence, reflection),
    )
    plan = ComplexRoutePlan(
        contour, matching, resolved, amplitudes, patch_count)
    metadata = (;
        backend=:principal_mst_logscaled,
        match_policy=:principal_mst_logscaled,
        match_x,
        xsplit=match_x,
        xm_match=match_x,
        split_mismatch=0.0,
        matching_condition=NaN,
        abel_error=0.0,
        patch_count,
        control_xm=match_x,
        control_rhom=nothing,
        control_N=order,
        control_tol=tolerance,
        control_small_frequency_expansion=true,
        control_large_frequency_expansion=false,
        amplitude_backend=:principal_mst_logscaled,
        amplitude_patch_count=patch_count,
    )
    evaluator = ComplexSFEvaluator(matched, complex_params)
    return ComplexRoute(
        branch_symbol,
        complex_params,
        requested,
        nothing,
        evaluator,
        metadata,
        plan,
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
end

# Route from the complete MST representation (backend=:mst).
function _complex_mst_route(
    s::Int,
    l::Int,
    m::Int,
    a,
    omega,
    branch;
    controls::Union{Nothing,ComplexControls}=nothing,
    lambda=nothing,
    nu=nothing,
    xm=nothing,
    rhom=nothing,
    N=nothing,
    tol=nothing,
    small_frequency_expansion=nothing,
    large_frequency_expansion=nothing,
)
    af, omegac = _validate_mode(s, l, m, a, omega)
    branch_symbol = _normalize_branch(branch)
    requested = controls === nothing ? complex_controls(
        xm=xm,
        rhom=rhom,
        N=N,
        tol=tol,
        small_frequency_expansion=small_frequency_expansion,
        large_frequency_expansion=large_frequency_expansion,
    ) : controls
    controls !== nothing && any(value !== nothing for value in
        (xm, rhom, N, tol, small_frequency_expansion, large_frequency_expansion)) && throw(ArgumentError(
            "pass either controls or individual complex-route control " *
            "keywords, not both."))
    requested.xm === nothing || throw(ArgumentError(
        "xm is not used by backend=:mst."))
    requested.rhom === nothing || throw(ArgumentError(
        "rhom is not used by backend=:mst."))
    requested.small_frequency_expansion === true && throw(ArgumentError(
        "small_frequency_expansion=true is incompatible with backend=:mst."))
    requested.large_frequency_expansion === true && throw(ArgumentError(
        "large_frequency_expansion=true is incompatible with backend=:mst."))

    order = requested.N === nothing ? 40 : requested.N
    tolerance = requested.tolerance === nothing ? 1.0e-14 :
        requested.tolerance
    params = gsn_parameters(
        s, l, m, af, omegac; lambda, nu)
    coefficient_controls = gsn_controls(
        params;
        N=order,
        xm=0.5,
        tol=tolerance,
        small_frequency_expansion=false,
        large_frequency_expansion=false,
    )
    coefficients = gsn_coefficients(
        params; controls=coefficient_controls)
    mst = mst_physical_plan(coefficients, branch_symbol)
    mst_amplitudes = mst_physical_amplitudes(mst)
    anchor = mst_physical_state(mst, 0.5)

    complex_params = complex_parameters(
        s, l, m, af, omegac, params.lambda)
    factors = _amplitude_factors(complex_params, branch_symbol)
    transmission = ComplexF64(1)
    incidence = ComplexF64(mst_amplitudes.incidence)
    reflection = ComplexF64(mst_amplitudes.reflection)
    teukolsky_transmission = ComplexF64(factors.transmission)
    teukolsky_incidence = ComplexF64(
        factors.incidence * incidence)
    teukolsky_reflection = ComplexF64(
        factors.reflection * reflection)
    all(_finite_complex, (
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
        anchor.X,
        anchor.dXdx,
    )) || error("nonfinite complete-MST complex route.")

    contour = (;
        model=:physical_mst_real_axis,
        omega=complex_params.omega,
        p=complex_params.p,
        horizon_direction=ComplexF64(1),
        infinity_direction=ComplexF64(1),
        horizon_rho=NaN,
        infinity_rho=NaN,
        rho_match=NaN,
        ray_angle=0.0,
    )
    matching = (;
        requested_x=nothing,
        allowed_min_x=0.0,
        allowed_max_x=1.0,
        candidate_x=missing,
        selected_x=missing,
        split_mismatch=anchor.estimated_relerr,
        condition=anchor.condition,
        policy=:direct_mst,
    )
    resolved = (;
        xm=nothing,
        rhom=nothing,
        N=order,
        tolerance,
        small_frequency_expansion=false,
        large_frequency_expansion=false,
        rescue_source=:none,
    )
    amplitude_summary = (;
        teukolsky=(
            transmission=teukolsky_transmission,
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        ),
        conversion_factors=factors,
        raw_gsn_transmission=transmission,
        gsn=(; transmission, incidence, reflection),
        mst=mst_amplitudes,
    )
    route_plan = ComplexRoutePlan(
        contour,
        matching,
        resolved,
        amplitude_summary,
        0,
    )
    metadata = (;
        backend=:direct_gsn_full_mst,
        match_policy=:direct_mst,
        match_x=missing,
        xsplit=missing,
        xm_match=missing,
        matching_condition=anchor.condition,
        split_mismatch=anchor.estimated_relerr,
        abel_error=0.0,
        patch_count=0,
        control_xm=nothing,
        control_rhom=nothing,
        control_N=order,
        control_tol=tolerance,
        control_small_frequency_expansion=false,
        control_large_frequency_expansion=false,
        amplitude_backend=:direct_mst,
        amplitude_patch_count=0,
        mst_nu=mst_amplitudes.nu,
        mst_nu_offset=mst_amplitudes.nu_offset,
        mst_representation=mst_amplitudes.representation,
        mst_amplitude_certificate=mst_amplitudes.certificate,
        mst_amplitude_nmax=mst_amplitudes.nmax,
        mst_amplitude_check_nmax=mst_amplitudes.check_nmax,
        mst_amplitude_precision_bits=mst_amplitudes.precision_bits,
        mst_nu_residual=mst_amplitudes.nu_residual,
        mst_amplitude_truncation=
            mst_amplitudes.truncation_agreement,
        mst_amplitude_spread=
            mst_amplitudes.representation_spread,
        mst_amplitude_nearest=
            mst_amplitudes.nearest_agreement,
        mst_amplitude_condition=mst_amplitudes.max_condition,
        mst_raw_teukolsky_logs=mst_amplitudes.raw_teukolsky_logs,
        mst_anchor_x=0.5,
        mst_anchor_error=anchor.estimated_relerr,
        mst_anchor_residual=anchor.residual,
        mst_anchor_tail=anchor.tail,
        mst_anchor_condition=anchor.condition,
        mst_anchor_nmin=anchor.nmin,
        mst_anchor_nmax=anchor.nmax,
        mst_anchor_representation=anchor.representation,
    )
    evaluator = ComplexMSTEvaluator(complex_params, mst)
    return ComplexRoute(
        branch_symbol,
        complex_params,
        requested,
        nothing,
        evaluator,
        metadata,
        route_plan,
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
end

# Two-ray rational route (default backend) wrapped as a ComplexRoute;
# with `return_retry_cache=true` also returns the data reused by retries.
function _complex_rational_route(
    s::Int,
    l::Int,
    m::Int,
    a,
    omega,
    branch;
    controls::Union{Nothing,ComplexControls}=nothing,
    lambda=nothing,
    nu=nothing,
    xm=nothing,
    rhom=nothing,
    N=nothing,
    tol=nothing,
    small_frequency_expansion=nothing,
    large_frequency_expansion=nothing,
    return_retry_cache::Bool=false,
)
    af, omegac = _validate_mode(s, l, m, a, omega)
    branch_symbol = _normalize_branch(branch)
    requested = controls === nothing ? complex_controls(
        xm=xm,
        rhom=rhom,
        N=N,
        tol=tol,
        small_frequency_expansion=small_frequency_expansion,
        large_frequency_expansion=large_frequency_expansion,
    ) : controls
    controls !== nothing && any(value !== nothing for value in
        (xm, rhom, N, tol, small_frequency_expansion, large_frequency_expansion)) &&
        throw(ArgumentError(
            "pass either controls or individual complex-route control keywords, not both."))
    requested.rhom === nothing || throw(ArgumentError(
        "rhom is specific to backend=:p_equation; the rational backend " *
        "selects both physical-ray endpoints adaptively."))
    requested.small_frequency_expansion === true && !_principal_small_frequency_expansion_axis(omegac) &&
        throw(ArgumentError(
            "complex small_frequency_expansion=true is restricted to the negative " *
            "pure-imaginary principal-MST axis."))
    requested.large_frequency_expansion === true && throw(ArgumentError(
        "the complex two-ray rational backend does not use the real-frequency large-frequency-expansion switch."))

    rational_order = requested.N === nothing &&
        _principal_small_frequency_expansion_axis(omegac) ? 56 : requested.N
    build = complex_rational_build(
        s,
        l,
        m,
        af,
        omegac,
        branch_symbol;
        lambda,
        nu,
        N=rational_order,
        tol=requested.tolerance,
        xm=requested.xm,
    )
    params = complex_parameters(
        s, l, m, af, omegac, build.params.lambda)
    factors = _amplitude_factors(params, branch_symbol)
    transmission = ComplexF64(build.transmission)
    incidence = ComplexF64(build.incidence)
    reflection = ComplexF64(build.reflection)
    teukolsky_transmission = ComplexF64(factors.transmission * transmission)
    teukolsky_incidence = ComplexF64(factors.incidence * incidence)
    teukolsky_reflection = ComplexF64(factors.reflection * reflection)
    amplitudes_finite = all(_finite_complex, (
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    ))
    (amplitudes_finite || _principal_small_frequency_expansion_axis(omegac)) ||
        error("nonfinite asymptotic amplitude in rational complex construction.")

    horizon_rho = build.metadata.horizon_rho
    infinity_rho = build.metadata.infinity_rho
    rho_match = max(
        isfinite(horizon_rho) ? abs(horizon_rho) : 0.0,
        abs(infinity_rho),
    )
    contour = (;
        model=:two_physical_rstar_rays,
        omega=params.omega,
        p=params.p,
        horizon_direction=ComplexF64(build.metadata.horizon_direction),
        infinity_direction=_physical_ray_direction(params.omega),
        horizon_rho,
        infinity_rho,
        rho_match,
        ray_angle=angle(_physical_ray_direction(params.omega)),
    )
    resolved = (;
        xm=build.metadata.match_x,
        rhom=rho_match,
        N=build.controls.ordinary_order,
        tolerance=build.controls.tolerance,
        small_frequency_expansion=false,
        large_frequency_expansion=false,
        rescue_source=:none,
    )
    matching = (;
        requested_x=requested.xm,
        allowed_min_x=0.0,
        allowed_max_x=1.0,
        candidate_x=build.metadata.match_x,
        selected_x=build.metadata.match_x,
        split_mismatch=build.metadata.abel_error,
        condition=build.metadata.matching_condition,
        policy=build.metadata.match_policy,
    )
    amplitude_summary = (;
        teukolsky=(
            transmission=teukolsky_transmission,
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        ),
        conversion_factors=factors,
        raw_gsn_transmission=transmission,
        gsn=(
            transmission=transmission,
            incidence=incidence,
            reflection=reflection,
        ),
    )
    plan = ComplexRoutePlan(
        contour,
        matching,
        resolved,
        amplitude_summary,
        build.metadata.patch_count,
    )
    metadata = (;
        build.metadata...,
        xsplit=build.metadata.match_x,
        xm_match=build.metadata.match_x,
        split_mismatch=build.metadata.abel_error,
        control_xm=build.metadata.match_x,
        control_rhom=rho_match,
        control_N=build.controls.ordinary_order,
        control_tol=build.controls.tolerance,
        control_small_frequency_expansion=false,
        control_large_frequency_expansion=false,
        amplitude_backend=:complex_rational,
        amplitude_patch_count=0,
    )
    route = ComplexRoute(
        branch_symbol,
        params,
        requested,
        nothing,
        build.evaluator,
        metadata,
        plan,
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
    return return_retry_cache ? (; route, retry_cache=build.retry_cache) : route
end

# Trigger and acceptance thresholds of the consensus retries and fallbacks below.
const _P_CONSENSUS_KAPPA_MAX = 0.15
const _P_CONSENSUS_CANCELLATION_MIN = 1.0e5
const _P_CONSENSUS_ROUTE_RESIDUAL_MAX = 1.0e-5
const _P_CONSENSUS_RESIDUAL_MAX = 1.0e-8
const _P_CONSENSUS_RESIDUAL_RATIO_MAX = 1.0e-2
const _P_CONSENSUS_REFLECTION_ERROR_MAX = 1.0e-7
const _P_CONSENSUS_SPLIT_MISMATCH_MAX = 1.0e-10
const _P_CONSENSUS_SCATTERING_CANCELLATION_MIN = 1.0e4
const _P_CONSENSUS_SCATTERING_RESIDUAL_MIN = 3.0e-8
const _P_CONSENSUS_SCATTERING_RESIDUAL_MAX = 2.0e-4
const _P_CONSENSUS_SCATTERING_RESIDUAL_RATIO_MIN = 0.5
const _P_CONSENSUS_SCATTERING_RESIDUAL_RATIO_MAX = 2.0
const _P_CONSENSUS_SCATTERING_REFLECTION_ERROR_MAX = 1.0e-8
const _P_CONSENSUS_UP_RESIDUAL_MIN = 0.5
const _P_CONSENSUS_UP_ORDER_DELTA = 4
const _P_CONSENSUS_UP_AGREEMENT_MAX = 1.0e-10
const _P_CONSENSUS_UP_SEPARATION_MIN = 1.0e-10
const _P_CONSENSUS_UP_SPLIT_MISMATCH_MAX = 1.0e-10
const _P_CONSENSUS_POLE_DUAL_RESIDUAL_RATIO_MAX = 5.0e-2
const _P_CONSENSUS_POLE_DUAL_AGREEMENT_MAX = 5.0e-9
const _P_CONSENSUS_POLE_DUAL_SEPARATION_MIN = 1.0e-8
const _P_CONSENSUS_POLE_DUAL_REFLECTION_ERROR_MAX = 1.0e-12
const _P_CONSENSUS_POLE_DUAL_SPLIT_MISMATCH_MAX = 2.0e-10
const _P_CONTOUR_CONDITION_MIN = 1.0e12
const _P_CONTOUR_SPLIT_MISMATCH_MIN = 1.0e-10
const _P_CONTOUR_ORDER_AGREEMENT_MAX = 1.0e-10
const _P_CONTOUR_SEPARATION_MIN = 1.0e-10
const _P_CONTOUR_SPLIT_MISMATCH_MAX = 1.0e-10
const _P_CONTOUR_USABLE_ORDER_AGREEMENT_MAX = 1.0e-8
const _P_CONTOUR_USABLE_SEPARATION_MIN = 1.0e-8
const _P_CONTOUR_USABLE_SPLIT_MISMATCH_MAX = 2.0e-8
const _P_NIA_REAL_RATIO_MAX = 1.0e-8
const _P_NIA_MAX_ABS_OMEGA = 0.1
const _P_NIA_MAX_SPHEROIDICITY = 64eps(Float64)
const _P_NIA_ORDER = 40
const _P_NIA_ORDER_DELTA = 4
const _P_NIA_ORDER_AGREEMENT_MAX = 1.0e-10
const _P_NIA_STATE_AGREEMENT_MAX = 1.0e-10
const _P_NIA_SEPARATION_MIN = 1.0e-8
const _P_NIA_BUILD_ORDER_AGREEMENT_MAX = 1.0e-9
const _P_NIA_BUILD_STATE_AGREEMENT_MAX = 1.0e-7
const _P_NIA_NONROTATING_MAX_ABS_A = 64eps(Float64)
const _P_NIA_PARTNER_ORDER_AGREEMENT_MAX = 1.0e-9
const _NIA_SIDE_AGREEMENT_MAX = 1.0e-8
const _NIA_SIDE_OFFSET_FACTOR = 4.0
const _PRINCIPAL_P_FALLBACK_ORDER = 40
const _PRINCIPAL_P_FALLBACK_ORDER_DELTA = 8
const _PRINCIPAL_P_FALLBACK_AGREEMENT_MAX = 1.0e-10
const _PRINCIPAL_P_FALLBACK_SPLIT_MAX = 1.0e-10
const _PRINCIPAL_P_FALLBACK_STATE_POINTS = (-2.0, 0.0, 2.0)
const _UP_INITIAL_MATCH_RESIDUAL_MAX = 1.0e-8
const _UP_INITIAL_MATCH_FAR_MID_AGREEMENT_MAX = 1.5e-8
const _UP_INITIAL_MATCH_MID_INITIAL_AGREEMENT_MAX = 5.0e-9
const _UP_INITIAL_MATCH_FAR_INITIAL_AGREEMENT_MAX = 1.5e-8
const _UP_INITIAL_MATCH_REFLECTION_AGREEMENT_MAX = 3.0e-11
const _UP_INITIAL_MATCH_SEPARATION_MIN = 1.0e-8
const _UP_INITIAL_MATCH_SPLIT_MISMATCH_MAX = 1.0e-10
const _UP_HORIZON_IN_KAPPA_MAX = 0.05
const _UP_HORIZON_IN_CANCELLATION_MIN = 1.1
const _UP_HORIZON_IN_ORDER_DELTA = 4
const _UP_HORIZON_IN_LOWER_BACKOFF = 12
const _UP_HORIZON_IN_AGREEMENT_MAX = 1.0e-10
const _UP_HORIZON_IN_SEPARATION_MIN = 1.0e-10
const _UP_HORIZON_IN_SPLIT_MISMATCH_MAX = 1.0e-10
const _POLE_RECIPROCAL_KAPPA_MAX = 0.05
const _POLE_RECIPROCAL_IN_RESIDUAL_MIN = 1.0e-10
const _POLE_RECIPROCAL_IN_RESIDUAL_MAX = 1.0e-7
const _POLE_RECIPROCAL_UP_RESIDUAL_MAX = 1.0e-6
const _POLE_RECIPROCAL_ORDER_DELTA = 4
const _POLE_RECIPROCAL_ORDER_AGREEMENT_MAX = 5.0e-10
const _POLE_RECIPROCAL_STATE_ERROR_MAX = 2.0e-5
const _POLE_RECIPROCAL_SEPARATION_MIN = 1.0e-10
const _POLE_RECIPROCAL_SPLIT_MISMATCH_MAX = 1.0e-10
const _OFFPOLE_MATCH_PLATEAU_KAPPA_MAX = 0.15
const _OFFPOLE_MATCH_PLATEAU_UP_AGREEMENT_MAX = 2.5e-10
const _OFFPOLE_MATCH_PLATEAU_IN_AGREEMENT_MAX = 2.5e-11
const _OFFPOLE_MATCH_PLATEAU_UP_SEPARATION_RATIO_MIN = 10.0
const _OFFPOLE_MATCH_PLATEAU_IN_SEPARATION_RATIO_MIN = 11.0
const _OFFPOLE_MATCH_PLATEAU_SPLIT_MISMATCH_MAX = 1.0e-10
const _OFFPOLE_MATCH_PLATEAU_CONTINUITY_MAX = 7.5e-10
const _OFFPOLE_MATCH_PLATEAU_MEMBER_AGREEMENT_MAX = 5.0e-12
const _OFFPOLE_MATCH_PLATEAU_P_OVERRIDE_RATIO_MAX = 1.0e5

@inline function _raw_incidence_residual(route::ComplexRoute)
    scale = max(
        abs(route.teukolsky_transmission),
        abs(route.teukolsky_incidence),
        abs(route.teukolsky_reflection),
        floatmin(Float64),
    )
    return abs(route.teukolsky_incidence) / scale
end

@inline function _endpoint_amplitude_match(route::ComplexRoute)
    policy = _metadata_value(
        route.metadata, :amplitude_match_policy, missing)
    return policy !== missing &&
        policy in (:scaled_infinity_endpoint, :scaled_horizon_endpoint)
end

@inline function _p_consensus_pretrigger(route::ComplexRoute)
    _endpoint_amplitude_match(route) && return false
    metadata = route.metadata
    hasproperty(metadata, :coefficient1_cancellation) || return false
    hasproperty(metadata, :coefficient2_cancellation) || return false
    coefficient1_cancellation =
        Float64(metadata.coefficient1_cancellation)
    coefficient2_cancellation =
        Float64(metadata.coefficient2_cancellation)
    cancellation = max(
        coefficient1_cancellation,
        coefficient2_cancellation,
    )
    direct_residual = _raw_incidence_residual(route)
    pole_trigger = cancellation >= _P_CONSENSUS_CANCELLATION_MIN &&
        direct_residual <= _P_CONSENSUS_ROUTE_RESIDUAL_MAX
    scattering_trigger = route.branch == :IN &&
        coefficient1_cancellation >=
            _P_CONSENSUS_SCATTERING_CANCELLATION_MIN &&
        _P_CONSENSUS_SCATTERING_RESIDUAL_MIN <= direct_residual <=
            _P_CONSENSUS_SCATTERING_RESIDUAL_MAX
    return route.params.kappa <= _P_CONSENSUS_KAPPA_MAX &&
        (pole_trigger || scattering_trigger)
end

@inline function _p_up_dual_pretrigger(route::ComplexRoute)
    _endpoint_amplitude_match(route) && return false
    metadata = route.metadata
    hasproperty(metadata, :match_policy) || return false
    return route.params.s == 2 &&
        route.branch == :UP &&
        route.params.kappa <= _P_CONSENSUS_KAPPA_MAX &&
        metadata.match_policy == :up_condition_retry &&
        _raw_incidence_residual(route) >= _P_CONSENSUS_UP_RESIDUAL_MIN
end

@inline function _p_pole_dual_pretrigger(route::ComplexRoute)
    return route.params.s == 2 &&
        route.branch == :UP &&
        _p_consensus_pretrigger(route)
end

@inline function _p_contour_pretrigger(route::ComplexRoute)
    _endpoint_amplitude_match(route) && return false
    condition_value = _metadata_value(
        route.metadata, :matching_condition, missing)
    condition = condition_value === missing ? 0.0 : Float64(condition_value)
    split_mismatch = _route_split_mismatch(route)
    return (isfinite(condition) && condition >= _P_CONTOUR_CONDITION_MIN) ||
        (isfinite(split_mismatch) &&
            split_mismatch >= _P_CONTOUR_SPLIT_MISMATCH_MIN)
end

@inline function _p_nia_pretrigger(route::ComplexRoute)
    omega = route.params.omega
    imag(omega) < 0 || return false
    !iszero(real(omega)) || return false
    abs(real(omega)) <= _P_NIA_REAL_RATIO_MAX * abs(imag(omega)) ||
        return false
    abs(omega) <= _P_NIA_MAX_ABS_OMEGA || return false
    return abs(route.params.a * omega) <= _P_NIA_MAX_SPHEROIDICITY
end

@inline function _p_nia_build_pretrigger(s, a, omega, branch)
    imag(omega) < 0 || return false
    !iszero(real(omega)) || return false
    abs(real(omega)) <= _P_NIA_REAL_RATIO_MAX * abs(imag(omega)) ||
        return false
    abs(a * omega) <= _P_NIA_MAX_SPHEROIDICITY || return false
    branch_symbol = _normalize_branch(branch)
    return s == 2 || (s == 0 && branch_symbol == :IN)
end

@inline function _p_nia_nonrotating(route::ComplexRoute)
    return abs(route.params.a) <= _P_NIA_NONROTATING_MAX_ABS_A
end

@inline function _p_nia_partner_pretrigger(s, a, omega)
    s < 0 || return false
    iszero(a) || return false
    imag(omega) < 0 || return false
    !iszero(real(omega)) || return false
    abs(real(omega)) <= _P_NIA_REAL_RATIO_MAX * abs(imag(omega)) ||
        return false
    return abs(omega) > _P_NIA_MAX_ABS_OMEGA
end

@inline function _nia_side_pretrigger(route::ComplexRoute)
    omega = route.params.omega
    route.branch == :UP || return false
    abs(route.params.a) > _P_NIA_NONROTATING_MAX_ABS_A || return false
    real(omega) < 0 || return false
    imag(omega) < 0 || return false
    abs(real(omega)) <= _P_NIA_REAL_RATIO_MAX * abs(imag(omega)) ||
        return false
    return abs(omega) <= _P_NIA_MAX_ABS_OMEGA
end

@inline function _up_initial_match_pretrigger(route::ComplexRoute)
    metadata = route.metadata
    hasproperty(metadata, :match_policy) || return false
    return route.params.s == 2 &&
        route.branch == :UP &&
        metadata.match_policy == :up_condition_retry &&
        _p_consensus_pretrigger(route)
end

@inline function _route_pair_errors(
    first::ComplexRoute,
    second::ComplexRoute,
)
    return _amplitude_pair_errors(first, second)
end

@inline function _unit_pair_errors(
    first::ComplexRoute,
    second::ComplexRoute,
)
    values = (
        first.transmission,
        first.incidence,
        first.reflection,
        second.transmission,
        second.incidence,
        second.reflection,
    )
    all(_finite_complex, values) || return (Inf, Inf, Inf)
    iszero(first.transmission) && return (Inf, Inf, Inf)
    iszero(second.transmission) && return (Inf, Inf, Inf)
    alignment = first.transmission / second.transmission
    aligned_incidence = alignment * second.incidence
    aligned_reflection = alignment * second.reflection
    scale = max(
        abs(first.incidence),
        abs(first.reflection),
        abs(aligned_incidence),
        abs(aligned_reflection),
        floatmin(Float64),
    )
    incidence_error = abs(
        first.incidence - aligned_incidence) / scale
    reflection_error = abs(
        first.reflection - aligned_reflection) / scale
    return incidence_error, reflection_error,
        max(incidence_error, reflection_error)
end

@inline function _amplitude_pair_errors(first, second)
    values = (
        first.incidence,
        first.reflection,
        second.incidence,
        second.reflection,
    )
    all(_finite_complex, values) || return (Inf, Inf, Inf)
    scale = max(maximum(abs, values), floatmin(Float64))
    incidence_error = abs(first.incidence - second.incidence) / scale
    reflection_error = abs(first.reflection - second.reflection) / scale
    return incidence_error, reflection_error,
        max(incidence_error, reflection_error)
end

@inline function _route_split_mismatch(route::ComplexRoute)
    value = route.plan.matching.split_mismatch
    return value === missing ? Inf : Float64(value)
end

@inline function _state_pair_error(first, second)
    values = (first[1], first[2], second[1], second[2])
    all(_finite_complex, values) || return Inf
    scale = max(maximum(abs, values), floatmin(Float64))
    return max(
        abs(first[1] - second[1]),
        abs(first[2] - second[2]),
    ) / scale
end

function _principal_p_state_agreement(
    first::ComplexRoute,
    second::ComplexRoute,
)
    agreement = 0.0
    for rstar in _PRINCIPAL_P_FALLBACK_STATE_POINTS
        current = _state_pair_error(
            first.state_evaluator(rstar),
            second.state_evaluator(rstar),
        )
        isfinite(current) || return Inf
        agreement = max(agreement, current)
    end
    return agreement
end

function _principal_p_unit_state_agreement(
    first::ComplexRoute,
    second::ComplexRoute,
)
    _finite_complex(first.transmission) || return Inf
    _finite_complex(second.transmission) || return Inf
    iszero(first.transmission) && return Inf
    iszero(second.transmission) && return Inf
    alignment = first.transmission / second.transmission
    agreement = 0.0
    for rstar in _PRINCIPAL_P_FALLBACK_STATE_POINTS
        first_state = first.state_evaluator(rstar)
        second_state = second.state_evaluator(rstar)
        aligned_second = (
            alignment * second_state[1],
            alignment * second_state[2],
        )
        current = _state_pair_error(first_state, aligned_second)
        isfinite(current) || return Inf
        agreement = max(agreement, current)
    end
    return agreement
end

function _tag_principal_build_fallback(
    route::ComplexRoute,
    fallback_backend::Symbol,
    build_error::String;
    certificate=missing,
)
    metadata = merge(route.metadata, (;
        fallback_from=:principal_mst_logscaled,
        fallback_trigger=:principal_small_frequency_expansion_build_failure,
        fallback_backend,
        principal_small_frequency_expansion_build_error=build_error,
        principal_fallback_certificate=certificate,
    ))
    return ComplexRoute(
        route.branch,
        route.params,
        route.controls,
        route.p_solution,
        route.state_evaluator,
        metadata,
        route.plan,
        route.transmission,
        route.incidence,
        route.reflection,
        route.teukolsky_transmission,
        route.teukolsky_incidence,
        route.teukolsky_reflection,
    )
end

# When the principal-MST route fails: use a dual-order P-equation route if its
# orders agree, otherwise the rational route with small_frequency_expansion=false.
function _principal_small_frequency_expansion_build_fallback(
    s,
    l,
    m,
    a,
    omega,
    branch,
    requested,
    build_error,
)
    first = try
        with_logger(NullLogger()) do
            _complex_p_route(
                s, l, m, a, omega, branch;
                N=_PRINCIPAL_P_FALLBACK_ORDER,
                tol=requested.tolerance,
            )
        end
    catch
        nothing
    end
    confirmation = first === nothing ? nothing : try
        with_logger(NullLogger()) do
            _complex_p_route(
                s, l, m, a, omega, branch;
                N=_PRINCIPAL_P_FALLBACK_ORDER +
                    _PRINCIPAL_P_FALLBACK_ORDER_DELTA,
                tol=requested.tolerance,
            )
        end
    catch
        nothing
    end
    if first !== nothing && confirmation !== nothing
        _, _, amplitude_agreement = _route_pair_errors(first, confirmation)
        state_agreement = try
            _principal_p_state_agreement(first, confirmation)
        catch
            Inf
        end
        split_mismatch = max(
            _route_split_mismatch(first),
            _route_split_mismatch(confirmation),
        )
        certificate = (;
            amplitude_agreement,
            state_agreement,
            split_mismatch,
            first_order=_PRINCIPAL_P_FALLBACK_ORDER,
            confirmation_order=_PRINCIPAL_P_FALLBACK_ORDER +
                _PRINCIPAL_P_FALLBACK_ORDER_DELTA,
        )
        if isfinite(amplitude_agreement) &&
                amplitude_agreement <=
                    _PRINCIPAL_P_FALLBACK_AGREEMENT_MAX &&
                isfinite(state_agreement) &&
                state_agreement <= _PRINCIPAL_P_FALLBACK_AGREEMENT_MAX &&
                isfinite(split_mismatch) &&
                split_mismatch <= _PRINCIPAL_P_FALLBACK_SPLIT_MAX
            return _tag_principal_build_fallback(
                confirmation,
                :p_equation_dual_order,
                sprint(showerror, build_error);
                certificate,
            )
        end
    end
    rational = complex_route(
        s, l, m, a, omega, branch; small_frequency_expansion=false)
    return _tag_principal_build_fallback(
        rational,
        :direct_rational_small_frequency_expansion_off,
        sprint(showerror, build_error),
    )
end

@inline function _raw_triplet_errors(
    first::ComplexRoute,
    second::ComplexRoute,
)
    values = (
        first.teukolsky_transmission,
        first.teukolsky_incidence,
        first.teukolsky_reflection,
        second.teukolsky_transmission,
        second.teukolsky_incidence,
        second.teukolsky_reflection,
    )
    all(_finite_complex, values) || return (Inf, Inf, Inf)
    iszero(second.teukolsky_transmission) && return (Inf, Inf, Inf)
    alignment = first.teukolsky_transmission /
        second.teukolsky_transmission
    aligned_incidence = alignment * second.teukolsky_incidence
    aligned_reflection = alignment * second.teukolsky_reflection
    scale = max(
        abs(first.teukolsky_transmission),
        abs(first.teukolsky_incidence),
        abs(first.teukolsky_reflection),
        abs(aligned_incidence),
        abs(aligned_reflection),
        floatmin(Float64),
    )
    incidence_error = abs(
        first.teukolsky_incidence - aligned_incidence) / scale
    reflection_error = abs(
        first.teukolsky_reflection - aligned_reflection) / scale
    return incidence_error, reflection_error,
        max(incidence_error, reflection_error)
end

@inline function _raw_triplet_errors(first::NTuple{3,ComplexF64},
    second::NTuple{3,ComplexF64})
    values = (first..., second...)
    all(_finite_complex, values) || return (Inf, Inf, Inf)
    iszero(second[1]) && return (Inf, Inf, Inf)
    alignment = first[1] / second[1]
    aligned_incidence = alignment * second[2]
    aligned_reflection = alignment * second[3]
    scale = max(
        maximum(abs, first),
        abs(aligned_incidence),
        abs(aligned_reflection),
        floatmin(Float64),
    )
    incidence_error = abs(first[2] - aligned_incidence) / scale
    reflection_error = abs(first[3] - aligned_reflection) / scale
    return incidence_error, reflection_error,
        max(incidence_error, reflection_error)
end

@inline function _route_raw_triplet(route::ComplexRoute)
    return (
        route.teukolsky_transmission,
        route.teukolsky_incidence,
        route.teukolsky_reflection,
    )
end

@inline function _candidate_raw_triplet(candidate, factors, transmission)
    return (
        ComplexF64(factors.transmission * transmission),
        ComplexF64(factors.incidence * candidate.incidence),
        ComplexF64(factors.reflection * candidate.reflection),
    )
end

@inline function _raw_triplet_residual(values::NTuple{3,ComplexF64})
    return abs(values[2]) /
        max(maximum(abs, values), floatmin(Float64))
end

function _up_initial_match_metrics(rational::ComplexRoute, candidates)
    factors = rational.plan.amplitudes.conversion_factors
    transmission = rational.transmission
    selected_raw = _route_raw_triplet(rational)
    far_raw = _candidate_raw_triplet(
        candidates.far, factors, transmission)
    mid_raw = _candidate_raw_triplet(
        candidates.mid, factors, transmission)
    initial_raw = _candidate_raw_triplet(
        candidates.initial, factors, transmission)
    _, far_mid_reflection, far_mid_agreement =
        _raw_triplet_errors(far_raw, mid_raw)
    _, mid_initial_reflection, mid_initial_agreement =
        _raw_triplet_errors(mid_raw, initial_raw)
    _, far_initial_reflection, far_initial_agreement =
        _raw_triplet_errors(far_raw, initial_raw)
    _, _, selected_initial_separation =
        _raw_triplet_errors(selected_raw, initial_raw)
    reflection_agreement = max(
        far_mid_reflection,
        mid_initial_reflection,
        far_initial_reflection,
    )
    initial_residual = _raw_triplet_residual(initial_raw)
    split_mismatch = max(
        candidates.far.abel_error,
        candidates.mid.abel_error,
        candidates.initial.abel_error,
    )
    accepted = isfinite(initial_residual) &&
        initial_residual <= _UP_INITIAL_MATCH_RESIDUAL_MAX &&
        isfinite(far_mid_agreement) &&
        far_mid_agreement <= _UP_INITIAL_MATCH_FAR_MID_AGREEMENT_MAX &&
        isfinite(mid_initial_agreement) &&
        mid_initial_agreement <=
            _UP_INITIAL_MATCH_MID_INITIAL_AGREEMENT_MAX &&
        isfinite(far_initial_agreement) &&
        far_initial_agreement <=
            _UP_INITIAL_MATCH_FAR_INITIAL_AGREEMENT_MAX &&
        isfinite(reflection_agreement) &&
        reflection_agreement <=
            _UP_INITIAL_MATCH_REFLECTION_AGREEMENT_MAX &&
        isfinite(selected_initial_separation) &&
        selected_initial_separation >= _UP_INITIAL_MATCH_SEPARATION_MIN &&
        isfinite(split_mismatch) &&
        split_mismatch <= _UP_INITIAL_MATCH_SPLIT_MISMATCH_MAX
    return (;
        accepted,
        initial_residual,
        far_mid_agreement,
        mid_initial_agreement,
        far_initial_agreement,
        reflection_agreement,
        selected_initial_separation,
        split_mismatch,
    )
end

# A route rebuilt from a full match result of the rational solver (same
# coefficients and settings) chosen by a consensus retry; `extra` holds the
# retry's own diagnostics, appended to the metadata.
function _route_from_match(rational::ComplexRoute, candidate, policy, extra)
    Base.@nospecialize rational candidate extra
    factors = rational.plan.amplitudes.conversion_factors
    transmission = rational.transmission
    incidence = ComplexF64(candidate.incidence)
    reflection = ComplexF64(candidate.reflection)
    teukolsky = _candidate_raw_triplet(candidate, factors, transmission)
    evaluator = ComplexRationalEvaluator(
        rational.state_evaluator.params,
        rational.state_evaluator.coefficients,
        rational.state_evaluator.settings,
        candidate.match_x,
        candidate.match_state,
    )
    horizon_rho = candidate.horizon_seed.path.rho
    infinity_rho = candidate.infinity_seed.path.rho
    rho_match = max(abs(horizon_rho), abs(infinity_rho))
    propagation_score = max(
        candidate.horizon.max_score,
        candidate.infinity.max_score,
    )
    metadata = merge(rational.metadata, (;
        match_policy=policy,
        match_rstar=candidate.match_rstar,
        match_x=candidate.match_x,
        xsplit=candidate.match_x,
        xm_match=candidate.match_x,
        control_xm=candidate.match_x,
        control_rhom=rho_match,
        matching_condition=candidate.matching_condition,
        coefficient1_cancellation=candidate.coefficient1_cancellation,
        coefficient2_cancellation=candidate.coefficient2_cancellation,
        horizon_path_kind=candidate.horizon_seed.path.kind,
        horizon_direction=candidate.horizon_seed.path.direction,
        horizon_angle_offset=candidate.horizon_seed.path.angle_offset,
        horizon_rho,
        infinity_rho,
        horizon_coordinate_steps=candidate.horizon_seed.path.accepted_steps,
        infinity_coordinate_steps=candidate.infinity_seed.path.accepted_steps,
        horizon_endpoint_score=candidate.horizon_seed.score,
        infinity_endpoint_score=candidate.infinity_seed.score,
        horizon_patches=candidate.horizon.patches,
        infinity_patches=candidate.infinity.patches,
        patch_count=candidate.patch_count,
        propagation_score,
        abel_ratio=candidate.abel_ratio,
        abel_raw_error=candidate.abel_raw_error,
        abel_error=candidate.abel_error,
        split_mismatch=candidate.abel_error,
        endpoint_us=candidate.endpoint_us,
        propagation_us=candidate.propagation_us,
        matching_us=candidate.matching_us,
        endpoint_states=candidate.endpoint_states,
        extra...,
    ))
    contour = merge(rational.plan.contour, (;
        horizon_direction=ComplexF64(candidate.horizon_seed.path.direction),
        horizon_rho,
        infinity_rho,
        rho_match,
    ))
    matching = merge(rational.plan.matching, (;
        candidate_x=candidate.match_x,
        selected_x=candidate.match_x,
        split_mismatch=candidate.abel_error,
        condition=candidate.matching_condition,
        policy=policy,
    ))
    controls = merge(rational.plan.controls, (;
        xm=candidate.match_x,
        rhom=rho_match,
    ))
    amplitudes = merge(rational.plan.amplitudes, (;
        teukolsky=(
            transmission=teukolsky[1],
            incidence=teukolsky[2],
            reflection=teukolsky[3],
        ),
        raw_gsn_transmission=transmission,
        gsn=(; transmission, incidence, reflection),
    ))
    plan = ComplexRoutePlan(
        contour,
        matching,
        controls,
        amplitudes,
        candidate.patch_count,
    )
    return ComplexRoute(
        rational.branch,
        rational.params,
        rational.controls,
        rational.p_solution,
        evaluator,
        metadata,
        plan,
        transmission,
        incidence,
        reflection,
        teukolsky[1],
        teukolsky[2],
        teukolsky[3],
    )
end

function _up_initial_match_route(
    rational::ComplexRoute,
    candidates,
    metrics,
    elapsed_us,
)
    Base.@nospecialize rational candidates metrics
    return _route_from_match(rational, candidates.initial,
        :up_initial_match_consensus, (;
            up_initial_match_scale=candidates.scale,
            up_initial_match_far_rstar=candidates.far_rstar,
            up_initial_match_mid_rstar=candidates.mid_rstar,
            up_initial_match_initial_residual=metrics.initial_residual,
            up_initial_match_far_mid_agreement=metrics.far_mid_agreement,
            up_initial_match_mid_initial_agreement=
                metrics.mid_initial_agreement,
            up_initial_match_far_initial_agreement=
                metrics.far_initial_agreement,
            up_initial_match_reflection_agreement=
                metrics.reflection_agreement,
            up_initial_match_selected_separation=
                metrics.selected_initial_separation,
            up_initial_match_split_mismatch=metrics.split_mismatch,
            up_initial_match_consensus_us=elapsed_us,
            rational_match_policy=rational.metadata.match_policy,
            rational_match_x=rational.metadata.match_x,
            rational_matching_condition=rational.metadata.matching_condition,
            rational_coefficient1_cancellation=
                rational.metadata.coefficient1_cancellation,
            rational_coefficient2_cancellation=
                rational.metadata.coefficient2_cancellation,
            rational_patch_count=rational.plan.patch_count,
        ))
end

# s = 2 UP after a condition retry: return to the r* = 0 junction when it agrees
# with two shifted junctions and has a small incidence residual.
function _try_up_initial_match_consensus(rational::ComplexRoute)
    Base.@nospecialize rational
    _up_initial_match_pretrigger(rational) || return rational
    started = time_ns()
    candidates = try
        complex_up_initial_match_candidates(rational.state_evaluator)
    catch
        return rational
    end
    metrics = _up_initial_match_metrics(rational, candidates)
    metrics.accepted || return rational
    elapsed_us = (time_ns() - started) / 1.0e3
    return _up_initial_match_route(
        rational, candidates, metrics, elapsed_us)
end

@inline function _offpole_match_plateau_pretrigger(
    rational::ComplexRoute,
    selected::ComplexRoute,
)
    _endpoint_amplitude_match(rational) && return false
    rational.params.kappa <= _OFFPOLE_MATCH_PLATEAU_KAPPA_MAX ||
        return false
    if rational.branch == :UP
        _p_up_dual_pretrigger(rational) && return false
        return selected === rational
    end
    rational.branch == :IN && rational.params.s < 0 || return false
    selected === rational && return true
    return _metadata_value(selected.metadata, :backend, nothing) ==
        :direct_gsn_p_consensus
end

function _offpole_match_plateau_metrics(rational, candidates)
    agreement = candidates.agreement
    agreement_max = rational.branch == :UP ?
        _OFFPOLE_MATCH_PLATEAU_UP_AGREEMENT_MAX :
        _OFFPOLE_MATCH_PLATEAU_IN_AGREEMENT_MAX
    separation_ratio_min = rational.branch == :UP ?
        _OFFPOLE_MATCH_PLATEAU_UP_SEPARATION_RATIO_MIN :
        _OFFPOLE_MATCH_PLATEAU_IN_SEPARATION_RATIO_MIN
    selected = candidates.selected
    continuity = try
        _offpole_match_plateau_continuity(rational, selected)
    catch
        (xerror=Inf, dxerror=Inf, error=Inf)
    end
    member_agreement = 0.0
    if rational.branch == :IN &&
            continuity.error > _OFFPOLE_MATCH_PLATEAU_CONTINUITY_MAX
        for candidate in candidates.all_candidates
            abs(candidate.match_rstar) < abs(selected.match_rstar) || continue
            _, _, candidate_agreement = _amplitude_pair_errors(
                selected, candidate)
            candidate_agreement <=
                _OFFPOLE_MATCH_PLATEAU_MEMBER_AGREEMENT_MAX || continue
            candidate_continuity = try
                _offpole_match_plateau_continuity(rational, candidate)
            catch
                (xerror=Inf, dxerror=Inf, error=Inf)
            end
            candidate_continuity.error <=
                _OFFPOLE_MATCH_PLATEAU_CONTINUITY_MAX || continue
            selected = candidate
            continuity = candidate_continuity
            member_agreement = candidate_agreement
        end
    end
    _, _, separation = _amplitude_pair_errors(rational, selected)
    separation_ratio = separation /
        max(agreement, floatmin(Float64))
    basic_accepted = isfinite(agreement) && agreement <= agreement_max &&
        isfinite(separation_ratio) && separation_ratio >=
            separation_ratio_min &&
        isfinite(candidates.split_mismatch) &&
        candidates.split_mismatch <=
            _OFFPOLE_MATCH_PLATEAU_SPLIT_MISMATCH_MAX
    accepted = basic_accepted && isfinite(continuity.error) &&
        continuity.error <= _OFFPOLE_MATCH_PLATEAU_CONTINUITY_MAX
    return (;
        accepted,
        selected,
        agreement,
        member_agreement,
        separation,
        separation_ratio,
        split_mismatch=candidates.split_mismatch,
        continuity_x_error=continuity.xerror,
        continuity_dx_error=continuity.dxerror,
        continuity_error=continuity.error,
    )
end

function _offpole_match_plateau_continuity(rational, candidate)
    evaluator = ComplexRationalEvaluator(
        rational.state_evaluator.params,
        rational.state_evaluator.coefficients,
        rational.state_evaluator.settings,
        candidate.match_x,
        candidate.match_state,
    )
    rstar = candidate.match_rstar
    epsilon = min(1.0e-3,
        max(1.0e-7, 1.0e-6 * max(1.0, abs(rstar))))
    left2 = evaluator(rstar - 2epsilon)
    left = evaluator(rstar - epsilon)
    center = evaluator(rstar)
    right = evaluator(rstar + epsilon)
    right2 = evaluator(rstar + 2epsilon)
    values = (
        left2[1], left2[2], left[1], left[2], center[1], center[2],
        right[1], right[2], right2[1], right2[2],
    )
    all(_finite_complex, values) ||
        return (xerror=Inf, dxerror=Inf, error=Inf)
    left_to_center = left[1] + epsilon * left[2]
    right_to_center = right[1] - epsilon * right[2]
    xscale = max(
        abs(left_to_center), abs(right_to_center), abs(center[1]),
        floatmin(Float64),
    )
    left_dx_to_center = 2left[2] - left2[2]
    right_dx_to_center = 2right[2] - right2[2]
    dxscale = max(
        abs(left_dx_to_center), abs(right_dx_to_center), abs(center[2]),
        floatmin(Float64),
    )
    xerror = max(
        abs(left_to_center - right_to_center),
        abs(center[1] - (left_to_center + right_to_center) / 2),
    ) / xscale
    dxerror = max(
        abs(left_dx_to_center - right_dx_to_center),
        abs(center[2] - (left_dx_to_center + right_dx_to_center) / 2),
    ) / dxscale
    return (; xerror, dxerror, error=max(xerror, dxerror))
end

function _offpole_match_plateau_route(
    rational::ComplexRoute,
    selected::ComplexRoute,
    candidates,
    metrics,
    elapsed_us,
)
    Base.@nospecialize rational selected candidates metrics
    return _route_from_match(rational, metrics.selected,
        :offpole_match_plateau, (;
            offpole_match_plateau_triggered=true,
            offpole_match_plateau_accepted=true,
            offpole_match_plateau_candidate_count=candidates.candidate_count,
            offpole_match_plateau_reused_candidate_count=
                candidates.reused_candidate_count,
            offpole_match_plateau_agreement=metrics.agreement,
            offpole_match_plateau_member_agreement=metrics.member_agreement,
            offpole_match_plateau_separation=metrics.separation,
            offpole_match_plateau_separation_ratio=metrics.separation_ratio,
            offpole_match_plateau_split_mismatch=metrics.split_mismatch,
            offpole_match_plateau_continuity_x_error=
                metrics.continuity_x_error,
            offpole_match_plateau_continuity_dx_error=
                metrics.continuity_dx_error,
            offpole_match_plateau_continuity_error=metrics.continuity_error,
            offpole_match_plateau_first_rstar=candidates.first_rstar,
            offpole_match_plateau_second_rstar=candidates.second_rstar,
            offpole_match_plateau_us=elapsed_us,
            offpole_match_plateau_replaced_backend=
                _metadata_value(selected.metadata, :backend, missing),
            offpole_match_plateau_rational_match_policy=
                _metadata_value(rational.metadata, :match_policy, missing),
        ))
end

# Re-match at a ladder of junctions along r* and adopt the best-agreeing match
# when it is continuous across the junction and clearly separated from the
# selected route.
function _try_offpole_match_plateau(
    rational::ComplexRoute,
    selected::ComplexRoute,
    retry_cache=nothing,
)
    Base.@nospecialize rational selected retry_cache
    _offpole_match_plateau_pretrigger(rational, selected) || return selected
    started = time_ns()
    candidates = try
        complex_offpole_match_plateau_candidate(
            rational.state_evaluator, rational.branch, retry_cache)
    catch
        return selected
    end
    metrics = _offpole_match_plateau_metrics(rational, candidates)
    metrics.accepted || return selected
    if selected !== rational && metrics.separation_ratio >
            _OFFPOLE_MATCH_PLATEAU_P_OVERRIDE_RATIO_MAX
        return selected
    end
    elapsed_us = (time_ns() - started) / 1.0e3
    return _offpole_match_plateau_route(
        rational, selected, candidates, metrics, elapsed_us)
end

@inline function _up_horizon_in_order_pretrigger(route::ComplexRoute)
    metadata = route.metadata
    hasproperty(metadata, :match_policy) || return false
    hasproperty(metadata, :coefficient1_cancellation) || return false
    cancellation = Float64(metadata.coefficient1_cancellation)
    return route.params.s == -2 &&
        route.branch == :UP &&
        route.params.kappa <= _UP_HORIZON_IN_KAPPA_MAX &&
        metadata.match_policy == :up_condition_retry &&
        cancellation >= _UP_HORIZON_IN_CANCELLATION_MIN
end

function _up_horizon_in_order_metrics(
    rational::ComplexRoute,
    candidate,
    confirmation,
    final_candidate,
    orders,
    ;
    window::Symbol=:upper,
)
    _, _, first_agreement = _amplitude_pair_errors(candidate, confirmation)
    _, _, second_agreement =
        _amplitude_pair_errors(confirmation, final_candidate)
    incidence_error, reflection_error, separation =
        _amplitude_pair_errors(rational, final_candidate)
    split_mismatch = max(
        candidate.split_mismatch,
        confirmation.split_mismatch,
        final_candidate.split_mismatch,
    )
    accepted = isfinite(first_agreement) &&
        first_agreement <= _UP_HORIZON_IN_AGREEMENT_MAX &&
        isfinite(second_agreement) &&
        second_agreement <= _UP_HORIZON_IN_AGREEMENT_MAX &&
        isfinite(separation) &&
        separation >= _UP_HORIZON_IN_SEPARATION_MIN &&
        isfinite(split_mismatch) &&
        split_mismatch <= _UP_HORIZON_IN_SPLIT_MISMATCH_MAX
    return (;
        accepted,
        orders,
        window,
        first_agreement,
        second_agreement,
        incidence_error,
        reflection_error,
        separation,
        split_mismatch,
    )
end

function _up_horizon_in_order_route(
    rational::ComplexRoute,
    confirmation::ComplexRoute,
    candidate,
    metrics,
    elapsed_us,
)
    Base.@nospecialize rational confirmation candidate metrics
    factors = confirmation.plan.amplitudes.conversion_factors
    transmission = confirmation.transmission
    incidence = ComplexF64(candidate.incidence)
    reflection = ComplexF64(candidate.reflection)
    teukolsky_transmission = confirmation.teukolsky_transmission
    teukolsky_incidence = ComplexF64(factors.incidence * incidence)
    teukolsky_reflection = ComplexF64(factors.reflection * reflection)
    endpoint_states = (;
        horizon_in=candidate.horizon_in,
        horizon_out=candidate.horizon_out,
        infinity_in=missing,
        infinity_out=candidate.infinity_out,
    )
    propagation_score = max(
        candidate.propagated.max_score,
        Float64(_metadata_value(
            confirmation.metadata, :propagation_score, 0.0)),
    )
    metadata = merge(confirmation.metadata, (;
        match_policy=:up_horizon_in_order_consensus,
        matching_condition=candidate.matching_condition,
        coefficient1_cancellation=candidate.coefficient1_cancellation,
        coefficient2_cancellation=candidate.coefficient2_cancellation,
        horizon_path_kind=candidate.seed.path.kind,
        horizon_direction=candidate.seed.path.direction,
        horizon_angle_offset=candidate.seed.path.angle_offset,
        horizon_rho=candidate.seed.path.rho,
        horizon_coordinate_steps=candidate.seed.path.accepted_steps,
        horizon_endpoint_score=candidate.seed.score,
        horizon_patches=candidate.propagated.patches,
        patch_count=candidate.patch_count,
        propagation_score,
        abel_ratio=candidate.normalization.ratio,
        abel_raw_error=candidate.normalization.raw_error,
        abel_error=candidate.normalization.corrected_error,
        split_mismatch=candidate.normalization.corrected_error,
        endpoint_states,
        up_horizon_in_order_orders=metrics.orders,
        up_horizon_in_order_window=metrics.window,
        up_horizon_in_order_first_agreement=metrics.first_agreement,
        up_horizon_in_order_second_agreement=metrics.second_agreement,
        up_horizon_in_order_incidence_separation=metrics.incidence_error,
        up_horizon_in_order_reflection_separation=metrics.reflection_error,
        up_horizon_in_order_separation=metrics.separation,
        up_horizon_in_order_split_mismatch=metrics.split_mismatch,
        up_horizon_in_order_us=elapsed_us,
        rational_patch_count=rational.plan.patch_count,
        rational_matching_condition=_metadata_value(
            rational.metadata, :matching_condition, missing),
        rational_coefficient1_cancellation=_metadata_value(
            rational.metadata, :coefficient1_cancellation, missing),
        rational_coefficient2_cancellation=_metadata_value(
            rational.metadata, :coefficient2_cancellation, missing),
    ))
    contour = merge(confirmation.plan.contour, (;
        horizon_direction=ComplexF64(candidate.seed.path.direction),
        horizon_rho=candidate.seed.path.rho,
    ))
    matching = merge(confirmation.plan.matching, (;
        split_mismatch=candidate.normalization.corrected_error,
        condition=candidate.matching_condition,
        policy=:up_horizon_in_order_consensus,
    ))
    amplitudes = merge(confirmation.plan.amplitudes, (;
        teukolsky=(
            transmission=teukolsky_transmission,
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        ),
        raw_gsn_transmission=transmission,
        gsn=(; transmission, incidence, reflection),
    ))
    plan = ComplexRoutePlan(
        contour,
        matching,
        confirmation.plan.controls,
        amplitudes,
        candidate.patch_count,
    )
    return ComplexRoute(
        confirmation.branch,
        confirmation.params,
        confirmation.controls,
        confirmation.p_solution,
        confirmation.state_evaluator,
        metadata,
        plan,
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
end

# s = -2 UP near extremal spin after a condition retry: recompute the amplitudes
# with a freshly propagated horizon-in solution at three rational orders
# (from the route order upwards, else from below it), accepted when
# consecutive orders agree.
function _try_up_horizon_in_order_consensus(
    rational::ComplexRoute,
    s,
    l,
    m,
    a,
    omega,
    branch,
)
    Base.@nospecialize rational
    _up_horizon_in_order_pretrigger(rational) || return rational
    started = time_ns()
    order = rational.state_evaluator.settings.order
    tolerance = rational.plan.controls.tolerance
    order1 = order + _UP_HORIZON_IN_ORDER_DELTA
    order2 = order1 + _UP_HORIZON_IN_ORDER_DELTA
    candidate = try
        complex_up_horizon_in_candidate(
            rational.state_evaluator, rational.metadata)
    catch
        return rational
    end
    params = rational.state_evaluator.params
    confirmation_route = try
        _complex_rational_route(
            s, l, m, a, omega, branch;
            lambda=params.lambda,
            nu=params.nu,
            N=order1,
            tol=tolerance,
        )
    catch
        return rational
    end
    confirmation = try
        complex_up_horizon_in_candidate(
            confirmation_route.state_evaluator, confirmation_route.metadata)
    catch
        return rational
    end
    final_route = try
        _complex_rational_route(
            s, l, m, a, omega, branch;
            lambda=params.lambda,
            nu=params.nu,
            N=order2,
            tol=tolerance,
        )
    catch
        return rational
    end
    final_candidate = try
        complex_up_horizon_in_candidate(
            final_route.state_evaluator, final_route.metadata)
    catch
        return rational
    end
    metrics = _up_horizon_in_order_metrics(
        rational,
        candidate,
        confirmation,
        final_candidate,
        (order, order1, order2),
    )
    if metrics.accepted
        elapsed_us = (time_ns() - started) / 1.0e3
        return _up_horizon_in_order_route(
            rational, final_route, final_candidate, metrics, elapsed_us)
    end

    lower0 = order - _UP_HORIZON_IN_LOWER_BACKOFF
    lower0 >= _UP_HORIZON_IN_ORDER_DELTA || return rational
    lower_orders = (
        lower0,
        lower0 + _UP_HORIZON_IN_ORDER_DELTA,
        lower0 + 2 * _UP_HORIZON_IN_ORDER_DELTA,
    )
    lower_routes = try
        map(lower_orders) do lower_order
            _complex_rational_route(
                s, l, m, a, omega, branch;
                lambda=params.lambda,
                nu=params.nu,
                N=lower_order,
                tol=tolerance,
            )
        end
    catch
        return rational
    end
    lower_candidates = try
        map(lower_routes) do lower_route
            complex_up_horizon_in_candidate(
                lower_route.state_evaluator, lower_route.metadata)
        end
    catch
        return rational
    end
    lower_metrics = _up_horizon_in_order_metrics(
        rational,
        lower_candidates...,
        lower_orders;
        window=:lower,
    )
    lower_metrics.accepted || return rational
    elapsed_us = (time_ns() - started) / 1.0e3
    return _up_horizon_in_order_route(
        rational,
        lower_routes[end],
        lower_candidates[end],
        lower_metrics,
        elapsed_us,
    )
end

function _p_up_dual_metrics(
    rational::ComplexRoute,
    candidate::ComplexRoute,
    confirmation::ComplexRoute,
    confirmation_order::Int,
)
    incidence_error, reflection_error, separation =
        _route_pair_errors(rational, confirmation)
    _, _, agreement = _route_pair_errors(candidate, confirmation)
    candidate_split_mismatch = _route_split_mismatch(candidate)
    confirmation_split_mismatch = _route_split_mismatch(confirmation)
    split_mismatch = max(
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
    direct_residual = _raw_incidence_residual(rational)
    candidate_residual = _raw_incidence_residual(confirmation)
    residual_ratio = candidate_residual /
        max(direct_residual, floatmin(Float64))
    accepted = isfinite(agreement) &&
        agreement <= _P_CONSENSUS_UP_AGREEMENT_MAX &&
        isfinite(separation) &&
        separation >= _P_CONSENSUS_UP_SEPARATION_MIN &&
        isfinite(split_mismatch) &&
        split_mismatch <= _P_CONSENSUS_UP_SPLIT_MISMATCH_MAX
    return (;
        accepted,
        acceptance_kind=:up_dual_order,
        direct_residual,
        candidate_residual,
        residual_ratio,
        incidence_error,
        reflection_error,
        split_mismatch,
        confirmation_order,
        confirmation_agreement=agreement,
        candidate_separation=separation,
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
end

function _p_contour_metrics(
    rational::ComplexRoute,
    candidate::ComplexRoute,
    confirmation::ComplexRoute,
    confirmation_order::Int,
)
    incidence_error, reflection_error, separation =
        _route_pair_errors(rational, confirmation)
    _, _, agreement = _route_pair_errors(candidate, confirmation)
    candidate_split_mismatch = _route_split_mismatch(candidate)
    confirmation_split_mismatch = _route_split_mismatch(confirmation)
    split_mismatch = max(
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
    direct_residual = _raw_incidence_residual(rational)
    candidate_residual = _raw_incidence_residual(confirmation)
    residual_ratio = candidate_residual /
        max(direct_residual, floatmin(Float64))
    strict_accepted = isfinite(agreement) &&
        agreement <= _P_CONTOUR_ORDER_AGREEMENT_MAX &&
        isfinite(separation) &&
        separation >= _P_CONTOUR_SEPARATION_MIN &&
        isfinite(split_mismatch) &&
        split_mismatch <= _P_CONTOUR_SPLIT_MISMATCH_MAX
    usable_accepted = !strict_accepted &&
        isfinite(agreement) &&
        agreement <= _P_CONTOUR_USABLE_ORDER_AGREEMENT_MAX &&
        isfinite(separation) &&
        separation >= _P_CONTOUR_USABLE_SEPARATION_MIN &&
        isfinite(split_mismatch) &&
        split_mismatch <= _P_CONTOUR_USABLE_SPLIT_MISMATCH_MAX
    accepted = strict_accepted || usable_accepted
    acceptance_kind = strict_accepted ?
        :contour_dual_order : :contour_dual_order_usable
    select_candidate = usable_accepted &&
        candidate_split_mismatch <= confirmation_split_mismatch
    return (;
        accepted,
        acceptance_kind,
        select_candidate,
        direct_residual,
        candidate_residual,
        residual_ratio,
        incidence_error,
        reflection_error,
        split_mismatch,
        confirmation_order,
        confirmation_agreement=agreement,
        candidate_separation=separation,
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
end

function _p_nia_metrics(
    rational::ComplexRoute,
    candidate::ComplexRoute,
    confirmation::ComplexRoute,
    confirmation_order::Int,
)
    incidence_error, reflection_error, separation =
        _unit_pair_errors(rational, confirmation)
    _, _, agreement = _unit_pair_errors(candidate, confirmation)
    state_agreement =
        _principal_p_unit_state_agreement(candidate, confirmation)
    candidate_split_mismatch = _route_split_mismatch(candidate)
    confirmation_split_mismatch = _route_split_mismatch(confirmation)
    split_mismatch = max(
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
    direct_residual = _raw_incidence_residual(rational)
    candidate_residual = _raw_incidence_residual(confirmation)
    residual_ratio = candidate_residual /
        max(direct_residual, floatmin(Float64))
    nonrotating = _p_nia_nonrotating(rational)
    state_accepted = nonrotating ||
        (isfinite(state_agreement) &&
            state_agreement <= _P_NIA_STATE_AGREEMENT_MAX)
    separation_accepted = nonrotating ||
        (isfinite(separation) && separation >= _P_NIA_SEPARATION_MIN)
    accepted = isfinite(agreement) &&
        agreement <= _P_NIA_ORDER_AGREEMENT_MAX &&
        state_accepted &&
        separation_accepted
    acceptance_kind = nonrotating ?
        :nia_dual_order_nonrotating : :nia_dual_order
    return (;
        accepted,
        acceptance_kind,
        select_candidate=false,
        direct_residual,
        candidate_residual,
        residual_ratio,
        incidence_error,
        reflection_error,
        split_mismatch,
        state_agreement,
        confirmation_order,
        confirmation_agreement=agreement,
        candidate_separation=separation,
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
end

function _p_nia_build_metrics(
    fallback::ComplexRoute,
    candidate::ComplexRoute,
    confirmation::ComplexRoute,
    confirmation_order::Int,
)
    incidence_error, reflection_error, separation =
        _unit_pair_errors(fallback, confirmation)
    _, _, agreement = _unit_pair_errors(candidate, confirmation)
    state_agreement =
        _principal_p_unit_state_agreement(candidate, confirmation)
    candidate_split_mismatch = _route_split_mismatch(candidate)
    confirmation_split_mismatch = _route_split_mismatch(confirmation)
    split_mismatch = max(
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
    direct_residual = _raw_incidence_residual(fallback)
    candidate_residual = _raw_incidence_residual(confirmation)
    residual_ratio = candidate_residual /
        max(direct_residual, floatmin(Float64))
    accepted = isfinite(agreement) &&
        agreement <= _P_NIA_BUILD_ORDER_AGREEMENT_MAX &&
        isfinite(state_agreement) &&
        state_agreement <= _P_NIA_BUILD_STATE_AGREEMENT_MAX &&
        isfinite(separation) &&
        separation >= _P_NIA_SEPARATION_MIN
    return (;
        accepted,
        acceptance_kind=:nia_build_dual_order,
        select_candidate=false,
        direct_residual,
        candidate_residual,
        residual_ratio,
        incidence_error,
        reflection_error,
        split_mismatch,
        state_agreement,
        confirmation_order,
        confirmation_agreement=agreement,
        candidate_separation=separation,
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
end

function _p_pole_dual_metrics(
    rational::ComplexRoute,
    candidate::ComplexRoute,
    confirmation::ComplexRoute,
    confirmation_order::Int,
)
    incidence_error, reflection_error, separation =
        _raw_triplet_errors(rational, confirmation)
    _, _, agreement = _raw_triplet_errors(candidate, confirmation)
    candidate_split_mismatch = _route_split_mismatch(candidate)
    confirmation_split_mismatch = _route_split_mismatch(confirmation)
    split_mismatch = max(
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
    direct_residual = _raw_incidence_residual(rational)
    candidate_residual = _raw_incidence_residual(confirmation)
    residual_ratio = candidate_residual /
        max(direct_residual, floatmin(Float64))
    accepted = isfinite(candidate_residual) &&
        candidate_residual <= _P_CONSENSUS_RESIDUAL_MAX &&
        isfinite(residual_ratio) &&
        residual_ratio <= _P_CONSENSUS_POLE_DUAL_RESIDUAL_RATIO_MAX &&
        isfinite(agreement) &&
        agreement <= _P_CONSENSUS_POLE_DUAL_AGREEMENT_MAX &&
        isfinite(separation) &&
        separation >= _P_CONSENSUS_POLE_DUAL_SEPARATION_MIN &&
        isfinite(reflection_error) &&
        reflection_error <=
            _P_CONSENSUS_POLE_DUAL_REFLECTION_ERROR_MAX &&
        isfinite(split_mismatch) &&
        split_mismatch <=
            _P_CONSENSUS_POLE_DUAL_SPLIT_MISMATCH_MAX
    return (;
        accepted,
        acceptance_kind=:pole_dual_order,
        direct_residual,
        candidate_residual,
        residual_ratio,
        incidence_error,
        reflection_error,
        split_mismatch,
        confirmation_order,
        confirmation_agreement=agreement,
        candidate_separation=separation,
        candidate_split_mismatch,
        confirmation_split_mismatch,
    )
end

function _p_consensus_metrics(
    rational::ComplexRoute,
    candidate::ComplexRoute,
)
    direct_residual = _raw_incidence_residual(rational)
    candidate_residual = _raw_incidence_residual(candidate)
    transmission = candidate.teukolsky_transmission
    _finite_complex(transmission) && !iszero(transmission) ||
        return (accepted=false, acceptance_kind=:rejected,
            direct_residual, candidate_residual, residual_ratio=Inf,
            incidence_error=Inf, reflection_error=Inf,
            split_mismatch=Inf)
    alignment = rational.teukolsky_transmission / transmission
    aligned_incidence = alignment * candidate.teukolsky_incidence
    aligned_reflection = alignment * candidate.teukolsky_reflection
    scale = max(
        abs(rational.teukolsky_transmission),
        abs(rational.teukolsky_reflection),
        abs(aligned_reflection),
        floatmin(Float64),
    )
    reflection_error = abs(
        rational.teukolsky_reflection - aligned_reflection) / scale
    incidence_error = abs(
        rational.teukolsky_incidence - aligned_incidence) /
        max(scale, abs(rational.teukolsky_incidence),
            abs(aligned_incidence))
    split_value = candidate.plan.matching.split_mismatch
    split_mismatch = split_value === missing ? Inf : Float64(split_value)
    residual_ratio = candidate_residual /
        max(direct_residual, floatmin(Float64))
    pole_accepted = isfinite(candidate_residual) &&
        candidate_residual <= _P_CONSENSUS_RESIDUAL_MAX &&
        candidate_residual <=
            _P_CONSENSUS_RESIDUAL_RATIO_MAX * direct_residual &&
        isfinite(reflection_error) &&
        reflection_error <= _P_CONSENSUS_REFLECTION_ERROR_MAX &&
        isfinite(split_mismatch) &&
        split_mismatch <= _P_CONSENSUS_SPLIT_MISMATCH_MAX
    scattering_accepted = rational.branch == :IN &&
        isfinite(candidate_residual) &&
        candidate_residual >= _P_CONSENSUS_SCATTERING_RESIDUAL_MIN &&
        isfinite(residual_ratio) &&
        _P_CONSENSUS_SCATTERING_RESIDUAL_RATIO_MIN <= residual_ratio <=
            _P_CONSENSUS_SCATTERING_RESIDUAL_RATIO_MAX &&
        isfinite(reflection_error) &&
        reflection_error <=
            _P_CONSENSUS_SCATTERING_REFLECTION_ERROR_MAX &&
        isfinite(split_mismatch) &&
        split_mismatch <= _P_CONSENSUS_SPLIT_MISMATCH_MAX
    accepted = pole_accepted || scattering_accepted
    acceptance_kind = pole_accepted ? :pole_residual :
        (scattering_accepted ? :scattering_cross_branch : :rejected)
    return (; accepted, acceptance_kind, direct_residual,
        candidate_residual, residual_ratio, incidence_error,
        reflection_error, split_mismatch)
end

# The adopted candidate's own error: its distance to the two-ray value plus
# that value's error, relative to the candidate's amplitude (both in
# Teukolsky normalisation, the candidate scaled to the same transmission).
function _p_consensus_errors(
    rational::ComplexRoute,
    candidate::ComplexRoute,
)
    alignment = rational.teukolsky_transmission /
        candidate.teukolsky_transmission
    aligned = alignment .* (candidate.teukolsky_incidence,
        candidate.teukolsky_reflection)
    two_ray = (rational.teukolsky_incidence, rational.teukolsky_reflection)
    return Float64.((abs.(two_ray .- aligned) .+
        amplitude_errors(rational) .* abs.(two_ray)) ./ abs.(aligned))
end

function _p_consensus_route(
    rational::ComplexRoute,
    candidate::ComplexRoute,
    metrics,
)
    Base.@nospecialize rational candidate metrics
    rational_metadata = rational.metadata
    errors = _p_consensus_errors(rational, candidate)
    metadata = merge(candidate.metadata, (;
        backend=:direct_gsn_p_consensus,
        incidence_error=errors[1],
        reflection_error=errors[2],
        consensus_source=:p_equation,
        consensus_acceptance_kind=metrics.acceptance_kind,
        consensus_direct_residual=metrics.direct_residual,
        consensus_candidate_residual=metrics.candidate_residual,
        consensus_residual_ratio=metrics.residual_ratio,
        consensus_incidence_error=metrics.incidence_error,
        consensus_reflection_error=metrics.reflection_error,
        consensus_split_mismatch=metrics.split_mismatch,
        consensus_confirmation_order=hasproperty(
            metrics, :confirmation_order) ?
            metrics.confirmation_order : missing,
        consensus_confirmation_agreement=hasproperty(
            metrics, :confirmation_agreement) ?
            metrics.confirmation_agreement : missing,
        consensus_state_agreement=hasproperty(
            metrics, :state_agreement) ?
            metrics.state_agreement : missing,
        consensus_candidate_separation=hasproperty(
            metrics, :candidate_separation) ?
            metrics.candidate_separation : missing,
        consensus_candidate_split_mismatch=hasproperty(
            metrics, :candidate_split_mismatch) ?
            metrics.candidate_split_mismatch : missing,
        consensus_confirmation_split_mismatch=hasproperty(
            metrics, :confirmation_split_mismatch) ?
            metrics.confirmation_split_mismatch : missing,
        rational_patch_count=rational.plan.patch_count,
        rational_matching_condition=hasproperty(
            rational_metadata, :matching_condition) ?
            rational_metadata.matching_condition : missing,
        rational_coefficient1_cancellation=hasproperty(
            rational_metadata, :coefficient1_cancellation) ?
            rational_metadata.coefficient1_cancellation : missing,
        rational_coefficient2_cancellation=hasproperty(
            rational_metadata, :coefficient2_cancellation) ?
            rational_metadata.coefficient2_cancellation : missing,
    ))
    return ComplexRoute(
        candidate.branch,
        candidate.params,
        candidate.controls,
        candidate.p_solution,
        candidate.state_evaluator,
        metadata,
        candidate.plan,
        candidate.transmission,
        candidate.incidence,
        candidate.reflection,
        candidate.teukolsky_transmission,
        candidate.teukolsky_incidence,
        candidate.teukolsky_reflection,
    )
end

# Rebuild the route with the P-equation backend when the rational match looks
# unreliable (amplitude cancellation near a pole, large matching condition or
# split mismatch, the negative imaginary axis); the P-equation route is adopted
# only when its own dual-order or residual checks pass.
function _try_p_consensus(
    rational::ComplexRoute,
    s,
    l,
    m,
    a,
    omega,
    branch,
)
    Base.@nospecialize rational
    regular_trigger = _p_consensus_pretrigger(rational)
    pole_dual_trigger = _p_pole_dual_pretrigger(rational)
    up_dual_trigger = _p_up_dual_pretrigger(rational)
    contour_trigger = _p_contour_pretrigger(rational)
    nia_trigger = _p_nia_pretrigger(rational)
    controls = rational.plan.controls
    if nia_trigger
        candidate = try
            with_logger(NullLogger()) do
                _complex_p_route(
                    s, l, m, a, omega, branch;
                    N=_P_NIA_ORDER,
                    tol=controls.tolerance,
                )
            end
        catch
            nothing
        end
        confirmation_order = _P_NIA_ORDER + _P_NIA_ORDER_DELTA
        confirmation = try
            with_logger(NullLogger()) do
                _complex_p_route(
                    s, l, m, a, omega, branch;
                    N=confirmation_order,
                    tol=controls.tolerance,
                )
            end
        catch
            nothing
        end
        if candidate !== nothing && confirmation !== nothing
            metrics = _p_nia_metrics(
                rational, candidate, confirmation, confirmation_order)
            if metrics.accepted
                return metrics.acceptance_kind ==
                    :nia_dual_order_nonrotating ?
                    _nia_p_overlay(
                        rational,
                        confirmation,
                        s,
                        metrics.confirmation_agreement,
                        :same_spin,
                    ) :
                    _p_consensus_route(rational, confirmation, metrics)
            end
        end
    end
    (regular_trigger || pole_dual_trigger || up_dual_trigger ||
        contour_trigger) || return rational
    candidate = try
        with_logger(NullLogger()) do
            _complex_p_route(
                s, l, m, a, omega, branch;
                N=controls.N,
                tol=controls.tolerance,
            )
        end
    catch
        return rational
    end
    if regular_trigger
        metrics = _p_consensus_metrics(rational, candidate)
        metrics.accepted &&
            return _p_consensus_route(rational, candidate, metrics)
    end
    if contour_trigger
        confirmation_order = controls.N + _P_CONSENSUS_UP_ORDER_DELTA
        confirmation = try
            with_logger(NullLogger()) do
                _complex_p_route(
                    s, l, m, a, omega, branch;
                    N=confirmation_order,
                    tol=controls.tolerance,
                )
            end
        catch
            nothing
        end
        if confirmation !== nothing
            metrics = _p_contour_metrics(
                rational, candidate, confirmation, confirmation_order)
            if metrics.accepted
                selected = metrics.select_candidate ?
                    candidate : confirmation
                return _p_consensus_route(rational, selected, metrics)
            end
        end
    end
    if pole_dual_trigger &&
        _route_split_mismatch(candidate) <=
            _P_CONSENSUS_POLE_DUAL_SPLIT_MISMATCH_MAX
        confirmation_order = controls.N + _P_CONSENSUS_UP_ORDER_DELTA
        confirmation = try
            with_logger(NullLogger()) do
                _complex_p_route(
                    s, l, m, a, omega, branch;
                    N=confirmation_order,
                    tol=controls.tolerance,
                )
            end
        catch
            nothing
        end
        if confirmation !== nothing
            metrics = _p_pole_dual_metrics(
                rational, candidate, confirmation, confirmation_order)
            metrics.accepted &&
                return _p_consensus_route(rational, confirmation, metrics)
        end
    end
    up_dual_trigger || return rational
    _route_split_mismatch(candidate) <=
        _P_CONSENSUS_UP_SPLIT_MISMATCH_MAX || return rational
    confirmation_order = controls.N + _P_CONSENSUS_UP_ORDER_DELTA
    confirmation = try
        with_logger(NullLogger()) do
            _complex_p_route(
                s, l, m, a, omega, branch;
                N=confirmation_order,
                tol=controls.tolerance,
            )
        end
    catch
        return rational
    end
    metrics = _p_up_dual_metrics(
        rational, candidate, confirmation, confirmation_order)
    return metrics.accepted ?
        _p_consensus_route(rational, confirmation, metrics) : rational
end

# After a rational build failure near the negative imaginary axis, adopt a
# dual-order P-equation route for the amplitudes when the two orders agree.
function _try_p_nia_build_consensus(
    fallback::ComplexRoute,
    s,
    l,
    m,
    a,
    omega,
    branch,
)
    _p_nia_build_pretrigger(s, a, omega, branch) || return fallback
    tolerance = fallback.plan.controls.tolerance
    candidate = try
        with_logger(NullLogger()) do
            _complex_p_route(
                s, l, m, a, omega, branch;
                N=_P_NIA_ORDER,
                tol=tolerance,
            )
        end
    catch
        return fallback
    end
    confirmation_order = _P_NIA_ORDER + _P_NIA_ORDER_DELTA
    confirmation = try
        with_logger(NullLogger()) do
            _complex_p_route(
                s, l, m, a, omega, branch;
                N=confirmation_order,
                tol=tolerance,
            )
        end
    catch
        return fallback
    end
    if _p_nia_nonrotating(fallback)
        _, _, agreement = _unit_pair_errors(candidate, confirmation)
        if isfinite(agreement) &&
                agreement <= _P_NIA_ORDER_AGREEMENT_MAX
            return _nia_p_overlay(
                fallback,
                confirmation,
                s,
                agreement,
                :same_spin,
            )
        end
    end
    metrics = _p_nia_build_metrics(
        fallback, candidate, confirmation, confirmation_order)
    return metrics.accepted ?
        _p_consensus_route(fallback, confirmation, metrics) : fallback
end

function _nia_p_overlay(
    route::ComplexRoute,
    source_route::ComplexRoute,
    source_spin::Int,
    agreement,
    source::Symbol,
)
    factors = route.plan.amplitudes.conversion_factors
    values = (
        route.transmission,
        route.teukolsky_transmission,
        factors.transmission,
        factors.incidence,
        factors.reflection,
        source_route.transmission,
        source_route.incidence,
        source_route.reflection,
    )
    all(_finite_complex, values) || return route
    iszero(route.transmission) && return route
    iszero(factors.transmission) && return route
    iszero(source_route.transmission) && return route

    incidence = ComplexF64(
        route.transmission *
        source_route.incidence / source_route.transmission)
    reflection = ComplexF64(
        route.transmission *
        source_route.reflection / source_route.transmission)
    teukolsky_scale = route.teukolsky_transmission /
        (factors.transmission * route.transmission)
    teukolsky_incidence =
        ComplexF64(factors.incidence * incidence * teukolsky_scale)
    teukolsky_reflection =
        ComplexF64(factors.reflection * reflection * teukolsky_scale)
    all(_finite_complex, (
        incidence,
        reflection,
        teukolsky_incidence,
        teukolsky_reflection,
    )) || return route

    metadata = merge(route.metadata, (;
        amplitude_backend=source == :spin_partner ?
            :nonrotating_spin_partner_p : :nonrotating_nia_p,
        nia_p_amplitude_source=source,
        nia_p_source_spin=source_spin,
        nia_p_orders=(
            _P_NIA_ORDER,
            _P_NIA_ORDER + _P_NIA_ORDER_DELTA,
        ),
        nia_p_order_agreement=agreement,
        nia_p_original_incidence=route.incidence,
        nia_p_original_reflection=route.reflection,
    ))
    amplitudes = merge(route.plan.amplitudes, (;
        teukolsky=merge(route.plan.amplitudes.teukolsky, (;
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        )),
        gsn=merge(route.plan.amplitudes.gsn, (;
            incidence,
            reflection,
        )),
    ))
    plan = ComplexRoutePlan(
        route.plan.contour,
        route.plan.matching,
        route.plan.controls,
        amplitudes,
        route.plan.patch_count,
    )
    return ComplexRoute(
        route.branch,
        route.params,
        route.controls,
        route.p_solution,
        route.state_evaluator,
        metadata,
        plan,
        route.transmission,
        incidence,
        reflection,
        route.teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
end

function _try_nia_partner(
    route::ComplexRoute,
    s,
    l,
    m,
    a,
    omega,
    branch,
)
    _p_nia_partner_pretrigger(s, a, omega) || return route
    # Function barrier: the triggered work is compiled only when it runs.
    return Base.inferencebarrier(_nia_partner_retry)(
        route, s, l, m, a, omega, branch)
end

# For s < 0 on the negative imaginary axis of a Schwarzschild hole, replace the
# unit-transmission amplitudes by those of the spin partner -s from two
# P-equation orders that agree.
function _nia_partner_retry(route, s, l, m, a, omega, branch)
    partner_spin = -s
    tolerance = route.plan.controls.tolerance
    candidate = try
        with_logger(NullLogger()) do
            _complex_p_route(
                partner_spin, l, m, a, omega, branch;
                N=_P_NIA_ORDER,
                tol=tolerance,
            )
        end
    catch
        return route
    end
    confirmation = try
        with_logger(NullLogger()) do
            _complex_p_route(
                partner_spin, l, m, a, omega, branch;
                N=_P_NIA_ORDER + _P_NIA_ORDER_DELTA,
                tol=tolerance,
            )
        end
    catch
        return route
    end
    _, _, agreement = _unit_pair_errors(candidate, confirmation)
    isfinite(agreement) &&
        agreement <= _P_NIA_PARTNER_ORDER_AGREEMENT_MAX || return route
    return _nia_p_overlay(
        route, confirmation, partner_spin, agreement, :spin_partner)
end

@inline function _teuk_unit_pair(route::ComplexRoute)
    values = (
        route.teukolsky_transmission,
        route.teukolsky_incidence,
        route.teukolsky_reflection,
    )
    all(_finite_complex, values) || return nothing
    iszero(route.teukolsky_transmission) && return nothing
    return (
        incidence=ComplexF64(
            route.teukolsky_incidence / route.teukolsky_transmission),
        reflection=ComplexF64(
            route.teukolsky_reflection / route.teukolsky_transmission),
    )
end

function _nia_side_pair(
    s,
    l,
    m,
    a,
    omega,
    q,
    target_factors,
    in_pair,
)
    up_route = with_logger(NullLogger()) do
        _complex_rational_route(
            s, l, m, a, omega, :UP)
    end
    up_pair = _teuk_unit_pair(up_route)
    up_pair === nothing && error("nonfinite right-side UP amplitudes")
    iszero(in_pair.incidence) &&
        error("zero right-side Teukolsky incidence")
    chi = ComplexF64(
        1 - im * q * in_pair.reflection / in_pair.incidence)
    teukolsky = (
        incidence=ComplexF64(chi * up_pair.incidence),
        reflection=ComplexF64(
            chi * up_pair.reflection + im * q / in_pair.incidence),
    )
    values = (
        target_factors.transmission,
        target_factors.incidence,
        target_factors.reflection,
        chi,
        teukolsky.incidence,
        teukolsky.reflection,
    )
    all(_finite_complex, values) ||
        error("nonfinite NIA lateral connection")
    iszero(target_factors.transmission) &&
        error("zero target transmission conversion factor")
    iszero(target_factors.incidence) &&
        error("zero target incidence conversion factor")
    iszero(target_factors.reflection) &&
        error("zero target reflection conversion factor")
    gsn = (
        incidence=ComplexF64(
            teukolsky.incidence *
            target_factors.transmission / target_factors.incidence),
        reflection=ComplexF64(
            teukolsky.reflection *
            target_factors.transmission / target_factors.reflection),
    )
    all(_finite_complex, (gsn.incidence, gsn.reflection)) ||
        error("nonfinite GSN NIA lateral amplitudes")
    return (
        gsn=gsn,
        teukolsky=teukolsky,
        chi=chi,
        patch_count=up_route.plan.patch_count,
    )
end

function _nia_in_pair(
    s,
    l,
    m,
    a,
    omega,
    tolerance,
)
    in_route = with_logger(NullLogger()) do
        _complex_p_route(
            s, l, m, a, omega, :IN;
            N=_P_NIA_ORDER,
            tol=tolerance,
        )
    end
    in_confirmation = with_logger(NullLogger()) do
        _complex_p_route(
            s, l, m, a, omega, :IN;
            N=_P_NIA_ORDER + _P_NIA_ORDER_DELTA,
            tol=tolerance,
        )
    end
    _, _, in_agreement = _unit_pair_errors(
        in_route, in_confirmation)
    isfinite(in_agreement) &&
        in_agreement <= _P_NIA_PARTNER_ORDER_AGREEMENT_MAX ||
        error("right-side IN dual-order certificate rejected")
    in_pair = _teuk_unit_pair(in_confirmation)
    in_pair === nothing && error("nonfinite right-side IN amplitudes")
    return (
        pair=in_pair,
        order_agreement=in_agreement,
        patch_count=in_confirmation.plan.patch_count,
        backend=:p_dual_order,
    )
end

function _nia_in_rational(
    s,
    l,
    m,
    a,
    omega,
)
    route = with_logger(NullLogger()) do
        _complex_rational_route(
            s, l, m, a, omega, :IN;
            N=_P_NIA_ORDER,
        )
    end
    confirmation = with_logger(NullLogger()) do
        _complex_rational_route(
            s, l, m, a, omega, :IN;
            N=_P_NIA_ORDER + _P_NIA_ORDER_DELTA,
        )
    end
    _, _, agreement = _unit_pair_errors(route, confirmation)
    isfinite(agreement) &&
        agreement <= _NIA_SIDE_AGREEMENT_MAX ||
        error("right-side rational IN dual-order certificate rejected")
    pair = _teuk_unit_pair(confirmation)
    pair === nothing && error("nonfinite rational right-side IN amplitudes")
    return (
        pair=pair,
        order_agreement=agreement,
        patch_count=confirmation.plan.patch_count,
        backend=:rational_dual_order,
    )
end

function _nia_side_overlay(
    route::ComplexRoute,
    mapped,
    confirmation,
    strength,
    offsets,
    agreement,
)
    factors = route.plan.amplitudes.conversion_factors
    iszero(route.transmission) && return route
    iszero(factors.transmission) && return route
    incidence = ComplexF64(
        route.transmission * mapped.gsn.incidence)
    reflection = ComplexF64(
        route.transmission * mapped.gsn.reflection)
    teukolsky_scale = route.teukolsky_transmission /
        (factors.transmission * route.transmission)
    teukolsky_incidence = ComplexF64(
        factors.incidence * incidence * teukolsky_scale)
    teukolsky_reflection = ComplexF64(
        factors.reflection * reflection * teukolsky_scale)
    all(_finite_complex, (
        incidence,
        reflection,
        teukolsky_incidence,
        teukolsky_reflection,
    )) || return route

    metadata = merge(route.metadata, (;
        amplitude_backend=:nia_lateral_connection,
        nia_lateral_offsets=offsets,
        nia_lateral_source_agreement=agreement,
        nia_lateral_q=ComplexF64(strength.branch_strength),
        nia_lateral_q_estimated_relerr=strength.estimated_relerr,
        nia_lateral_q_adjacent_error=strength.adjacent_error,
        nia_lateral_chi=mapped.chi,
        nia_lateral_confirmation_chi=confirmation.chi,
        nia_lateral_source_patch_counts=(
            mapped.patch_count,
            confirmation.patch_count,
        ),
        nia_lateral_original_incidence=route.incidence,
        nia_lateral_original_reflection=route.reflection,
    ))
    amplitudes = merge(route.plan.amplitudes, (;
        teukolsky=merge(route.plan.amplitudes.teukolsky, (;
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        )),
        gsn=merge(route.plan.amplitudes.gsn, (;
            incidence,
            reflection,
        )),
    ))
    plan = ComplexRoutePlan(
        route.plan.contour,
        route.plan.matching,
        route.plan.controls,
        amplitudes,
        route.plan.patch_count,
    )
    return ComplexRoute(
        route.branch,
        route.params,
        route.controls,
        route.p_solution,
        route.state_evaluator,
        metadata,
        plan,
        route.transmission,
        incidence,
        reflection,
        route.teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
end

function _try_nia_side(
    route::ComplexRoute,
    s,
    l,
    m,
    a,
    omega,
)
    _nia_side_pretrigger(route) || return route
    # Function barrier: the triggered work is compiled only when it runs.
    return Base.inferencebarrier(_nia_side_retry)(route, s, l, m, a, omega)
end

# UP amplitudes just left of the negative imaginary axis (Re omega < 0) of a
# rotating hole from the NIA lateral connection: the MST branch strength q on
# the axis combined with the IN amplitudes at a small positive real offset and
# the UP amplitudes at two such offsets, accepted when the two agree.
function _nia_side_retry(route, s, l, m, a, omega)
    omegac = ComplexF64(omega)
    axis_omega = ComplexF64(0, imag(omegac))
    strength = try
        axis_params = gsn_parameters(
            s, l, m, a, axis_omega)
        nia_strength_fast(
            s, l, m, a, axis_omega, axis_params.lambda)
    catch
        return route
    end
    strength.certificate_accepted || return route
    q = ComplexF64(strength.branch_strength)
    epsilon = abs(real(omegac))
    offsets = (epsilon, _NIA_SIDE_OFFSET_FACTOR * epsilon)
    factors = route.plan.amplitudes.conversion_factors
    tolerance = route.plan.controls.tolerance
    in_source = try
        _nia_in_pair(
            s, l, m, a,
            ComplexF64(offsets[1], imag(omegac)),
            tolerance,
        )
    catch
        try
            _nia_in_rational(
                s, l, m, a,
                ComplexF64(offsets[1], imag(omegac)),
            )
        catch
            return route
        end
    end
    first = try
        _nia_side_pair(
            s, l, m, a,
            ComplexF64(offsets[1], imag(omegac)),
            q,
            factors,
            in_source.pair,
        )
    catch
        return route
    end
    second = try
        _nia_side_pair(
            s, l, m, a,
            ComplexF64(offsets[2], imag(omegac)),
            q,
            factors,
            in_source.pair,
        )
    catch
        return route
    end
    _, _, agreement = _amplitude_pair_errors(
        first.gsn, second.gsn)
    isfinite(agreement) &&
        agreement <= _NIA_SIDE_AGREEMENT_MAX || return route
    overlaid = _nia_side_overlay(
        route, first, second, strength, offsets, agreement)
    overlaid === route && return route
    metadata = merge(overlaid.metadata, (;
        nia_lateral_in_order_agreement=in_source.order_agreement,
        nia_lateral_in_patch_count=in_source.patch_count,
        nia_lateral_in_backend=in_source.backend,
    ))
    return ComplexRoute(
        overlaid.branch,
        overlaid.params,
        overlaid.controls,
        overlaid.p_solution,
        overlaid.state_evaluator,
        metadata,
        overlaid.plan,
        overlaid.transmission,
        overlaid.incidence,
        overlaid.reflection,
        overlaid.teukolsky_transmission,
        overlaid.teukolsky_incidence,
        overlaid.teukolsky_reflection,
    )
end

@inline function _pole_reciprocal_pretrigger(route::ComplexRoute)
    route.branch == :IN || return false
    route.params.s > 0 || return false
    route.params.a * route.params.m > 0 || return false
    route.params.kappa <= _POLE_RECIPROCAL_KAPPA_MAX || return false
    route.state_evaluator isa ComplexRationalEvaluator || return false
    _metadata_value(route.metadata, :match_policy, missing) == :rstar_zero ||
        return false
    residual = _raw_incidence_residual(route)
    return _POLE_RECIPROCAL_IN_RESIDUAL_MIN <= residual <=
        _POLE_RECIPROCAL_IN_RESIDUAL_MAX
end

@inline function _relative_complex_error(first, second)
    scale = max(abs(first), abs(second), floatmin(Float64))
    return abs(first - second) / scale
end

function _pole_reciprocal_state_error(
    in_route::ComplexRoute,
    up_route::ComplexRoute,
    scale,
)
    in_x, in_dx, _ = in_route.state_evaluator(0.0)
    up_x, up_dx, _ = up_route.state_evaluator(0.0)
    scaled_up = (scale * up_x, scale * up_dx)
    in_values = (in_x, in_dx)
    state_scale = max(
        maximum(abs, in_values),
        maximum(abs, scaled_up),
        floatmin(Float64),
    )
    return maximum(abs.(in_values .- scaled_up)) / state_scale
end

@inline function _scaled_rational_state(state, scale)
    return typeof(state)(scale * state.X, scale * state.dXdx)
end

function _pole_reciprocal_route(
    original::ComplexRoute,
    up_route::ComplexRoute,
    confirmation::ComplexRoute,
    metrics,
    elapsed_us,
)
    Base.@nospecialize original up_route confirmation metrics
    reflection = ComplexF64(metrics.reflection)
    transmission = ComplexF64(1)
    incidence = ComplexF64(0)
    factors = _amplitude_factors(original.params, :IN)
    teukolsky_transmission = ComplexF64(factors.transmission)
    teukolsky_incidence = ComplexF64(0)
    teukolsky_reflection = ComplexF64(factors.reflection * reflection)

    source_evaluator = up_route.state_evaluator
    match_state = _scaled_rational_state(
        source_evaluator.match_state, reflection)
    evaluator = ComplexRationalEvaluator(
        source_evaluator.params,
        source_evaluator.coefficients,
        source_evaluator.settings,
        source_evaluator.match_x,
        match_state,
    )
    source_endpoints = up_route.metadata.endpoint_states
    infinity_out = _scaled_rational_state(
        source_endpoints.infinity_out, reflection)
    endpoint_states = (;
        horizon_in=missing,
        horizon_out=missing,
        infinity_in=missing,
        infinity_out,
    )

    metadata = merge(up_route.metadata, (;
        match_policy=:exact_pole_reciprocal_up_consensus,
        endpoint_states,
        pole_normalization=:exact_incidence_zero,
        pole_incidence_before=original.incidence,
        pole_teukolsky_incidence_before=original.teukolsky_incidence,
        pole_incidence_residual_before=metrics.in_residual,
        pole_reciprocal_source=:up_dual_order,
        pole_reciprocal_orders=metrics.orders,
        pole_reciprocal_up_residual=metrics.up_residual,
        pole_reciprocal_confirmation_residual=
            metrics.confirmation_residual,
        pole_reciprocal_order_agreement=metrics.order_agreement,
        pole_reciprocal_state_error=metrics.state_error,
        pole_reciprocal_separation=metrics.separation,
        pole_reciprocal_split_mismatch=metrics.split_mismatch,
        pole_reciprocal_us=elapsed_us,
        pole_reciprocal_original_match_policy=
            _metadata_value(original.metadata, :match_policy, missing),
        pole_reciprocal_original_matching_condition=
            _metadata_value(original.metadata, :matching_condition, missing),
        pole_reciprocal_original_reflection=original.reflection,
        pole_reciprocal_source_endpoint_states=source_endpoints,
        pole_reciprocal_confirmation_match_policy=
            _metadata_value(confirmation.metadata, :match_policy, missing),
    ))
    matching = merge(up_route.plan.matching, (;
        policy=:exact_pole_reciprocal_up_consensus,
    ))
    amplitudes = (;
        teukolsky=(
            transmission=teukolsky_transmission,
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        ),
        conversion_factors=factors,
        raw_gsn_transmission=transmission,
        gsn=(; transmission, incidence, reflection),
    )
    plan = ComplexRoutePlan(
        up_route.plan.contour,
        matching,
        up_route.plan.controls,
        amplitudes,
        up_route.plan.patch_count,
    )
    return ComplexRoute(
        :IN,
        original.params,
        original.controls,
        nothing,
        evaluator,
        metadata,
        plan,
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )
end

# With pole normalization at an IN pole (s > 0, a*m > 0, near-extremal spin),
# take the reflection amplitude as the reciprocal of the UP reflection, accepted
# when two UP orders agree and the UP route's state matches the IN state.
function _try_exact_pole_reciprocal(
    route::ComplexRoute,
    s,
    l,
    m,
    a,
    omega,
)
    Base.@nospecialize route
    _pole_reciprocal_pretrigger(route) || return route
    started = time_ns()
    evaluator = route.state_evaluator
    order = evaluator.settings.order
    tolerance = route.plan.controls.tolerance
    up_route = try
        _complex_rational_route(
            s, l, m, a, omega, :UP;
            lambda=route.params.lambda,
            nu=evaluator.params.nu,
            N=order,
            tol=tolerance,
        )
    catch
        return route
    end
    confirmation_order = order + _POLE_RECIPROCAL_ORDER_DELTA
    confirmation = try
        _complex_rational_route(
            s, l, m, a, omega, :UP;
            lambda=route.params.lambda,
            nu=evaluator.params.nu,
            N=confirmation_order,
            tol=tolerance,
        )
    catch
        return route
    end
    iszero(up_route.reflection) && return route
    reflection = ComplexF64(inv(up_route.reflection))
    all(_finite_complex, (
        reflection,
        up_route.reflection,
        confirmation.reflection,
    )) || return route
    in_residual = _raw_incidence_residual(route)
    up_residual = _raw_incidence_residual(up_route)
    confirmation_residual = _raw_incidence_residual(confirmation)
    order_agreement = _relative_complex_error(
        up_route.reflection, confirmation.reflection)
    state_error = _pole_reciprocal_state_error(
        route, up_route, reflection)
    separation = _relative_complex_error(route.reflection, reflection)
    split_mismatch = max(
        _route_split_mismatch(up_route),
        _route_split_mismatch(confirmation),
    )
    accepted = up_residual <= _POLE_RECIPROCAL_UP_RESIDUAL_MAX &&
        confirmation_residual <= _POLE_RECIPROCAL_UP_RESIDUAL_MAX &&
        order_agreement <= _POLE_RECIPROCAL_ORDER_AGREEMENT_MAX &&
        state_error <= _POLE_RECIPROCAL_STATE_ERROR_MAX &&
        separation >= _POLE_RECIPROCAL_SEPARATION_MIN &&
        split_mismatch <= _POLE_RECIPROCAL_SPLIT_MISMATCH_MAX
    accepted || return route
    metrics = (;
        reflection,
        orders=(order, confirmation_order),
        in_residual,
        up_residual,
        confirmation_residual,
        order_agreement,
        state_error,
        separation,
        split_mismatch,
    )
    elapsed_us = (time_ns() - started) / 1.0e3
    return _pole_reciprocal_route(
        route, up_route, confirmation, metrics, elapsed_us)
end

function _apply_exact_pole_normalization(route::ComplexRoute)
    Base.@nospecialize route
    incidence = zero(route.incidence)
    teukolsky_incidence = zero(route.teukolsky_incidence)
    metadata = merge(route.metadata, (;
        pole_normalization=:exact_incidence_zero,
        pole_incidence_before=route.incidence,
        pole_teukolsky_incidence_before=route.teukolsky_incidence,
        pole_incidence_residual_before=_raw_incidence_residual(route),
    ))
    amplitudes = merge(route.plan.amplitudes, (;
        teukolsky=merge(route.plan.amplitudes.teukolsky, (;
            incidence=teukolsky_incidence,
        )),
        gsn=merge(route.plan.amplitudes.gsn, (;
            incidence,
        )),
    ))
    plan = ComplexRoutePlan(
        route.plan.contour,
        route.plan.matching,
        route.plan.controls,
        amplitudes,
        route.plan.patch_count,
    )
    return ComplexRoute(
        route.branch,
        route.params,
        route.controls,
        route.p_solution,
        route.state_evaluator,
        metadata,
        plan,
        route.transmission,
        incidence,
        route.reflection,
        route.teukolsky_transmission,
        teukolsky_incidence,
        route.teukolsky_reflection,
    )
end

@inline function _apply_exact_pole_normalization(
    route::ComplexRoute,
    requested::Bool,
)
    return requested ? _apply_exact_pole_normalization(route) : route
end

# Consensus retries for the two-ray rational route, in priority order: the UP
# horizon-in order retry, then a pole dual-order P-equation retry; otherwise
# the UP initial-match retry is preferred over any other accepted P-equation
# retry. A retry is built only when no earlier one has already been returned.
function _select_consensus_retry(rational, s, l, m, a, omega, branch)
    horizon_retry = _try_up_horizon_in_order_consensus(
        rational, s, l, m, a, omega, branch)
    horizon_retry !== rational && return horizon_retry
    p_retry = _try_p_consensus(
        rational, s, l, m, a, omega, branch)
    p_retry !== rational && _metadata_value(
        p_retry.metadata, :consensus_acceptance_kind, nothing) ==
        :pole_dual_order && return p_retry
    initial_match_retry = _try_up_initial_match_consensus(rational)
    return initial_match_retry !== rational ? initial_match_retry : p_retry
end

function _complex_route_impl(
    s::Int,
    l::Int,
    m::Int,
    a,
    omega,
    branch;
    backend::Symbol=:direct_rational,
    controls::Union{Nothing,ComplexControls}=nothing,
    lambda=nothing,
    nu=nothing,
    xm=nothing,
    rhom=nothing,
    N=nothing,
    tol=nothing,
    small_frequency_expansion=nothing,
    large_frequency_expansion=nothing,
    pole_normalization::Bool=false,
)
    if backend == :mst
        selected = _complex_mst_route(
            s, l, m, a, omega, branch;
            controls,
            lambda,
            nu,
            xm,
            rhom,
            N,
            tol,
            small_frequency_expansion,
            large_frequency_expansion,
        )
        return _apply_exact_pole_normalization(
            selected, pole_normalization)
    end
    # Controls that fix the route itself (matching point, order, eigenvalue,
    # expansions) are not overridden by retries or fallbacks; a tolerance alone
    # sets accuracy targets, which the retries and fallbacks carry along.
    pinned_controls = controls !== nothing || lambda !== nothing ||
        nu !== nothing || xm !== nothing || rhom !== nothing ||
        N !== nothing || small_frequency_expansion !== nothing ||
        large_frequency_expansion !== nothing
    if backend == :direct_rational
        _, omegac = _validate_mode(s, l, m, a, omega)
        principal_requested = controls === nothing ?
            small_frequency_expansion !== false : controls.small_frequency_expansion !== false
        if _principal_small_frequency_expansion_axis(omegac) && principal_requested
            requested = controls === nothing ? complex_controls(
                xm=xm,
                rhom=rhom,
                N=N,
                tol=tol,
                small_frequency_expansion=small_frequency_expansion,
                large_frequency_expansion=large_frequency_expansion,
            ) : controls
            controls !== nothing && any(value !== nothing for value in
                (xm, rhom, N, tol, small_frequency_expansion, large_frequency_expansion)) &&
                throw(ArgumentError(
                    "pass either controls or individual complex-route " *
                    "control keywords, not both."))
            requested.rhom === nothing || throw(ArgumentError(
                "rhom is not used by the principal-MST complex route."))
            requested.large_frequency_expansion === true && throw(ArgumentError(
                "complex large_frequency_expansion=true is incompatible with principal-MST small-frequency-expansion."))
            selected = try
                _complex_principal_small_frequency_expansion_route(
                    s, l, m, a, omegac, branch, requested; lambda, nu)
            catch error
                pinned_controls && rethrow()
                _principal_small_frequency_expansion_build_fallback(
                    s, l, m, a, omegac, branch, requested, error)
            end
            return _apply_exact_pole_normalization(
                selected, pole_normalization)
        end
    end
    if backend == :direct_rational
        rational_build = try
            _complex_rational_route(
                s, l, m, a, omega, branch;
                controls,
                lambda,
                nu,
                xm,
                rhom,
                N,
                tol,
                small_frequency_expansion,
                large_frequency_expansion,
                return_retry_cache=true,
            )
        catch error
            fallback_trigger = _rational_build_fallback_trigger(error)
            principal_fallback = !pinned_controls &&
                _principal_small_frequency_expansion_axis(ComplexF64(omega))
            if pinned_controls ||
                    (fallback_trigger === nothing && !principal_fallback)
                rethrow()
            end
            fallback_trigger === nothing &&
                (fallback_trigger = :principal_mst_evaluator_fallback)
            fallback = _complex_p_route(
                s, l, m, a, omega, branch; tol)
            fallback = _tag_p_build_fallback(
                fallback, fallback_trigger, sprint(showerror, error))
            fallback = _principal_small_frequency_expansion_overlay(
                fallback, fallback.controls)
            fallback = _try_p_nia_build_consensus(
                fallback, s, l, m, a, omega, branch)
            fallback = _try_nia_partner(
                fallback, s, l, m, a, omega, branch)
            fallback = _try_nia_side(
                fallback, s, l, m, a, omega)
            return _apply_exact_pole_normalization(
                fallback, pole_normalization)
        end
        rational = _principal_small_frequency_expansion_overlay(
            rational_build.route, rational_build.route.controls)
        _metadata_value(
            rational.metadata, :amplitude_backend, :complex_rational) ==
            :principal_mst_logscaled &&
            return _apply_exact_pole_normalization(
                rational, pole_normalization)
        pinned_controls && return _apply_exact_pole_normalization(
            rational, pole_normalization)
        selected = _select_consensus_retry(
            rational, s, l, m, a, omega, branch)
        pole_normalization || (selected =
            _try_offpole_match_plateau(
                rational, selected, rational_build.retry_cache))
        selected = _try_nia_partner(
            selected, s, l, m, a, omega, branch)
        selected = _try_nia_side(
            selected, s, l, m, a, omega)
        if pole_normalization
            reciprocal = _try_exact_pole_reciprocal(
                selected, s, l, m, a, omega)
            reciprocal !== selected && return reciprocal
        end
        return _apply_exact_pole_normalization(
            selected, pole_normalization)
    elseif backend == :p_equation
        lambda === nothing || throw(ArgumentError(
            "lambda is not accepted by the legacy backend=:p_equation builder."))
        nu === nothing || throw(ArgumentError(
            "nu is not accepted by the legacy backend=:p_equation builder."))
        selected = _complex_p_route(
            s, l, m, a, omega, branch;
            controls,
            xm,
            rhom,
            N,
            tol,
            small_frequency_expansion,
            large_frequency_expansion,
        )
        return _apply_exact_pole_normalization(
            selected, pole_normalization)
    end
    throw(ArgumentError(
        "complex-frequency backend must be :direct_rational, :mst, or " *
        ":p_equation."))
end

# Builder hooks used by `complex_route`: a special-case builder tried
# first (returns `nothing` to decline), the default builder, and a step that
# post-processes the built route. NIAAmplitudes registers the first and
# the last.
const _SPECIAL_ROUTE_BUILDER = Ref{Any}(
    (args...; kwargs...) -> nothing)
const _DEFAULT_ROUTE_BUILDER = Ref{Any}(
    _complex_route_impl)
const _POST_ROUTE_BUILDER = Ref{Any}(
    (route, args...; kwargs...) -> route)

function _register_special_route_builder!(builder)
    _SPECIAL_ROUTE_BUILDER[] = builder
    return nothing
end

function _register_post_route_builder!(builder)
    _POST_ROUTE_BUILDER[] = builder
    return nothing
end

"""
    complex_route(s, l, m, a, omega, branch; backend=:direct_rational, kwargs...)

Build the IN or UP GSN solution at complex `omega` (nonzero imaginary part).
`backend = :direct_rational` (default) integrates along two physical r* rays
and matches at a real junction, with consensus retries and fallbacks where the
match is ill-conditioned; `:mst` uses the complete MST representation and
`:p_equation` the P-equation builder. Passing route-fixing controls (`xm`, `N`,
`lambda`, `nu`, `small_frequency_expansion`, ...) disables the retries and fallbacks.
`pole_normalization=true` sets the incidence amplitude exactly to zero (for
use at a mode frequency).
"""
function complex_route(
    s::Int,
    l::Int,
    m::Int,
    a,
    omega,
    branch;
    backend::Symbol=:direct_rational,
    controls::Union{Nothing,ComplexControls}=nothing,
    lambda=nothing,
    nu=nothing,
    xm=nothing,
    rhom=nothing,
    N=nothing,
    tol=nothing,
    small_frequency_expansion=nothing,
    large_frequency_expansion=nothing,
    pole_normalization::Bool=false,
)
    special_builder =
        Base.inferencebarrier(_SPECIAL_ROUTE_BUILDER[])
    special = special_builder(
        s, l, m, a, omega, branch;
        backend,
        controls,
        lambda,
        nu,
        xm,
        rhom,
        N,
        tol,
        small_frequency_expansion,
        large_frequency_expansion,
        pole_normalization,
    )
    special === nothing || return special
    default_builder =
        Base.inferencebarrier(_DEFAULT_ROUTE_BUILDER[])
    route = default_builder(
        s, l, m, a, omega, branch;
        backend,
        controls,
        lambda,
        nu,
        xm,
        rhom,
        N,
        tol,
        small_frequency_expansion,
        large_frequency_expansion,
        pole_normalization,
    )
    post_builder =
        Base.inferencebarrier(_POST_ROUTE_BUILDER[])
    return post_builder(
        route, s, l, m, a, omega, branch;
        backend,
        controls,
        lambda,
        nu,
        xm,
        rhom,
        N,
        tol,
        small_frequency_expansion,
        large_frequency_expansion,
        pole_normalization,
    )
end

complex_route_plan_of(route::ComplexRoute) = route.plan
complex_route_patch_count(route::ComplexRoute) = route.plan.patch_count

complex_route_plan_of(
    route::ConjugatedRoute{R},
) where {R<:ComplexRoute} = complex_route_plan_of(route.route)

complex_route_patch_count(
    route::ConjugatedRoute{R},
) where {R<:ComplexRoute} = complex_route_patch_count(route.route)

"""
    complex_amplitudes(route)

GSN asymptotic amplitudes of a complex-frequency route in unit-transmission
normalization, with the names of the incidence and reflection amplitudes
(`Binc`/`Bref` for IN, `Cinc`/`Cref` for UP).
"""
function complex_amplitudes(route::ComplexRoute)
    names = route.branch == :IN ? (:Binc, :Bref) : (:Cinc, :Cref)
    return (
        transmission=route.transmission,
        incidence=route.incidence,
        reflection=route.reflection,
        incidence_name=names[1],
        reflection_name=names[2],
    )
end

function complex_amplitudes(
    route::ConjugatedRoute{R},
) where {R<:ComplexRoute}
    amplitudes = complex_amplitudes(route.route)
    return merge(amplitudes, (
        transmission=conj(amplitudes.transmission),
        incidence=conj(amplitudes.incidence),
        reflection=conj(amplitudes.reflection),
    ))
end

"""
    complex_nia_jump(route)

Analytic (MST) jump of the amplitudes across the negative imaginary frequency
axis for an IN route, from its unit-transmission Teukolsky amplitudes.
"""
function complex_nia_jump(route::ComplexRoute)
    route.branch == :IN ||
        throw(ArgumentError(
            "the analytic NIA jump requires an IN route."))
    transmission = ComplexF64(route.teukolsky_transmission)
    _finite_complex(transmission) && !iszero(transmission) ||
        throw(ArgumentError(
            "the analytic NIA jump requires finite Teukolsky transmission."))
    return mst_nia_jump(
        route.params.s,
        route.params.l,
        route.params.m,
        route.params.a,
        route.params.omega,
        route.params.lambda,
        route.teukolsky_incidence / transmission,
        route.teukolsky_reflection / transmission,
    )
end

"""
    complex_endpoint_states(route)

Endpoint states of a complex-frequency route: the horizon and infinity
solutions used in the matching (`missing` where the backend did not record one).
"""
function complex_endpoint_states(route::ComplexRoute)
    m = route.metadata
    if _metadata_value(m, :backend, nothing) == :direct_gsn_two_ray_rational
        states = _metadata_value(m, :endpoint_states, missing)
        states === missing || return route.branch == :IN ? (
            horizon_in=states.horizon_in,
            infinity_in=states.infinity_in,
            infinity_out=states.infinity_out,
        ) : (
            horizon_in=states.horizon_in,
            horizon_out=states.horizon_out,
            infinity_out=states.infinity_out,
        )
    end
    if route.branch == :IN
        return (
            horizon_in=_metadata_value(m, :P_zero_in, missing),
            infinity_in=_metadata_value(m, :P_inf_in, missing),
            infinity_out=_metadata_value(m, :P_inf_out, missing),
        )
    end
    return (
        horizon_in=_metadata_value(m, :P_zero_in, missing),
        horizon_out=_metadata_value(m, :P_zero_out, missing),
        infinity_out=_metadata_value(m, :P_inf_out, missing),
    )
end

complex_route_state_r(route::ComplexRoute, r::Real) =
    _state_at_r(route.state_evaluator, r)

function complex_route_state_r(
    route::ConjugatedRoute{R},
    r::Real,
) where {R<:ComplexRoute}
    state = complex_route_state_r(route.route, r)
    return (
        X=conj(state.X),
        dXdrstar=conj(state.dXdrstar),
        error=state.error,
    )
end

function complex_route_state_rstar(route::ComplexRoute, rstar::Real)
    X, dXdrstar, error = route.state_evaluator(rstar)
    return (X=X, dXdrstar=dXdrstar, error=error)
end

function complex_route_state_rstar(
    route::ConjugatedRoute{R},
    rstar::Real,
) where {R<:ComplexRoute}
    state = complex_route_state_rstar(route.route, rstar)
    return (
        X=conj(state.X),
        dXdrstar=conj(state.dXdrstar),
        error=state.error,
    )
end

complex_route_evaluate_r(route::ComplexRoute, r::Real) =
    complex_route_state_r(route, r).X
complex_route_evaluate_rstar(route::ComplexRoute, rstar::Real) =
    complex_route_state_rstar(route, rstar).X
complex_route_evaluate_r(
    route::ConjugatedRoute{R},
    r::Real,
) where {R<:ComplexRoute} = complex_route_state_r(route, r).X
complex_route_evaluate_rstar(
    route::ConjugatedRoute{R},
    rstar::Real,
) where {R<:ComplexRoute} =
    complex_route_state_rstar(route, rstar).X

"""
    evaluate_complex_route_on_real_axis(route, coordinates; coordinate=:rstar, state=false)

Evaluate the GSN function of a complex-frequency route at real `coordinates`
(`coordinate = :rstar` or `:r`). Returns the values of X, or with `state=true`
named tuples `(X, dXdrstar, error)`.
"""
function evaluate_complex_route_on_real_axis(
    route::ComplexRoute,
    coordinates;
    coordinate::Symbol=:rstar,
    state::Bool=false,
)
    if coordinate == :rstar && coordinates isa AbstractVector{<:Real}
        values = route.state_evaluator(coordinates)
        return state ? [
            (
                X=values[1][index],
                dXdrstar=values[2][index],
                error=values[3][index],
            )
            for index in eachindex(values[1])
        ] : values[1]
    end
    evaluator = if coordinate == :rstar
        state ? (value -> complex_route_state_rstar(route, value)) :
            (value -> complex_route_evaluate_rstar(route, value))
    elseif coordinate == :r
        state ? (value -> complex_route_state_r(route, value)) :
            (value -> complex_route_evaluate_r(route, value))
    else
        throw(ArgumentError("coordinate must be :rstar or :r."))
    end
    return map(evaluator, coordinates)
end

complex_gsn_solution_rstar(route::ComplexRoute) = route.state_evaluator

struct ComplexTeukolskyPair{R,C,T}
    route::R
    converter::C
    inv_scale::T
end

@inline function _teukolsky_state(evaluator::ComplexTeukolskyPair, r::Real)
    route = evaluator.route
    state = complex_route_state_r(route, r)
    p = route.params
    x = (Float64(r) - p.rplus) / (Float64(r) - p.rminus)
    dxdrstar = dx_drstar(p, x)
    iszero(dxdrstar) && throw(DomainError(
        r, "ISEM complex Teukolsky evaluation requires a finite exterior radius."))
    dXdx = state.dXdrstar / dxdrstar
    R, Rp = evaluator.converter(r, state.X, dXdx)
    # Near the horizon the conversion can cancel strongly (for s = -2 the
    # ingoing R ~ Delta^2 comes from terms of size one), multiplying the state's
    # error; there the horizon series with the route's Teukolsky amplitudes is
    # used when its error (the amplitudes' errors) is the smaller.
    value_part = evaluator.converter(r, state.X, zero(dXdx))
    derivative_part = evaluator.converter(r, zero(state.X), dXdx)
    cancellation = max(
        (abs(value_part[1]) + abs(derivative_part[1])) / abs(R),
        (abs(value_part[2]) + abs(derivative_part[2])) / abs(Rp))
    conversion_error = state.error * cancellation + eps(Float64) * cancellation
    if Float64(r) - p.rplus < 2 * p.kappa && conversion_error > eps(Float64)
        series = _teukolsky_horizon_state(evaluator, r)
        series !== nothing && series.error < conversion_error && return series
    end
    return (
        R=evaluator.inv_scale * R,
        Rp=evaluator.inv_scale * Rp,
        error=isfinite(conversion_error) ? conversion_error : Inf,
    )
end

function _teukolsky_horizon_state(evaluator::ComplexTeukolskyPair, r::Real)
    route = evaluator.route
    p = route.params
    Rin, dRin, Rout, dRout, tail = teukolsky_horizon_pair(
        p.s, p.m, p.a, p.omega, p.lambda, r)
    isfinite(tail) || return nothing
    if route.branch == :IN
        A = route.teukolsky_transmission
        R = evaluator.inv_scale * A * Rin
        Rp = evaluator.inv_scale * A * dRin
        return (R=R, Rp=Rp, error=tail)
    end
    # The horizon coefficients of the solution the state evaluator represents
    # (they differ from the route's amplitudes when an analytic overlay
    # replaced those), in the route's Teukolsky normalization.
    A = route.teukolsky_reflection
    B = route.teukolsky_incidence
    incidence_error, reflection_error = amplitude_errors(route)
    evaluator_state = route.state_evaluator
    if evaluator_state isa ComplexRationalEvaluator &&
            evaluator_state.amplitudes[] !== nothing &&
            !iszero(route.reflection) && !iszero(route.incidence)
        a, b, da, db = evaluator_state.amplitudes[]
        A *= a / route.reflection
        B *= b / route.incidence
        reflection_error = da / abs(a)
        incidence_error = db / abs(b)
    end
    R = evaluator.inv_scale * (A * Rin + B * Rout)
    Rp = evaluator.inv_scale * (A * dRin + B * dRout)
    error = max(
        (abs(A * Rin) * reflection_error + abs(B * Rout) * incidence_error) /
            abs(R / evaluator.inv_scale),
        (abs(A * dRin) * reflection_error + abs(B * dRout) * incidence_error) /
            abs(Rp / evaluator.inv_scale)) + tail
    return (R=R, Rp=Rp, error=isfinite(error) ? error : Inf)
end

function (evaluator::ComplexTeukolskyPair)(r::Real)
    state = _teukolsky_state(evaluator, r)
    return state.R, state.Rp
end


struct ComplexTeukolskyTuple{P,E,F}
    params::P
    pair::E
    d2R::F
end

function (evaluator::ComplexTeukolskyTuple)(r::Real)
    state = _teukolsky_state(evaluator.pair, r)
    p = evaluator.params
    Rpp = evaluator.d2R(r, state.R, state.Rp,
        p.s, p.a, p.omega, p.m, p.lambda)
    return state.R, state.Rp, Rpp, state.error
end

function complex_gsn_solution_rstar(
    route::ConjugatedRoute{R},
) where {R<:ComplexRoute}
    source = complex_gsn_solution_rstar(route.route)
    return rstar -> begin
        X, dXdrstar, error = source(rstar)
        return conj(X), conj(dXdrstar), error
    end
end

function _public_match_rstar(route::ComplexRoute, backend)
    backend == :direct_gsn_two_ray_rational &&
        return _metadata_value(route.metadata, :match_rstar, missing)
    selected_x = route.plan.matching.selected_x
    selected_x === missing && return missing
    p = route.params
    candidate_r = p.rplus - 2 * p.kappa * selected_x
    iszero(imag(candidate_r)) && real(candidate_r) > p.rplus ||
        return missing
    return rstar_from_r(p.a, real(candidate_r))
end

function _complex_horizon_resonance(params::ComplexParameters)
    iszero(params.p) && return (order=0, threshold=true)
    difference = 2im * (params.rplus / params.kappa) * params.p
    order = round(Int, real(difference))
    tolerance = 2048eps(Float64) * max(1.0, abs(difference))
    resonant = order >= 1 &&
        abs(difference - order) <= tolerance
    return (order=resonant ? order : 0, threshold=false)
end

function _complex_amplitude_observable(route::ComplexRoute)
    explicit = _metadata_value(
        route.metadata, :amplitude_observable, nothing)
    explicit === nothing || return explicit
    resonance = _complex_horizon_resonance(route.params)
    resonance.threshold &&
        return :threshold_log_connection_coefficients
    resonance.order > 0 &&
        return :resonant_horizon_connection_coefficients
    return :plane_wave_amplitudes
end

function _complex_gsn_wrapper(route::ComplexRoute)
    root = _root_module()
    p = route.params
    backend = _metadata_value(
        route.metadata, :backend, :p_equation_complex_contour)
    rsmp = _public_match_rstar(route, backend)
    N = route.plan.controls.N
    order = N === nothing ? missing : Int(N)
    horizon_resonance = _complex_horizon_resonance(route.params)
    matching = route.plan.matching
    selected_x = _metadata_value(matching, :selected_x, missing)
    match_policy = _metadata_value(
        matching, :policy,
        _metadata_value(route.metadata, :match_policy, :unknown))
    small_frequency_expansion = _metadata_value(route.plan.controls, :small_frequency_expansion, false)
    large_frequency_expansion = _metadata_value(route.plan.controls, :large_frequency_expansion, false)
    horizon_patches = Int(_metadata_value(
        route.metadata, :horizon_patches, 0))
    infinity_patches = Int(_metadata_value(
        route.metadata, :infinity_patches, 0))
    ordinary_patches = max(
        0, route.plan.patch_count - horizon_patches - infinity_patches)
    anchor_backend = backend == :direct_gsn_mst_anchor
    full_mst_backend = backend == :direct_gsn_full_mst
    spin_regime = p.kappa <= 0.02 ? :near_extreme : :regular
    numerical_metadata = (
        backend,
        amplitude_backend=_metadata_value(
            route.metadata, :amplitude_backend, :complex_rational),
        branch=route.branch,
        route_plan=route.plan,
        xm=selected_x,
        patch_count=route.plan.patch_count,
        small_frequency_expansion,
        large_frequency_expansion,
        frequency_regime=full_mst_backend ? :mst_full :
            anchor_backend ? :mst_anchor : :complex,
        frequency_reason=backend,
        spin_regime,
        spin_reason=spin_regime,
        endpoint_representation=full_mst_backend ? :physical_mst :
            anchor_backend ? :mst : backend,
        propagation_representation=full_mst_backend ? :direct_mst :
            anchor_backend ? :ordinary_taylor : :complex_rational,
        matching_representation=match_policy,
        control_tolerance=route.plan.controls.tolerance,
        horizon_endpoint_score=_metadata_value(
            route.metadata, :horizon_endpoint_score, missing),
        infinity_endpoint_score=_metadata_value(
            route.metadata, :infinity_endpoint_score, missing),
        propagation_score=_metadata_value(
            route.metadata, :propagation_score, missing),
        horizon_path_kind=_metadata_value(
            route.metadata, :horizon_path_kind, missing),
        horizon_coordinate_steps=_metadata_value(
            route.metadata, :horizon_coordinate_steps, missing),
        infinity_horizon_winding=_metadata_value(
            route.metadata, :infinity_horizon_winding, 0),
        horizon_monodromy=_metadata_value(
            route.metadata, :horizon_monodromy, one(ComplexF64)),
        eikonal_candidate=false,
        horizon_order=order,
        ordinary_order=order,
        infinity_out_order=order,
        infinity_in_order=order,
        ordinary_patches,
        eikonal_patches=0,
        endpoint_patches=0,
        horizon_patches,
        infinity_patches,
        solution_scale=1.0,
        unit_pair_available=true,
        unit_state_available=true,
        normalization=:unit_gsn_transmission,
        horizon_in_representation=_metadata_value(
            route.metadata, :horizon_in_representation, missing),
        horizon_out_representation=_metadata_value(
            route.metadata, :horizon_out_representation, missing),
        horizon_resonance_order=_metadata_value(
            route.metadata, :horizon_resonance_order,
            horizon_resonance.order),
        amplitude_observable=_complex_amplitude_observable(route),
        fallback_from=_metadata_value(
            route.metadata, :fallback_from, missing),
        fallback_trigger=_metadata_value(
            route.metadata, :fallback_trigger, missing),
        rational_build_error=_metadata_value(
            route.metadata, :rational_build_error, missing),
        mst_nu=_metadata_value(route.metadata, :mst_nu, missing),
        mst_nu_offset=_metadata_value(
            route.metadata, :mst_nu_offset, missing),
        mst_representation=_metadata_value(
            route.metadata, :mst_representation, missing),
        mst_amplitude_certificate=_metadata_value(
            route.metadata, :mst_amplitude_certificate, missing),
        mst_amplitude_precision_bits=_metadata_value(
            route.metadata, :mst_amplitude_precision_bits, missing),
        mst_nu_residual=_metadata_value(
            route.metadata, :mst_nu_residual, missing),
        mst_amplitude_truncation=_metadata_value(
            route.metadata, :mst_amplitude_truncation, missing),
        mst_amplitude_condition=_metadata_value(
            route.metadata, :mst_amplitude_condition, missing),
        mst_anchor_error=_metadata_value(
            route.metadata, :mst_anchor_error, missing),
        mst_anchor_residual=_metadata_value(
            route.metadata, :mst_anchor_residual, missing),
        mst_anchor_tail=_metadata_value(
            route.metadata, :mst_anchor_tail, missing),
        derivative_coordinate=:rstar,
        stored_normalization=:unit_gsn_transmission,
    )
    mode = getfield(root, :Mode)(p.s, p.l, p.m, p.a, p.omega, p.lambda)
    return getfield(root, :GSNRadialFunction)(
        mode,
        _boundary_condition(route.branch),
        missing,
        missing,
        rsmp,
        order,
        order,
        route.transmission,
        route.incidence,
        route.reflection,
        numerical_metadata,
        missing,
        complex_gsn_solution_rstar(route),
        getfield(root, :UNIT_GSN_TRANS),
        "ISEM",
    )
end

function complex_gsn_radial_function(route::ComplexRoute)
    return _complex_gsn_wrapper(route)
end

function _teukolsky_scale(route::ComplexRoute, normalize)
    root = _root_module()
    if normalize == :teukolsky
        scale = route.teukolsky_transmission
        _finite_complex(scale) && !iszero(scale) || throw(ArgumentError(
            "unit-Teukolsky normalization requires finite nonzero transmission."))
        return scale, getfield(root, :UNIT_TEUKOLSKY_TRANS)
    elseif normalize == :gsn
        return one(route.transmission), getfield(root, :UNIT_GSN_TRANS)
    end
    throw(ArgumentError("normalize must be :teukolsky or :gsn."))
end

function route_teukolsky_solution(
    route::ComplexRoute;
    normalize=:teukolsky,
    scale=nothing,
    converter=nothing,
)
    p = route.params
    scale_value = scale === nothing ? first(_teukolsky_scale(route, normalize)) : scale
    pair_converter = converter === nothing ?
        _teukolsky_converter(p) : converter
    pair = ComplexTeukolskyPair(
        route, pair_converter, inv(scale_value))
    transformation = _matching_transformation_module()
    tuple_value = ComplexTeukolskyTuple(
        p, pair, getfield(transformation, :d2R))
    radial_value = TeukolskyRadialEvaluator(pair)
    return TeukolskySolution(tuple_value, pair, radial_value)
end

function teukolsky_radial_function(
    route::ComplexRoute;
    normalize=:teukolsky,
)
    root = _root_module()
    scale, convention = _teukolsky_scale(route, normalize)
    solution = route_teukolsky_solution(
        route; normalize, scale)
    p = route.params
    mode = getfield(root, :Mode)(
        p.s, p.l, p.m, p.a, p.omega, p.lambda)
    return getfield(root, :TeukolskyRadialFunction)(
        mode,
        _boundary_condition(route.branch),
        route.teukolsky_transmission / scale,
        route.teukolsky_incidence / scale,
        route.teukolsky_reflection / scale,
        route.p_solution,
        complex_gsn_radial_function(route),
        solution,
        convention,
    )
end

y_branch_supported(route::ComplexRoute) =
    (route.params.s == -2 && route.branch == :IN) ||
    (route.params.s == 2 && route.branch == :UP)

function _check_y_branch(route::ComplexRoute)
    y_branch_supported(route) && return nothing
    error("ISEM Y transformation supports only s = -2 with IN and s = +2 with UP.")
end

function _y_amplitudes(route::ComplexRoute)
    _check_y_branch(route)
    return (
        trans=route.teukolsky_transmission,
        incidence=route.teukolsky_incidence,
        reflection=route.teukolsky_reflection,
    )
end

function route_y_solution(route::ComplexRoute; amplitudes=nothing)
    _check_y_branch(route)
    p = route.params
    iszero(p.omega) && error("ISEM Y transformation requires nonzero omega.")
    amps = amplitudes === nothing ? _y_amplitudes(route) : amplitudes
    transformation = _matching_transformation_module()
    matrix = r ->
        Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(
            p.s, p.m, p.a, p.omega, p.lambda, r)
    return getfield(transformation, :GSN_to_Y_solution_from_matrix)(
        r -> rstar_from_r(p.a, r),
        matrix,
        complex_gsn_solution_rstar(route),
        p.s,
        p.m,
        p.a,
        p.omega,
        p.lambda,
        amps.trans,
    )
end

function y_radial_function(route::ComplexRoute)
    _check_y_branch(route)
    root = _root_module()
    isem = _isem_module()
    y_type = getfield(isem, :YRadialFunction)
    amps = _y_amplitudes(route)
    y_solution = route_y_solution(route; amplitudes=amps)
    x_state = complex_gsn_solution_rstar(route)
    x_solution = r -> x_state(rstar_from_r(route.params.a, r))[1]
    y_scalar = r -> y_solution(r)[1]
    p = route.params
    mode = getfield(root, :Mode)(
        p.s, p.l, p.m, p.a, p.omega, p.lambda)
    p_solution = route.p_solution isa Function ? route.p_solution : missing
    return y_type(
        mode,
        _boundary_condition(route.branch),
        amps.trans,
        amps.incidence,
        amps.reflection,
        p_solution,
        missing,
        x_solution,
        y_scalar,
        y_solution,
        getfield(root, :UNIT_TEUKOLSKY_TRANS),
    )
end

complex_gsn_radial_function(
    route::ConjugatedRoute{R},
) where {R<:ComplexRoute} = gsn_radial_function(route)

"""
    complex_radial_route(s, l, m, a, omega, branch; kwargs...)

Build a complex-frequency route with `complex_route` (keywords
are passed through) and return it as a `GSNRadialFunction`.
"""
function complex_radial_route(s::Int, l::Int, m::Int, a, omega, branch; kwargs...)
    return complex_gsn_radial_function(
        complex_route(s, l, m, a, omega, branch; kwargs...),
    )
end

function radial_route(
    s::Integer,
    l::Integer,
    m::Integer,
    a,
    omega::Complex,
    branch=:IN;
    lambda=nothing,
    nu=nothing,
    controls=nothing,
    small_frequency_expansion=:auto,
    large_frequency_expansion=:auto,
    kwargs...,
)
    iszero(imag(omega)) && return radial_route(
        s, l, m, a, real(omega), branch;
        lambda, nu, controls, small_frequency_expansion, large_frequency_expansion, kwargs...)
    complex_small_frequency_expansion = small_frequency_expansion === :auto ? nothing : small_frequency_expansion
    complex_large_frequency_expansion = large_frequency_expansion === :auto ? nothing : large_frequency_expansion
    omegac = ComplexF64(omega)
    default_controls = lambda === nothing && nu === nothing &&
        controls === nothing && small_frequency_expansion === :auto && large_frequency_expansion === :auto &&
        all(value === nothing for value in values(kwargs))
    # Re omega < 0: the exact reflection R_m(omega) = conj R_{-m}(-conj(omega))
    # (the MST closed forms hold only for Re epsilon > 0); the sign of the
    # zero of Re omega picks the side of the negative imaginary axis.
    if default_controls && signbit(real(omegac))
        source_omega = -conj(omegac)
        source = complex_route(
            Int(s), Int(l), -Int(m), a, source_omega, branch;
            small_frequency_expansion=complex_small_frequency_expansion,
            large_frequency_expansion=complex_large_frequency_expansion,
        )
        target_params = complex_parameters(
            Int(s),
            Int(l),
            Int(m),
            a,
            omegac,
            conj(source.params.lambda),
        )
        return ConjugatedRoute(source, target_params)
    end
    return complex_route(
        Int(s), Int(l), Int(m), a, omegac, branch;
        lambda,
        nu,
        controls,
        small_frequency_expansion=complex_small_frequency_expansion,
        large_frequency_expansion=complex_large_frequency_expansion,
        kwargs...,
    )
end

gsn_radial_function(route::ComplexRoute) =
    complex_gsn_radial_function(route)


# Relative errors (incidence, reflection) of the route's amplitudes: those
# recorded by the two-ray build, or, when the analytic (MST) amplitudes
# replaced them, the analytic estimate (condition and the agreement of their
# representations, rescaled to each amplitude, plus the running error bound of
# their sums where recorded); Inf when neither is recorded.
function amplitude_errors(route::ComplexRoute)
    metadata = route.metadata
    if _metadata_value(metadata, :amplitude_backend, nothing) ==
            :complex_mst_analytic
        analytic = max(
            eps(Float64) * _metadata_value(
                metadata, :mst_amplitude_max_condition, Inf),
            _metadata_value(metadata, :mst_amplitude_nearest_agreement, Inf),
            _metadata_value(metadata,
                :mst_amplitude_truncation_agreement, Inf))
        bound = _metadata_value(metadata, :mst_amplitude_errors, (NaN, NaN))
        own = all(isfinite, bound) ? bound : (0.0, 0.0)
        # Each amplitude's own relative distance from the other evaluations
        # where recorded; otherwise the pair distance (units of the larger
        # amplitude) rescaled to it. Plus the running error of its own sums.
        components = _metadata_value(metadata,
            :mst_amplitude_component_agreement, (NaN, NaN))
        all(isfinite, components) && return Float64.(own .+ max.(
            eps(Float64) * _metadata_value(
                metadata, :mst_amplitude_max_condition, Inf), components))
        amplitudes = abs.((route.incidence, route.reflection))
        return Float64.(own .+ analytic .* maximum(amplitudes) ./ amplitudes)
    end
    return (Float64(_metadata_value(metadata, :incidence_error, Inf)),
        Float64(_metadata_value(metadata, :reflection_error, Inf)))
end

end
