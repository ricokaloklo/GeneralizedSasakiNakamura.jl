module Matching

using ..CoefficientTables: CoefficientSet, gsn_coefficients
using ..Iteration:
    Basis,
    Patch,
    ScaledPatch,
    LogScaledPatch,
    LogScaledState,
    LogScaledBasis,
    IterationScratch,
    Truncation,
    basis_patch_count,
    basis_state,
    basis_value,
    endpoint_seed_scale,
    infinity_exponents,
    iterate_scaled_y,
    iterate_logscaled_from_state,
    iterate_from_state,
    iterate_from_infinity,
    iterate_pair_from_zero,
    iterate_from_zero,
    iterate_pair_from_infinity,
    infinity_seed_wronskian,
    select_infinity_endpoint,
    logscaled_state,
    logscaled_basis_state,
    materialize_logscaled_state,
    _SMALL_FREQUENCY_EXPANSION_BRIDGE_TARGET_Y,
    _endpoint_pair,
    _endpoint_state,
    _ENDPOINT_SEED_TARGET,
    _scaled_y_infinity_selection,
    _INFINITY_ENDPOINT_TARGET
using ..Eikonal: EikonalPatch, eikonal_preflight
using ..NearExtreme: near_extreme_prepare, near_extreme_selected
using ..ManoSuzukiTakasugiInfinity: ManoSuzukiTakasugiCertificateError, abel_denominator
using ..ManoSuzukiTakasugiInfinity: build_mano_suzuki_takasugi_eval_plan, mano_suzuki_takasugi_infinity_basis
using ..ManoSuzukiTakasugiInfinity: mano_suzuki_takasugi_infinity_pair, build_mano_suzuki_takasugi_plan, mano_suzuki_takasugi_route_state
using ..ManoSuzukiTakasugiInfinity: mano_suzuki_takasugi_logscaled_seed, mano_suzuki_takasugi_pin_state
using ..ManoSuzukiTakasugiInfinity: mano_suzuki_takasugi_scaled_state
using ..ManoSuzukiTakasugiInfinity: _pin_unit_scale
using ..Parameters:
    GSNControls,
    GSNParameters,
    SMALL_FREQUENCY_EXPANSION_FREQUENCY_LIMIT,
    SMALL_FREQUENCY_EXPANSION_TRANSITION_LIMIT,
    _ordinary_fallback,
    gsn_controls,
    gsn_parameters,
    horizon_tail,
    select_spin_regime

export Route, ConjugatedRoute, RoutePlan
export match_route, radial_route
export evaluate_route, route_state, route_patch_count, route_plan_of, route_truncations

const DENOMINATOR_REL_MIN = 1.0e-10
const LOW_FREQUENCY_MIN_ORDER = 40
const LOW_FREQUENCY_TAIL_HORIZON_ORDER = 80
const LOW_FREQUENCY_MATCH_X = 0.9
const LOW_FREQUENCY_AXISYMMETRIC_MATCH_X = 0.99
const ORDINARY_IN_MATCH_X = 0.7
const ORDINARY_UP_MATCH_X = 0.3
const BRIDGE_CONDITION_LIMIT = 1.0e12
# Diagnostic switch: ISEM_PHYSICAL_BRIDGE=0 disables the physical bridge (_needs_bridge)
# and the small-frequency endpoint bridge, so routes evaluate the matched bases directly.
const PHYSICAL_BRIDGE = get(ENV, "ISEM_PHYSICAL_BRIDGE", "1") != "0"
const PUBLIC_EVALUATION_ENDPOINT_X = 0.999995
const MATCHING_SCALE_LOG_MARGIN = 32.0
const SMALL_FREQUENCY_EXPANSION_ENDPOINT_HANDOFF_TOLERANCE = 1.0e-10
const SMALL_FREQUENCY_EXPANSION_UP_TRIAL_STEPS = 128

struct EndpointPlan{S}
    selection::S
    scale::ComplexF64
    seed_x::Float64
end

struct ScaledBasis{B}
    basis::B
    scale::ComplexF64
end

"""
    Route

Matched ISEM GSN solution for one branch (`:IN` or `:UP`): the horizon and
infinity bases, the incidence and reflection coefficients at `match_x`, a common
`solution_scale`, and optional bridges, endpoint plan and MST evaluation plan.
Evaluate with `evaluate_route`, `route_state`, or by calling the route on `x`.
"""
struct Route{P,C,HI,HO,IO,II,B,EB,EP,M}
    branch::Symbol
    params::P
    controls::GSNControls
    coefficients::C
    match_x::Float64
    infinity_endpoint_y::Float64
    infinity_out_order::Int
    infinity_in_order::Int
    incidence::ComplexF64
    reflection::ComplexF64
    solution_scale::Float64
    horizon_in::HI
    horizon_out::HO
    infinity_out::IO
    infinity_in::II
    bridge::B
    endpoint_bridge::EB
    endpoint_plan::EP
    mst_plan::M
end

"""
    ConjugatedRoute

Route for a real negative frequency: wraps the route at (-m, -omega) and returns
complex-conjugated values; `params` holds the requested parameters.
"""
struct ConjugatedRoute{R,P}
    route::R
    params::P
end

function Base.getproperty(route::ConjugatedRoute, name::Symbol)
    name === :route && return getfield(route, :route)
    name === :params && return getfield(route, :params)
    return getproperty(getfield(route, :route), name)
end

"""
    RoutePlan

Diagnostic summary of a route returned by `route_plan_of`.
"""
struct RoutePlan
    frequency::Symbol
    frequency_reason::Symbol
    spin::Symbol
    spin_reason::Symbol
    endpoint::Symbol
    propagation::Symbol
    matching::Symbol
    match_x::Float64
    horizon_order::Int
    ordinary_order::Int
    infinity_order::Int
    eikonal_candidate::Bool
    ordinary_patches::Int
    eikonal_patches::Int
    endpoint_patches::Int
    horizon_patches::Int
    infinity_patches::Int
    angular_compression::Float64
    phase_compression::Float64
end

function _normalized_route_branch(branch)
    text = uppercase(String(Symbol(branch)))
    text == "IN" && return :IN
    text == "UP" && return :UP
    throw(ArgumentError("ISEM GSN route branch must be :IN or :UP."))
end

_basis_tuple(basis::Basis) = (X=basis.match_X, dXdx=basis.match_dXdx)
_logscaled_match_state(basis::LogScaledBasis) = basis.match_state
_logscaled_match_state(basis::Basis) = logscaled_state(
    basis.match_X, basis.match_dXdx)
_logscaled_match_state(state::LogScaledState) = state
_logscaled_match_state(state) = logscaled_state(state.X, state.dXdx)
_scaled_state(state, scale) = (
    X=ComplexF64(scale * state.X),
    dXdx=ComplexF64(scale * state.dXdx),
)

function _required_basis(basis, name)
    basis === nothing && error("missing ISEM GSN basis $name")
    return basis
end

# Controls for small_frequency_expansion routes: in the transition band only the default IN match
# point moves; below SMALL_FREQUENCY_EXPANSION_FREQUENCY_LIMIT the default match point moves out, the
# orders are raised to at least LOW_FREQUENCY_MIN_ORDER, and the scaled-y
# infinity basis is selected.
function _with_low_frequency_defaults(ctrl::GSNControls, params, route_branch::Symbol)
    ctrl.small_frequency_expansion || return ctrl
    omega = abs(getproperty(params, :omega))
    0 < omega <= SMALL_FREQUENCY_EXPANSION_TRANSITION_LIMIT || return ctrl
    if omega > SMALL_FREQUENCY_EXPANSION_FREQUENCY_LIMIT
        route_branch == :IN || return ctrl
        getproperty(ctrl, :source) == :default || return ctrl
        return GSNControls(
            ctrl.horizon_order,
            ctrl.ordinary_order,
            ctrl.infinity_order,
            LOW_FREQUENCY_MATCH_X,
            ctrl.tolerance,
            ctrl.basis,
            ctrl.endpoint_basis,
            ctrl.ordinary_basis,
            ctrl.small_frequency_expansion,
            ctrl.large_frequency_expansion,
            ctrl.frequency_selection,
            :small_frequency_expansion_transition_match,
        )
    end
    match_x = if getproperty(ctrl, :source) == :default
        spin = select_spin_regime(params)
        spin.regime == :near_extremal && iszero(getproperty(params, :m)) ?
            LOW_FREQUENCY_AXISYMMETRIC_MATCH_X : LOW_FREQUENCY_MATCH_X
    else
        ctrl.match_x
    end
    0 < match_x < 1 || return ctrl
    horizon_order = horizon_tail(params) ?
        LOW_FREQUENCY_TAIL_HORIZON_ORDER : LOW_FREQUENCY_MIN_ORDER
    return GSNControls(
        max(ctrl.horizon_order, horizon_order),
        max(ctrl.ordinary_order, LOW_FREQUENCY_MIN_ORDER),
        ctrl.infinity_order,
        match_x,
        ctrl.tolerance,
        :low_frequency_scaled_y,
        ctrl.endpoint_basis,
        ctrl.ordinary_basis,
        ctrl.small_frequency_expansion,
        ctrl.large_frequency_expansion,
        ctrl.frequency_selection,
        :low_frequency_default,
    )
end

# Default controls at an ordinary real frequency: match at ORDINARY_IN_MATCH_X
# for IN and ORDINARY_UP_MATCH_X for UP.
function _with_ordinary_defaults(
    ctrl::GSNControls,
    params,
    route_branch::Symbol,
)
    getproperty(ctrl, :source) == :default || return ctrl
    ctrl.frequency_selection.regime == :ordinary || return ctrl
    isreal(getproperty(params, :omega)) || return ctrl
    match_x = route_branch == :IN ?
        ORDINARY_IN_MATCH_X : ORDINARY_UP_MATCH_X
    return GSNControls(
        ctrl.horizon_order,
        ctrl.ordinary_order,
        ctrl.infinity_order,
        match_x,
        ctrl.tolerance,
        ctrl.basis,
        ctrl.endpoint_basis,
        ctrl.ordinary_basis,
        ctrl.small_frequency_expansion,
        ctrl.large_frequency_expansion,
        ctrl.frequency_selection,
        :ordinary_branch_match,
    )
end

@inline _finite_complex(value) =
    isfinite(real(value)) && isfinite(imag(value))

function _logabs(value)
    re = abs(real(value))
    im = abs(imag(value))
    scale = max(re, im)
    iszero(scale) && return -Inf
    isfinite(scale) || return Inf
    return log(scale) + 0.5 * log((re / scale)^2 + (im / scale)^2)
end

function _logadd(a, b)
    top = max(a, b)
    top == -Inf && return -Inf
    isfinite(top) || return top
    return top + log(exp(a - top) + exp(b - top))
end

# Log of the cancellation in c1 b1 + c2 b2 = target, for value and derivative.
function _eval_condition(target, coefficient1, basis1, coefficient2, basis2)
    value_terms = _logadd(
        _logabs(coefficient1) + _logabs(basis1.X),
        _logabs(coefficient2) + _logabs(basis2.X),
    )
    derivative_terms = _logadd(
        _logabs(coefficient1) + _logabs(basis1.dXdx),
        _logabs(coefficient2) + _logabs(basis2.dXdx),
    )
    value_scale = _logabs(target.X)
    derivative_scale = _logabs(target.dXdx)
    value_condition = value_scale == -Inf ? Inf : value_terms - value_scale
    derivative_condition = derivative_scale == -Inf ? Inf :
        derivative_terms - derivative_scale
    return max(value_condition, derivative_condition)
end

function _needs_bridge(target, coefficient1, basis1, coefficient2, basis2;
    tolerance=nothing)
    PHYSICAL_BRIDGE || return false
    condition = _eval_condition(
        target, coefficient1, basis1, coefficient2, basis2)
    log_limit = log(BRIDGE_CONDITION_LIMIT)
    if tolerance !== nothing
        # This chooses an already-supported matched-state bridge for the
        # low-frequency UP route; it does not relax the MST sum certificate.
        log_limit = min(log_limit, log(tolerance) - log(eps(Float64)))
    end
    return condition > log_limit
end

# When the matched combination loses more than BRIDGE_CONDITION_LIMIT to
# cancellation, propagate the matched state itself from match_x to target_x
# and use it there instead of the combination; nothing otherwise.
function _build_bridge(
    coefficients,
    ctrl,
    scratch,
    branch,
    match_x,
    target_x,
    target,
    coefficient1,
    basis1,
    coefficient2,
    basis2,
)
    _needs_bridge(target, coefficient1, basis1, coefficient2, basis2;
        tolerance=(branch == :UP && ctrl.small_frequency_expansion ?
            ctrl.tolerance : nothing)) || return nothing
    abs(target_x - match_x) > 100eps(Float64) || return nothing
    kind = branch == :IN ? :horizon_in : :infinity_out
    try
        return iterate_from_state(
            coefficients,
            kind,
            match_x,
            target.X,
            target.dXdx,
            target_x;
            controls=ctrl,
            scratch=scratch,
            allow_eikonal=true,
            allow_partial=true,
        )
    catch err
        message = sprint(showerror, err)
        startswith(message, "ISEM GSN hybrid recurrence failed") || rethrow()
        return nothing
    end
end

# Normalized infinity endpoint series of one branch at x, or at the endpoint
# coordinate computed from `radius` when given.
function _endpoint_branch_state(coefficients, selection, branch, x; radius=nothing)
    kind = branch == :in ? :infinity_in : :infinity_out
    coeffs = branch == :in ? selection.in_coeffs : selection.out_coeffs
    order = branch == :in ? selection.in_order : selection.out_order
    rho, sigma = infinity_exponents(coefficients.params, branch)
    value, derivative = if radius === nothing
        _endpoint_state(kind, @view(coeffs[1:(order + 1)]), rho,
            sigma / selection.variable_scale, x, selection.variable_scale)
    else
        kappa = coefficients.params.kappa
        y = 2kappa / (radius - 1 + kappa)
        _endpoint_pair(kind, @view(coeffs[1:(order + 1)]), rho,
            sigma / selection.variable_scale, y, selection.variable_scale)
    end
    scale = endpoint_seed_scale(kind, coefficients.params)
    return (X=scale * value, dXdx=scale * derivative)
end

function _physical_endpoint_state(
    coefficients,
    route_branch,
    selection,
    x,
    incidence,
    reflection,
    solution_scale,
    scale=1.0 + 0.0im;
    radius=nothing,
)
    outgoing = _endpoint_branch_state(coefficients, selection, :out, x; radius)
    if route_branch == :UP
        return (
            X=scale * solution_scale * outgoing.X,
            dXdx=scale * solution_scale * outgoing.dXdx,
        )
    end
    incoming = _endpoint_branch_state(coefficients, selection, :in, x; radius)
    return (
        X=scale * (incidence * incoming.X + reflection * outgoing.X),
        dXdx=scale * (incidence * incoming.dXdx +
            reflection * outgoing.dXdx),
    )
end

# Least-squares complex scale c minimizing |reference - c * candidate| over
# (X, dXdx); nothing if either state is zero or not finite.
function _fit_state_scale(reference, candidate)
    reference_scale = max(abs(reference.X), abs(reference.dXdx))
    candidate_scale = max(abs(candidate.X), abs(candidate.dXdx))
    reference_scale > 0 && candidate_scale > 0 || return nothing
    isfinite(reference_scale) && isfinite(candidate_scale) || return nothing
    ref_x = reference.X / reference_scale
    ref_dx = reference.dXdx / reference_scale
    candidate_x = candidate.X / candidate_scale
    candidate_dx = candidate.dXdx / candidate_scale
    denominator = abs2(candidate_x) + abs2(candidate_dx)
    denominator > 0 && isfinite(denominator) || return nothing
    unit_scale = (conj(candidate_x) * ref_x +
        conj(candidate_dx) * ref_dx) / denominator
    scale = (reference_scale / candidate_scale) * unit_scale
    return _finite_complex(scale) ? ComplexF64(scale) : nothing
end

function _state_error(left, right)
    scale = max(abs(left.X), abs(left.dXdx), abs(right.X), abs(right.dXdx),
        floatmin(Float64))
    isfinite(scale) || return Inf
    return max(abs(left.X - right.X), abs(left.dXdx - right.dXdx)) / scale
end

function _small_frequency_expansion_endpoint_selection(coefficients, ctrl, match_x)
    try
        return select_infinity_endpoint(
            coefficients, match_x; controls=ctrl)
    catch err
        message = sprint(showerror, err)
        startswith(message, "ISEM GSN infinity") || rethrow()
        return nothing
    end
end

# Small-frequency bridge: propagate the matched state outward in scaled y and
# stop once it agrees, up to a fitted scale, with the infinity endpoint series.
# Returns (bridge or nothing, endpoint plan or nothing).
function _build_small_frequency_expansion_endpoint_bridge(
    coefficients,
    ctrl,
    kind,
    route_branch,
    match_x,
    target,
    incidence,
    reflection,
    solution_scale,
)
    PHYSICAL_BRIDGE && ctrl.small_frequency_expansion || return nothing, nothing
    match_x < PUBLIC_EVALUATION_ENDPOINT_X || return nothing, nothing
    omega = abs(getproperty(coefficients.params, :omega))
    omega > 0 || return nothing, nothing
    target_y = (1.0 - PUBLIC_EVALUATION_ENDPOINT_X) / omega
    selection = _small_frequency_expansion_endpoint_selection(coefficients, ctrl, match_x)
    plan = Ref{Any}(nothing)

    function stop_check(x, X, dXdx)
        selection === nothing && return false
        1.0 - x <= selection.endpoint_y * (1 + 20eps(Float64)) ||
            return false
        bridge_state = (X=X, dXdx=dXdx)
        endpoint = _physical_endpoint_state(
            coefficients,
            route_branch,
            selection,
            x,
            incidence,
            reflection,
            solution_scale,
        )
        endpoint_scale = _fit_state_scale(bridge_state, endpoint)
        endpoint_scale === nothing && return false
        aligned = _physical_endpoint_state(
            coefficients,
            route_branch,
            selection,
            x,
            incidence,
            reflection,
            solution_scale,
            endpoint_scale,
        )
        _state_error(bridge_state, aligned) <=
            SMALL_FREQUENCY_EXPANSION_ENDPOINT_HANDOFF_TOLERANCE || return false
        plan[] = EndpointPlan(selection, endpoint_scale, Float64(x))
        return true
    end

    bridge = iterate_scaled_y(
        coefficients,
        kind,
        match_x,
        target.X,
        target.dXdx;
        target_y,
        stop_check,
    )
    return bridge, plan[]
end

# Small-frequency UP: try an inward infinity basis (step budget
# SMALL_FREQUENCY_EXPANSION_UP_TRIAL_STEPS) aligned to the matched state at match_x. Returns
# (scaled basis or nothing, endpoint plan or nothing).
function _small_frequency_expansion_up_basis(coefficients, ctrl, scratch, match_x, target)
    omega = abs(getproperty(coefficients.params, :omega))
    omega > 0 || return nothing, nothing
    target_y = (1.0 - PUBLIC_EVALUATION_ENDPOINT_X) / omega
    extended = target_y > _SMALL_FREQUENCY_EXPANSION_BRIDGE_TARGET_Y
    selection = _small_frequency_expansion_endpoint_selection(coefficients, ctrl, match_x)
    selection === nothing && return nothing, nothing
    # The UP boundary fixes this normalization, even when inward propagation
    # cannot reach the matching point. Use it only in the certified endpoint domain.
    order = max(ctrl.infinity_order, ctrl.ordinary_order)
    boundary = order > selection.requested_order ?
        _scaled_y_infinity_selection(coefficients, order) : selection
    plan = boundary !== nothing && boundary.score <= _INFINITY_ENDPOINT_TARGET ?
        EndpointPlan(boundary, 1.0 + 0.0im, 1.0-boundary.endpoint_y) : nothing
    extended && eikonal_preflight(coefficients.params, :small_frequency_expansion).candidate &&
        return nothing, plan
    basis = try
        iterate_from_infinity(
            coefficients,
            :out;
            controls=ctrl,
            match_x,
            selection,
            scratch,
            certified_steps=extended,
            step_budget=SMALL_FREQUENCY_EXPANSION_UP_TRIAL_STEPS,
        )
    catch err
        startswith(sprint(showerror, err), "ISEM GSN") || rethrow()
        return nothing, plan
    end
    basis === nothing && return nothing, plan
    scale = _fit_state_scale(target, _basis_tuple(basis))
    scale === nothing && return nothing, plan
    aligned = _scaled_state(_basis_tuple(basis), scale)
    _state_error(target, aligned) <= SMALL_FREQUENCY_EXPANSION_ENDPOINT_HANDOFF_TOLERANCE ||
        return nothing, plan
    return ScaledBasis(basis, scale), plan
end

# (a b - c d) / det, recomputed from scaled factors when the direct product
# overflows; errors if the result is outside the Float64 range.
function _scaled_det_ratio(a, b, c, d, det)
    direct = (a * b - c * d) / det
    _finite_complex(direct) && return direct

    left_scale = max(abs(a), abs(c))
    right_scale = max(abs(b), abs(d))
    (iszero(left_scale) || iszero(right_scale)) && return 0.0 + 0.0im

    mantissa = (a / left_scale) * (b / right_scale) -
        (c / left_scale) * (d / right_scale)
    iszero(mantissa) && return 0.0 + 0.0im

    det_scale = abs(det)
    (!isfinite(det_scale) || iszero(det_scale)) &&
        error("singular ISEM GSN scaled matching denominator")
    log_magnitude = log(abs(mantissa)) + log(left_scale) +
        log(right_scale) - log(det_scale)
    log_magnitude <= log(floatmax(Float64)) ||
        error("ISEM GSN matching coefficient exceeds Float64 range")
    phase = angle(mantissa) - angle(det)
    value = exp(log_magnitude) * cis(phase)
    _finite_complex(value) ||
        error("nonfinite ISEM GSN scaled matching coefficient")
    return ComplexF64(value)
end

function _scaled_det_ratio_polar(a, b, c, d, det)
    left_scale = max(abs(a), abs(c))
    right_scale = max(abs(b), abs(d))
    (iszero(left_scale) || iszero(right_scale)) && return (-Inf, 0.0)
    isfinite(left_scale) && isfinite(right_scale) ||
        error("nonfinite ISEM GSN matching numerator scale")

    mantissa = (a / left_scale) * (b / right_scale) -
        (c / left_scale) * (d / right_scale)
    iszero(mantissa) && return (-Inf, 0.0)
    det_scale = abs(det)
    isfinite(det_scale) && !iszero(det_scale) ||
        error("singular ISEM GSN scaled matching denominator")
    log_magnitude = log(abs(mantissa)) + log(left_scale) +
        log(right_scale) - log(det_scale)
    phase = angle(mantissa) - angle(det)
    return log_magnitude, phase
end

function _scaled_polar_value(log_magnitude, phase, log_scale)
    log_magnitude == -Inf && return 0.0 + 0.0im
    value = exp(log_magnitude + log_scale) * cis(phase)
    _finite_complex(value) ||
        error("nonfinite ISEM GSN jointly scaled matching coefficient")
    return ComplexF64(value)
end

# Cramer solution with a given denominator `det`. If a coefficient would come
# near the Float64 range, both are scaled down by a common factor, returned as
# the third value (solution_scale) and applied to the target state.
function _connection_coefficients(target, basis1, basis2, det)
    (!isfinite(real(det)) || !isfinite(imag(det)) || iszero(det)) &&
        error("singular ISEM GSN Abel matching denominator")
    target_log = log(floatmax(Float64)) - MATCHING_SCALE_LOG_MARGIN
    try
        c1 = _scaled_det_ratio(
            target.X, basis2.dXdx, basis2.X, target.dXdx, det)
        c2 = _scaled_det_ratio(
            basis1.X, target.dXdx, target.X, basis1.dXdx, det)
        max(_logabs(c1), _logabs(c2)) <= target_log &&
            return c1, c2, 1.0
    catch err
        sprint(showerror, err) ==
            "ISEM GSN matching coefficient exceeds Float64 range" || rethrow()
    end

    log1, phase1 = _scaled_det_ratio_polar(
        target.X, basis2.dXdx, basis2.X, target.dXdx, det)
    log2, phase2 = _scaled_det_ratio_polar(
        basis1.X, target.dXdx, target.X, basis1.dXdx, det)
    max_log = max(log1, log2)
    isfinite(max_log) || error("nonfinite ISEM GSN matching scale")
    log_scale = min(0.0, target_log - max_log)
    solution_scale = exp(log_scale)
    isfinite(solution_scale) && solution_scale > 0 ||
        error("ISEM GSN matching normalization is outside Float64 range")
    c1 = _scaled_polar_value(log1, phase1, log_scale)
    c2 = _scaled_polar_value(log2, phase2, log_scale)
    return c1, c2, solution_scale
end

function _connection_coefficients_logscaled(
    target::LogScaledState,
    basis1::LogScaledState,
    basis2::LogScaledState,
    det,
)
    _finite_complex(det) && !iszero(det) ||
        error("singular ISEM GSN Abel matching denominator")
    numerator1 = target.X * basis2.dXdx - basis2.X * target.dXdx
    numerator2 = basis1.X * target.dXdx - target.X * basis1.dXdx
    _finite_complex(numerator1) && _finite_complex(numerator2) ||
        error("nonfinite ISEM GSN log-scaled matching numerator")
    iszero(numerator1) && iszero(numerator2) &&
        error("singular ISEM GSN log-scaled matching numerators")

    log_det = log(abs(det))
    phase_det = angle(det)
    log1 = iszero(numerator1) ? -Inf :
        target.log_scale + basis2.log_scale +
        log(abs(numerator1)) - log_det
    log2 = iszero(numerator2) ? -Inf :
        basis1.log_scale + target.log_scale +
        log(abs(numerator2)) - log_det
    phase1 = iszero(numerator1) ? 0.0 : angle(numerator1) - phase_det
    phase2 = iszero(numerator2) ? 0.0 : angle(numerator2) - phase_det
    max_log = max(log1, log2)
    isfinite(max_log) || error("nonfinite ISEM GSN log-scaled matching scale")

    target_log = log(floatmax(Float64)) - MATCHING_SCALE_LOG_MARGIN
    joint_log_scale = min(0.0, target_log - max_log)
    solution_scale = exp(joint_log_scale)
    isfinite(solution_scale) && solution_scale > 0 ||
        error("ISEM GSN log-scaled matching normalization is outside Float64 range")
    coefficient1 = _scaled_polar_value(log1, phase1, joint_log_scale)
    coefficient2 = _scaled_polar_value(log2, phase2, joint_log_scale)
    return coefficient1, coefficient2, solution_scale
end

function _denom_ok(fallback, scale)
    _finite = isfinite(real(fallback)) && isfinite(imag(fallback)) && !iszero(fallback)
    _finite || return false
    return abs(fallback) / max(scale * scale, eps(Float64)) >= DENOMINATOR_REL_MIN
end

function _basis_scale(basis1, basis2)
    return max(
        abs(basis1.X),
        abs(basis1.dXdx),
        abs(basis2.X),
        abs(basis2.dXdx),
        eps(Float64),
    )
end

@inline function _small_frequency_expansion_auto_fallback(ctrl::GSNControls, err)
    return ctrl.frequency_selection.small_frequency_expansion_request == :auto &&
        err isa ManoSuzukiTakasugiCertificateError
end

# Wronskian of the two bases at match_x: the Abel denominator, or for the
# horizon pair the numerical `fallback` when it is well scaled and more
# accurate than the quadrature, or when the Abel evaluation fails.
function _match_denominator(
    coefficients,
    endpoint::Symbol,
    match_x::Float64,
    ctrl::GSNControls,
    fallback;
    fallback_scale=NaN,
)
    denom = abel_denominator(coefficients, endpoint, match_x)
    direct_ok = endpoint == :horizon && _denom_ok(fallback, Float64(fallback_scale))
    if direct_ok && denom.method == :rational_quadrature_depth
        direct_rel = abs(fallback) /
            max(Float64(fallback_scale)^2, eps(Float64))
        direct_error = 1024eps(Float64) / max(direct_rel, eps(Float64))
        direct_error < denom.tail && return fallback
    end
    if denom.status != "OK"
        direct_ok && return fallback
        error("ISEM GSN $(endpoint) Abel denominator failed: $(denom.status)")
    end
    return denom.denominator
end

# Relative difference between the numerical Wronskian of the two infinity seeds
# and the closed-form Abel denominator at the seed point (Inf if unavailable).
function _infinity_abel_certificate(
    coefficients,
    endpoint_y,
    out_coeffs,
    out_order,
    in_coeffs,
    in_order,
    variable_scale,
)
    closed = abel_denominator(
        coefficients, :infinity, 1.0 - endpoint_y)
    closed.status == "OK" || return Inf
    exact = closed.denominator
    iszero(exact) && return Inf
    numerical = infinity_seed_wronskian(
        coefficients,
        endpoint_y,
        out_coeffs,
        out_order,
        in_coeffs,
        in_order,
        variable_scale,
    )
    ratio = numerical / exact
    isfinite(real(ratio)) && isfinite(imag(ratio)) || return Inf
    return abs(ratio - 1)
end

function _infinity_basis(
    coefficients,
    branch::Symbol,
    selection,
    match_x::Float64,
    ctrl::GSNControls,
    scratch,
    mst_plan=nothing,
)
    if ctrl.small_frequency_expansion
        return mano_suzuki_takasugi_infinity_basis(
            coefficients,
            branch,
            match_x;
            controls=ctrl,
            scratch=scratch,
            plan=mst_plan,
        )
    end
    return iterate_from_infinity(coefficients, branch; controls=ctrl,
        match_x=match_x, selection=selection, scratch=scratch)
end

function _logscaled_mano_suzuki_takasugi_basis(
    coefficients,
    branch::Symbol,
    match_x::Float64,
    ctrl,
    scratch,
    mst_plan,
)
    normalized = branch == :in ? :in : :out
    kind = normalized == :in ? :infinity_in : :infinity_out
    seed, _ = mano_suzuki_takasugi_logscaled_seed(
        coefficients, mst_plan, normalized)
    return iterate_logscaled_from_state(
        coefficients,
        kind,
        mst_plan.seed_x,
        seed,
        match_x;
        controls=ctrl,
        scratch=scratch,
    )
end

function _logscaled_mano_suzuki_takasugi_pair(
    coefficients,
    match_x::Float64,
    ctrl,
    scratch,
    mst_plan,
)
    incoming = _logscaled_mano_suzuki_takasugi_basis(
        coefficients, :in, match_x, ctrl, scratch, mst_plan)
    outgoing = _logscaled_mano_suzuki_takasugi_basis(
        coefficients, :out, match_x, ctrl, scratch, mst_plan)
    return incoming, outgoing
end

function _physical_handoff_state(
    coefficients,
    route_branch,
    bridge,
    endpoint_plan,
    x,
    incidence,
    reflection,
    solution_scale,
)
    if bridge !== nothing && _bridge_covers(bridge, x)
        value, derivative = _bridge_state(bridge, x)
        return (X=value, dXdx=derivative)
    end
    if endpoint_plan !== nothing && _endpoint_plan_covers(endpoint_plan, x)
        return _physical_endpoint_state(
            coefficients,
            route_branch,
            endpoint_plan.selection,
            x,
            incidence,
            reflection,
            solution_scale,
            endpoint_plan.scale,
        )
    end
    error("ISEM GSN log-scaled route lacks a physical MST handoff state")
end

# MST evaluation plan for a log-scaled route: fix the MST pin (IN) or outgoing
# (UP) normalization against the physical state at a handoff point and check
# the agreement to SMALL_FREQUENCY_EXPANSION_ENDPOINT_HANDOFF_TOLERANCE.
function _logscaled_mano_suzuki_takasugi_eval_plan(
    coefficients,
    route_branch,
    mst_plan,
    bridge,
    endpoint_plan,
    incidence,
    reflection,
    solution_scale,
    logscaled_basis=nothing,
)
    use_logscaled_basis = route_branch == :UP &&
        logscaled_basis isa LogScaledBasis
    desired_x = min(Float64(mst_plan.seed_x), PUBLIC_EVALUATION_ENDPOINT_X)
    handoff_x = if use_logscaled_basis
        Float64(mst_plan.seed_x)
    elseif bridge !== nothing && _bridge_covers(bridge, desired_x)
        desired_x
    elseif endpoint_plan !== nothing &&
            _endpoint_plan_covers(endpoint_plan, desired_x)
        desired_x
    elseif bridge !== nothing && !isempty(_route_basis(bridge).patches)
        Float64(_route_basis(bridge).patches[end].next_x)
    else
        desired_x
    end
    physical = if use_logscaled_basis
        materialize_logscaled_state(
            logscaled_basis_state(logscaled_basis, handoff_x);
            scale=solution_scale,
        )
    else
        _physical_handoff_state(
            coefficients,
            route_branch,
            bridge,
            endpoint_plan,
            handoff_x,
            incidence,
            reflection,
            solution_scale,
        )
    end
    if route_branch == :IN
        candidate = mano_suzuki_takasugi_pin_state(
            coefficients, mst_plan, handoff_x; scale=1)
        pin_scale = _fit_state_scale(physical, candidate)
        pin_scale === nothing &&
            error("ISEM GSN physical-in MST handoff scale is unavailable")
        aligned = mano_suzuki_takasugi_pin_state(
            coefficients, mst_plan, handoff_x; scale=pin_scale)
        _state_error(physical, aligned) <= SMALL_FREQUENCY_EXPANSION_ENDPOINT_HANDOFF_TOLERANCE ||
            error("ISEM GSN physical-in MST handoff failed validation")
        return build_mano_suzuki_takasugi_eval_plan(
            mst_plan;
            pin_scale=pin_scale,
            logscaled=true,
            handoff_x=handoff_x,
        )
    end

    aligned = if use_logscaled_basis
        seed, estimated_relerr = mano_suzuki_takasugi_logscaled_seed(
            coefficients, mst_plan, :out)
        materialized = materialize_logscaled_state(
            seed; scale=solution_scale)
        (X=materialized.X, dXdx=materialized.dXdx,
            estimated_relerr=estimated_relerr)
    else
        mano_suzuki_takasugi_scaled_state(
            coefficients,
            mst_plan,
            :out,
            handoff_x;
            scale=solution_scale,
        )
    end
    _state_error(physical, aligned) <= SMALL_FREQUENCY_EXPANSION_ENDPOINT_HANDOFF_TOLERANCE ||
        error("ISEM GSN outgoing MST handoff failed validation")
    return build_mano_suzuki_takasugi_eval_plan(
        mst_plan;
        out_scale=solution_scale,
        logscaled=true,
        handoff_x=handoff_x,
    )
end

function _public_exterior_covered(bridge, endpoint_plan)
    bridge !== nothing &&
        _bridge_covers(bridge, PUBLIC_EVALUATION_ENDPOINT_X) && return true
    return endpoint_plan !== nothing &&
        _endpoint_plan_covers(endpoint_plan, PUBLIC_EVALUATION_ENDPOINT_X)
end

"""
    match_route(coefficients, branch; controls=nothing)

Build the horizon and infinity bases for `branch` (`:IN` or `:UP`), match them
at `controls.match_x` with the Abel (Wronskian) denominator, and return the
`Route` with its incidence and reflection coefficients.
"""
function match_route(
    coefficients::CoefficientSet,
    branch;
    controls=nothing,
)
    return _match_impl(coefficients, branch, controls)
end

function _match_impl(
    coefficients::CoefficientSet,
    branch,
    controls,
)
    base_ctrl = controls === nothing ? gsn_controls(coefficients.params) : controls
    route_branch = _normalized_route_branch(branch)
    ctrl = _with_ordinary_defaults(
        _with_low_frequency_defaults(
            base_ctrl, coefficients.params, route_branch),
        coefficients.params,
        route_branch,
    )
    coefficients = near_extreme_prepare(coefficients, ctrl)
    # The small-frequency-expansion (MST) machinery is large; with small_frequency_expansion as a type parameter an ordinary
    # solve compiles only the branches it can reach.
    return _match_body(
        coefficients, route_branch, base_ctrl, ctrl, Val(ctrl.small_frequency_expansion))
end

# Core of match_route: propagates the horizon and infinity bases to match_x, solves
# the connection system, attaches the bridge, endpoint and MST evaluation plans, and
# re-enters _match_impl with ordinary controls when an MST certificate fails under
# sfe=:auto.
function _match_body(
    coefficients::CoefficientSet,
    route_branch,
    base_ctrl,
    ctrl,
    ::Val{UseSmallFrequencyExpansion},
) where {UseSmallFrequencyExpansion}
    match_x = ctrl.match_x
    selection = UseSmallFrequencyExpansion ? nothing : select_infinity_endpoint(
        coefficients,
        match_x;
        controls=ctrl,
        certificate=_infinity_abel_certificate,
    )
    scratch = IterationScratch(coefficients,
        getproperty(ctrl, :ordinary_order);
        force_eikonal=ctrl.large_frequency_expansion, large_frequency_expansion=ctrl.large_frequency_expansion, small_frequency_expansion=UseSmallFrequencyExpansion,
        route_branch=route_branch)
    mst_branches = route_branch == :IN ? (:in, :out) : (:out,)
    mst_plan = if UseSmallFrequencyExpansion
        try
            build_mano_suzuki_takasugi_plan(coefficients, match_x; branches=mst_branches)
        catch err
            _small_frequency_expansion_auto_fallback(base_ctrl, err) || rethrow()
            fallback = _ordinary_fallback(
                base_ctrl, :small_frequency_expansion_infinity_in_tail_fallback)
            return _match_impl(coefficients, route_branch, fallback)
        end
    else
        nothing
    end
    if mst_plan !== nothing &&
            mst_plan.selection.status == "NO_SAFE_SEED" &&
            base_ctrl.frequency_selection.small_frequency_expansion_request == :auto
        fallback = _ordinary_fallback(
            base_ctrl, :small_frequency_expansion_no_safe_seed_fallback)
        return _match_impl(coefficients, route_branch, fallback)
    end
    endpoint_y = mst_plan === nothing ? selection.endpoint_y :
        1.0 - mst_plan.seed_x

    if route_branch == :IN
        h_in = iterate_from_zero(coefficients, :in; controls=ctrl, match_x=match_x, scratch=scratch)
        i_in, i_out, logscaled_mst = if UseSmallFrequencyExpansion
            try
                incoming, outgoing = mano_suzuki_takasugi_infinity_pair(
                    coefficients,
                    :in,
                    :out,
                    match_x;
                    controls=ctrl,
                    scratch=scratch,
                    plan=mst_plan,
                )
                (incoming, outgoing, false)
            catch err
                if err isa ManoSuzukiTakasugiCertificateError
                    _small_frequency_expansion_auto_fallback(base_ctrl, err) || rethrow()
                    fallback = _ordinary_fallback(
                        base_ctrl, :small_frequency_expansion_infinity_in_tail_fallback)
                    return _match_impl(
                        coefficients, route_branch, fallback)
                end
                try
                    incoming, outgoing = _logscaled_mano_suzuki_takasugi_pair(
                        coefficients, match_x, ctrl, scratch, mst_plan)
                    (incoming, outgoing, true)
                catch robust_err
                    _small_frequency_expansion_auto_fallback(base_ctrl, robust_err) ||
                        throw(robust_err)
                    fallback = _ordinary_fallback(
                        base_ctrl, :small_frequency_expansion_infinity_in_tail_fallback)
                    return _match_impl(
                        coefficients, route_branch, fallback)
                end
            end
        else
            incoming, outgoing = iterate_pair_from_infinity(
                coefficients, :in, :out;
                controls=ctrl, selection, match_x, scratch)
            (incoming, outgoing, false)
        end
        fallback = logscaled_mst ? 0.0 + 0.0im :
            _basis_tuple(i_in).X * _basis_tuple(i_out).dXdx -
            _basis_tuple(i_out).X * _basis_tuple(i_in).dXdx
        denom = _match_denominator(coefficients, :infinity, match_x, ctrl, fallback)
        incidence, reflection, solution_scale = if logscaled_mst
            _connection_coefficients_logscaled(
                _logscaled_match_state(h_in),
                _logscaled_match_state(i_in),
                _logscaled_match_state(i_out),
                denom,
            )
        else
            _connection_coefficients(
                _basis_tuple(h_in),
                _basis_tuple(i_in),
                _basis_tuple(i_out),
                denom,
            )
        end
        target = _scaled_state(_basis_tuple(h_in), solution_scale)
        bridge, endpoint_plan = if UseSmallFrequencyExpansion
            _build_small_frequency_expansion_endpoint_bridge(
                coefficients,
                ctrl,
                :horizon_in,
                :IN,
                match_x,
                target,
                incidence,
                reflection,
                solution_scale,
            )
        else
            (_build_bridge(
                coefficients,
                ctrl,
                scratch,
                :IN,
                match_x,
                1.0 - endpoint_y,
                target,
                incidence,
                _basis_tuple(i_in),
                reflection,
                _basis_tuple(i_out),
            ), nothing)
        end
        mst_eval_plan = if mst_plan === nothing
            nothing
        else
            try
                logscaled_mst ? _logscaled_mano_suzuki_takasugi_eval_plan(
                    coefficients,
                    :IN,
                    mst_plan,
                    bridge,
                    endpoint_plan,
                    incidence,
                    reflection,
                    solution_scale,
                ) : build_mano_suzuki_takasugi_eval_plan(mst_plan)
            catch err
                if logscaled_mst &&
                        _public_exterior_covered(bridge, endpoint_plan)
                    nothing
                else
                    _small_frequency_expansion_auto_fallback(base_ctrl, err) || rethrow()
                    fallback = _ordinary_fallback(
                        base_ctrl, :small_frequency_expansion_evaluation_plan_fallback)
                    return _match_impl(
                        coefficients, route_branch, fallback)
                end
            end
        end
        return Route(
            :IN,
            coefficients.params,
            ctrl,
            coefficients,
            match_x,
            endpoint_y,
            selection === nothing ?
                (endpoint_plan === nothing ? 0 : endpoint_plan.selection.out_order) :
                selection.out_order,
            selection === nothing ?
                (endpoint_plan === nothing ? 0 : endpoint_plan.selection.in_order) :
                selection.in_order,
            ComplexF64(incidence),
            ComplexF64(reflection),
            solution_scale,
            h_in,
            nothing,
            i_out,
            i_in,
            bridge,
            nothing,
            endpoint_plan,
            mst_eval_plan,
        )
    end

    i_out, logscaled_mst = if UseSmallFrequencyExpansion
        try
            (_infinity_basis(
                coefficients, :out, selection, match_x, ctrl, scratch, mst_plan),
                false,
            )
        catch err
            if err isa ManoSuzukiTakasugiCertificateError
                _small_frequency_expansion_auto_fallback(base_ctrl, err) || rethrow()
                fallback = _ordinary_fallback(
                    base_ctrl, :small_frequency_expansion_infinity_in_tail_fallback)
                return _match_impl(
                    coefficients, route_branch, fallback)
            end
            try
                (_logscaled_mano_suzuki_takasugi_basis(
                    coefficients, :out, match_x, ctrl, scratch, mst_plan),
                    true,
                )
            catch robust_err
                _small_frequency_expansion_auto_fallback(base_ctrl, robust_err) ||
                    throw(robust_err)
                fallback = _ordinary_fallback(
                    base_ctrl, :small_frequency_expansion_infinity_in_tail_fallback)
                return _match_impl(
                    coefficients, route_branch, fallback)
            end
        end
    else
        (iterate_from_infinity(coefficients, :out;
            controls=ctrl, match_x, selection, scratch),
            false,
        )
    end
    h_in, h_out = iterate_pair_from_zero(
        coefficients,
        :in,
        :out;
        controls=ctrl,
        match_x=match_x,
        scratch=scratch,
    )
    fallback = _basis_tuple(h_in).X * _basis_tuple(h_out).dXdx -
        _basis_tuple(h_out).X * _basis_tuple(h_in).dXdx
    denom = _match_denominator(
        coefficients,
        :horizon,
        match_x,
        ctrl,
        fallback;
        fallback_scale=_basis_scale(_basis_tuple(h_in), _basis_tuple(h_out)),
    )
    reflection, incidence, solution_scale = if logscaled_mst
        _connection_coefficients_logscaled(
            _logscaled_match_state(i_out),
            _logscaled_match_state(h_in),
            _logscaled_match_state(h_out),
            denom,
        )
    else
        _connection_coefficients(
            _basis_tuple(i_out),
            _basis_tuple(h_in),
            _basis_tuple(h_out),
            denom,
        )
    end
    target = logscaled_mst ? materialize_logscaled_state(
        _logscaled_match_state(i_out); scale=solution_scale) :
        _scaled_state(_basis_tuple(i_out), solution_scale)
    bridge = _build_bridge(
        coefficients,
        ctrl,
        scratch,
        :UP,
        match_x,
        min(h_in.seed_x, h_out.seed_x),
        target,
        reflection,
        _basis_tuple(h_in),
        incidence,
        _basis_tuple(h_out),
    )
    endpoint_basis, boundary_plan = UseSmallFrequencyExpansion && !logscaled_mst ?
        _small_frequency_expansion_up_basis(coefficients, ctrl, scratch, match_x, target) : (nothing,nothing)
    endpoint_bridge, endpoint_plan = if !UseSmallFrequencyExpansion
        (nothing, nothing)
    elseif endpoint_basis === nothing && boundary_plan === nothing
        _build_small_frequency_expansion_endpoint_bridge(
            coefficients,
            ctrl,
            :infinity_out,
            :UP,
            match_x,
            target,
            incidence,
            reflection,
            solution_scale,
        )
    else
        (endpoint_basis, boundary_plan)
    end
    mst_eval_plan = if mst_plan === nothing
        nothing
    else
        try
            logscaled_mst ? _logscaled_mano_suzuki_takasugi_eval_plan(
                coefficients,
                :UP,
                mst_plan,
                endpoint_bridge,
                endpoint_plan,
                incidence,
                reflection,
                solution_scale,
                i_out,
            ) : build_mano_suzuki_takasugi_eval_plan(mst_plan)
        catch err
            if logscaled_mst
                _public_exterior_covered(endpoint_bridge, endpoint_plan) ?
                    nothing : rethrow()
            else
                _small_frequency_expansion_auto_fallback(base_ctrl, err) || rethrow()
                fallback = _ordinary_fallback(
                    base_ctrl, :small_frequency_expansion_evaluation_plan_fallback)
                return _match_impl(coefficients, route_branch, fallback)
            end
        end
    end
    return Route(
        :UP,
        coefficients.params,
        ctrl,
        coefficients,
        match_x,
        endpoint_y,
        selection === nothing ?
            (endpoint_plan === nothing ? 0 : endpoint_plan.selection.out_order) :
            selection.out_order,
        selection === nothing ?
            (endpoint_plan === nothing ? 0 : endpoint_plan.selection.in_order) :
            selection.in_order,
        ComplexF64(incidence),
        ComplexF64(reflection),
        solution_scale,
        h_in,
        h_out,
        i_out,
        nothing,
        bridge,
        endpoint_bridge,
        endpoint_plan,
        mst_eval_plan,
    )
end

"""
    radial_route(s, l, m, a, omega, branch=:IN; lambda=nothing, nu=nothing, controls=nothing, kwargs...)
    radial_route(params::GSNParameters, branch=:IN; controls=nothing, coefficients=nothing, kwargs...)

Build the ISEM GSN route for the `:IN` or `:UP` solution. Extra keywords are
passed to `gsn_controls` when `controls` is not given. A real negative
frequency is solved at (-m, -omega) and returned as a `ConjugatedRoute`.
"""
function radial_route(
    s::Integer,
    l::Integer,
    m::Integer,
    a,
    omega,
    branch=:IN;
    lambda=nothing,
    nu=nothing,
    controls=nothing,
    kwargs...,
)
    params = gsn_parameters(s, l, m, a, omega; lambda=lambda, nu=nu)
    ctrl = controls === nothing ? gsn_controls(params; kwargs...) : controls
    if isreal(params.omega) && params.omega < 0
        positive_params = GSNParameters(
            params.s,
            params.l,
            -params.m,
            params.a,
            -params.omega,
            params.lambda,
            params.nu,
            params.kappa,
        )
        # The reflected route is the positive-frequency route itself, so the
        # two signs agree by construction.
        positive_route = radial_route(
            positive_params,
            branch;
            controls=ctrl,
        )
        return ConjugatedRoute(positive_route, params)
    end
    coefficients = gsn_coefficients(params; controls=ctrl)
    return match_route(coefficients, branch; controls=ctrl)
end

function radial_route(
    params::GSNParameters,
    branch=:IN;
    controls=nothing,
    coefficients=nothing,
    kwargs...,
)
    ctrl = controls === nothing ? gsn_controls(params; kwargs...) : controls
    coeffs = coefficients === nothing ? gsn_coefficients(params; controls=ctrl) : coefficients
    return match_route(coeffs, branch; controls=ctrl)
end

function _bridge_covers(basis::Basis, x)
    if basis.endpoint_valid
        horizon = basis.kind in (:horizon_in, :horizon_out)
        ((horizon && x <= basis.seed_x) || (!horizon && x >= basis.seed_x)) &&
            return true
    end
    isempty(basis.patches) && return x == basis.seed_x
    target_x = basis.patches[end].next_x
    lo = min(basis.seed_x, target_x) - 20eps(Float64)
    hi = max(basis.seed_x, target_x) + 20eps(Float64)
    return lo <= x <= hi
end

_bridge_covers(basis::ScaledBasis, x) =
    _bridge_covers(basis.basis, x)

function _bridge_covers(basis::LogScaledBasis, x)
    isempty(basis.patches) && return x == basis.seed_x
    target_x = basis.patches[end].next_x
    lo = min(basis.seed_x, target_x) - 20eps(Float64)
    hi = max(basis.seed_x, target_x) + 20eps(Float64)
    return lo <= x <= hi
end

_bridge_value(basis::Basis, x) = basis_value(basis, x)
_bridge_value(basis::ScaledBasis, x) =
    basis.scale * basis_value(basis.basis, x)

_bridge_state(basis::Basis, x) = basis_state(basis, x)
function _bridge_state(basis::ScaledBasis, x)
    value, derivative = basis_state(basis.basis, x)
    return basis.scale * value, basis.scale * derivative
end

function _physical_logscaled_basis_state(basis, x, scale)
    return materialize_logscaled_state(
        logscaled_basis_state(basis, x); scale=scale)
end

_route_basis(basis::Basis) = basis
_route_basis(basis::ScaledBasis) = basis.basis
_route_basis(basis::LogScaledBasis) = basis

_endpoint_plan_covers(plan::EndpointPlan, x) =
    x >= plan.seed_x - 20eps(Float64)

function _endpoint_plan_state(route::Route, x; radius=nothing)
    plan = route.endpoint_plan
    plan === nothing && error("missing ISEM GSN endpoint evaluation plan")
    return _physical_endpoint_state(
        route.coefficients,
        route.branch,
        plan.selection,
        x,
        route.incidence,
        route.reflection,
        route.solution_scale,
        plan.scale;
        radius,
    )
end

@inline _logscaled_mano_suzuki_takasugi_eval(plan) = plan !== nothing &&
    hasproperty(plan, :logscaled) && plan.logscaled

function _physical_mano_suzuki_takasugi_out_state(route::Route, x; radius=nothing)
    if _logscaled_mano_suzuki_takasugi_eval(route.mst_plan)
        route.mst_plan.out_scale === nothing &&
            error("missing ISEM GSN log-scaled outgoing normalization")
        return mano_suzuki_takasugi_scaled_state(
            route.coefficients,
            route.mst_plan,
            :out,
            x;
            scale=route.mst_plan.out_scale,
            radius,
        )
    end
    state = mano_suzuki_takasugi_route_state(
        route.coefficients, route.mst_plan, :out, x; radius)
    return (
        X=route.solution_scale * state.X,
        dXdx=route.solution_scale * state.dXdx,
        estimated_relerr=state.estimated_relerr,
    )
end

function _mano_suzuki_takasugi_combined_in_state(route::Route, x; radius=nothing)
    if x <= route.mst_plan.seed_x
        return mano_suzuki_takasugi_pin_state(route.coefficients, route.mst_plan, x;
            scale=route.solution_scale * _pin_unit_scale(route.params), radius)
    end
    incoming = mano_suzuki_takasugi_route_state(route.coefficients, route.mst_plan, :in, x; radius)
    outgoing = mano_suzuki_takasugi_route_state(route.coefficients, route.mst_plan, :out, x; radius)
    xi, xo = route.incidence * incoming.X, route.reflection * outgoing.X
    di, do_ = route.incidence * incoming.dXdx, route.reflection * outgoing.dXdx
    value, derivative = xi + xo, di + do_
    # A near-zone IN solution must not be recovered by subtracting two
    # much larger infinity solutions when that subtraction loses precision.
    if eps(Float64) * (abs(xi) + abs(xo)) > route.controls.tolerance * abs(value) ||
            eps(Float64) * (abs(di) + abs(do_)) > route.controls.tolerance * abs(derivative)
        return mano_suzuki_takasugi_pin_state(route.coefficients, route.mst_plan, x;
            scale=route.solution_scale * _pin_unit_scale(route.params), radius)
    end
    return (X=value, dXdx=derivative)
end

# Whether a point of a route with an MST plan is evaluated from the MST
# representation instead of the propagated bases, bridges or endpoint plan.
function _mano_suzuki_takasugi_evaluation_region(route::Route, x)
    route.mst_plan === nothing && return false
    if route.branch == :UP
        if x > route.match_x
            # Outward propagation of the subdominant UP solution amplifies
            # its growing companion. An inward asymptotic basis is stable.
            route.endpoint_bridge isa ScaledBasis &&
                _bridge_covers(route.endpoint_bridge, x) && return false
            return !(route.endpoint_plan !== nothing && _endpoint_plan_covers(route.endpoint_plan, x))
        end
        x <= min(route.horizon_in.seed_x, route.horizon_out.seed_x) && return false
        if _bridge_covers(route.horizon_in, x) && _bridge_covers(route.horizon_out, x)
            vin, din = basis_state(route.horizon_in, x)
            vout, dout = basis_state(route.horizon_out, x)
            xi, xo = route.reflection * vin, route.incidence * vout
            di, do_ = route.reflection * din, route.incidence * dout
            eps(Float64) * (abs(xi) + abs(xo)) <=
                route.controls.tolerance * abs(xi + xo) &&
                eps(Float64) * (abs(di) + abs(do_)) <=
                route.controls.tolerance * abs(di + do_) && return false
        end
        return !(route.bridge !== nothing && _bridge_covers(route.bridge, x))
    end
    return x > route.match_x &&
        !(route.bridge !== nothing && _bridge_covers(route.bridge, x)) &&
        !(route.endpoint_plan !== nothing && _endpoint_plan_covers(route.endpoint_plan, x))
end

"""
    evaluate_route(route, x)

Value of the route solution at the compactified coordinate `x`, taken from the
basis, bridge, endpoint plan or MST evaluation that covers `x`.
"""
function evaluate_route(route::Route, x)
    x = Float64(x)
    if route.branch == :IN
        if x <= route.match_x
            return route.solution_scale *
                basis_value(_required_basis(route.horizon_in, "horizon_in"), x)
        end
        route.bridge !== nothing && _bridge_covers(route.bridge, x) &&
            return _bridge_value(route.bridge, x)
        if route.endpoint_plan !== nothing &&
                _endpoint_plan_covers(route.endpoint_plan, x)
            return _endpoint_plan_state(route, x).X
        end
        if route.mst_plan !== nothing
            if _logscaled_mano_suzuki_takasugi_eval(route.mst_plan)
                return mano_suzuki_takasugi_pin_state(
                    route.coefficients, route.mst_plan, x).X
            end
            return _mano_suzuki_takasugi_combined_in_state(route, x).X
        end
        return route.incidence * basis_value(_required_basis(route.infinity_in, "infinity_in"), x) +
            route.reflection * basis_value(_required_basis(route.infinity_out, "infinity_out"), x)
    end
    if _mano_suzuki_takasugi_evaluation_region(route, x)
        return _physical_mano_suzuki_takasugi_out_state(route, x).X
    end
    if x <= route.match_x
        route.bridge !== nothing && _bridge_covers(route.bridge, x) &&
            return _bridge_value(route.bridge, x)
        return route.reflection * basis_value(_required_basis(route.horizon_in, "horizon_in"), x) +
            route.incidence * basis_value(_required_basis(route.horizon_out, "horizon_out"), x)
    end
    if route.endpoint_plan !== nothing &&
            _endpoint_plan_covers(route.endpoint_plan, x)
        return _endpoint_plan_state(route, x).X
    end
    route.endpoint_bridge !== nothing &&
        _bridge_covers(route.endpoint_bridge, x) &&
        return _bridge_value(route.endpoint_bridge, x)
    if route.infinity_out isa LogScaledBasis &&
            _bridge_covers(route.infinity_out, x)
        return _physical_logscaled_basis_state(
            route.infinity_out, x, route.solution_scale).X
    end
    return route.solution_scale *
        basis_value(_required_basis(route.infinity_out, "infinity_out"), x)
end

evaluate_route(route::ConjugatedRoute, x) =
    conj(evaluate_route(route.route, x))

"""
    route_state(route, x; radius=nothing)

Value and x-derivative `(X, dXdx)` of the route solution at the compactified
coordinate `x`. Near infinity, `radius` (if given) is used to form the
endpoint coordinate directly, avoiding cancellation in 1 - x.
"""
function route_state(route::Route, x; radius=nothing)
    x = Float64(x)
    if route.branch == :IN
        if x <= route.match_x
            value, derivative = basis_state(_required_basis(route.horizon_in, "horizon_in"), x)
            return (
                X=route.solution_scale * value,
                dXdx=route.solution_scale * derivative,
            )
        end
        if route.bridge !== nothing && _bridge_covers(route.bridge, x)
            value, derivative = _bridge_state(route.bridge, x)
            return (X=value, dXdx=derivative)
        end
        if route.endpoint_plan !== nothing &&
                _endpoint_plan_covers(route.endpoint_plan, x)
            return _endpoint_plan_state(route, x; radius)
        end
        if route.mst_plan !== nothing
            if _logscaled_mano_suzuki_takasugi_eval(route.mst_plan)
                state = mano_suzuki_takasugi_pin_state(
                    route.coefficients, route.mst_plan, x; radius)
                return (X=state.X, dXdx=state.dXdx)
            end
            return _mano_suzuki_takasugi_combined_in_state(route, x; radius)
        end
        vin, din = basis_state(_required_basis(route.infinity_in, "infinity_in"), x)
        vout, dout = basis_state(_required_basis(route.infinity_out, "infinity_out"), x)
        return (
            X=route.incidence * vin + route.reflection * vout,
            dXdx=route.incidence * din + route.reflection * dout,
        )
    end
    if _mano_suzuki_takasugi_evaluation_region(route, x)
        return _physical_mano_suzuki_takasugi_out_state(route, x; radius)
    end
    if x <= route.match_x
        if route.bridge !== nothing && _bridge_covers(route.bridge, x)
            value, derivative = _bridge_state(route.bridge, x)
            return (X=value, dXdx=derivative)
        end
        vin, din = basis_state(_required_basis(route.horizon_in, "horizon_in"), x)
        vout, dout = basis_state(_required_basis(route.horizon_out, "horizon_out"), x)
        return (
            X=route.reflection * vin + route.incidence * vout,
            dXdx=route.reflection * din + route.incidence * dout,
        )
    end
    if route.endpoint_plan !== nothing &&
            _endpoint_plan_covers(route.endpoint_plan, x)
        return _endpoint_plan_state(route, x; radius)
    end
    if route.endpoint_bridge !== nothing &&
            _bridge_covers(route.endpoint_bridge, x)
        basis = route.endpoint_bridge
        value, derivative = if radius !== nothing &&
                basis isa ScaledBasis{<:Basis{Patch}}
            kappa = route.params.kappa
            y = 2kappa / (radius - 1 + kappa)
            basis_state(basis.basis, x, y) .* basis.scale
        else
            _bridge_state(basis, x)
        end
        return (X=value, dXdx=derivative)
    end
    if route.infinity_out isa LogScaledBasis &&
            _bridge_covers(route.infinity_out, x)
        return _physical_logscaled_basis_state(
            route.infinity_out, x, route.solution_scale)
    end
    value, derivative = basis_state(_required_basis(route.infinity_out, "infinity_out"), x)
    return (
        X=route.solution_scale * value,
        dXdx=route.solution_scale * derivative,
    )
end

function route_state(route::ConjugatedRoute, x)
    state = route_state(route.route, x)
    return (X=conj(state.X), dXdx=conj(state.dXdx))
end

(route::Route)(x) = evaluate_route(route, x)
(route::ConjugatedRoute)(x) = evaluate_route(route, x)

"""
    route_patch_count(route)

Total number of propagation patches in all bases of a route.
"""
function route_patch_count(route::Route)
    count = 0
    for basis in (route.horizon_in, route.horizon_out, route.infinity_out,
            route.infinity_in, route.bridge, route.endpoint_bridge)
        count += basis === nothing ? 0 :
            basis_patch_count(_route_basis(basis))
    end
    return count
end

route_patch_count(route::ConjugatedRoute) =
    route_patch_count(route.route)

function _patch_counts(route::Route)
    ordinary = 0
    eikonal = 0
    endpoint = 0
    horizon = 0
    infinity = 0
    for (location, basis) in ((:horizon, route.horizon_in),
            (:horizon, route.horizon_out),
            (:infinity, route.infinity_out),
            (:infinity, route.infinity_in),
            (route.branch == :IN ? :infinity : :horizon, route.bridge),
            (:infinity, route.endpoint_bridge))
        basis === nothing && continue
        basis = _route_basis(basis)
        count = length(basis.patches)
        location == :horizon ? (horizon += count) : (infinity += count)
        for patch in basis.patches
            if patch isa Patch || patch isa ScaledPatch ||
                    patch isa LogScaledPatch
                ordinary += 1
            elseif patch isa EikonalPatch
                eikonal += 1
            else
                endpoint += 1
            end
        end
    end
    return ordinary, eikonal, endpoint, horizon, infinity
end

"""
    route_plan_of(route)

Summary (`RoutePlan`) of how a route was built: frequency and spin
regimes, endpoint and propagation methods, orders, and patch counts by kind.
"""
function route_plan_of(route::Route)
    controls = route.controls
    frequency = controls.frequency_selection
    spin = select_spin_regime(route.params)
    eikonal = eikonal_preflight(
        route.params, frequency.regime, route.branch)
    ordinary_patches, eikonal_patches, endpoint_patches,
        horizon_patches, infinity_patches = _patch_counts(route)
    logscaled_mst = route.infinity_out isa LogScaledBasis ||
        route.infinity_in isa LogScaledBasis
    endpoint = route.endpoint_bridge isa ScaledBasis ?
        :direct_endpoint_basis :
        logscaled_mst ? :mst_logscaled :
        route.endpoint_plan !== nothing ? :mst_scaled_y_handoff :
        controls.small_frequency_expansion ? :mst : controls.large_frequency_expansion ? :large_frequency_expansion_carrier :
        near_extreme_selected(route.coefficients) ? :near_extremal : :taylor
    propagation = eikonal_patches > 0 ? :eikonal : :taylor
    return RoutePlan(
        frequency.regime,
        frequency.reason,
        spin.regime,
        spin.reason,
        endpoint,
        propagation,
        :abel,
        route.match_x,
        controls.horizon_order,
        controls.ordinary_order,
        controls.infinity_order,
        eikonal.candidate,
        ordinary_patches,
        eikonal_patches,
        endpoint_patches,
        horizon_patches,
        infinity_patches,
        spin.angular_compression,
        spin.phase_compression,
    )
end

route_plan_of(route::ConjugatedRoute) =
    route_plan_of(route.route)

# Relative errors (incidence, reflection) of a real-frequency route. Each basis
# is propagated with local errors below the route tolerance per step from a
# seed certified to max(tolerance, the endpoint seed target); a local relative
# error e perturbs both terms of the matched state by about e times the larger
# one, so the smaller term c_j g_j loses e |larger| / |c_j g_j|. On the real
# axis both solutions keep their size along the path, so the term ratio at the
# match holds along it. The rounding of the 2x2 solve adds eps times the
# cancellation of each coefficient's numerator.
function amplitude_errors(route::Route)
    target, first, second = route.branch == :IN ?
        (route.horizon_in, route.infinity_in, route.infinity_out) :
        (route.infinity_out, route.horizon_in, route.horizon_out)
    c1, c2 = route.branch == :IN ? (route.incidence, route.reflection) :
        (route.reflection, route.incidence)
    tolerance = route.controls.tolerance
    seed = max(tolerance, _ENDPOINT_SEED_TARGET)
    local_error = sum(basis -> basis.step_count * tolerance + seed,
        (target, first, second))
    t = _logscaled_match_state(target)
    b1 = _logscaled_match_state(first)
    b2 = _logscaled_match_state(second)
    log1 = _logabs(c1) + _logabs(b1.X) + b1.log_scale
    log2 = _logabs(c2) + _logabs(b2.X) + b2.log_scale
    larger = max(log1, log2)
    # Numerators t.X b2' - b2.X t' and b1.X t' - t.X b1': cancellation from the
    # ratio of their two products (scales cancel).
    cancellation(ratio) = (1 + abs(ratio)) / abs(1 - ratio)
    cancel1 = cancellation((b2.X * t.dXdx) / (t.X * b2.dXdx))
    cancel2 = cancellation((t.X * b1.dXdx) / (b1.X * t.dXdx))
    e1 = local_error * exp(larger - log1) + eps(Float64) * cancel1
    e2 = local_error * exp(larger - log2) + eps(Float64) * cancel2
    finite(value) = isfinite(value) ? Float64(value) : Inf
    return route.branch == :IN ? (finite(e1), finite(e2)) :
        (finite(e2), finite(e1))
end

amplitude_errors(route::ConjugatedRoute) =
    amplitude_errors(route.route)

"""
    route_truncations(route)

All series truncations (`Truncation`) recorded by the bases of a route.
"""
function route_truncations(route::Route)
    records = Truncation[]
    for basis in (route.horizon_in, route.horizon_out, route.infinity_out,
            route.infinity_in, route.bridge, route.endpoint_bridge)
        basis === nothing || append!(records, _route_basis(basis).truncations)
    end
    return records
end

route_truncations(route::ConjugatedRoute) =
    route_truncations(route.route)

end
