module NIAAmplitudes

import ..MSTInfinity
using ..MSTInfinity:
    mst_anchor,
    mst_nia_amplitudes,
    mst_nu_complex,
    mst_principal_amplitudes
using ..CoefficientTables: gsn_coefficients
using ..Parameters:
    gsn_controls,
    gsn_parameters
using ..Iteration:
    iterate_from_zero
import ..ComplexFrequency
import ..ComplexRational
using ..ComplexFrequency:
    ComplexRoute,
    ComplexRoutePlan,
    ComplexSFEvaluator,
    ComplexStateEvaluator
using ..ComplexRational:
    ComplexRationalEvaluator,
    ComplexRationalState

const _REAL_RATIO_MAX = 1.0e-8
const _NONROTATING_MAX_ABS_A = 64eps(Float64)
const _MIN_ABS_IMAG_OMEGA = 1.8
const _MAX_ABS_IMAG_OMEGA = 2.2
const _NIA_LATERAL_REAL_RATIO_MAX = 1.0e-6
const _STATE_ANCHOR_X = 0.5
const _STATE_ANCHOR_MAX_ERROR = 1.0e-10
const _REAL_AXIS_ANCHOR_MIN_DEPARTURE = 1.0e-8
const _STATE_ANCHOR_SUM_TOLERANCE = 1.0e-11
const _STATE_ANCHOR_RELAXED_TOLERANCE = 1.0e-9
const _STATE_ANCHOR_RELAXED_REAL_RATIO_MAX = 1.0e-3
const _AMPLITUDE_AGREEMENT_MAX = 1.0e-10
const _MST_STATE_NU_DELTA_MIN = 1.0e-10
const _FAST_ANCHOR_MIN_OMEGA = 0.025
const _FAST_ANCHOR_MAX_OMEGA = 2.0
const _FAST_ANCHOR_RATIO_LOW = 0.2
const _FAST_ANCHOR_RATIO_HIGH = 5.0
const _FAST_ANCHOR_ORDER = 40
const _FAST_ANCHOR_TOLERANCE = 1.0e-14

# State evaluator that multiplies X and dX/dr* of another evaluator by a
# constant complex scale.
struct ScaledEvaluator{E}
    evaluator::E
    scale::ComplexF64
end

@inline _finite(value) =
    isfinite(real(value)) && isfinite(imag(value))

function ComplexFrequency._state_at_r(
    evaluator::ScaledEvaluator,
    r::Real,
)
    state = ComplexFrequency._state_at_r(evaluator.evaluator, r)
    return (
        X=ComplexF64(evaluator.scale * state.X),
        dXdrstar=ComplexF64(evaluator.scale * state.dXdrstar),
        error=state.error,
    )
end

function (evaluator::ScaledEvaluator)(rstar::Real)
    X, dXdrstar, error = evaluator.evaluator(rstar)
    return (
        ComplexF64(evaluator.scale * X),
        ComplexF64(evaluator.scale * dXdrstar),
        error,
    )
end

function (evaluator::ScaledEvaluator)(rstars::AbstractVector{<:Real})
    states = map(evaluator, rstars)
    return (
        getindex.(states, 1),
        getindex.(states, 2),
        getindex.(states, 3),
    )
end

@inline _rational_core(evaluator::ComplexRationalEvaluator) =
    (evaluator, ComplexF64(1))

function _rational_core(evaluator::ScaledEvaluator)
    core = _rational_core(evaluator.evaluator)
    core === nothing && return nothing
    return (core[1], ComplexF64(evaluator.scale * core[2]))
end

@inline _rational_core(evaluator) = nothing

@inline function _pair_error(first, second)
    scale = max(
        abs(first[1]), abs(first[2]),
        abs(second[1]), abs(second[2]),
        floatmin(Float64),
    )
    return max(
        abs(first[1] - second[1]),
        abs(first[2] - second[2]),
    ) / scale
end

# Least-squares complex scale c minimizing |target - c * candidate| over the
# (value, derivative) pair; nothing if either pair is zero or not finite.
function _fit_pair_scale(target, candidate)
    target_scale = max(abs(target[1]), abs(target[2]))
    candidate_scale = max(abs(candidate[1]), abs(candidate[2]))
    target_scale > 0 && candidate_scale > 0 || return nothing
    isfinite(target_scale) && isfinite(candidate_scale) || return nothing
    target_value = target[1] / target_scale
    target_derivative = target[2] / target_scale
    candidate_value = candidate[1] / candidate_scale
    candidate_derivative = candidate[2] / candidate_scale
    denominator =
        abs2(candidate_value) + abs2(candidate_derivative)
    denominator > 0 && isfinite(denominator) || return nothing
    unit_scale = (
        conj(candidate_value) * target_value +
        conj(candidate_derivative) * target_derivative
    ) / denominator
    scale = (target_scale / candidate_scale) * unit_scale
    return _finite(scale) ? ComplexF64(scale) : nothing
end

# For an IN route with a rational evaluator and omega in the near-real sector,
# compare the state at the match point with a real-axis propagation from the
# horizon and rescale the evaluator when they differ by a pure scale.
# Returns (evaluator, metadata).
function _in_real_axis_anchor(route, evaluator)
    rejected(status; projective_error=Inf, scale=ComplexF64(NaN, NaN),
            departure=Inf, steps=-1) = (
        evaluator,
        metadata=(;
            real_axis_state_anchor_status=status,
            real_axis_state_anchor_projective_error=
                Float64(projective_error),
            real_axis_state_anchor_scale=ComplexF64(scale),
            real_axis_state_anchor_departure=Float64(departure),
            real_axis_state_anchor_steps=Int(steps),
        ),
    )
    route.branch == :IN || return rejected(:not_applicable)
    core = _rational_core(evaluator)
    core === nothing && return rejected(:not_rational)
    rational, existing_scale = core
    omega = ComplexF64(rational.params.omega)
    abs(real(omega)) >=
        _FAST_ANCHOR_RATIO_LOW * abs(imag(omega)) ||
        return rejected(:outside_real_axis_sector)
    controls = gsn_controls(
        rational.params;
        N=rational.settings.order,
        xm=rational.match_x,
        tol=route.plan.controls.tolerance,
        small_frequency_expansion=false,
        large_frequency_expansion=false,
    )
    basis = try
        iterate_from_zero(
            rational.coefficients,
            :in;
            controls,
            match_x=rational.match_x,
        )
    catch
        return rejected(:construction_error)
    end
    basis.endpoint_valid ||
        return rejected(:invalid_endpoint; steps=basis.step_count)
    target = (basis.match_X, basis.match_dXdx)
    current = (
        existing_scale * rational.match_state.X,
        existing_scale * rational.match_state.dXdx,
    )
    all(_finite, (target..., current...)) ||
        return rejected(:nonfinite_state; steps=basis.step_count)
    correction = _fit_pair_scale(target, current)
    correction === nothing &&
        return rejected(:rejected_scale; steps=basis.step_count)
    corrected = (
        correction * current[1],
        correction * current[2],
    )
    projective_error = Float64(_pair_error(corrected, target))
    departure = Float64(abs(correction - 1))
    projective_error <= _STATE_ANCHOR_MAX_ERROR || return rejected(
        :rejected_projective;
        projective_error,
        scale=correction,
        departure,
        steps=basis.step_count,
    )
    departure > _REAL_AXIS_ANCHOR_MIN_DEPARTURE || return rejected(
        :retained;
        projective_error,
        scale=correction,
        departure,
        steps=basis.step_count,
    )
    # A projective match cannot certify the normalization of a dominant state.
    # Require the reference propagation to converge before using its scale.
    confirmation_controls = gsn_controls(
        rational.params;
        N=rational.settings.order,
        horizon_order=ComplexRational._horizon_series_order(rational.settings.order),
        xm=rational.match_x,
        tol=route.plan.controls.tolerance,
        small_frequency_expansion=false,
        large_frequency_expansion=false,
    )
    confirmation = iterate_from_zero(rational.coefficients, :in;
        controls=confirmation_controls, match_x=rational.match_x)
    reference_error = _pair_error(target,
        (confirmation.match_X,confirmation.match_dXdx))
    reference_error <= _STATE_ANCHOR_MAX_ERROR || return rejected(
        :uncertified_reference;
        projective_error,
        scale=correction,
        departure,
        steps=basis.step_count,
    )
    return (
        evaluator=ScaledEvaluator(evaluator, correction),
        metadata=(;
            real_axis_state_anchor_status=:rescaled,
            real_axis_state_anchor_projective_error=projective_error,
            real_axis_state_anchor_scale=correction,
            real_axis_state_anchor_departure=departure,
            real_axis_state_anchor_steps=basis.step_count,
        ),
    )
end

# Whether `build` overlays analytic MST amplitudes: default controls, |s| = 2,
# a != 0, and omega within _REAL_RATIO_MAX of the negative imaginary axis with
# _MIN_ABS_IMAG_OMEGA <= |Im(omega)| <= _MAX_ABS_IMAG_OMEGA.
@inline function _pretrigger(
    s,
    a,
    omega;
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
    backend == :direct_rational || return false
    controls === nothing || return false
    # A tolerance alone sets accuracy targets and does not fix the route.
    all(value === nothing for value in (
        lambda,
        nu,
        xm,
        rhom,
        N,
        small_frequency_expansion,
        large_frequency_expansion,
    )) || return false
    pole_normalization && return false
    omegac = ComplexF64(omega)
    abs(s) == 2 || return false
    abs(a) > _NONROTATING_MAX_ABS_A || return false
    imag(omegac) < 0 || return false
    !iszero(real(omegac)) || return false
    abs(real(omegac)) <= _REAL_RATIO_MAX * abs(imag(omegac)) ||
        return false
    return _MIN_ABS_IMAG_OMEGA <= abs(imag(omegac)) <=
        _MAX_ABS_IMAG_OMEGA
end

# Replace the route's incidence and reflection amplitudes (GSN and Teukolsky)
# by transmission times the analytic ratios in `result`, and record the
# analytic certificate in the metadata. Returns the route unchanged if any
# value is not finite.
function _overlay(
    route::ComplexRoute,
    result;
    amplitude_backend::Symbol=:nia_mst_analytic,
    state_evaluator=route.state_evaluator,
    state_metadata=(;),
)
    values = (
        route.transmission,
        route.teukolsky_transmission,
        result.gsn...,
        result.teuk...,
    )
    all(_finite, values) || return route
    iszero(route.transmission) && return route
    iszero(route.teukolsky_transmission) && return route

    incidence = ComplexF64(route.transmission * result.gsn[1])
    reflection = ComplexF64(route.transmission * result.gsn[2])
    teukolsky_incidence =
        ComplexF64(route.teukolsky_transmission * result.teuk[1])
    teukolsky_reflection =
        ComplexF64(route.teukolsky_transmission * result.teuk[2])
    all(_finite, (
        incidence,
        reflection,
        teukolsky_incidence,
        teukolsky_reflection,
    )) || return route
    # The real-axis state is not rebuilt from the reported amplitudes: a
    # change of the state at the match point that is within the amplitudes'
    # error grows along the real axis (weak-damping UP solutions), so the
    # state keeps its own construction and certificate.

    metadata = merge(route.metadata, (;
        amplitude_backend,
        mst_amplitude_certificate=result.certificate_kind,
        mst_amplitude_representation_spread=
            result.representation_spread,
        mst_amplitude_nearest_agreement=result.nearest_agreement,
        mst_amplitude_truncation_agreement=
            result.truncation_agreement,
        mst_amplitude_max_condition=result.max_condition,
        mst_amplitude_errors=hasproperty(result, :gsn_error) ?
            Float64.(result.gsn_error) : (NaN, NaN),
        mst_amplitude_component_agreement=
            hasproperty(result, :component_agreement) ?
            Float64.(result.component_agreement) : (NaN, NaN),
        mst_amplitude_coefficient_residual=
            hasproperty(result, :coefficient_residual) ?
            Float64(result.coefficient_residual) : NaN,
        mst_amplitude_coefficient_source_shift=
            hasproperty(result, :coefficient_source_shift) ?
            Int(result.coefficient_source_shift) : 0,
    ), state_metadata)
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
        state_evaluator,
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

function _state_overlay(route::ComplexRoute, anchor)
    anchor.evaluator === route.state_evaluator &&
        isempty(anchor.metadata) && return route
    return ComplexRoute(
        route.branch,
        route.params,
        route.controls,
        route.p_solution,
        anchor.evaluator,
        merge(route.metadata, anchor.metadata),
        route.plan,
        route.transmission,
        route.incidence,
        route.reflection,
        route.teukolsky_transmission,
        route.teukolsky_incidence,
        route.teukolsky_reflection,
    )
end

function _with_amplitude_errors(route, errors)
    metadata = merge(route.metadata, (;
        incidence_error=Float64(errors[1]),
        reflection_error=Float64(errors[2])))
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

function _retain_analytic_agreement(route, result, agreement)
    metadata = merge(route.metadata, (;
        amplitude_backend=agreement <= _AMPLITUDE_AGREEMENT_MAX ?
            :complex_rational_certified : :complex_rational,
        mst_amplitude_certificate=result.certificate_kind,
        mst_amplitude_agreement=Float64(agreement),
        mst_amplitude_representation_spread=
            result.representation_spread,
        mst_amplitude_nearest_agreement=result.nearest_agreement,
        mst_amplitude_truncation_agreement=
            result.truncation_agreement,
        mst_amplitude_max_condition=result.max_condition,
        mst_amplitude_errors=hasproperty(result, :gsn_error) ?
            Float64.(result.gsn_error) : (NaN, NaN),
        mst_amplitude_component_agreement=
            hasproperty(result, :component_agreement) ?
            Float64.(result.component_agreement) : (NaN, NaN),
        mst_amplitude_coefficient_residual=
            hasproperty(result, :coefficient_residual) ?
            Float64(result.coefficient_residual) : NaN,
        mst_amplitude_coefficient_source_shift=
            hasproperty(result, :coefficient_source_shift) ?
            Int(result.coefficient_source_shift) : 0,
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

@inline function _nu_delta(nu)
    value = real(ComplexF64(nu))
    return Float64(value - round(value))
end

@inline function _allows_uncertified_anchor(omega)
    value = ComplexF64(omega)
    denominator = abs(imag(value))
    denominator > 0 && isfinite(denominator) || return false
    ratio = abs(real(value)) / denominator
    return isfinite(ratio) && ratio < _FAST_ANCHOR_RATIO_HIGH
end

@inline function _anchor_accepted(anchor, omega, require_certificate)
    anchor.certified && return true
    require_certificate && return false
    _allows_uncertified_anchor(omega) && return true
    return isfinite(anchor.condition) &&
        anchor.condition <= MSTInfinity.MST_COEFF_CONDITION_MAX
end

@inline function _allows_relaxed_state_anchor(omega)
    value = ComplexF64(omega)
    imag(value) < 0 || return false
    !iszero(real(value)) || return false
    return abs(real(value)) <=
        _STATE_ANCHOR_RELAXED_REAL_RATIO_MAX *
        max(abs(imag(value)), 1.0)
end

# Outgoing MST state at x = _STATE_ANCHOR_X. Uses the monodromy state at the
# complex nu when nu is away from an integer and the sum converges, else the
# principal-branch anchor; `certified` and `source` record which was used.
function _mst_out_anchor(
    source;
    nu=nothing,
    allow_relaxed::Bool=false,
)
    params = hasproperty(source, :params) ? source.params : source
    selected_nu = if nu === nothing
        try
            mst_nu_complex(
                params.s,
                params.l,
                params.m,
                params.a,
                params.omega,
                params.lambda,
            )
        catch
            nothing
        end
    else
        ComplexF64(nu)
    end
    if selected_nu !== nothing
        selected_nu = ComplexF64(selected_nu)
        nu_delta = _nu_delta(selected_nu)
        if nu_delta > _MST_STATE_NU_DELTA_MIN
            selected = try
                MSTInfinity.mst_monodromy_state(
                    params,
                    _STATE_ANCHOR_X;
                    tolerance=_STATE_ANCHOR_SUM_TOLERANCE,
                    nu=selected_nu,
                )
            catch
                nothing
            end
            if selected !== nothing &&
                    selected.estimated_relerr <= _STATE_ANCHOR_MAX_ERROR &&
                    all(_finite, (selected.X, selected.dXdx))
                return merge(selected, (
                    source=:complex_nu,
                    nu_delta,
                    selector_nu=selected_nu,
                    certified=true,
                    relaxed_sum=false,
                ))
            end
        end
    end

    fallback_tolerance =
        allow_relaxed && _allows_relaxed_state_anchor(params.omega) ?
        _STATE_ANCHOR_RELAXED_TOLERANCE : nothing
    principal = mst_anchor(
        source,
        :out,
        _STATE_ANCHOR_X;
        tolerance=_STATE_ANCHOR_SUM_TOLERANCE,
        fallback_tolerance,
    )
    relaxed_sum = hasproperty(principal, :relaxed_sum) &&
        principal.relaxed_sum
    principal_nu = ComplexF64(principal.nu)
    selected_nu === nothing && return merge(principal, (
        source=relaxed_sum ?
            :principal_relaxed_fallback : :principal_fallback,
        nu_delta=_nu_delta(principal_nu),
        selector_nu=principal_nu,
        certified=false,
        relaxed_sum,
    ))
    nu_delta = _nu_delta(selected_nu)
    nu_delta <= _MST_STATE_NU_DELTA_MIN &&
        return merge(principal, (
            source=relaxed_sum ? :principal_relaxed : :principal,
            nu_delta,
            selector_nu=selected_nu,
            certified=!relaxed_sum,
            relaxed_sum,
        ))
    return merge(principal, (
        source=relaxed_sum ?
            :principal_relaxed_fallback : :principal_fallback,
        nu_delta,
        selector_nu=selected_nu,
        certified=false,
        relaxed_sum,
    ))
end

function _p_state_anchor(
    route,
    evaluator;
    nu=nothing,
    require_certificate::Bool=false,
)
    params = evaluator.params
    mst_params = gsn_parameters(
        params.s,
        params.l,
        params.m,
        params.a,
        params.omega;
        lambda=params.lambda,
    )
    anchor = _mst_out_anchor(mst_params; nu)
    source = anchor.source
    !_anchor_accepted(anchor, params.omega, require_certificate) && return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:uncertified,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=Inf,
            mst_state_anchor_projective_error=Inf,
            mst_state_anchor_scale=ComplexF64(NaN, NaN),
        ),
    )
    anchor.estimated_relerr <= _STATE_ANCHOR_MAX_ERROR || return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:rejected_error,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=Inf,
            mst_state_anchor_projective_error=Inf,
            mst_state_anchor_scale=ComplexF64(NaN, NaN),
        ),
    )
    r = MSTInfinity._x_to_r(params, _STATE_ANCHOR_X)
    current_state = ComplexFrequency._state_at_r(evaluator, r)
    current = (current_state.X, current_state.dXdrstar)
    target = (
        anchor.X,
        anchor.dXdx * MSTInfinity._dx_drstar(
            params, _STATE_ANCHOR_X),
    )
    disagreement = all(_finite, current) ?
        Float64(_pair_error(current, target)) : Inf
    disagreement <= _STATE_ANCHOR_MAX_ERROR && return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:retained,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
            mst_state_anchor_projective_error=disagreement,
            mst_state_anchor_scale=ComplexF64(1),
        ),
    )
    correction = _fit_pair_scale(target, current)
    correction === nothing && return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:rejected_scale,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
            mst_state_anchor_projective_error=Inf,
            mst_state_anchor_scale=ComplexF64(NaN, NaN),
        ),
    )
    corrected = (
        correction * current[1],
        correction * current[2],
    )
    projective_error = Float64(_pair_error(corrected, target))
    projective_error <= _STATE_ANCHOR_MAX_ERROR || return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:rejected_projective,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
            mst_state_anchor_projective_error=projective_error,
            mst_state_anchor_scale=correction,
        ),
    )
    anchored = ComplexStateEvaluator(
        evaluator.params,
        evaluator.p_solution,
        evaluator.coefficients,
        ComplexF64(evaluator.scale * correction),
    )
    return (
        evaluator=anchored,
        metadata=(;
            mst_state_anchor_status=:rescaled,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
            mst_state_anchor_projective_error=projective_error,
            mst_state_anchor_scale=correction,
        ),
    )
end

function _small_frequency_expansion_state_anchor(
    route,
    evaluator;
    nu=nothing,
    require_certificate::Bool=false,
)
    coefficients = evaluator.route.coefficients
    anchor = _mst_out_anchor(coefficients; nu)
    source = anchor.source
    !_anchor_accepted(
        anchor,
        coefficients.params.omega,
        require_certificate,
    ) && return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:uncertified,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=Inf,
            mst_state_anchor_projective_error=Inf,
            mst_state_anchor_scale=ComplexF64(NaN, NaN),
        ),
    )
    anchor.estimated_relerr <= _STATE_ANCHOR_MAX_ERROR || return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:rejected_error,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=Inf,
            mst_state_anchor_projective_error=Inf,
            mst_state_anchor_scale=ComplexF64(NaN, NaN),
        ),
    )
    params = evaluator.route.params
    r = MSTInfinity._x_to_r(params, _STATE_ANCHOR_X)
    current_state = ComplexFrequency._state_at_r(evaluator, r)
    dxdrstar =
        MSTInfinity._dx_drstar(params, _STATE_ANCHOR_X)
    current = (current_state.X, current_state.dXdrstar)
    target = (anchor.X, anchor.dXdx * dxdrstar)
    disagreement = all(_finite, current) ?
        Float64(_pair_error(current, target)) : Inf
    disagreement <= _STATE_ANCHOR_MAX_ERROR && return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:retained,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
            mst_state_anchor_projective_error=disagreement,
            mst_state_anchor_scale=ComplexF64(1),
        ),
    )
    correction = _fit_pair_scale(target, current)
    correction === nothing && return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:rejected_scale,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
            mst_state_anchor_projective_error=Inf,
            mst_state_anchor_scale=ComplexF64(NaN, NaN),
        ),
    )
    corrected = (
        correction * current[1],
        correction * current[2],
    )
    projective_error = Float64(_pair_error(corrected, target))
    projective_error <= _STATE_ANCHOR_MAX_ERROR || return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:rejected_projective,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
            mst_state_anchor_projective_error=projective_error,
            mst_state_anchor_scale=correction,
        ),
    )
    return (
        evaluator=ScaledEvaluator(evaluator, correction),
        metadata=(;
            mst_state_anchor_status=:rescaled,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
            mst_state_anchor_projective_error=projective_error,
            mst_state_anchor_scale=correction,
        ),
    )
end

# Check the route state against an independent reference and rescale or
# replace its evaluator when they disagree: the real-axis propagation for IN,
# the MST outgoing anchor for UP. Returns (evaluator, metadata).
function _state_anchor(
    route;
    nu=nothing,
    require_certificate::Bool=false,
)
    route.branch == :IN &&
        return _in_real_axis_anchor(route, route.state_evaluator)
    route.branch == :UP || return (
        evaluator=route.state_evaluator,
        metadata=(;),
    )
    evaluator = route.state_evaluator
    evaluator isa ComplexSFEvaluator &&
        return _small_frequency_expansion_state_anchor(
            route,
            evaluator;
            nu,
            require_certificate,
        )
    evaluator isa ComplexStateEvaluator &&
        return _p_state_anchor(
            route,
            evaluator;
            nu,
            require_certificate,
        )
    evaluator isa ComplexRationalEvaluator || return (
        evaluator,
        metadata=(;),
    )
    anchor = _mst_out_anchor(
        evaluator.coefficients;
        nu,
        allow_relaxed=true,
    )
    source = anchor.source
    relaxed_sum = hasproperty(anchor, :relaxed_sum) &&
        anchor.relaxed_sum
    !_anchor_accepted(
        anchor,
        evaluator.coefficients.params.omega,
        require_certificate,
    ) && return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:uncertified,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=Inf,
        ),
    )
    anchor_error_limit = relaxed_sum ?
        _STATE_ANCHOR_RELAXED_TOLERANCE :
        _STATE_ANCHOR_MAX_ERROR
    anchor.estimated_relerr <= anchor_error_limit || return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:rejected_error,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=Inf,
        ),
    )
    current = try
        state = ComplexRational._state_x(
            evaluator, _STATE_ANCHOR_X)
        (state.X, state.dXdx)
    catch
        (ComplexF64(NaN, NaN), ComplexF64(NaN, NaN))
    end
    target = (anchor.X, anchor.dXdx)
    disagreement = all(_finite, current) ?
        Float64(_pair_error(current, target)) : Inf
    if relaxed_sum
        correction = _fit_pair_scale(target, current)
        correction === nothing && return (
            evaluator,
            metadata=(;
                mst_state_anchor_status=:rejected_scale,
                mst_state_anchor_source=source,
                mst_state_anchor_x=_STATE_ANCHOR_X,
                mst_state_anchor_error=anchor.estimated_relerr,
                mst_state_anchor_disagreement=disagreement,
                mst_state_anchor_projective_error=Inf,
                mst_state_anchor_scale=ComplexF64(NaN, NaN),
                mst_state_anchor_relaxed_sum=true,
            ),
        )
        corrected = (
            correction * current[1],
            correction * current[2],
        )
        projective_error =
            Float64(_pair_error(corrected, target))
        projective_error <= _STATE_ANCHOR_MAX_ERROR || return (
            evaluator,
            metadata=(;
                mst_state_anchor_status=:rejected_projective,
                mst_state_anchor_source=source,
                mst_state_anchor_x=_STATE_ANCHOR_X,
                mst_state_anchor_error=anchor.estimated_relerr,
                mst_state_anchor_disagreement=disagreement,
                mst_state_anchor_projective_error=projective_error,
                mst_state_anchor_scale=correction,
                mst_state_anchor_relaxed_sum=true,
            ),
        )
        return (
            evaluator=ScaledEvaluator(evaluator, correction),
            metadata=(;
                mst_state_anchor_status=:rescaled_relaxed,
                mst_state_anchor_source=source,
                mst_state_anchor_x=_STATE_ANCHOR_X,
                mst_state_anchor_error=anchor.estimated_relerr,
                mst_state_anchor_disagreement=disagreement,
                mst_state_anchor_projective_error=projective_error,
                mst_state_anchor_scale=correction,
                mst_state_anchor_relaxed_sum=true,
            ),
        )
    end
    disagreement <= _STATE_ANCHOR_MAX_ERROR && return (
        evaluator,
        metadata=(;
            mst_state_anchor_status=:retained,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
        ),
    )
    anchored = ComplexRationalEvaluator(
        evaluator.params,
        evaluator.coefficients,
        evaluator.settings,
        _STATE_ANCHOR_X,
        ComplexRationalState(anchor.X, anchor.dXdx),
    )
    return (
        evaluator=anchored,
        metadata=(;
            mst_state_anchor_status=:replaced,
            mst_state_anchor_source=source,
            mst_state_anchor_x=_STATE_ANCHOR_X,
            mst_state_anchor_error=anchor.estimated_relerr,
            mst_state_anchor_disagreement=disagreement,
        ),
    )
end

# Analytic MST amplitudes for the route: principal amplitudes for a rational
# evaluator, else refined amplitudes at the complex nu when certified. Returns
# (result, nu, kind) or nothing.
function _analytic_result(route)
    p = route.params
    evaluator = route.state_evaluator
    if evaluator isa ComplexRationalEvaluator
        principal = try
            mst_principal_amplitudes(
                evaluator.coefficients, route.branch)
        catch
            nothing
        end
        principal !== nothing &&
            return principal, nothing, :principal
    end
    refined_nu = try
        mst_nu_complex(
            p.s, p.l, p.m, p.a, p.omega, p.lambda)
    catch
        nothing
    end
    if refined_nu !== nothing
        refined = try
            mst_nia_amplitudes(
                p.s,
                p.l,
                p.m,
                p.a,
                p.omega,
                p.lambda,
                route.branch;
                nu=refined_nu,
            )
        catch
            nothing
        end
        refined !== nothing && refined.certificate_accepted &&
            return refined, refined.nu, :complex_refined
    end
    return nothing
end

# Preconditions of `_fast_route`: default controls, IN or UP, kappa above the
# P-consensus limit, Im(omega) < 0 with Re(omega) != 0, |omega| in
# (_FAST_ANCHOR_MIN_OMEGA, max_omega], and omega close to either axis.
@inline function _fast_anchor_pretrigger(
    a,
    omega,
    branch;
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
    max_omega::Float64=_FAST_ANCHOR_MAX_OMEGA,
)
    backend == :direct_rational || return false
    controls === nothing || return false
    # A tolerance alone sets accuracy targets and does not fix the route.
    all(value === nothing for value in (
        lambda,
        nu,
        xm,
        rhom,
        N,
        small_frequency_expansion,
        large_frequency_expansion,
    )) || return false
    pole_normalization && return false
    branch_symbol = try
        ComplexFrequency._normalize_branch(branch)
    catch
        return false
    end
    branch_symbol in (:IN, :UP) || return false
    af = try
        Float64(a)
    catch
        return false
    end
    kappa_squared = (1.0 - af) * (1.0 + af)
    kappa_squared > 0.0 || return false
    sqrt(kappa_squared) >
        ComplexFrequency._P_CONSENSUS_KAPPA_MAX || return false
    omegac = ComplexF64(omega)
    _finite(omegac) && imag(omegac) < 0 &&
        !iszero(real(omegac)) || return false
    omega_abs = abs(omegac)
    _FAST_ANCHOR_MIN_OMEGA < omega_abs <=
        max_omega || return false
    ratio = abs(real(omegac)) / abs(imag(omegac))
    return ratio <= _FAST_ANCHOR_RATIO_LOW ||
        ratio >= _FAST_ANCHOR_RATIO_HIGH
end

# Route built directly from certified principal (or refined) MST amplitudes and
# an MST state anchor at x = _STATE_ANCHOR_X, without a two-ray propagation.
# Returns nothing whenever a precondition, certificate or tolerance fails.
function _fast_route(
    s,
    l,
    m,
    a,
    omega,
    branch;
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
    max_omega::Float64=_FAST_ANCHOR_MAX_OMEGA,
    rescue_source::Symbol=:prebuild,
)
    _fast_anchor_pretrigger(
        a,
        omega,
        branch;
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
        max_omega,
    ) || return nothing
    branch_symbol = ComplexFrequency._normalize_branch(branch)

    params = try
        gsn_parameters(s, l, m, a, omega)
    catch
        return nothing
    end
    data = try
        MSTInfinity._mst_data(params)
    catch
        return nothing
    end
    result = try
        MSTInfinity._principal_amplitudes(data, branch_symbol)
    catch error
        error isa MSTInfinity.MSTCertificateError &&
            error.kind == :nia_amplitude_certificate ||
            return nothing
        try
            MSTInfinity._refined_amplitudes(
                data, branch_symbol)
        catch
            return nothing
        end
    end
    result.certificate_accepted || return nothing
    # The analytic amplitudes stand alone here, with nothing to compare them
    # with, so they are taken only when their error bound meets the tolerance;
    # otherwise the two-ray route is built and the overlay adjudicates.
    maximum(result.gsn_error) <=
        (tol === nothing ? _FAST_ANCHOR_TOLERANCE : tol) || return nothing
    state_data = if hasproperty(result, :nu)
        refined_nu = ComplexF64(result.nu)
        MSTInfinity._nia_series_data(
            MSTInfinity.MSTParams(
                data.params.s,
                data.params.l,
                data.params.m,
                data.params.a,
                data.params.omega,
                data.params.lambda,
                refined_nu,
            ),
            MSTInfinity.MST_NIA_TAIL,
        )
    else
        data
    end
    direct_controls = gsn_controls(
        params;
        N=_FAST_ANCHOR_ORDER,
        xm=_STATE_ANCHOR_X,
        tol=_FAST_ANCHOR_TOLERANCE,
        small_frequency_expansion=false,
        large_frequency_expansion=false,
    )
    coefficients = try
        gsn_coefficients(params; controls=direct_controls)
    catch
        return nothing
    end
    state = try
        if branch_symbol == :UP
            anchor = _mst_out_anchor(
                coefficients;
                nu=hasproperty(result, :nu) ? result.nu : nothing,
                allow_relaxed=true,
            )
            (
                X=ComplexF64(anchor.X),
                dXdx=ComplexF64(anchor.dXdx),
                estimated_relerr=Float64(anchor.estimated_relerr),
                nu=ComplexF64(anchor.nu),
                condition=Float64(anchor.condition),
                source=anchor.source,
                certified=Bool(anchor.certified),
            )
        else
            pin = if hasproperty(result, :nu_offset_dd)
                MSTInfinity._dd_pin_state(
                    coefficients,
                    state_data,
                    result.nu_offset_dd,
                    _STATE_ANCHOR_X;
                    scale=MSTInfinity._pin_unit_scale(params),
                )
            else
                plan = (
                    data=state_data,
                    converter=MSTInfinity._p_converter(params),
                    pin_norm=MSTInfinity._incoming_raw_factor(
                        state_data),
                    pin_scale=nothing,
                )
                MSTInfinity.mst_pin_state(
                    coefficients,
                    plan,
                    _STATE_ANCHOR_X;
                    scale=MSTInfinity._pin_unit_scale(params),
                )
            end
            (
                X=ComplexF64(pin.X),
                dXdx=ComplexF64(pin.dXdx),
                estimated_relerr=Float64(max(
                    pin.estimated_relerr,
                    pin.residual,
                )),
                nu=ComplexF64(NaN, NaN),
                condition=hasproperty(pin, :condition) ?
                    Float64(pin.condition) : NaN,
                source=:physical_in,
                representation=hasproperty(pin, :representation) ?
                    pin.representation : :native_type1,
                certified=true,
            )
        end
    catch
        return nothing
    end
    all(_finite, (state.X, state.dXdx)) || return nothing
    state.estimated_relerr <= _STATE_ANCHOR_MAX_ERROR ||
        return nothing
    _anchor_accepted(
        state,
        params.omega,
        false,
    ) || return nothing

    settings = ComplexRational.ComplexRationalSettings(
        _FAST_ANCHOR_ORDER, _FAST_ANCHOR_TOLERANCE)
    evaluator = ComplexRationalEvaluator(
        params,
        coefficients,
        settings,
        _STATE_ANCHOR_X,
        ComplexRationalState(state.X, state.dXdx),
    )

    complex_params = try
        ComplexFrequency.complex_parameters(
            params.s,
            params.l,
            params.m,
            params.a,
            params.omega,
            params.lambda,
        )
    catch
        return nothing
    end
    factors =
        ComplexFrequency._amplitude_factors(
            complex_params, branch_symbol)
    transmission = ComplexF64(1)
    incidence = ComplexF64(result.gsn[1])
    reflection = ComplexF64(result.gsn[2])
    teukolsky_transmission = ComplexF64(factors.transmission)
    teukolsky_incidence =
        ComplexF64(teukolsky_transmission * result.teuk[1])
    teukolsky_reflection =
        ComplexF64(teukolsky_transmission * result.teuk[2])
    all(_finite, (
        transmission,
        incidence,
        reflection,
        teukolsky_transmission,
        teukolsky_incidence,
        teukolsky_reflection,
    )) || return nothing

    requested = ComplexFrequency.complex_controls()
    resolved = (
        xm=_STATE_ANCHOR_X,
        rhom=nothing,
        N=_FAST_ANCHOR_ORDER,
        tolerance=_FAST_ANCHOR_TOLERANCE,
        small_frequency_expansion=false,
        large_frequency_expansion=false,
        rescue_source=:none,
    )
    contour = ComplexFrequency.complex_contour_parameters(
        complex_params; rho_match=1.0)
    matching = (
        requested_x=nothing,
        allowed_min_x=0.0,
        allowed_max_x=1.0,
        candidate_x=_STATE_ANCHOR_X,
        selected_x=_STATE_ANCHOR_X,
        split_mismatch=state.estimated_relerr,
        condition=1.0,
        policy=:mst_anchor,
    )
    amplitude_summary = (
        teukolsky=(
            transmission=teukolsky_transmission,
            incidence=teukolsky_incidence,
            reflection=teukolsky_reflection,
        ),
        conversion_factors=factors,
        raw_gsn_transmission=transmission,
        gsn=(;
            transmission,
            incidence,
            reflection,
        ),
    )
    plan = ComplexRoutePlan(
        contour,
        matching,
        resolved,
        amplitude_summary,
        0,
    )
    metadata = (
        backend=:direct_gsn_mst_anchor,
        match_policy=:mst_anchor,
        match_x=_STATE_ANCHOR_X,
        xsplit=_STATE_ANCHOR_X,
        xm_match=_STATE_ANCHOR_X,
        matching_condition=1.0,
        split_mismatch=state.estimated_relerr,
        abel_error=state.estimated_relerr,
        patch_count=0,
        control_xm=_STATE_ANCHOR_X,
        control_rhom=nothing,
        control_N=_FAST_ANCHOR_ORDER,
        control_tol=_FAST_ANCHOR_TOLERANCE,
        control_small_frequency_expansion=false,
        control_large_frequency_expansion=false,
        amplitude_backend=:principal_mst_anchor,
        amplitude_patch_count=0,
        mst_amplitude_certificate=result.certificate_kind,
        mst_amplitude_representation_spread=
            result.representation_spread,
        mst_amplitude_nearest_agreement=result.nearest_agreement,
        mst_amplitude_truncation_agreement=
            result.truncation_agreement,
        mst_amplitude_max_condition=result.max_condition,
        mst_amplitude_errors=hasproperty(result, :gsn_error) ?
            Float64.(result.gsn_error) : (NaN, NaN),
        mst_amplitude_component_agreement=
            hasproperty(result, :component_agreement) ?
            Float64.(result.component_agreement) : (NaN, NaN),
        mst_amplitude_coefficient_residual=
            hasproperty(result, :coefficient_residual) ?
            Float64(result.coefficient_residual) : NaN,
        mst_amplitude_coefficient_source_shift=
            hasproperty(result, :coefficient_source_shift) ?
            Int(result.coefficient_source_shift) : 0,
        mst_amplitude_nu_residual=
            hasproperty(result, :nu_residual) ?
            result.nu_residual : NaN,
        mst_state_anchor_status=:accepted,
        mst_state_anchor_source=state.source,
        mst_state_anchor_x=_STATE_ANCHOR_X,
        mst_state_anchor_error=state.estimated_relerr,
        mst_state_anchor_nu=state.nu,
        mst_state_anchor_condition=state.condition,
        mst_state_anchor_certified=state.certified,
        mst_state_anchor_representation=hasproperty(
            state, :representation) ?
            state.representation : :unknown,
        mst_fast_rescue_source=rescue_source,
    )
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

# Whether an IN route from the P-route fallback or consensus backend, with
# _FAST_ANCHOR_MAX_OMEGA < |omega| <= max_omega, is retried with `_fast_route`.
@inline function _fast_p_pretrigger(
    route::ComplexRoute;
    max_omega::Float64=_MAX_ABS_IMAG_OMEGA,
)
    route.branch == :IN || return false
    route_backend = hasproperty(route.metadata, :backend) ?
        route.metadata.backend : :unknown
    route_backend in (
        :direct_gsn_p_build_fallback,
        :direct_gsn_p_consensus,
    ) || return false
    omega_abs = abs(ComplexF64(route.params.omega))
    return _FAST_ANCHOR_MAX_OMEGA < omega_abs <= max_omega
end

function _fast_p_rescue(
    route::ComplexRoute,
    s,
    l,
    m,
    a,
    omega,
    branch;
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
    max_omega::Float64=_MAX_ABS_IMAG_OMEGA,
)
    _fast_p_pretrigger(route; max_omega) || return route
    route_backend = route.metadata.backend
    candidate = _fast_route(
        s,
        l,
        m,
        a,
        omega,
        branch;
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
        max_omega,
        rescue_source=route_backend,
    )
    return candidate === nothing ? route : candidate
end

# Whether the analytic-amplitude comparison runs after a route was built:
# default controls, Im(omega) < 0, and either omega near the negative imaginary
# axis or a route whose amplitudes came from the complex-rational or
# log-scaled MST backends.
@inline function _offaxis_pretrigger(
    route::ComplexRoute,
    ;
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
    backend == :direct_rational || return false
    controls === nothing || return false
    # A tolerance alone sets accuracy targets and does not fix the route.
    all(value === nothing for value in (
        lambda,
        nu,
        xm,
        rhom,
        N,
        small_frequency_expansion,
        large_frequency_expansion,
    )) || return false
    pole_normalization && return false
    omega = ComplexF64(route.params.omega)
    imag(omega) < 0 || return false
    !iszero(real(omega)) || return false
    near_nia = abs(real(omega)) <=
        _NIA_LATERAL_REAL_RATIO_MAX * max(abs(imag(omega)), 1.0)
    near_nia && return true
    amplitude_backend = hasproperty(route.metadata, :amplitude_backend) ?
        route.metadata.amplitude_backend : :unknown
    match_policy = hasproperty(route.metadata, :match_policy) ?
        route.metadata.match_policy : :unknown
    match_policy == :branch_shift_consensus && return false
    return amplitude_backend in (
        :complex_rational,
        :principal_mst_logscaled,
    )
end

@inline function _state_pretrigger(
    route::ComplexRoute,
    ;
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
    backend == :direct_rational || return false
    controls === nothing || return false
    # A tolerance alone sets accuracy targets and does not fix the route.
    all(value === nothing for value in (
        lambda,
        nu,
        xm,
        rhom,
        N,
        small_frequency_expansion,
        large_frequency_expansion,
    )) || return false
    pole_normalization && return false
    route.branch == :UP || return false
    route.state_evaluator isa ComplexStateEvaluator || return false
    omega = ComplexF64(route.params.omega)
    imag(omega) < 0 || return false
    !iszero(real(omega)) || return false
    route_backend = hasproperty(route.metadata, :backend) ?
        route.metadata.backend : :unknown
    return route_backend == :direct_gsn_p_build_fallback
end

# Post-route hook registered with ComplexFrequency: tries the fast
# MST-anchor rescue, then (when triggered) compares the route amplitudes with
# the analytic MST ones and re-anchors the UP/IN state. Returns the route
# unchanged when nothing applies.
@noinline function _post_overlay(
    route::ComplexRoute,
    s,
    l,
    m,
    a,
    omega,
    branch;
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
    Base.@nospecialize route
    rescued = _fast_p_rescue(
        route,
        s,
        l,
        m,
        a,
        omega,
        branch;
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
        max_omega=_MAX_ABS_IMAG_OMEGA,
    )
    rescued === route || return rescued
    amplitude_trigger = _offaxis_pretrigger(
        route;
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
    state_trigger = _state_pretrigger(
        route;
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
    (amplitude_trigger || state_trigger) || return route
    # Function barrier: the triggered work is compiled only when it runs.
    return Base.inferencebarrier(_post_overlay_triggered)(
        route, amplitude_trigger, state_trigger)
end

# Error of the analytic amplitudes from their condition and the agreement of
# their nearest and truncated representations, in units of the larger
# amplitude (the representation distances are max norms of the pair).
_analytic_error(result) = max(eps(Float64) * result.max_condition,
    result.nearest_agreement, result.truncation_agreement)

# Relative error bound of each analytic amplitude: its own relative distance
# from the other evaluations it was compared with (another representation and
# the truncated sums), plus the running error of its own sums (a distance to
# another evaluation bounds the error only together with that evaluation's own
# error). Where only the pair distance is recorded, it is in units of the
# larger amplitude and is rescaled to each amplitude.
function _analytic_bound(result)
    own = hasproperty(result, :gsn_error) ? result.gsn_error : (0.0, 0.0)
    hasproperty(result, :component_agreement) && return own .+ max.(
        eps(Float64) * result.max_condition, result.component_agreement)
    scale = maximum(abs, result.gsn)
    return own .+ _analytic_error(result) .* scale ./ abs.(result.gsn)
end

# Which of the two-ray amplitudes u the analytic ones v replace, and the
# errors to report for the kept ones. v replaces u only where it proves u
# wrong: |v - u| > 2 E(v) gives |u - true| > E(v) >= |v - true|, with E(v) the
# running error bound of the analytic sums and their combinations (eps times
# the largest condition alone misses the per-term rounding and the
# cancellation of the final combinations). The reflection coefficient is
# judged on its own relative scale. The incidence coefficient, which vanishes
# at a mode frequency where the analytic bound is not reliable in relative
# terms, is judged with the pair as a vector (max norm, units of the larger
# amplitude). A kept u has an error of at least |v - u| - E(v) and at most
# |v - u| + E(v).
function _rational_preferred(route, result)
    metadata = route.metadata
    hasproperty(metadata, :incidence_error) &&
        hasproperty(metadata, :reflection_error) || return nothing
    _finite(route.transmission) && !iszero(route.transmission) ||
        return nothing
    bound = _analytic_bound(result)
    current = (route.incidence / route.transmission,
        route.reflection / route.transmission)
    differences = abs.(current .- result.gsn) ./ abs.(result.gsn)
    scale = maximum(abs, result.gsn)
    pair = maximum(differences .* abs.(result.gsn)) / scale >
        2 * maximum(bound .* abs.(result.gsn)) / scale
    replace = (pair, differences[2] > 2 * bound[2])
    # A kept value is within its distance to the analytic value plus the
    # analytic value's own error. An analytic value less precise than the
    # kept value's own estimate can only contradict that estimate, when the
    # two disagree by more than both errors.
    own = (metadata.incidence_error, metadata.reflection_error)
    relative = abs.(result.gsn) ./ abs.(current)
    analytic = bound .* relative
    distance = differences .* relative
    reported = ifelse.(
        (analytic .< own) .| (distance .> own .+ analytic),
        distance .+ analytic, own)
    return (; replace, reported, bound)
end

function _post_overlay_triggered(route, amplitude_trigger, state_trigger)
    Base.@nospecialize route
    analytic = amplitude_trigger ? _analytic_result(route) : nothing
    p = route.params
    near_nia = abs(real(p.omega)) <=
        _NIA_LATERAL_REAL_RATIO_MAX * max(abs(imag(p.omega)), 1.0)
    if analytic === nothing
        (near_nia || state_trigger) || return route
        anchor = try
            _state_anchor(route)
        catch
            (evaluator=route.state_evaluator, metadata=(;))
        end
        return _state_overlay(route, anchor)
    end
    result, nu_value, _ = analytic
    amplitude_agreement = if _finite(route.transmission) &&
            !iszero(route.transmission)
        current = (
            route.incidence / route.transmission,
            route.reflection / route.transmission,
        )
        _pair_error(current, result.gsn)
    else
        Inf
    end
    decision = _rational_preferred(route, result)
    # The per-coefficient decision governs where there is one: a pair that
    # agrees in units of the larger amplitude can still be wrong in the
    # smaller one.
    retain = decision === nothing ?
        amplitude_agreement <= _AMPLITUDE_AGREEMENT_MAX :
        !any(decision.replace)
    if retain
        decision === nothing || (route = _with_amplitude_errors(route,
            decision.reported))
        anchor = try
            _state_anchor(
                route;
                nu=nu_value,
                require_certificate=true,
            )
        catch
            (evaluator=route.state_evaluator, metadata=(;))
        end
        return _retain_analytic_agreement(
            _state_overlay(route, anchor),
            result,
            amplitude_agreement,
        )
    end
    # The real-axis state is chosen by its own certificate, as on the retain
    # path; which amplitude values are reported does not change it.
    anchor = try
        _state_anchor(route; nu=nu_value, require_certificate=true)
    catch
        (evaluator=route.state_evaluator, metadata=(;))
    end
    anchored_route = _state_overlay(route, anchor)
    near_nia && !(result.certificate_kind in (
        :global_medoid,
        :principal_single_nu,
    )) &&
        return anchored_route
    decision === nothing || all(decision.replace) ||
        return _overlay_part(anchored_route, route, result, decision, anchor)
    return _overlay(
        anchored_route,
        result;
        amplitude_backend=:complex_mst_analytic,
        state_evaluator=anchored_route.state_evaluator,
        state_metadata=anchor.metadata,
    )
end

# The analytic value of one coefficient with the two-ray value of the other
# (each kept with its own error), applied as an overlay; the state keeps the
# anchored evaluator and is not rebuilt from the pair.
function _overlay_part(anchored_route, route, result, decision, anchor)
    two_ray = (route.incidence, route.reflection) ./ route.transmission
    two_ray_teukolsky = (route.teukolsky_incidence,
        route.teukolsky_reflection) ./ route.teukolsky_transmission
    pick(analytic, own) = ifelse.(decision.replace, analytic, own)
    part = merge(result, (;
        gsn=pick(result.gsn, two_ray),
        teuk=pick(result.teuk, two_ray_teukolsky),
    ))
    errors = pick(decision.bound, decision.reported)
    overlaid = _overlay(
        anchored_route,
        part;
        amplitude_backend=:complex_mst_partial,
        state_evaluator=anchored_route.state_evaluator,
        state_metadata=anchor.metadata,
    )
    return _with_amplitude_errors(overlaid, errors)
end

# Special route builder registered with ComplexFrequency: returns a route
# from the fast MST-anchor path, or for default controls with omega near the
# negative imaginary axis the two-ray route with the analytic MST amplitudes
# overlaid; nothing lets the caller build its usual route.
@noinline function build(
    s,
    l,
    m,
    a,
    omega,
    branch;
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
    fast = _fast_route(
        s,
        l,
        m,
        a,
        omega,
        branch;
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
        max_omega=_MAX_ABS_IMAG_OMEGA,
    )
    fast === nothing || return fast
    _pretrigger(
        s,
        a,
        omega;
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
    ) || return nothing
    # Function barrier: the triggered work is compiled only when it runs.
    return Base.inferencebarrier(_build_triggered)(s, l, m, a, omega, branch)
end

function _build_triggered(s, l, m, a, omega, branch)
    route = try
        ComplexFrequency._complex_p_route(
            s, l, m, a, omega, branch)
    catch
        return nothing
    end
    result = try
        mst_nia_amplitudes(
            route.params.s,
            route.params.l,
            route.params.m,
            route.params.a,
            route.params.omega,
            route.params.lambda,
            route.branch,
        )
    catch
        return nothing
    end
    result.certificate_accepted || return nothing
    anchor = try
        _state_anchor(route)
    catch
        (evaluator=route.state_evaluator, metadata=(;))
    end
    return _overlay(
        route,
        result;
        state_evaluator=anchor.evaluator,
        state_metadata=anchor.metadata,
    )
end

ComplexFrequency._register_special_route_builder!(build)
ComplexFrequency._register_post_route_builder!(_post_overlay)

end
