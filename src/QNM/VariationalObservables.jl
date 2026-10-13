# Keep the fast route within a bounded mode family, with per-call checks. No root
# finding, mirror construction or extremal-family classification happens here.
function _qnm_variational_scope(root)
    m = root.mode
    return m.s in -2:2 && abs(m.s) <= m.l <= 4 && 0 <= m.n <= 3 &&
        abs(root.a) <= 1 && isfinite(root.omega) &&
        !iszero(real(root.omega)) && imag(root.omega) < 0 &&
        !(get(root.provenance, :extremal_endpoint, false) &&
          !get(root.provenance, :simple_pole, true))
end

function _qnm_use_variational(root, options)
    backend = get(options, :backend, :auto)
    backend in (:auto, :variational, :richardson) ||
        throw(ArgumentError("excitation backend must be :auto, :variational or :richardson"))
    compatible = get(options, :h, nothing) === nothing &&
        get(options, :direction_tolerance, nothing) === nothing &&
        get(options, :radial_solver, GSN_radial) in (GSN_radial, _qnm_default_radial_solver) &&
        get(options, :teukolsky_solver, Teukolsky_radial) in (nothing, Teukolsky_radial)
    if backend === :variational && (!compatible || !_qnm_variational_scope(root))
        throw(ArgumentError("variational backend supports default solvers, no explicit stencil/direction check, and integer s=-2:2,abs(s)<=l<=4,n=0:3 simple poles"))
    end
    return backend !== :richardson && compatible && _qnm_variational_scope(root)
end

_without_excitation_backend(options) =
    (; (k => v for (k,v) in pairs(options) if k !== :backend)...)

function _qnm_checked_variational(root; h=nothing, radial_solver=GSN_radial,
        teukolsky_solver=Teukolsky_radial, incidence_tolerance=1e-8,
        step_tolerance=1e-7, wronskian_tolerance=1e-8,
        direction_tolerance=nothing, bridge_tolerance=nothing)
    target = min(Float64(step_tolerance), 1e-8)
    all(t -> isfinite(t) && t > 0,
        (incidence_tolerance, target, wronskian_tolerance)) ||
        throw(ArgumentError("excitation tolerances must be positive and finite"))
    inputs = _radial_evaluation_inputs(root)
    # Start farther from the regular singular point with more series terms;
    # near extremality this reduces near-horizon ODE step rejection.
    extreme = abs(root.a) == 1
    af,wf=Float64(inputs.a),ComplexF64(inputs.omega)
    kappa=sqrt(max(0.,1-af*af))
    p=wf-root.mode.m*af/(2*(1+kappa))
    complex_horizon=!extreme && abs(af)>.9 && abs(real(p))>.1 &&
        abs(imag(wf))/kappa>.5
    # Spin-field endpoint conversion is ill-conditioned at extremality.
    # Only the small boundary calculation uses higher precision, not the IVP.
    bits = extreme && root.mode.s != 0 ? 128 : 53
    coarse = _qnm_variational_excitation(root; horizon_scale=0.4,
        horizon_order=24, infinity_radius=80.0, horizon_ray=extreme,
        extremal_scale=24.0, extremal_radius=0.2, extremal_order=40,
        extremal_precision=bits,horizon_tortoise=complex_horizon)
    fine = _qnm_variational_excitation(root; tol=3e-12, horizon_scale=0.8,
        horizon_order=40, infinity_order=40, infinity_radius=100.0,
        horizon_ray=extreme, extremal_scale=32.0, extremal_radius=0.15,
        extremal_order=64, extremal_precision=bits,horizon_tortoise=complex_horizon)
    drift(x,y) = abs(x-y)/max(abs(x), abs(y), floatmin(Float64))
    refinement = maximum((drift(coarse.alpha,fine.alpha),
        drift(coarse.reflection,fine.reflection),drift(coarse.B_gsn,fine.B_gsn)))
    scaled = abs(coarse.incidence)/max(1,abs(coarse.reflection))
    fine_scaled = abs(fine.incidence)/max(1,abs(fine.reflection))
    wronskian = max(coarse.wronskian_drift,fine.wronskian_drift)
    finite = all(isfinite,(coarse.incidence,coarse.reflection,coarse.alpha,
        coarse.B_gsn,fine.incidence,fine.reflection,fine.alpha,fine.B_gsn,
        refinement,wronskian))
    valid = finite && max(scaled,fine_scaled) <= incidence_tolerance &&
        wronskian <= wronskian_tolerance
    reason = !finite ? :nonfinite_variational_excitation :
        max(scaled,fine_scaled) > incidence_tolerance ? :scaled_incidence_gate :
        wronskian > wronskian_tolerance ? :wronskian_constancy_gate : :accepted
    validation = ISEMValidationResult(root.omega,coarse.incidence,
        coarse.reflection,one(coarse.incidence),abs(coarse.incidence),scaled,
        :variational,valid ? :accepted : :failed,reason,
        (incidence_tolerance,wronskian_tolerance,wronskian_drift=wronskian,
         incidence_backend=:variational_connection,
         inferred_root_shift=abs(coarse.incidence/coarse.alpha),
         evaluation_a=inputs.a,evaluation_omega=inputs.omega,
         evaluation_policy=inputs.policy,frequency_projection_drift=inputs.projection_drift))
    accepted = valid && refinement <= target
    return (omega=root.omega,alpha=coarse.alpha,reflection=coarse.reflection,
        B_gsn=coarse.B_gsn,step_drift=refinement,
        status=accepted ? :accepted : :failed,
        stop_reason=!valid ? reason : accepted ? :accepted : :variational_refinement_gate,
        metadata=(validation,backend=:variational,
            derivative_direction=:holomorphic_frequency_tangent,
            precision_basis=:same_frequency_resolution_indicator_not_error_bound,
            step_tolerance=target,maximum_wronskian_drift=wronskian,wronskian_tolerance,
            incidence_backend=:variational_connection,
            evaluation_policy=inputs.policy,frequency_projection_drift=inputs.projection_drift,
            evaluation_a=inputs.a,evaluation_omega=inputs.omega,
            radial_sample_requests=2,radial_unique_samples=1,
            augmented_ivps=extreme ? 6 : 4,shifted_frequency_solves=0,
            horizon_endpoint_precision_bits=bits,
            horizon_coordinate=complex_horizon ? :complex_tortoise : :radius,
            refinement_alpha=drift(coarse.alpha,fine.alpha),
            refinement_reflection=drift(coarse.reflection,fine.reflection),
            refinement_B=drift(coarse.B_gsn,fine.B_gsn),
            lambda=coarse.lambda,dlambda=coarse.dlambda,
            angular_seconds=coarse.angular_s+fine.angular_s,
            endpoint_seconds=coarse.endpoint_s+fine.endpoint_s,
            propagation_seconds=coarse.propagation_s+fine.propagation_s))
end

# The detailed result shares precisely the compact GSN calculation. The
# Teukolsky bridge is analytic, not a second independent radial validation.
function _qnm_variational_teukolsky(root, gsn; bridge_tolerance=1e-10, kwargs...)
    a,w = gsn.metadata.evaluation_a,gsn.metadata.evaluation_omega
    s,m = root.mode.s,root.mode.m
    wd = _frequency_dual(w)
    ld = _frequency_dual(gsn.metadata.lambda,gsn.metadata.dlambda)
    bt = ConversionFactors.Btrans(s,m,a,wd,ld)
    ci = ConversionFactors.Binc(s,m,a,wd,ld)/bt
    cr = ConversionFactors.Bref(s,m,a,wd,ld)/bt
    # Preserve c'(omega)*Ain away from an exact zero.
    alpha_t = _primal(ci)*gsn.alpha + _tangent(ci)*gsn.metadata.validation.incidence
    reflection_t = _primal(cr)*gsn.reflection
    B_t = reflection_t/(2im*w*alpha_t)
    pole_B = _primal(cr)/(im*_primal(ci))*gsn.B_gsn
    residual = abs(B_t-pole_B)/max(abs(B_t),abs(pole_B),floatmin(Float64))
    accepted = gsn.status === :accepted && isfinite(residual) && residual <= bridge_tolerance
    return ExcitationFactorResult(root.omega,gsn.alpha,alpha_t,gsn.reflection,
        reflection_t,gsn.B_gsn,B_t,NaN,gsn.step_drift,residual,
        accepted ? :accepted : :failed,
        gsn.status !== :accepted ? gsn.stop_reason :
            accepted ? :accepted : :normalization_bridge_gate,
        merge(gsn.metadata,(bridge_tolerance,
            bridge_method=:analytic_conversion_with_off_root_product_rule,
            bridge_independence=:shared_gsn_solution_not_independent_solver,
            direction_check=:not_sampled_holomorphic_tangent,
            convention_gsn=:Aout_over_2omega_alpha,
            convention_teukolsky=:Aout_over_2iomega_alpha)))
end

function _qnm_selected_excitation(root, options; detailed=false)
    legacy = _without_excitation_backend(options)
    if _qnm_use_variational(root, options)
        candidate = try
            _qnm_checked_variational(root; legacy...)
        catch err
            err isa InterruptException && rethrow()
            get(options,:backend,:auto) === :variational && rethrow()
            nothing
        end
        if candidate !== nothing && candidate.status === :accepted
            return detailed ? _qnm_variational_teukolsky(root,candidate;legacy...) : candidate
        elseif get(options,:backend,:auto) === :variational
            return detailed ? _qnm_variational_teukolsky(root,candidate;legacy...) : candidate
        end
        # Failed numerical checks never promote a fast result. The unchanged
        # legacy algorithm remains available; record the rejected attempt.
        result = detailed ? _qnm_richardson_full_excitation_factor(root;legacy...) :
            _qnm_richardson_excitation_factor(root;legacy...)
        reason = candidate === nothing ? :variational_exception : candidate.stop_reason
        if !detailed
            return merge(result,(metadata=merge(result.metadata,
                (variational_fallback_reason=reason,)),))
        end
        return ExcitationFactorResult(result.omega,result.alpha,result.alpha_teukolsky,
            result.reflection,result.reflection_teukolsky,result.B_gsn,result.B_teukolsky,
            result.direction_drift,result.step_drift,result.bridge_residual,result.status,
            result.stop_reason,merge(result.metadata,(variational_fallback_reason=reason,)))
    end
    return detailed ? _qnm_richardson_full_excitation_factor(root;legacy...) :
        _qnm_richardson_excitation_factor(root;legacy...)
end

_qnm_gsn_excitation_factor(root::LeaverResult;kwargs...) =
    _qnm_selected_excitation(root,(;kwargs...))

"""
    qnm_excitation_factor(root::LeaverResult; backend=:auto, kwargs...)

Use checked frequency sensitivities for integer spins -2 through 2, l<=4,
n<=3 simple poles, or the original Richardson algorithm otherwise. Failed
refinement or contour checks fall back to Richardson. Force either with
`backend=:variational` or `:richardson`. Explicit stencils/custom solvers
retain Richardson. The variational Teukolsky bridge includes the off-root
product rule and does not claim an independent radial or direction check.
"""
qnm_excitation_factor(root::LeaverResult;kwargs...) =
    _qnm_selected_excitation(root,(;kwargs...);detailed=true)
