function _radial_metadata(solution)
    numerical = solution.numerical_GSN_solution
    numerical_summary = if numerical isa NamedTuple
        propertynames(numerical)
    elseif numerical === missing || numerical === nothing
        ()
    else
        (type=string(typeof(numerical)),)
    end
    numerical_value(name, default=missing) =
        numerical isa NamedTuple && hasproperty(numerical, name) ?
        getproperty(numerical, name) : default
    return (
        method=solution.method,
        normalization=solution.normalization_convention,
        lambda=solution.mode.lambda,
        numerical_summary=numerical_summary,
        control_tolerance=numerical_value(:control_tolerance),
        horizon_endpoint_score=numerical_value(:horizon_endpoint_score),
        infinity_endpoint_score=numerical_value(:infinity_endpoint_score),
        propagation_score=numerical_value(:propagation_score),
        horizon_path_kind=numerical_value(:horizon_path_kind),
        horizon_coordinate_steps=numerical_value(:horizon_coordinate_steps),
        solution_type=string(typeof(solution)),
    )
end

function _qnm_conditioned_radial_tolerance(a)
    af = Float64(real(a))
    isfinite(af) && abs(af) < 1 || return nothing
    kappa = sqrt(max(0.0, 1 - af^2))
    kappa <= 0.02 || return nothing
    conditioned = 64eps(Float64) / (1000 * max(kappa^2, eps(Float64)))
    return clamp(conditioned, 1.0e-14, 1.0e-11)
end

# Radial solver used for QNM amplitudes: GSN-ISEM with a fixed collocation order
# away from extremality; the solver default with a kappa-dependent tolerance for
# kappa <= 0.02, and with its own defaults at abs(a) = 1.
function _qnm_default_radial_solver(s, l, m, a, omega, boundary;
        solver=GSN_radial)
    tolerance = _qnm_conditioned_radial_tolerance(a)
    if tolerance === nothing && abs(a) < 1
        # Keep the connection normalization fixed across the derivative
        # stencil; the automatic MST amplitude overlay uses another route.
        order = abs(imag(omega)) > 1 ? 96 : 48
        solution = try
            solver(s, l, m, a, omega, boundary;
                method="GSN-ISEM", N=order, tol=1.0e-14)
        catch error
            # A straight tortoise-coordinate ray can hit a coordinate
            # branch point. Move its real junction, without changing the ODE.
            occursin("two-ray coordinate path", sprint(showerror,error)) || rethrow()
            return solver(s, l, m, a, omega, boundary;
                method="GSN-ISEM", N=order, tol=1.0e-14, xm=0.3)
        end
        gsn = solver === Teukolsky_radial ? solution.GSN_solution : solution
        numerical = gsn.numerical_GSN_solution
        if hasproperty(numerical, :infinity_horizon_winding) &&
                numerical.infinity_horizon_winding != 0
            # Removing monodromy cannot recover digits already lost on
            # a winding path. Prefer an unwound real junction when available.
            alternate = solver(s, l, m, a, omega, boundary;
                method="GSN-ISEM", N=order, tol=1.0e-14, xm=0.3)
            alternate_gsn = solver === Teukolsky_radial ? alternate.GSN_solution : alternate
            alternate_gsn.numerical_GSN_solution.infinity_horizon_winding == 0 &&
                return alternate
        end
        return solution
    end
    return tolerance === nothing ?
        solver(s, l, m, a, omega, boundary) :
        solver(s, l, m, a, omega, boundary; tol=tolerance)
end

function _radial_evaluation_inputs(result::LeaverResult, omega=result.omega)
    if result.a isa BigFloat || real(omega) isa BigFloat
        evaluation_a = Float64(result.a)
        evaluation_omega = ComplexF64(omega)
        T = typeof(real(omega))
        projected = Complex{T}(
            T(real(evaluation_omega)), T(imag(evaluation_omega)))
        drift = T(abs(projected - omega) /
            max(one(T), abs(projected), abs(omega)))
        return (
            a=evaluation_a,
            omega=evaluation_omega,
            policy=:explicit_binary64_projection,
            projection_drift=drift,
            source_precision_bits=result.precision_bits,
            evaluation_precision_bits=53,
        )
    end
    T = typeof(real(omega))
    return (
        a=result.a,
        omega=omega,
        policy=:native_input_precision,
        projection_drift=zero(T),
        source_precision_bits=result.precision_bits,
        evaluation_precision_bits=precision(real(omega)),
    )
end

# Use the IN/UP Wronskian for the incidence amplitude at exact extremality and
# for conditioned (estimated) roots with abs(a) >= 0.9999 and n >= 4.
function _use_wronskian_incidence(result::LeaverResult, omega=result.omega)
    result.status == :estimated && abs(result.a) >= oftype(result.a, 0.9999) &&
        result.mode.n >= 4 && return true
    return abs(result.a) == one(result.a) && !isreal(omega)
end

function _is_extremal_damped_simple_pole(result::LeaverResult)
    provenance = result.provenance
    return abs(result.a) == one(result.a) &&
        hasproperty(provenance, :extremal_endpoint) &&
        provenance.extremal_endpoint &&
        hasproperty(provenance, :simple_pole) && provenance.simple_pole &&
        hasproperty(provenance, :family) && provenance.family == :DM
end

# At a = +-1 the horizon is an irregular singular point. The horizon state
# e^{alpha/h} h^beta sum c_n h^n is asymptotic in h = r - 1 with
# alpha = i(2 omega - m a); the integration starts at h = 2|alpha|/40 so that the
# series is used inside its region of usefulness for every detuning, instead of
# a fixed r* that fails when 2 omega is close to m. For damped frequencies with
# Re(alpha) > 0 the outgoing solution still grows against the ingoing one by
# e^{2 Re(alpha)/h}; in Float64 this limits the incidence closure for corotating
# damped modes, which the incidence gate reports.
function _extremal_qnm_horizon_start(m, a, omega; series_terms=40,
        minimum_h=1.0e-6, maximum_h=0.05)
    detuning = abs(2omega - m * sign(a))
    h = clamp(2detuning / series_terms, minimum_h, maximum_h)
    return rstar_from_r(a, one(h) + h)
end

function _extremal_qnm_radial_solver(s, l, m, a, omega, boundary)
    rsin = _extremal_qnm_horizon_start(m, a, ComplexF64(omega))
    return GSN_radial(
        s, l, m, a, omega, boundary, rsin, 60.0;
        method="GSN-ISEM",
        tolerance=1.0e-14,
        horizon_expansion_order=64,
        infinity_expansion_order=32,
    )
end

# Incidence amplitude from the IN/UP Wronskian, averaged over `evaluation_points`
# in r*; their spread is returned as `wronskian_drift`.
function _wronskian_incidence(result::LeaverResult, omega=result.omega;
        radial_solver=GSN_radial, evaluation_points=(0.0, 1.0))
    inputs = _radial_evaluation_inputs(result, omega)
    mode = result.mode
    Xin = radial_solver(
        mode.s, mode.l, mode.m, inputs.a, inputs.omega, IN)
    Xup = radial_solver(
        mode.s, mode.l, mode.m, inputs.a, inputs.omega, UP)
    incidence_values = map(evaluation_points) do rstar
        xin = Xin.GSN_solution(rstar)
        xup = Xup.GSN_solution(rstar)
        raw_wronskian = xin[1] * xup[2] - xup[1] * xin[2]
        radius = r_from_rstar(inputs.a, rstar)
        eta = Transformation.eta(
            mode.s, mode.m, inputs.a, inputs.omega,
            Xin.mode.lambda, radius)
        eta_infinity = Transformation.eta_coefficient(
            mode.s, mode.m, inputs.a, inputs.omega,
            Xin.mode.lambda, 0)
        eta_infinity * raw_wronskian /
            (2im * inputs.omega * eta)
    end
    incidence = sum(incidence_values) / length(incidence_values)
    scale = max(one(abs(incidence)), abs(Xin.reflection_amplitude),
        abs(Xin.transmission_amplitude), maximum(abs, incidence_values))
    wronskian_drift = maximum(
        abs(value - incidence) / scale for value in incidence_values)
    return (
        incidence=incidence,
        reflection=Xin.reflection_amplitude,
        transmission=Xin.transmission_amplitude,
        lambda=Xin.mode.lambda,
        method=Xin.method,
        normalization=Xin.normalization_convention,
        solution_in=Xin,
        solution_up=Xup,
        incidence_values=incidence_values,
        evaluation_points=evaluation_points,
        wronskian_drift=wronskian_drift,
        inputs=inputs,
    )
end

function _scaled_wronskian_incidence(pair)
    absolute_incidence = abs(pair.incidence)
    scale = max(one(absolute_incidence), abs(pair.reflection),
        abs(pair.transmission), maximum(abs, pair.incidence_values))
    return absolute_incidence / scale
end

function _replace_extremal_damped_root_frequency(
        result::LeaverResult, omega, polish_metadata)
    provenance = result.provenance
    angular_order = hasproperty(provenance, :angular_order) ?
        provenance.angular_order : 80
    sheet_id = hasproperty(provenance, :angular_sheet) ?
        Symbol(provenance.angular_sheet, :_gsn_endpoint_polish) :
        :extremal_gsn_endpoint_polish
    angular = angular_branch(result.mode, result.a, omega;
        truncation_order=angular_order, sheet_id)
    angular.status in (
        :continued, :predictor_corrected, :spherical_anchor,
        :high_precision_refined) || return result
    T = typeof(result.a)
    return LeaverResult(
        result.mode, result.convention, result.overtone_index,
        result.inversion_index, result.a, omega, angular.angular_A,
        angular.lambda, angular.mixing, result.cf_value, result.cf_error,
        result.cf_iterations, result.root_residual, T(angular.residual),
        result.precision_bits, result.status, result.stop_reason,
        merge(provenance, polish_metadata, (
            angular_route=angular.metadata.route,
            angular_sheet=angular.sheet_id,
            angular_order=angular.truncation_order,
        )),
    )
end

# At a = +-1, check (and, for roots not from the exact-extremal recurrence,
# Newton-correct) a damped root against the GSN Wronskian incidence.
function _polish_extremal_damped_root_with_gsn(result::LeaverResult;
        incidence_tolerance=1.0e-8,
        wronskian_tolerance=1.0e-8,
        radial_solver=GSN_radial)
    _is_extremal_damped_simple_pole(result) || return result
    T = typeof(real(result.omega))
    tolerance = T(incidence_tolerance)
    wronskian_gate = T(wronskian_tolerance)
    initial = _wronskian_incidence(
        result, result.omega; radial_solver)
    initial_scaled = T(_scaled_wronskian_incidence(initial))
    if get(result.provenance, :root_equation, :none) ==
            :richartz_extremal_recurrence
        # The exact-extremal recurrence fixes the frequency; the radial
        # incidence is only a closure check of the radial basis.
        return _replace_extremal_damped_root_frequency(
            result, result.omega, (
                exact_gsn_root_polish=initial_scaled <= tolerance &&
                    initial.wronskian_drift <= wronskian_gate ?
                    :not_needed : :radial_closure_failed,
                exact_gsn_initial_scaled_incidence=initial_scaled,
                exact_gsn_selected_scaled_incidence=initial_scaled,
                exact_gsn_root_shift=zero(T),
                exact_gsn_polish_rows=NamedTuple[],
            ))
    end
    if initial_scaled <= tolerance &&
            initial.wronskian_drift <= wronskian_gate
        return _replace_extremal_damped_root_frequency(
            result, result.omega, (
                exact_gsn_root_polish=:not_needed,
                exact_gsn_initial_scaled_incidence=initial_scaled,
                exact_gsn_selected_scaled_incidence=initial_scaled,
                exact_gsn_root_shift=zero(T),
                exact_gsn_polish_rows=NamedTuple[],
            ))
    end

    base_step = T(1.0e-4) * max(one(T), abs(result.omega))
    best_omega = result.omega
    best_scaled = initial_scaled
    rows = NamedTuple[]
    for factor in T[one(T), T(2), T(0.5)]
        step = base_step * factor
        plus = _wronskian_incidence(
            result, result.omega + step; radial_solver)
        minus = _wronskian_incidence(
            result, result.omega - step; radial_solver)
        derivative = (plus.incidence - minus.incidence) / (2step)
        finite_derivative = isfinite(derivative) &&
            abs(derivative) > sqrt(eps(T))
        candidate_omega = finite_derivative ?
            result.omega - initial.incidence / derivative : result.omega
        root_shift = abs(candidate_omega - result.omega) /
            max(one(T), abs(result.omega))
        candidate = if finite_derivative && root_shift <= T(1.0e-5)
            _wronskian_incidence(
                result, candidate_omega; radial_solver)
        else
            initial
        end
        candidate_scaled = T(_scaled_wronskian_incidence(candidate))
        candidate_healthy = finite_derivative &&
            root_shift <= T(1.0e-5) &&
            candidate.wronskian_drift <= wronskian_gate &&
            isfinite(candidate_scaled)
        push!(rows, (
            factor,
            step,
            derivative,
            candidate_omega,
            root_shift=T(root_shift),
            scaled_incidence=candidate_scaled,
            wronskian_drift=T(candidate.wronskian_drift),
            healthy=candidate_healthy,
        ))
        if candidate_healthy && candidate_scaled < best_scaled
            best_omega = candidate_omega
            best_scaled = candidate_scaled
        end
        best_scaled <= tolerance && break
    end
    status = best_scaled <= tolerance ? :accepted : :failed
    return _replace_extremal_damped_root_frequency(
        result, best_omega, (
            exact_gsn_root_polish=status,
            exact_gsn_initial_scaled_incidence=initial_scaled,
            exact_gsn_selected_scaled_incidence=best_scaled,
            exact_gsn_root_shift=T(abs(best_omega - result.omega) /
                max(one(T), abs(result.omega))),
            exact_gsn_incidence_tolerance=tolerance,
            exact_gsn_wronskian_tolerance=wronskian_gate,
            exact_gsn_polish_rows=rows,
        ))
end

"""
    validate_qnm_with_isem(result::LeaverResult; kwargs...)

Validate an accepted or estimated Leaver root with the independent public call
`GSN_radial(s,l,m,a,omega,IN)` (default `method="auto"`, which uses GSN-ISEM).
At exact extremality, and for conditioned roots with `abs(a) >= 0.9999` and
`n >= 4`, the incidence amplitude comes from the IN/UP Wronskian instead. The
scaled gate is `abs(Ain)/max(abs(Aout),abs(Atrans),1) <= incidence_tolerance`;
this result never feeds back into the Leaver root finder.
"""
function validate_qnm_with_isem(result::LeaverResult;
        incidence_tolerance=1.0e-8,
        wronskian_tolerance=1.0e-8,
        radial_solver=GSN_radial)
    result.status in (:accepted, :estimated) || return ISEMValidationResult(
        result.omega, zero(result.omega), zero(result.omega),
        zero(result.omega), oftype(real(result.omega), Inf),
        oftype(real(result.omega), Inf), :not_run, :failed,
        :leaver_root_not_accepted, (leaver_status=result.status,))
    use_wronskian = _use_wronskian_incidence(result)
    if use_wronskian
        pair = _wronskian_incidence(result; radial_solver)
        inputs = pair.inputs
        solution = pair.solution_in
        incidence = pair.incidence
        reflection = pair.reflection
        transmission = pair.transmission
    else
        inputs = _radial_evaluation_inputs(result)
        solution = radial_solver(
            result.mode.s, result.mode.l, result.mode.m,
            inputs.a, inputs.omega, IN)
        incidence = solution.incidence_amplitude
        reflection = solution.reflection_amplitude
        transmission = solution.transmission_amplitude
    end
    absolute_incidence = abs(incidence)
    scale = max(one(absolute_incidence), abs(reflection), abs(transmission))
    scaled_incidence = absolute_incidence / scale
    tolerance = oftype(scaled_incidence, incidence_tolerance)
    wronskian_tolerance_value = oftype(
        scaled_incidence, wronskian_tolerance)
    wronskian_accepted = !use_wronskian ||
        pair.wronskian_drift <= wronskian_tolerance_value
    accepted = isfinite(scaled_incidence) &&
        scaled_incidence <= tolerance && wronskian_accepted
    return ISEMValidationResult(
        result.omega, incidence, reflection, transmission,
        absolute_incidence, scaled_incidence, Symbol(solution.method),
        accepted ? :accepted : :failed,
        accepted ? :accepted :
            !wronskian_accepted ? :wronskian_constancy_gate :
            :scaled_incidence_gate,
        merge(_radial_metadata(solution), (
            incidence_tolerance=tolerance,
            wronskian_tolerance=wronskian_tolerance_value,
            scale_definition=:max_abs_reflection_abs_transmission_one,
            route=use_wronskian ?
                :public_default_GSN_radial_IN_UP_wronskian :
                :public_default_GSN_radial_IN,
            incidence_backend=use_wronskian ?
                :direct_real_axis_wronskian : :direct_asymptotic_amplitude,
            wronskian_evaluation_points=use_wronskian ?
                pair.evaluation_points : (),
            wronskian_incidence_values=use_wronskian ?
                pair.incidence_values : (),
            wronskian_drift=use_wronskian ?
                pair.wronskian_drift : zero(scaled_incidence),
            evaluation_policy=inputs.policy,
            frequency_projection_drift=inputs.projection_drift,
            source_precision_bits=inputs.source_precision_bits,
            evaluation_precision_bits=inputs.evaluation_precision_bits,
            evaluation_a=inputs.a,
            evaluation_omega=inputs.omega,
            leaver_root_status=result.status,
            leaver_root_stop_reason=result.stop_reason,
        )),
    )
end
