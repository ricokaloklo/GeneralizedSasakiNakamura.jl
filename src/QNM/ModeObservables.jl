function _solve_qnm_root(mode, a, guess, root_options;
        convention=:overtone, inversion_index=:auto, segment=0)
    options = merge(root_options,
        (; guess, convention, inversion_index, segment))
    return qnm_frequency(mode, a; options...)
end

_qnm_root_segment(root) = hasproperty(root, :provenance) ?
    get(root.provenance, :sequence_segment, 0) : 0

function _qnm_radial_solutions(root)
    inputs = _radial_evaluation_inputs(root)
    mode = root.mode
    radial_boundary = IN
    if abs(inputs.a) == one(inputs.a)
        X = GSN_radial(
            mode.s, mode.l, mode.m, inputs.a, inputs.omega, IN;
            method="GSN-ISEM")
        R = Teukolsky_radial(
            mode.s, mode.l, mode.m, inputs.a, inputs.omega, IN;
            method="GSN-ISEM")
        y_boundary = mode.s == -2 ? IN : mode.s == 2 ? UP : nothing
        Y = y_boundary === nothing ? missing : Y_radial(
            mode.s, mode.l, mode.m, inputs.a, inputs.omega, y_boundary;
            method="GSN-ISEM")
        return (
            X=X,
            Y=Y,
            R=R,
            boundary_condition=radial_boundary,
            Y_boundary_condition=y_boundary,
            evaluation_a=inputs.a,
            evaluation_omega=inputs.omega,
            evaluation_policy=inputs.policy,
            frequency_projection_drift=inputs.projection_drift,
            source_precision_bits=inputs.source_precision_bits,
            evaluation_precision_bits=inputs.evaluation_precision_bits,
            direct_route=missing,
            Y_direct_route=missing,
            pole_normalization=:exact_extremal_irregular_horizon_unit_transmission,
            coordinate_conventions=(X=:rstar, Y=:boyer_lindquist_r,
                R=:boyer_lindquist_r),
        )
    end
    direct_route = ISEM.DirectGSN.direct_complex_route(
        mode.s, mode.l, mode.m, inputs.a, inputs.omega, :IN;
        pole_normalization=true)
    X = ISEM.DirectGSN.direct_complex_gsn_radial_function(direct_route)
    R = ISEM.DirectGSN.direct_teukolsky_radial_function(direct_route)

    y_boundary = mode.s == -2 ? IN : mode.s == 2 ? UP : nothing
    y_route = if y_boundary === nothing
        nothing
    elseif y_boundary == IN
        direct_route
    else
        ISEM.DirectGSN.direct_complex_route(
            mode.s, mode.l, mode.m, inputs.a, inputs.omega, :UP;
            pole_normalization=true)
    end
    Y = y_route === nothing ? missing :
        ISEM.DirectGSN.direct_y_radial_function(y_route)
    return (
        X=X,
        Y=Y,
        R=R,
        boundary_condition=radial_boundary,
        Y_boundary_condition=y_boundary,
        evaluation_a=inputs.a,
        evaluation_omega=inputs.omega,
        evaluation_policy=inputs.policy,
        frequency_projection_drift=inputs.projection_drift,
        source_precision_bits=inputs.source_precision_bits,
        evaluation_precision_bits=inputs.evaluation_precision_bits,
        direct_route=direct_route,
        Y_direct_route=y_route,
        pole_normalization=:exact_incidence_zero,
        coordinate_conventions=(X=:rstar, Y=:boyer_lindquist_r,
            R=:boyer_lindquist_r),
    )
end

function _finite_qnm_observables(validation, excitation)
    return all(isfinite, (
        validation.incidence,
        validation.reflection,
        validation.transmission,
        excitation.alpha,
        excitation.reflection,
        excitation.B_gsn,
    ))
end

function _estimated_decimal_digits(relative_drift)
    T = typeof(relative_drift)
    iszero(relative_drift) && return T(Inf)
    isfinite(relative_drift) || return zero(T)
    return max(zero(T), -log10(relative_drift))
end

function _qnm_branch_observables(root, excitation_options;
        input_a=root.a, detailed::Bool=false)
    if !(root.status in (:accepted, :estimated))
        return QNMFailure((
            a=input_a,
            mode=root.mode,
            frequency=root.omega,
            omega=root.omega,
            convention=root.convention,
            overtone_index=root.overtone_index,
            inversion_index=root.inversion_index,
            segment=_qnm_root_segment(root),
            root=root,
            validation=missing,
            excitation=missing,
            stage=:root,
            status=:failed,
            stop_reason=root.stop_reason,
            formalism=:GSN,
        ))
    end
    if get(root.provenance, :extremal_endpoint, false) &&
            !get(root.provenance, :simple_pole, true)
        return QNMEndpointResult((
            a=input_a,
            mode=root.mode,
            frequency=root.omega,
            omega=root.omega,
            angular_A=root.angular_A,
            lambda=root.lambda,
            root=root,
            convention=root.convention,
            overtone_index=root.overtone_index,
            inversion_index=root.inversion_index,
            segment=_qnm_root_segment(root),
            spectrum_object=root.provenance.spectrum_object,
            radial_observables=root.provenance.radial_observables,
            endpoint_method=root.provenance.endpoint_method,
            simple_pole=false,
            coalesced_overtone_labels=get(
                root.provenance, :coalesced_overtone_labels, false),
            requested_overtone_label=root.mode.n,
            overtone_interpretation=
                :near_extremal_family_label_not_distinct_endpoint_pole,
            incidence_amplitude=missing,
            reflection_amplitude=missing,
            incidence_derivative=missing,
            excitation_factor=missing,
            X=missing,
            Y=missing,
            R=missing,
            status=:accepted,
            stop_reason=root.stop_reason,
            formalism=:GSN_Teukolsky_endpoint,
        ))
    end
    if !get(root.provenance, :gsn_amplitude_limit_implemented, true)
        return QNMFailure((
            a=input_a,
            mode=root.mode,
            frequency=root.omega,
            omega=root.omega,
            angular_A=root.angular_A,
            lambda=root.lambda,
            convention=root.convention,
            overtone_index=root.overtone_index,
            inversion_index=root.inversion_index,
            segment=_qnm_root_segment(root),
            root=root,
            validation=missing,
            excitation=missing,
            stage=:pole_skipping_endpoint,
            status=:failed,
            stop_reason=:matsubara_qnm_pole_skipping_no_simple_gsn_residue,
            criterion=:gsn_incident_amplitude_simple_zero,
            formalism=:GSN,
        ))
    end
    if _is_extremal_damped_simple_pole(root)
        radial_solver = get(excitation_options, :radial_solver,
            _extremal_qnm_radial_solver)
        incidence_tolerance = get(excitation_options, :incidence_tolerance, 1.0e-8)
        wronskian_tolerance = get(excitation_options, :wronskian_tolerance, 1.0e-8)
        root = _polish_extremal_damped_root_with_gsn(root;
            incidence_tolerance, wronskian_tolerance, radial_solver)
        excitation_options = merge(
            excitation_options, (radial_solver=radial_solver,))
    end
    if !hasproperty(excitation_options, :radial_solver)
        excitation_options = merge(
            excitation_options, (radial_solver=_qnm_default_radial_solver,))
    end
    excitation = try
        detailed ? qnm_excitation_factor(root; excitation_options...) :
            _qnm_gsn_excitation_factor(root; excitation_options...)
    catch error
        return QNMFailure((
            a=input_a,
            mode=root.mode,
            frequency=root.omega,
            omega=root.omega,
            convention=root.convention,
            overtone_index=root.overtone_index,
            inversion_index=root.inversion_index,
            segment=_qnm_root_segment(root),
            root=root,
            validation=missing,
            excitation=missing,
            stage=:radial_evaluation,
            status=:failed,
            stop_reason=:radial_evaluation_error,
            error_type=Symbol(nameof(typeof(error))),
            error_message=sprint(showerror, error),
            formalism=:GSN,
        ))
    end
    validation = excitation.metadata.validation
    conditioned_near_extremal = validation.status == :accepted &&
        hasproperty(validation.metadata,
            :conditioned_near_extremal_validation) &&
        validation.metadata.conditioned_near_extremal_validation
    high_damping_diagnostic =
        _finite_high_damping_diagnostic_candidate(root, validation) &&
        hasproperty(excitation.metadata, :high_damping_diagnostic) &&
        excitation.metadata.high_damping_diagnostic
    accepted = root.status == :accepted &&
        validation.status == :accepted && excitation.status == :accepted &&
        !conditioned_near_extremal
    stop_reason = root.status != :accepted ? root.stop_reason :
        validation.status != :accepted ? validation.stop_reason :
        conditioned_near_extremal ?
            :near_extremal_scaled_incidence_estimate :
            excitation.stop_reason
    finite_observables = _finite_qnm_observables(validation, excitation)
    estimated = !accepted && root.status in (:accepted, :estimated) &&
        finite_observables &&
        ((validation.status == :accepted &&
          (excitation.stop_reason == :richardson_step_gate ||
           conditioned_near_extremal || root.status == :estimated)) ||
         (high_damping_diagnostic &&
          excitation.stop_reason == :scaled_incidence_gate))
    status = accepted ? :accepted : estimated ? :estimated : :failed
    achieved_relative_precision = excitation.step_drift
    target_relative_precision = hasproperty(
        excitation.metadata, :step_tolerance) ?
        excitation.metadata.step_tolerance :
        oftype(achieved_relative_precision, 1.0e-7)
    limiting_gate = !estimated ? :none : high_damping_diagnostic ?
        :scaled_incidence_gate : root.status == :estimated ?
        :conditioned_root_gate : conditioned_near_extremal ?
        excitation.stop_reason == :richardson_step_gate ?
            :richardson_step_and_scaled_incidence_gates :
            :scaled_incidence_gate : :richardson_step_gate
    gsn = (
        incidence_amplitude=validation.incidence,
        reflection_amplitude=excitation.reflection,
        transmission_amplitude=validation.transmission,
        incidence_derivative=excitation.alpha,
        excitation_factor=excitation.B_gsn,
        normalization=:gsn_unit_transmission,
    )
    compact = (
        a=input_a,
        mode=root.mode,
        frequency=root.omega,
        omega=root.omega,
        angular_A=root.angular_A,
        lambda=root.lambda,
        incidence_amplitude=gsn.incidence_amplitude,
        reflection_amplitude=gsn.reflection_amplitude,
        incidence_derivative=gsn.incidence_derivative,
        excitation_factor=gsn.excitation_factor,
        gsn=gsn,
        root=root,
        convention=root.convention,
        overtone_index=root.overtone_index,
        inversion_index=root.inversion_index,
        segment=_qnm_root_segment(root),
        validation=validation,
        excitation=excitation,
        status=status,
        stop_reason=stop_reason,
        scientific_acceptance=accepted,
        achieved_relative_precision=achieved_relative_precision,
        target_relative_precision=target_relative_precision,
        scaled_incidence_residual=validation.scaled_incidence,
        incidence_target=validation.metadata.incidence_tolerance,
        inferred_root_shift=hasproperty(
            validation.metadata, :inferred_root_shift) ?
            validation.metadata.inferred_root_shift : missing,
        estimated_decimal_digits=
            _estimated_decimal_digits(achieved_relative_precision),
        limiting_gate=limiting_gate,
        precision_basis=
            high_damping_diagnostic ?
            :richardson_drift_with_unclosed_scaled_incidence_not_error_bound :
            conditioned_near_extremal ?
            :richardson_drift_and_scaled_incidence_indicators_not_error_bounds :
            :richardson_step_drift_convergence_indicator_not_error_bound,
        primary_convention=:gsn,
        formalism=:GSN,
        detailed=false,
    )
    if status == :failed
        stage = validation.status != :accepted ?
            :GSN_ISEM_validation : :excitation_factor
        return QNMFailure(merge(compact, (stage=stage,)))
    end
    ResultType = status == :estimated ? QNMEstimate : QNMResult
    detailed || return ResultType(compact)

    solutions = _qnm_radial_solutions(root)
    teukolsky = (
        reflection_amplitude=excitation.reflection_teukolsky,
        incidence_derivative=excitation.alpha_teukolsky,
        excitation_factor=excitation.B_teukolsky,
        normalization=:teukolsky_unit_transmission,
    )
    return ResultType(merge(compact, (
        X=solutions.X,
        Y=solutions.Y,
        R=solutions.R,
        solutions=solutions,
        teukolsky=teukolsky,
        detailed=true,
    )))
end

function _qnm_selector(branch::QNMBranch)
    return branch
end

function _qnm_selector(branch::Union{Symbol,AbstractString})
    normalized = lowercase(strip(String(branch)))
    normalized in ("ordinary", "positive_real") && return Ordinary
    normalized in ("mirror", "negative_real") && return Mirror
    throw(ArgumentError(
        "branch must be Ordinary or Mirror, not $(repr(branch))."))
end

function _ordinary_qnm_root(a, s, l, m, n, primary_guess, root_options;
        convention=:overtone, inversion_index=:auto, segment=0)
    mode = QNMMode(Int(s), Int(l), Int(m), Int(n), :positive_real)
    return _solve_qnm_root(mode, a, primary_guess, root_options;
        convention, inversion_index, segment)
end

function _mirror_qnm_root(a, s, l, m, n, mirror_guess, root_options;
        ordinary_root=nothing, convention=:overtone,
        inversion_index=:auto, segment=0)
    mode = QNMMode(Int(s), Int(l), Int(m), Int(n), :negative_real)
    if abs(a) == one(a) && _extremal_zdm_endpoint_scope(mode, a)
        root = _solve_qnm_root(
            mode, a, mirror_guess, root_options;
            convention, inversion_index, segment)
        return root, root.omega, nothing
    end
    endpoint_contract = _schwarzschild_label_contract(
        mode, _qnm_convention(convention))
    if mirror_guess === nothing && iszero(a) &&
            endpoint_contract.near_as_role != :none
        root = _solve_qnm_root(mode, a, nothing, root_options;
            convention, inversion_index, segment)
        return root, root.omega, nothing
    end
    partner_root = if mirror_guess === nothing
        if m == 0 && ordinary_root !== nothing
            ordinary_root
        else
            partner_mode = QNMMode(
                Int(s), Int(l), -Int(m), Int(n), :positive_real)
            _solve_qnm_root(partner_mode, a, nothing, root_options;
                convention, inversion_index, segment)
        end
    else
        nothing
    end
    # The mirror of segment k is -conj of the ordinary segment k; it inherits
    # the partner's segment label.
    mirror_segment = partner_root !== nothing ?
        _qnm_root_segment(partner_root) : segment === :auto ? 0 : segment
    function solve_mirror()
        seed = mirror_guess === nothing ?
            -conj(partner_root.omega) : mirror_guess
        root = _solve_qnm_root(mode, a, seed, root_options;
            convention, inversion_index, segment=mirror_segment)
        return root, seed, partner_root
    end
    # Keep an automatic high-precision fallback at its actual precision;
    # unary arithmetic must not silently promote it to the global 256 bits.
    if partner_root !== nothing && partner_root.omega isa Complex{BigFloat}
        return setprecision(solve_mirror, BigFloat, partner_root.precision_bits)
    end
    return solve_mirror()
end

"""
    qnm(a, s, l, m, n, branch=ordinary; kwargs...)

Compute one complete QNM branch in `M=1` units using five numerical labels.
The optional sixth positional argument is `ordinary` or `mirror`; it defaults
to `ordinary`. The default result computes only the frequency, angular
eigenvalue, directly evaluated GSN scattering amplitudes, real-direction
Richardson incidence derivative, and GSN excitation factor shown by its
compact display.

Set `detailed=true` to additionally compute the Teukolsky excitation-factor
bridge and callable radial solution objects `X`, `Y`, and `R`.

`X(rstar)` evaluates the GSN solution, while `Y(r)` and `R(r)` use the
Boyer-Lindquist radius. `X` and `R` use the `IN` QNM basis. The direct public
`Y_radial` transformation supports the matching `IN` branch for `s=-2`; for
`s=2` it returns the supported `UP` branch, and for other spins `Y` is
`missing`. The exact boundary and coordinate conventions are recorded in
`branch.solutions`.

For `mirror`, the seed is `-conj(omega[s,l,-m,n])`, but the returned root is
re-polished on the negative-real branch. Amplitudes and frequency derivatives
are evaluated at that branch; the radial solver may use exact conjugation
symmetry rather than solve the reflected equation again.

`segment` selects the segment of a sequence that ends on the negative
imaginary axis (see [`qnm_frequency`](@ref)); the mirror of segment `k` is
seeded from the ordinary segment `k`.

Use [`qnm_pair`](@ref) only when both branches are required in one result.
`root_options` and `excitation_options` are optional NamedTuples forwarded to
`qnm_frequency` and `qnm_excitation_factor`. The no-keyword call uses package
defaults. If the frequency, incidence-zero, normalization, and finiteness
checks pass but only the Richardson step-drift target is missed, the call
returns `QNMEstimate` with the finite observables and measured convergence
indicator. It does not relax the production acceptance threshold.
"""
function qnm(a, s::Integer, l::Integer, m::Integer, n::Integer,
        branch::QNMBranch=ordinary;
        convention=:overtone,
        inversion_index=:auto,
        primary_guess=nothing,
        mirror_guess=nothing,
        segment=0,
        root_options::NamedTuple=NamedTuple(),
        excitation_options::NamedTuple=NamedTuple(),
        detailed::Bool=false)
    root = if branch == ordinary
        _ordinary_qnm_root(
            a, s, l, m, n, primary_guess, root_options;
            convention, inversion_index, segment)
    else
        first(_mirror_qnm_root(
            a, s, l, m, n, mirror_guess, root_options;
            convention, inversion_index, segment))
    end
    return _qnm_branch_observables(
        root, excitation_options; input_a=a, detailed)
end

function qnm(a, s::Integer, l::Integer, m::Integer, n::Integer,
        branch::Union{Symbol,AbstractString}; kwargs...)
    return qnm(a, s, l, m, n, _qnm_selector(branch); kwargs...)
end

"""
    qnm_pair(a, s, l, m, n; kwargs...)

Compute and return both complete QNM branches. Prefer [`qnm`](@ref) when only
one branch is needed so the other branch is not evaluated.
"""
function qnm_pair(a, s::Integer, l::Integer, m::Integer, n::Integer;
        convention=:overtone,
        inversion_index=:auto,
        primary_guess=nothing,
        mirror_guess=nothing,
        segment=0,
        root_options::NamedTuple=NamedTuple(),
        excitation_options::NamedTuple=NamedTuple(),
        detailed::Bool=false)
    ordinary_root = _ordinary_qnm_root(
        a, s, l, m, n, primary_guess, root_options;
        convention, inversion_index, segment)
    mirror_root, seed, partner_root = _mirror_qnm_root(
        a, s, l, m, n, mirror_guess, root_options;
        ordinary_root, convention, inversion_index, segment)

    ordinary = _qnm_branch_observables(
        ordinary_root, excitation_options; input_a=a, detailed)
    mirror = _qnm_branch_observables(
        mirror_root, excitation_options; input_a=a, detailed)
    branch_failed = ordinary.status == :failed || mirror.status == :failed
    branch_estimated = ordinary.status == :estimated ||
        mirror.status == :estimated
    status = branch_failed ? :failed : branch_estimated ?
        :estimated : :accepted
    stop_reason = ordinary.status != :accepted ?
        Symbol("ordinary_", ordinary.stop_reason) :
        mirror.status != :accepted ?
        Symbol("mirror_", mirror.stop_reason) : :accepted
    partner_summary = partner_root === nothing ? nothing : (
        mode=partner_root.mode,
        frequency=partner_root.omega,
        status=partner_root.status,
        root_residual=partner_root.root_residual,
    )
    return QNMPairResult((
        a=a,
        s=Int(s),
        l=Int(l),
        m=Int(m),
        n=Int(n),
        ordinary=ordinary,
        mirror=mirror,
        status=status,
        stop_reason=stop_reason,
        scientific_acceptance=status == :accepted,
        metadata=(
            mirror_seed=seed,
            mirror_seed_policy=mirror_guess === nothing ?
                :homogeneous_spectral_symmetry_direct_repolish :
                :explicit_mirror_guess_direct_repolish,
            mirror_partner=partner_summary,
            amplitude_fill_policy=:direct_no_conjugate_fill,
            convention=ordinary_root.convention,
            overtone_index=ordinary_root.overtone_index,
            inversion_index=ordinary_root.inversion_index,
            segment=_qnm_root_segment(ordinary_root),
            units=:M_equals_one,
        ),
    ))
end
