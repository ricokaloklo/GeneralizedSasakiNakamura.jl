function _solve_qnm_root(mode, a, guess, root_options;
        convention=:overtone, inversion_index=:auto, segment=0)
    options = merge(root_options,
        (; guess, convention, inversion_index, segment))
    return qnm_frequency(mode, a; options...)
end

_qnm_root_segment(root) = hasproperty(root, :provenance) ?
    get(root.provenance, :sequence_segment, 0) : 0

# Radial solutions X, Y, R at an accepted root with unit-transmission pole normalization: the public ISEM
# solvers at exact extremality, otherwise the complex-frequency ISEM route.
function _qnm_radial_solutions(root)
    inputs = _radial_evaluation_inputs(root)
    mode = root.mode
    radial_boundary = IN
    if abs(inputs.a) == one(inputs.a)
        X = GSN_radial(
            mode.s, mode.l, mode.m, inputs.a, inputs.omega, IN;
            method="ISEM")
        R = Teukolsky_radial(
            mode.s, mode.l, mode.m, inputs.a, inputs.omega, IN;
            method="ISEM")
        y_boundary = mode.s == -2 ? IN : mode.s == 2 ? UP : nothing
        Y = y_boundary === nothing ? missing : Y_radial(
            mode.s, mode.l, mode.m, inputs.a, inputs.omega, y_boundary;
            method="ISEM")
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
            isem_route=missing,
            Y_isem_route=missing,
            pole_normalization=:exact_extremal_irregular_horizon_unit_transmission,
            coordinate_conventions=(X=:rstar, Y=:boyer_lindquist_r,
                R=:boyer_lindquist_r),
        )
    end
    isem_route = ISEM.complex_route(
        mode.s, mode.l, mode.m, inputs.a, inputs.omega, :IN;
        pole_normalization=true)
    X = ISEM.complex_gsn_radial_function(isem_route)
    R = ISEM.teukolsky_radial_function(isem_route)

    y_boundary = mode.s == -2 ? IN : mode.s == 2 ? UP : nothing
    y_route = if y_boundary === nothing
        nothing
    elseif y_boundary == IN
        isem_route
    else
        ISEM.complex_route(
            mode.s, mode.l, mode.m, inputs.a, inputs.omega, :UP;
            pole_normalization=true)
    end
    Y = y_route === nothing ? missing :
        ISEM.y_radial_function(y_route)
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
        isem_route=isem_route,
        Y_isem_route=y_route,
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

# Turn one `LeaverResult` into the public result: validate it with ISEM,
# compute the excitation factor, and wrap it as QNMResult, QNMEstimate,
# QNMEndpointResult or QNMFailure (with radial solutions when `detailed`).
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
    if _is_extremal_damped_simple_pole(root) &&
            !_qnm_use_variational(root, excitation_options)
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
            hasproperty(excitation.metadata, :precision_basis) ?
            excitation.metadata.precision_basis : high_damping_diagnostic ?
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
            :ISEM_validation : :excitation_factor
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

# Solve the negative-real (mirror) root. Without `mirror_guess` it is the reflection -conj(omega) of the
# accepted (s, l, -m, n) ordinary root at the same spin (or of the supplied `ordinary_root` for m = 0), returned as
# the partner; a = 0 near-algebraically-special modes are solved directly. With `mirror_guess`, or when the partner
# is not accepted, the root is polished on the negative-real branch from the guess or from -conj of the partner.
# Returns (root, seed, partner_root).
function _mirror_qnm_root(a, s, l, m, n, mirror_guess, root_options;
        ordinary_root=nothing, convention=:overtone,
        inversion_index=:auto, segment=0)
    mode = QNMMode(Int(s), Int(l), Int(m), Int(n), :negative_real)
    endpoint_contract = _schwarzschild_label_contract(
        mode, _qnm_convention(convention))
    if mirror_guess === nothing && iszero(a) &&
            endpoint_contract.near_as_role != :none
        root = _solve_qnm_root(mode, a, nothing, root_options;
            convention, inversion_index, segment)
        return root, root.omega, nothing
    end
    partner_root = nothing
    if mirror_guess === nothing
        partner_root = m == 0 && ordinary_root !== nothing ? ordinary_root :
            _solve_qnm_root(QNMMode(Int(s), Int(l), -Int(m), Int(n), :positive_real),
                a, nothing, root_options; convention, inversion_index, segment)
        # A partner sequence that ended on the negative imaginary axis ends the mirror
        # sequence at the same spin: its reflection is returned with the same status.
        if partner_root.status == :accepted ||
                partner_root.stop_reason == :negative_imaginary_axis_endpoint
            reflect() = _mirror_reflected_result(mode, partner_root;
                angular_order=get(root_options, :angular_order, 32),
                sheet_id=get(root_options, :sheet_id, :qnm_straight_from_spherical))
            root = partner_root.omega isa Complex{BigFloat} ?
                setprecision(reflect, BigFloat, partner_root.precision_bits) : reflect()
            return root, root.omega, partner_root
        end
    end
    # Explicit guess, or a partner that was not accepted: polish on the negative-real
    # branch from the guess or from the partner's reflected frequency.
    # A reflected partner frequency on the positive-real side of the axis is moved to
    # Re(omega) = -max(|Re|, 1e-3) so that the polish starts on the mirror branch.
    seed = mirror_guess !== nothing ? mirror_guess :
        let reflected = -conj(partner_root.omega)
            _branch_consistent(mode, reflected) ? reflected :
                complex(-max(abs(real(reflected)), oftype(real(reflected), 1.0e-3)),
                    imag(reflected))
        end
    mirror_segment = partner_root !== nothing ? _qnm_root_segment(partner_root) :
        segment === :auto ? 0 : segment
    solve_mirror() = _solve_qnm_root(mode, a, seed, root_options;
        convention, inversion_index, segment=mirror_segment)
    root = partner_root !== nothing && partner_root.omega isa Complex{BigFloat} ?
        setprecision(solve_mirror, BigFloat, partner_root.precision_bits) :
        solve_mirror()
    return root, seed, partner_root
end

"""
    qnm(a, s, l, m, n, branch=ordinary; kwargs...)

Compute one QNM in `M=1` units: its frequency, angular eigenvalue, GSN
incidence and reflection amplitudes, the frequency derivative of the incidence
amplitude and the GSN excitation factor. The optional sixth positional argument
is `ordinary` (default) or `mirror`.

For `s = -2, …, 2`, `l ≤ 15` and `n ≤ 10` the frequency derivative comes from a
single radial solve that also integrates the derivative with respect to `ω`;
when its accuracy checks fail, and for all other modes, it is
Richardson-extrapolated from frequency stencils. Set
`excitation_options=(backend=:richardson,)` to always use the stencils.

Set `detailed=true` to also compute the excitation factor in the Teukolsky
normalization and the callable radial solutions `X`, `Y` and `R`.

`X(rstar)` evaluates the GSN solution, while `Y(r)` and `R(r)` use the
Boyer-Lindquist radius. `X` and `R` use the `IN` QNM basis. The direct public
`Y_radial` transformation supports the matching `IN` branch for `s=-2`; for
`s=2` it returns the supported `UP` branch, and for other spins `Y` is
`missing`. The exact boundary and coordinate conventions are recorded in
`branch.solutions`.

For `mirror`, the root is the exact reflection `-conj(ω)` of the ordinary
`(s, l, -m, n)` root at the same spin, with the angular data recomputed; with
`mirror_guess` it is polished on the negative-real branch from that guess. The
amplitudes and the frequency derivative are evaluated at the mirror root.

`convention` (`:overtone`, `:complete_spectrum` or `:leaver`), `inversion_index`
and `segment` (`0`, `1` or `:auto`) have the meaning described in
[`qnm_frequency`](@ref); `segment` selects the segment of a sequence that ends
on the negative imaginary axis, and the mirror of segment `k` is seeded from
the ordinary segment `k`. `primary_guess` and `mirror_guess` are optional
explicit frequency guesses for the `ordinary` and `mirror` roots
(`convention=:leaver` requires one).

Use [`qnm_pair`](@ref) only when both branches are required in one result.
`root_options` and `excitation_options` are optional NamedTuples forwarded to
`qnm_frequency` and `qnm_excitation_factor`. The no-keyword call uses package
defaults. The result is a `QNMResult` when every check passes, a
`QNMEstimate` when the observables are finite but one check misses its target
(named in `limiting_gate`), a `QNMFailure` when a step fails, and a
`QNMEndpointResult` for a zero-damping mode at `a = ±1`.
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
    frequency_start = time_ns()
    root = if branch == ordinary
        _ordinary_qnm_root(
            a, s, l, m, n, primary_guess, root_options;
            convention, inversion_index, segment)
    else
        first(_mirror_qnm_root(
            a, s, l, m, n, mirror_guess, root_options;
            convention, inversion_index, segment))
    end
    frequency_seconds = (time_ns() - frequency_start) / 1e9
    excitation_start = time_ns()
    result = _qnm_branch_observables(
        root, excitation_options; input_a=a, detailed)
    excitation_seconds = (time_ns() - excitation_start) / 1e9
    ResultType = result isa QNMResult ? QNMResult :
        result isa QNMEstimate ? QNMEstimate :
        result isa QNMEndpointResult ? QNMEndpointResult : QNMFailure
    return ResultType(merge(result.data,
        (timing=(frequency_seconds=frequency_seconds,
            excitation_seconds=excitation_seconds,
            includes_detailed_radial_objects=detailed),)))
end

function qnm(a, s::Integer, l::Integer, m::Integer, n::Integer,
        branch::Union{Symbol,AbstractString}; kwargs...)
    return qnm(a, s, l, m, n, _qnm_selector(branch); kwargs...)
end

"""
    qnm_pair(a, s, l, m, n; kwargs...)

Compute and return both complete QNM branches as a `QNMPairResult` with fields
`ordinary` and `mirror`. Keywords are those of [`qnm`](@ref). Without
`mirror_guess`, the mirror root is the reflection `-conj(ω)` of the ordinary
`(s, l, -m, n)` root at the same spin. Prefer [`qnm`](@ref) when only
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
                :homogeneous_spectral_symmetry_reflection :
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
