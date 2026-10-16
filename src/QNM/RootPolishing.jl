const _SCHWARZSCHILD_SEEDS = Dict{NTuple{3,Int},ComplexF64}(
    (-2, 2, 0) => 0.37367168441804185 - 0.08896231568893566im,
    (-2, 2, 1) => 0.34671099687916350 - 0.27391487529123500im,
    (-2, 2, 2) => 0.30105345461236650 - 0.47827698322307200im,
    (-2, 2, 3) => 0.25150496218564372 - 0.70514820243348242im,
    (-2, 2, 4) => 0.20751457981301125 - 0.94684489086636203im,
    (-2, 2, 5) => 0.16929940309302113 - 1.19560805413578164im,
    (-2, 2, 6) => 0.13325234024514371 - 1.44791062616211419im,
    (-2, 2, 7) => 0.09282233367019177 - 1.70384117220627562im,
    # n = 8 is the algebraically special frequency -2i, returned exactly. These are the next roots.
    (-2, 2, 9) => 0.063263505125603553 - 2.3026447651585391im,
    (-2, 2, 10) => 0.07655346288598433 - 2.5608266173815060im,
    (-2, 3, 0) => 0.59944328843749000 - 0.09270304794494700im,
    (-2, 3, 1) => 0.58264380303330066 - 0.28129811343504524im,
    (-2, 3, 2) => 0.55168490077844956 - 0.47909275096696286im,
    (-2, 4, 0) => 0.80917837753223900 - 0.09416396098892300im,
    (-2, 4, 1) => 0.79663153203450299 - 0.28433434940484159im,
    (-2, 4, 2) => 0.77270953260670861 - 0.47990817512116291im,
)

const _KERR_L2_M2_NEAR_AS_ANCHOR_SPIN = 0.05
const _KERR_L2_M2_MANO_SUZUKI_TAKASUGI_SWITCH_SPIN = 0.015
const _KERR_L2_M2_NEAR_AS_ANCHORS = Dict{Symbol,ComplexF64}(
    :unconventional => 0.06177258790671765 - 1.8876918643227598im,
    :algebraically_special =>
        0.060299316438925596 - 2.0888402014191692im,
)
# Counterrotating n=8 at a=0.5, Table I of arXiv:2609.09531v1 (M=1).
# This rounded value is only a seed; the continued fraction is solved anew.
const _KERR_L2_COUNTER_N8_ANCHOR = 0.098356 - 2.078506im

# Branch anchors from the Python qnm 0.4.4 (commit f3abd18) for the counterrotating
# (s,l,|m|)=(-2,2,2) family.  They select the physical overtone branch only;
# angular data and the radial continued fraction are recomputed at the target.
const _COUNTERROTATING_NEAR_EXTREMAL_ANCHORS = Dict(
    (0.9999, 4) => 0.10966391100861249 - 1.0484565784875364im,
    (0.9999, 5) => 0.09622818654671349 - 1.3188614163368266im,
    (0.9999, 6) => 0.08552024051137216 - 1.5859422845532654im,
    (0.9999, 7) => 0.07630883605751265 - 1.8541542337264556im,
    (0.9999, 9) => -0.13788097978067113 - 2.2230115092378702im,
    (0.9999, 10) => 0.08172665916503848 - 2.674770803913548im,
    (0.99999, 4) => 0.10967190514278498 - 1.0484624076661069im,
    (0.99999, 5) => 0.09620454978969084 - 1.3188953723824268im,
    (0.99999, 6) => 0.08546785861495902 - 1.5860447304425451im,
    (0.99999, 7) => 0.07632391674711030 - 1.8542317885405424im,
    (0.99999, 9) => -0.14229956759043466 - 2.2191402886345815im,
    (0.99999, 10) => 0.08165609544784828 - 2.6748151109278746im,
)
const _COUNTERROTATING_ANCHOR_SPINS = (0.9999, 0.99999)

function _working_type(a, guess)
    if a isa BigFloat || guess isa Complex{BigFloat} || guess isa BigFloat
        return BigFloat
    end
    return Float64
end

function _complex_of_type(::Type{T}, value) where {T<:AbstractFloat}
    return Complex{T}(T(real(value)), T(imag(value)))
end

# The l=2, |m|=2 complete-spectrum n=8/n=9 trajectories cross the NIA. Their
# label is fixed by continuation from the finite-spin anchor, not by sign(Re w).
_negative_imaginary_axis_crossing_exempt(mode::QNMMode) = mode.branch == :positive_real &&
    mode.s == -2 && mode.l == 2 && abs(mode.m) == 2 && mode.n in (8, 9)

_branch_sign(mode::QNMMode) = mode.branch == :positive_real ? 1 : -1

function _branch_consistent(mode::QNMMode, omega)
    _negative_imaginary_axis_crossing_exempt(mode) && return true
    if mode.branch == :positive_real
        return real(omega) > 0
    end
    return real(omega) < 0
end

# Label metadata for a mode under `convention`: near-algebraically-special role,
# endpoint kind, branch label and the Schwarzschild sequence index to continue from.
function _schwarzschild_label_contract(mode::QNMMode, convention::Symbol, a=0;
        segment::Integer=0)
    # Negative public spins are evaluated through (a,m) -> (-a,-m), so the
    # same physical near-AS branch appears internally with m=-2.
    l2m2 = mode.s == -2 && mode.l == 2 && abs(mode.m) == 2
    kerr_overtone = convention == :overtone && l2m2 && !iszero(a)
    corotating = a * mode.m * (mode.branch == :positive_real ? 1 : -1) > 0
    complete = convention == :complete_spectrum || (kerr_overtone && corotating)
    counter_n8 = kerr_overtone && !corotating && mode.n == 8
    # At a = 0 the spectrum does not depend on m, and s = +2 shares it, so every
    # |s| = 2, l = 2 overtone 8 is the algebraically special frequency -2i.
    standard_as = convention != :leaver && abs(mode.s) == 2 && mode.l == 2 &&
        mode.n == 8 && iszero(a) && !(complete && l2m2)
    complete_unconventional = complete && l2m2 &&
        mode.n == 8
    complete_as = complete && l2m2 && mode.n == 9
    near_as_role = counter_n8 ? :counterrotating :
        standard_as || complete_as ? :algebraically_special :
        complete_unconventional ? :unconventional : :none
    is_algebraically_special_pair = near_as_role != :none
    endpoint_kind = near_as_role == :none ? :regular : near_as_role
    sequence_index = convention == :leaver ? missing :
        complete && l2m2 && mode.n >= 10 ?
            mode.n - 1 :
        complete_as ? missing : mode.n
    branch_label = convention == :leaver ? :unassigned_leaver_root :
        near_as_role == :unconventional ? :kerr_branch_8_0 :
        near_as_role == :algebraically_special ? :kerr_branch_8_1 :
        _segment_branch_label(mode.n, segment)
    return (
        label_convention=convention == :overtone ? :standard_overtone :
            convention == :complete_spectrum ?
                :chen_complete_spectrum : :leaver_inversion,
        branch_label,
        endpoint_kind,
        near_as_role,
        algebraically_special_pair=is_algebraically_special_pair,
        schwarzschild_sequence_index=sequence_index,
        schwarzschild_multiplicity=1,
    )
end

function _schwarzschild_algebraically_special_frequency(
        mode::QNMMode, ::Type{T}) where {T<:AbstractFloat}
    coefficient = T((mode.l - 1) * mode.l * (mode.l + 1) *
        (mode.l + 2)) / T(12)
    return complex(zero(T), -coefficient)
end

function _algebraically_special_endpoint_result(mode::QNMMode, a::T;
        angular_order::Int, sheet_id::Symbol,
        requested_inversion_index,
        convention::Symbol=:overtone) where {T<:AbstractFloat}
    omega = _schwarzschild_algebraically_special_frequency(mode, T)
    angular = angular_branch(mode, a, omega;
        truncation_order=angular_order, sheet_id)
    nan_complex = complex(T(NaN), T(NaN))
    contract = _schwarzschild_label_contract(mode, convention)
    return LeaverResult(
        mode, convention, mode.n, missing, a, omega,
        angular.angular_A, angular.lambda, angular.mixing,
        nan_complex, T(NaN), 0, zero(T), T(angular.residual),
        _precision_bits(T), :accepted, :exact_algebraically_special_endpoint,
        (
            input_type=string(typeof(a)),
            working_type=string(T),
            guess_source=:exact_algebraically_special_formula,
            label_convention=contract.label_convention,
            branch_label=contract.branch_label,
            endpoint_kind=contract.endpoint_kind,
            schwarzschild_sequence_index=missing,
            schwarzschild_multiplicity=1,
            endpoint_formula=
                :minus_i_lminus1_l_lplus1_lplus2_over_12,
            endpoint_frequency=omega,
            root_equation=:exact_algebraically_special_formula,
            radial_cf_applicable=false,
            requested_inversion_index,
            inversion_index=missing,
            angular_route=angular.metadata.route,
            angular_sheet=angular.sheet_id,
            angular_order=angular.truncation_order,
            root_drift=zero(T),
            family=:schwarzschild_algebraically_special,
            matsubara_index=8,
            matsubara_frequency=omega,
            pole_skipping=true,
            qnm_criterion=:incident_amplitude_pole_zero_collision,
            gsn_amplitude_limit_implemented=false,
            claim_boundary=
                :exact_frequency_and_pole_skipping_endpoint_not_simple_gsn_residue,
        ),
    )
end

function _resolve_continuation_coordinate(value, a, mode::QNMMode)
    coordinate = Symbol(lowercase(String(value)))
    if coordinate == :auto
        low_overtone_zdm = mode.n <= 4 &&
            _near_extremal_zdm_seed_scope(mode, a)
        # A nu path is efficient for a single near-horizon ZDM endpoint, but
        # it does not by itself preserve the overtone label once damped and
        # zero-damping branches crowd together. Track n >= 5 from the
        # Schwarzschild sequence in spin unless the caller explicitly asks
        # for a nu continuation.
        return low_overtone_zdm ? :nu : :spin
    end
    coordinate in (:spin, :nu) || throw(ArgumentError(
        "continuation_coordinate must be :auto, :spin, or :nu."))
    return coordinate
end

# Classify a near-extremal root as zero-damping (:ZDM), damped (:DM) or :unclassified.
function _near_extremal_family(mode::QNMMode, a, omega)
    abs(a) >= oftype(a, 0.99) || return :unclassified
    kappa = sqrt(max(zero(a), one(a) - a^2))
    rplus = one(a) + kappa
    horizon_frequency = a / (2rplus)
    detuning = abs(omega - mode.m * horizon_frequency)
    damping = abs(imag(omega))
    near_horizon_scale = max(oftype(a, 50) * kappa,
        oftype(a, 256) * eps(typeof(a)))
    if mode.m * a > zero(a) &&
            detuning <= near_horizon_scale &&
            damping <= near_horizon_scale
        return :ZDM
    elseif damping > oftype(a, 0.05) &&
            max(detuning, damping) > near_horizon_scale
        return :DM
    end
    return :unclassified
end

function _extremal_zdm_endpoint_scope(mode::QNMMode, a)
    # Equatorial corotating families have the unambiguous synchronous ZDM
    # endpoint. Other families can contain both zero-damping and damped
    # branches, so their endpoint is classified by continuation instead.
    abs(mode.m) == mode.l || return false
    synchronous_frequency = mode.m * a / 2
    iszero(synchronous_frequency) && return false
    return mode.branch == :positive_real ?
        real(synchronous_frequency) > 0 : real(synchronous_frequency) < 0
end

function _near_extremal_zdm_seed_scope(mode::QNMMode, a)
    abs(a) >= oftype(a, 0.999) || return false
    abs(mode.m) == mode.l || return false
    return _extremal_zdm_endpoint_scope(mode, sign(a))
end

function _near_extremal_zdm_seed(mode::QNMMode, a::T,
        inversion_index::Int) where {T<:AbstractFloat}
    horizon_gap = sqrt(max(zero(T), one(T) - a^2))
    rplus = one(T) + horizon_gap
    horizon_frequency = a / (2rplus)
    surface_gravity = horizon_gap / (2rplus)
    omega = complex(mode.m * horizon_frequency,
        -surface_gravity * (T(mode.n) + T(0.5)))
    return (
        seed=omega,
        seed_spin=a,
        source=:near_extremal_zdm_asymptotic_seed,
        seed_inversion_index=inversion_index,
        seed_predictor=:synchronous_frequency_surface_gravity,
        seed_duplicate_distance=T(Inf),
    )
end

function _extremal_zdm_endpoint_result(mode::QNMMode, a::T;
        angular_order::Int, sheet_id::Symbol, convention::Symbol,
        requested_inversion_index) where {T<:AbstractFloat}
    omega = complex(T(mode.m) * a / 2, zero(T))
    angular = angular_branch(mode, a, omega;
        truncation_order=angular_order, sheet_id)
    contract = _schwarzschild_label_contract(mode, convention)
    nan_complex = complex(T(NaN), T(NaN))
    return LeaverResult(
        mode, convention, convention == :leaver ? missing : mode.n,
        missing, a, omega, angular.angular_A, angular.lambda,
        angular.mixing, nan_complex, T(NaN), 0, zero(T),
        T(angular.residual), _precision_bits(T), :accepted,
        :exact_extremal_zdm_accumulation_endpoint,
        (
            input_type=string(typeof(a)),
            working_type=string(T),
            label_convention=contract.label_convention,
            branch_label=contract.branch_label,
            endpoint_kind=contract.endpoint_kind,
            extremal_endpoint=true,
            endpoint_method=:analytic_synchronous_frequency,
            spectrum_object=:zdm_branch_point_accumulation,
            simple_pole=false,
            coalesced_overtone_labels=true,
            requested_overtone_label=mode.n,
            distinct_overtone_frequency=false,
            radial_observables=:not_applicable_to_isolated_simple_pole,
            gsn_amplitude_limit_implemented=false,
            horizon_frequency=a / 2,
            synchronous_frequency=omega,
            requested_inversion_index,
            inversion_index=missing,
            angular_route=angular.metadata.route,
            angular_sheet=angular.sheet_id,
            angular_order=angular.truncation_order,
            root_equation=:analytic_extremal_synchronous_limit,
            family=:ZDM,
            claim_boundary=
                :frequency_endpoint_and_branch_point_not_simple_pole_residue,
        ),
    )
end

function _spin_sequence_anchor(mode::QNMMode, ::Type{T}) where {
        T<:AbstractFloat}
    if mode.m < 0 && mode.n >= 5
        # Counterrotating damped branches remain separated from the ZDM
        # accumulation but become poorly conditioned close to extremality.
        # Continue to the last supported finite-spin guard instead of making
        # one cubic jump from a=0.9995.
        return T(0.99999)
    end
    return T(0.9995)
end

function _default_cf_maximum_iterations(mode::QNMMode, a=nothing,
        guess=nothing)
    # This is a ceiling, not a forced iteration count. Near-extremal
    # retrograde modes can require the same depth as high overtones.
    baseline = mode.n >= 8 ? 48000 : 12000
    high_damping_near_extremal = a !== nothing && guess !== nothing &&
        abs(a) >= 0.99 && abs(imag(guess)) >= 0.3
    high_damping_near_extremal || return baseline
    ceiling = mode.branch == :negative_real ? 240000 : 120000
    return max(baseline, ceiling)
end

function _root_derivative_step(::Type{Float64}, omega)
    return eps(Float64)^(1 / 5) * max(1.0, abs(omega))
end

function _root_derivative_step(::Type{BigFloat}, omega)
    return eps(BigFloat)^(BigFloat(1) / 5) *
        max(one(BigFloat), abs(omega))
end

# Real step of the four-point stencil omega +- h, omega +- 2h. Outside the
# NIA-crossing exemption it is clamped to a quarter of the distance to the
# axis, so all four points stay on the physical side of the Leaver tail cut.
function _stencil_step(mode::QNMMode, omega::Complex{T},
        derivative_step_scale=one(T)) where {T<:AbstractFloat}
    h = T(_root_derivative_step(T, omega)) * T(derivative_step_scale)
    _negative_imaginary_axis_crossing_exempt(mode) && return h
    return min(h, _branch_sign(mode) * real(omega) / 4)
end

function _near_extremal_derivative_scale(mode::QNMMode, a::T) where {
        T<:AbstractFloat}
    abs(a) >= T(0.99) && mode.n >= 5 || return one(T)
    kappa = sqrt(max(zero(T), one(T) - a^2))
    return clamp(T(2) * kappa, T(0.002), one(T))
end

function _conditioned_root_residual_limit(mode::QNMMode, a::T,
        root_tolerance::T) where {T<:AbstractFloat}
    abs(a) >= T(0.999) && mode.n >= 5 || return root_tolerance
    kappa = sqrt(max(zero(T), one(T) - a^2))
    conditioning = clamp(T(0.01) / max(kappa, sqrt(eps(T))), one(T), T(10))
    return max(root_tolerance, T(1.0e-8) * conditioning)
end

function _conditioned_root_stability_limit(mode::QNMMode, a::T,
        stability_tolerance::T) where {T<:AbstractFloat}
    abs(a) >= T(0.999) && mode.n >= 5 || return stability_tolerance
    return max(stability_tolerance, T(1.0e-6))
end

# Seed for counterrotating (-2, 2, |2|) overtones at abs(a) >= 0.9999: an anchor
# value, or linear interpolation in the surface-gravity variable kappa.
function _counterrotating_near_extremal_seed(
        mode::QNMMode, a::T) where {T<:AbstractFloat}
    mode.s == -2 && mode.l == 2 && abs(mode.m) == 2 || return nothing
    mode.branch == :positive_real || return nothing
    mode.m * a < zero(T) || return nothing
    mode.n in (4, 5, 6, 7, 9, 10) || return nothing
    abs(a) >= T(first(_COUNTERROTATING_ANCHOR_SPINS)) || return nothing

    spin = abs(a)
    spin1, spin2 = T.(_COUNTERROTATING_ANCHOR_SPINS)
    omega1 = _complex_of_type(T,
        _COUNTERROTATING_NEAR_EXTREMAL_ANCHORS[(Float64(spin1), mode.n)])
    omega2 = _complex_of_type(T,
        _COUNTERROTATING_NEAR_EXTREMAL_ANCHORS[(Float64(spin2), mode.n)])
    if spin == spin1
        return (seed=omega1, source=:frozen_qnm_0p4p4_exact_anchor,
            anchor_spins=(spin1, spin2), interpolation_weight=zero(T))
    elseif spin == spin2
        return (seed=omega2, source=:frozen_qnm_0p4p4_exact_anchor,
            anchor_spins=(spin1, spin2), interpolation_weight=one(T))
    end

    kappa = sqrt(max(zero(T), one(T) - spin^2))
    kappa1 = sqrt(max(zero(T), one(T) - spin1^2))
    kappa2 = sqrt(max(zero(T), one(T) - spin2^2))
    weight = (kappa - kappa1) / (kappa2 - kappa1)
    return (
        seed=omega1 + weight * (omega2 - omega1),
        source=:frozen_qnm_0p4p4_kappa_predictor,
        anchor_spins=(spin1, spin2),
        interpolation_weight=weight,
    )
end

# Evaluate the counterrotating anchor seed directly at the target spin and gate
# it as accepted, estimated (conditioned) or failed.
function _counterrotating_conditioned_result(
        mode::QNMMode, a::T, seed;
        convention::Symbol, overtone_index,
        inversion_index::Int, angular_order::Int, sheet_id::Symbol,
        cf_tolerance::T, cf_minimum_iterations::Int,
        cf_maximum_iterations::Int, root_tolerance::T,
        stability_tolerance::T) where {T<:AbstractFloat}
    primary = leaver_radial_residual(mode, a, seed.seed;
        angular_order, sheet_id, inversion_index, cf_tolerance,
        cf_minimum_iterations, cf_maximum_iterations,
        require_cf_convergence=false)
    guard = leaver_radial_residual(mode, a, seed.seed;
        angular_order=angular_order + 8,
        sheet_id=Symbol(sheet_id, :_counterrotating_guard),
        inversion_index, cf_tolerance,
        cf_minimum_iterations=cf_minimum_iterations + 100,
        cf_maximum_iterations, require_cf_convergence=false)
    root_residual = T(abs(guard.value))
    angular_drift = T(abs(guard.angular.angular_A -
        primary.angular.angular_A) / max(one(T),
        abs(guard.angular.angular_A), abs(primary.angular.angular_A)))
    strict = root_residual <= root_tolerance &&
        guard.error <= cf_tolerance && guard.cf_converged &&
        angular_drift <= stability_tolerance
    conditioned_limit = max(
        _conditioned_root_residual_limit(mode, a, root_tolerance), T(1.0e-6))
    conditioned = !strict && root_residual <= conditioned_limit &&
        guard.error <= max(cf_tolerance, T(1.0e-10)) &&
        isfinite(guard.error) && angular_drift <= stability_tolerance
    status = strict ? :accepted : conditioned ? :estimated : :failed
    stop_reason = strict ? :accepted : conditioned ?
        :conditioned_near_extremal_reference_seed :
        root_residual > conditioned_limit ? :reference_seed_residual_gate :
        angular_drift > stability_tolerance ? :angular_stability_gate :
        :lentz_error_gate
    contract = _schwarzschild_label_contract(mode, convention)
    return LeaverResult(
        mode, convention, overtone_index, inversion_index, a, seed.seed,
        guard.angular.angular_A, guard.angular.lambda, guard.angular.mixing,
        guard.value, T(guard.error), guard.iterations, root_residual,
        T(guard.angular.residual), _precision_bits(T), status, stop_reason,
        (
            input_type=string(typeof(a)),
            working_type=string(T),
            guess_source=seed.source,
            reference_role=:branch_selection_only,
            target_recomputation=:angular_and_leaver_cf,
            upstream_qnm_version="0.4.4",
            upstream_qnm_commit=
                "f3abd18e59828e7e7d75d07f20c7cbc87925edfa",
            anchor_spins=seed.anchor_spins,
            interpolation_weight=seed.interpolation_weight,
            convention,
            overtone_index,
            inversion_index,
            label_convention=contract.label_convention,
            branch_label=contract.branch_label,
            endpoint_kind=contract.endpoint_kind,
            near_as_role=contract.near_as_role,
            schwarzschild_sequence_index=contract.schwarzschild_sequence_index,
            schwarzschild_multiplicity=contract.schwarzschild_multiplicity,
            family=_near_extremal_family(mode, a, seed.seed),
            continuation_coordinate=:frozen_kappa_predictor,
            angular_route=guard.angular.metadata.route,
            angular_sheet=guard.angular.sheet_id,
            angular_order=guard.angular.truncation_order,
            angular_order_drift=angular_drift,
            cf_minimum_iterations=cf_minimum_iterations + 100,
            cf_maximum_iterations,
            cf_tolerance,
            root_tolerance,
            conditioned_root_limit=conditioned_limit,
            conditioned_root=conditioned,
            stability_tolerance,
            root_drift=zero(T),
            root_equation=:leaver_continued_fraction,
            acceptance_basis=strict ? :direct_target_cf :
                :conditioned_reference_branch_target_cf,
            cf_tail_gate_passed=guard.error <= cf_tolerance,
            radial_cf_applicable=true,
            scientific_acceptance=status == :accepted,
        ),
    )
end

function _polynomial_predictor(xs, ys, x)
    length(xs) == length(ys) || throw(ArgumentError(
        "Continuation abscissas and values must have equal lengths."))
    isempty(xs) && throw(ArgumentError(
        "The continuation predictor requires at least one point."))
    value = zero(first(ys))
    for index in eachindex(xs)
        weight = one(x)
        for other in eachindex(xs)
            index == other && continue
            weight *= (x - xs[other]) / (xs[index] - xs[other])
        end
        value += weight * ys[index]
    end
    return value
end

function _continuation_curvature(spins, roots)
    length(spins) >= 3 || return zero(real(first(roots)))
    x1, x2, x3 = spins[(end - 2):end]
    y1, y2, y3 = roots[(end - 2):end]
    slope12 = (y2 - y1) / (x2 - x1)
    slope23 = (y3 - y2) / (x3 - x2)
    return 2abs((slope23 - slope12) / (x3 - x1))
end

function _continuation_predictor(spins, roots, next_spin)
    length(spins) >= 3 && return _polynomial_predictor(
        spins[(end - 2):end], roots[(end - 2):end], next_spin)
    if length(spins) == 2
        slope = (roots[end] - roots[end - 1]) /
            (spins[end] - spins[end - 1])
        return roots[end] + slope * (next_spin - spins[end])
    end
    return roots[end]
end


# Not-a-knot cubic spline through the accepted path, evaluated at `target_spin`
# (interpolating polynomial for fewer than four points).
function _terminal_extrapolation_predictor(spins, roots, target_spin)
    count = length(spins)
    count >= 4 || return _polynomial_predictor(
        spins, roots, target_spin)
    T = typeof(first(spins))
    matrix = zeros(T, count, count)
    rhs = zeros(eltype(roots), count)
    steps = diff(spins)

    matrix[1, 1] = -steps[2]
    matrix[1, 2] = steps[1] + steps[2]
    matrix[1, 3] = -steps[1]
    for index in 2:(count - 1)
        left_step = steps[index - 1]
        right_step = steps[index]
        matrix[index, index - 1] = left_step
        matrix[index, index] = 2(left_step + right_step)
        matrix[index, index + 1] = right_step
        rhs[index] = 6((roots[index + 1] - roots[index]) / right_step -
            (roots[index] - roots[index - 1]) / left_step)
    end
    previous_step = steps[end - 1]
    last_step = steps[end]
    matrix[end, end - 2] = last_step
    matrix[end, end - 1] = -(last_step + previous_step)
    matrix[end, end] = previous_step
    second_derivatives = matrix \ rhs

    left = spins[end - 1]
    right = spins[end]
    width = right - left
    left_weight = (right - target_spin) / width
    right_weight = (target_spin - left) / width
    return second_derivatives[end - 1] *
        (right - target_spin)^3 / (6width) +
        second_derivatives[end] *
        (target_spin - left)^3 / (6width) +
        (roots[end - 1] - second_derivatives[end - 1] * width^2 / 6) *
        left_weight +
        (roots[end] - second_derivatives[end] * width^2 / 6) *
        right_weight
end

# Accept a continuation candidate only if it lies within an allowance of the
# predictor set by the recent root motion.
function _continuation_branch_guard(spins, roots, next_spin, predictor,
        candidate, root_tolerance)
    T = typeof(real(candidate))
    length(spins) < 2 && return (accepted=true, discrepancy=zero(T),
        allowance=T(Inf))
    previous_step = abs(spins[end] - spins[end - 1])
    step_ratio = abs(next_spin - spins[end]) / previous_step
    history_move = abs(roots[end] - roots[end - 1]) * step_ratio
    predicted_move = abs(predictor - roots[end])
    scale = max(one(T), abs(candidate), abs(predictor))
    allowance = max(
        T(32) * sqrt(T(root_tolerance)) * scale,
        T(0.2) * max(history_move, predicted_move),
    )
    discrepancy = abs(candidate - predictor)
    return (accepted=discrepancy <= allowance,
        discrepancy=discrepancy, allowance=allowance)
end

# Damped Newton polish of the Leaver residual with a four-point derivative,
# optional deflation by `excluded_roots`, and an optional position certificate
# (`position_tolerance`) for final roots next to the NIA.
function _polish_leaver_root(mode::QNMMode, a::T, initial::Complex{T};
        angular_order::Int,
        sheet_id::Symbol,
        inversion_index::Int=mode.n,
        excluded_roots=Complex{T}[],
        cf_tolerance,
        cf_minimum_iterations::Int,
        cf_maximum_iterations::Int,
        root_tolerance::T,
        maximum_root_iterations::Int,
        trust_radius::T,
        derivative_step_scale::T=one(T),
        position_tolerance=nothing) where {T<:AbstractFloat}
    omega = initial
    # The Leaver tail (ContinuedFractions.jl, sqrt(-(D1+D2+D3))) has a branch
    # cut on the NIA, so a root on the wrong side of it is not this branch.
    clamp_stencil = !_negative_imaginary_axis_crossing_exempt(mode)
    if clamp_stencil && _branch_sign(mode) * real(omega) <= zero(T)
        return (
            omega=omega, evaluation=nothing, history=NamedTuple[],
            converged=false, conditioned=false,
            stop_reason=:initial_guess_off_branch)
    end
    evaluation = leaver_radial_residual(mode, a, omega;
        angular_order, sheet_id, inversion_index, cf_tolerance,
        cf_minimum_iterations, cf_maximum_iterations,
        require_cf_convergence=false)
    history = NamedTuple[]

    function deflated_value(value, frequency)
        isempty(excluded_roots) && return value
        separation = minimum(abs(frequency - root) for root in excluded_roots)
        scale = max(one(T), abs(frequency),
            maximum(abs, excluded_roots))
        separation <= T(64) * eps(T) * scale &&
            return complex(T(Inf), zero(T))
        denominator = prod(frequency - root for root in excluded_roots)
        return value / denominator
    end

    for iteration in 1:maximum_root_iterations
        residual = T(abs(evaluation.value))
        deflated_residual = T(abs(deflated_value(evaluation.value, omega)))
        push!(history, (
            iteration=iteration,
            omega=omega,
            residual=residual,
            deflated_residual=deflated_residual,
            cf_error=evaluation.error,
            cf_iterations=evaluation.iterations,
            angular_residual=evaluation.angular.residual,
        ))
        residual <= root_tolerance && position_tolerance === nothing && return (
            omega=omega, evaluation=evaluation, history=history,
            converged=true, conditioned=false,
            stop_reason=:root_tolerance)

        h = _stencil_step(mode, omega, derivative_step_scale)
        if position_tolerance !== nothing
            isempty(excluded_roots) || throw(ArgumentError(
                "Position certification requires the undeflated final root."))
            h = min(h, _branch_sign(mode) * real(omega) / 16)
            if !isfinite(h) || h <= zero(T) || omega + h/2 == omega ||
                    omega - h/2 == omega || omega + im*h == omega
                return (omega, evaluation, history, converged=false,
                    conditioned=false, stop_reason=:near_axis_position_stencil_unverified)
            end
        end
        xplus = leaver_radial_residual(mode, a, omega + h;
            angular_order, sheet_id, inversion_index, cf_tolerance,
            cf_minimum_iterations, cf_maximum_iterations,
            require_cf_convergence=false)
        xminus = leaver_radial_residual(mode, a, omega - h;
            angular_order, sheet_id, inversion_index, cf_tolerance,
            cf_minimum_iterations, cf_maximum_iterations,
            require_cf_convergence=false)
        xplus2 = leaver_radial_residual(mode, a, omega + 2h;
            angular_order, sheet_id, inversion_index, cf_tolerance,
            cf_minimum_iterations, cf_maximum_iterations,
            require_cf_convergence=false)
        xminus2 = leaver_radial_residual(mode, a, omega - 2h;
            angular_order, sheet_id, inversion_index, cf_tolerance,
            cf_minimum_iterations, cf_maximum_iterations,
            require_cf_convergence=false)
        derivative = (-deflated_value(xplus2.value, omega + 2h) +
            8 * deflated_value(xplus.value, omega + h) -
            8 * deflated_value(xminus.value, omega - h) +
            deflated_value(xminus2.value, omega - 2h)) / (12h)
        # Deflation can scale a regular derivative arbitrarily small.
        # Root acceptance is checked against the undeflated residual.
        if !isfinite(derivative) || iszero(derivative)
            return (
                omega=omega, evaluation=evaluation, history=history,
                converged=false, conditioned=false,
                stop_reason=:singular_root_jacobian)
        end
        delta = -deflated_value(evaluation.value, omega) / derivative
        if position_tolerance !== nothing
            # The first stencil is the ordinary one; at most two halved
            # stencils resolve fourth-order truncation error.
            for stencil_attempt in 0:2
                if stencil_attempt > 0
                    h /= 2
                    if !isfinite(h) || h <= zero(T) || omega + h/2 == omega ||
                            omega - h/2 == omega || omega + im*h == omega
                        return (omega=omega, evaluation=evaluation, history=history,
                            converged=false, conditioned=false,
                            stop_reason=:near_axis_position_stencil_unverified)
                    end
                    xplus, xminus, xplus2, xminus2 = map((h, -h, 2h, -2h)) do offset
                        leaver_radial_residual(mode, a, omega + offset;
                            angular_order, sheet_id, inversion_index, cf_tolerance,
                            cf_minimum_iterations, cf_maximum_iterations,
                            require_cf_convergence=false)
                    end
                    derivative = (-xplus2.value + 8 * xplus.value -
                        8 * xminus.value + xminus2.value) / (12h)
                end
                # Final BigFloat roots near the NIA need a position certificate:
                # the residual can be tiny while |F/F'| is still material.
                cf_valid(e) = isfinite(e.value) && isfinite(e.error) &&
                    e.cf_converged && e.error <= cf_tolerance &&
                    isfinite(e.angular.residual)
                extra = map((h/2, -h/2, im*h, -im*h, 2im*h, -2im*h)) do offset
                    leaver_radial_residual(mode, a, omega + offset;
                        angular_order, sheet_id, inversion_index, cf_tolerance,
                        cf_minimum_iterations, cf_maximum_iterations,
                        require_cf_convergence=false)
                end
                if !all(cf_valid, (evaluation, xplus, xminus, xplus2, xminus2, extra...))
                    push!(history, (iteration=iteration, omega=omega, position_verified=false,
                        stop_reason=:near_axis_position_cf_unverified))
                    return (omega, evaluation, history, converged=false,
                        conditioned=false, stop_reason=:near_axis_position_cf_unverified)
                end
                hp, hm, ip, iminus, ip2, iminus2 = extra
                half_derivative = (-xplus.value + 8hp.value - 8hm.value +
                    xminus.value) / (6h)
                imaginary_derivative = (-ip2.value + 8ip.value -
                    8iminus.value + iminus2.value) / (12im*h)
                finite_derivatives = isfinite(half_derivative) &&
                    isfinite(imaginary_derivative) && !iszero(half_derivative)
                derivative_scale = max(abs(half_derivative), eps(T))
                half_discrepancy = abs(derivative - half_derivative) / derivative_scale
                imaginary_discrepancy = abs(imaginary_derivative - half_derivative) /
                    derivative_scale
                derivative_verified = finite_derivatives &&
                    max(half_discrepancy, imaginary_discrepancy) <= position_tolerance
                delta = -evaluation.value / half_derivative
                scaled_correction = abs(delta) / max(one(T), abs(omega))
                position_verified = derivative_verified && isfinite(delta) &&
                    scaled_correction <= position_tolerance && residual <= root_tolerance
                position_check = (verified=position_verified,
                    precision_bits=_precision_bits(T), omega, h, stencil_attempt,
                    residual, cf_error=evaluation.error,
                    cf_iterations=evaluation.iterations, cf_converged=evaluation.cf_converged,
                    all_samples_cf_verified=true,
                    sample_offsets=(zero(omega), h, -h, 2h, -2h,
                        h/2, -h/2, im*h, -im*h, 2im*h, -2im*h),
                    sample_cf_errors=map(e -> e.error,
                        (evaluation, xplus, xminus, xplus2, xminus2, extra...)),
                    sample_cf_iterations=map(e -> e.iterations,
                        (evaluation, xplus, xminus, xplus2, xminus2, extra...)),
                    sample_angular_residuals=map(e -> e.angular.residual,
                        (evaluation, xplus, xminus, xplus2, xminus2, extra...)),
                    derivative, half_derivative, imaginary_derivative,
                    half_discrepancy, imaginary_discrepancy,
                    uncapped_correction=delta, scaled_correction,
                    root_tolerance, cf_tolerance, position_tolerance,
                    certificate=:local_newton_position_estimate)
                push!(history, (iteration=iteration, omega=omega, position_check=position_check))
                if derivative_verified && isfinite(delta)
                    position_verified && return (omega=omega, evaluation=evaluation,
                        history=history, converged=true, conditioned=false,
                        stop_reason=:root_position_tolerance, position_check=position_check)
                    break
                elseif stencil_attempt == 2
                    return (omega=omega, evaluation=evaluation, history=history,
                        converged=false, conditioned=false,
                        stop_reason=:near_axis_position_derivative_unverified,
                        position_check=position_check)
                end
            end
        end
        if abs(delta) > trust_radius
            delta *= trust_radius / abs(delta)
        end

        accepted = false
        for backtrack in 0:10
            scale = T(2)^(-backtrack)
            candidate = omega + scale * delta
            _branch_consistent(mode, candidate) || continue
            candidate_evaluation = leaver_radial_residual(
                mode, a, candidate;
                angular_order, sheet_id, inversion_index, cf_tolerance,
                cf_minimum_iterations, cf_maximum_iterations,
                require_cf_convergence=false)
            candidate_cf_valid = position_tolerance === nothing ||
                (isfinite(candidate_evaluation.value) &&
                 isfinite(candidate_evaluation.error) &&
                 candidate_evaluation.cf_converged &&
                 candidate_evaluation.error <= cf_tolerance)
            if candidate_cf_valid &&
                    abs(deflated_value(candidate_evaluation.value, candidate)) <
                    abs(deflated_value(evaluation.value, omega))
                omega = candidate
                evaluation = candidate_evaluation
                accepted = true
                break
            end
        end
        if !accepted
            conditioned_limit = _conditioned_root_residual_limit(
                mode, a, root_tolerance)
            conditioned = position_tolerance === nothing &&
                residual <= conditioned_limit &&
                evaluation.cf_converged && evaluation.error <= cf_tolerance
            return (
                omega=omega, evaluation=evaluation, history=history,
                converged=conditioned, conditioned,
                stop_reason=conditioned ? :conditioned_cf_floor :
                    :root_line_search_stagnation)
        end
    end
    final_residual = T(abs(evaluation.value))
    conditioned_limit = _conditioned_root_residual_limit(
        mode, a, root_tolerance)
    # With position checking enabled, no residual-only iteration-limit exit
    # may bypass the explicit derivative and uncapped-correction certificate.
    strict = position_tolerance === nothing && final_residual <= root_tolerance
    conditioned = position_tolerance === nothing && !strict &&
        final_residual <= conditioned_limit &&
        evaluation.cf_converged && evaluation.error <= cf_tolerance
    return (
        omega=omega, evaluation=evaluation, history=history,
        converged=strict || conditioned, conditioned,
        stop_reason=strict ? :root_tolerance : conditioned ?
            :conditioned_cf_floor : :root_iteration_limit)
end

_near_as_complete_scope(mode::QNMMode, near_as_role::Symbol) =
    mode.s == -2 && mode.l == 2 && abs(mode.m) == 2 &&
    near_as_role in (:unconventional, :algebraically_special) &&
    mode.branch == :positive_real

function _near_as_mano_suzuki_takasugi_required(mode::QNMMode, a, ::Type{T},
        near_as_role::Symbol) where {T<:AbstractFloat}
    return T === Float64 && _near_as_complete_scope(mode, near_as_role) &&
        !iszero(a) && abs(a) < T(_KERR_L2_M2_MANO_SUZUKI_TAKASUGI_SWITCH_SPIN)
end

function _near_as_mano_suzuki_takasugi_finalization_required(
        mode::QNMMode, a, ::Type{T}, near_as_role::Symbol) where {
        T<:AbstractFloat}
    return T === Float64 &&
        _near_as_complete_scope(mode, near_as_role) &&
        near_as_role == :algebraically_special &&
        abs(a) >= T(_KERR_L2_M2_MANO_SUZUKI_TAKASUGI_SWITCH_SPIN) &&
        abs(a) < T(_KERR_L2_M2_NEAR_AS_ANCHOR_SPIN)
end

function _mano_suzuki_takasugi_incidence_evaluation(mode::QNMMode, a::Float64,
        omega::ComplexF64; angular_order::Int, sheet_id::Symbol)
    route = ISEM.complex_route(
        mode.s, mode.l, mode.m, a, omega, :IN; backend=:mst)
    scale = max(1.0, abs(route.transmission), abs(route.incidence),
        abs(route.reflection))
    metadata = route.metadata
    representation_error = maximum((
        abs(metadata.mst_amplitude_truncation),
        abs(metadata.mst_amplitude_spread),
        abs(metadata.mst_amplitude_nearest),
        abs(metadata.mst_nu_residual),
        abs(metadata.mst_anchor_error),
        abs(metadata.mst_anchor_residual),
    ))
    angular = angular_branch(mode, a, omega;
        truncation_order=angular_order, sheet_id)
    return (
        value=route.incidence,
        scaled_residual=abs(route.incidence) / scale,
        error=Float64(representation_error),
        iterations=0,
        angular,
        route,
        scale,
        root_equation=:gsn_mst_incidence,
    )
end

# Newton polish of the GSN incidence amplitude from the MST backend; used for
# the near-AS pair at small nonzero spin.
function _polish_mano_suzuki_takasugi_incidence_root(mode::QNMMode, a::Float64,
        initial::ComplexF64; angular_order::Int, sheet_id::Symbol,
        root_tolerance::Float64, maximum_root_iterations::Int,
        trust_radius::Float64, derivative_step_scale::Float64=1.0)
    omega = initial
    evaluation = _mano_suzuki_takasugi_incidence_evaluation(
        mode, a, omega; angular_order, sheet_id)
    history = NamedTuple[]
    for iteration in 1:maximum_root_iterations
        residual = evaluation.scaled_residual
        push!(history, (
            iteration,
            omega,
            residual,
            deflated_residual=residual,
            cf_error=evaluation.error,
            cf_iterations=0,
            angular_residual=evaluation.angular.residual,
            root_equation=:gsn_mst_incidence,
        ))
        residual <= root_tolerance && return (
            omega, evaluation, history, converged=true,
            stop_reason=:root_tolerance)

        h = derivative_step_scale * 2.0e-5 * max(1.0, abs(omega))
        minus2 = _mano_suzuki_takasugi_incidence_evaluation(
            mode, a, omega - 2h; angular_order, sheet_id)
        minus1 = _mano_suzuki_takasugi_incidence_evaluation(
            mode, a, omega - h; angular_order, sheet_id)
        plus1 = _mano_suzuki_takasugi_incidence_evaluation(
            mode, a, omega + h; angular_order, sheet_id)
        plus2 = _mano_suzuki_takasugi_incidence_evaluation(
            mode, a, omega + 2h; angular_order, sheet_id)
        derivative = (minus2.value - 8minus1.value + 8plus1.value -
            plus2.value) / (12h)
        if !isfinite(derivative) || abs(derivative) <= sqrt(eps(Float64))
            return (
                omega, evaluation, history, converged=false,
                stop_reason=:singular_root_jacobian)
        end
        delta = -evaluation.value / derivative
        if abs(delta) > trust_radius
            delta *= trust_radius / abs(delta)
        end

        accepted = false
        for backtrack in 0:10
            candidate = omega + delta / 2.0^backtrack
            candidate_evaluation = _mano_suzuki_takasugi_incidence_evaluation(
                mode, a, candidate; angular_order, sheet_id)
            if candidate_evaluation.scaled_residual < residual
                omega = candidate
                evaluation = candidate_evaluation
                accepted = true
                break
            end
        end
        accepted || return (
            omega, evaluation, history, converged=false,
            stop_reason=:root_line_search_stagnation)
    end
    return (
        omega, evaluation, history,
        converged=evaluation.scaled_residual <= root_tolerance,
        stop_reason=evaluation.scaled_residual <= root_tolerance ?
            :root_tolerance : :root_iteration_limit)
end

function _polish_frequency_root(mode::QNMMode, a::T, initial::Complex{T};
        angular_order::Int, sheet_id::Symbol, inversion_index::Int,
        excluded_roots=Complex{T}[],
        cf_tolerance, cf_minimum_iterations::Int,
        cf_maximum_iterations::Int, root_tolerance::T,
        maximum_root_iterations::Int, trust_radius::T,
        derivative_step_scale=one(T),
        near_as_role::Symbol=:none,
        position_tolerance=nothing) where {T<:AbstractFloat}
    if _near_as_mano_suzuki_takasugi_required(mode, a, T, near_as_role)
        return _polish_mano_suzuki_takasugi_incidence_root(
            mode, Float64(a), ComplexF64(initial);
            angular_order, sheet_id,
            root_tolerance=Float64(max(root_tolerance, T(2.0e-11))),
            maximum_root_iterations,
            trust_radius=Float64(trust_radius),
            derivative_step_scale=Float64(derivative_step_scale),
        )
    end
    return _polish_leaver_root(
        mode, a, initial;
        angular_order, sheet_id, inversion_index, excluded_roots,
        cf_tolerance, cf_minimum_iterations, cf_maximum_iterations,
        root_tolerance, maximum_root_iterations, trust_radius,
        derivative_step_scale=T(derivative_step_scale) *
            _near_extremal_derivative_scale(mode, a), position_tolerance)
end

_root_equation(evaluation) = hasproperty(evaluation, :root_equation) ?
    evaluation.root_equation : :leaver_continued_fraction
_root_residual(evaluation) = hasproperty(evaluation, :scaled_residual) ?
    evaluation.scaled_residual : abs(evaluation.value)

function _near_as_anchor_seed(mode::QNMMode, target_spin::T,
        near_as_role::Symbol) where {T<:AbstractFloat}
    if near_as_role == :counterrotating
        seed = _complex_of_type(T, _KERR_L2_COUNTER_N8_ANCHOR)
        return (
            seed=mode.branch == :positive_real ? seed : -conj(seed),
            seed_spin=T(0.5),
            source=:published_finite_spin_anchor_repolished_locally,
            seed_inversion_index=8,
            seed_predictor=:counterrotating_a0p5,
            seed_duplicate_distance=T(Inf),
        )
    end
    _near_as_complete_scope(mode, near_as_role) || throw(ArgumentError(
        "The automatic near-AS anchor is available only for " *
        "the explicitly selected (s,l,abs(m))=(-2,2,2) near-AS branch."))
    anchor_inversion = near_as_role == :algebraically_special ? 9 : 8
    positive_seed = _complex_of_type(T,
        _KERR_L2_M2_NEAR_AS_ANCHORS[near_as_role])
    effective_spin = mode.m < 0 ? -target_spin : target_spin
    seed = effective_spin < zero(T) ? -conj(positive_seed) : positive_seed
    seed_spin = copysign(T(_KERR_L2_M2_NEAR_AS_ANCHOR_SPIN), target_spin)
    return (
        seed,
        seed_spin,
        source=:published_finite_spin_anchor_repolished_locally,
        seed_inversion_index=anchor_inversion,
        seed_predictor=mode.m < 0 ?
            :complete_spectrum_a0p05_spin_reflected :
            :complete_spectrum_a0p05,
        seed_duplicate_distance=T(Inf),
    )
end

function _failed_result(mode::QNMMode, a::T, omega::Complex{T}, reason::Symbol,
        provenance; convention::Symbol=:overtone,
        overtone_index::Union{Int,Missing}=mode.n,
        inversion_index::Int=mode.n) where {T<:AbstractFloat}
    nan_complex = complex(T(NaN), T(NaN))
    return LeaverResult(
        mode, convention, overtone_index, inversion_index,
        a, omega, nan_complex, nan_complex, Complex{T}[],
        nan_complex, T(Inf), 0, T(Inf), T(Inf), _precision_bits(T),
        :failed, reason, provenance)
end

function _schwarzschild_eikonal_guess(mode::QNMMode, ::Type{T}) where {
        T<:AbstractFloat}
    scale = T(3) * sqrt(T(3))
    return complex((T(mode.l) + T(0.5)) / scale,
        -(T(mode.n) + T(0.5)) / scale)
end

function _schwarzschild_overtone_predictor(
        roots, mode::QNMMode, ::Type{T}) where {T<:AbstractFloat}
    overtone = mode.n
    if overtone >= 2 && haskey(roots, overtone - 1) &&
            haskey(roots, overtone - 2)
        predictor = 2 * roots[overtone - 1] - roots[overtone - 2]
        if real(predictor) <= zero(T)
            predictor = complex(abs(real(predictor)), imag(predictor))
        end
        return predictor
    end
    return _schwarzschild_eikonal_guess(mode, T)
end

function _schwarzschild_inversion_candidates(overtone::Int)
    candidates = Int[]
    for candidate in (overtone, overtone - 1, overtone + 1,
            overtone - 2, overtone + 2)
        candidate >= 0 && candidate ∉ candidates && push!(candidates, candidate)
    end
    return candidates
end

# Schwarzschild (a = 0) seed for `mode`: a tabulated root, or the overtone
# sequence generated root by root with deflation of the lower overtones.
function _schwarzschild_seed(mode::QNMMode, ::Type{T};
        angular_order::Int,
        sheet_id::Symbol,
        cf_tolerance,
        cf_minimum_iterations::Int,
        cf_maximum_iterations::Int,
        root_tolerance::T,
        maximum_root_iterations::Int,
        trust_radius::T) where {T<:AbstractFloat}
    mode.branch == :positive_real || throw(ArgumentError(
        "A negative-real QNM requires an explicit directly sourced guess; " *
        "it is not inferred by conjugation."))
    target_key = (mode.s, mode.l, mode.n)
    if haskey(_SCHWARZSCHILD_SEEDS, target_key)
        return (
            seed=_complex_of_type(T, _SCHWARZSCHILD_SEEDS[target_key]),
            seed_spin=zero(T),
            source=:frozen_schwarzschild_seed,
            seed_inversion_index=mode.n,
            seed_predictor=:frozen,
            seed_duplicate_distance=T(Inf),
        )
    end

    return let roots = Dict{Int,Complex{T}}(),
            generated_metadata = Dict{Int,NamedTuple}()
        for (key, value) in _SCHWARZSCHILD_SEEDS
            key[1] == mode.s && key[2] == mode.l || continue
            roots[key[3]] = _complex_of_type(T, value)
        end
        for overtone in 0:mode.n
            haskey(roots, overtone) && continue
            if mode.s == -2 && mode.l == 2 && abs(mode.m) == 2 && overtone == 8
                # The complete-spectrum n=8/n=9 endpoints are not ordinary
                # seeds in the distinct Leaver-root sequence. Their separate
                # unconventional and pole-skipping contracts handle them.
                continue
            end
            overtone_mode = QNMMode(
                mode.s, mode.l, mode.m, overtone, :positive_real)
            predictor = _schwarzschild_overtone_predictor(
                roots, overtone_mode, T)
            eikonal = _schwarzschild_eikonal_guess(overtone_mode, T)
            predictors = predictor == eikonal ? Complex{T}[predictor] :
                Complex{T}[predictor, eikonal]
            excluded_roots = Complex{T}[
                roots[index] for index in sort(collect(keys(roots)))
                if index < overtone]
            candidates = NamedTuple[]
            trial_specs = Tuple{Int,Int,Complex{T}}[]
            for (predictor_index, trial_predictor) in pairs(predictors)
                for inversion in _schwarzschild_inversion_candidates(overtone)
                    push!(trial_specs,
                        (predictor_index, inversion, trial_predictor))
                end
            end
            for (predictor_index, inversion, trial_predictor) in trial_specs
                    polished = _polish_frequency_root(
                        overtone_mode, zero(T), trial_predictor;
                        angular_order,
                        sheet_id=Symbol(sheet_id,
                            :_schwarzschild_overtone_, inversion, :_,
                            predictor_index),
                        inversion_index=inversion,
                        excluded_roots,
                        cf_tolerance,
                        cf_minimum_iterations,
                        cf_maximum_iterations,
                        root_tolerance,
                        maximum_root_iterations=max(
                            maximum_root_iterations, 60),
                        trust_radius,
                    )
                    polished.converged || continue
                    duplicate_distance = isempty(excluded_roots) ? T(Inf) :
                        minimum(abs(polished.omega - root)
                            for root in excluded_roots)
                    duplicate_tolerance = T(1.0e-6) *
                        max(one(T), abs(polished.omega))
                    duplicate_distance > duplicate_tolerance || continue
                    overtone > 0 &&
                        imag(polished.omega) >= imag(roots[overtone - 1]) &&
                        continue
                    seed_cf_limit = max(
                        T(100) * T(cf_tolerance), sqrt(eps(T)))
                    polished.evaluation.error <= seed_cf_limit || continue
                    push!(candidates, (
                        polished,
                        inversion,
                        predictor=trial_predictor,
                        predictor_distance=abs(
                            polished.omega - trial_predictor),
                        duplicate_distance,
                    ))
                    # The physical-overtone predictor with n_inv=n is the
                    # least ambiguous route.  Once it passes all guards there
                    # is no benefit in solving the remaining fallback roots.
                    predictor_index == 1 && inversion == overtone && break
            end
            isempty(candidates) && throw(ErrorException(
                "Automatic Schwarzschild overtone continuation failed for " *
                "(s=$(mode.s), l=$(mode.l), n=$overtone): no unique " *
                "continued-fraction root survived the inversion and " *
                "duplicate-root gates."))
            selected = candidates[argmin(
                candidate.predictor_distance for candidate in candidates)]
            roots[overtone] = selected.polished.omega
            generated_metadata[overtone] = (
                inversion_index=selected.inversion,
                predictor=selected.predictor,
                duplicate_distance=Float64(selected.duplicate_distance),
                root_residual=Float64(abs(
                    selected.polished.evaluation.value)),
                cf_error=Float64(selected.polished.evaluation.error),
            )
        end
        metadata = generated_metadata[mode.n]
        return (
            seed=roots[mode.n],
            seed_spin=zero(T),
            source=:generated_schwarzschild_overtone_sequence,
            seed_inversion_index=metadata.inversion_index,
            seed_predictor=metadata.predictor,
            seed_duplicate_distance=T(metadata.duplicate_distance),
        )
    end
end

# Re-solve a failed Float64 root in BigFloat with the same acceptance gates, and
# check that the BigFloat root stays on the Float64 branch.
function _automatic_precision_fallback(mode::QNMMode, a, omega;
        angular_order, sheet_id, cf_minimum_iterations,
        cf_maximum_iterations, maximum_root_iterations, trust_radius,
        fallback_precision_bits, float64_provenance,
        cf_tolerance, root_tolerance, stability_tolerance,
        continuation_coordinate,
        convention::Symbol, overtone_index::Union{Int,Missing},
        inversion_index::Int, segment::Int=0)
    fallback_precision_bits >= 96 || throw(ArgumentError(
        "fallback_precision_bits must be at least 96."))
    high_damping_near_extremal = abs(a) >= 0.99 &&
        abs(imag(omega)) >= 0.3
    fallback_cf_ceiling = high_damping_near_extremal ?
        (mode.branch == :negative_real ? 240000 : 120000) : 24000
    fallback_cf_maximum_iterations = max(
        cf_maximum_iterations, fallback_cf_ceiling)
    return setprecision(BigFloat, fallback_precision_bits) do
        result = _qnm_frequency_labeled(
            mode, BigFloat(a);
            guess=Complex{BigFloat}(
                BigFloat(real(omega)), BigFloat(imag(omega))),
            angular_order,
            sheet_id=Symbol(sheet_id, :_automatic_bigfloat),
            # A precision fallback must re-evaluate the same scientific
            # contract. Tightening the acceptance gates here makes the result
            # depend on whether Float64 happened to request the fallback.
            cf_tolerance=BigFloat(cf_tolerance),
            cf_minimum_iterations=max(cf_minimum_iterations, 600),
            cf_maximum_iterations=fallback_cf_maximum_iterations,
            root_tolerance=BigFloat(root_tolerance),
            stability_tolerance=BigFloat(stability_tolerance),
            maximum_root_iterations,
            maximum_spin_step=BigFloat("0.004"),
            minimum_spin_step=BigFloat("1e-5"),
            trust_radius=BigFloat(trust_radius),
            automatic_precision_fallback=false,
            fallback_precision_bits,
            continuation_coordinate,
            convention,
            inversion_index,
            segment,
        )
        return setprecision(BigFloat, max(fallback_precision_bits, result.precision_bits)) do
            precision_drift = abs(result.omega - omega) /
                max(one(real(result.omega)), abs(result.omega), abs(omega))
            branch_reference = _complex_of_type(BigFloat, omega)
            raw_branch_allowance = zero(BigFloat)
            if hasproperty(float64_provenance, :terminal_predictor)
                branch_reference = _complex_of_type(
                    BigFloat, float64_provenance.terminal_predictor)
            elseif hasproperty(float64_provenance, :terminal_extrapolation_rows) &&
                    !isempty(float64_provenance.terminal_extrapolation_rows)
                terminal_row = last(float64_provenance.terminal_extrapolation_rows)
                branch_reference = _complex_of_type(BigFloat, terminal_row.predictor)
                if hasproperty(terminal_row, :branch_allowance)
                    raw_branch_allowance = BigFloat(terminal_row.branch_allowance)
                end
            end
            branch_scale = max(one(BigFloat), abs(result.omega),
                abs(branch_reference))
            precision_candidate_allowance =
                BigFloat(64) * sqrt(BigFloat(stability_tolerance))
            precision_predictor_allowance = max(
                precision_candidate_allowance,
                raw_branch_allowance / branch_scale,
            )
            precision_predictor_drift = abs(result.omega - branch_reference) /
                branch_scale
            precision_branch_drift = min(
                BigFloat(precision_drift), precision_predictor_drift)
            precision_branch_allowance = precision_candidate_allowance
            precision_branch_guard_passed =
                BigFloat(precision_drift) <= precision_candidate_allowance ||
                precision_predictor_drift <= precision_predictor_allowance
            provenance = merge(result.provenance, (
                precision_policy=:automatic_bigfloat_fallback,
                high_damping_near_extremal,
                fallback_cf_maximum_iterations,
                float64_candidate=omega,
                float64_bigfloat_drift=precision_drift,
                precision_branch_reference=branch_reference,
                precision_candidate_drift=precision_drift,
                precision_candidate_allowance,
                precision_predictor_drift,
                precision_predictor_allowance,
                precision_branch_drift,
                precision_branch_allowance,
                precision_branch_guard_passed,
                float64_provenance=float64_provenance,
                continuation_axis_touches=get(float64_provenance,
                    :continuation_axis_touches,
                    get(result.provenance, :continuation_axis_touches, 0)),
                nia_restart_rows=get(float64_provenance, :nia_restart_rows,
                    get(result.provenance, :nia_restart_rows, NamedTuple[])),
            ))
            status = result.status in (:accepted, :estimated) &&
                precision_branch_guard_passed ? result.status : :failed
            stop_reason = !(result.status in (:accepted, :estimated)) ?
                result.stop_reason :
                precision_branch_guard_passed ? result.stop_reason :
                :precision_fallback_branch_guard
            return LeaverResult(
                result.mode, result.convention, result.overtone_index,
                result.inversion_index, result.a, result.omega, result.angular_A,
                result.lambda, result.mixing, result.cf_value, result.cf_error,
                result.cf_iterations, result.root_residual,
                result.angular_residual, result.precision_bits, status,
                stop_reason, provenance)
        end
    end
end

# Seed for an ordinary root at 0 <= a <= _SPIN_TABLE_MAXIMUM_SPIN from the packaged
# spin table (SpinRootTable.jl): cubic interpolation in u = sqrt(1 - a) on the four
# nearest nodes, and the distance to the quadratic interpolant on the three nearest as
# an error estimate. Returns nothing outside the table or when a node of the stencil is
# missing or the stencil changes Leaver inversion index.
function _tabulated_spin_seed(mode::QNMMode, a)
    mode.branch == :positive_real || return nothing
    a isa Union{Float64,Integer} || return nothing
    spin = Float64(a)
    _SPIN_TABLE_MAXIMUM_SPIN < spin < one(spin) &&
        return _two_end_spin_seed(mode, spin)
    zero(spin) <= spin <= _SPIN_TABLE_MAXIMUM_SPIN || return nothing
    entry = get(_SPIN_ROOT_TABLE, (mode.s, mode.l, mode.m, mode.n), nothing)
    entry === nothing && return nothing
    roots, inversions = entry
    us = _SPIN_TABLE_U
    u = sqrt(one(spin) - spin)
    nodes = sortperm(abs.(us .- u))[1:4]
    sort!(nodes)
    inversion = inversions[nodes[1]]
    (inversion >= 0 && all(inversions[k] == inversion for k in nodes) &&
        all(isfinite, roots[nodes])) || return nothing
    lagrange(ks) = sum(roots[j] * prod((u - us[i]) / (us[j] - us[i])
        for i in ks if i != j) for j in ks)
    cubic = lagrange(nodes)
    nearest = sort(sortperm(abs.(us .- u))[1:3])
    quadratic = lagrange(nearest)
    error = abs(cubic - quadratic)
    return (seed=cubic, inversion=inversion, error,
        radius=max(10 * error, 1.0e-6 * abs(cubic)),
        source=:packaged_spin_table,
        predictor=:packaged_table_cubic_interpolation)
end

# Seed for 0.9995 < a < 1, anchored at both ends: the packaged root at a = 1
# (ExtremalRootTable.jl; the synchronous endpoint m/2 for a ZDM entry) and the two
# spin-table nodes nearest the extremal end. From the nearest node the root is marched in
# kappa = sqrt(1 - a^2) by factors of 0.8, each step predicted by the quadratic in kappa
# through the a = 1 root and the last two roots, and kept only within a quarter of the
# surface gravity sqrt(1 - a^2) / (2 r_+) of its prediction, the level-spacing fraction
# that also bounds the terminal extrapolation. Returns nothing when a step fails.
function _two_end_spin_seed(mode::QNMMode, spin::Float64)
    key = (mode.s, mode.l, mode.m, mode.n)
    entry = get(_SPIN_ROOT_TABLE, key, nothing)
    extremal = get(_EXTREMAL_ROOT_TABLE, key, nothing)
    (entry === nothing || extremal === nothing) && return nothing
    roots, inversions = entry
    nodes = findall(u -> one(u) - u^2 <= _SPIN_TABLE_MAXIMUM_SPIN, _SPIN_TABLE_U)
    near, far = nodes[end], nodes[end-1]
    inversion = inversions[near]
    (inversion >= 0 && inversions[far] == inversion &&
        isfinite(roots[near]) && isfinite(roots[far]) &&
        isfinite(extremal[2])) || return nothing
    nodal_kappa(u) = sqrt(one(u) - (one(u) - u^2)^2)
    radius(kappa) = kappa / (8 * (one(kappa) + kappa))
    target = sqrt(one(spin) - spin^2)
    kappas = [nodal_kappa(_SPIN_TABLE_U[far]), nodal_kappa(_SPIN_TABLE_U[near])]
    omegas = ComplexF64[roots[far], roots[near]]
    predictor = omegas[end]
    while kappas[end] > target
        kappa = max(target, 0.8 * kappas[end])
        k0, k1, k2 = zero(spin), kappas[end-1], kappas[end]
        w0, w1, w2 = extremal[2], omegas[end-1], omegas[end]
        predictor = w0 * (kappa - k1) * (kappa - k2) / ((k0 - k1) * (k0 - k2)) +
            w1 * (kappa - k0) * (kappa - k2) / ((k1 - k0) * (k1 - k2)) +
            w2 * (kappa - k0) * (kappa - k1) / ((k2 - k0) * (k2 - k1))
        step_spin = kappa == target ? spin : sqrt(one(spin) - kappa^2)
        polished = _polish_frequency_root(mode, step_spin, predictor;
            angular_order=32, sheet_id=:qnm_straight_from_spherical,
            inversion_index=inversion,
            cf_tolerance=mode.n >= 8 ? 1.0e-12 : 1.0e-13,
            cf_minimum_iterations=300,
            cf_maximum_iterations=_default_cf_maximum_iterations(
                mode, step_spin, nothing),
            root_tolerance=5.0e-10, maximum_root_iterations=30,
            trust_radius=0.2)
        (polished.converged &&
            abs(polished.omega - predictor) <= radius(kappa)) || return nothing
        push!(kappas, kappa)
        push!(omegas, polished.omega)
    end
    return (seed=omegas[end], inversion, error=abs(omegas[end] - predictor),
        radius=radius(target),
        source=:packaged_two_end_anchor,
        predictor=:packaged_table_kappa_anchored_march)
end

# Ordinary DM root for 0.999 < a < 1 from the near-extremal expansion
# omega = c_0 + c_1 x + ... + c_8 x^8 in x = 1 - a^2. The c_k are solved order by order from the
# Teukolsky equation at a = 1, computed on first use (0.03-0.3 s, then cached) from the
# packaged exact-extremal damped root; no coefficient is stored.
# The series misses only the exponentially small term exp(-2 pi Re(2 omega - a m) / sqrt(x)), so
# it is used for x <= x_max, where that term is below 1e-16; root_residual carries the larger of
# that term and the first omitted power. Here the series is more accurate than a Float64 Leaver
# solve, whose continued fraction cancels terms of size 1/sqrt(x). Computed coefficients must
# reproduce the packaged spin-table root nearest below a = 0.999 to 1e-8, which ties the a = 1
# root to the same overtone label. Returns nothing for any other mode or spin.
const _NEAR_EXTREMAL_EXPANSION_MINIMUM_SPIN = 0.999

function _near_extremal_expansion_data(mode::QNMMode)
    key = (mode.s, mode.l, mode.m, mode.n)
    extremal = get(_EXTREMAL_ROOT_TABLE, key, nothing)
    (extremal === nothing || extremal[1] != :DM) && return nothing
    spin = get(_SPIN_ROOT_TABLE, key, nothing)
    spin === nothing && return nothing
    roots, inversions = spin
    node = findlast(u -> one(u) - u^2 <= _NEAR_EXTREMAL_EXPANSION_MINIMUM_SPIN, _SPIN_TABLE_U)
    node === nothing && return nothing
    (isfinite(roots[node]) && inversions[node] >= 0) || return nothing
    coefficients = _cached_near_extremal_coefficients(mode, extremal[2])
    coefficients === nothing && return nothing
    xmax = (2pi * max(real(2 * coefficients[1] - mode.m), 0.0) / 36.8)^2
    unode = _SPIN_TABLE_U[node]
    xnode = (unode^2) * (2 - unode^2)
    xnode <= xmax || return nothing
    abs(evalpoly(xnode, coefficients) - roots[node]) <=
        1.0e-8 * abs(roots[node]) || return nothing
    return (inversions[node], coefficients, xmax)
end

function _kappa_squared_table_root(mode::QNMMode, a::Float64)
    mode.branch == :positive_real &&
        _NEAR_EXTREMAL_EXPANSION_MINIMUM_SPIN < a < 1 || return nothing
    entry = _near_extremal_expansion_data(mode)
    entry === nothing && return nothing
    inversion, coefficients, xmax = entry
    x = (1 - a) * (1 + a)
    x <= xmax || return nothing
    omega = evalpoly(x, coefficients)
    detuning = real(2 * coefficients[1] - mode.m)
    bound = max(abs(coefficients[end]) * x^length(coefficients),
        exp(-2pi * detuning / sqrt(x))) / abs(omega)
    angular_order = 32
    angular = angular_branch(mode, a, omega; truncation_order=angular_order,
        sheet_id=:qnm_kappa_squared_table)
    angular.status in (:continued, :predictor_corrected, :spherical_anchor,
        :high_precision_refined) || return nothing
    contract = _schwarzschild_label_contract(mode, :overtone)
    return LeaverResult(
        mode, :overtone, mode.n, inversion, a, omega,
        angular.angular_A, angular.lambda, angular.mixing,
        complex(NaN, NaN), NaN, 0, bound, Float64(angular.residual), 53,
        :accepted, :accepted_near_extremal_expansion,
        (
            input_type=string(typeof(a)),
            working_type="Float64",
            guess_source=:packaged_kappa_squared_table,
            convention=:overtone,
            overtone_index=mode.n,
            inversion_index=inversion,
            label_convention=contract.label_convention,
            branch_label=contract.branch_label,
            endpoint_kind=contract.endpoint_kind,
            near_as_role=contract.near_as_role,
            schwarzschild_sequence_index=contract.schwarzschild_sequence_index,
            schwarzschild_multiplicity=contract.schwarzschild_multiplicity,
            sequence_segment=0,
            family=:DM,
            angular_route=angular.metadata.route,
            angular_sheet=angular.sheet_id,
            angular_order,
            root_drift=0.0,
            root_equation=:near_extremal_expansion_in_x,
            acceptance_basis=:near_extremal_expansion_validity,
            kappa_squared=x,
            expansion_validity_x_max=xmax,
            relative_frequency_bound=bound,
            claim_boundary=:near_extremal_expansion_root,
        ),
    )
end

# Root for 0.9995 < |a| < 1 with a two-end anchored seed. A Float64 solve from the seed
# is tried first; when it is not kept, the default solve is returned if it agrees with the
# seed, and otherwise the seeded solve is repeated with the automatic BigFloat fallback,
# so that the seed fixes the label while the default path keeps its cost where it agrees.
function _two_end_anchored_frequency(mode::QNMMode, a, spin_seed; kwargs...)
    fallback = get(kwargs, :automatic_precision_fallback, true)
    quick = _qnm_frequency_labeled(mode, a; kwargs..., spin_seed,
        automatic_precision_fallback=false)
    # A failed solve records no guess source; the seed was then used.
    get(quick.provenance, :guess_source, spin_seed.source) == spin_seed.source ||
        return quick.status == :accepted || !fallback ? quick :
            _qnm_frequency_labeled(mode, a; kwargs...)
    _spin_seed_consistent(quick, spin_seed) && return quick
    default = _qnm_frequency_labeled(mode, a; kwargs...)
    (_spin_seed_consistent(default, spin_seed) || !fallback) && return default
    seeded = _qnm_frequency_labeled(mode, a; kwargs..., spin_seed)
    return _spin_seed_consistent(seeded, spin_seed) ? seeded : default
end

# A tabulated seed is kept when the polished root is accepted with the tabulated
# inversion index and lies within the seed's acceptance radius: ten error estimates (at
# least 1e-6 relative) for an interpolated node seed, a quarter of the surface gravity
# for a two-end anchored one.
_spin_seed_consistent(result, seed) = result.status == :accepted &&
    result.inversion_index == seed.inversion &&
    abs(result.omega - seed.seed) <= seed.radius

# Root at a = 1 from the packaged table of exact-extremal roots (ExtremalRootTable.jl),
# which were computed with the kappa ladder below. A ZDM entry returns the synchronous
# endpoint; a DM entry seeds the exact-extremal recurrence, and its root is kept only
# when the recurrence is accepted and reproduces the tabulated frequency to 1e-10.
# Otherwise (or for modes, conventions and guesses outside the table) returns nothing
# and the ladder runs.
function _tabulated_extremal_root(mode::QNMMode, a::T; guess, convention::Symbol,
        inversion_index, angular_order::Int, sheet_id::Symbol) where {T<:AbstractFloat}
    (T === Float64 && a == one(T) && guess === nothing &&
        mode.branch == :positive_real && convention == :overtone &&
        inversion_index === :auto) || return nothing
    entry = get(_EXTREMAL_ROOT_TABLE, (mode.s, mode.l, mode.m, mode.n), nothing)
    entry === nothing && return nothing
    family, seed, cf_error, cf_iterations = entry
    if family == :ZDM
        endpoint = _extremal_zdm_endpoint_result(mode, a; angular_order,
            sheet_id, convention, requested_inversion_index=inversion_index)
        return LeaverResult(endpoint.mode, endpoint.convention,
            endpoint.overtone_index, endpoint.inversion_index, endpoint.a,
            endpoint.omega, endpoint.angular_A, endpoint.lambda,
            endpoint.mixing, endpoint.cf_value, endpoint.cf_error,
            endpoint.cf_iterations, endpoint.root_residual,
            endpoint.angular_residual, endpoint.precision_bits,
            endpoint.status, endpoint.stop_reason,
            merge(endpoint.provenance, (family_classification=:ZDM,
                extremal_seed_source=:packaged_table)))
    end
    exact = extremal_qnm_root(mode, complex(T(real(seed)), T(imag(seed)));
        angular_order=angular_order + 8,
        sheet_id=Symbol(sheet_id, :_exact_extremal_root))
    angular = exact.angular
    angular_usable = angular !== nothing && angular.status in (
        :continued, :predictor_corrected, :spherical_anchor,
        :high_precision_refined)
    (exact.status == :accepted && angular_usable &&
        abs(exact.omega - seed) <= T(1e-10) * max(one(T), abs(seed))) ||
        return nothing
    contract = _schwarzschild_label_contract(mode, convention)
    nan_complex = complex(T(NaN), T(NaN))
    return LeaverResult(
        mode, convention, mode.n, mode.n, a, exact.omega,
        angular.angular_A, angular.lambda, angular.mixing, nan_complex,
        T(cf_error), cf_iterations, exact.residual, T(angular.residual),
        _precision_bits(T),
        :accepted, :accepted_exact_extremal_root,
        (
            family_classification=:DM,
            input_type=string(typeof(a)),
            working_type=string(T),
            label_convention=contract.label_convention,
            branch_label=contract.branch_label,
            endpoint_kind=contract.endpoint_kind,
            extremal_endpoint=true,
            endpoint_method=:exact_extremal_recurrence,
            spectrum_object=:discrete_damped_pole,
            simple_pole=true,
            radial_observables=:exact_extremal_gsn_scattering_basis,
            gsn_amplitude_limit_implemented=true,
            exact_extremal_depth=exact.depth,
            exact_extremal_depth_drift=exact.depth_drift,
            exact_extremal_depth_rows=exact.depth_rows,
            exact_extremal_seed_shift=exact.seed_shift,
            extremal_seed_source=:packaged_table,
            angular_route=angular.metadata.route,
            angular_sheet=angular.sheet_id,
            angular_order=angular_order + 8,
            root_equation=:richartz_extremal_recurrence,
            family=:DM,
            claim_boundary=:exact_extremal_frequency_root,
        ),
    )
end

# Root at abs(a) = 1: solve a near-extremal kappa ladder, classify the family,
# and return the ZDM endpoint or the exact-extremal damped root.
# Damped extremal roots from moderate spins. A damped frequency is analytic in
# x = 1 - a^2 at a = 1, so roots at |a| = 0.9 ... 0.99 (one continuation to 0.9,
# then polishes seeded by the polynomial predictor) extrapolate to x = 0 well
# enough to seed the exact-extremal recurrence. None of these spins needs the
# deep continued fractions of the near-extremal kappa ladder. Returns `nothing`
# (and the caller runs the ladder) unless every source root is accepted, the
# exact root lies within 1e-4 of the extrapolation, and Re(2 omega - m a) is
# clearly positive (zero-damping families approach m/2 like sqrt(x), so their
# extrapolation in x does not land on a damped root).
function _extremal_moderate_spin_result(mode::QNMMode, a::T;
        angular_order::Int, sheet_id::Symbol, cf_tolerance,
        cf_minimum_iterations::Int, cf_maximum_iterations::Int,
        root_tolerance, stability_tolerance, maximum_root_iterations::Int,
        trust_radius, automatic_precision_fallback::Bool,
        fallback_precision_bits::Int) where {T<:AbstractFloat}
    spins = T[copysign(T(s), a) for s in (0.9, 0.94, 0.96, 0.98, 0.99)]
    xs = T[one(T) - spin^2 for spin in spins]
    roots = Complex{T}[]
    rows = NamedTuple[]
    selected_inversion = :auto
    for (index, spin) in pairs(spins)
        source_guess = isempty(roots) ? nothing :
            _polynomial_predictor(xs[1:length(roots)], roots, xs[index])
        result = _qnm_frequency_labeled(
            mode, spin;
            guess=source_guess,
            convention=:overtone,
            inversion_index=selected_inversion,
            angular_order,
            sheet_id=Symbol(sheet_id, :_extremal_moderate_, index),
            cf_tolerance,
            cf_minimum_iterations,
            cf_maximum_iterations,
            root_tolerance,
            stability_tolerance,
            maximum_root_iterations,
            trust_radius,
            automatic_precision_fallback,
            fallback_precision_bits,
        )
        result.status == :accepted || return nothing
        omega = _complex_of_type(T, result.omega)
        push!(roots, omega)
        push!(rows, (
            a=spin, one_minus_a_squared=xs[index], omega,
            root_residual=T(result.root_residual),
            cf_error=T(result.cf_error),
            cf_iterations=result.cf_iterations,
            precision_bits=result.precision_bits,
            inversion_index=result.inversion_index,
            status=result.status,
            stop_reason=result.stop_reason,
        ))
        if selected_inversion === :auto && result.inversion_index isa Integer
            selected_inversion = result.inversion_index
        end
    end
    seed = _polynomial_predictor(xs, roots, zero(T))
    lower = _polynomial_predictor(xs[2:end], roots[2:end], zero(T))
    exact = extremal_qnm_root(mode, seed;
        angular_order=angular_order + 8,
        sheet_id=Symbol(sheet_id, :_exact_extremal_root))
    angular = exact.angular
    angular_usable = angular !== nothing && angular.status in (
        :continued, :predictor_corrected, :spherical_anchor,
        :high_precision_refined)
    detuning = real(2 * exact.omega - T(mode.m) * sign(a))
    (exact.status == :accepted && angular_usable &&
        exact.seed_shift <= T(1e-4) * max(one(T), abs(exact.omega)) &&
        detuning > T(0.02) && abs(real(exact.omega)) > T(1e-6)) ||
        return nothing
    contract = _schwarzschild_label_contract(mode, :overtone)
    return LeaverResult(
        mode, :overtone, mode.n, selected_inversion, a, exact.omega,
        angular.angular_A, angular.lambda, angular.mixing,
        complex(T(NaN), T(NaN)),
        maximum(row.cf_error for row in rows),
        maximum(row.cf_iterations for row in rows), exact.residual,
        T(angular.residual), _precision_bits(T),
        :accepted, :accepted_exact_extremal_root,
        (
            source_rows=rows,
            extrapolation_coordinate=:one_minus_a_squared,
            maximum_source_precision_bits=
                maximum(row.precision_bits for row in rows),
            quartic_limit=seed,
            cubic_limit=lower,
            limit_drift=abs(seed - lower) / max(one(T), abs(seed)),
            all_sources_accepted=true,
            family_classification=:DM,
            extremal_seed_source=:moderate_spin_extrapolation,
            input_type=string(typeof(a)),
            working_type=string(T),
            label_convention=contract.label_convention,
            branch_label=contract.branch_label,
            endpoint_kind=contract.endpoint_kind,
            extremal_endpoint=true,
            endpoint_method=:exact_extremal_recurrence,
            spectrum_object=:discrete_damped_pole,
            simple_pole=true,
            radial_observables=:exact_extremal_gsn_scattering_basis,
            gsn_amplitude_limit_implemented=true,
            exact_extremal_depth=exact.depth,
            exact_extremal_depth_drift=exact.depth_drift,
            exact_extremal_depth_rows=exact.depth_rows,
            exact_extremal_seed_shift=exact.seed_shift,
            exact_extremal_detuning=detuning,
            angular_route=angular.metadata.route,
            angular_sheet=angular.sheet_id,
            angular_order=angular_order + 8,
            root_equation=:richartz_extremal_recurrence,
            family=:DM,
            claim_boundary=:exact_extremal_frequency_root,
        ),
    )
end

function _extremal_damped_limit_result(mode::QNMMode, a::T;
        guess,
        convention::Symbol, inversion_index,
        angular_order::Int, sheet_id::Symbol,
        cf_tolerance, cf_minimum_iterations::Int,
        cf_maximum_iterations::Int, root_tolerance,
        stability_tolerance, maximum_root_iterations::Int,
        trust_radius, automatic_precision_fallback::Bool,
        fallback_precision_bits::Int) where {T<:AbstractFloat}
    tabulated = _tabulated_extremal_root(mode, a; guess, convention,
        inversion_index, angular_order, sheet_id)
    tabulated === nothing || return tabulated
    if guess === nothing && convention == :overtone && inversion_index === :auto
        fast = _extremal_moderate_spin_result(mode, a; angular_order,
            sheet_id, cf_tolerance, cf_minimum_iterations,
            cf_maximum_iterations, root_tolerance, stability_tolerance,
            maximum_root_iterations, trust_radius,
            automatic_precision_fallback, fallback_precision_bits)
        fast === nothing || return fast
    end
    kappas = T[T(0.04), T(0.02), T(0.01), T(0.005), T(0.0025)]
    spins = T[copysign(sqrt(one(T) - kappa^2), a)
        for kappa in kappas]
    spin_defects = T[one(T) - abs(spin) for spin in spins]
    roots = Complex{T}[]
    rows = NamedTuple[]
    previous = guess === nothing ? nothing : _complex_of_type(T, guess)
    selected_inversion = inversion_index

    for (index, spin) in pairs(spins)
        # Without a caller guess every ladder spin is reached by its own
        # continuation in nu from the Schwarzschild label. Chaining polished
        # predictions between ladder spins hops between the crowded
        # zero-damping levels.
        source_guess = guess === nothing ? nothing :
            length(roots) >= 2 ? _polynomial_predictor(
                spin_defects[max(1, length(roots) - 2):length(roots)],
                roots[max(1, length(roots) - 2):end],
                spin_defects[index]) : previous
        result = _qnm_frequency_labeled(
            mode, spin;
            guess=source_guess,
            convention,
            inversion_index=selected_inversion,
            continuation_coordinate=:nu,
            angular_order,
            sheet_id=Symbol(sheet_id, :_extremal_limit_, index),
            cf_tolerance,
            cf_minimum_iterations,
            cf_maximum_iterations,
            root_tolerance,
            stability_tolerance,
            maximum_root_iterations,
            trust_radius,
            automatic_precision_fallback,
            fallback_precision_bits,
        )
        initial_failure = result.status == :accepted ? nothing : result
        source_precision_fallback = false
        source_usable = result.status == :accepted ||
            (mode.n >= 5 && result.status == :estimated)
        # A source sequence that ended on the negative imaginary axis is not revived by
        # higher precision (as for finite spins, segment 0 is not restarted).
        if !source_usable && T === Float64 &&
                result.stop_reason != :negative_imaginary_axis_endpoint
            bits = max(192, fallback_precision_bits)
            source_precision_fallback = true
            # The continued fraction of the high-precision fallback may use up to 600/kappa
            # terms (the tail converges more slowly as kappa -> 0); its tolerance is unchanged.
            # The high-precision fallback first polishes from the extrapolation of
            # the ladder roots already accepted, instead of repeating the
            # continuation from the Schwarzschild label at the deepest
            # continued-fraction depth. The polished root is kept only within a
            # quarter of the last ladder step of that prediction (the same level);
            # otherwise the continuation below runs as before.
            ladder_prediction = guess === nothing && length(roots) >= 2 ?
                _polynomial_predictor(
                    spin_defects[max(1, length(roots) - 2):length(roots)],
                    roots[max(1, length(roots) - 2):end],
                    spin_defects[index]) : nothing
            result = nothing
            if ladder_prediction !== nothing
                predicted = setprecision(bits) do
                    kappa_big = BigFloat(kappas[index])
                    spin_big = copysign(
                        sqrt(one(BigFloat) - kappa_big^2), BigFloat(a))
                    _qnm_frequency_labeled(
                        mode, spin_big;
                        guess=_complex_of_type(BigFloat, ladder_prediction),
                        convention,
                        inversion_index=selected_inversion,
                        continuation_coordinate=:nu,
                        angular_order=angular_order + 8,
                        sheet_id=Symbol(
                            sheet_id, :_extremal_limit_high_precision_, index),
                        cf_tolerance=BigFloat(cf_tolerance),
                        cf_minimum_iterations=max(cf_minimum_iterations, 600),
                        cf_maximum_iterations=max(cf_maximum_iterations, 60000,
                            ceil(Int, 600 / kappas[index])),
                        root_tolerance=min(BigFloat(root_tolerance), big"1e-14"),
                        stability_tolerance=min(
                            BigFloat(stability_tolerance), big"1e-14"),
                        maximum_root_iterations=max(maximum_root_iterations, 60),
                        trust_radius=min(BigFloat(trust_radius), big"0.02"),
                        automatic_precision_fallback=false,
                        fallback_precision_bits=bits,
                    )
                end
                ladder_step = abs(roots[end] - roots[end-1])
                same_level = abs(_complex_of_type(T, predicted.omega) -
                    ladder_prediction) <= ladder_step / 4
                predicted_usable = predicted.status == :accepted ||
                    (mode.n >= 5 && predicted.status == :estimated)
                predicted_usable && same_level && (result = predicted)
            end
            result === nothing && (result = setprecision(bits) do
                kappa_big = BigFloat(kappas[index])
                spin_big = copysign(
                    sqrt(one(BigFloat) - kappa_big^2), BigFloat(a))
                guess_big = source_guess === nothing ? nothing :
                    _complex_of_type(BigFloat, source_guess)
                _qnm_frequency_labeled(
                    mode, spin_big;
                    guess=guess_big,
                    convention,
                    inversion_index=selected_inversion,
                    continuation_coordinate=:nu,
                    angular_order=angular_order + 8,
                    sheet_id=Symbol(
                        sheet_id, :_extremal_limit_high_precision_, index),
                    # Preserve the caller's continued-fraction tail gate.  The
                    # endpoint fallback raises arithmetic precision so that
                    # the root can be polished; it must not silently replace
                    # a 1e-13 tail contract by an unattainable 1e-18 contract.
                    cf_tolerance=BigFloat(cf_tolerance),
                    cf_minimum_iterations=max(cf_minimum_iterations, 600),
                    cf_maximum_iterations=max(cf_maximum_iterations, 60000,
                            ceil(Int, 600 / kappas[index])),
                    root_tolerance=min(BigFloat(root_tolerance), big"1e-14"),
                    stability_tolerance=min(
                        BigFloat(stability_tolerance), big"1e-14"),
                    maximum_root_iterations=max(maximum_root_iterations, 60),
                    trust_radius=min(BigFloat(trust_radius), big"0.02"),
                    automatic_precision_fallback=false,
                    fallback_precision_bits=bits,
                )
            end)
        end
        if !(result.status == :accepted ||
                (mode.n >= 5 && result.status == :estimated)) &&
                _extremal_zdm_endpoint_scope(mode, a)
            endpoint = _extremal_zdm_endpoint_result(mode, a; angular_order,
                sheet_id, convention,
                requested_inversion_index=inversion_index)
            return LeaverResult(endpoint.mode, endpoint.convention,
                endpoint.overtone_index, endpoint.inversion_index,
                endpoint.a, endpoint.omega, endpoint.angular_A,
                endpoint.lambda, endpoint.mixing, endpoint.cf_value,
                endpoint.cf_error, endpoint.cf_iterations,
                endpoint.root_residual, endpoint.angular_residual,
                endpoint.precision_bits, :estimated,
                :extremal_family_unverified,
                merge(endpoint.provenance, (
                    family_classification=:unverified_source_failure,
                    failed_spin=spin,
                    failed_source_reason=result.stop_reason,
                    source_rows=rows,
                )))
        end
        (result.status == :accepted ||
            (mode.n >= 5 && result.status == :estimated)) ||
            return _failed_result(
            mode, a, isempty(roots) ? complex(T(NaN), T(NaN)) : roots[end],
            :extremal_limit_source_failure,
            (
                extremal_endpoint=true,
                endpoint_method=:spin_defect_richardson_extrapolation,
                failed_spin=spin,
                failed_source_reason=result.stop_reason,
                initial_source_reason=initial_failure === nothing ?
                    :not_applicable : initial_failure.stop_reason,
                source_rows=rows,
            );
            convention,
            overtone_index=convention == :leaver ? missing : mode.n,
            inversion_index=selected_inversion === :auto ? mode.n :
                Int(selected_inversion),
        )
        omega = _complex_of_type(T, result.omega)
        push!(roots, omega)
        push!(rows, (
            kappa=kappas[index],
            spin_defect=spin_defects[index],
            a=spin,
            omega,
            root_residual=T(result.root_residual),
            cf_error=T(result.cf_error),
            cf_iterations=result.cf_iterations,
            precision_bits=result.precision_bits,
            inversion_index=result.inversion_index,
            status=result.status,
            stop_reason=result.stop_reason,
            source_precision_fallback,
        ))
        previous = omega
        if selected_inversion === :auto &&
                result.inversion_index isa Integer
            selected_inversion = result.inversion_index
        end
    end

    cubic = _polynomial_predictor(
        spin_defects[2:end], roots[2:end], zero(T))
    quadratic = _polynomial_predictor(
        spin_defects[end-2:end], roots[end-2:end], zero(T))
    drift = abs(cubic - quadratic) /
        max(one(T), abs(cubic), abs(quadratic))
    contract = _schwarzschild_label_contract(mode, convention)
    selected_output_inversion = selected_inversion === :auto ? mode.n :
        Int(selected_inversion)
    all_sources_accepted = all(row.status == :accepted for row in rows)
    family = _extremal_ladder_family(mode, a, roots)
    ladder = (
        source_rows=rows,
        extrapolation_coordinate=:one_minus_abs_a,
        maximum_source_precision_bits=
            maximum(row.precision_bits for row in rows),
        quadratic_limit=quadratic,
        cubic_limit=cubic,
        limit_drift=drift,
        all_sources_accepted,
        family_classification=family.family,
        family_ratios=family.ratios,
    )

    if family.family == :ZDM
        endpoint = _extremal_zdm_endpoint_result(mode, a; angular_order,
            sheet_id, convention,
            requested_inversion_index=inversion_index)
        return LeaverResult(endpoint.mode, endpoint.convention,
            endpoint.overtone_index, endpoint.inversion_index, endpoint.a,
            endpoint.omega, endpoint.angular_A, endpoint.lambda,
            endpoint.mixing, endpoint.cf_value, endpoint.cf_error,
            endpoint.cf_iterations, endpoint.root_residual,
            endpoint.angular_residual, endpoint.precision_bits,
            endpoint.status, endpoint.stop_reason,
            merge(endpoint.provenance, ladder))
    end

    # Damped family: solve the exact-extremal condition at a = 1, seeded by
    # the analytic-in-(1 - |a|) limit of the near-extremal ladder.
    exact = extremal_qnm_root(mode, quadratic;
        angular_order=angular_order + 8,
        sheet_id=Symbol(sheet_id, :_exact_extremal_root))
    # Same branch as the ladder: the remaining distance from the smallest
    # spin defect must not exceed the last ladder step.
    last_step = abs(roots[end] - roots[end-1])
    continuity = abs(exact.omega - roots[end]) <=
        2 * last_step + sqrt(eps(T))
    angular = exact.angular
    angular_usable = angular !== nothing && angular.status in (
        :continued, :predictor_corrected, :spherical_anchor,
        :high_precision_refined)
    exact_ok = exact.status == :accepted && angular_usable
    accepted = exact_ok && continuity && family.family == :DM &&
        all(row.status in (:accepted, :estimated) for row in rows)
    estimated = !accepted && angular_usable && isfinite(exact.omega)
    nan_complex = complex(T(NaN), T(NaN))
    omega = angular_usable ? exact.omega : quadratic
    return LeaverResult(
        mode, convention, convention == :leaver ? missing : mode.n,
        selected_output_inversion, a, omega,
        angular_usable ? angular.angular_A : nan_complex,
        angular_usable ? angular.lambda : nan_complex,
        angular_usable ? angular.mixing : nothing, nan_complex,
        maximum(row.cf_error for row in rows),
        maximum(row.cf_iterations for row in rows), exact.residual,
        angular_usable ? T(angular.residual) : T(Inf), _precision_bits(T),
        accepted ? :accepted : estimated ? :estimated : :failed,
        accepted ? :accepted_exact_extremal_root :
            !continuity ? :extremal_root_off_ladder_branch :
            family.family != :DM ? :extremal_family_ambiguous :
            exact.stop_reason,
        merge(ladder, (
            input_type=string(typeof(a)),
            working_type=string(T),
            label_convention=contract.label_convention,
            branch_label=contract.branch_label,
            endpoint_kind=contract.endpoint_kind,
            extremal_endpoint=true,
            endpoint_method=:exact_extremal_recurrence,
            spectrum_object=:discrete_damped_pole,
            simple_pole=true,
            radial_observables=:exact_extremal_gsn_scattering_basis,
            gsn_amplitude_limit_implemented=true,
            exact_extremal_depth=exact.depth,
            exact_extremal_depth_drift=exact.depth_drift,
            exact_extremal_depth_rows=exact.depth_rows,
            exact_extremal_seed_shift=exact.seed_shift,
            ladder_continuity=continuity,
            ladder_last_step=last_step,
            angular_route=angular_usable ? angular.metadata.route : :failed,
            angular_sheet=angular_usable ? angular.sheet_id : :failed,
            angular_order=angular_order + 8,
            root_equation=:richartz_extremal_recurrence,
            family=:DM,
            claim_boundary=:exact_extremal_frequency_root,
        )),
    )
end

function _extremal_ladder_family(mode::QNMMode, a, roots)
    # Zero-damping branches approach the synchronous frequency m a / 2 linearly
    # in the surface gravity, so their distance halves along the kappa ladder;
    # damped branches converge to a separate limit.
    target = mode.m * a / 2
    distances = [abs(root - target) for root in roots[end-2:end]]
    ratios = (distances[2] / distances[1], distances[3] / distances[2])
    synchronous = mode.m * a > 0 && (mode.branch == :positive_real ?
        real(target) > 0 : real(target) < 0)
    family = !synchronous ? :DM :
        all(r -> r <= oftype(r, 0.7), ratios) ? :ZDM :
        all(r -> r >= oftype(r, 0.8), ratios) ? :DM : :ambiguous
    return (family=family, ratios=ratios)
end

# Map a root solved at (-a, -m) to the requested (a, m); the frequency is
# unchanged and the angular data are recomputed.
function _spin_reflected_result(
        requested_mode::QNMMode, requested_a, source::LeaverResult;
        angular_order::Int, sheet_id::Symbol)
    if hasproperty(source.provenance, :near_axis_position_recovery) &&
            source.omega isa Complex{BigFloat} &&
            precision(BigFloat) < source.precision_bits
        return setprecision(BigFloat, source.precision_bits) do
            _spin_reflected_result(requested_mode, requested_a, source;
                angular_order, sheet_id)
        end
    end
    mapped_a = -source.a
    angular = angular_branch(
        requested_mode, mapped_a, source.omega;
        truncation_order=angular_order,
        sheet_id=Symbol(sheet_id, :_spin_reflection_guard),
    )
    scale = max(one(real(abs(angular.angular_A))),
        abs(angular.angular_A), abs(source.angular_A))
    angular_drift = abs(angular.angular_A - source.angular_A) / scale
    return LeaverResult(
        requested_mode,
        source.convention,
        source.overtone_index,
        source.inversion_index,
        mapped_a,
        source.omega,
        angular.angular_A,
        angular.lambda,
        angular.mixing,
        source.cf_value,
        source.cf_error,
        source.cf_iterations,
        source.root_residual,
        typeof(mapped_a)(angular.residual),
        source.precision_bits,
        source.status,
        source.stop_reason,
        merge(source.provenance, (
            spin_reflection_applied=true,
            spin_reflection_map=:a_m_to_minus_a_minus_m,
            spin_reflection_source_a=source.a,
            spin_reflection_source_m=source.mode.m,
            spin_reflection_target_a=mapped_a,
            spin_reflection_target_m=requested_mode.m,
            spin_reflection_frequency_map=:identity,
            spin_reflection_angular_drift=angular_drift,
            angular_route=angular.metadata.route,
            angular_sheet=angular.sheet_id,
            angular_order=angular.truncation_order,
        )),
    )
end

# Map an ordinary root of (s, l, -m, n) at spin a to the mirror root of (s, l, m, n)
# at the same spin: omega -> -conj(omega) (the Teukolsky equation is invariant under
# m -> -m with complex conjugation); the angular data are recomputed at the mapped
# frequency. A zero-damping endpoint maps to the mirror endpoint omega = m a / 2.
function _mirror_reflected_result(
        mirror_mode::QNMMode, source::LeaverResult;
        angular_order::Int, sheet_id::Symbol)
    omega = -conj(source.omega)
    angular = angular_branch(
        mirror_mode, source.a, omega;
        truncation_order=angular_order,
        sheet_id=Symbol(sheet_id, :_mirror_reflection),
    )
    endpoint_fields = hasproperty(source.provenance, :synchronous_frequency) ?
        (synchronous_frequency=omega,) : NamedTuple()
    return LeaverResult(
        mirror_mode,
        source.convention,
        source.overtone_index,
        source.inversion_index,
        source.a,
        omega,
        angular.angular_A,
        angular.lambda,
        angular.mixing,
        conj(source.cf_value),
        source.cf_error,
        source.cf_iterations,
        source.root_residual,
        typeof(real(omega))(angular.residual),
        source.precision_bits,
        source.status,
        source.stop_reason,
        merge(source.provenance, endpoint_fields, (
            mirror_reflection_applied=true,
            mirror_reflection_map=:omega_to_minus_conj_omega_m_to_minus_m,
            mirror_reflection_source_m=source.mode.m,
            angular_route=angular.metadata.route,
            angular_sheet=angular.sheet_id,
            angular_order=angular.truncation_order,
        )),
    )
end

"""
    qnm_frequency(a, s, l, m, n, branch=ordinary; convention=:overtone,
        segment=0, guess=nothing, kwargs...)
    qnm_frequency(mode::QNMMode, a; kwargs...)

Return the QNM frequency of the mode `(s, l, m, n)` at spin `a` as a
`LeaverResult` (root only, no amplitudes). `branch` is `ordinary` (default) or
`mirror`. Without `guess`, the mirror root is the reflection `-conj(ω)` of the
ordinary `(s, l, -m, n)` root at the same spin, with the angular data
recomputed and the status of that ordinary root. Only when that ordinary root
is not accepted (and did not stop at the negative imaginary axis) is the mirror
root polished on the negative-real branch, starting from the reflection. With
`guess`, the mirror root is polished from `guess`. The `QNMMode` form is also
accepted; there a mirror (`:negative_real`) mode needs an explicit `guess`.

The root solves the coupled angular and radial Leaver equations in `M=1` units.
Without an explicit guess it starts from the packaged root table or, for a mode
or spin the table lacks, from the labeled Schwarzschild root continued in spin;
the final root is re-polished with a higher angular truncation as a check. For
some damped modes at `0.999 < |a| < 1` the root is instead evaluated from a
near-extremal expansion, and at `|a| = 1` it is solved with an exact-extremal
recurrence (see the QNM manual page). `status` is `:accepted`, `:estimated` or
`:failed`, with the reason in `stop_reason`.

With the default `convention=:overtone`, overtone labels follow the Python
`qnm` package and Cook & Zalutskiy: overtone `n` is the Kerr continuation of
the `n`-th Schwarzschild root. For `s=-2`, `l=2`, `|m|=2` co-rotating modes the
algebraically special frequency `-2i` splits into a pair (Cook's 8₀ and 8₁);
`n=8` returns the less damped 8₀, and 8₁ is `n=9` with
`convention=:complete_spectrum`, which numbers every root of that family in
order of damping (overtone `n >= 9` is complete-spectrum `n+1`).
`convention=:leaver` only selects the continued-fraction inversion index
(`inversion_index`); it does not label a mode, requires an explicit `guess`,
and returns `overtone_index = missing`.

Some sequences end on the negative imaginary axis (NIA) at finite spin and
re-emerge from it at larger spin (Cook & Zalutskiy's `n₀` and `n₁`). The
default `segment=0` continues the sequence attached to `a=0`; when it ends on
the NIA the result fails with `stop_reason=:negative_imaginary_axis_endpoint`
and records the extrapolated `endpoint_spin`. `segment=1` returns the
re-emerged `n₁` (implemented for `m=0` only), and `segment=:auto` returns the first
segment that exists at `a`. See the QNM manual page for the stop reasons.
"""
function qnm_frequency(mode::QNMMode, a; segment=0, kwargs...)
    selected_segment = _qnm_segment(segment)
    convention = _qnm_convention(get(kwargs, :convention, :overtone))
    internal_n = _overtone_internal_index(mode, a, convention)
    internal_mode = internal_n == mode.n ? mode :
        QNMMode(mode.s, mode.l, mode.m, internal_n, mode.branch)
    # keywords at their default values (as passed on by `qnm`) keep the near-extremal expansion
    default_keyword(k) = k === :convention ||
        (k === :guess && kwargs[k] === nothing) ||
        (k === :inversion_index && kwargs[k] === :auto)
    table_root = selected_segment === 0 && a isa Float64 &&
        all(default_keyword, keys(kwargs)) &&
        _qnm_convention(get(kwargs, :convention, :overtone)) == :overtone ?
        _kappa_squared_table_root(a < 0 ? QNMMode(internal_mode.s,
            internal_mode.l, -internal_mode.m, internal_mode.n,
            internal_mode.branch) : internal_mode, abs(a)) : nothing
    result = if table_root !== nothing
        a < 0 ? _spin_reflected_result(internal_mode, a, table_root;
            angular_order=32, sheet_id=:qnm_kappa_squared_table) : table_root
    elseif selected_segment === 0
        spin_seed = get(kwargs, :guess, nothing) === nothing ?
            _tabulated_spin_seed(a isa Real && a < 0 ?
                QNMMode(internal_mode.s, internal_mode.l, -internal_mode.m,
                    internal_mode.n, internal_mode.branch) : internal_mode,
                a isa Real ? abs(a) : a) : nothing
        if spin_seed !== nothing && spin_seed.source == :packaged_two_end_anchor
            _two_end_anchored_frequency(internal_mode, a, spin_seed; kwargs...)
        else
            seeded = _qnm_frequency_labeled(internal_mode, a; spin_seed, kwargs...)
            spin_seed === nothing ||
                get(seeded.provenance, :guess_source, nothing) != spin_seed.source ||
                _spin_seed_consistent(seeded, spin_seed) ? seeded :
                _qnm_frequency_labeled(internal_mode, a; kwargs...)
        end
    else
        _qnm_frequency_segmented(internal_mode, a, selected_segment;
            kwargs...)
    end
    result = _with_sequence_segment(result, 0)
    internal_n == mode.n && return result
    return _relabel_result(result, mode)
end

function qnm_frequency(a::Real, s::Integer, l::Integer, m::Integer,
        n::Integer, branch=ordinary; guess=nothing, kwargs...)
    if _qnm_selector(branch) == mirror && guess === nothing
        options = Base.structdiff(values(kwargs),
            NamedTuple{(:convention, :inversion_index, :segment)})
        return first(_mirror_qnm_root(a, s, l, m, n, nothing, options;
            convention=get(kwargs, :convention, :overtone),
            inversion_index=get(kwargs, :inversion_index, :auto),
            segment=get(kwargs, :segment, 0)))
    end
    return qnm_frequency(_qnm_flat_mode(s, l, m, n, branch), a; guess,
        kwargs...)
end

_qnm_flat_mode(s, l, m, n, branch) = QNMMode(Int(s), Int(l), Int(m), Int(n),
    _qnm_selector(branch) == mirror ? :negative_real : :positive_real)

function _qnm_segment(segment)
    if segment isa Union{Symbol,AbstractString}
        Symbol(lowercase(String(segment))) == :auto && return :auto
        throw(ArgumentError("segment must be 0, 1, or :auto."))
    end
    segment isa Integer || throw(ArgumentError(
        "segment must be 0, 1, or :auto."))
    segment >= 0 || throw(ArgumentError("segment must be nonnegative."))
    segment <= 1 || throw(ArgumentError(
        "Only segments 0 and 1 are implemented."))
    return Int(segment)
end

# Paths that never enter the spin continuation (Schwarzschild endpoints,
# extremal limits) are segment 0 unless they say otherwise.
function _with_sequence_segment(result::LeaverResult, segment::Int)
    hasproperty(result.provenance, :sequence_segment) && return result
    return _with_provenance(result, (sequence_segment=segment,))
end
_with_sequence_segment(result, segment::Int) = result

function _with_provenance(result::LeaverResult, extra::NamedTuple;
        status=result.status, stop_reason=result.stop_reason)
    return LeaverResult(result.mode, result.convention,
        result.overtone_index, result.inversion_index, result.a,
        result.omega, result.angular_A, result.lambda, result.mixing,
        result.cf_value, result.cf_error, result.cf_iterations,
        result.root_residual, result.angular_residual,
        result.precision_bits, status, stop_reason,
        merge(result.provenance, extra))
end

# Public overtone label -> label of the internal complete-spectrum ordering. Only
# s=-2, l=2, |m|=2 co-rotating Kerr modes differ: internally the pair split from
# -2i is 8 (8_0) and 9 (8_1), and the Schwarzschild continuation n >= 9 is n+1.
function _overtone_internal_index(mode::QNMMode, a, convention::Symbol)
    l2m2 = mode.s == -2 && mode.l == 2 && abs(mode.m) == 2
    corotating = a * mode.m * (mode.branch == :positive_real ? 1 : -1) > 0
    (convention == :overtone && l2m2 && corotating && !iszero(a)) || return mode.n
    return mode.n >= 9 ? mode.n + 1 : mode.n
end

# Cook & Zalutskiy label the segment of a sequence that re-emerges from the NIA
# as n_1; segment 0 keeps the plain overtone label.
_segment_branch_label(n::Integer, segment::Integer) = segment == 0 ?
    Symbol(:kerr_overtone_n, n) : Symbol(:kerr_overtone_n, n, :_, segment)

function _relabel_result(result::LeaverResult, mode::QNMMode)
    label = mode.n == 8 ? result.provenance.branch_label :
        _segment_branch_label(mode.n,
            get(result.provenance, :sequence_segment, 0))
    provenance = hasproperty(result.provenance, :branch_label) ?
        merge(result.provenance, (branch_label=label,)) : result.provenance
    return LeaverResult(mode, result.convention,
        ismissing(result.overtone_index) ? missing : mode.n,
        result.inversion_index, result.a, result.omega, result.angular_A,
        result.lambda, result.mixing, result.cf_value, result.cf_error,
        result.cf_iterations, result.root_residual, result.angular_residual,
        result.precision_bits, result.status, result.stop_reason, provenance)
end
_relabel_result(result, mode::QNMMode) = result

# Spin continuation of one labeled root from the last accepted point to
# `sequence_target`. `state` holds the accepted path and the row logs and is
# extended in place. Returns `completed=false` with the last polish when the
# step falls below `minimum_step`.
function _continue_in_spin!(state, mode::QNMMode, sequence_target::T;
        selected_coordinate::Symbol, initial_step::T, minimum_step::T,
        maximum_step::T, continuation_root_tolerance::T,
        polish_options) where {T<:AbstractFloat}
    (; accepted_spins, accepted_coordinates, accepted_roots,
        accepted_angular_A, continuation_rows,
        rejected_continuation_rows) = state
    direction = sign(sequence_target - accepted_coordinates[end])
    proposed_step = min(initial_step,
        abs(sequence_target - accepted_coordinates[end]))
    while direction *
            (sequence_target - accepted_coordinates[end]) > zero(T)
        if length(accepted_coordinates) < 3
            proposed_step = initial_step
        else
            curvature = max(
                _continuation_curvature(
                    accepted_coordinates, accepted_roots),
                _continuation_curvature(
                    accepted_coordinates, accepted_angular_A),
            )
            if curvature > zero(T) && isfinite(curvature)
                proposed_step = clamp(T(0.05) / sqrt(curvature),
                    minimum_step, maximum_step)
            else
                proposed_step = min(proposed_step, maximum_step)
            end
        end
        remaining = abs(sequence_target - accepted_coordinates[end])
        attempted_step = min(proposed_step, remaining)
        attempt = 0
        accepted_step = false

        while !accepted_step
            attempt += 1
            next_coordinate = accepted_coordinates[end] +
                direction * attempted_step
            next_spin = selected_coordinate == :nu ?
                sin(next_coordinate) : next_coordinate
            predictor = _continuation_predictor(
                accepted_coordinates, accepted_roots, next_coordinate)
            polished = _polish_frequency_root(
                mode, next_spin, predictor;
                polish_options...)
            append!(continuation_rows, [merge(row, (
                a=next_spin, continuation_attempt=attempt,
                continuation_coordinate=selected_coordinate,
                coordinate_value=next_coordinate,
                predictor=predictor)) for row in polished.history])
            state.last_polish[] = polished
            branch_guard = polished.converged ?
                _continuation_branch_guard(
                    accepted_coordinates, accepted_roots,
                    next_coordinate,
                    predictor, polished.omega,
                    continuation_root_tolerance) :
                (accepted=false, discrepancy=T(Inf),
                    allowance=zero(T))

            # A converged root beyond an NIA crossing belongs to another
            # sequence even when the branch guard admits it.
            axis_crossing = polished.converged && branch_guard.accepted &&
                _negative_imaginary_axis_step_crosses_axis(mode, accepted_coordinates,
                    accepted_roots, next_coordinate)
            if polished.converged && branch_guard.accepted && !axis_crossing
                push!(accepted_spins, next_spin)
                push!(accepted_coordinates, next_coordinate)
                push!(accepted_roots, polished.omega)
                push!(accepted_angular_A,
                    polished.evaluation.angular.angular_A)
                proposed_step = attempted_step
                accepted_step = true
                continue
            end

            push!(rejected_continuation_rows, (
                a=next_spin,
                continuation_coordinate=selected_coordinate,
                coordinate_value=next_coordinate,
                attempted_step=attempted_step,
                attempt=attempt,
                predictor=predictor,
                candidate=polished.omega,
                converged=polished.converged,
                stop_reason=polished.stop_reason,
                branch_discrepancy=branch_guard.discrepancy,
                branch_allowance=branch_guard.allowance,
                axis_crossing,
            ))
            attempted_step /= 2
            if attempted_step < minimum_step
                reason = axis_crossing ? :nia_crossing_step :
                    polished.converged ? :continuation_branch_guard :
                    polished.stop_reason
                return (completed=false, polished, reason,
                    failed_spin=next_spin)
            end
        end
    end

    return (completed=true, polished=nothing, reason=:completed,
        failed_spin=nothing)
end

# When the 128-bit position certificate of a final root fails on its
# derivative, solve again at the same spin in 192-bit arithmetic (explicit
# guess, no continuation). The 192-bit solve cannot re-enter this recovery.
function _near_axis_position_precision_recovery(mode::QNMMode, a, omega;
        options, context, failed_primary, failed_guard)
    lower_precision_bits = _precision_bits(BigFloat)
    lower_precision_bits == 128 || throw(ArgumentError(
        "Near-axis position recovery requires a failed 128-bit finalization."))
    return setprecision(BigFloat, 192) do
        B = BigFloat
        precise_a = B(a)
        precise_seed = _complex_of_type(B, omega)
        result = _qnm_frequency_labeled(mode, precise_a;
            options..., guess=precise_seed, automatic_precision_fallback=false,
            fallback_precision_bits=192)
        result.precision_bits == 192 && precision(real(result.omega)) == 192 &&
            precision(imag(result.omega)) == 192 || throw(ErrorException(
                "Position recovery did not retain its actual 192-bit precision."))

        # The explicit-guess solve has no continuation history; the branch
        # guard uses the parent's path when it has at least two points.
        available = length(context.coordinates) >= 2
        branch_guard = available ? merge(_continuation_branch_guard(
            B.(context.coordinates), _complex_of_type.(Ref(B), context.roots),
            B(last(context.coordinates)), _complex_of_type(B, last(context.roots)),
            result.omega, B(context.continuation_root_tolerance)),
            (available=true, source=:pre_recovery_GSN_path)) :
            (available=false, accepted=missing,
                reason=:explicit_guess_has_no_continuation_history)
        child_position = get(result.provenance, :near_axis_position_check, nothing)
        child_position_verified = child_position !== nothing &&
            get(child_position, :verified, false) === true &&
            get(child_position, :precision_bits, 0) == 192
        status = result.status
        stop_reason = result.stop_reason
        if status in (:accepted, :estimated) && !child_position_verified
            status = :failed
            stop_reason = :near_axis_position_recovery_unverified
        end
        terminal_label_uncertain = false
        if status in (:accepted, :estimated) && available && !branch_guard.accepted
            status = :failed
            stop_reason = :near_axis_position_branch_guard
        end
        # Apply the parent's terminal-extrapolation label gate; a failed
        # 192-bit solve is never upgraded.
        if status == :accepted && !isempty(context.terminal_rows)
            terminal_row = last(context.terminal_rows)
            kappa = sqrt(max(zero(B), one(B) - precise_a^2)) /
                (2 * (one(B) + sqrt(max(zero(B), one(B) - precise_a^2))))
            jump = abs(terminal_row.candidate - context.roots[end-1])
            if jump > kappa / 4
                status = :estimated
                stop_reason = :terminal_extrapolation_label_uncertain
                terminal_label_uncertain = true
            end
        end
        p = result.provenance
        recovery = (
            policy=:near_axis_finalization_precision_recovery,
            trigger=failed_guard.stop_reason, source_precision_bits=lower_precision_bits,
            working_precision_bits=192, source_spin=a, exact_promoted_spin=precise_a,
            source_frequency=omega, promoted_seed=precise_seed,
            recovered_frequency=result.omega, recovery_options=options,
            continuation_recomputed=false, child_position_verified=child_position_verified,
            previous_attempt=(precision_bits=lower_precision_bits,
                final_primary=failed_primary, guard=failed_guard,
                provenance=context.provenance),
            child_final_primary_rows=get(p, :final_primary_rows, NamedTuple[]),
            child_guard_rows=get(p, :guard_rows, NamedTuple[]),
            branch_guard=branch_guard,
            terminal_label_uncertain=terminal_label_uncertain,
        )
        # Report the 192-bit solve's metrics and histories, with the parent's
        # seed and continuation path.
        merged = merge(context.provenance, p, (
            input_type=context.provenance.input_type,
            guess_source=context.provenance.guess_source,
            seed_inversion_index=context.provenance.seed_inversion_index,
            seed_spin=context.seed_spin,
            seed_predictor=context.provenance.seed_predictor,
            seed_duplicate_distance=context.provenance.seed_duplicate_distance,
            spin_path=context.provenance.spin_path,
            coordinate_path=context.provenance.coordinate_path,
            continuation_coordinate=context.provenance.continuation_coordinate,
            continuation_rows=context.provenance.continuation_rows,
            rejected_continuation_rows=context.provenance.rejected_continuation_rows,
            terminal_extrapolation_rows=context.terminal_rows,
            terminal_extrapolation=context.provenance.terminal_extrapolation,
            continuation_axis_touches=context.provenance.continuation_axis_touches,
            nia_restart_rows=context.provenance.nia_restart_rows,
            final_primary_rows=get(p, :final_primary_rows, NamedTuple[]),
            guard_rows=get(p, :guard_rows, NamedTuple[]),
            final_position_branch_guard=branch_guard,
            near_axis_position_recovery=recovery,
            conditioned_root=terminal_label_uncertain || get(p, :conditioned_root, false),
        ))
        return _with_provenance(result, merged; status, stop_reason)
    end
end

# Segment-0 (or seeded segment-1) solver behind `qnm_frequency`: seed, spin
# continuation, final polish with an independent higher-order guard, and gates.
function _qnm_frequency_labeled(mode::QNMMode, a; guess=nothing,
        convention=:overtone,
        inversion_index=:auto,
        continuation_coordinate=:auto,
        angular_order::Int=32,
        sheet_id::Symbol=:qnm_straight_from_spherical,
        cf_tolerance=nothing,
        cf_minimum_iterations::Int=300,
        cf_maximum_iterations=nothing,
        root_tolerance=nothing,
        stability_tolerance=nothing,
        maximum_root_iterations::Int=30,
        initial_spin_step=0.005,
        maximum_spin_step=0.02,
        minimum_spin_step=1.0e-6,
        trust_radius=0.2,
        automatic_precision_fallback::Bool=true,
        fallback_precision_bits::Int=128,
        segment::Int=0,
        segment_seed=nothing,
        spin_seed=nothing)
    _require_qnm_angular_backend()
    selected_convention = _qnm_convention(convention)
    T = _working_type(a, guess)
    aT = T(a)
    abs(aT) > one(T) && throw(DomainError(a,
        "A Kerr spin must satisfy abs(a) <= 1."))
    if aT < zero(T)
        source_mode = QNMMode(
            mode.s, mode.l, -mode.m, mode.n, mode.branch)
        source = _qnm_frequency_labeled(
            source_mode, -aT;
            guess,
            convention=selected_convention,
            inversion_index,
            continuation_coordinate,
            angular_order,
            sheet_id=Symbol(sheet_id, :_positive_spin_source),
            cf_tolerance,
            cf_minimum_iterations,
            cf_maximum_iterations,
            root_tolerance,
            stability_tolerance,
            maximum_root_iterations,
            initial_spin_step,
            maximum_spin_step,
            minimum_spin_step,
            trust_radius,
            automatic_precision_fallback,
            fallback_precision_bits,
            segment,
            spin_seed,
        )
        mapped_order = hasproperty(source.provenance, :angular_order) ?
            source.provenance.angular_order : angular_order
        return _spin_reflected_result(
            mode, aT, source; angular_order=mapped_order, sheet_id)
    end
    label_contract = _schwarzschild_label_contract(
        mode, selected_convention, aT; segment)
    near_as_role = label_contract.near_as_role
    automatic_inversion = inversion_index === :auto
    default_inversion = selected_convention == :leaver ||
        label_contract.algebraically_special_pair ? mode.n :
        label_contract.schwarzschild_sequence_index
    selected_inversion = automatic_inversion ? default_inversion :
        Int(inversion_index)
    selected_inversion >= 0 || throw(ArgumentError(
        "inversion_index must be nonnegative."))
    overtone_index = selected_convention == :leaver ? missing : mode.n
    selected_convention == :leaver && guess === nothing &&
        throw(ArgumentError(
            "convention=:leaver requires an explicit guess because a " *
            "continued-fraction inversion index does not uniquely label a root."))
    # At a = 0 both co-rotating branches 8_0 and 8_1 end at the algebraically special frequency -2i,
    # so convention=:complete_spectrum n = 8 and n = 9 return the same endpoint as convention=:overtone n = 8.
    if label_contract.algebraically_special_pair && iszero(aT)
        return _algebraically_special_endpoint_result(
            mode, aT; angular_order, sheet_id,
            requested_inversion_index=inversion_index,
            convention=selected_convention)
    end
    selected_coordinate = _resolve_continuation_coordinate(
        continuation_coordinate, aT, mode)
    tolerance_cf = cf_tolerance === nothing ?
        T(mode.n >= 8 ? 1.0e-12 : 1.0e-13) : T(cf_tolerance)
    maximum_cf_iterations = cf_maximum_iterations === nothing ?
        _default_cf_maximum_iterations(mode, aT, guess) :
        Int(cf_maximum_iterations)
    maximum_cf_iterations >= cf_minimum_iterations || throw(ArgumentError(
        "cf_maximum_iterations must not be below cf_minimum_iterations."))
    tolerance_root = root_tolerance === nothing ? T(1.0e-12) :
        T(root_tolerance)
    tolerance_stability = stability_tolerance === nothing ? T(1.0e-10) :
        T(stability_tolerance)
    continuation_root_tolerance = T === Float64 ?
        max(tolerance_root, T(5.0e-10)) : tolerance_root
    radius = T(trust_radius)
    initial_step = T(initial_spin_step)
    maximum_step = T(maximum_spin_step)
    minimum_step = T(minimum_spin_step)
    initial_step > zero(T) || throw(ArgumentError(
        "initial_spin_step must be positive."))
    maximum_step > zero(T) || throw(ArgumentError(
        "maximum_spin_step must be positive."))
    minimum_step > zero(T) || throw(ArgumentError(
        "minimum_spin_step must be positive."))
    minimum_step <= maximum_step || throw(ArgumentError(
        "minimum_spin_step must not exceed maximum_spin_step."))

    if abs(aT) == one(T) && segment != 0
        return _failed_result(mode, aT, complex(T(NaN), T(NaN)),
            :segment_unsupported_at_extremal, (
                input_type=string(typeof(a)),
                label_convention=label_contract.label_convention,
                branch_label=label_contract.branch_label,
                sequence_segment=segment,
            );
            convention=selected_convention, overtone_index,
            inversion_index=selected_inversion)
    end
    if abs(aT) == one(T)
        return _extremal_damped_limit_result(
            mode, aT;
            guess,
            convention=selected_convention,
            inversion_index,
            angular_order,
            sheet_id,
            cf_tolerance=tolerance_cf,
            cf_minimum_iterations,
            cf_maximum_iterations=maximum_cf_iterations,
            root_tolerance=tolerance_root,
            stability_tolerance=tolerance_stability,
            maximum_root_iterations,
            trust_radius=radius,
            automatic_precision_fallback,
            fallback_precision_bits,
        )
    end

    label_provenance = (
        label_convention=label_contract.label_convention,
        branch_label=label_contract.branch_label,
        endpoint_kind=label_contract.endpoint_kind,
        near_as_role,
        schwarzschild_sequence_index=
            label_contract.schwarzschild_sequence_index,
        schwarzschild_multiplicity=
            label_contract.schwarzschild_multiplicity,
        sequence_segment=segment,
    )
    counterrotating_seed = guess === nothing && segment_seed === nothing ?
        _counterrotating_near_extremal_seed(mode, aT) : nothing
    if counterrotating_seed !== nothing
        return _counterrotating_conditioned_result(
            mode, aT, counterrotating_seed;
            convention=selected_convention,
            overtone_index,
            inversion_index=selected_inversion,
            angular_order,
            sheet_id,
            cf_tolerance=tolerance_cf,
            cf_minimum_iterations,
            cf_maximum_iterations=maximum_cf_iterations,
            root_tolerance=tolerance_root,
            stability_tolerance=tolerance_stability,
        )
    end
    failed_result(omega, reason, provenance) = _failed_result(
        mode, aT, omega, reason, merge(label_provenance, provenance);
        convention=selected_convention,
        overtone_index,
        inversion_index=selected_inversion)

    seed_mode = guess === nothing && selected_convention != :leaver &&
        !ismissing(label_contract.schwarzschild_sequence_index) &&
        label_contract.schwarzschild_sequence_index != mode.n ?
        QNMMode(mode.s, mode.l, mode.m,
            label_contract.schwarzschild_sequence_index, mode.branch) : mode
    automatic_coordinate = Symbol(lowercase(String(
        continuation_coordinate))) == :auto
    seed_result = segment_seed !== nothing ? segment_seed :
        automatic_coordinate && guess === nothing &&
            selected_convention != :leaver &&
            mode.n <= 4 && _near_extremal_zdm_seed_scope(mode, aT) ?
        _near_extremal_zdm_seed(mode, aT, selected_inversion) :
        guess === nothing &&
            label_contract.algebraically_special_pair ?
        _near_as_anchor_seed(mode, aT, near_as_role) :
        guess === nothing && spin_seed !== nothing && segment == 0 &&
            automatic_coordinate && selected_convention == :overtone &&
            automatic_inversion && T === Float64 ? (
            seed=complex(T(real(spin_seed.seed)), T(imag(spin_seed.seed))),
            seed_spin=aT,
            source=spin_seed.source,
            seed_inversion_index=spin_seed.inversion,
            seed_predictor=spin_seed.predictor,
            seed_duplicate_distance=T(Inf),
        ) : guess === nothing ?
        _schwarzschild_seed(
            seed_mode, T;
            angular_order,
            sheet_id,
            cf_tolerance=tolerance_cf,
            cf_minimum_iterations,
            cf_maximum_iterations=maximum_cf_iterations,
            root_tolerance=continuation_root_tolerance,
            maximum_root_iterations,
            trust_radius=radius,
        ) : (
        seed=_complex_of_type(T, guess),
        seed_spin=aT,
        source=:explicit_guess,
        seed_inversion_index=selected_inversion,
        seed_predictor=:explicit,
        seed_duplicate_distance=T(Inf),
    )
    if automatic_inversion && selected_convention != :leaver &&
            guess === nothing
        selected_inversion = seed_result.seed_inversion_index
    end
    inversion_selection_policy = automatic_inversion ?
        (selected_convention != :leaver && guess === nothing ?
            (label_contract.algebraically_special_pair ?
                :finite_spin_complete_spectrum_anchor :
                :schwarzschild_seed_stability_selection) :
            :input_n_leaver_inversion) : :explicit
    initial = seed_result.seed
    guess_source = seed_result.source
    _branch_consistent(mode, initial) || throw(ArgumentError(
        "The supplied guess is inconsistent with mode.branch."))
    current = initial
    continuation_rows = NamedTuple[]
    rejected_continuation_rows = NamedTuple[]
    terminal_extrapolation_rows = NamedTuple[]
    accepted_spins = T[]
    accepted_coordinates = T[]
    accepted_roots = Complex{T}[]
    accepted_angular_A = Complex{T}[]
    last_polish = nothing

    polish_options = (;
        angular_order, sheet_id, inversion_index=selected_inversion,
        cf_tolerance=tolerance_cf, cf_minimum_iterations,
        cf_maximum_iterations=maximum_cf_iterations,
        root_tolerance=continuation_root_tolerance,
        maximum_root_iterations, trust_radius=radius, near_as_role)
    fallback_options = (;
        angular_order, sheet_id, cf_minimum_iterations,
        cf_maximum_iterations=maximum_cf_iterations,
        maximum_root_iterations, trust_radius=radius,
        fallback_precision_bits, cf_tolerance=tolerance_cf,
        root_tolerance=tolerance_root, stability_tolerance=tolerance_stability,
        continuation_coordinate=selected_coordinate,
        convention=selected_convention, overtone_index,
        inversion_index=selected_inversion, segment)

    try
        if guess === nothing
            spin = seed_result.seed_spin
            coordinate = selected_coordinate == :nu ? asin(spin) : spin
            polished = _polish_frequency_root(mode, spin, current;
                polish_options...)
            append!(continuation_rows, [merge(row, (a=spin,
                continuation_attempt=1, predictor=current))
                for row in polished.history])
            last_polish = polished
            if !polished.converged
                return failed_result(polished.omega,
                    polished.stop_reason, (
                        input_type=string(typeof(a)),
                        guess_source=guess_source,
                        spin_path=accepted_spins,
                        continuation_rows=continuation_rows,
                        rejected_continuation_rows=rejected_continuation_rows,
                    ))
            end
            current = polished.omega
            push!(accepted_spins, spin)
            push!(accepted_coordinates, coordinate)
            push!(accepted_roots, current)
            push!(accepted_angular_A,
                polished.evaluation.angular.angular_A)

            target_coordinate = selected_coordinate == :nu ? asin(aT) : aT
            sequence_anchor = _spin_sequence_anchor(mode, T)
            # A two-end anchored seed is polished at the target itself.
            sequence_target = selected_coordinate == :nu ||
                guess_source == :packaged_two_end_anchor ?
                target_coordinate : abs(aT) > sequence_anchor ?
                copysign(sequence_anchor, aT) : aT
            continuation_state = (; accepted_spins, accepted_coordinates,
                accepted_roots, accepted_angular_A, continuation_rows,
                rejected_continuation_rows,
                last_polish=Ref{Any}(last_polish),
                nia_restart_rows=NamedTuple[])
            outcome = _continue_with_negative_imaginary_axis_restarts!(continuation_state, mode,
                sequence_target; selected_coordinate, initial_step,
                minimum_step, maximum_step, continuation_root_tolerance,
                polish_options, sequence_segment=segment)
            last_polish = continuation_state.last_polish[]
            current = accepted_roots[end]
            if !outcome.completed
                classified = _classify_continuation_failure(mode,
                    outcome.reason, outcome.polished, outcome.failed_spin,
                    accepted_spins, accepted_roots, accepted_angular_A,
                    rejected_continuation_rows; sequence_segment=segment,
                    after_birth=segment >= 1)
                failure_provenance = merge((
                        input_type=string(typeof(a)),
                        guess_source=guess_source,
                        spin_path=accepted_spins,
                        coordinate_path=accepted_coordinates,
                        continuation_coordinate=selected_coordinate,
                        continuation_rows=continuation_rows,
                        rejected_continuation_rows=
                            rejected_continuation_rows,
                        failed_spin=outcome.failed_spin,
                        continuation_axis_touches=count(r -> r.accepted,
                            continuation_state.nia_restart_rows),
                        nia_restart_rows=continuation_state.nia_restart_rows,
                        continuation_settings=(; polish_options,
                            initial_step, minimum_step, maximum_step,
                            continuation_root_tolerance),
                    ), classified.provenance)
                return failed_result(outcome.polished.omega,
                    classified.reason, failure_provenance)
            end

            if selected_coordinate == :spin && sequence_target != aT
                predictor = _terminal_extrapolation_predictor(
                    accepted_coordinates, accepted_roots, aT)
                polished = _polish_frequency_root(mode, aT, predictor;
                    polish_options...)
                append!(continuation_rows, [merge(row, (
                    a=aT, continuation_attempt=1,
                    predictor=predictor,
                    terminal_extrapolation=true))
                    for row in polished.history])
                last_polish = polished
                if !polished.converged
                    terminal_provenance = (
                        input_type=string(typeof(a)),
                        guess_source=:terminal_cubic_extrapolation,
                        spin_path=accepted_spins,
                        coordinate_path=accepted_coordinates,
                        continuation_coordinate=selected_coordinate,
                        continuation_rows=continuation_rows,
                        rejected_continuation_rows=
                            rejected_continuation_rows,
                        terminal_predictor=predictor,
                        terminal_stop_reason=polished.stop_reason,
                        continuation_axis_touches=count(r -> r.accepted,
                            continuation_state.nia_restart_rows),
                        nia_restart_rows=continuation_state.nia_restart_rows,
                    )
                    if T === Float64 && automatic_precision_fallback
                        return _automatic_precision_fallback(
                            mode, aT, polished.omega;
                            fallback_options...,
                            float64_provenance=terminal_provenance)
                    end
                    return failed_result(polished.omega,
                        polished.stop_reason, terminal_provenance)
                end
                branch_guard = _continuation_branch_guard(
                    accepted_coordinates, accepted_roots, aT, predictor,
                    polished.omega, continuation_root_tolerance)
                push!(terminal_extrapolation_rows, (
                    a=aT,
                    attempted_step=abs(aT - sequence_target),
                    attempt=1,
                    predictor=predictor,
                    candidate=polished.omega,
                    converged=polished.converged,
                    stop_reason=:terminal_extrapolation_diagnostic,
                    branch_discrepancy=branch_guard.discrepancy,
                    branch_allowance=branch_guard.allowance,
                ))
                current = polished.omega
                push!(accepted_spins, aT)
                push!(accepted_coordinates, aT)
                push!(accepted_roots, current)
                push!(accepted_angular_A,
                    polished.evaluation.angular.angular_A)
            end
        else
            polished = _polish_frequency_root(mode, aT, current;
                polish_options...)
            append!(continuation_rows, [merge(row, (a=aT,
                continuation_attempt=1, predictor=current))
                for row in polished.history])
            last_polish = polished
            if !polished.converged
                classified = _classify_continuation_failure(mode,
                    polished.stop_reason, polished, aT, T[], Complex{T}[],
                    Complex{T}[], NamedTuple[(predictor=current,)];
                    sequence_segment=segment, after_birth=segment >= 1)
                return failed_result(polished.omega,
                    classified.reason, merge((
                        input_type=string(typeof(a)),
                        guess_source=:explicit_guess,
                        spin_path=T[aT],
                        continuation_rows=continuation_rows,
                        rejected_continuation_rows=rejected_continuation_rows,
                    ), classified.provenance))
            end
            current = polished.omega
            push!(accepted_spins, aT)
            push!(accepted_coordinates,
                selected_coordinate == :nu ? asin(aT) : aT)
            push!(accepted_roots, current)
            push!(accepted_angular_A,
                polished.evaluation.angular.angular_A)
        end

        mst_primary = nothing
        final_primary = nothing
        # Use the Float64 risk geometry for BigFloat finalizations (automatic
        # precision fallback): BigFloat epsilon would shrink this NIA region.
        final_position_tolerance = segment >= 1 && mode.m == 0 &&
            T === BigFloat && _precision_bits(T) >= 128 &&
            zero(T) < _branch_sign(mode) * real(current) <
                2 * T(_root_derivative_step(Float64, ComplexF64(current))) *
                    T(_near_extremal_derivative_scale(mode, Float64(aT))) ?
            tolerance_stability : nothing
        if _near_as_mano_suzuki_takasugi_finalization_required(
                mode, aT, T, near_as_role)
            mst_primary = _polish_mano_suzuki_takasugi_incidence_root(
                mode, Float64(aT), ComplexF64(current);
                angular_order,
                sheet_id=Symbol(sheet_id, :_mst_final),
                root_tolerance=Float64(max(tolerance_root, T(2.0e-11))),
                maximum_root_iterations,
                trust_radius=Float64(radius),
            )
            guard = mst_primary.converged ? _polish_mano_suzuki_takasugi_incidence_root(
                mode, Float64(aT), ComplexF64(mst_primary.omega);
                angular_order=angular_order + 8,
                sheet_id=Symbol(sheet_id, :_mst_guard),
                root_tolerance=Float64(max(tolerance_root, T(2.0e-11))),
                maximum_root_iterations,
                trust_radius=Float64(radius),
                derivative_step_scale=0.5,
            ) : mst_primary
            current = mst_primary.omega
        else
            # Continuation deliberately uses a looser Float64 tolerance so that
            # an ill-conditioned high-overtone root can be tracked cheaply.
            # Polish at the final tolerance before measuring the independent
            # angular-order drift; otherwise ordinary Newton correction is
            # misclassified as a branch jump.
            final_primary = _polish_frequency_root(mode, aT, current;
                polish_options...,
                sheet_id=Symbol(sheet_id, :_final_primary),
                root_tolerance=tolerance_root,
                position_tolerance=final_position_tolerance)
            current = final_primary.omega
            guard = final_primary.converged ? _polish_frequency_root(
                mode, aT, current;
                polish_options...,
                angular_order=angular_order + 8,
                sheet_id=Symbol(sheet_id, :_guard),
                cf_minimum_iterations=cf_minimum_iterations + 100,
                root_tolerance=tolerance_root,
                derivative_step_scale=T(0.5),
                position_tolerance=final_position_tolerance) : final_primary
        end
        final_position_branch_guard = nothing
        if final_position_tolerance !== nothing && guard.converged
            # Apply the continuation's branch allowance to the final
            # correction from the GSN path, without limiting that correction
            # itself by the much smaller convergence/stability tolerance.
            available = length(accepted_coordinates) >= 2
            final_position_branch_guard = available ? merge(
                _continuation_branch_guard(accepted_coordinates, accepted_roots,
                    last(accepted_coordinates), last(accepted_roots), guard.omega,
                    continuation_root_tolerance), (available=true,)) :
                (available=false, accepted=missing,
                    reason=:explicit_guess_has_no_continuation_history)
            if available && !final_position_branch_guard.accepted
                guard = merge(guard, (converged=false, conditioned=false,
                    stop_reason=:near_axis_position_branch_guard))
            end
        end
        float64_provenance = (
            input_type=string(typeof(a)),
            working_type=string(T),
            guess_source=guess_source,
            inversion_selection_policy,
            seed_inversion_index=seed_result.seed_inversion_index,
            seed_predictor=seed_result.seed_predictor,
            seed_duplicate_distance=seed_result.seed_duplicate_distance,
            spin_path=accepted_spins,
            coordinate_path=accepted_coordinates,
            continuation_coordinate=selected_coordinate,
            spin_sequence_anchor=_spin_sequence_anchor(mode, T),
            terminal_extrapolation=selected_coordinate == :spin &&
                guess === nothing &&
                abs(aT) > _spin_sequence_anchor(mode, T),
            initial_spin_step=initial_step,
            maximum_spin_step=maximum_step,
            minimum_spin_step=minimum_step,
            continuation_rows=continuation_rows,
            rejected_continuation_rows=rejected_continuation_rows,
            terminal_extrapolation_rows=terminal_extrapolation_rows,
            continuation_axis_touches=guess === nothing ?
                count(r -> r.accepted, continuation_state.nia_restart_rows) : 0,
            nia_restart_rows=guess === nothing ?
                continuation_state.nia_restart_rows : NamedTuple[],
            near_as_mst_finalization=mst_primary !== nothing,
            mst_primary_rows=mst_primary === nothing ?
                NamedTuple[] : mst_primary.history,
            final_primary_rows=final_primary === nothing ?
                NamedTuple[] : final_primary.history,
            guard_rows=guard.history,
        )
        guard_root_equation = _root_equation(guard.evaluation)
        if !guard.converged
            if final_position_tolerance !== nothing && T === BigFloat &&
                    _precision_bits(T) == 128 &&
                    guard.stop_reason == :near_axis_position_derivative_unverified
                return _near_axis_position_precision_recovery(mode, aT, guard.omega;
                    options=(angular_order=angular_order, sheet_id=sheet_id,
                        convention=selected_convention, inversion_index=selected_inversion,
                        segment=segment, cf_tolerance=tolerance_cf,
                        cf_minimum_iterations=cf_minimum_iterations,
                        cf_maximum_iterations=maximum_cf_iterations,
                        root_tolerance=tolerance_root, stability_tolerance=tolerance_stability,
                        maximum_root_iterations=maximum_root_iterations, trust_radius=radius,
                        continuation_coordinate=selected_coordinate,
                        initial_spin_step=initial_step, minimum_spin_step=minimum_step,
                        maximum_spin_step=maximum_step),
                    context=(coordinates=accepted_coordinates, roots=accepted_roots,
                        continuation_root_tolerance=continuation_root_tolerance,
                        terminal_rows=terminal_extrapolation_rows,
                        seed_spin=seed_result.seed_spin, provenance=float64_provenance),
                    failed_primary=final_primary, failed_guard=guard)
            end
            if T === Float64 && automatic_precision_fallback &&
                    guard_root_equation == :leaver_continued_fraction
                return _automatic_precision_fallback(
                    mode, aT, guard.omega;
                    fallback_options...,
                    float64_provenance=merge(float64_provenance, (
                        float64_stop_reason=guard.stop_reason,)))
            end
            return failed_result(
                guard.omega, guard.stop_reason, (
                input_type=string(typeof(a)), spin_path=accepted_spins,
                coordinate_path=accepted_coordinates,
                continuation_coordinate=selected_coordinate,
                continuation_rows=continuation_rows,
                rejected_continuation_rows=rejected_continuation_rows,
                guard_rows=guard.history,
            ))
        end
        root_drift = T(abs(guard.omega - current) /
            max(one(T), abs(guard.omega), abs(current)))
        final = guard.evaluation
        root_equation = _root_equation(final)
        final_root_residual = T(_root_residual(final))
        effective_root_tolerance = root_equation == :gsn_mst_incidence ?
            max(tolerance_root, T(2.0e-11)) : tolerance_root
        cf_tail_gate_passed = final.error <= tolerance_cf
        near_as_cf_tail_ceiling = T(1.0e-8)
        near_as_root_position_certificate =
            root_equation == :leaver_continued_fraction &&
            _near_as_complete_scope(mode, near_as_role) &&
            near_as_role == :unconventional &&
            zero(T) < abs(aT) < T(_KERR_L2_M2_NEAR_AS_ANCHOR_SPIN) &&
            root_drift <= tolerance_stability &&
            final_root_residual <= effective_root_tolerance &&
            final.error <= near_as_cf_tail_ceiling
        accepted = root_drift <= tolerance_stability &&
            final_root_residual <= effective_root_tolerance &&
            (cf_tail_gate_passed || near_as_root_position_certificate)
        conditioned_root_limit = _conditioned_root_residual_limit(
            mode, aT, tolerance_root)
        conditioned_stability_limit = _conditioned_root_stability_limit(
            mode, aT, tolerance_stability)
        conditioned_root = root_equation == :leaver_continued_fraction &&
            !accepted && root_drift <= conditioned_stability_limit &&
            final_root_residual <= conditioned_root_limit &&
            cf_tail_gate_passed &&
            ((hasproperty(guard, :conditioned) && guard.conditioned) ||
             (final_primary !== nothing &&
              hasproperty(final_primary, :conditioned) &&
              final_primary.conditioned))
        terminal_label_uncertain = false
        if !accepted && !conditioned_root && T === Float64 &&
                automatic_precision_fallback &&
                root_equation == :leaver_continued_fraction
            return _automatic_precision_fallback(
                mode, aT, guard.omega;
                fallback_options...,
                float64_provenance=merge(float64_provenance, (
                    float64_root_drift=root_drift,
                    float64_root_residual=final_root_residual,
                    float64_cf_error=final.error,
                )))
        end
        # A cubic jump from the sequence anchor cannot resolve the overtone
        # label once the predicted move exceeds a fraction of the local level
        # spacing, which near extremality is set by the surface gravity.
        if accepted && !isempty(terminal_extrapolation_rows)
            terminal_row = last(terminal_extrapolation_rows)
            kappa = sqrt(max(zero(T), one(T) - aT^2)) /
                (2 * (one(T) + sqrt(max(zero(T), one(T) - aT^2))))
            jump = abs(terminal_row.candidate - accepted_roots[end-1])
            if jump > kappa / 4
                accepted = false
                conditioned_root = true
                terminal_label_uncertain = true
            end
        end
        # Next to the NIA cut the Float64 continued fraction is poorly
        # conditioned. A re-emerged segment's Float64 root within 2h of the
        # axis is kept only when a BigFloat re-polish reproduces it to the
        # stability gate.
        near_axis_check = nothing
        if segment >= 1 && T === Float64 && (accepted || conditioned_root) &&
                root_equation == :leaver_continued_fraction &&
                _branch_sign(mode) * real(guard.omega) <
                2 * T(_root_derivative_step(T, guard.omega)) *
                    _near_extremal_derivative_scale(mode, aT)
            near_axis_check = _near_axis_precision_check(
                mode, aT, guard.omega, polish_options;
                precision_bits=fallback_precision_bits,
                root_tolerance=tolerance_root,
                stability_tolerance=tolerance_stability)
            if !near_axis_check.verified
                accepted = false
                conditioned_root = false
            end
        end
        stop_reason = near_axis_check !== nothing &&
            !near_axis_check.verified ? :near_axis_precision_unverified :
            terminal_label_uncertain ?
            :terminal_extrapolation_label_uncertain :
            accepted ? :accepted : conditioned_root ?
            :conditioned_near_extremal_cf_floor :
            (root_drift > tolerance_stability ? :root_stability_gate :
             final_root_residual > effective_root_tolerance ?
                (root_equation == :gsn_mst_incidence ?
                    :mst_incidence_residual_gate : :raw_cf_residual_gate) :
             root_equation == :gsn_mst_incidence ?
                :mst_representation_error_gate : :lentz_error_gate)
        output_inversion = root_equation == :gsn_mst_incidence ?
            missing : selected_inversion
        return LeaverResult(
            mode, selected_convention, overtone_index,
            output_inversion, aT, guard.omega,
            final.angular.angular_A, final.angular.lambda,
            final.angular.mixing, final.value, T(final.error),
            final.iterations, final_root_residual,
            T(final.angular.residual), _precision_bits(T),
            accepted ? :accepted : conditioned_root ? :estimated : :failed,
            stop_reason,
            (
                input_type=string(typeof(a)),
                working_type=string(T),
                guess_source=guess_source,
                inversion_selection_policy,
                seed_inversion_index=seed_result.seed_inversion_index,
                seed_spin=seed_result.seed_spin,
                seed_predictor=seed_result.seed_predictor,
                seed_duplicate_distance=seed_result.seed_duplicate_distance,
                convention=selected_convention,
                overtone_index=overtone_index,
                inversion_index=selected_inversion,
                label_convention=label_contract.label_convention,
                branch_label=label_contract.branch_label,
                endpoint_kind=label_contract.endpoint_kind,
                near_as_role,
                schwarzschild_sequence_index=
                    label_contract.schwarzschild_sequence_index,
                schwarzschild_multiplicity=
                    label_contract.schwarzschild_multiplicity,
                sequence_segment=segment,
                near_axis_precision_check=near_axis_check,
                near_axis_position_check=hasproperty(guard, :position_check) ?
                    guard.position_check : nothing,
                final_position_branch_guard,
                continuation_axis_touches=guess === nothing ?
                    count(r -> r.accepted, continuation_state.nia_restart_rows) : 0,
                nia_restart_rows=guess === nothing ?
                    continuation_state.nia_restart_rows : NamedTuple[],
                family=_near_extremal_family(mode, aT, guard.omega),
                spin_path=accepted_spins,
                coordinate_path=accepted_coordinates,
                continuation_coordinate=selected_coordinate,
                maximum_spin_step=maximum_step,
                minimum_spin_step=minimum_step,
                angular_route=final.angular.metadata.route,
                angular_sheet=final.angular.sheet_id,
                angular_order=final.angular.truncation_order,
                cf_minimum_iterations=cf_minimum_iterations + 100,
                cf_maximum_iterations=maximum_cf_iterations,
                cf_tolerance=tolerance_cf,
                root_tolerance=tolerance_root,
                effective_root_tolerance,
                conditioned_root_limit,
                conditioned_stability_limit,
                conditioned_root,
                stability_tolerance=tolerance_stability,
                root_drift=root_drift,
                root_equation,
                acceptance_basis=near_as_root_position_certificate ?
                    :near_as_root_position_stability : root_equation,
                cf_tail_gate_passed,
                near_as_cf_tail_ceiling,
                radial_cf_applicable=
                    root_equation == :leaver_continued_fraction,
                near_as_mst_finalization=mst_primary !== nothing,
                mst_primary_rows=mst_primary === nothing ?
                    NamedTuple[] : mst_primary.history,
                final_primary_rows=final_primary === nothing ?
                    NamedTuple[] : final_primary.history,
                continuation_rows=continuation_rows,
                rejected_continuation_rows=rejected_continuation_rows,
                guard_rows=guard.history,
                upstream_qnm_commit=
                    "f3abd18e59828e7e7d75d07f20c7cbc87925edfa",
            ),
        )
    catch error
        return failed_result(current, :evaluation_error, (
            input_type=string(typeof(a)),
            spin_path=accepted_spins,
            continuation_rows=continuation_rows,
            rejected_continuation_rows=rejected_continuation_rows,
            error=sprint(showerror, error),
            last_polish=last_polish,
        ))
    end
end
