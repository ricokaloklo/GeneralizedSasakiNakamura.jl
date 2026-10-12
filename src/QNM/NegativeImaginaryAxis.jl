# Sequences that end on, and re-emerge from, the negative imaginary axis (NIA).
#
# On the NIA the Leaver continued fraction has no minimal solution; a mode
# there is a confluent Heun polynomial solution (Cook & Zalutskiy, PRD 94,
# 104074 (2016), Sec. III).  For s=-2, l=2, m=0 and 9 <= n <= 26 the sequence
# attached to a=0 (Cook's n_0) ends on the axis at finite spin, and after a gap
# a second segment n_1 leaves the axis at larger spin.  Segments are ordered by
# spin; values are never interpolated across the gap.
#
# Three design constants of the re-emergence search:
#   - edge line Re(omega) = 8h, with h the unclamped derivative step at the
#     endpoint: the stencil clamp keeps every Newton stencil on the physical
#     side, and 8h puts the scan line outside the region where a Newton start
#     would be clamped to a few-h stencil;
#   - window gamma_end +- 1/4 in gamma = -Im(omega): the polynomial frequencies
#     -i N/4 are spaced by 1/4, so the window covers one spacing on each side;
#   - candidate cap Re(omega) <= 1/4: a re-emerged root is born on the axis,
#     so at the first spin rung after the birth it is still close to it.

const _NIA_EDGE_STEP_FACTOR = 8
const _NIA_SEARCH_HALF_WINDOW = 1 // 4
const _NIA_CANDIDATE_REAL_CAP = 1 // 4
const _NIA_EDGE_SCAN_POINTS = 33
const _NIA_FIRST_SEARCH_SPAN = 1 // 10

"""
    _nia_endpoint_detected(mode, last_root, last_predictor, last_stop_reason;
        derivative_step_scale=1)

Classify a stalled continuation as reaching the NIA: the last accepted root is
within two unclamped derivative steps of the axis, or the last predictor lies
on or beyond it, or the last polish started off the branch.
"""
function _nia_endpoint_detected(mode::QNMMode, last_root, last_predictor,
        last_stop_reason::Symbol; derivative_step_scale=1)
    _nia_crossing_exempt(mode) && return false
    last_stop_reason == :initial_guess_off_branch && return true
    sign = _branch_sign(mode)
    if last_root !== nothing && isfinite(last_root)
        T = typeof(float(real(last_root)))
        h = T(_root_derivative_step(T, last_root)) * derivative_step_scale
        sign * real(last_root) <= 2h && return true
    end
    return last_predictor !== nothing && isfinite(last_predictor) &&
        sign * real(last_predictor) <= 0
end

# Extrapolate the last (at most three) accepted roots to Re(omega)=0. Near a
# collision with the mirror root Re(omega)^2 is regular in a, so the zero of
# its interpolant is used; it must lie within four bracket widths.
function _nia_endpoint_extrapolation(mode::QNMMode, spins, roots, failed_spin)
    isempty(spins) && return (spin=missing, frequency=missing,
        method=:unavailable)
    T = eltype(spins)
    count = min(3, length(spins))
    xs = spins[(end - count + 1):end]
    ys = roots[(end - count + 1):end]
    sign = _branch_sign(mode)
    rho2 = [(sign * real(y))^2 for y in ys]
    last_spin = xs[end]
    width = abs(T(failed_spin) - last_spin)
    direction = _sign_or_one(T(failed_spin) - last_spin)
    candidates = count == 1 ? T[] : _interpolant_zeros(xs, rho2)
    filter!(x -> isfinite(x) && direction * (x - last_spin) >= 0 &&
        abs(x - last_spin) <= 4width, candidates)
    spin, method = isempty(candidates) ?
        (last_spin + direction * width / 2, :bracket_midpoint) :
        (candidates[argmin(abs.(candidates .- last_spin))],
            :real_part_squared_interpolant)
    gamma = _polynomial_predictor(xs, imag.(ys), spin)
    return (spin=spin, frequency=complex(zero(T), gamma), method=method)
end

_sign_or_one(x) = x < zero(x) ? -one(x) : one(x)

# Real zeros of the polynomial through (xs[i], ys[i]), length 2 or 3.
function _interpolant_zeros(xs, ys)
    T = eltype(xs)
    found = T[]
    if length(xs) == 2
        slope = (ys[2] - ys[1]) / (xs[2] - xs[1])
        iszero(slope) || push!(found, xs[2] - ys[2] / slope)
    elseif length(xs) == 3
        d1 = (ys[2] - ys[1]) / (xs[2] - xs[1])
        d2 = (ys[3] - ys[2]) / (xs[3] - xs[2])
        c2 = (d2 - d1) / (xs[3] - xs[1])
        b = d2 + c2 * (xs[3] - xs[2])
        if iszero(c2)
            iszero(b) || push!(found, xs[3] - ys[3] / b)
        else
            discriminant = b^2 - 4 * c2 * ys[3]
            if discriminant >= 0
                root = sqrt(discriminant)
                push!(found, xs[3] + (-b + root) / (2c2),
                    xs[3] + (-b - root) / (2c2))
            end
        end
    end
    return filter(isfinite, found)
end

"""
    _nia_step_crosses_axis(mode, coordinates, roots, next_coordinate)

True when the interpolant of Re(omega)^2 or of Re(omega) through the last (at
most three) accepted points has a zero between the last accepted coordinate and
`next_coordinate`: the path reaches the NIA inside the step, so a root found at
`next_coordinate` is not the continuation of this sequence. Re(omega)^2 is
regular where the root meets its mirror (Re(omega) ~ sqrt(a_end - a)); Re(omega)
itself is regular where the sequence ends on a continued-fraction pole -i N/4
(Re(omega) ~ a_end - a, e.g. l=2, m=0, n=13).
"""
function _nia_step_crosses_axis(mode::QNMMode, coordinates, roots,
        next_coordinate)
    _nia_crossing_exempt(mode) && return false
    count = min(3, length(coordinates))
    count >= 2 || return false
    xs = coordinates[(end - count + 1):end]
    rho = [real(root) for root in roots[(end - count + 1):end]]
    last = xs[end]
    direction = _sign_or_one(next_coordinate - last)
    inside(x) = zero(x) < direction * (x - last) <=
        direction * (next_coordinate - last)
    return any(inside, _interpolant_zeros(xs, rho .^ 2)) ||
        any(inside, _interpolant_zeros(xs, rho))
end

# A later segment can touch the NIA without ending. Before returning an
# endpoint failure, restart on the physical side at Re(omega)=+/-1e-3 and
# continue in the original direction. Segment 0 keeps its original path.
function _continue_with_nia_restarts!(state, mode::QNMMode, target::T;
        sequence_segment::Int=0, selected_coordinate::Symbol,
        initial_step::T, minimum_step::T, maximum_step::T,
        continuation_root_tolerance::T, polish_options) where {T<:AbstractFloat}
    settings = (; selected_coordinate, initial_step, minimum_step,
        maximum_step, continuation_root_tolerance, polish_options)
    outcome = _continue_in_spin!(state, mode, target; settings...)
    while !outcome.completed && sequence_segment > 0
        classified = _classify_continuation_failure(mode, outcome.reason,
            outcome.polished, outcome.failed_spin, state.accepted_spins,
            state.accepted_roots, state.accepted_angular_A,
            state.rejected_continuation_rows; sequence_segment)
        # The continuation already established an interpolated axis crossing.
        # At higher precision its accepted root and final predictor can both
        # remain farther from the axis than the derivative-step heuristic.
        # Preserve that explicit crossing evidence when dispatching a later
        # segment's guarded restart; no root-acceptance condition changes.
        (outcome.reason == :nia_crossing_step ||
            classified.reason == :negative_imaginary_axis_endpoint) || break
        direction = sign(target - state.accepted_coordinates[end])
        failed_coordinate = selected_coordinate == :nu ?
            asin(T(outcome.failed_spin)) : T(outcome.failed_spin)
        remaining = abs(target - state.accepted_coordinates[end])
        step = min(initial_step, remaining)
        recovered = false
        while step >= minimum_step
            coordinate = state.accepted_coordinates[end] + direction * step
            # Reconnecting on the near side is not a crossing. Do not count
            # it or repeatedly restart without passing the failed bracket.
            direction * (coordinate - failed_coordinate) > zero(T) || break
            spin = selected_coordinate == :nu ? sin(coordinate) : coordinate
            # Restart at the last accepted imaginary frequency. Extrapolating
            # a tiny-step polynomial across the touch can jump to another root.
            seed = complex(T(_branch_sign(mode)) * T(1.0e-3),
                imag(state.accepted_roots[end]))
            # The off-axis Newton seed is an artificial displacement, not
            # a prediction of this branch. Compare the polished root with
            # the last accepted physical root; reducing the spin step still
            # enforces the unchanged continuation guard across the touch.
            predictor = state.accepted_roots[end]
            polished = _polish_frequency_root(mode, spin, seed;
                polish_options...)
            guard = polished.converged ? _continuation_branch_guard(
                state.accepted_coordinates, state.accepted_roots,
                coordinate, predictor, polished.omega,
                continuation_root_tolerance) :
                (accepted=false, discrepancy=T(Inf), allowance=zero(T))
            accepted = polished.converged && guard.accepted &&
                _branch_consistent(mode, polished.omega)
            push!(state.nia_restart_rows, (a=spin, coordinate,
                seed, predictor, candidate=polished.omega, accepted,
                branch_predictor_kind=:last_accepted_root,
                stop_reason=polished.stop_reason,
                original_failed_spin=outcome.failed_spin,
                original_stop_reason=outcome.reason,
                branch_discrepancy=guard.discrepancy,
                branch_allowance=guard.allowance))
            append!(state.continuation_rows, [merge(row, (a=spin,
                continuation_coordinate=selected_coordinate,
                coordinate_value=coordinate, predictor=predictor,
                nia_restart_seed=seed, nia_restart=true))
                for row in polished.history])
            if accepted
                # Resume with a fresh local predictor history on the far
                # side, retaining the complete pre-touch path and failures.
                resumed = (accepted_spins=T[spin],
                    accepted_coordinates=T[coordinate],
                    accepted_roots=Complex{T}[polished.omega],
                    accepted_angular_A=Complex{T}[
                        polished.evaluation.angular.angular_A],
                    continuation_rows=NamedTuple[],
                    rejected_continuation_rows=NamedTuple[],
                    last_polish=Ref{Any}(polished))
                outcome = _continue_in_spin!(resumed, mode, target;
                    settings...)
                for key in (:accepted_spins, :accepted_coordinates,
                        :accepted_roots, :accepted_angular_A,
                        :continuation_rows, :rejected_continuation_rows)
                    append!(getproperty(state, key), getproperty(resumed, key))
                end
                state.last_polish[] = resumed.last_polish[]
                recovered = true
                break
            end
            step /= 2
        end
        recovered || break
    end
    return outcome
end

"""
    _classify_continuation_failure(mode, reason, polished, failed_spin,
        spins, roots, angular_A, rejected_rows; sequence_segment=0)

Attach the accepted root path to a stalled continuation and, when
[`_nia_endpoint_detected`](@ref) holds, reclassify the failure as
`:negative_imaginary_axis_endpoint` with the extrapolated endpoint.
"""
function _classify_continuation_failure(mode::QNMMode, reason::Symbol,
        polished, failed_spin, spins, roots, angular_A, rejected_rows;
        sequence_segment::Int=0, after_birth::Bool=false)
    root_path = copy(roots)
    angular_path = copy(angular_A)
    last_root = isempty(roots) ? polished.omega : roots[end]
    last_predictor = isempty(rejected_rows) ? nothing :
        rejected_rows[end].predictor
    T = typeof(real(last_root))
    scale = _near_extremal_derivative_scale(mode, T(failed_spin))
    detected = _nia_endpoint_detected(mode, last_root, last_predictor,
        polished.stop_reason; derivative_step_scale=scale)
    detected || return (reason=reason, provenance=(root_path=root_path,
        angular_path=angular_path))
    endpoint = _nia_endpoint_extrapolation(mode, spins, roots, failed_spin)
    # A re-emerged segment continues to a -> 1 while touching the axis ever
    # more often (Cook & Zalutskiy); a stall there after the tangent restarts
    # is not an end of the sequence.
    after_birth && return (reason=:negative_imaginary_axis_tangent_accumulation,
        provenance=(
        endpoint_kind=:negative_imaginary_axis_tangent_accumulation,
        stall_spin=endpoint.spin,
        stall_frequency=endpoint.frequency,
        stall_spin_bracket=(isempty(spins) ? missing : spins[end],
            failed_spin),
        root_path=root_path,
        angular_path=angular_path,
        continuation_stop_reason=reason,
        sequence_segment=sequence_segment,
        next_segment_hint=nothing,
    ))
    return (reason=:negative_imaginary_axis_endpoint, provenance=(
        endpoint_kind=:negative_imaginary_axis,
        endpoint_spin=endpoint.spin,
        endpoint_frequency=endpoint.frequency,
        endpoint_extrapolation=endpoint.method,
        endpoint_spin_bracket=(isempty(spins) ? missing : spins[end],
            failed_spin),
        root_path=root_path,
        angular_path=angular_path,
        continuation_stop_reason=reason,
        sequence_segment=sequence_segment,
        next_segment_hint=sequence_segment < 1 ? sequence_segment + 1 :
            nothing,
    ))
end

_segment_succeeded(result) = result.status in (:accepted, :estimated)

function _segment_failure(mode::QNMMode, a::T, omega, reason::Symbol,
        provenance::NamedTuple; convention::Symbol,
        inversion_index::Int, segment) where {T<:AbstractFloat}
    contract = _schwarzschild_label_contract(mode, convention, a;
        segment=segment === :auto ? 0 : segment)
    return _failed_result(mode, a, _complex_of_type(T, omega), reason,
        merge((
            label_convention=contract.label_convention,
            branch_label=contract.branch_label,
            endpoint_kind=contract.endpoint_kind,
            sequence_segment=segment,
        ), provenance);
        convention,
        overtone_index=convention == :leaver ? missing : mode.n,
        inversion_index)
end

"""
    _qnm_frequency_segmented(mode, a, segment; kwargs...)

Solve segment `segment` (`1` or `:auto`) of a sequence. Segment 0 is solved
first; only its NIA endpoint starts the re-emergence search.
"""
function _qnm_frequency_segmented(mode::QNMMode, a, segment;
        guess=nothing, kwargs...)
    convention = _qnm_convention(get(kwargs, :convention, :overtone))
    contract = _schwarzschild_label_contract(mode, convention, a)
    contract.near_as_role in (:unconventional, :algebraically_special) &&
        throw(ArgumentError(
            "segment applies to ordinary overtone sequences, not to the " *
            "8_0/8_1 near-algebraically-special pair."))
    # An explicit guess fixes the root; the segment only labels it.
    guess === nothing || return _qnm_frequency_labeled(mode, a;
        guess, segment=segment === :auto ? 0 : segment, kwargs...)
    T = _working_type(a, nothing)
    aT = T(a)
    abs(aT) > one(T) && throw(DomainError(a,
        "A Kerr spin must satisfy abs(a) <= 1."))
    if abs(aT) == one(T)
        segment === :auto && return _qnm_frequency_labeled(mode, aT;
            kwargs...)
        return _qnm_frequency_labeled(mode, aT; segment, kwargs...)
    end
    if aT < zero(T)
        sheet_id = get(kwargs, :sheet_id, :qnm_straight_from_spherical)
        source = _qnm_frequency_segmented(
            QNMMode(mode.s, mode.l, -mode.m, mode.n, mode.branch), -aT,
            segment; kwargs...,
            sheet_id=Symbol(sheet_id, :_positive_spin_source))
        angular_order = get(source.provenance, :angular_order,
            get(kwargs, :angular_order, 32))
        return _spin_reflected_result(mode, aT, source;
            angular_order, sheet_id)
    end

    first_segment = _qnm_frequency_labeled(mode, aT; kwargs...)
    failure(omega, reason, provenance; requested=segment) =
        _segment_failure(mode, aT, omega, reason, provenance;
            convention, inversion_index=first_segment.inversion_index,
            segment=requested)
    if _segment_succeeded(first_segment)
        segment === :auto && return first_segment
        return failure(first_segment.omega, :no_further_segment, (
            segment0_status=first_segment.status,
            segment0_stop_reason=first_segment.stop_reason,
            segment0_omega=first_segment.omega,
        ))
    end
    first_segment.stop_reason == :negative_imaginary_axis_endpoint ||
        return first_segment
    endpoint = (spin=first_segment.provenance.endpoint_spin,
        frequency=first_segment.provenance.endpoint_frequency)
    if mode.m != 0
        segment === :auto && return first_segment
        search = _nia_reemergence_seed(mode, first_segment, aT;
            attribute=false, kwargs...)
        return failure(first_segment.omega,
            :segment_attribution_unavailable, (
                segment0_endpoint=endpoint,
                reemergence_events=search.events,
                searched_spin_range=search.searched_spin_range,
                edge_scan_rows=search.edge_scan_rows,
            ))
    end

    search = _nia_reemergence_seed(mode, first_segment, aT; kwargs...)
    search_provenance = (
        segment0_endpoint=endpoint,
        reemergence_events=search.events,
        searched_spin_range=search.searched_spin_range,
        edge_scan_rows=search.edge_scan_rows,
        attribution_window=search.attribution_window,
    )
    if search.event === nothing
        segment === :auto && return _with_provenance(first_segment,
            merge(search_provenance,
                (segment_search=:negative_imaginary_axis_gap,)))
        return failure(first_segment.omega, :negative_imaginary_axis_gap,
            search_provenance)
    end
    second_segment = _qnm_frequency_labeled(mode, aT;
        kwargs..., segment=1, segment_seed=search.seed)
    event = search.event
    segment_provenance = merge(search_provenance, (
        reemergence=(spin=event.spin, frequency=event.frequency),
        segment_endpoints=(segment0_end=endpoint,
            segment1_start=(spin=event.spin, frequency=event.frequency)),
    ))
    return _with_provenance(second_segment, segment_provenance)
end

# Roots on the edge line Re(omega) = rho_edge at spin `a`, polished and
# restricted to the candidate strip 0 < Re(omega) <= 1/4.
function _nia_edge_scan(mode::QNMMode, a::T, gamma_center, rho_edge,
        polish_options) where {T<:AbstractFloat}
    sign = _branch_sign(mode)
    gammas = range(gamma_center - T(_NIA_SEARCH_HALF_WINDOW),
        gamma_center + T(_NIA_SEARCH_HALF_WINDOW);
        length=_NIA_EDGE_SCAN_POINTS)
    residual_options = (; polish_options.angular_order,
        polish_options.sheet_id, polish_options.inversion_index,
        polish_options.cf_tolerance, polish_options.cf_minimum_iterations,
        polish_options.cf_maximum_iterations)
    points = [complex(sign * rho_edge, -gamma) for gamma in gammas]
    # The continued fraction has poles at the polynomial frequencies -i N/4.
    # A root leaving the axis next to one is masked by it, so the scan uses
    # F * prod(omega - omega_N) over the poles in the window, and every pole
    # also gets Newton starts on the edge line beside it.
    poles = [complex(zero(T), -T(N) / 4) for N in
        ceil(Int, 4 * first(gammas)):floor(Int, 4 * last(gammas))]
    magnitudes = map(points) do omega
        value = try
            abs(leaver_radial_residual(mode, a, omega;
                residual_options..., require_cf_convergence=false).value) *
                prod(abs(omega - pole) for pole in poles; init=one(T))
        catch
            T(Inf)
        end
        isfinite(value) ? T(value) : T(Inf)
    end
    starts = Complex{T}[points[index] for index in 2:(length(points) - 1)
        if magnitudes[index] < magnitudes[index - 1] &&
            magnitudes[index] < magnitudes[index + 1]]
    # A root born from a pole can stay within ~1e-5 of it well after the
    # birth (l=2, m=0, n=19 at a=0.56), so the starts beside each pole close in
    # geometrically from the edge line, rho_edge / 8^k for k = 0..3.
    for pole in poles, k in 0:3
        distance = rho_edge / T(8)^k
        for offset in (-distance, zero(T), distance)
            push!(starts, complex(sign * distance, imag(pole) + offset))
        end
    end
    candidates = NamedTuple[]
    for start in starts
        polished = _polish_frequency_root(mode, a, start;
            polish_options...)
        polished.converged || continue
        rho = sign * real(polished.omega)
        zero(T) < rho <= T(_NIA_CANDIDATE_REAL_CAP) || continue
        duplicate = any(candidates) do candidate
            abs(candidate.omega - polished.omega) <=
                T(1.0e-6) * max(one(T), abs(polished.omega))
        end
        duplicate && continue
        push!(candidates, (omega=polished.omega,
            angular_A=polished.evaluation.angular.angular_A,
            scan_start=start))
    end
    return (candidates=candidates, rows=(a=a, gammas=collect(gammas),
        magnitudes=magnitudes, poles=poles, starts=starts,
        roots=[c.omega for c in candidates]))
end

# Continue a candidate root at `a_start` down to `a_end`. A candidate that
# reaches a_end without touching the axis belongs to another sequence; one
# that ends on the NIA above a_end is a re-emergence event.
function _nia_downward_event(mode::QNMMode, a_start::T, candidate, a_end::T,
        settings) where {T<:AbstractFloat}
    state = (accepted_spins=T[a_start], accepted_coordinates=T[a_start],
        accepted_roots=Complex{T}[candidate.omega],
        accepted_angular_A=Complex{T}[candidate.angular_A],
        continuation_rows=NamedTuple[],
        rejected_continuation_rows=NamedTuple[],
        last_polish=Ref{Any}(nothing))
    outcome = _continue_in_spin!(state, mode, a_end;
        selected_coordinate=:spin, settings.initial_step,
        settings.minimum_step, settings.maximum_step,
        settings.continuation_root_tolerance, settings.polish_options)
    outcome.completed && return (event=false, kind=:continues_below_endpoint,
        lowest_spin=state.accepted_spins[end])
    classified = _classify_continuation_failure(mode, outcome.reason,
        outcome.polished, outcome.failed_spin, state.accepted_spins,
        state.accepted_roots, state.accepted_angular_A,
        state.rejected_continuation_rows; sequence_segment=1)
    classified.reason == :negative_imaginary_axis_endpoint ||
        return (event=false, kind=classified.reason,
            lowest_spin=state.accepted_spins[end])
    spin = classified.provenance.endpoint_spin
    a_end < spin <= a_start || return (event=false,
        kind=:endpoint_outside_search_interval, lowest_spin=spin)
    return (event=true, kind=:negative_imaginary_axis_endpoint, spin=spin,
        frequency=classified.provenance.endpoint_frequency,
        gamma=-imag(classified.provenance.endpoint_frequency),
        lowest_spin=state.accepted_spins[end])
end

"""
    _nia_reemergence_seed(mode, first_segment, target; attribute=true,
        kwargs...)

Search for the segment that re-emerges from the NIA after `first_segment`
ended there. Spin rungs `a_end + j*maximum_spin_step` up to `target` are
scanned on the edge line (first within `a_end + 0.1`, then up to the
sequence anchor); each candidate is continued down to `a_end`. An event is
attributed to `mode` (m = 0) when its axis frequency lies between the
segment-0 endpoints of overtones `n` and `n+1`; the earliest such event seeds
the segment-1 continuation at its rung.
"""
function _nia_reemergence_seed(mode::QNMMode, first_segment, target::T;
        attribute::Bool=true, kwargs...) where {T<:AbstractFloat}
    provenance = first_segment.provenance
    settings = provenance.continuation_settings
    a_end = T(provenance.endpoint_spin)
    omega_end = provenance.endpoint_frequency
    gamma_end = -imag(omega_end)
    rho_edge = T(_NIA_EDGE_STEP_FACTOR) *
        T(_root_derivative_step(T, omega_end))
    step = settings.maximum_step
    anchor = _spin_sequence_anchor(mode, T)
    first_upper = min(target, a_end + T(_NIA_FIRST_SEARCH_SPAN), anchor)
    final_upper = min(target, anchor)
    rungs = T[]
    j = 1
    while a_end + j * step <= final_upper
        push!(rungs, a_end + j * step)
        j += 1
    end
    final_upper > a_end && (isempty(rungs) || rungs[end] < final_upper) &&
        push!(rungs, final_upper)
    upper_gamma = Ref{Union{Nothing,T}}(nothing)
    function next_endpoint_gamma()
        upper_gamma[] === nothing || return upper_gamma[]
        next = _qnm_frequency_labeled(
            QNMMode(mode.s, mode.l, mode.m, mode.n + 1, mode.branch),
            anchor; kwargs...)
        upper_gamma[] = next.stop_reason ==
            :negative_imaginary_axis_endpoint ?
            T(-imag(next.provenance.endpoint_frequency)) : T(Inf)
        return upper_gamma[]
    end
    events = NamedTuple[]
    rows = NamedTuple[]
    searched = (a_end, a_end)
    for (pass, upper) in ((1, first_upper), (2, final_upper))
        for rung in rungs
            rung <= upper || break
            pass == 2 && rung <= first_upper && continue
            searched = (a_end, rung)
            scan = _nia_edge_scan(mode, rung, gamma_end, rho_edge,
                settings.polish_options)
            push!(rows, scan.rows)
            rung_events = NamedTuple[]
            for candidate in scan.candidates
                outcome = _nia_downward_event(mode, rung, candidate, a_end,
                    settings)
                record = merge(outcome, (rung=rung,
                    rung_frequency=candidate.omega))
                push!(events, record)
                outcome.event && push!(rung_events, merge(record,
                    (seed=candidate,)))
            end
            attribute || continue
            attributed = filter(rung_events) do event
                gamma_end < event.gamma < next_endpoint_gamma()
            end
            isempty(attributed) && continue
            event = attributed[argmin([e.spin for e in attributed])]
            seed = (
                seed=event.seed.omega,
                seed_spin=rung,
                source=:nia_reemergence_edge_scan,
                seed_inversion_index=first_segment.inversion_index,
                seed_predictor=:edge_scan_polished_root,
                seed_duplicate_distance=T(Inf),
            )
            return (event=event, seed=seed, events=events,
                searched_spin_range=searched, edge_scan_rows=rows,
                attribution_window=(gamma_end, upper_gamma[]))
        end
    end
    return (event=nothing, seed=nothing, events=events,
        searched_spin_range=searched, edge_scan_rows=rows,
        attribution_window=(gamma_end, upper_gamma[]))
end


"""
    _near_axis_precision_check(mode, a, omega, polish_options; precision_bits,
        root_tolerance, stability_tolerance)

Re-polish a Float64 root next to the NIA in BigFloat at the same spin, starting
from `omega`, with the same angular order and continued-fraction settings (the
ceiling raised as for the automatic precision fallback). The root is verified
when the BigFloat polish converges and moves it by at most the stability gate
relative to max(1, |omega|).
"""
function _near_axis_precision_check(mode::QNMMode, a::Float64,
        omega::ComplexF64, polish_options; precision_bits::Int,
        root_tolerance, stability_tolerance)
    return setprecision(BigFloat, precision_bits) do
        B = BigFloat
        options = merge(polish_options, (
            sheet_id=Symbol(polish_options.sheet_id, :_near_axis_check),
            cf_tolerance=B(polish_options.cf_tolerance),
            cf_minimum_iterations=max(polish_options.cf_minimum_iterations,
                600),
            cf_maximum_iterations=max(polish_options.cf_maximum_iterations,
                24000),
            root_tolerance=B(root_tolerance),
            trust_radius=B(polish_options.trust_radius),
        ))
        polished = try
            _polish_frequency_root(mode, B(a), _complex_of_type(B, omega);
                options...)
        catch
            nothing
        end
        polished === nothing && return (verified=false,
            precision_bits, converged=false, stop_reason=:evaluation_error,
            bigfloat_omega=missing, drift=Inf)
        drift = Float64(abs(polished.omega - omega) /
            max(one(B), abs(polished.omega)))
        return (
            verified=polished.converged && drift <= stability_tolerance,
            precision_bits,
            converged=polished.converged,
            stop_reason=polished.stop_reason,
            bigfloat_omega=ComplexF64(polished.omega),
            drift,
        )
    end
end
