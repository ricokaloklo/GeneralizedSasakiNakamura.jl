module OrdinaryPointExpansion

using ..CoefficientTables:
    CoefficientSet,
    EndpointCoefficientSet,
    RationalCoefficients,
    evaluate_ordinary_coefficients!,
    evaluate_ordinary_variable!

export PQLayout
export pq_layout, split_starts, ab_term_counts
export rational_series_from_coefficients!
export endpoint_ab_series!, endpoint_ab_series
export ordinary_ab_series!, ordinary_ab_series
export ordinary_scaled_y_ab_series!
export shifted_ab!, shifted_infinity_ab!
export ordinary_point_coeffs!
export endpoint_power_lookup
export poly_value, poly_pair, poly_triple
export series_issue, denom_ok

const _SERIES_MAG_LIMIT = sqrt(floatmax(Float64))
const _SERIES_RATIO_LIMIT = inv(sqrt(eps(Float64)))
const _DENOM_REL_LIMIT = 64eps(Float64)
const _SHIFT_EPS_LIMIT = 1e-13

@inline _finite_complex(z) = isfinite(real(z)) && isfinite(imag(z))

# Check the newest series coefficient for non-finite values, over/underflow, or
# two consecutive coefficient ratios above 1/sqrt(eps); returns :none if it is usable.
function series_issue(coeffs, index::Integer, scale::Real=1.0)
    coeff = coeffs[index]
    _finite_complex(coeff) || return :nonfinite_coefficient

    magnitude = abs(coeff) / max(abs(scale), floatmin(Float64))
    !isfinite(magnitude) && return :coefficient_too_large
    magnitude > _SERIES_MAG_LIMIT && return :coefficient_too_large
    !iszero(coeff) && iszero(magnitude) && return :coefficient_too_small
    0 < magnitude < floatmin(Float64) && return :coefficient_too_small

    if index >= 3
        previous2 = abs(coeffs[index - 2])
        previous1 = abs(coeffs[index - 1])
        current = abs(coeff)
        if previous2 > 0 && previous1 > 0 && current > 0
            ratio1 = previous1 / previous2
            ratio2 = current / previous1
            growth1 = !isfinite(ratio1) || ratio1 > _SERIES_RATIO_LIMIT
            growth2 = !isfinite(ratio2) || ratio2 > _SERIES_RATIO_LIMIT
            if growth1 && growth2
                return :consecutive_ratio_growth
            end
        end
    end
    return :none
end

# A recurrence denominator is usable if finite and not small relative to `scale`.
function denom_ok(denominator, scale)
    _finite_complex(denominator) || return false
    denom_scale = max(abs(scale), floatmin(Float64))
    return abs(denominator) > _DENOM_REL_LIMIT * denom_scale
end

function _q0_ok(q0)
    _finite_complex(q0) && !iszero(q0) || return false
    return _finite_complex(inv(ComplexF64(q0)))
end

# Lengths of the A and B numerator/denominator polynomials (ap, aq, bp, bq) in a packed P/Q vector.
struct PQLayout
    ap_len::Int
    aq_len::Int
    bp_len::Int
    bq_len::Int
end

# Block lengths of the packed A/B numerator and denominator coefficients for the
# horizon (:H), infinity (:I) and ordinary-point (:O) tables.
function pq_layout(s::Integer, group::Symbol)
    if s == 0
        group == :H && return PQLayout(11, 11, 11, 13)
        group == :I && return PQLayout(10, 10, 9, 9)
        group == :O && return PQLayout(10, 11, 9, 13)
    else
        group == :H && return PQLayout(11, 11, 17, 19)
        group == :I && return PQLayout(10, 10, 15, 15)
        group == :O && return PQLayout(10, 11, 15, 19)
    end
    throw(ArgumentError("unsupported ISEM GSN coefficient group $group."))
end

# 1-based start offsets of the four polynomials in a packed P/Q vector.
function split_starts(layout::PQLayout)
    ap = 1
    aq = ap + layout.ap_len
    bp = aq + layout.aq_len
    bq = bp + layout.bp_len
    return ap, aq, bp, bq
end

# Number of Laurent/Taylor terms of A and B needed for a solution series of
# `order`; endpoint series start at x^-1 and x^-2 (horizon) or y^-1 and y^-4 (infinity).
function ab_term_counts(group::Symbol, order::Integer)
    order >= 0 || throw(ArgumentError("series order must be non-negative."))
    group == :H && return order + 2, order + 3
    group == :I && return order + 2, order + 5
    return order + 1, order + 1
end

function _denominator_leading_shift(values)
    @inbounds for k in 0:(length(values) - 1)
        !iszero(values[k + 1]) && return k
    end
    error("all-zero denominator in ISEM GSN rational P/Q recurrence")
end

function _check_regular_leading_terms(numerator, shift::Integer)
    @inbounds for k in 0:(shift - 1)
        if k < length(numerator) && !iszero(numerator[k + 1])
            error("non-regular ISEM GSN rational P/Q leading terms")
        end
    end
    return nothing
end

function _ordinary_rational_series_from_pq_len!(
    out::Vector{ComplexF64},
    out_len::Integer,
    pq::Vector{ComplexF64},
    p_start::Integer,
    p_len::Integer,
    q_start::Integer,
    q_len::Integer,
)
    out_len <= length(out) || throw(ArgumentError("output buffer is too short."))
    q0 = pq[q_start]
    _q0_ok(q0) || error("noninvertible q0 in ISEM GSN ordinary P/Q recurrence")
    inv_q0 = inv(q0)
    order = out_len - 1
    @inbounds for n in 0:order
        acc = n < p_len ? pq[p_start + n] : 0.0 + 0.0im
        for k in 1:min(n, q_len - 1)
            acc -= pq[q_start + k] * out[n - k + 1]
        end
        out[n + 1] = acc * inv_q0
    end
    return out
end

@inline function _shift_source(source, pole_order, index)
    index <= pole_order && return 0.0 + 0.0im
    return ComplexF64(source[index - pole_order])
end

# Taylor shift of a polynomial (times x^pole_order) to the point `center`,
# written into `work` starting at `start`; returns the number of coefficients.
function _shift_poly!(
    work::Vector{<:Complex},
    start::Integer,
    source,
    pole_order::Integer,
    center,
)
    count = length(source) + pole_order
    start + count - 1 <= length(work) ||
        throw(ArgumentError("shift workspace is too short."))
    work[start] = _shift_source(source, pole_order, count)
    active = 1
    @inbounds for source_index in (count - 1):-1:1
        work[start + active] = work[start + active - 1]
        for index in active:-1:2
            slot = start + index - 1
            work[slot] = work[slot - 1] + center * work[slot]
        end
        work[start] = _shift_source(source, pole_order, source_index) +
            center * work[start]
        active += 1
    end
    return count
end

# Condition number of the Taylor shift: sum of absolute contributions over the
# magnitude of each shifted coefficient.
function _shift_condition(
    work::Vector{<:Complex},
    start::Integer,
    source,
    pole_order::Integer,
    center,
)
    count = length(source) + pole_order
    degree = count - 1
    maximum_condition = 1.0
    @inbounds for n in 0:degree
        scale = 0.0
        weight = 1.0
        for k in n:degree
            scale += abs(_shift_source(source, pole_order, k + 1)) * weight
            k < degree &&
                (weight *= abs(center) * (k + 1) / (k + 1 - n))
        end
        shifted = work[start + n]
        condition = scale / max(abs(shifted), floatmin(Float64))
        isfinite(condition) || return Inf
        maximum_condition = max(maximum_condition, condition)
    end
    return maximum_condition
end

# Series division of the shifted numerator by the shifted denominator; returns
# the largest cancellation ratio met in the recurrence.
function _shift_series!(
    out::Vector{<:Complex},
    out_len::Integer,
    work::Vector{<:Complex},
    p_start::Integer,
    p_len::Integer,
    q_start::Integer,
    q_len::Integer,
)
    q0 = work[q_start]
    _q0_ok(q0) || return Inf
    inv_q0 = inv(q0)
    maximum_condition = 1.0
    @inbounds for n in 0:(out_len - 1)
        acc = n < p_len ? work[p_start + n] : 0.0 + 0.0im
        scale = abs(acc)
        for k in 1:min(n, q_len - 1)
            term = work[q_start + k] * out[n - k + 1]
            acc -= term
            scale += abs(term)
        end
        out[n + 1] = acc * inv_q0
        condition = scale / max(abs(acc), floatmin(Float64))
        isfinite(condition) || return Inf
        maximum_condition = max(maximum_condition, condition)
    end
    return maximum_condition
end

function _shift_rational!(
    out::Vector{<:Complex},
    out_len::Integer,
    work::Vector{<:Complex},
    rational::RationalCoefficients,
    pole_order::Integer,
    center,
)
    p_start = 1
    p_len = _shift_poly!(work, p_start, rational.numerator, 0, center)
    p_condition = _shift_condition(
        work, p_start, rational.numerator, 0, center)
    q_start = p_start + p_len
    q_len = _shift_poly!(work, q_start, rational.denominator,
        pole_order, center)
    q_condition = _shift_condition(
        work, q_start, rational.denominator, pole_order, center)
    recurrence_condition = _shift_series!(
        out, out_len, work, p_start, p_len, q_start, q_len)
    return max(p_condition, q_condition, recurrence_condition)
end

# Shift in binary64, falling back to 128-bit BigFloat when eps * condition
# exceeds `eps_limit`; returns the estimated relative error.
function _accurate_shift_rational!(out, terms, work, rational, pole_order, center, eps_limit)
    condition = _shift_rational!(out, terms, work, rational, pole_order, center)
    error = eps(Float64) * condition
    error <= eps_limit && return error
    # Only the polynomial shift and series division need extra precision.
    # The returned local coefficients and propagated solution remain binary64.
    return setprecision(BigFloat, 128) do
        big_out = Vector{Complex{BigFloat}}(undef, terms)
        big_work = Vector{Complex{BigFloat}}(undef, length(work))
        condition = _shift_rational!(big_out, terms, big_work,
            rational, pole_order, Complex{BigFloat}(center))
        copyto!(out, 1, big_out, 1, terms)
        Float64(eps(BigFloat) * condition + eps(Float64))
    end
end

@inline function _valid_shift_center(center)
    finite = isfinite(real(center)) && isfinite(imag(center))
    finite || return false
    if center isa Real
        return 0 < center < 1
    end
    # A rational Taylor shift is not restricted to the unit disk. Its
    # denominator and conditioning checks below determine local validity.
    return !iszero(center) && !isone(center)
end

# Taylor coefficients of A and B about x0 obtained by shifting the endpoint
# rational functions. Returns (accepted, A terms, B terms, error estimate).
function shifted_ab!(
    a_out::Vector{ComplexF64},
    b_out::Vector{ComplexF64},
    work::Vector{ComplexF64},
    endpoint::EndpointCoefficientSet,
    x0,
    order::Integer;
    eps_limit::Real=_SHIFT_EPS_LIMIT,
)
    center = x0 isa Real ? Float64(x0) : ComplexF64(x0)
    _valid_shift_center(center) ||
        return false, order + 1, order + 1, Inf
    a_terms, b_terms = ab_term_counts(:O, order)
    length(a_out) >= a_terms ||
        throw(ArgumentError("A-series output buffer is too short."))
    length(b_out) >= b_terms ||
        throw(ArgumentError("B-series output buffer is too short."))
    a_error = _accurate_shift_rational!(
        a_out, a_terms, work, endpoint.A, 1, center, eps_limit)
    b_error = _accurate_shift_rational!(
        b_out, b_terms, work, endpoint.B, 2, center, eps_limit)
    eps_condition = max(a_error, b_error)
    accepted = isfinite(eps_condition) && eps_condition <= eps_limit
    return accepted, a_terms, b_terms, eps_condition
end

# As shifted_ab!, for the infinity rational functions about y0; A is negated
# so that the result is for derivatives in y.
function shifted_infinity_ab!(
    a_out::Vector{ComplexF64},
    b_out::Vector{ComplexF64},
    work::Vector{ComplexF64},
    endpoint::EndpointCoefficientSet,
    y0,
    order::Integer;
    eps_limit::Real=_SHIFT_EPS_LIMIT,
)
    center = y0 isa Real ? Float64(y0) : ComplexF64(y0)
    _valid_shift_center(center) ||
        return false, order + 1, order + 1, Inf
    a_terms, b_terms = ab_term_counts(:O, order)
    length(a_out) >= a_terms ||
        throw(ArgumentError("A-series output buffer is too short."))
    length(b_out) >= b_terms ||
        throw(ArgumentError("B-series output buffer is too short."))
    a_error = _accurate_shift_rational!(
        a_out, a_terms, work, endpoint.A, 1, center, eps_limit)
    b_error = _accurate_shift_rational!(
        b_out, b_terms, work, endpoint.B, 4, center, eps_limit)
    eps_condition = max(a_error, b_error)
    accepted = isfinite(eps_condition) && eps_condition <= eps_limit
    if accepted
        @inbounds for index in 1:a_terms
            a_out[index] = -a_out[index]
        end
    end
    return accepted, a_terms, b_terms, eps_condition
end

# Power series of numerator/denominator by recursive division into the first
# `out_len` entries of `out`, after removing a common leading power of the
# argument from the denominator.
function rational_series_from_coefficients!(
    out::Vector{ComplexF64},
    out_len::Integer,
    numerator,
    denominator,
)
    out_len <= length(out) || throw(ArgumentError("output buffer is too short."))
    shift = _denominator_leading_shift(denominator)
    _check_regular_leading_terms(numerator, shift)
    q0 = denominator[shift + 1]
    _q0_ok(q0) || error("noninvertible q0 in ISEM GSN rational P/Q recurrence")
    inv_q0 = inv(q0)
    order = out_len - 1
    @inbounds for n in 0:order
        p_index = n + shift
        acc = p_index < length(numerator) ? numerator[p_index + 1] : 0.0 + 0.0im
        for k in 1:min(n, length(denominator) - shift - 1)
            acc -= denominator[shift + k + 1] * out[n - k + 1]
        end
        out[n + 1] = acc * inv_q0
    end
    return out
end

# Laurent series of A and B about the horizon (:H, in x) or infinity (:I, in y),
# written into `a_out` and `b_out`; returns the numbers of terms.
function endpoint_ab_series!(
    a_out::Vector{ComplexF64},
    b_out::Vector{ComplexF64},
    endpoint::EndpointCoefficientSet,
    group::Symbol,
    order::Integer,
)
    a_terms, b_terms = ab_term_counts(group, order)
    length(a_out) >= a_terms || throw(ArgumentError("A-series output buffer is too short."))
    length(b_out) >= b_terms || throw(ArgumentError("B-series output buffer is too short."))
    rational_series_from_coefficients!(a_out, a_terms, endpoint.A.numerator, endpoint.A.denominator)
    rational_series_from_coefficients!(b_out, b_terms, endpoint.B.numerator, endpoint.B.denominator)
    return a_terms, b_terms
end

# Allocating form of endpoint_ab_series!; returns (A series, B series).
function endpoint_ab_series(coefficients::CoefficientSet, group::Symbol, order::Integer)
    group in (:H, :I) || throw(ArgumentError("endpoint group must be :H or :I."))
    endpoint = group == :H ? coefficients.horizon : coefficients.infinity
    a_terms, b_terms = ab_term_counts(group, order)
    a_coeffs = Vector{ComplexF64}(undef, a_terms)
    b_coeffs = Vector{ComplexF64}(undef, b_terms)
    endpoint_ab_series!(a_coeffs, b_coeffs, endpoint, group, order)
    return a_coeffs, b_coeffs
end

# Taylor series of A and B from P/Q values already evaluated at the expansion point.
function ordinary_ab_series!(
    a_out::Vector{ComplexF64},
    b_out::Vector{ComplexF64},
    pq::Vector{ComplexF64},
    s::Integer,
    order::Integer,
)
    layout = pq_layout(s, :O)
    ap, aq, bp, bq = split_starts(layout)
    a_terms, b_terms = ab_term_counts(:O, order)
    length(a_out) >= a_terms || throw(ArgumentError("A-series output buffer is too short."))
    length(b_out) >= b_terms || throw(ArgumentError("B-series output buffer is too short."))
    _ordinary_rational_series_from_pq_len!(a_out, a_terms, pq, ap, layout.ap_len, aq, layout.aq_len)
    _ordinary_rational_series_from_pq_len!(b_out, b_terms, pq, bp, layout.bp_len, bq, layout.bq_len)
    return a_terms, b_terms
end

# A and B series in the scaled variable u = (1 - x) / omega_scale about u = y0
# (y0 is given in scaled units); the P/Q blocks are rescaled before division.
function ordinary_scaled_y_ab_series!(a_out, b_out, pq,
        coefficients::CoefficientSet, y0, omega_scale, order::Integer)
    variable = coefficients.ordinary.variable == :y ?
        omega_scale * y0 : 1.0 - omega_scale * y0
    evaluate_ordinary_variable!(pq, coefficients.ordinary, variable)
    layout = pq_layout(coefficients.params.s, :O)
    starts = split_starts(layout)
    lengths = (layout.ap_len, layout.aq_len, layout.bp_len, layout.bq_len)
    # Rescale the rational polynomials before division; the unscaled Taylor
    # coefficients can overflow even when the scaled equation is regular.
    for (start, count) in zip(starts, lengths)
        power = 1.0
        for j in 0:(count-1)
            pq[start+j] *= power
            power *= -omega_scale
        end
    end
    na, nb = ordinary_ab_series!(
        a_out, b_out, pq, coefficients.params.s, order)
    for j in 1:na
        a_out[j] *= -omega_scale
    end
    for j in 1:nb
        b_out[j] *= omega_scale^2
    end
    return na, nb
end

# Taylor series of A and B about the ordinary point x0: by shifting the horizon
# rational functions when `metadata.horizon_shift` is set and the shift is
# accepted, otherwise from the ordinary-point tables.
function ordinary_ab_series!(
    a_out::Vector{ComplexF64},
    b_out::Vector{ComplexF64},
    pq::Vector{ComplexF64},
    coefficients::CoefficientSet,
    x0,
    order::Integer,
)
    use_shift = hasproperty(coefficients.metadata, :horizon_shift) &&
        getproperty(coefficients.metadata, :horizon_shift)
    if use_shift
        accepted, a_terms, b_terms, _ = shifted_ab!(
            a_out, b_out, pq, coefficients.horizon, x0, order)
        accepted && return a_terms, b_terms
    end
    evaluate_ordinary_coefficients!(pq, coefficients.ordinary, x0)
    return ordinary_ab_series!(a_out, b_out, pq, coefficients.params.s, order)
end

# Allocating form of ordinary_ab_series!; returns (A series, B series, P/Q work vector).
function ordinary_ab_series(coefficients::CoefficientSet, x0, order::Integer)
    a_terms, b_terms = ab_term_counts(:O, order)
    pq = Vector{ComplexF64}(undef, coefficients.ordinary.value_count)
    a_coeffs = Vector{ComplexF64}(undef, a_terms)
    b_coeffs = Vector{ComplexF64}(undef, b_terms)
    ordinary_ab_series!(a_coeffs, b_coeffs, pq, coefficients, x0, order)
    return a_coeffs, b_coeffs, pq
end

# Taylor coefficients of a solution about an ordinary point from the A and B
# series and the initial value c0 and slope c1.
function ordinary_point_coeffs!(
    coeffs::Vector{ComplexF64},
    a_coeffs,
    b_coeffs,
    c0,
    c1,
    order::Integer,
)
    order >= 0 || throw(ArgumentError("series order must be non-negative."))
    length(coeffs) >= order + 1 || throw(ArgumentError("solution coefficient buffer is too short."))
    coeffs[1] = ComplexF64(c0)
    order >= 1 && (coeffs[2] = ComplexF64(c1))
    order < 2 && return coeffs
    @inbounds for n in 0:(order - 2)
        numerator = 0.0 + 0.0im
        for j in 0:n
            numerator += a_coeffs[j + 1] * (n - j + 1) * coeffs[n - j + 2]
            numerator += b_coeffs[j + 1] * coeffs[n - j + 1]
        end
        coeffs[n + 3] = -numerator / ((n + 2) * (n + 1))
    end
    return coeffs
end

# Closure p -> coefficient of power p, for a series whose first stored power is
# `first_power`; zero outside the stored range.
function endpoint_power_lookup(coeffs, first_power::Integer)
    return p -> begin
        idx = p - first_power + 1
        1 <= idx <= length(coeffs) ? coeffs[idx] : 0.0 + 0.0im
    end
end

# Horner evaluation of a polynomial, with first (and second) derivatives in the
# _pair and _triple versions.
function poly_value(coeffs, z)
    isempty(coeffs) && return 0.0 + 0.0im
    value = coeffs[end]
    @inbounds for i in (length(coeffs) - 1):-1:1
        value = value * z + coeffs[i]
    end
    return value
end

function poly_pair(coeffs, z)
    n = length(coeffs) - 1
    n < 0 && return 0.0 + 0.0im, 0.0 + 0.0im
    n == 0 && return coeffs[1], 0.0 + 0.0im
    value = coeffs[n + 1]
    derivative = 0.0 + 0.0im
    @inbounds for i in n:-1:1
        derivative = derivative * z + value
        value = value * z + coeffs[i]
    end
    return value, derivative
end

function poly_triple(coeffs, z)
    n = length(coeffs) - 1
    n < 0 && return 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im
    n == 0 && return coeffs[1], 0.0 + 0.0im, 0.0 + 0.0im
    value = coeffs[n + 1]
    derivative = 0.0 + 0.0im
    second = 0.0 + 0.0im
    @inbounds for i in n:-1:1
        second = second * z + 2 * derivative
        derivative = derivative * z + value
        value = value * z + coeffs[i]
    end
    return value, derivative, second
end

end
