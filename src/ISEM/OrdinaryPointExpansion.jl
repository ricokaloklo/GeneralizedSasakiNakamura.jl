module OrdinaryPointExpansion

using ..CoefficientTables:
    CoefficientSet,
    EndpointCoefficientSet,
    RationalCoefficients,
    evaluate_ordinary_coefficients!,
    evaluate_ordinary_variable!
using ...Potentials: sF, sU
using TaylorSeries: Taylor1, getcoeff, derivative

export PQLayout
export pq_layout, split_starts, ab_term_counts
export rational_series_from_coefficients!
export endpoint_ab_series!, endpoint_ab_series
export ordinary_ab_series!, ordinary_ab_series
export ordinary_scaled_y_ab_series!
export shifted_ab!, shifted_infinity_ab!
export ordinary_point_coeffs!
export DIRECT_AB_KEY, direct_ab_series!, direct_ab_value
export endpoint_power_lookup
export poly_value, poly_pair, poly_triple
export series_issue, denom_ok

const _SERIES_MAG_LIMIT = sqrt(floatmax(Float64))
const _SERIES_RATIO_LIMIT = inv(sqrt(eps(Float64)))
const _DENOM_REL_LIMIT = 64eps(Float64)
const _SHIFT_EPS_LIMIT = 1e-13

# The same guard limits in the working floating type; Float64 uses the constants above.
@inline _series_mag_limit(::Type{Float64}) = _SERIES_MAG_LIMIT
@inline _series_mag_limit(::Type{T}) where {T<:AbstractFloat} = sqrt(floatmax(T))
@inline _series_ratio_limit(::Type{Float64}) = _SERIES_RATIO_LIMIT
@inline _series_ratio_limit(::Type{T}) where {T<:AbstractFloat} = inv(sqrt(eps(T)))
@inline _denom_rel_limit(::Type{Float64}) = _DENOM_REL_LIMIT
@inline _denom_rel_limit(::Type{T}) where {T<:AbstractFloat} = 64eps(T)
@inline _working_float(value) = typeof(float(real(value)))

@inline _finite_complex(z) = isfinite(real(z)) && isfinite(imag(z))

# Check the newest series coefficient for non-finite values, over/underflow, or
# two consecutive coefficient ratios above 1/sqrt(eps); returns :none if it is usable.
function series_issue(coeffs, index::Integer, scale::Real=1.0)
    coeff = coeffs[index]
    _finite_complex(coeff) || return :nonfinite_coefficient
    T = _working_float(coeff)

    magnitude = abs(coeff) / max(abs(scale), floatmin(T))
    !isfinite(magnitude) && return :coefficient_too_large
    magnitude > _series_mag_limit(T) && return :coefficient_too_large
    !iszero(coeff) && iszero(magnitude) && return :coefficient_too_small
    0 < magnitude < floatmin(T) && return :coefficient_too_small

    if index >= 3
        previous2 = abs(coeffs[index - 2])
        previous1 = abs(coeffs[index - 1])
        current = abs(coeff)
        if previous2 > 0 && previous1 > 0 && current > 0
            ratio1 = previous1 / previous2
            ratio2 = current / previous1
            growth1 = !isfinite(ratio1) || ratio1 > _series_ratio_limit(T)
            growth2 = !isfinite(ratio2) || ratio2 > _series_ratio_limit(T)
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
    T = _working_float(denominator)
    denom_scale = max(abs(scale), floatmin(T))
    return abs(denominator) > _denom_rel_limit(T) * denom_scale
end

function _q0_ok(q0::Union{Real,Complex})
    _finite_complex(q0) && !iszero(q0) || return false
    return _finite_complex(inv(complex(float(q0))))
end

# Double-double work values of the ill-conditioned shifts: the binary64 inverse is checked.
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
    out::AbstractVector{<:Complex},
    out_len::Integer,
    pq::AbstractVector{<:Complex},
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
        acc = n < p_len ? pq[p_start + n] : zero(eltype(out))
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
    work::Vector{<:Number},
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
    work::Vector{<:Number},
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
    out::Vector{<:Number},
    out_len::Integer,
    work::Vector{<:Number},
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
    out::Vector{<:Number},
    out_len::Integer,
    work::Vector{<:Number},
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

# Complex double-double number (about 106 bits) for the ill-conditioned shifts:
# each part is an unevaluated sum high + low of binary64 values.
struct _DoubleDoubleComplex <: Number
    re_high::Float64
    re_low::Float64
    im_high::Float64
    im_low::Float64
end

const _DDC = _DoubleDoubleComplex
const _DD_EPS = 2.0^-104

@inline function _two_sum(a::Float64, b::Float64)
    s = a + b
    v = s - a
    return s, (a - (s - v)) + (b - v)
end
@inline function _quick_two_sum(a::Float64, b::Float64)
    s = a + b
    return s, b - (s - a)
end
@inline function _dd_add(ah, al, bh, bl)
    s, e = _two_sum(ah, bh)
    return _quick_two_sum(s, e + al + bl)
end
@inline function _dd_mul(ah, al, bh, bl)
    p = ah * bh
    e = fma(ah, bh, -p)
    return _quick_two_sum(p, e + ah * bl + al * bh)
end
@inline function _dd_div(ah, al, bh, bl)
    q1 = ah / bh
    ph, pl = _dd_mul(q1, 0.0, bh, bl)
    rh, rl = _dd_add(ah, al, -ph, -pl)
    q2 = rh / bh
    ph, pl = _dd_mul(q2, 0.0, bh, bl)
    rh, rl = _dd_add(rh, rl, -ph, -pl)
    return _dd_add(q1, 0.0, q2, rh / bh)
end

_DDC(z::Complex) = _DDC(Float64(real(z)), 0.0, Float64(imag(z)), 0.0)
_DDC(x::Real) = _DDC(Float64(x), 0.0, 0.0, 0.0)
Base.convert(::Type{_DDC}, z::Number) = _DDC(complex(z))
Base.convert(::Type{_DDC}, z::_DDC) = z
Base.promote_rule(::Type{_DDC}, ::Type{<:Number}) = _DDC
Base.ComplexF64(z::_DDC) = complex(z.re_high + z.re_low, z.im_high + z.im_low)
Base.convert(::Type{ComplexF64}, z::_DDC) = ComplexF64(z)
Base.real(z::_DDC) = z.re_high
Base.imag(z::_DDC) = z.im_high
Base.zero(::Type{_DDC}) = _DDC(0.0, 0.0, 0.0, 0.0)
Base.iszero(z::_DDC) = iszero(z.re_high) && iszero(z.im_high)
Base.abs(z::_DDC) = hypot(z.re_high, z.im_high)
Base.:+(a::_DDC, b::_DDC) = _DDC(_dd_add(a.re_high, a.re_low, b.re_high, b.re_low)...,
    _dd_add(a.im_high, a.im_low, b.im_high, b.im_low)...)
Base.:-(a::_DDC) = _DDC(-a.re_high, -a.re_low, -a.im_high, -a.im_low)
Base.:-(a::_DDC, b::_DDC) = a + (-b)
function Base.:*(a::_DDC, b::_DDC)
    rh, rl = _dd_mul(a.re_high, a.re_low, b.re_high, b.re_low)
    sh, sl = _dd_mul(a.im_high, a.im_low, b.im_high, b.im_low)
    th, tl = _dd_mul(a.re_high, a.re_low, b.im_high, b.im_low)
    uh, ul = _dd_mul(a.im_high, a.im_low, b.re_high, b.re_low)
    return _DDC(_dd_add(rh, rl, -sh, -sl)..., _dd_add(th, tl, uh, ul)...)
end
function Base.inv(b::_DDC)
    nh, nl = _dd_add(_dd_mul(b.re_high, b.re_low, b.re_high, b.re_low)...,
        _dd_mul(b.im_high, b.im_low, b.im_high, b.im_low)...)
    return _DDC(_dd_div(b.re_high, b.re_low, nh, nl)...,
        _dd_div(-b.im_high, -b.im_low, nh, nl)...)
end

# Shift in binary64, falling back to double-double and then to 128-bit BigFloat
# when eps * condition exceeds `eps_limit`; returns the estimated relative error.
function _accurate_shift_rational!(out, terms, work, rational, pole_order, center, eps_limit)
    condition = _shift_rational!(out, terms, work, rational, pole_order, center)
    error = eps(Float64) * condition
    error <= eps_limit && return error
    dd_out = Vector{_DDC}(undef, terms)
    dd_work = Vector{_DDC}(undef, length(work))
    condition = _shift_rational!(dd_out, terms, dd_work,
        rational, pole_order, _DDC(complex(center)))
    error = _DD_EPS * condition + eps(Float64)
    if error <= eps_limit
        copyto!(out, 1, dd_out, 1, terms)
        return error
    end
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
    out::AbstractVector{<:Complex},
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
        acc = p_index < length(numerator) ? numerator[p_index + 1] : zero(eltype(out))
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
    a_out::AbstractVector{<:Complex},
    b_out::AbstractVector{<:Complex},
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
    C = eltype(endpoint.A.numerator)
    a_coeffs = Vector{C}(undef, a_terms)
    b_coeffs = Vector{C}(undef, b_terms)
    endpoint_ab_series!(a_coeffs, b_coeffs, endpoint, group, order)
    return a_coeffs, b_coeffs
end

# Taylor series of A and B from P/Q values already evaluated at the expansion point.
function ordinary_ab_series!(
    a_out::AbstractVector{<:Complex},
    b_out::AbstractVector{<:Complex},
    pq::AbstractVector{<:Complex},
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

# Below this x the ordinary-point tables lose accuracy: they form B from O(1)
# terms that cancel to O(x^2) near the horizon, so B carries a relative error
# of about 1e-12 / x^2 (4e-8 at x = 0.0044, 9e-3 at x = 1e-5), while the Taylor
# shift of the horizon rational functions stays at about 2e-13 there. The
# horizon series of Y starts at a small x when a root of eta lies near the
# horizon or the spin is near extremal, and its first patches then took the
# table values: Y(r) of (s, l, m) = (-2, 2, 2) was off by 1.9e-5 at
# (a, omega) = (0.8, 0.9) and 1.0e-6 at (0.99, 0.76) without any warning. The
# shift loses accuracy towards x = 1 (relative B error 2e-10 at 0.5, 1e-6 at
# 0.9), so it is used only below this point; up to 0.25 it made Y of
# (0.95, 0.9) jump by 4e-9 between neighbouring omega (1e-10 with the tables).
const HORIZON_SHIFT_MAX_X = 0.1

# Task-local switch (set by the Y self-check re-solve): real points between
# HORIZON_SHIFT_MAX_X and DIRECT_AB_MAX_X take A and B from direct_ab_series!.
const DIRECT_AB_KEY = :GSN_ISEM_DIRECT_AB
const DIRECT_AB_MAX_X = 0.7

# r as a truncated series with its powers r^2 ... r^n computed once: the closed-form potentials take
# powers of r dozens of times (sU of s = -2, order 41: 135 us per call with the powers taken each time,
# 62 us cached).
struct _PowerCachedSeries{T}
    powers::Vector{T}
end

function _PowerCachedSeries(r, n::Integer)
    powers = [r]
    for _ in 2:n
        push!(powers, powers[end] * r)
    end
    return _PowerCachedSeries(powers)
end

Base.literal_pow(::typeof(^), c::_PowerCachedSeries, ::Val{k}) where {k} =
    k <= length(c.powers) ? c.powers[k] : c.powers[1]^k
Base.:-(c::_PowerCachedSeries) = -c.powers[1]
for op in (:+, :-, :*, :/)
    @eval Base.$op(c::_PowerCachedSeries, y) = $op(c.powers[1], y)
    @eval Base.$op(y, c::_PowerCachedSeries) = $op(y, c.powers[1])
    @eval Base.$op(c::_PowerCachedSeries, d::_PowerCachedSeries) = $op(c.powers[1], d.powers[1])
end

# Taylor series of A and B about a real x0 from the closed-form GSN potentials:
# with g = dx/dr*, X'' + A X' + B X = 0 in x has A = (dg/dx - F) / g and
# B = -U / g^2, for X_{r*r*} - F X_{r*} - U X = 0. Its B is good to about 2e-13
# where the tables are off by 1e-12 to 2e-12 (x = 0.2 - 0.4 for s = -2,
# l = 3, a = 0.9, omega = 2.35), at about 0.1 ms per call.
function direct_ab_series!(a_out, b_out, params, x0::Real, order::Integer)
    s, m, a, omega, lambda = params.s, params.m, params.a, params.omega, params.lambda
    kappa = sqrt(1 - a^2)
    rplus, rminus = 1 + kappa, 1 - kappa
    x = x0 + Taylor1(ComplexF64, order + 1)
    r = (rplus - rminus * x) / (1 - x)
    delta = r^2 - 2r + a^2
    g = (rplus - rminus) * delta / ((r - rminus)^2 * (r^2 + a^2))
    cached = _PowerCachedSeries(r, 7)
    A = (derivative(g) - sF(s, m, a, omega, lambda, cached)) / g
    B = -sU(s, m, a, omega, lambda, cached) / g^2
    a_terms, b_terms = ab_term_counts(:O, order)
    for k in 1:a_terms
        a_out[k] = getcoeff(A, k - 1)
    end
    for k in 1:b_terms
        b_out[k] = getcoeff(B, k - 1)
    end
    return a_terms, b_terms
end

# A and B at a (complex) point x from the closed-form potentials, as in direct_ab_series!.
function direct_ab_value(params, x)
    s, m, a, omega, lambda = params.s, params.m, params.a, params.omega, params.lambda
    kappa = sqrt(1 - a^2)
    rplus, rminus = 1 + kappa, 1 - kappa
    r = (rplus - rminus * x) / (1 - x)
    drdx = (rplus - rminus) / (1 - x)^2
    delta = r^2 - 2r + a^2
    den = (r - rminus)^2 * (r^2 + a^2)
    g = (rplus - rminus) * delta / den
    dgdr = (rplus - rminus) * ((2r - 2) - delta * (2 / (r - rminus) + 2r / (r^2 + a^2))) / den
    return (dgdr * drdx - sF(s, m, a, omega, lambda, r)) / g, -sU(s, m, a, omega, lambda, r) / g^2
end

# Taylor series of A and B about the ordinary point x0: by shifting the horizon
# rational functions when `metadata.horizon_shift` is set or x0 is a real point
# below HORIZON_SHIFT_MAX_X, and the shift is accepted; from the closed-form
# potentials when DIRECT_AB_KEY is set and x0 is real and below DIRECT_AB_MAX_X;
# otherwise from the ordinary-point tables.
function ordinary_ab_series!(
    a_out::AbstractVector{<:Complex},
    b_out::AbstractVector{<:Complex},
    pq::AbstractVector{<:Complex},
    coefficients::CoefficientSet,
    x0,
    order::Integer,
)
    if x0 isa Real && HORIZON_SHIFT_MAX_X <= x0 < DIRECT_AB_MAX_X &&
            coefficients.params.omega isa Real &&
            get(task_local_storage(), DIRECT_AB_KEY, false)
        return direct_ab_series!(a_out, b_out, coefficients.params, x0, order)
    end
    use_shift = (hasproperty(coefficients.metadata, :horizon_shift) &&
        getproperty(coefficients.metadata, :horizon_shift)) ||
        (x0 isa Real && 0 < x0 < HORIZON_SHIFT_MAX_X)
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
    C = eltype(coefficients.ordinary.coefficients)
    pq = Vector{C}(undef, coefficients.ordinary.value_count)
    a_coeffs = Vector{C}(undef, a_terms)
    b_coeffs = Vector{C}(undef, b_terms)
    ordinary_ab_series!(a_coeffs, b_coeffs, pq, coefficients, x0, order)
    return a_coeffs, b_coeffs, pq
end

# Taylor coefficients of a solution about an ordinary point from the A and B
# series and the initial value c0 and slope c1.
function ordinary_point_coeffs!(
    coeffs::AbstractVector{<:Complex},
    a_coeffs,
    b_coeffs,
    c0,
    c1,
    order::Integer,
)
    order >= 0 || throw(ArgumentError("series order must be non-negative."))
    length(coeffs) >= order + 1 || throw(ArgumentError("solution coefficient buffer is too short."))
    C = eltype(coeffs)
    coeffs[1] = convert(C, c0)
    order >= 1 && (coeffs[2] = convert(C, c1))
    order < 2 && return coeffs
    @inbounds for n in 0:(order - 2)
        numerator = zero(C)
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
        1 <= idx <= length(coeffs) ? coeffs[idx] : zero(eltype(coeffs))
    end
end

# Horner evaluation of a polynomial, with first (and second) derivatives in the
# _pair and _triple versions.
function poly_value(coeffs, z)
    isempty(coeffs) && return zero(eltype(coeffs))
    value = coeffs[end]
    @inbounds for i in (length(coeffs) - 1):-1:1
        value = value * z + coeffs[i]
    end
    return value
end

function poly_pair(coeffs, z)
    n = length(coeffs) - 1
    C = eltype(coeffs)
    n < 0 && return zero(C), zero(C)
    n == 0 && return coeffs[1], zero(C)
    value = coeffs[n + 1]
    derivative = zero(C)
    @inbounds for i in n:-1:1
        derivative = derivative * z + value
        value = value * z + coeffs[i]
    end
    return value, derivative
end

function poly_triple(coeffs, z)
    n = length(coeffs) - 1
    C = eltype(coeffs)
    n < 0 && return zero(C), zero(C), zero(C)
    n == 0 && return coeffs[1], zero(C), zero(C)
    value = coeffs[n + 1]
    derivative = zero(C)
    second = zero(C)
    @inbounds for i in n:-1:1
        second = second * z + 2 * derivative
        derivative = derivative * z + value
        value = value * z + coeffs[i]
    end
    return value, derivative, second
end

end
