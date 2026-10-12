module DirectLFE

export LFEPlan, lfe_plan, lfe_ab!, lfe_ab_fast!
export lfe_carrier_ab!, lfe_log_series!
export lfe_carrier_phase, lfe_carrier_rate, lfe_table_stats

# Large-frequency (LFE) ordinary-point coefficients. A binary table stores each
# coefficient as a weighted sum of monomials in (lambda, m, omega, sin nu, cos nu);
# these are summed in double-double arithmetic to limit cancellation at large omega.
# DIRECT_GSN_LFE_FAST_RATIONAL and DIRECT_GSN_LFE_FUSED_DD_HORNER (read at load
# time, default "1") select the binary64 series division and the fused DD Horner step.
const TABLE_PATH = joinpath(
    @__DIR__, "Generated", "direct_gsn_lfe_tables_v1.bin")
const TABLE_MAGIC = UInt8[0x44, 0x47, 0x4e, 0x4c, 0x46, 0x45, 0x31, 0x00]
const LFE_FAST_RATIONAL = get(
    ENV, "DIRECT_GSN_LFE_FAST_RATIONAL", "1") == "1"
const LFE_FUSED_DD_HORNER = get(
    ENV, "DIRECT_GSN_LFE_FUSED_DD_HORNER", "1") == "1"

# Double-double numbers: value = high + low with |low| <= ulp(high) / 2.
struct DDReal
    high::Float64
    low::Float64
end

DDReal(value::Real) = DDReal(Float64(value), 0.0)

struct DDComplex
    high::ComplexF64
    low::ComplexF64
end

DDComplex(value::Number) = DDComplex(ComplexF64(value), 0.0 + 0.0im)

struct LFETable
    s::Int
    basis::Symbol
    exponents::Vector{NTuple{5,UInt8}}
    coefficient_offsets::Vector{Int32}
    monomial_indices::Vector{UInt32}
    weights::Vector{ComplexF64}
    polynomial_offsets::Vector{Int32}
    polynomial_indices::Vector{Int32}
    coefficient_count::Int
    output_count::Int
end

mutable struct LFEPlan
    table::LFETable
    coefficients::Vector{DDComplex}
    pq::Vector{DDComplex}
    fast_pq::Vector{ComplexF64}
    avec::Vector{DDComplex}
    bvec::Vector{DDComplex}
    carrier::Vector{DDComplex}
    transformed_a::Vector{DDComplex}
    transformed_b::Vector{DDComplex}
    wvec::Vector{DDComplex}
    ellvec::Vector{DDComplex}
end

# Error-free sum (Knuth TwoSum): value + error == left + right exactly.
@inline function two_sum(left::Float64, right::Float64)
    value = left + right
    right_virtual = value - left
    error = (left - (value - right_virtual)) + (right - right_virtual)
    return value, error
end

@inline function dd_normalize(high::Float64, low::Float64)
    value, error = two_sum(high, low)
    return DDReal(value, error)
end

@inline function dd_add(left::DDReal, right::DDReal)
    high, error = two_sum(left.high, right.high)
    return dd_normalize(high, left.low + right.low + error)
end

@inline function dd_mul(left::DDReal, right::DDReal)
    high = left.high * right.high
    error = fma(left.high, right.high, -high)
    low = error + left.high * right.low + left.low * right.high
    return dd_normalize(high, low)
end

@inline function dc_normalize(high::ComplexF64, low::ComplexF64)
    real_high, real_low = two_sum(real(high), real(low))
    imag_high, imag_low = two_sum(imag(high), imag(low))
    return DDComplex(
        ComplexF64(real_high, imag_high),
        ComplexF64(real_low, imag_low),
    )
end

@inline function dc_add(left::DDComplex, right::DDComplex)
    real_high, real_error = two_sum(real(left.high), real(right.high))
    imag_high, imag_error = two_sum(imag(left.high), imag(right.high))
    low = left.low + right.low + ComplexF64(real_error, imag_error)
    return dc_normalize(ComplexF64(real_high, imag_high), low)
end

@inline dc_neg(value::DDComplex) = DDComplex(-value.high, -value.low)
@inline dc_sub(left::DDComplex, right::DDComplex) =
    dc_add(left, dc_neg(right))
@inline dc_imul(value::DDComplex) = DDComplex(
    ComplexF64(-imag(value.high), real(value.high)),
    ComplexF64(-imag(value.low), real(value.low)),
)

# Complex product with its rounding error, using fma for the partial products.
@inline function complex_product(left::ComplexF64, right::ComplexF64)
    rr = real(left) * real(right)
    ii = imag(left) * imag(right)
    ri = real(left) * imag(right)
    ir = imag(left) * real(right)
    real_high, real_sum_error = two_sum(rr, -ii)
    imag_high, imag_sum_error = two_sum(ri, ir)
    real_error = fma(real(left), real(right), -rr) -
        fma(imag(left), imag(right), -ii) + real_sum_error
    imag_error = fma(real(left), imag(right), -ri) +
        fma(imag(left), real(right), -ir) + imag_sum_error
    return ComplexF64(real_high, imag_high),
        ComplexF64(real_error, imag_error)
end

@inline function dc_mul(left::DDComplex, right::DDComplex)
    high, error = complex_product(left.high, right.high)
    low = error + left.high * right.low + left.low * right.high
    return dc_normalize(high, low)
end

@inline function dc_scale(value::DDComplex, factor::Float64)
    high_re = real(value.high) * factor
    high_im = imag(value.high) * factor
    error = ComplexF64(
        fma(real(value.high), factor, -high_re),
        fma(imag(value.high), factor, -high_im),
    )
    return dc_normalize(
        ComplexF64(high_re, high_im), error + value.low * factor)
end

@inline function dc_madd(
        value::DDComplex, factor::Float64, addend::DDComplex)
    product_re = real(value.high) * factor
    product_im = imag(value.high) * factor
    product_error = ComplexF64(
        fma(real(value.high), factor, -product_re),
        fma(imag(value.high), factor, -product_im),
    )
    high_re, sum_error_re = two_sum(product_re, real(addend.high))
    high_im, sum_error_im = two_sum(product_im, imag(addend.high))
    low = product_error + value.low * factor + addend.low +
        ComplexF64(sum_error_re, sum_error_im)
    return dc_normalize(ComplexF64(high_re, high_im), low)
end

@inline function dc_weight(value::DDReal, weight::ComplexF64)
    high_re = value.high * real(weight)
    high_im = value.high * imag(weight)
    error = ComplexF64(
        fma(value.high, real(weight), -high_re),
        fma(value.high, imag(weight), -high_im),
    )
    low = error + value.low * weight
    return dc_normalize(ComplexF64(high_re, high_im), low)
end

# Double-double division: binary64 quotient refined by two residual corrections.
@inline function dc_div(numerator::DDComplex, denominator::DDComplex)
    quotient = DDComplex(numerator.high / denominator.high)
    for _ in 1:2
        residual = dc_sub(numerator, dc_mul(denominator, quotient))
        quotient = dc_add(quotient,
            DDComplex(residual.high / denominator.high))
    end
    return quotient
end

@inline dc_value(value::DDComplex) = value.high + value.low

function read_vector(io, ::Type{T}, count::Integer) where {T}
    values = Vector{T}(undef, count)
    read!(io, values)
    return values
end

# One (s, basis) table: header, monomial exponents, per-coefficient term lists
# and weights, and per-output Horner coefficient lists.
function read_table(io)
    s = Int(read(io, Int8))
    basis_code = read(io, UInt8)
    read(io, UInt16)
    counts = read_vector(io, UInt32, 5)
    coefficient_count, output_count, monomial_count, term_count,
        polynomial_count = Int.(counts)
    exponents = Vector{NTuple{5,UInt8}}(undef, monomial_count)
    for index in eachindex(exponents)
        values = read_vector(io, UInt8, 5)
        exponents[index] = Tuple(values)
    end
    coefficient_offsets = read_vector(
        io, Int32, coefficient_count + 1)
    monomial_indices = read_vector(io, UInt32, term_count)
    weights = Vector{ComplexF64}(undef, term_count)
    for index in eachindex(weights)
        weights[index] = ComplexF64(read(io, Float64), read(io, Float64))
    end
    polynomial_offsets = read_vector(io, Int32, output_count + 1)
    polynomial_indices = read_vector(io, Int32, polynomial_count)
    basis = basis_code == 1 ? :yc :
        basis_code == 2 ? :ys : error("invalid LFE basis code $basis_code")
    return LFETable(
        s, basis, exponents, coefficient_offsets, monomial_indices,
        weights, polynomial_offsets, polynomial_indices,
        coefficient_count, output_count)
end

function load_tables(path)
    isfile(path) || error("Direct GSN LFE table file is missing: $path")
    tables = Dict{Tuple{Int,Symbol},LFETable}()
    open(path, "r") do io
        read_vector(io, UInt8, length(TABLE_MAGIC)) == TABLE_MAGIC ||
            error("invalid Direct GSN LFE table header")
        version = read(io, UInt32)
        version == 1 || error("unsupported Direct GSN LFE table version $version")
        count = Int(read(io, UInt32))
        for _ in 1:count
            table = read_table(io)
            tables[(table.s, table.basis)] = table
        end
        eof(io) || error("trailing bytes in Direct GSN LFE table file")
    end
    return tables
end

const TABLES = load_tables(TABLE_PATH)

function lfe_table_stats()
    return map(sort!(collect(keys(TABLES)))) do key
        table = TABLES[key]
        (
            s=table.s,
            basis=table.basis,
            coefficients=table.coefficient_count,
            monomials=length(table.exponents),
            terms=length(table.weights),
            outputs=table.output_count,
        )
    end
end

function real_powers(value, maximum_power)
    powers = Vector{DDReal}(undef, maximum_power + 1)
    powers[1] = DDReal(1.0)
    maximum_power == 0 && return powers
    powers[2] = DDReal(value)
    @inbounds for power in 2:maximum_power
        powers[power + 1] = dd_mul(powers[power], powers[2])
    end
    return powers
end

# Evaluate the table coefficients for `params` in double-double and allocate
# the series work buffers for `order`.
function lfe_plan(params, basis::Symbol, order::Integer)
    table = get(TABLES, (Int(params.s), basis), nothing)
    table === nothing && error(
        "no Direct GSN LFE coefficient table for s=$(params.s), basis=$basis")
    sine = sin(Float64(params.nu))
    cosine = basis == :yc ?
        sqrt((1.0 - sine) * (1.0 + sine)) :
        sqrt(max(0.0, 1.0 - sine * sine))
    bases = (
        Float64(params.lambda),
        Float64(params.m),
        Float64(params.omega),
        sine,
        cosine,
    )
    maxima = ntuple(variable ->
        maximum(exponents[variable] for exponents in table.exponents), 5)
    powers = ntuple(variable ->
        real_powers(bases[variable], maxima[variable]), 5)
    monomials = Vector{DDReal}(undef, length(table.exponents))
    @inbounds for index in eachindex(monomials)
        value = DDReal(1.0)
        exponents = table.exponents[index]
        for variable in 1:5
            exponent = Int(exponents[variable])
            iszero(exponent) ||
                (value = dd_mul(value, powers[variable][exponent + 1]))
        end
        monomials[index] = value
    end
    coefficients = Vector{DDComplex}(undef, table.coefficient_count)
    @inbounds for coefficient_index in eachindex(coefficients)
        value = DDComplex(0)
        first_term = Int(table.coefficient_offsets[coefficient_index])
        last_term = Int(table.coefficient_offsets[coefficient_index + 1]) - 1
        for term_index in first_term:last_term
            monomial = monomials[Int(table.monomial_indices[term_index])]
            value = dc_add(value,
                dc_weight(monomial, table.weights[term_index]))
        end
        coefficients[coefficient_index] = value
    end
    order_int = Int(order)
    return LFEPlan(
        table,
        coefficients,
        Vector{DDComplex}(undef, table.output_count),
        Vector{ComplexF64}(undef, table.output_count),
        Vector{DDComplex}(undef, order_int + 1),
        Vector{DDComplex}(undef, order_int + 1),
        Vector{DDComplex}(undef, order_int + 2),
        Vector{DDComplex}(undef, order_int + 1),
        Vector{DDComplex}(undef, order_int + 1),
        Vector{DDComplex}(undef, order_int + 1),
        Vector{DDComplex}(undef, order_int + 2),
    )
end

# Packed P/Q values at x0, by Horner evaluation in y0 = 1 - x0.
function evaluate_pq!(plan::LFEPlan, x0)
    table = plan.table
    factor = 1.0 - Float64(x0)
    @inbounds for output_index in 1:table.output_count
        first_index = Int(table.polynomial_offsets[output_index])
        last_index = Int(table.polynomial_offsets[output_index + 1]) - 1
        if first_index > last_index
            plan.pq[output_index] = DDComplex(0)
            continue
        end
        coefficient_index = Int(table.polynomial_indices[first_index])
        value = plan.coefficients[coefficient_index]
        for position in (first_index + 1):last_index
            coefficient_index = Int(table.polynomial_indices[position])
            @static if LFE_FUSED_DD_HORNER
                value = dc_madd(value, factor,
                    plan.coefficients[coefficient_index])
            else
                value = dc_add(dc_scale(value, factor),
                    plan.coefficients[coefficient_index])
            end
        end
        plan.pq[output_index] = value
    end
    return plan.pq
end

# Power series of numerator / denominator by recursive division (double-double
# here, binary64 in rational_fast!).
function rational_series!(out, numerator, denominator, order)
    q0 = denominator[1]
    qinv = dc_div(DDComplex(1), q0)
    @inbounds for n in 0:order
        value = n < length(numerator) ? numerator[n + 1] : DDComplex(0)
        for k in 1:min(n, length(denominator) - 1)
            value = dc_sub(value,
                dc_mul(denominator[k + 1], out[n - k + 1]))
        end
        out[n + 1] = dc_mul(value, qinv)
    end
    return out
end

function rational_fast!(out, values, numerator_start, numerator_length,
        denominator_start, denominator_length, order)
    q0 = values[denominator_start]
    @inbounds for n in 0:order
        value = n < numerator_length ? values[numerator_start + n] :
            ComplexF64(0)
        for k in 1:min(n, denominator_length - 1)
            value -= values[denominator_start + k] * out[n - k + 1]
        end
        out[n + 1] = value / q0
    end
    return out
end

function lfe_layout(s::Integer)
    s == 0 && return (10, 11, 9, 13)
    s in (-2, -1, 1, 2) && return (10, 11, 15, 19)
    error("unsupported Direct GSN LFE spin $s")
end

# Taylor coefficients of A and B about x0 up to `order`.
function lfe_ab!(a_out, b_out, plan::LFEPlan, x0, order::Integer)
    LFE_FAST_RATIONAL && return lfe_ab_fast!(
        a_out, b_out, plan, x0, order)
    order_int = Int(order)
    order_int + 1 <= length(plan.avec) ||
        throw(ArgumentError("LFE plan order is too short"))
    order_int + 1 <= length(a_out) ||
        throw(ArgumentError("A output is too short"))
    order_int + 1 <= length(b_out) ||
        throw(ArgumentError("B output is too short"))
    pq = evaluate_pq!(plan, x0)
    ap_len, aq_len, bp_len, bq_len = lfe_layout(plan.table.s)
    aq_start = ap_len + 1
    bp_start = aq_start + aq_len
    bq_start = bp_start + bp_len
    rational_series!(plan.avec,
        @view(pq[1:ap_len]), @view(pq[aq_start:(aq_start + aq_len - 1)]),
        order_int)
    rational_series!(plan.bvec,
        @view(pq[bp_start:(bp_start + bp_len - 1)]),
        @view(pq[bq_start:(bq_start + bq_len - 1)]), order_int)
    @inbounds for index in 1:(order_int + 1)
        a_out[index] = dc_value(plan.avec[index])
        b_out[index] = dc_value(plan.bvec[index])
    end
    return order_int + 1, order_int + 1
end

function lfe_ab_fast!(a_out, b_out, plan::LFEPlan, x0, order::Integer)
    order_int = Int(order)
    order_int + 1 <= length(a_out) ||
        throw(ArgumentError("A output is too short"))
    order_int + 1 <= length(b_out) ||
        throw(ArgumentError("B output is too short"))
    pq = evaluate_pq!(plan, x0)
    @inbounds for index in eachindex(plan.fast_pq)
        plan.fast_pq[index] = dc_value(pq[index])
    end
    ap_len, aq_len, bp_len, bq_len = lfe_layout(plan.table.s)
    aq_start = ap_len + 1
    bp_start = aq_start + aq_len
    bq_start = bp_start + bp_len
    rational_fast!(a_out, plan.fast_pq,
        1, ap_len, aq_start, aq_len, order_int)
    rational_fast!(b_out, plan.fast_pq,
        bp_start, bp_len, bq_start, bq_len, order_int)
    return order_int + 1, order_int + 1
end

@inline function _log_ratio(numerator, denominator)
    delta = (numerator - denominator) / denominator
    return abs(delta) <= 0.5 ? log1p(delta) :
        log(numerator) - log(denominator)
end

# Integral of lfe_carrier_rate from start_x to destination_x, in closed form.
function lfe_carrier_phase(params, sign, start_x, destination_x)
    kappa = params.kappa
    rplus = 1 + kappa
    horizon = (2 * params.omega * rplus - params.a * params.m) /
        (2 * kappa)
    delta_r = 2 * kappa * (destination_x - start_x) /
        ((1 - destination_x) * (1 - start_x))
    return sign * (
        params.omega * delta_r +
        horizon * _log_ratio(destination_x, start_x) -
        2 * params.omega * _log_ratio(
            1 - destination_x, 1 - start_x)
    )
end

# WKB phase rate in x: sign * (2 kappa omega / y^2 + 2 omega / y + nu_h / x),
# with y = 1 - x and nu_h the horizon exponent scale.
function lfe_carrier_rate(params, sign, x)
    y = 1 - x
    horizon = (2 * params.omega * (1 + params.kappa) -
        params.a * params.m) / (2 * params.kappa)
    return sign * (
        2 * params.kappa * params.omega / y^2 +
        2 * params.omega / y + horizon / x)
end

# Taylor coefficients of lfe_carrier_rate about x0 (one order higher than A, B).
function carrier_series!(plan::LFEPlan, params, x0, order, sign)
    order + 2 <= length(plan.carrier) ||
        throw(ArgumentError("LFE carrier buffer is too short"))
    inverse_x = DDReal(inv(Float64(x0)))
    inverse_y = DDReal(inv(1.0 - Float64(x0)))
    xpower = inverse_x
    ypower = inverse_y
    first_scale = 2 * Float64(params.kappa) * Float64(params.omega)
    second_scale = 2 * Float64(params.omega)
    horizon = (2 * Float64(params.omega) *
        (1 + Float64(params.kappa)) -
        Float64(params.a) * Float64(params.m)) /
        (2 * Float64(params.kappa))
    @inbounds for n in 0:(order + 1)
        ypower2 = dd_mul(ypower, inverse_y)
        first = dd_mul(DDReal(first_scale * (n + 1)), ypower2)
        second = dd_mul(DDReal(second_scale), ypower)
        alternating = isodd(n) ? -horizon : horizon
        third = dd_mul(DDReal(alternating), xpower)
        value = dd_add(dd_add(first, second), third)
        plan.carrier[n + 1] = DDComplex(
            ComplexF64(sign * value.high),
            ComplexF64(sign * value.low),
        )
        xpower = dd_mul(xpower, inverse_x)
        ypower = dd_mul(ypower, inverse_y)
    end
    return plan.carrier
end

# A and B series for the envelope W in X = exp(i int c dx) W, with c the carrier
# rate: A + 2 i c and B + i c' - c^2 + i A c.
function carrier_coefficients!(plan::LFEPlan, params, x0, order, sign)
    order_int = Int(order)
    order_int + 1 <= length(plan.avec) ||
        throw(ArgumentError("LFE plan order is too short"))
    pq = evaluate_pq!(plan, x0)
    ap_len, aq_len, bp_len, bq_len = lfe_layout(plan.table.s)
    aq_start = ap_len + 1
    bp_start = aq_start + aq_len
    bq_start = bp_start + bp_len
    rational_series!(plan.avec,
        @view(pq[1:ap_len]), @view(pq[aq_start:(aq_start + aq_len - 1)]),
        order_int)
    rational_series!(plan.bvec,
        @view(pq[bp_start:(bp_start + bp_len - 1)]),
        @view(pq[bq_start:(bq_start + bq_len - 1)]), order_int)
    carrier = carrier_series!(plan, params, x0, order_int, sign)
    @inbounds for n in 0:order_int
        square = DDComplex(0)
        product = DDComplex(0)
        for k in 0:n
            square = dc_add(square,
                dc_mul(carrier[k + 1], carrier[n - k + 1]))
            product = dc_add(product,
                dc_mul(plan.avec[k + 1], carrier[n - k + 1]))
        end
        derivative = dc_scale(carrier[n + 2], Float64(n + 1))
        transformed_b = dc_add(plan.bvec[n + 1], dc_imul(derivative))
        transformed_b = dc_sub(transformed_b, square)
        transformed_b = dc_add(transformed_b, dc_imul(product))
        plan.transformed_b[n + 1] = transformed_b
    end
    @inbounds for n in 0:order_int
        plan.transformed_a[n + 1] = dc_add(plan.avec[n + 1],
            dc_scale(dc_imul(carrier[n + 1]), 2.0))
    end
    return plan.transformed_a, plan.transformed_b
end

function lfe_carrier_ab!(a_out, b_out, plan::LFEPlan, params,
        x0, order::Integer, sign)
    order_int = Int(order)
    order_int + 1 <= length(a_out) ||
        throw(ArgumentError("A output is too short"))
    order_int + 1 <= length(b_out) ||
        throw(ArgumentError("B output is too short"))
    transformed_a, transformed_b = carrier_coefficients!(
        plan, params, x0, order_int, sign)
    @inbounds for index in 1:(order_int + 1)
        a_out[index] = dc_value(transformed_a[index])
        b_out[index] = dc_value(transformed_b[index])
    end
    return order_int + 1, order_int + 1
end

# Series of the Riccati variable w = W'/W of the carrier envelope (from w0) and of
# its integral ell = log W - log W(x0).
function lfe_log_series!(w_out, ell_out, plan::LFEPlan, params,
        x0, order::Integer, sign, w0)
    order_int = Int(order)
    order_int + 1 <= length(w_out) ||
        throw(ArgumentError("log-derivative output is too short"))
    order_int + 2 <= length(ell_out) ||
        throw(ArgumentError("log-envelope output is too short"))
    transformed_a, transformed_b = carrier_coefficients!(
        plan, params, x0, order_int, sign)
    wcoeffs = plan.wvec
    ellcoeffs = plan.ellvec
    wcoeffs[1] = DDComplex(w0)
    @inbounds for n in 0:(order_int - 1)
        rhs = transformed_b[n + 1]
        for k in 0:n
            rhs = dc_add(rhs,
                dc_mul(wcoeffs[k + 1], wcoeffs[n - k + 1]))
            rhs = dc_add(rhs,
                dc_mul(transformed_a[k + 1], wcoeffs[n - k + 1]))
        end
        wcoeffs[n + 2] = dc_scale(rhs, -1.0 / Float64(n + 1))
    end
    ellcoeffs[1] = DDComplex(0)
    @inbounds for n in 0:order_int
        ellcoeffs[n + 2] = dc_scale(
            wcoeffs[n + 1], 1.0 / Float64(n + 1))
        w_out[n + 1] = dc_value(wcoeffs[n + 1])
        ell_out[n + 2] = dc_value(ellcoeffs[n + 2])
    end
    ell_out[1] = 0
    return order_int + 1, order_int + 2
end

end
