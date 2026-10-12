module NearExtreme

using ..CoefficientTables:
    CoefficientSet,
    EndpointCoefficientSet,
    RationalCoefficients
using ..LargeFrequencyExpansion:
    DDComplex,
    DDReal,
    dc_add,
    dc_div,
    dc_mul,
    dc_scale,
    dc_sub,
    dc_value,
    dc_weight,
    dd_mul,
    real_powers
using ..OrdinaryPointExpansion:
    ab_term_counts,
    endpoint_ab_series!
using ..Parameters: select_spin_regime

export factor_values_complex_dd
export near_extreme_prepare, near_extreme_selected

# Near-extremal horizon coefficients. Per-spin binary tables store the horizon
# A and B rational functions with common powers of kappa = cos(nu) factored out of
# numerator and denominator, avoiding the cancellation of the generic tables as kappa -> 0.
# ISEM_NEAR_EXTREME_COMPENSATED (default "1", read at load time) enables this.
const TABLE_DIR = joinpath(@__DIR__, "Generated")
const TABLE_MAGIC = UInt8[0x44, 0x47, 0x4e, 0x46, 0x43, 0x54, 0x31, 0x00]
const CERTIFICATE_ORDER = 24
const CERTIFICATE_LIMITS = (0.001, 0.001, 3e-5, 0.001, 0.001)
const FULL_COMPENSATION_LIMITS = (0.01, 0.5, 3e-5, 0.5, 0.012)
const ORDINARY_COMPRESSION_MIN = 128.0
const SMALL_FREQUENCY_EXPANSION_CERTIFICATE_LIMIT = 0.0032
const ORDINARY_VALIDATED_SPINS = (-2, -1, 0, 1, 2)
const SMALL_FREQUENCY_EXPANSION_VALIDATED_SPINS = (2,)
const COMPENSATED_ENABLED = get(
    ENV, "ISEM_NEAR_EXTREME_COMPENSATED", "1") == "1"

struct FactorTable
    s::Int
    powers::NTuple{4,Int}
    exponents::Vector{NTuple{5,UInt8}}
    offsets::Vector{Int32}
    indices::Vector{UInt32}
    weights::Vector{ComplexF64}
    maxima::NTuple{5,Int}
    outputs::Int
end

function spin_tag(spin)
    spin == -2 && return "sminus2"
    spin == -1 && return "sminus1"
    spin == 0 && return "s0"
    spin == 1 && return "splus1"
    spin == 2 && return "splus2"
    throw(ArgumentError("unsupported spin $spin"))
end

factor_table_path(spin) = joinpath(TABLE_DIR,
    "gsn_endpoint_factor_$(spin_tag(spin))_v1.bin")

function read_vector(io, ::Type{T}, count::Integer) where {T}
    values = Vector{T}(undef, count)
    read!(io, values)
    return values
end

function load_factor_table(path)
    isfile(path) || error("ISEM endpoint factor table is missing: $path")
    return open(path, "r") do io
        read_vector(io, UInt8, length(TABLE_MAGIC)) == TABLE_MAGIC ||
            error("invalid ISEM endpoint factor table header")
        version = read(io, UInt32)
        version == 1 || error(
            "unsupported ISEM endpoint factor table version $version")
        s = Int(read(io, Int8))
        powers = Tuple(Int.(read_vector(io, UInt8, 4)))
        read(io, UInt8)
        output_count = Int(read(io, UInt32))
        monomial_count = Int(read(io, UInt32))
        term_count = Int(read(io, UInt32))
        exponents = Vector{NTuple{5,UInt8}}(undef, monomial_count)
        for index in eachindex(exponents)
            exponents[index] = Tuple(read_vector(io, UInt8, 5))
        end
        offsets = read_vector(io, Int32, output_count + 1)
        indices = read_vector(io, UInt32, term_count)
        weights = Vector{ComplexF64}(undef, term_count)
        for index in eachindex(weights)
            weights[index] = ComplexF64(
                read(io, Float64), read(io, Float64))
        end
        eof(io) || error("trailing bytes in ISEM endpoint factor table")
        maxima = ntuple(variable -> maximum(
            Int(item[variable]) for item in exponents), 5)
        FactorTable(s, powers, exponents, offsets, indices,
            weights, maxima, output_count)
    end
end

const FACTOR_TABLES = ntuple(index ->
    load_factor_table(factor_table_path(index - 3)), 5)

@inline factor_table(spin) = FACTOR_TABLES[spin + 3]

function float_powers(value, maximum_power)
    powers = Vector{Float64}(undef, maximum_power + 1)
    powers[1] = 1.0
    maximum_power == 0 && return powers
    powers[2] = Float64(value)
    @inbounds for power in 2:maximum_power
        powers[power + 1] = powers[power] * powers[2]
    end
    return powers
end

# Factored horizon coefficients for real parameters in binary64, with the
# kappa powers restored by scale_factor!.
function factor_values(params)
    table = factor_table(params.s)
    sine = sin(params.nu)
    cosine = sqrt((1.0 - sine) * (1.0 + sine))
    bases = (params.lambda, Float64(params.m), params.omega, sine, cosine)
    powers = ntuple(variable ->
        float_powers(bases[variable], table.maxima[variable]), 5)
    monomials = Vector{Float64}(undef, length(table.exponents))
    @inbounds for index in eachindex(monomials)
        value = 1.0
        exponents = table.exponents[index]
        for variable in 1:5
            exponent = Int(exponents[variable])
            iszero(exponent) ||
                (value *= powers[variable][exponent + 1])
        end
        monomials[index] = value
    end

    values = Vector{ComplexF64}(undef, table.outputs)
    @inbounds for slot in eachindex(values)
        real_value = 0.0
        imag_value = 0.0
        first_term = Int(table.offsets[slot])
        last_term = Int(table.offsets[slot + 1]) - 1
        for term in first_term:last_term
            monomial = monomials[Int(table.indices[term])]
            weight = table.weights[term]
            real_value = muladd(monomial, real(weight), real_value)
            imag_value = muladd(monomial, imag(weight), imag_value)
        end
        values[slot] = ComplexF64(real_value, imag_value)
    end
    scale_factor!(values, cosine, table.powers, params.s)
    return values
end

# As factor_values, with the monomial sums in double-double.
function factor_values_dd(params)
    table = factor_table(params.s)
    sine = sin(params.nu)
    cosine = sqrt((1.0 - sine) * (1.0 + sine))
    bases = (params.lambda, Float64(params.m), params.omega, sine, cosine)
    powers = ntuple(variable ->
        real_powers(bases[variable], table.maxima[variable]), 5)
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

    values = Vector{ComplexF64}(undef, table.outputs)
    @inbounds for slot in eachindex(values)
        value = DDComplex(0)
        first_term = Int(table.offsets[slot])
        last_term = Int(table.offsets[slot + 1]) - 1
        for term in first_term:last_term
            value = dc_add(value, dc_weight(
                monomials[Int(table.indices[term])], table.weights[term]))
        end
        values[slot] = dc_value(value)
    end
    scale_factor!(values, cosine, table.powers, params.s)
    return values
end

function complex_powers_dd(value::DDComplex, maximum_power)
    powers = Vector{DDComplex}(undef, maximum_power + 1)
    powers[1] = DDComplex(1)
    maximum_power == 0 && return powers
    powers[2] = value
    @inbounds for power in 2:maximum_power
        powers[power + 1] = dc_mul(powers[power], powers[2])
    end
    return powers
end

# Double-double square root by three Newton steps from `initial`.
function sqrt_dd(value::DDComplex, initial::Float64)
    root = DDComplex(initial)
    @inbounds for _ in 1:3
        root = dc_scale(dc_add(root, dc_div(value, root)), 0.5)
    end
    return root
end

function pow_dd(value::DDComplex, exponent::Int)
    iszero(exponent) && return DDComplex(1)
    exponent < 0 && return dc_div(
        DDComplex(1), pow_dd(value, -exponent))
    result = DDComplex(1)
    base = value
    power = exponent
    while power > 0
        isodd(power) && (result = dc_mul(result, base))
        power >>= 1
        power > 0 && (base = dc_mul(base, base))
    end
    return result
end

# Factored horizon coefficients in double-double for complex omega and lambda,
# returned unrounded as DDComplex values.
function factor_values_complex_dd(params)
    table = factor_table(params.s)
    one = DDComplex(1)
    sine = DDComplex(params.a)
    cosine = sqrt_dd(
        dc_mul(dc_sub(one, sine), dc_add(one, sine)),
        params.kappa,
    )
    bases = (
        DDComplex(params.lambda),
        DDComplex(Float64(params.m)),
        DDComplex(params.omega),
        sine,
        cosine,
    )
    powers = ntuple(variable ->
        complex_powers_dd(bases[variable], table.maxima[variable]), 5)
    monomials = Vector{DDComplex}(undef, length(table.exponents))
    @inbounds for index in eachindex(monomials)
        value = DDComplex(1)
        exponents = table.exponents[index]
        for variable in 1:5
            exponent = Int(exponents[variable])
            iszero(exponent) ||
                (value = dc_mul(
                    value, powers[variable][exponent + 1]))
        end
        monomials[index] = value
    end

    values = Vector{DDComplex}(undef, table.outputs)
    @inbounds for slot in eachindex(values)
        value = DDComplex(0)
        first_term = Int(table.offsets[slot])
        last_term = Int(table.offsets[slot + 1]) - 1
        for term in first_term:last_term
            value = dc_add(value, dc_mul(
                monomials[Int(table.indices[term])],
                DDComplex(table.weights[term]),
            ))
        end
        values[slot] = value
    end

    ap, aq, bp, bq = table.powers
    a_scale = pow_dd(cosine, ap - aq)
    b_scale = pow_dd(cosine, bp - bq)
    @inbounds for index in 1:11
        values[index] = dc_mul(values[index], a_scale)
    end
    b_end = params.s == 0 ? 33 : 39
    @inbounds for index in 23:b_end
        values[index] = dc_mul(values[index], b_scale)
    end
    return values
end

# Multiply the A and B numerator blocks by cos(nu)^(numerator power - denominator power).
function scale_factor!(values, cosine, powers, spin)
    ap, aq, bp, bq = powers
    a_scale = cosine^(ap - aq)
    b_scale = cosine^(bp - bq)
    a_scale == 1 || @views values[1:11] .*= a_scale
    b_end = spin == 0 ? 33 : 39
    b_scale == 1 || @views values[23:b_end] .*= b_scale
    return values
end

# Wrap packed factored values as a horizon (:hc) endpoint coefficient set.
function factor_horizon(values, spin)
    b_length = spin == 0 ? 11 : 17
    q_start = 23 + b_length
    expected = spin == 0 ? 46 : 58
    length(values) == expected || error(
        "invalid s=$spin horizon factor layout")
    A = RationalCoefficients(
        @view(values[1:11]), @view(values[12:22]))
    B = RationalCoefficients(
        @view(values[23:(q_start - 1)]), @view(values[q_start:expected]))
    return EndpointCoefficientSet(A, B, :hc)
end

@inline finite_complex(value) =
    isfinite(real(value)) && isfinite(imag(value))

function scaled_delta(left, right)
    scale = max(maximum(abs, left), maximum(abs, right),
        floatmin(Float64))
    return maximum(abs.(left .- right)) / scale
end

# Relative difference of two horizon coefficient sets, compared through their
# A and B series to order CERTIFICATE_ORDER.
function endpoint_delta(left, right)
    a_terms, b_terms = ab_term_counts(:H, CERTIFICATE_ORDER)
    left_a = Vector{ComplexF64}(undef, a_terms)
    left_b = Vector{ComplexF64}(undef, b_terms)
    right_a = similar(left_a)
    right_b = similar(left_b)
    endpoint_ab_series!(left_a, left_b, left,
        :H, CERTIFICATE_ORDER)
    endpoint_ab_series!(right_a, right_b, right,
        :H, CERTIFICATE_ORDER)
    return max(scaled_delta(left_a, right_a),
        scaled_delta(left_b, right_b))
end

function near_extreme_selected(coefficients::CoefficientSet)
    metadata = coefficients.metadata
    return hasproperty(metadata, :near_extreme_compensated) &&
        getproperty(metadata, :near_extreme_compensated)
end

# For validated spins near extremality, compare the generic horizon coefficients
# with the factored ones; if they differ by more than the certificate limit (the
# generic table has lost accuracy), substitute the factored double-double set.
# `near_extreme_compensated` in the metadata marks full compensation.
function near_extreme_prepare(coefficients::CoefficientSet, controls)
    COMPENSATED_ENABLED || return coefficients
    getproperty(controls, :large_frequency_expansion) && return coefficients
    params = coefficients.params
    use_small_frequency_expansion = getproperty(controls, :small_frequency_expansion)
    validated_spins = use_small_frequency_expansion ? SMALL_FREQUENCY_EXPANSION_VALIDATED_SPINS : ORDINARY_VALIDATED_SPINS
    params.s in validated_spins || return coefficients
    spin = select_spin_regime(params)
    if use_small_frequency_expansion
        spin.regime == :near_extremal || return coefficients
    else
        spin.regime == :exact_extremal && return coefficients
        max(spin.angular_compression, spin.phase_compression, spin.tau_abs) >=
            ORDINARY_COMPRESSION_MIN || return coefficients
    end

    factorized = factor_horizon(factor_values(params), params.s)
    certificate = endpoint_delta(coefficients.horizon, factorized)
    certificate_limit = use_small_frequency_expansion ? SMALL_FREQUENCY_EXPANSION_CERTIFICATE_LIMIT :
        CERTIFICATE_LIMITS[params.s + 3]
    isfinite(certificate) && certificate >= certificate_limit ||
        return coefficients

    values = factor_values_dd(params)
    all(finite_complex, values) || return coefficients
    horizon = factor_horizon(values, params.s)
    full_compensation = use_small_frequency_expansion || (
        spin.regime == :near_extremal &&
        certificate >= FULL_COMPENSATION_LIMITS[params.s + 3]
    )
    metadata = merge(coefficients.metadata, (
        near_extreme_compensated=full_compensation,
        near_extreme_horizon_compensated=true,
        near_extreme_certificate=certificate,
    ))
    return CoefficientSet(params, horizon, coefficients.infinity,
        coefficients.ordinary, metadata)
end

end
