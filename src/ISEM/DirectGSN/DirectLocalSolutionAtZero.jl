module DirectLocalSolutionAtZero

using ..DirectCoefficientTables: DirectCoefficientSet
using ..DirectOrdinaryPointExpansion:
    direct_endpoint_ab_series,
    direct_endpoint_power_lookup,
    direct_denom_ok,
    direct_poly_pair,
    direct_poly_triple,
    direct_poly_value,
    direct_series_issue

export direct_horizon_frequency, direct_horizon_exponent
export direct_horizon_coefficients!, direct_horizon_coefficients
export direct_zero_local_solution, direct_zero_log_solution
export direct_zero_in_coeffs, direct_zero_out_coeffs
export direct_horizon_series_value, direct_horizon_series_pair, direct_horizon_series_triple

_param_a(params) = hasproperty(params, :a) ? getproperty(params, :a) : sin(getproperty(params, :nu))
_param_kappa(params) = hasproperty(params, :kappa) ? getproperty(params, :kappa) : cos(getproperty(params, :nu))

# Horizon exponent scale (2 r_+ omega - a m) / (2 kappa); the local solutions
# behave as x^(-i nu_h) (in) and x^(+i nu_h) (out) at x = 0.
function direct_horizon_frequency(params)
    a = _param_a(params)
    kappa = _param_kappa(params)
    r_plus = 1 + kappa
    return (2 * r_plus * getproperty(params, :omega) - a * getproperty(params, :m)) / (2 * kappa)
end

function direct_horizon_exponent(params, branch::Symbol)
    nu_h = direct_horizon_frequency(params)
    branch in (:in, :ingoing, :IN) && return -1im * nu_h
    branch in (:out, :outgoing, :OUT) && return 1im * nu_h
    throw(ArgumentError("horizon branch must be :in or :out."))
end

# Frobenius recurrence at x = 0 for X'' + A X' + B X = 0 with A ~ x^-1, B ~ x^-2
# and exponent mu. With `guarded`, stops at the first ill-conditioned denominator
# or suspicious coefficient and returns (effective order, reason).
function _horizon_coefficients!(
    coeffs::Vector{ComplexF64},
    a_power,
    b_power,
    mu,
    order::Integer,
    guarded::Bool,
)
    order >= 0 || throw(ArgumentError("series order must be non-negative."))
    length(coeffs) >= order + 1 || throw(ArgumentError("coefficient buffer is too short."))
    mu = ComplexF64(mu)
    coeffs[1] = 1.0 + 0.0im
    a_minus_one = a_power(-1)
    b_minus_two = b_power(-2)
    @inbounds for n in 1:order
        numerator = 0.0 + 0.0im
        for j in 0:(n - 1)
            numerator += a_power(j) * (mu + n - j - 1) * coeffs[n - j]
        end
        for j in -1:(n - 2)
            numerator += b_power(j) * coeffs[n - j - 1]
        end
        denominator = (mu + n) * (mu + n - 1) + a_minus_one * (mu + n) + b_minus_two
        denominator_scale = abs((mu + n) * (mu + n - 1)) +
            abs(a_minus_one * (mu + n)) + abs(b_minus_two)
        if !direct_denom_ok(denominator, denominator_scale)
            guarded || error("ill-conditioned denominator in direct GSN horizon recurrence at n=$n")
            return n - 1, :small_recurrence_denominator
        end
        coeffs[n + 1] = -numerator / denominator
        if guarded
            issue = direct_series_issue(coeffs, n + 1)
            issue == :none || return n - 1, issue
        end
    end
    return order, :none
end

function direct_horizon_coefficients!(
    coeffs::Vector{ComplexF64},
    a_power,
    b_power,
    mu,
    order::Integer,
)
    _horizon_coefficients!(coeffs, a_power, b_power, mu, order, false)
    return coeffs
end

function direct_horizon_coefficients(a_power, b_power, mu, order::Integer)
    coeffs = Vector{ComplexF64}(undef, order + 1)
    return direct_horizon_coefficients!(coeffs, a_power, b_power, mu, order)
end

function direct_horizon_coefficients(
    a_coeffs::Vector{ComplexF64},
    b_coeffs::Vector{ComplexF64},
    mu,
    order::Integer;
    a_first_power::Integer=-1,
    b_first_power::Integer=-2,
)
    a_power = direct_endpoint_power_lookup(a_coeffs, a_first_power)
    b_power = direct_endpoint_power_lookup(b_coeffs, b_first_power)
    return direct_horizon_coefficients(a_power, b_power, mu, order)
end

# Guarded horizon Frobenius series for the :in or :out branch, with metadata.
function direct_zero_local_solution(
    coefficients::DirectCoefficientSet,
    branch::Symbol,
    order::Integer,
)
    a_coeffs, b_coeffs = direct_endpoint_ab_series(coefficients, :H, order)
    mu = direct_horizon_exponent(coefficients.params, branch)
    a_power = direct_endpoint_power_lookup(a_coeffs, -1)
    b_power = direct_endpoint_power_lookup(b_coeffs, -2)
    solution_coeffs = Vector{ComplexF64}(undef, order + 1)
    effective_order, reason = _horizon_coefficients!(
        solution_coeffs,
        a_power,
        b_power,
        mu,
        order,
        true,
    )
    minimum_order = min(2, Int(order))
    effective_order >= minimum_order || error(
        "direct GSN horizon recurrence lost Float64 credibility before order $minimum_order: $reason")
    effective_order < order && resize!(solution_coeffs, effective_order + 1)
    return (
        branch=branch,
        exponent=mu,
        coefficients=solution_coeffs,
        A=a_coeffs,
        B=b_coeffs,
        variable=:x,
        requested_order=Int(order),
        effective_order=effective_order,
        truncation_reason=reason,
    )
end

function _horizon_numerator(coeffs, a_power, b_power, mu, n::Int)
    numerator = zero(ComplexF64)
    @inbounds for j in 0:(n - 1)
        numerator +=
            a_power(j) * (mu + n - j - 1) * coeffs[n - j]
    end
    @inbounds for j in -1:(n - 2)
        numerator += b_power(j) * coeffs[n - j - 1]
    end
    return ComplexF64(numerator)
end

function _horizon_log_forcing(coeffs, a_power, mu, n::Int)
    forcing = (2 * (mu + n) - 1 + a_power(-1)) * coeffs[n + 1]
    @inbounds for j in 0:(n - 1)
        forcing += a_power(j) * coeffs[n - j]
    end
    return ComplexF64(forcing)
end

# At the threshold frequency (nu_h = 0) the two horizon exponents coincide; this
# returns the logarithmic companion X = log(x) X_in + x^mu sum c_n x^n.
function direct_zero_log_solution(
    coefficients::DirectCoefficientSet,
    branch::Symbol,
    order::Integer,
)
    nu_h = direct_horizon_frequency(coefficients.params)
    iszero(nu_h) ||
        throw(ArgumentError(
            "logarithmic horizon solution requires the exact threshold frequency."))
    normalized = branch in (:out, :outgoing, :OUT) ? :out :
        throw(ArgumentError(
            "the threshold logarithmic companion is the outgoing horizon branch."))
    base = direct_zero_local_solution(coefficients, :in, order)
    a_power = direct_endpoint_power_lookup(base.A, -1)
    b_power = direct_endpoint_power_lookup(base.B, -2)
    mu = ComplexF64(base.exponent)
    companion = Vector{ComplexF64}(undef, order + 1)
    companion[1] = zero(ComplexF64)
    a_minus_one = a_power(-1)
    b_minus_two = b_power(-2)
    @inbounds for n in 1:order
        numerator = _horizon_numerator(
            companion, a_power, b_power, mu, n)
        numerator += _horizon_log_forcing(
            base.coefficients, a_power, mu, n)
        denominator = (mu + n) * (mu + n - 1) +
            a_minus_one * (mu + n) + b_minus_two
        denominator_scale = abs((mu + n) * (mu + n - 1)) +
            abs(a_minus_one * (mu + n)) + abs(b_minus_two)
        direct_denom_ok(denominator, denominator_scale) ||
            error("threshold logarithmic horizon recurrence is ill-conditioned at n=$n")
        companion[n + 1] = -numerator / denominator
        direct_series_issue(companion, n + 1) == :none ||
            error("threshold logarithmic horizon recurrence lost Float64 credibility at n=$n")
    end
    return (
        branch=normalized,
        exponent=mu,
        coefficients=companion,
        A=base.A,
        B=base.B,
        variable=:x,
        requested_order=Int(order),
        effective_order=Int(order),
        truncation_reason=:threshold_log_companion,
        representation=:threshold_log_frobenius,
        log_exponent=mu,
        log_coefficient=one(ComplexF64),
        log_coefficients=base.coefficients,
    )
end

direct_zero_in_coeffs(coefficients::DirectCoefficientSet, order::Integer) =
    direct_zero_local_solution(coefficients, :in, order).coefficients

direct_zero_out_coeffs(coefficients::DirectCoefficientSet, order::Integer) =
    direct_zero_local_solution(coefficients, :out, order).coefficients

# x^mu times the power series, with first and second x-derivatives in the
# _pair and _triple versions.
function direct_horizon_series_value(coeffs, mu, x)
    return x^mu * direct_poly_value(coeffs, x)
end

function direct_horizon_series_pair(coeffs, mu, x)
    series, dseries = direct_poly_pair(coeffs, x)
    prefactor = x^mu
    value = prefactor * series
    derivative = prefactor * (dseries + mu * series / x)
    return value, derivative
end

function direct_horizon_series_triple(coeffs, mu, x)
    series, dseries, ddseries = direct_poly_triple(coeffs, x)
    prefactor = x^mu
    log_derivative = mu / x
    log_second = -mu / x^2
    value = prefactor * series
    derivative = prefactor * (dseries + log_derivative * series)
    second = prefactor * (
        ddseries + 2 * log_derivative * dseries +
        (log_derivative^2 + log_second) * series
    )
    return value, derivative, second
end

end
