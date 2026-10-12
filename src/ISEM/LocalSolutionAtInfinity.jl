module LocalSolutionAtInfinity

using ..CoefficientTables: CoefficientSet
using ..OrdinaryPointExpansion:
    endpoint_ab_series,
    endpoint_power_lookup,
    denom_ok,
    poly_pair,
    poly_triple,
    poly_value,
    series_issue

export infinity_exponents
export infinity_coefficients!, infinity_coefficients
export infinity_guarded
export infinity_local_solution, inf_in_coeffs, inf_out_coeffs
export infinity_series_value, infinity_series_pair, infinity_series_triple

_param_kappa(params) = hasproperty(params, :kappa) ? getproperty(params, :kappa) : cos(getproperty(params, :nu))

# (rho, sigma) of the asymptotic form exp(sigma / y) y^rho at y = 1 - x = 0.
function infinity_exponents(params, branch::Symbol)
    omega = getproperty(params, :omega)
    kappa = _param_kappa(params)
    if branch in (:out, :outgoing, :up, :UP)
        return -2im * omega, 2im * kappa * omega
    elseif branch in (:in, :ingoing, :down, :DOWN)
        return 2im * omega, -2im * kappa * omega
    end
    throw(ArgumentError("infinity branch must be :out/:up or :in/:down."))
end

# Irregular-singular-point recurrence at y = 0 for the series multiplying
# exp(sigma / y) y^rho. With `guarded`, stops at the first ill-conditioned
# denominator or suspicious coefficient and returns (effective order, reason).
function _infinity_coefficients!(
    coeffs::Vector{ComplexF64},
    a_power,
    b_power,
    rho,
    sigma,
    order::Integer,
    guarded::Bool,
)
    order >= 0 || throw(ArgumentError("series order must be non-negative."))
    length(coeffs) >= order + 1 || throw(ArgumentError("coefficient buffer is too short."))
    rho = ComplexF64(rho)
    sigma = ComplexF64(sigma)
    p(j) = j == 0 ? -2 * sigma : (j == 1 ? 2 * rho - a_power(-1) : -a_power(j - 2))
    q(j) = b_power(j - 2) + sigma * a_power(j) - rho * a_power(j - 1) +
        (j == 0 ? rho^2 - rho : 0.0 + 0.0im)
    coeffs[1] = 1.0 + 0.0im
    p0 = p(0)
    if !isfinite(real(p0)) || !isfinite(imag(p0)) || abs(p0) <= floatmin(Float64)
        guarded || error("invalid p0 in ISEM GSN infinity recurrence")
        return 0, :small_recurrence_denominator
    end
    @inbounds for n in 0:(order - 1)
        numerator = 0.0 + 0.0im
        if n >= 1
            numerator += (n - 1) * n * coeffs[n + 1]
        end
        for j in 1:(n + 1)
            numerator += (n + 1 - j) * p(j) * coeffs[n - j + 2]
        end
        for j in 0:n
            numerator += q(j) * coeffs[n - j + 1]
        end
        denominator = (n + 1) * p0
        if !denom_ok(denominator, abs(n + 1) * abs(p0))
            guarded || error("ill-conditioned denominator in ISEM GSN infinity recurrence at n=$n")
            return n, :small_recurrence_denominator
        end
        coeffs[n + 2] = -numerator / denominator
        if guarded
            issue = series_issue(coeffs, n + 2)
            issue == :none || return n, issue
        end
    end
    return order, :none
end

function infinity_coefficients!(
    coeffs::Vector{ComplexF64},
    a_power,
    b_power,
    rho,
    sigma,
    order::Integer,
)
    _infinity_coefficients!(coeffs, a_power, b_power, rho, sigma, order, false)
    return coeffs
end

function infinity_coefficients(a_power, b_power, rho, sigma, order::Integer)
    coeffs = Vector{ComplexF64}(undef, order + 1)
    return infinity_coefficients!(coeffs, a_power, b_power, rho, sigma, order)
end

function infinity_coefficients(
    a_coeffs::Vector{ComplexF64},
    b_coeffs::Vector{ComplexF64},
    rho,
    sigma,
    order::Integer;
    a_first_power::Integer=-1,
    b_first_power::Integer=-4,
)
    a_power = endpoint_power_lookup(a_coeffs, a_first_power)
    b_power = endpoint_power_lookup(b_coeffs, b_first_power)
    return infinity_coefficients(a_power, b_power, rho, sigma, order)
end

function infinity_guarded(
    a_coeffs::Vector{ComplexF64},
    b_coeffs::Vector{ComplexF64},
    rho,
    sigma,
    order::Integer;
    a_first_power::Integer=-1,
    b_first_power::Integer=-4,
)
    a_power = endpoint_power_lookup(a_coeffs, a_first_power)
    b_power = endpoint_power_lookup(b_coeffs, b_first_power)
    coeffs = Vector{ComplexF64}(undef, order + 1)
    effective_order, reason = _infinity_coefficients!(
        coeffs,
        a_power,
        b_power,
        rho,
        sigma,
        order,
        true,
    )
    minimum_order = min(2, Int(order))
    effective_order >= minimum_order || error(
        "ISEM GSN infinity recurrence lost Float64 credibility before order $minimum_order: $reason")
    effective_order < order && resize!(coeffs, effective_order + 1)
    return (
        coefficients=coeffs,
        requested_order=Int(order),
        effective_order=effective_order,
        truncation_reason=reason,
    )
end

# Laurent coefficients of A and B rewritten for the scaled variable y / variable_scale.
function _scale_infinity_ab(a_coeffs, b_coeffs, variable_scale)
    scale = ComplexF64(variable_scale)
    iszero(scale) && throw(DomainError(
        variable_scale, "scaled infinity variable requires a nonzero scale."))
    a_scaled = similar(a_coeffs)
    b_scaled = similar(b_coeffs)
    @inbounds for i in eachindex(a_coeffs)
        power = i - 2
        a_scaled[i] = a_coeffs[i] * scale^(power + 1)
    end
    @inbounds for i in eachindex(b_coeffs)
        power = i - 5
        b_scaled[i] = b_coeffs[i] * scale^(power + 2)
    end
    return a_scaled, b_scaled
end

# Guarded asymptotic series at infinity for the :out/:up or :in/:down branch,
# optionally in a scaled variable.
function infinity_local_solution(
    coefficients::CoefficientSet,
    branch::Symbol,
    order::Integer,
    ;
    variable_scale=nothing,
)
    a_coeffs, b_coeffs = endpoint_ab_series(coefficients, :I, order)
    rho, sigma = infinity_exponents(coefficients.params, branch)
    scale = variable_scale === nothing ?
        1.0 + 0.0im : ComplexF64(variable_scale)
    recurrence_a, recurrence_b, recurrence_sigma = if isone(scale)
        a_coeffs, b_coeffs, sigma
    else
        scaled_a, scaled_b = _scale_infinity_ab(
            a_coeffs, b_coeffs, scale)
        scaled_a, scaled_b, sigma / scale
    end
    guarded = infinity_guarded(
        recurrence_a, recurrence_b, rho, recurrence_sigma, order)
    return (
        branch=branch,
        rho=rho,
        sigma=sigma,
        recurrence_sigma=recurrence_sigma,
        coefficients=guarded.coefficients,
        A=recurrence_a,
        B=recurrence_b,
        variable=isone(scale) ? :y : :scaled_y,
        variable_scale=scale,
        requested_order=guarded.requested_order,
        effective_order=guarded.effective_order,
        truncation_reason=guarded.truncation_reason,
    )
end

inf_out_coeffs(coefficients::CoefficientSet, order::Integer) =
    infinity_local_solution(coefficients, :out, order).coefficients

inf_in_coeffs(coefficients::CoefficientSet, order::Integer) =
    infinity_local_solution(coefficients, :in, order).coefficients

# exp(sigma / y) y^rho times the series. The _pair and _triple versions return
# derivatives with respect to x = 1 - y.
function infinity_series_value(coeffs, rho, sigma, y)
    return exp(sigma / y) * y^rho * poly_value(coeffs, y)
end

function infinity_series_pair(coeffs, rho, sigma, y)
    series, dseries = poly_pair(coeffs, y)
    prefactor = exp(sigma / y) * y^rho
    value = prefactor * series
    derivative_y = prefactor * (dseries + (-sigma / y^2 + rho / y) * series)
    return value, -derivative_y
end

function infinity_series_triple(coeffs, rho, sigma, y)
    series, dseries, ddseries = poly_triple(coeffs, y)
    prefactor = exp(sigma / y) * y^rho
    log_derivative = -sigma / y^2 + rho / y
    log_second = 2 * sigma / y^3 - rho / y^2
    value = prefactor * series
    derivative_y = prefactor * (dseries + log_derivative * series)
    second_y = prefactor * (
        ddseries + 2 * log_derivative * dseries +
        (log_derivative^2 + log_second) * series
    )
    return value, -derivative_y, second_y
end

end
