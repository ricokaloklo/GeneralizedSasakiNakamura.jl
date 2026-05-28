module LocalSolutionAtZero

include("TeukolskyTransformation.jl")
using .TeukolskyTransformation

export an_zero, bn_zero, zero_in_coeffs, zero_out_coeffs,
       convergence_radius_an_zero, convergence_radius_bn_zero,
       convergence_radius_zero_in, convergence_radius_zero_out,
       evaluate_zero_in, evaluate_zero_out, error_zero_in, error_zero_out,
       evaluate_zero_in_rho, evaluate_zero_out_rho, error_zero_in_rho, error_zero_out_rho

function an_zero(s, epsilon, tau, kappa, lambda, N)
    a = zeros(ComplexF64, N + 1)
    a[1] = 1.0 + 0.0im # a_0 = 1

    if N == 0
        return a
    end

    for n = 1:N
        term1_num = n^2 - s^2 + epsilon^2 + 1im * epsilon * kappa + s * (2im * epsilon * kappa - 1) -
                    lambda + tau * (1im - tau) - 1im * n * (2 * tau + 2 * epsilon * kappa - 1im)

        term2_num = - 2 * epsilon * kappa * (1im - 1im * n + 1im * s + epsilon - tau)

        denom = n * (n - s - 1im * (tau + epsilon))

        if isapprox(abs(denom), 0.0, atol=eps(Float64))
            error("Denominator in a_n recurrence is zero for n=$n. This implies superradiance, which requires special handling.")
        end

        if n == 1
            # a_1 calculation (a_{-1} = 0)
            a[n+1] = (term1_num / denom) * a[n] # a[n+1] corresponds to a_n
        else
            # a_n calculation (n >= 2)
            a[n+1] = (term1_num / denom) * a[n] + (term2_num / denom) * a[n-1]
        end
    end

    return a
end

function bn_zero(s, epsilon, tau, kappa, lambda, N)
    b = zeros(ComplexF64, N + 1)
    b[1] = 1.0 + 0.0im # b_0 = 1

    if N == 0
        return b
    end

    for n = 1:N
        term1_num = n^2 - 1im * epsilon + 2im * s * (1im + epsilon) + n * (2 * s - 1 - 2im * epsilon * (kappa - 1)) +
                    epsilon * kappa * (1im + 2 * epsilon + 2 * tau) - lambda

        term2_num = 2im * kappa * epsilon * (n + 2im * epsilon - 1)

        denom = n * (n + s + 1im * (tau + epsilon))

        if isapprox(denom, 0.0 + 0.0im, atol=eps(Float64))
            error("Denominator in b_n recurrence is zero for n=$n. This implies -s - i(epsilon + tau) is an integer, which requires special handling (logarithmic term).")
        end

        if n == 1
            # b_1 calculation (b_{-1} = 0)
            b[n+1] = (term1_num / denom) * b[n] # b[n+1] corresponds to b_n
        else
            # b_n calculation (n >= 2)
            b[n+1] = (term1_num / denom) * b[n] + (term2_num / denom) * b[n-1]
        end
    end

    return b
end

function convergence_radius_an_zero(coeffs, tol = 1e-15)
    N = length(coeffs) - 1
    if N <= 0
        return 0.0 # No coefficients to determine radius
    end
    a_N = coeffs[N+1] # coeffs is 1-indexed, so a_N is at index N+1
    if abs(a_N) == 0.0
        return Inf # If the last coefficient is zero, it implies infinite radius or a problem
    end
    return (tol / abs(a_N))^(1.0 / N)
end

function convergence_radius_bn_zero(coeffs, tol = 1e-15)
    N = length(coeffs) - 1
    if N <= 0
        return 0.0
    end
    b_N = coeffs[N+1]
    if abs(b_N) == 0.0
        return 0.6
    end
    return min((tol / abs(b_N))^(1.0 / N), 0.6)
end

function zero_in_coeffs(s, epsilon, tau, kappa, lambda, z, N::Integer; tol = 1e-15)
    coeffs = Vector{ComplexF64}(undef, N + 1)
    coeffs[1] = 1.0 + 0.0im
    N == 0 && return coeffs, 0.95

    for n in 1:N
        prev1 = -epsilon^2 + lambda + n - n^2 + 1im * epsilon * kappa * (-1 + 2 * n - 2 * s) +
                s + s^2 - 1im * tau + 2im * n * tau + tau^2
        prev2 = -2 * epsilon * kappa * (1im + epsilon - 1im * n + 1im * s - tau)
        denom = n * (epsilon + 1im * n - 1im * s + tau) * z^2
        coeffs[n + 1] = 1im * (prev2 * (n >= 2 ? coeffs[n - 1] : 0.0 + 0.0im) +
                               prev1 * z * coeffs[n]) / denom
    end

    R = convergence_radius_zero_in(coeffs, tol)
    return coeffs, R
end

function zero_out_coeffs(s, epsilon, tau, kappa, lambda, z, N::Integer; tol = 1e-15)
    coeffs = Vector{ComplexF64}(undef, N + 1)
    coeffs[1] = 1.0 + 0.0im
    N == 0 && return coeffs, 0.95

    for n in 1:N
        prev1 = -2 * epsilon^2 * kappa + lambda - (n - 1) * (n + 2 * s) +
                1im * epsilon * (1 - 2 * n - 2 * s + kappa * (-1 + 2 * n + 2im * tau))
        prev2 = -2 * epsilon * kappa * (2 * epsilon - 1im * (n - 1))
        denom = n * (1im * epsilon + n + s + 1im * tau) * z^2
        coeffs[n + 1] = (prev2 * (n >= 2 ? coeffs[n - 1] : 0.0 + 0.0im) +
                         prev1 * z * coeffs[n]) / denom
    end

    R = convergence_radius_zero_out(coeffs, tol)
    return coeffs, R
end

function convergence_radius_zero_in(coeffs, tol = 1e-15)
    N = length(coeffs) - 1
    N <= 0 && return 0.0
    last = coeffs[end]
    abs(last) == 0 && return 0.6
    return min((tol / abs(last))^(1.0 / N), 0.6)
end

function convergence_radius_zero_out(coeffs, tol = 1e-15)
    N = length(coeffs) - 1
    N <= 0 && return 0.0
    last = coeffs[end]
    abs(last) == 0 && return 0.6
    return min((tol / abs(last))^(1.0 / N), 0.6)
end

function evaluate_series(x, coeffs)
    n = length(coeffs)
    n == 0 && return (0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)
    val  = coeffs[end]
    dval = zero(val)
    ddval = zero(val)
    @inbounds for i in (n-1):-1:1
        ddval = muladd(ddval, x, 2dval)
        dval  = muladd(dval,  x, val)
        val   = muladd(val,   x, coeffs[i])
    end
    return (val, dval, ddval)
end

function evaluate_zero_in(x, coeffs)
    return evaluate_series(x, coeffs)
end

function evaluate_zero_out(x, s, epsilon, tau, coeffs)
    # P0_out(x) = (-x)^(s+i*(epsilon+tau)) * sum(b_n * x^n)
    P, P_prime, P_primeprime = evaluate_series(x, coeffs)
    prefactor = (-x)^(s + 1im * (epsilon + tau))
    prefactor_prime = (s + 1im * (epsilon + tau)) * prefactor / x
    prefactor_primeprime = (s + 1im * (epsilon + tau) - 1) * prefactor_prime / x
    P0out = prefactor * P
    P0out_prime = prefactor_prime * P + prefactor * P_prime
    P0out_primeprime = prefactor_primeprime * P + 2 * prefactor_prime * P_prime + prefactor * P_primeprime
    return (P0out, P0out_prime, P0out_primeprime)
end

function evaluate_zero_out_rho(rho, s, epsilon, tau, coeffs, z)
    sigma = s + 1im * (epsilon + tau)
    y, yp, ypp = evaluate_frobenius_series(rho, sigma, coeffs)
    pref = z^(-sigma)
    return (pref * y, pref * yp, pref * ypp)
end

function error_zero_in(x, s, epsilon, tau, kappa, lambda, coeffs)
    P, P_prime, P_primeprime = evaluate_zero_in(x, coeffs)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * P_primeprime
    term2 = q1 * P_prime
    term3 = q2 * P
    return _relative_residual(term1, term2, term3)
end

function error_zero_in((x, P, P_prime, P_primeprime), s, epsilon, tau, kappa, lambda)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * P_primeprime
    term2 = q1 * P_prime
    term3 = q2 * P
    return _relative_residual(term1, term2, term3)
end

function error_zero_out(x, s, epsilon, tau, kappa, lambda, coeffs)
    P, P_prime, P_primeprime = evaluate_zero_out(x, s, epsilon, tau, coeffs)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * P_primeprime
    term2 = q1 * P_prime
    term3 = q2 * P
    return _relative_residual(term1, term2, term3)
end

function error_zero_out((x, P, P_prime, P_primeprime), s, epsilon, tau, kappa, lambda)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * P_primeprime
    term2 = q1 * P_prime
    term3 = q2 * P
    return _relative_residual(term1, term2, term3)
end

function evaluate_zero_in_rho(rho, coeffs)
    return evaluate_frobenius_series(rho, coeffs)
end

function error_zero_in_rho(rho, s, epsilon, tau, kappa, lambda, z, coeffs)
    P, Pp, Ppp = evaluate_zero_in_rho(rho, coeffs)
    term1 = rho * (z + rho) * Ppp
    term2 = contour_q1(rho, s, epsilon, tau, kappa, z) * Pp
    term3 = contour_q2(rho, s, epsilon, tau, kappa, lambda, z) * P
    return _relative_residual(term1, term2, term3)
end

function error_zero_out_rho(rho, s, epsilon, tau, kappa, lambda, z, coeffs)
    P, Pp, Ppp = evaluate_zero_out_rho(rho, s, epsilon, tau, coeffs, z)
    term1 = rho * (z + rho) * Ppp
    term2 = contour_q1(rho, s, epsilon, tau, kappa, z) * Pp
    term3 = contour_q2(rho, s, epsilon, tau, kappa, lambda, z) * P
    return _relative_residual(term1, term2, term3)
end

end
