module OrdinaryPointExpansion

include("TeukolskyTransformation.jl")
using .TeukolskyTransformation

export ordinary_point_coeffs, evaluate_ordinary_point_series, error_ordinary,
       convergence_radius_x0_left, convergence_radius_x0_right,
       convergence_radius_rho_left, convergence_radius_rho_right

const _TOLERANCE = 1e-15
const C0 = 0.0 + 0.0im
const C1 = 1.0 + 0.0im
const I1 = 1.0im


function ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, x0, P, P_prime, N::Integer)
    a = Vector{ComplexF64}(undef, N + 1)
    a[1] = P
    N == 0 && return a

    a[2] = P_prime
    N == 1 && return a

    x0m1 = x0 - 1
    abs(x0 * x0m1) ≤ eps(Float64) && error("Denominator in ordinary point recurrence is zero because x0=0 or x0=1. x0 must be an ordinary point.")

    ek = epsilon * kappa
    x02 = x0 * x0
    two_x0_minus_1 = 2 * x0 - 1
    common_denom = x0m1 * x0

    @inbounds for n in 2:N
        nc = Float64(n)
        nm1 = nc - 1
        denom = nm1 * nc * common_denom

        term_an_minus_1 = - nm1 * (
            s + nc * two_x0_minus_1 +
            I1 * (-I1 + epsilon + 2 * x02 * ek + tau - 2 * x0 * (ek + tau - I1))
        )

        term_an_minus_2 = (
            -2 - nc^2 + s + s^2 - epsilon^2 - 3I1 * ek - 2I1 * s * ek +
            6I1 * x0 * ek + 2I1 * s * x0 * ek + 2 * x0 * epsilon * ek +
            lambda + nc * (3 - 2I1 * two_x0_minus_1 * ek + 2I1 * tau) -
            3I1 * tau - 2 * x0 * ek * tau + tau^2
        )

        if n == 2
            a[3] = (term_an_minus_1 * a[2] + term_an_minus_2 * a[1]) / denom
        else
            term_an_minus_3 = 2 * ek * (2I1 - I1 * nc + I1 * s - tau + epsilon)
            a[n + 1] = (term_an_minus_1 * a[n] + term_an_minus_2 * a[n - 1] + term_an_minus_3 * a[n - 2]) / denom
        end
    end

    return a
end

function ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, z, rho0, P, P_prime, N::Integer)
    coeffs = Vector{ComplexF64}(undef, N + 1)
    coeffs[1] = P
    N == 0 && return coeffs
    coeffs[2] = P_prime
    N == 1 && return coeffs

    A0 = rho0 * (z + rho0)
    A1 = z + 2 * rho0
    A2 = 1.0 + 0.0im
    B0 = contour_q1(rho0, s, epsilon, tau, kappa, z)
    B1 = 2 * (1 - 1im * epsilon * kappa - 1im * tau) - 4im * epsilon * kappa * rho0 / z
    B2 = -2im * epsilon * kappa / z
    C0 = contour_q2(rho0, s, epsilon, tau, kappa, lambda, z)
    C1 = 2 * epsilon * kappa * (-1im + epsilon + 1im * s - tau) / z

    abs(A0) == 0 && error("ordinary-point expansion requires rho0 != 0 and rho0 != -z.")

    for n in 0:(N - 2)
        term1 = (n + 1) * (n * A1 + B0) * coeffs[n + 2]
        term2 = (n * (n - 1) * A2 + n * B1 + C0) * coeffs[n + 1]
        term3 = n >= 1 ? ((n - 1) * B2 + C1) * coeffs[n] : 0.0 + 0.0im
        coeffs[n + 3] = -(term1 + term2 + term3) / ((n + 1) * (n + 2) * A0)
    end

    return coeffs
end

function evaluate_ordinary_point_series(x, x0, coeffs::AbstractVector)
    n = length(coeffs)
    n == 0 && return (C0, C0, C0)

    Δx = x - x0

    # p(Δx) = sum_{j=0}^{N} coeffs[j+1] * Δx^j
    # simultaneous Horner for value / first derivative / second derivative
    val = coeffs[n]
    dval = zero(val)
    ddval = zero(val)

    @inbounds for i in (n - 1):-1:1
        ddval = muladd(ddval, Δx, 2 * dval)
        dval  = muladd(dval,  Δx, val)
        val   = muladd(val,   Δx, coeffs[i])
    end

    return (val, dval, ddval)
end

function error_ordinary(x, x0, s, epsilon, tau, kappa, lambda, coeffs)
    P0, P0_prime, P0_primeprime = evaluate_ordinary_point_series(x, x0, coeffs)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * P0_primeprime
    term2 = q1 * P0_prime
    term3 = q2 * P0
    return _relative_residual(term1, term2, term3)
end

function error_ordinary((x, P0, P0_prime, P0_primeprime), s, epsilon, tau, kappa, lambda)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * P0_primeprime
    term2 = q1 * P0_prime
    term3 = q2 * P0
    return _relative_residual(term1, term2, term3)
end

function error_ordinary(rho, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
    P0, P0_prime, P0_primeprime = evaluate_ordinary_point_series(rho, rho0, coeffs)
    term1 = rho * (z + rho) * P0_primeprime
    term2 = contour_q1(rho, s, epsilon, tau, kappa, z) * P0_prime
    term3 = contour_q2(rho, s, epsilon, tau, kappa, lambda, z) * P0
    return _relative_residual(term1, term2, term3)
end

function convergence_radius_x0_right(x0, s, epsilon, tau, kappa, lambda, coeffs, tol = _TOLERANCE)
    R = 0.5 * x0
    err = error_ordinary(x0 - R, x0, s, epsilon, tau, kappa, lambda, coeffs)
    while err > tol
        R = 0.5 * R
        err = error_ordinary(x0 - R, x0, s, epsilon, tau, kappa, lambda, coeffs)
        if abs(R) < 1e3 * tol
            break
        end
    end
    return - R
end

function convergence_radius_x0_left(x0, s, epsilon, tau, kappa, lambda, coeffs, tol = _TOLERANCE)
    R = 0.5 * x0
    err = error_ordinary(x0 + R, x0, s, epsilon, tau, kappa, lambda, coeffs)
    while err > tol
        R = 0.5 * R
        err = error_ordinary(x0 + R, x0, s, epsilon, tau, kappa, lambda, coeffs)
        if abs(R) < 1e3 * tol
            break
        end
    end
    return - R
end

function convergence_radius_rho_right(rho0, s, epsilon, tau, kappa, lambda, z, coeffs, tol = _TOLERANCE)
    R = max(0.5 * min(abs(rho0), abs(rho0 + z)), 1e-6)
    err = error_ordinary(rho0 + R, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
    while err > tol
        R *= 0.5
        err = error_ordinary(rho0 + R, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
    end
    return R
end

function convergence_radius_rho_left(rho0, s, epsilon, tau, kappa, lambda, z, coeffs, tol = _TOLERANCE)
    R = max(0.5 * min(abs(rho0), abs(rho0 + z)), 1e-6)
    err = error_ordinary(rho0 - R, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
    while err > tol
        R *= 0.5
        err = error_ordinary(rho0 - R, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
    end
    return R
end

end
