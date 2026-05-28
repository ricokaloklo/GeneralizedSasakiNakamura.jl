module SmallFrequencyExpansion

include("TeukolskyTransformation.jl")
using .TeukolskyTransformation
using HypergeometricFunctions
using Roots
using SpecialFunctions

export sfe_in, sfe_out

_TruncatioN = 20
_TOLERANCE = 1e-14

function mst_abc(nu, n, s, epsilon, tau, kappa, lambda)
    alpha = 1im * epsilon * kappa * (n + nu + 1 + s + 1im * epsilon) * (n + nu + 1 + s - 1im * epsilon) * (n +
            nu + 1 + 1im * tau) / (n + nu + 1) / (2 * n + 2 * nu + 3)
    beta = - lambda - s * (s + 1) + (n + nu) * (n + nu + 1) + epsilon * (epsilon + tau * kappa) + epsilon * tau *
            kappa * (s^2 + epsilon^2) / (n + nu) / (n + nu + 1)
    gamma = - 1im * epsilon * kappa * (n + nu - s + 1im * epsilon) * (n + nu - s - 1im * epsilon) * (n +
            nu - 1im * tau) / (n + nu) / (2 * n + 2 * nu - 1)
    return alpha, beta, gamma
end

function backward_rec(nu, s, epsilon, tau, kappa, lambda, N = _TruncatioN)
    f = zeros(ComplexF64, N + 1)
    if N == 0 || N == 1
        @error("N should be larger than 1.")
    end
    f[N] = 1.0 + 0.0im

    for n = N-1:-1:1
        alpha, beta, gamma = mst_abc(nu, n, s, epsilon, tau, kappa, lambda)
        f[n] = - (beta * f[n + 1] + alpha * f[n + 2]) / gamma
    end
    return f, f[2] / f[1]
end

function forward_rec(nu, s, epsilon, tau, kappa, lambda, N = _TruncatioN)
    f = zeros(ComplexF64, N + 1)
    if N == 0 || N == 1
        @error("N should be larger than 1.")
    end
    f[2] = 1.0 + 0.0im

    for n = 3:N+1
        alpha, beta, gamma = mst_abc(nu, n - N - 2, s, epsilon, tau, kappa, lambda)
        f[n] = - (beta * f[n - 1] + gamma * f[n - 2]) / alpha
    end
    return f, f[N] / f[N+1]
end

function mst_series(s, epsilon, tau, kappa, lambda, N = _TruncatioN, l = nothing)
    function g_nu(nu)
        f_positive, R = backward_rec(nu, s, epsilon, tau, kappa, lambda, N)
        f_negative, L = forward_rec(nu, s, epsilon, tau, kappa, lambda, N)
        alpha0, beta0, gamma0 = mst_abc(nu, 0, s, epsilon, tau, kappa, lambda)
        f = zeros(ComplexF64, 2 * N - 1)
        for n = 1:N
            f[n] = f_negative[n + 1] / f_negative[N + 1]
            f[n + N - 1] = f_positive[n] / f_positive[1]
        end
        g = beta0 + alpha0 * R + gamma0 * L
        return real(g), f
    end
    if l == nothing
        l = sqrt(1 + 4 * lambda * 4 * s * (s + 1) + 4 * epsilon * (epsilon - tau * kappa)) / 2 - 0.5
    end
    nu_attempt = l + (-2 - s^2 / (l * (l + 1)) + ((l + 1)^2 - s^2)^2 / ((2 * l + 1) * (2 * l + 2) * (2 * l + 3)) -
                (l^2 - s^2)^2 / ((2 * l - 1) * 2 * l * (2 * l + 1))) * epsilon^2 / (2 * l + 1)
    nu_conv = find_zero(g_nu, nu_attempt)
    return nu_conv, g_nu(nu_conv)[2]
end

function sfe_in(s, epsilon, tau, kappa, lambda, N = _TruncatioN, l = nothing)
    if epsilon < 0.0
        nu, coeffs, P, R = sfe_in(s, -epsilon, -tau, kappa, lambda, N, l)
        eval_P = x -> conj.(P(x))
        return nu, coeffs, eval_P, R
    end
    N_check = Int64(floor(150 / log10(1 / abs(epsilon * kappa))))
    N = min(N_check, N)
    nu, coeffs = mst_series(s, epsilon, tau, kappa, lambda, N, l)
    F = sum(coeffs)
    function evaluate_P(x)
        z = 2im * epsilon * kappa * (1 - x)
        prefactor = (1 - x)^(-1im * epsilon + s + 1im * tau - 1)
        prefactorPrime = (-1im * epsilon + s + 1im * tau - 1) * prefactor / (x - 1)
        prefactorPrimePrime = (-1im * epsilon + s + 1im * tau - 1) * (-1im * epsilon +
                                 s + 1im * tau - 2) * prefactor / (x - 1)^2
        S = zero(ComplexF64)
        Sprime = zero(ComplexF64)
        Sprimeprime = zero(ComplexF64)

        U = zeros(ComplexF64, 2 * N - 1)
        U1 = zeros(ComplexF64, 2 * N - 1)
        U2 = zeros(ComplexF64, 2 * N - 1)
        for i in 1:2*N-1
            n = i - N
            a = n + nu + 1 - s + 1im * epsilon
            b = ComplexF64(2 * n + 2 * nu + 2)
            c = (1 + a - b) / (b - 2)
            U[i] = HypergeometricFunctions.U(a, b, z)
            if i == 1
                U1[i] = c * HypergeometricFunctions.U(a - 1, b - 2, z) / z - (c + (b - 1) / z) * U[i]
                U2[i] = c * ((2 + a - b) / (b - 4) * HypergeometricFunctions.U(a - 2, b - 4, z) / z - ((2 + a - b) / (b - 4) + (b - 3) / z) *
                            HypergeometricFunctions.U(a - 1, b - 2, z)) / z - (c + (b - 1) / z) * U1[i] - c *
                            HypergeometricFunctions.U(a - 1, b - 2, z) / z^2 + (b - 1) * U[i] / z^2
            else
                U1[i] = c * U[i - 1] / z - (c + (b - 1) / z) * U[i]
                U2[i] = c * U1[i - 1] / z - (c + (b - 1) / z) * U1[i] - c * U[i - 1] / z^2 + (b - 1) * U[i] / z^2
            end
        end

        for i in eachindex(coeffs)
            n = i - N
            a = n + nu + 1 - s + 1im * epsilon
            S += coeffs[i] * z^a * U[i]
            Sprime += (- 2im * epsilon * kappa) * coeffs[i] * (a * z^(a - 1) * U[i] + z^a * U1[i])
            Sprimeprime += (- 2im * epsilon * kappa)^2 * coeffs[i] * (a * (a - 1) * z^(a - 2) * U[i] +  2 * a * z^(a - 1) * U1[i] + z^a * U2[i])
        end
        P = prefactor * S / F
        Pprime = (prefactor * Sprime + prefactorPrime * S) / F
        Pprimeprime = (prefactor * Sprimeprime + 2 * prefactorPrime * Sprime + prefactorPrimePrime * S) / F
        q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
        term1 = (1 - x) * x * Pprimeprime
        term2 = q1 * Pprime
        term3 = q2 * P
        error = abs(term1 + term2 + term3) / max(abs(term1), abs(term2), abs(term3))
        return P, Pprime, Pprimeprime, error
    end
    R = N / abs(10im * epsilon * kappa)
    return nu, coeffs, evaluate_P, R
end


function sfe_out(s, epsilon, tau, kappa, lambda, N = _TruncatioN, l = nothing)
    if epsilon < 0.0
        nu, coeffs, P, R = sfe_out(s, -epsilon, -tau, kappa, lambda, N, l)
        eval_P = x -> conj.(P(x))
        return nu, coeffs, eval_P, R
    end
    N_check = Int64(floor(150 / log10(1 / abs(epsilon * kappa))))
    N = min(N_check, N)
    nu, coeffs = mst_series(s, epsilon, tau, kappa, lambda, N, l)
    function evaluate_P(x)
        F = 0.0 + 0.0im
        z = 2im * epsilon * kappa * (1 - x)
        prefactor = exp(- 2im * epsilon * kappa * x) * (1 - x)^(1im * epsilon - s + 1im * tau - 1)
        prefactorPrime = - (1 + s + 1im * epsilon * (2 * (x - 1) * kappa - 1) - 1im * tau) * prefactor / (x - 1)
        prefactorPrimePrime = (2 + s^2 - epsilon^2 * (1 - 2 * (x - 1) * kappa)^2 + s * (3 + 1im * epsilon *
                                (4 * (x - 1) * kappa - 2) - 2im * tau) - 3im * tau - tau^2 + epsilon * (- 3im -
                                2 * tau + 4 * kappa * (x - 1) * (1im + tau))) * prefactor / (x - 1)^2
        S = zero(ComplexF64)
        Sprime = zero(ComplexF64)
        Sprimeprime = zero(ComplexF64)

        U = zeros(ComplexF64, 2 * N - 1)
        U1 = zeros(ComplexF64, 2 * N - 1)
        U2 = zeros(ComplexF64, 2 * N - 1)
        for i in 1:2*N-1
            n = i - N
            a = n + nu + 1 + s - 1im * epsilon
            b = ComplexF64(2 * n + 2 * nu + 2)
            c = (1 + a - b) / (b - 2)
            U[i] = HypergeometricFunctions.U(a, b, - z)
            if i == 1
                U1[i] = c * HypergeometricFunctions.U(a - 1, b - 2, - z) / z + (c - (b - 1) / z) * U[i]
                U2[i] = c * ((2 + a - b) / (b - 4) * HypergeometricFunctions.U(a - 2, b - 4, - z) / z + ((2 + a - b) / (b - 4) - (b - 3) / z) *
                            HypergeometricFunctions.U(a - 1, b - 2, - z)) / z + (c - (b - 1) / z) * U1[i] - c *
                            HypergeometricFunctions.U(a - 1, b - 2, - z) / z^2 + (b - 1) * U[i] / z^2
            else
                U1[i] = c * U[i - 1] / z + (c - (b - 1) / z) * U[i]
                U2[i] = c * U1[i - 1] / z + (c - (b - 1) / z) * U1[i] - c * U[i - 1] / z^2 + (b - 1) * U[i] / z^2
            end
        end
        for i in eachindex(coeffs)
            n = i - N
            a = n + nu + 1 + s - 1im * epsilon
            factor = HypergeometricFunctions.pochhammer(nu + 1 + s - 1im * epsilon, n) / HypergeometricFunctions.pochhammer(nu + 1 - s + 1im * epsilon, n)
            S += factor * coeffs[i] * z^a * U[i]
            Sprime += factor * (- 2im * epsilon * kappa) * coeffs[i] * (a * z^(a - 1) * U[i] + z^a * U1[i])
            Sprimeprime += factor * (- 2im * epsilon * kappa)^2 * coeffs[i] * (a * (a - 1) * z^(a - 2) * U[i] +  2 * a * z^(a - 1) * U1[i] + z^a * U2[i])
            F += coeffs[i] * exp(1im * π * (n + nu + 1 + s - 1im * epsilon)) * factor
        end
        P = prefactor * S / F
        Pprime = (prefactor * Sprime + prefactorPrime * S) / F
        Pprimeprime = (prefactor * Sprimeprime + 2 * prefactorPrime * Sprime + prefactorPrimePrime * S) / F
        q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
        term1 = (1 - x) * x * Pprimeprime
        term2 = q1 * Pprime
        term3 = q2 * P
        error = abs(term1 + term2 + term3) / max(abs(term1), abs(term2), abs(term3))
        return P, Pprime, Pprimeprime, error
    end
    R = N / abs(2im * epsilon * kappa)
    return nu, coeffs, evaluate_P, R
end

end
