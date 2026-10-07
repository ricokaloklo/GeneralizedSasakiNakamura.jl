module SmallFrequencyExpansion

include("TeukolskyTransformation.jl")
using .TeukolskyTransformation
using HypergeometricFunctions
using Roots
using SpecialFunctions

export sfe_in, sfe_out

_TruncatioN = 20
_TOLERANCE = 1e-13
# Below this |epsilon| (= 2|omega|) the renormalized angular momentum is taken from
# its O(epsilon^2) series instead of the root finder: the residual g(nu) is flat to
# round-off there and the secant solver fails, while the series error is O(epsilon^4).
const _NU_SERIES_EPS = 1e-6

"""
    NuSplit(l, delta)

The MST renormalized angular momentum nu = l + delta stored as an integer part `l` and
a residual `delta`.  For small frequencies delta = O(epsilon^2) falls below the
floating-point spacing of `l` (e.g. 2e-22 next to 2.0), and the MST recursion
coefficients contain factors (n + nu) and (n + nu + 1) that vanish exactly at
n = -l and n = -l - 1.  Carrying delta separately keeps those factors equal to
delta instead of rounding them to zero, which is what produced Inf/NaN below
omega ~ 1e-8.  Use `nuplus(nu, c)` for every expression of the form nu + c.
"""
struct NuSplit
    l::Int
    delta::Float64
end
@inline nuplus(nu::NuSplit, c) = (nu.l + c) + nu.delta
@inline nuplus(nu::Number, c) = nu + c
Base.:+(nu::NuSplit, c::Number) = nuplus(nu, c)
Base.:+(c::Number, nu::NuSplit) = nuplus(nu, c)
Base.:-(nu::NuSplit, c::Number) = nuplus(nu, -c)
Base.Float64(nu::NuSplit) = nu.l + nu.delta
Base.show(io::IO, nu::NuSplit) = print(io, "NuSplit(", nu.l, " + ", nu.delta, ")")

function mst_abc(nu, n, s, epsilon, tau, kappa, lambda)
    # nu + n and nu + n + 1 are built with nuplus so that the exact zeros at
    # n = -l and n = -l - 1 become the (tiny, non-zero) residual delta.
    nun = nuplus(nu, n)          # n + nu
    nun1 = nuplus(nu, n + 1)     # n + nu + 1
    alpha = 1im * epsilon * kappa * (nun1 + s + 1im * epsilon) * (nun1 + s - 1im * epsilon) * (nun1 + 1im * tau) /
            nun1 / (2 * nun1 + 1)
    beta = - lambda - s * (s + 1) + nun * nun1 + epsilon * (epsilon + tau * kappa) + epsilon * tau *
            kappa * (s^2 + epsilon^2) / nun / nun1
    gamma = - 1im * epsilon * kappa * (nun - s + 1im * epsilon) * (nun - s - 1im * epsilon) * (nun - 1im * tau) /
            nun / (2 * nun - 1)
    return alpha, beta, gamma
end

function backward_rec(nu, s, epsilon, tau, kappa, lambda, N = _TruncatioN)
    f = zeros(ComplexF64, N + 1)
    if N == 0 || N == 1
        error("N should be larger than 1.")
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
        error("N should be larger than 1.")
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
        # leading-order MST relation l(l + 1) = lambda + s(s + 1) + O(epsilon); l is an integer
        l = round(Int, sqrt(1 + 4 * (lambda + s * (s + 1)) + 4 * epsilon * (epsilon - tau * kappa)) / 2 - 0.5)
    end
    l = Int(l)
    # O(epsilon^2) term of nu - l (Mano, Suzuki & Takasugi), kept as a separate residual
    delta2 = (-2 - s^2 / (l * (l + 1)) + ((l + 1)^2 - s^2)^2 / ((2 * l + 1) * (2 * l + 2) * (2 * l + 3)) -
              (l^2 - s^2)^2 / ((2 * l - 1) * 2 * l * (2 * l + 1))) * epsilon^2 / (2 * l + 1)
    nu_series = NuSplit(l, delta2)
    g_split(nu_float) = g_nu(NuSplit(l, nu_float - l))[1]
    nu_conv = if abs(epsilon) < _NU_SERIES_EPS
        nu_series
    else
        try
            NuSplit(l, find_zero(g_split, l + delta2) - l)
        catch err
            err isa Roots.ConvergenceFailed || rethrow()
            nu_series
        end
    end
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
            a = nuplus(nu, n + 1 - s) + 1im * epsilon
            b = ComplexF64(2 * nuplus(nu, n + 1))
            # d/dz U(a, b, z) = -a U(a + 1, b + 1, z); the former three-term recurrence
            # divided by (b - 2) = 2 (n + nu), which vanishes at n = -l for integer nu.
            U[i] = HypergeometricFunctions.U(a, b, z)
            U1[i] = - a * HypergeometricFunctions.U(a + 1, b + 1, z)
            U2[i] = a * (a + 1) * HypergeometricFunctions.U(a + 2, b + 2, z)
        end

        for i in eachindex(coeffs)
            n = i - N
            a = nuplus(nu, n + 1 - s) + 1im * epsilon
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
            a = nuplus(nu, n + 1 + s) - 1im * epsilon
            b = ComplexF64(2 * nuplus(nu, n + 1))
            # U is evaluated at -z: d/dz U(a, b, -z) = a U(a + 1, b + 1, -z)
            U[i] = HypergeometricFunctions.U(a, b, - z)
            U1[i] = a * HypergeometricFunctions.U(a + 1, b + 1, - z)
            U2[i] = a * (a + 1) * HypergeometricFunctions.U(a + 2, b + 2, - z)
        end
        for i in eachindex(coeffs)
            n = i - N
            a = nuplus(nu, n + 1 + s) - 1im * epsilon
            factor = HypergeometricFunctions.pochhammer(nuplus(nu, 1 + s) - 1im * epsilon, n) / HypergeometricFunctions.pochhammer(nuplus(nu, 1 - s) + 1im * epsilon, n)
            S += factor * coeffs[i] * z^a * U[i]
            Sprime += factor * (- 2im * epsilon * kappa) * coeffs[i] * (a * z^(a - 1) * U[i] + z^a * U1[i])
            Sprimeprime += factor * (- 2im * epsilon * kappa)^2 * coeffs[i] * (a * (a - 1) * z^(a - 2) * U[i] +  2 * a * z^(a - 1) * U1[i] + z^a * U2[i])
            F += coeffs[i] * exp(1im * π * (nuplus(nu, n + 1 + s) - 1im * epsilon)) * factor
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
