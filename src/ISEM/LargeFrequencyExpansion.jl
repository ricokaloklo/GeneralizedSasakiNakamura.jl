module LargeFrequencyExpansion

using TaylorSeries

export lfe_f, lfe_gamma, _eta, _psi, A_series, B_series, _chi, _zeta, lfe_in_zero, lfe_out_zero, lfe_ordinary_point, X_to_P

_TruncatioN = 20
_TOLERANCE = 1e-13

function lfe_f(s, epsilon, tau, kappa, lambda)
    α0 = (s + 1im * (epsilon + tau) + 1) * (s + 1im * (epsilon + tau) - 1) / 4
    α1 = epsilon * (1 + kappa) * (epsilon + tau) + 1im * s * (1im * s + 1im + epsilon * kappa - tau) - lambda
    α2 = s * (s + 1) - 3im * s * epsilon * kappa - epsilon^2 * (1 + 3 * kappa + kappa^2) + lambda - epsilon * kappa * tau
    α3 = 2 * epsilon * kappa * (1im * s + epsilon + epsilon * kappa)
    α4 = - epsilon^2 * kappa^2
    f = x -> (α0 + α1 * x + α2 * x^2 + α3 * x^3 + α4 * x^4) / (1 - x)^2 / x^2

    return f
end

function lfe_gamma(s, epsilon, tau, kappa, lambda)
    α0 = (s + 1im * (epsilon + tau) + 1) * (s + 1im * (epsilon + tau) - 1) / 4
    α1 = epsilon * (1 + kappa) * (epsilon + tau) + 1im * s * (1im * s + 1im + epsilon * kappa - tau) - lambda
    α2 = s * (s + 1) - 3im * s * epsilon * kappa - epsilon^2 * (1 + 3 * kappa + kappa^2) + lambda - epsilon * kappa * tau
    α3 = 2 * epsilon * kappa * (1im * s + epsilon + epsilon * kappa)
    α4 = - epsilon^2 * kappa^2
    Λ = (s + 1im * (epsilon + tau) + 1) / 2
    return (α1/α0, α2/α0, α3/α0, α4/α0, Λ)
end

function _eta(ft)
    sq_ft = sqrt(ft)
    return integrate(sq_ft)
end

function _psi(γ, N = 20)
    γ1, γ2, γ3, γ4, Λ = γ
    xf_sq = taylor_expand(x -> sqrt(1 + γ1 * x + γ2 * x^2 + γ3 * x^3 + γ4 * x^4) / (1 - x), 0.0, order = N + 2)
    coeffs = zeros(ComplexF64, N + 3)
    coeffs[2] = 1.0 + 0.0im
    x = Taylor1(coeffs)
    y = taylor_expand(x -> (8 * (3 - 10 * x + 10 * x^2 + x * (1 - 4 * x + 6 * x^2) * γ1 + 3 * x^4 * γ2 + 2 * x^4 * γ3 + x^5 *
         γ3 + x^4 * γ4 + 2 * x^5 * γ4) / (1 + γ1 * x + γ2 * x^2 + γ3 * x^3 + γ4 * x^4)^2 - 5 * (-2 + 4 * x + x *
         (-1 + 3 * x) * γ1 + 2 * x^3 * γ2 + x^3 * γ3 + x^4 * γ3 + 2 * x^4 * γ4)^2 / (1 + γ1 * x + γ2 * x^2 +
         γ3 * x^3 + γ4 * x^4)^3 - 4) / 16, 0.0, order = N+2)
    ψ0 = Taylor1(zeros(ComplexF64, N+2))
    for i in 1:N+1
        ψ0[i] = y[i+1]
    end
    g = Taylor1(zeros(ComplexF64, N + 3))
    for i in 1:N+2
        g[i] = xf_sq[i] / i
    end
    return - x / xf_sq, - ψ0 * xf_sq, exp(-g)
end

function A_series(x0, f, N = 20)
    A = Vector{}(undef, N + 1)
    ft = taylor_expand(x -> (f(x) / f(x0)), x0, order = 2 * N + 2)
    fp = differentiate(ft)
    fpp = differentiate(fp)
    coeffs = zeros(ComplexF64,  2*N+1)
    coeffs[1] = 1.0 + 0.0im
    A[1] = Taylor1(coeffs)
    ϕ = (4 * ft * fpp - 5 * fp^2) / (16 * ft^3)
    for i in 2:N+1
        A[i] = - differentiate(A[i - 1]) / sqrt(ft) / 2 + integrate(ϕ * sqrt(ft) * A[i - 1] / 2)
    end
    return A, ft
end

function B_series(γ, N = 20)
    B = Vector{}(undef, N + 1)
    term1, term2, g = _psi(γ, 2 * N)
    coeffs = zeros(ComplexF64, 2*N+1)
    coeffs[1] = 1.0 + 0.0im
    B[1] = Taylor1(coeffs)
    for i in 2:N+1
        B[i] = B[i - 1] - differentiate(B[i - 1]) * term1 + integrate(B[i - 1] * term2)
    end
    return B, g
end

function _chi(β, η, A)
    N = size(A)[1]
    dη = differentiate(η)
    d2η = differentiate(dη)
    function chi1(xt)
        val1 = 0.0 + 0.0im
        dval1 = 0.0 + 0.0im
        d2val1 = 0.0 + 0.0im
        for i in 1:N
            Ai = A[i]
            dAi = differentiate(Ai)
            d2Ai = differentiate(dAi)
            val1 += Ai(xt) / β^(i - 1)
            dval1 += dAi(xt) / β^(i - 1)
            d2val1 += d2Ai(xt) / β^(i - 1)
        end
        χ1 = exp(β * η(xt)) * val1
        dχ1 = exp(β * η(xt)) * (dval1 + β * dη(xt) * val1)
        d2χ1 = exp(β * η(xt)) * (d2val1 + 2 * β * dη(xt) * dval1 + β * (d2η(xt) + β * dη(xt)^2) * val1)
        return (χ1, dχ1, d2χ1)
    end
    function chi2(xt)
        val2 = 0.0 + 0.0im
        dval2 = 0.0 + 0.0im
        d2val2 = 0.0 + 0.0im
        for i in 1:N
            Ai = A[i]
            dAi = differentiate(Ai)
            d2Ai = differentiate(dAi)
            val2 += Ai(xt) / (-β)^(i - 1)
            dval2 += dAi(xt) / (-β)^(i - 1)
            d2val2 += d2Ai(xt) / (-β)^(i - 1)
        end
        χ2 = exp(-β * η(xt)) * val2
        dχ2 = exp(-β * η(xt)) * (dval2 - β * dη(xt) * val2)
        d2χ2 = exp(-β * η(xt)) * (d2val2 - 2 * β * dη(xt) * dval2 - β * (d2η(xt) - β * dη(xt)^2) * val2)
        return (χ2, dχ2, d2χ2)
    end
    return chi1, chi2
end

function _zeta(Λ, g, B)
    N = size(B)[1]
    dg = differentiate(g)
    d2g = differentiate(dg)
    ξ = x -> - g(x) / x
    dξ = x -> - dg(x) / x + g(x) / x^2
    d2ξ = x -> - d2g(x) / x + 2 * dg(x) / x^2 - 2 * g(x) / x^3
    function zeta1(x)
        val1 = 0.0 + 0.0im
        dval1 = 0.0 + 0.0im
        d2val1 = 0.0 + 0.0im
        for i in 1:N
            Bi = B[i]
            dBi = differentiate(Bi)
            d2Bi = differentiate(dBi)
            val1 += Bi(x) / (2 * Λ)^(i - 1)
            dval1 += dBi(x) / (2 * Λ)^(i - 1)
            d2val1 += d2Bi(x) / (2 * Λ)^(i - 1)
        end
        ζ1 = ξ(x)^Λ * val1
        dζ1 = ξ(x)^(Λ - 1) * (ξ(x) * dval1 + Λ * dξ(x) * val1)
        d2ζ1 = ξ(x)^(Λ - 2) * (Λ * (Λ - 1) * dξ(x)^2 * val1 + ξ(x)^2 * d2val1 + Λ * ξ(x) * (2 * dval1 * dξ(x) + val1 * d2ξ(x)))
        return (ζ1, dζ1, d2ζ1)
    end
    function zeta2(x)
        val2 = 0.0 + 0.0im
        dval2 = 0.0 + 0.0im
        d2val2 = 0.0 + 0.0im
        for i in 1:N
            Bi = B[i]
            dBi = differentiate(Bi)
            d2Bi = differentiate(dBi)
            val2 += Bi(x) / (2 * (1 - Λ))^(i - 1)
            dval2 += dBi(x) / (2 * (1 - Λ))^(i - 1)
            d2val2 += d2Bi(x) / (2 * (1 - Λ))^(i - 1)
        end
        ζ2 = ξ(x)^(1-Λ) * val2
        dζ2 = ξ(x)^(-Λ) * (ξ(x) * dval2 - (Λ - 1) * dξ(x) * val2)
        d2ζ2 = ξ(x)^(-Λ - 1) * (Λ * (Λ - 1) * dξ(x)^2 * val2 + ξ(x)^2 * d2val2 - (Λ - 1) * ξ(x) * (2 * dval2 * dξ(x) + val2 * d2ξ(x)))
        return (ζ2, dζ2, d2ζ2)
    end
    return zeta1, zeta2
end

function lfe_in_zero(s, epsilon, tau, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE)
    γ = lfe_gamma(s, epsilon, tau, kappa, lambda)
    B, g = B_series(γ, N)
    Λ = γ[5]
    ζ1, ζ2 = _zeta(Λ, g, B)
    gt = Taylor1(zeros(ComplexF64, size(g)[1]))
    f = lfe_f(s, epsilon, tau, kappa, lambda)
    for i in 0:size(g)[1]-1
        gt[i] = (1 - i) * g[i]
    end
    coes = zeros(ComplexF64, size(g)[1])
    coes[2] += 1.0
    x =  Taylor1(coes)
    Z = x / sqrt(gt)
    dZ = differentiate(Z)
    d2Z = differentiate(dZ)
    prefactor =  (1 - 2 * Λ) / (2 * Λ)
    function X_in(x)
        ζ, dζ, d2ζ = ζ1(x)
        Xin = prefactor * Z(x) * ζ
        dXin = prefactor * (dZ(x) * ζ + Z(x) * dζ)
        d2Xin = prefactor * (d2Z(x) * ζ + 2 * dZ(x) * dζ + Z(x) * d2ζ)
        fx = f(x)
        error = abs(d2Xin - fx * Xin) / max(abs(d2Xin), abs(fx * Xin))
        return (Xin, dXin, d2Xin, error)
    end
    R = 0.6
    err = X_in(-R)[4]
    while err > tol
        R *= 0.9
        err = X_in(-R)[4]
    end
    return X_in, R
end

function lfe_out_zero(s, epsilon, tau, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE)
    γ = lfe_gamma(s, epsilon, tau, kappa, lambda)
    B, g = B_series(γ, N)
    Λ = γ[5]
    ζ1, ζ2 = _zeta(Λ, g, B)
    gt = Taylor1(zeros(ComplexF64, size(g)[1]))
    f = lfe_f(s, epsilon, tau, kappa, lambda)
    for i in 0:size(g)[1]-1
        gt[i] = (1 - i) * g[i]
    end
    coes = zeros(ComplexF64, size(g)[1])
    coes[2] += 1.0
    x = Taylor1(coes)
    Z = x / sqrt(gt)
    dZ = differentiate(Z)
    d2Z = differentiate(dZ)
    prefactor = (1 - 2 * Λ) / (2 * Λ - 2)
    function X_out(x)
        ζ, dζ, d2ζ = ζ2(x)
        Xout = prefactor * Z(x) * ζ
        dXout = prefactor * (dZ(x) * ζ + Z(x) * dζ)
        d2Xout = prefactor * (d2Z(x) * ζ + 2 * dZ(x) * dζ + Z(x) * d2ζ)
        fx = f(x)
        error = abs(d2Xout - fx * Xout) / max(abs(d2Xout), abs(fx * Xout))
        return (Xout, dXout, d2Xout, error)
    end
    R = 0.6
    err = X_out(-R)[4]
    while err > tol
        R *= 0.9
        err = X_out(-R)[4]
    end
    return X_out, R
end

function lfe_ordinary_point(s, epsilon, tau, kappa, lambda, x0, X0, dX0, N = _TruncatioN, tol = _TOLERANCE)
    N = min(N, _TruncatioN)
    f = lfe_f(s, epsilon, tau, kappa, lambda)
    β = sqrt(f(x0))
    A, ft = A_series(x0, f, N)
    η = _eta(ft)
    χ1, χ2 = _chi(β, η, A)
    Z = ft^(-1/4)
    dZ = differentiate(Z)
    d2Z = differentiate(dZ)
    Z0 = Z(0.0)
    dZ0 = dZ(0.0)
    χ1x0, dχ1x0, d2χ1x0 = χ1(0.0)
    χ2x0, dχ2x0, d2χ2x0 = χ2(0.0)
    X10 = Z0 * χ1x0
    dX10 = (dZ0 * χ1x0 + Z0 * dχ1x0)
    X20 = Z0 * χ2x0
    dX20 = (dZ0 * χ2x0 + Z0 * dχ2x0)
    C1 = - (dX0 * X20 - X0 * dX20) / (X10 * dX20 - X20 * dX10)
    C2 = - (dX10 * X0 - X10 * dX0) / (X10 * dX20 - X20 * dX10)
    function X(x)
        xt = x - x0
        χ10, dχ10, d2χ10 = χ1(xt)
        χ20, dχ20, d2χ20 = χ2(xt)
        X1 = Z(xt) * χ10
        dX1 = (dZ(xt) * χ10 + Z(xt) * dχ10)
        d2X1 = (d2Z(xt) * χ10 + 2 * dZ(xt) * dχ10 + Z(xt) * d2χ10)
        X2 = Z(xt) * χ20
        dX2 = (dZ(xt) * χ20 + Z(xt) * dχ20)
        d2X2 = (d2Z(xt) * χ20 + 2 * dZ(xt) * dχ20 + Z(xt) * d2χ20)
        Xx = C1 * X1 + C2 * X2
        dXx = C1 * dX1 + C2 * dX2
        d2Xx = C1 * d2X1 + C2 * d2X2
        fx = f(x)
        error = abs(d2Xx - fx * Xx) / max(abs(d2Xx), abs(fx * Xx))
        return (Xx, dXx, d2Xx, error)
    end
    R = - 0.5 * x0
    err = X(x0 - R)[4]
    while err > tol
        R *= 0.5
        err = X(x0 - R)[4]
    end
    return X, R
end

function X_to_P(X, s, epsilon, tau, kappa)
    f = x -> exp(- 1im * x * epsilon * kappa) * (-x)^((s + 1im * (1im + epsilon + tau)) / 2) *
        (1 - x)^((- 1 - s - 1im * epsilon + 1im * tau) / 2)
    fp = x -> f(x) * (-(s + 1im * (1im + epsilon + 2 * x^2 * epsilon * kappa + tau - 2 * x * (1im + epsilon * kappa +
        tau))) / (2 * x * (x - 1)))
    fpp = x -> f(x) * ((s^2 - 4 * x^4 * epsilon^2 * kappa^2 - (1im + epsilon + tau) * (3im + epsilon + tau) +
        8 * x^3 * epsilon * kappa * (1im + epsilon * kappa + tau) + 4 * x * (1im + epsilon + tau) * (2im +
        epsilon * kappa + tau) - 4 * x^2 * (- 2 + epsilon * kappa * (3im + epsilon + epsilon * kappa) + 3im * tau
        + 3 * epsilon * kappa * tau + tau^2) + 2im * s * (2im + epsilon + 2 * x^2 * epsilon * kappa + tau - 2 * x *
        (2im + epsilon * kappa + tau))) / (4 * x^2 * (x - 1)^2))
    function P(x)
        Xx, dXx, d2Xx, error = X(x)
        prefactor = f(x)
        prefactor_prime = fp(x)
        prefactor_primeprime = fpp(x)
        Px = prefactor * Xx
        dPx = prefactor_prime * Xx + prefactor * dXx
        d2Px = prefactor_primeprime * Xx + 2 * prefactor_prime * dXx + prefactor * d2Xx
        return (Px, dPx, d2Px, error)
    end
    return P
end

end
