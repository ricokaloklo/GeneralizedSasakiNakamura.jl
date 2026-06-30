module LocalSolutionAtInfinity

include("TeukolskyTransformation.jl")
using .TeukolskyTransformation

export an_inf, bn_inf, inf_in_coeffs, inf_out_coeffs,
       evaluate_inf_in, evaluate_inf_out, error_inf_in, error_inf_out,
       convergence_radius_an_inf, convergence_radius_bn_inf

const C0 = 0.0 + 0.0im
const C1 = 1.0 + 0.0im
const I1 = 1.0im

function an_inf(s, epsilon, tau, kappa, lambda, tol=1e-15, N=500, patience=20)
    a = Vector{ComplexF64}(undef, N+1)
    a[1] = C1
    d = 2*kappa*epsilon
    abs(d) ≤ eps(Float64) && error("singular")
    Nopt = 0; Ropt = Inf; bad = 0
    for n in 1:N
        nc = Float64(n)
        den = d*nc
        t1 = -I1*nc^2 + nc*(I1+2I1*s+2*epsilon*(1+kappa)) + I1*(2*epsilon^2*kappa+lambda+I1*epsilon*(1+2s+kappa+2I1*kappa*tau))
        if n==1
            a[2] = (t1/den)*a[1]
        else
            t2 = (nc+2I1*epsilon-1)*(I1*nc-I1-I1*s-epsilon+tau)
            a[n+1] = (t1*a[n]+t2*a[n-1])/den
        end
        Rn = abs(a[n+1]/tol)^(1/n)
        if Rn < Ropt
            Ropt = Rn; Nopt = n; bad = 0
        else
            bad += 1
        end
        bad ≥ patience && n > Nopt && break
    end
    return a[1:Nopt+1], Nopt, Ropt
end

function bn_inf(s, epsilon, tau, kappa, lambda, tol=1e-15, N=500, patience=20)
    b = Vector{ComplexF64}(undef, N+1)
    b[1] = C1
    d = 2*kappa*epsilon
    abs(d) ≤ eps(Float64) && error("singular")
    Nopt = 0; Ropt = Inf; bad = 0
    for n in 1:N
        nc = Float64(n)
        den = d*nc
        t1 = I1*nc^2 - epsilon - epsilon*kappa + 2s*(epsilon+epsilon*kappa-I1) + nc*(2I1*s-I1+2*epsilon*(1+kappa)) - I1*lambda
        if n==1
            b[2] = (t1/den)*b[1]
        else
            t2 = -I1*(nc+2s-1)*(nc+s-I1*(epsilon+tau-I1))
            b[n+1] = (t1*b[n]+t2*b[n-1])/den
        end
        Rn = abs(b[n+1]/tol)^(1/n)
        if Rn < Ropt
            Ropt = Rn; Nopt = n; bad = 0
        else
            bad += 1
        end
        bad ≥ patience && n > Nopt && break
    end
    return b[1:Nopt+1], Nopt, Ropt
end

function evaluate_series_inverse_x(x, coeffs::AbstractVector)
    n = length(coeffs)
    n == 0 && return (C0, C0, C0)

    u = inv(x)

    # g(u) = sum_{j=0}^{N} coeffs[j+1] * u^j
    g = coeffs[n]
    gu = zero(g)
    guu = zero(g)

    @inbounds for i in (n - 1):-1:1
        guu = muladd(guu, u, 2 * gu)
        gu  = muladd(gu,  u, g)
        g   = muladd(g,   u, coeffs[i])
    end

    # u = 1/x
    up = - u^2
    upp = 2 * u^3

    fx = g
    fpx = gu * up
    fppx = guu * up^2 + gu * upp

    return (fx, fpx, fppx)
end

function evaluate_inf_in(x, s, epsilon, tau, coeffs)
    P, Pp, Ppp = evaluate_series_inverse_x(x, coeffs)
    α = -I1 * epsilon + s + I1 * tau - 1
    pref = (-x)^α
    pref_p = α * pref / x
    pref_pp = (α - 1) * pref_p / x
    y = pref * P
    yp = pref_p * P + pref * Pp
    ypp = pref_pp * P + 2 * pref_p * Pp + pref * Ppp
    return (y, yp, ypp)
end

function evaluate_inf_out(x, s, epsilon, tau, kappa, coeffs)
    P, Pp, Ppp = evaluate_series_inverse_x(x, coeffs)
    pref = exp(-2I1 * epsilon * kappa * x) * (-x)^(I1 * epsilon - s + I1 * tau - 1)
    pref_p = -(1 + s + I1 * epsilon * (2 * x * kappa - 1) - I1 * tau) * pref / x
    pref_pp = (2 + s^2 - epsilon^2 * (1 - 2 * x * kappa)^2 +
        s * (3 + 2I1 * epsilon * (2 * x * kappa - 1) - 2I1 * tau) - 3I1 * tau - tau^2 +
        epsilon * (-3I1 - 2 * tau + 4 * x * kappa * (I1 + tau))) * pref / x^2
    y = pref * P
    yp = pref_p * P + pref * Pp
    ypp = pref_pp * P + 2 * pref_p * Pp + pref * Ppp
    return (y, yp, ypp)
end

function error_inf_in(x, s, epsilon, tau, kappa, lambda, coeffs)
    P, Pp, Ppp = evaluate_series_inverse_x(x, coeffs)
    q3, q4 = calculate_q3_q4(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * Ppp
    term2 = q3 * Pp
    term3 = q4 * P
    return _relative_residual(term1, term2, term3)
end

function error_inf_in((x, P, Pp, Ppp), s, epsilon, tau, kappa, lambda)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * Ppp
    term2 = q1 * Pp
    term3 = q2 * P
    return _relative_residual(term1, term2, term3)
end

function error_inf_out(x, s, epsilon, tau, kappa, lambda, coeffs)
    P, Pp, Ppp = evaluate_series_inverse_x(x, coeffs)
    q5, q6 = calculate_q5_q6(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * Ppp
    term2 = q5 * Pp
    term3 = q6 * P
    return _relative_residual(term1, term2, term3)
end

function error_inf_out((x, P, Pp, Ppp), s, epsilon, tau, kappa, lambda)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    term1 = (1 - x) * x * Ppp
    term2 = q1 * Pp
    term3 = q2 * P
    return _relative_residual(term1, term2, term3)
end

function convergence_radius_an_inf(R, s, epsilon, tau, kappa, lambda, coeffs, tol = 1e-15)
    err = error_inf_in(- R, s, epsilon, tau, kappa, lambda, coeffs)
    while err > tol
        R *= 1.1
        err = error_inf_in(- R, s, epsilon, tau, kappa, lambda, coeffs)
    end
    return R
end

function convergence_radius_bn_inf(R, s, epsilon, tau, kappa, lambda, coeffs, tol = 1e-15)
    err = error_inf_out(- R, s, epsilon, tau, kappa, lambda, coeffs)
    while err > tol
        R *= 1.1
        err = error_inf_out(- R, s, epsilon, tau, kappa, lambda, coeffs)
    end
    return R
end

function inf_in_coeffs(s, epsilon, tau, kappa, lambda, z, Nmax::Integer = 500; tol = 1e-15, patience = 20)
    coeffs = Vector{ComplexF64}(undef, Nmax + 1)
    coeffs[1] = 1.0 + 0.0im
    d = 2 * epsilon * kappa
    isapprox(abs(d), 0.0, atol = eps(Float64)) && error("singular infinity recurrence: epsilon*kappa ≈ 0")
    Nopt = 0
    Ropt = Inf
    bad = 0
    for n in 1:Nmax
        nc = Float64(n)
        term1 = (-1im) * nc^2 +
                nc * (1im + 2im * s + 2 * epsilon * (1 + kappa)) +
                1im * (2 * epsilon^2 * kappa + lambda + 1im * epsilon * (1 + 2 * s + kappa + 2im * kappa * tau))
        denom = d * nc
        if n == 1
            coeffs[2] = (term1 / denom) * coeffs[1]
        else
            term2 = (nc + 2im * epsilon - 1) * (1im * nc - 1im - 1im * s - epsilon + tau)
            coeffs[n + 1] = (term1 * coeffs[n] + term2 * coeffs[n - 1]) / denom
        end
        Rn = abs(coeffs[n + 1] / tol)^(1 / n)
        if Rn < Ropt
            Ropt = Rn
            Nopt = n
            bad = 0
        else
            bad += 1
        end
        bad >= patience && n > Nopt && break
    end
    return coeffs[1:Nopt + 1], Nopt, Ropt
end

function inf_out_coeffs(s, epsilon, tau, kappa, lambda, z, Nmax::Integer = 500; tol = 1e-15, patience = 20)
    coeffs = Vector{ComplexF64}(undef, Nmax + 1)
    coeffs[1] = 1.0 + 0.0im
    d = 2 * epsilon * kappa
    isapprox(abs(d), 0.0, atol = eps(Float64)) && error("singular infinity recurrence: epsilon*kappa ≈ 0")
    Nopt = 0
    Ropt = Inf
    bad = 0
    for n in 1:Nmax
        nc = Float64(n)
        term1 = 1im * nc^2 +
                (-epsilon - epsilon * kappa + 2 * s * (epsilon + epsilon * kappa - 1im)) +
                nc * (2im * s - 1im + 2 * epsilon * (1 + kappa)) -
                1im * lambda
        denom = d * nc
        if n == 1
            coeffs[2] = (term1 / denom) * coeffs[1]
        else
            term2 = -1im * (nc + 2 * s - 1) * (nc + s - 1im * (-1im + epsilon + tau))
            coeffs[n + 1] = (term1 * coeffs[n] + term2 * coeffs[n - 1]) / denom
        end
        Rn = abs(coeffs[n + 1] / tol)^(1 / n)
        if Rn < Ropt
            Ropt = Rn
            Nopt = n
            bad = 0
        else
            bad += 1
        end
        bad >= patience && n > Nopt && break
    end
    return coeffs[1:Nopt + 1], Nopt, Ropt
end

function evaluate_inf_in(rho, s, epsilon, tau, coeffs, z)
    sigma = s - 1 + 1im * (tau - epsilon)
    u = -rho / z
    g, gu, guu = evaluate_inverse_power_series(u, coeffs)
    du = -1 / z
    g1 = gu * du
    g2 = guu * du^2
    pref = (rho / z)^sigma
    prefp = sigma * pref / rho
    prefpp = (sigma - 1) * prefp / rho
    y = pref * g
    yp = prefp * g + pref * g1
    ypp = prefpp * g + 2 * prefp * g1 + pref * g2
    return (y, yp, ypp)
end

function evaluate_inf_out(rho, s, epsilon, tau, kappa, coeffs, z)
    sigma = -s - 1 + 1im * (tau + epsilon)
    u = -rho / z
    g, gu, guu = evaluate_inverse_power_series(u, coeffs)
    du = -1 / z
    gp = gu * du
    gpp = guu * du^2
    pref = exp(2im * epsilon * kappa * rho / z) * (rho / z)^sigma
    prefp = pref * (2im * epsilon * kappa / z + sigma / rho)
    prefpp = pref * ((2im * epsilon * kappa / z + sigma / rho)^2 - sigma / rho^2)
    y = pref * g
    yp = prefp * g + pref * gp
    ypp = prefpp * g + 2 * prefp * gp + pref * gpp
    return (y, yp, ypp)
end

function error_inf_in(rho, s, epsilon, tau, kappa, lambda, coeffs, z)
    P, Pp, Ppp = evaluate_inf_in(rho, s, epsilon, tau, coeffs, z)
    term1 = rho * (z + rho) * Ppp
    term2 = contour_q1(rho, s, epsilon, tau, kappa, z) * Pp
    term3 = contour_q2(rho, s, epsilon, tau, kappa, lambda, z) * P
    return _relative_residual(term1, term2, term3)
end

function error_inf_out(rho, s, epsilon, tau, kappa, lambda, coeffs, z)
    P, Pp, Ppp = evaluate_inf_out(rho, s, epsilon, tau, kappa, coeffs, z)
    term1 = rho * (z + rho) * Ppp
    term2 = contour_q1(rho, s, epsilon, tau, kappa, z) * Pp
    term3 = contour_q2(rho, s, epsilon, tau, kappa, lambda, z) * P
    return _relative_residual(term1, term2, term3)
end

@inline evaluate_inf_in_rho(rho, s, epsilon, tau, coeffs, z) = evaluate_inf_in(rho, s, epsilon, tau, coeffs, z)
@inline evaluate_inf_out_rho(rho, s, epsilon, tau, kappa, coeffs, z) = evaluate_inf_out(rho, s, epsilon, tau, kappa, coeffs, z)
@inline error_inf_in_rho(rho, s, epsilon, tau, kappa, lambda, coeffs, z) = error_inf_in(rho, s, epsilon, tau, kappa, lambda, coeffs, z)
@inline error_inf_out_rho(rho, s, epsilon, tau, kappa, lambda, coeffs, z) = error_inf_out(rho, s, epsilon, tau, kappa, lambda, coeffs, z)

end
