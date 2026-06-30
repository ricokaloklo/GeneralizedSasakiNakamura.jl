module LocalSolutionAtZero

include("TeukolskyTransformation.jl")
using .TeukolskyTransformation

export an_zero, bn_zero, an_zero_auto, bn_zero_auto, zero_in_coeffs, zero_out_coeffs,
       convergence_radius_an_zero, convergence_radius_bn_zero,
       convergence_radius_zero_in, convergence_radius_zero_out,
       evaluate_zero_in, evaluate_zero_out, error_zero_in, error_zero_out,
       evaluate_zero_in_rho, evaluate_zero_out_rho, error_zero_in_rho, error_zero_out_rho,
       zero_out_resonant_log_coeffs, evaluate_zero_out_resonant_log, error_zero_out_resonant_log,
       zero_in_resonant_log_coeffs, evaluate_zero_in_resonant_log, error_zero_in_resonant_log

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

@inline function _threshold_zero_polynomial_coeffs(s, epsilon, kappa, lambda)
    B0 = 1 - s
    B1 = -2 - 2im * epsilon * (1 - kappa)
    B2 = -2im * epsilon * kappa
    C0 = lambda + s + s^2 + 1im * epsilon * (-1 + kappa - 2 * kappa * s)
    C1 = 4 * epsilon^2 * kappa + 2im * epsilon * kappa * (-1 + s)
    return B0, B1, B2, C0, C1
end

@inline function _resonant_regular_source_coeff(j, h, B0, B1, B2)
    hj = j <= length(h) - 1 ? h[j + 1] : 0.0 + 0.0im
    hjm1 = j >= 1 ? h[j] : 0.0 + 0.0im
    hjm2 = j >= 2 ? h[j - 1] : 0.0 + 0.0im
    return (2 * j - 1 + B0) * hj + (-2 * (j - 1) + 1 + B1) * hjm1 + B2 * hjm2
end

function zero_out_resonant_log_coeffs(s, epsilon, kappa, lambda, N)
    s < 0 || error("zero_out_resonant_log_coeffs requires a negative spin weight. Got s=$s.")
    h = an_zero(s, epsilon, -epsilon, kappa, lambda, N)
    b = zeros(ComplexF64, N + 1)
    b[1] = 1.0 + 0.0im
    log_coefficient = 0.0 + 0.0im
    nres = -s
    B0, B1, B2, C0, C1 = _threshold_zero_polynomial_coeffs(s, epsilon, kappa, lambda)

    for k in 1:N
        D = (s + k) * k
        E = -(s + k - 1) * (s + k - 2) + B1 * (s + k - 1) + C0
        F = B2 * (s + k - 2) + C1
        prev1 = b[k]
        prev2 = k >= 2 ? b[k - 1] : 0.0 + 0.0im
        source = k >= nres ? _resonant_regular_source_coeff(k - nres, h, B0, B1, B2) : 0.0 + 0.0im

        if abs(D) < 1e-13
            abs(source) < 1e-13 && error("Resonant logarithmic source vanished at k=$k.")
            log_coefficient = -(E * prev1 + F * prev2) / source
        else
            b[k + 1] = -(E * prev1 + F * prev2 + log_coefficient * source) / D
        end
    end

    return h, b, log_coefficient
end

@inline function _resonant_source_coeff_exponent_s(j, h, s, B1, B2)
    hj = j <= length(h) - 1 ? h[j + 1] : 0.0 + 0.0im
    hjm1 = j >= 1 ? h[j] : 0.0 + 0.0im
    hjm2 = j >= 2 ? h[j - 1] : 0.0 + 0.0im
    return (s + 2j) * hj + (-2 * (s + j - 1) + 1 + B1) * hjm1 + B2 * hjm2
end

function zero_in_resonant_log_coeffs(s, epsilon, kappa, lambda, N)
    s > 0 || error("zero_in_resonant_log_coeffs requires a positive spin weight. Got s=$s.")
    b = bn_zero(s, epsilon, -epsilon, kappa, lambda, N)
    a = zeros(ComplexF64, N + 1)
    a[1] = 1.0 + 0.0im
    log_coefficient = 0.0 + 0.0im
    nres = s
    B0, B1, B2, C0, C1 = _threshold_zero_polynomial_coeffs(s, epsilon, kappa, lambda)
    phase = (-1.0 + 0.0im)^s

    for k in 1:N
        D = k * (k - s)
        E = -(k - 1) * (k - 2) + B1 * (k - 1) + C0
        F = B2 * (k - 2) + C1
        prev1 = a[k]
        prev2 = k >= 2 ? a[k - 1] : 0.0 + 0.0im
        source = k >= nres ? phase * _resonant_source_coeff_exponent_s(k - nres, b, s, B1, B2) : 0.0 + 0.0im

        if abs(D) < 1e-13
            abs(source) < 1e-13 && error("Resonant logarithmic source vanished at k=$k.")
            log_coefficient = -(E * prev1 + F * prev2) / source
        else
            a[k + 1] = -(E * prev1 + F * prev2 + log_coefficient * source) / D
        end
    end

    return a, b, log_coefficient
end

function _zero_best_radius(coeff, n, tol)
    cn = abs(coeff)
    return iszero(cn) ? Inf : (tol / cn)^(1 / n)
end

function _zero_optimal_prefix(coeffs, tol, patience; tail = 4, min_order = 10)
    Nmax = length(coeffs) - 1
    Nmax <= 0 && return coeffs, 0, 0.0
    Npeak = 1
    peak = abs(coeffs[2])
    bad = 0
    for n in 2:Nmax
        cn = abs(coeffs[n + 1])
        if cn > peak
            peak = cn
            Npeak = n
            bad = 0
        else
            bad += 1
        end
        bad >= patience && n > Npeak && break
    end
    Nuse = min(Nmax, max(min_order, Npeak + tail))
    return coeffs[1:Nuse + 1], Nuse, _zero_best_radius(coeffs[Nuse + 1], Nuse, tol)
end

const _ZERO_AUTO_TARGET_RADIUS = 0.1

function an_zero_auto(s, epsilon, tau, kappa, lambda, Nmax; tol = 1e-15, patience = 12, min_order = 10, target_radius = _ZERO_AUTO_TARGET_RADIUS)
    coeffs = ComplexF64[1.0 + 0.0im]
    Nopt = min_order
    Ropt = 0.0
    bad = 0
    for n in 1:Nmax
        term1_num = n^2 - s^2 + epsilon^2 + 1im * epsilon * kappa + s * (2im * epsilon * kappa - 1) -
                    lambda + tau * (1im - tau) - 1im * n * (2 * tau + 2 * epsilon * kappa - 1im)
        term2_num = - 2 * epsilon * kappa * (1im - 1im * n + 1im * s + epsilon - tau)
        denom = n * (n - s - 1im * (tau + epsilon))
        isapprox(abs(denom), 0.0, atol=eps(Float64)) && error("Denominator in a_n recurrence is zero for n=$n.")
        cn = n == 1 ? (term1_num / denom) * coeffs[n] :
             (term1_num / denom) * coeffs[n] + (term2_num / denom) * coeffs[n - 1]
        push!(coeffs, cn)
        if n >= min_order
            Rn = convergence_radius_an_zero(@view(coeffs[1:n + 1]), tol)
            if isfinite(Rn) && Rn >= target_radius
                return coeffs[1:n + 1], n, Rn
            end
            if isfinite(Rn) && Rn > Ropt
                Ropt = Rn
                Nopt = n
                bad = 0
            else
                bad += 1
            end
        elseif n == Nmax
            Nopt = n
            Ropt = convergence_radius_an_zero(coeffs, tol)
        end
        if n >= min_order && bad >= patience && n > Nopt
            break
        end
    end
    Nuse = min(length(coeffs) - 1, max(min_order, Nopt))
    return coeffs[1:Nuse + 1], Nuse, convergence_radius_an_zero(@view(coeffs[1:Nuse + 1]), tol)
end

function bn_zero_auto(s, epsilon, tau, kappa, lambda, Nmax; tol = 1e-15, patience = 12, min_order = 10, target_radius = _ZERO_AUTO_TARGET_RADIUS)
    coeffs = ComplexF64[1.0 + 0.0im]
    Nopt = min_order
    Ropt = 0.0
    bad = 0
    for n in 1:Nmax
        term1_num = n^2 - 1im * epsilon + 2im * s * (1im + epsilon) + n * (2 * s - 1 - 2im * epsilon * (kappa - 1)) +
                    epsilon * kappa * (1im + 2 * epsilon + 2 * tau) - lambda
        term2_num = 2im * kappa * epsilon * (n + 2im * epsilon - 1)
        denom = n * (n + s + 1im * (tau + epsilon))
        isapprox(denom, 0.0 + 0.0im, atol=eps(Float64)) && error("Denominator in b_n recurrence is zero for n=$n.")
        cn = n == 1 ? (term1_num / denom) * coeffs[n] :
             (term1_num / denom) * coeffs[n] + (term2_num / denom) * coeffs[n - 1]
        push!(coeffs, cn)
        if n >= min_order
            Rn = convergence_radius_bn_zero(@view(coeffs[1:n + 1]), tol)
            if isfinite(Rn) && Rn >= target_radius
                return coeffs[1:n + 1], n, Rn
            end
            if isfinite(Rn) && Rn > Ropt
                Ropt = Rn
                Nopt = n
                bad = 0
            else
                bad += 1
            end
        elseif n == Nmax
            Nopt = n
            Ropt = convergence_radius_bn_zero(coeffs, tol)
        end
        if n >= min_order && bad >= patience && n > Nopt
            break
        end
    end
    Nuse = min(length(coeffs) - 1, max(min_order, Nopt))
    return coeffs[1:Nuse + 1], Nuse, convergence_radius_bn_zero(@view(coeffs[1:Nuse + 1]), tol)
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

function evaluate_zero_out_resonant_log(x, s, h, b, log_coefficient)
    H, Hp, Hpp = evaluate_series(x, h)
    B, Bp, Bpp = evaluate_series(x, b)
    prefactor = (-x)^s
    prefactor_prime = s * prefactor / x
    prefactor_primeprime = (s - 1) * prefactor_prime / x
    logx = log(-x)

    Psing = prefactor * B
    Psing_prime = prefactor_prime * B + prefactor * Bp
    Psing_primeprime = prefactor_primeprime * B + 2 * prefactor_prime * Bp + prefactor * Bpp

    Plog = log_coefficient * H * logx
    Plog_prime = log_coefficient * (Hp * logx + H / x)
    Plog_primeprime = log_coefficient * (Hpp * logx + 2 * Hp / x - H / x^2)

    return Psing + Plog, Psing_prime + Plog_prime, Psing_primeprime + Plog_primeprime
end

function evaluate_zero_in_resonant_log(x, s, a, b, log_coefficient)
    A, Ap, App = evaluate_series(x, a)
    H, Hp, Hpp = evaluate_zero_out(x, s, 0.0, 0.0, b)
    logx = log(-x)

    Plog = log_coefficient * H * logx
    Plog_prime = log_coefficient * (Hp * logx + H / x)
    Plog_primeprime = log_coefficient * (Hpp * logx + 2 * Hp / x - H / x^2)

    return A + Plog, Ap + Plog_prime, App + Plog_primeprime
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

function error_zero_out_resonant_log(x, s, epsilon, kappa, lambda, h, b, log_coefficient)
    P, P_prime, P_primeprime = evaluate_zero_out_resonant_log(x, s, h, b, log_coefficient)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, -epsilon, lambda)
    term1 = (1 - x) * x * P_primeprime
    term2 = q1 * P_prime
    term3 = q2 * P
    return _relative_residual(term1, term2, term3)
end

function error_zero_in_resonant_log(x, s, epsilon, kappa, lambda, a, b, log_coefficient)
    P, P_prime, P_primeprime = evaluate_zero_in_resonant_log(x, s, a, b, log_coefficient)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, -epsilon, lambda)
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
