# This file translates the radial continued fraction in qnm/radial.py,
# commit f3abd18e59828e7e7d75d07f20c7cbc87925edfa, lines 12-405.
# qnm is Copyright (c) 2019 Leo C. Stein and is distributed under the MIT
# License. The translated formulas retain that attribution and license notice.

@inline _real_type(::Type{Complex{T}}) where {T<:AbstractFloat} = T
@inline _real_type(::Type{T}) where {T<:AbstractFloat} = T
@inline _precision_bits(::Type{Float16}) = 11
@inline _precision_bits(::Type{Float32}) = 24
@inline _precision_bits(::Type{Float64}) = 53
@inline _precision_bits(::Type{BigFloat}) = precision(BigFloat)

# Horizon and infinity exponents of Leaver's series for nonextremal Kerr (M = 1).
function leaver_characteristic_exponents(omega, a, s::Int, m::Int)
    T = promote_type(typeof(real(omega)), typeof(real(a)))
    one_T = one(T)
    root = sqrt(one_T - a * a)
    iszero(root) && throw(DomainError(a,
        "Leaver's nonextremal recurrence requires abs(a) < 1."))
    rplus = one_T + root
    rminus = one_T - root
    sigma_plus = (2 * omega * rplus - m * a) / (2 * root)
    sigma_minus = (2 * omega * rminus - m * a) / (2 * root)
    zeta = complex(zero(T), one_T) * omega
    xi = -s - complex(zero(T), one_T) * sigma_plus
    eta = -complex(zero(T), one_T) * sigma_minus
    return (zeta=zeta, xi=xi, eta=eta, rplus=rplus, rminus=rminus,
        sigma_plus=sigma_plus, sigma_minus=sigma_minus)
end

# The five D coefficients of Leaver's three-term radial recurrence (qnm/radial.py).
function leaver_D_coefficients(omega, a, s::Int, m::Int, angular_A)
    exponents = leaver_characteristic_exponents(omega, a, s, m)
    T = promote_type(typeof(real(omega)), typeof(real(a)),
        typeof(real(angular_A)))
    root = sqrt(one(T) - a * a)
    zeta = exponents.zeta
    xi = exponents.xi
    eta = exponents.eta
    p = root * zeta
    alpha = one(T) + s + xi + eta - 2 * zeta + s
    gamma = one(T) + s + 2 * eta
    delta = one(T) + s + 2 * xi
    sigma = angular_A + a * a * omega * omega - 8 * omega * omega +
        p * (2 * alpha + gamma - delta) +
        (one(T) + s - (gamma + delta) / 2) *
        (s + (gamma + delta) / 2)
    return (
        delta,
        4 * p - 2 * alpha + gamma - delta - 2,
        2 * alpha - gamma + 2,
        alpha * (4 * p - delta) - sigma,
        alpha * (alpha - gamma + one(T)),
    )
end

# Recurrence coefficients alpha_n, beta_n, gamma_n built from D.
@inline function leaver_recurrence_coefficients(n::Integer, D)
    nT = n * one(real(D[1]))
    alpha_n = nT^2 + (D[1] + 1) * nT + D[1]
    beta_n = -2 * nT^2 + (D[2] + 2) * nT + D[4]
    gamma_n = nT^2 + (D[3] - 3) * nT + D[5] - D[3] + 2
    return (alpha=alpha_n, beta=beta_n, gamma=gamma_n)
end

# Continued fraction with inversion index `n_inv`, truncated at `depth` and
# closed by the ratio `tail`.
function leaver_cf_truncated(omega, a, s::Int, m::Int, angular_A,
        n_inv::Int; depth::Int=300, tail=one(omega))
    depth > n_inv || throw(ArgumentError("depth must exceed n_inv."))
    D = leaver_D_coefficients(omega, a, s, m, angular_A)
    lower = zero(omega)
    for index in 0:(n_inv - 1)
        coefficients = leaver_recurrence_coefficients(index, D)
        denominator = coefficients.beta - coefficients.gamma * lower
        iszero(denominator) && throw(DomainError(denominator,
            "The finite terminating side encountered a zero denominator."))
        lower = coefficients.alpha / denominator
    end
    upper = -tail
    for index in depth:-1:(n_inv + 1)
        coefficients = leaver_recurrence_coefficients(index, D)
        denominator = coefficients.beta - coefficients.alpha * upper
        iszero(denominator) && throw(DomainError(denominator,
            "The truncated infinite side encountered a zero denominator."))
        upper = coefficients.gamma / denominator
    end
    pivot = leaver_recurrence_coefficients(n_inv, D)
    return pivot.beta - pivot.gamma * lower - pivot.alpha * upper
end

# Coefficients of the large-n expansion of the continued-fraction tail in
# powers of 1/sqrt(n+1), used to close the truncated fraction.
function _leaver_tail_coefficients(D, order::Int=12)
    C = typeof(D[1])
    T = typeof(real(D[1]))
    coefficients = zeros(C, order + 3)
    shifted = similar(coefficients)
    coefficients[1] = -1
    # Principal sqrt: its cut lies on the negative imaginary frequency axis.
    # _polish_leaver_root clamps its difference stencil to the physical side.
    coefficients[2] = sqrt(-(D[1] + D[2] + D[3]))
    iszero(coefficients[2]) && throw(DomainError(coefficients[2],
        "The half-integer asymptotic tail requires nonzero frequency."))
    coefficient(k) = k < 0 ? zero(C) : coefficients[k + 1]
    product(k) = k < 0 ? zero(C) :
        sum(coefficients[i + 1] * shifted[k - i + 1] for i in 0:k)
    for j in 2:order
        # U(t/sqrt(1+t^2)): each t^p acquires (1+t^2)^(-p/2).
        fill!(shifted, zero(C))
        shifted[1] = coefficients[1]
        for p in 1:(j - 1)
            weight = one(T)
            for k in p:2:(j + 1)
                shifted[k + 1] += weight * coefficients[p + 1]
                r = (k - p) ÷ 2
                weight *= (-T(p) / 2 - r) / (r + 1)
            end
        end
        k = j + 1
        defect = product(k) + (D[1] + 1)*product(k - 2) +
            D[1]*product(k - 4) + 2coefficient(k) -
            (D[2] + 2)*coefficient(k - 2) - D[4]*coefficient(k - 4)
        k == 4 && (defect += D[5] - D[3] + 2)
        coefficients[j + 1] = -defect / (2coefficients[2])
    end
    return coefficients[1:(order + 1)]
end

# Truncated fraction with the asymptotic tail; the depth is doubled until two
# successive relative changes are below `tolerance`.
function _leaver_cf_asymptotic(omega, a, s, m, angular_A, n_inv;
        tolerance, minimum_iterations, maximum_iterations)
    T = typeof(real(omega))
    D = leaver_D_coefficients(omega, a, s, m, angular_A)
    coefficients = _leaver_tail_coefficients(D)
    depth = min(maximum_iterations, max(minimum_iterations, n_inv + 32))
    evaluate(N) = leaver_cf_truncated(omega, a, s, m, angular_A, n_inv;
        depth=N, tail=-evalpoly(inv(sqrt(T(N + 1))), coefficients))
    previous = evaluate(depth)
    previous_error = T(Inf)
    error = T(Inf)
    convergence_error = T(Inf)
    value = previous
    while depth < maximum_iterations
        depth = min(2depth, maximum_iterations)
        value = evaluate(depth)
        error = abs(value - previous) / max(one(T), abs(value))
        convergence_error = max(error, previous_error)
        if convergence_error <= tolerance
            return (value=value, error=convergence_error,
                iterations=depth, converged=true, D=D)
        end
        previous = value
        previous_error = error
    end
    return (value=value, error=convergence_error, iterations=depth,
        converged=false, D=D)
end

# Leaver continued fraction with inversion index `n_inv`; its zero in omega is a
# QNM. Uses the asymptotic tail for omega != 0 (`algorithm=:auto`) and the
# modified Lentz method otherwise. Returns (value, error, iterations, converged, D).
function leaver_cf_inversion(omega, a, s::Int, m::Int, angular_A,
        n_inv::Int; tolerance=nothing, minimum_iterations::Int=300,
        maximum_iterations::Int=4000, require_convergence::Bool=true,
        algorithm::Symbol=:auto)
    n_inv >= 0 || throw(ArgumentError("n_inv must be nonnegative."))
    maximum_iterations > minimum_iterations || throw(ArgumentError(
        "maximum_iterations must exceed minimum_iterations."))
    T = promote_type(typeof(real(omega)), typeof(real(a)),
        typeof(real(angular_A)))
    tol = tolerance === nothing ? T(1.0e-13) : T(tolerance)
    use_tail = algorithm === :asymptotic ||
        (algorithm === :auto && !iszero(omega))
    if use_tail
        result = _leaver_cf_asymptotic(omega, a, s, m, angular_A, n_inv;
            tolerance=tol, minimum_iterations, maximum_iterations)
        if require_convergence && !result.converged
            throw(LentzConvergenceError(result.error, result.iterations,
                maximum_iterations))
        end
        return result
    end
    D = leaver_D_coefficients(omega, a, s, m, angular_A)

    lower = zero(omega)
    for index in 0:(n_inv - 1)
        coefficients = leaver_recurrence_coefficients(index, D)
        denominator = coefficients.beta - coefficients.gamma * lower
        iszero(denominator) && throw(DomainError(denominator,
            "The finite terminating side encountered a zero denominator."))
        lower = coefficients.alpha / denominator
    end

    tiny_real = max(T(1.0e-30), eps(T)^2)
    tiny = complex(tiny_real, zero(T))
    fraction = tiny
    C = tiny
    D_lentz = zero(tiny)
    delta = complex(T(Inf), zero(T))
    converged = false
    iterations = 0
    n = n_inv

    for iteration in 1:maximum_iterations
        current = leaver_recurrence_coefficients(n, D)
        denominator = current.gamma
        iszero(denominator) && throw(DomainError(denominator,
            "The infinite continued fraction encountered gamma_n = 0."))
        a_term = -current.alpha / denominator
        n += 1
        next = leaver_recurrence_coefficients(n, D)
        iszero(next.gamma) && throw(DomainError(next.gamma,
            "The infinite continued fraction encountered gamma_(n+1) = 0."))
        b_term = next.beta / next.gamma

        D_new = b_term + a_term * D_lentz
        iszero(D_new) && (D_new = tiny)
        C_new = b_term + a_term / C
        iszero(C_new) && (C_new = tiny)
        D_new = inv(D_new)
        delta = C_new * D_new
        fraction *= delta
        iterations = iteration

        D_lentz = D_new
        C = C_new
        if iteration > minimum_iterations && abs(delta - one(delta)) < tol
            converged = true
            break
        end
    end

    error_estimate = T(abs(delta - one(delta)))
    if require_convergence && !converged
        throw(LentzConvergenceError(
            error_estimate, iterations, maximum_iterations))
    end
    pivot = leaver_recurrence_coefficients(n_inv, D)
    value = pivot.beta - pivot.gamma * lower + pivot.gamma * fraction
    return (value=value, error=error_estimate, iterations=iterations,
        converged=converged, D=D)
end
