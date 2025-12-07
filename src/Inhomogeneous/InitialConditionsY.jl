module InitialConditionsY

using SpecialFunctions
using HypergeometricFunctions

using ..Coordinates: r_from_rstar
using ..AsymptoticExpansionCoefficientsY

export Y_initial_values_inf_m2, Y_initial_values_hor_m2
export Y_and_Yp_ingoing_hor_m2, Y_and_Yp_outgoing_inf_m2
export Y_initial_values_inf_p2, Y_initial_values_hor_p2
export Y_and_Yp_ingoing_hor_p2, Y_and_Yp_outgoing_inf_p2
export Y_initial_values_inf_up_m2, Y_initial_values_inf_up_p2
export Y_and_Yp_outgoing_hor_m2, Y_and_Yp_outgoing_hor_p2

_DEFAULT_infinity_expansion_order = 6
_DEFAULT_horizon_expansion_order = 3
_DEFAULTDATATYPE = ComplexF64

# Coefficients ck for E-term
function compute_ck(a, b, kmax::Int; data_type=_DEFAULTDATATYPE)
    ν = data_type(- b - 0.5)
    c = zeros(data_type, kmax + 1)
    c[1] = data_type(1)  # c₀

    for k in 1:kmax
        s = 0.0 + 0.0im
        for m in 0:k - 1
            term1 = (1 - ν - 2 * b + m)
            term2 = (1 - ν - 2 * (a + 1) + m)
            term3 = (m - 1 - ν)

            e_km = (gamma(term1 + 2 + k - m) * (a - b)) / ((a + 1 - b) * (1 - b) * gamma(term1)) +
                   (gamma(term2 + 2 + k - m)) / ((b - a - 1) * a * gamma(term2)) -
                   (gamma(term3 + 2 + k - m) * (a - 1)) / ((1 - b) * a * gamma(term3))

            s += c[m + 1] * e_km
        end
        c[k + 1] = - s / (4 * k)
    end
    return c
end

# E-term: oscillatory term

function E_asym(a, b, z; data_type=_DEFAULTDATATYPE, reltol=1e-12)
    if z <= 0
        error("Asymptotic expansion is valid only for large negative z")
    end
    ν = - b - 0.5
    factor = (2 * π)^(-1/2) * 2^(- ν - 1/2)

    exp_phase1 = exp(2im * sqrt(z))
    exp_phase2 = exp(-2im * sqrt(z))

    s1 = 0.0 + 0.0im
    s2 = 0.0 + 0.0im
    kmax = 20
    ck = compute_ck(a, b, kmax)

    for k in 0:kmax
        s1 += ck[k + 1] * (2im * sqrt(z))^(ν - k)
        s2 += ck[k + 1] * (-2im * sqrt(z))^(ν - k)
        if abs(ck[k + 1] * (2im * sqrt(z))^(ν - k)) / abs(s1) < reltol || !isfinite(abs(ck[k + 1] * (2im * sqrt(z))^(ν - k)))
            break
        end
    end
    return factor * (exp_phase1 * s1 + exp_phase2 * s2)
end

# H-term: algebraic term
function H_asym(a, b, z; reltol=1e-12)
    if z <= 0
        error("Asymptotic expansion is valid only for large negative z")
    end
    s = 0.0 + 0.0im
    k = 0
    term = 1.0 + 0.0im
    while abs(term) > abs(s) * reltol || k == 0
        num = (-1)^k * gamma(a + k)
        den = factorial(k) * gamma(b - a - k) * gamma(1 - k)
        term = num / den * z^(- a - k)
        s += term
        k += 1
    end
    return s
end

# Full asymptotic form
function my_1F2_asym(a, b, z; reltol=1e-12)
    if z >= 0
        error("Asymptotic expansion is valid only for large negative z")
    end
    c = a + 1
    pre = gamma(b) * a
    threshold = -100.0

    if z > threshold
        # Use direct hypergeometric series for small |z|
        return pFq([a], [b, c], z)
    else
        H = H_asym(a, b, -z; reltol=reltol)
        E = E_asym(a, b, -z; reltol=reltol)

        return pre * (H + E)
    end
end

function y_n_p(omega, r, order::Int; data_type=_DEFAULTDATATYPE)
    n = order
    if n < 2
        error("The order must be at least 2.")
    end

    hyper_term1 = my_1F2_asym((1-n)/2+2im*omega, 1/2, -omega^2*r^2)
    hyper_term2 = my_1F2_asym((2-n)/2+2im*omega, 3/2, -omega^2*r^2)
    sin_term = sin(π*(n-4im*omega)/2)
    cos_term = cos(π*(n-4im*omega)/2)

    y = r^(1-n+4im*omega) * (hyper_term1/(n-1-4im*omega) + 2im*omega*r*hyper_term2/(n-2-4im*omega)) + 
        (sin_term + 1im*cos_term * sign(omega)) * gamma(1-n+4im*omega) * abs(2*omega)^(n-1-4im*omega)
    return y
end

function y_n_m(omega, r, order::Int; data_type=_DEFAULTDATATYPE)
    n = order
    if n < 2
        error("The order must be at least 2.")
    end

    hyper_term1 = my_1F2_asym((1-n)/2-2im*omega, 1/2, -omega^2*r^2)
    hyper_term2 = my_1F2_asym((2-n)/2-2im*omega, 3/2, -omega^2*r^2)
    sin_term = sin(π*(n+4im*omega)/2)
    cos_term = cos(π*(n+4im*omega)/2)

    y = r^(1-n-4im*omega) * (hyper_term1/(n-1+4im*omega) - 2im*omega*r*hyper_term2/(n-2+4im*omega)) + 
        (sin_term - 1im*cos_term * sign(omega)) * gamma(1-n-4im*omega) * abs(2*omega)^(n-1+4im*omega)
    return y
end


function Yp_outgoing_inf_m2(func, r, order)
    ans = 0.0 + 0.0im
    for n in 0:order
        ans += - func(n) / ((n + 2) * r^(n + 2))
    end
    return ans
end

function Yp_outgoing_inf_p2(func, omega, r, order)
    ans = 0.0 + 0.0im
    for n in 0:order
        ans += - func(n) * y_n_p(omega, r, n+3)
    end
    return ans
end

function Y_outgoing_inf_m2(func, r, order)
    ans = 0.0 + 0.0im
    for n in 0:order
        ans += func(n) / ((n + 2) * (n + 1) * r^(n + 1))
    end
    return ans
end

function Y_outgoing_inf_p2(func, omega, r, order)
    ans = 0.0 + 0.0im
    for n in 0:order
        ans += func(n) * (y_n_p(omega, r, n+2) - r * y_n_p(omega, r, n+3))
    end
    return ans
end

function Yp_ingoing_inf_m2(func, omega, r, order)
    ans = 0.0 + 0.0im
    for n in 0:order
        ans += - func(n) * y_n_m(omega, r, n+3)
    end
    return ans
end

function Yp_ingoing_inf_p2(func, r, order)
    ans = 0.0 + 0.0im
    for n in 0:order
        ans += - func(n) / ((n + 2) * r^(n + 2))
    end
    return ans
end

function Y_ingoing_inf_m2(func, omega, r, order)
    ans = 0.0 + 0.0im
    for n in 0:order
        ans += func(n) * (y_n_m(omega, r, n+2) - r * y_n_m(omega, r, n+3))
    end
    return ans
end

function Y_ingoing_inf_p2(func, r, order)
    ans = 0.0 + 0.0im
    for n in 0:order
        ans += func(n) / ((n + 2) * (n + 1) * r^(n + 1))
    end
    return ans
end

function Yp_outgoing_hor_m2(func, a, r, order)
    ans = 0.0 + 0.0im
    rp = 1 + sqrt(1 - a^2)
    for n in 0:order
        ans += func(n) * (r - rp) ^ (n + 1) / (n + 1)
    end
    return ans
end

function Yp_outgoing_hor_p2(func, a, omega, m, r, order)
    ans = 0.0 + 0.0im
    rp = 1 + sqrt(1 - a^2)
    q = (a * rp * m + 2 * a^2 * omega - 4 * rp * omega) / (rp * sqrt(1 - a^2))
    for n in 0:order
        ans += func(n) * (r - rp) ^ (n + 1 - im * q) / (n + 1 - im * q)
    end
    return ans
end

function Y_outgoing_hor_m2(func, a, r, order)
    ans = 0.0 + 0.0im
    rp = 1 + sqrt(1 - a^2)
    for n in 0:order
        ans += func(n) * (r - rp) ^ (n + 2) / ((n + 1) * (n + 2))
    end
    return ans
end

function Y_outgoing_hor_p2(func, a, omega, m, r, order)
    ans = 0.0 + 0.0im
    rp = 1 + sqrt(1 - a^2)
    q = (a * rp * m + 2 * a^2 * omega - 4 * rp * omega) / (rp * sqrt(1 - a^2))

    for n in 0:order
        ans += func(n) * (r - rp) ^ (n + 2 - im * q) / ((n + 1 - im * q) * (n + 2 - im * q))
    end
    return ans
end

function Yp_ingoing_hor_m2(func, a, omega, m, r, order)
    ans = 0.0 + 0.0im
    rp = 1 + sqrt(1 - a^2)
    q = (a * rp * m + 2 * a^2 * omega - 4 * rp * omega) / (rp * sqrt(1 - a^2))
    for n in 0:order
        ans += func(n) * (r - rp) ^ (n + 1 + im * q) / (n + 1 + im * q)
    end
    return ans
end

function Yp_ingoing_hor_p2(func, a, r, order)
    ans = 0.0 + 0.0im
    rp = 1 + sqrt(1 - a^2)
    for n in 0:order
        ans += func(n) * (r - rp) ^ (n + 1) / (n + 1)
    end
    return ans
end

function Y_ingoing_hor_m2(func, a, omega, m, r, order)
    ans = 0.0 + 0.0im
    rp = 1 + sqrt(1 - a^2)
    q = (a * rp * m + 2 * a^2 * omega - 4 * rp * omega) / (rp * sqrt(1 - a^2))

    for n in 0:order
        ans += func(n) * (r - rp) ^ (n + 2 + im * q) / ((n + 1 + im * q) * (n + 2 + im * q))
    end

    return ans
end

function Y_ingoing_hor_p2(func, a, r, order)
    ans = 0.0 + 0.0im
    rp = 1 + sqrt(1 - a^2)

    for n in 0:order
        ans += func(n) * (r - rp) ^ (n + 2) / ((n + 1) * (n + 2))
    end
    return ans
end

function Y_initial_values_inf_m2(m::Int, a, omega, lambda, Binc, Bref, rsout; order::Int=-1)
    #=
    The default order is set to 3 since we have derived the analytical expressions of the them.
    If you want to use a higher order (>3), you can set it by passing the `order` keyword argument.
    =#
    _default_order = _DEFAULT_infinity_expansion_order
    order = (order == -1 ? _default_order : order)

    Ypp_out(ord) = Ypp_outgoing_inf_m2(m, a, omega, lambda, ord)
    Ypp_in(ord) = Ypp_ingoing_inf_m2(m, a, omega, lambda, ord)

    rout = r_from_rstar(a, rsout)
    Yp_out(r) = Yp_outgoing_inf_m2(Ypp_out, r, order)
    Y_out(r) = Y_outgoing_inf_m2(Ypp_out, r, order)
    Yp_in(r) = Yp_ingoing_inf_m2(Ypp_in, omega, r, order)
    Y_in(r) = Y_ingoing_inf_m2(Ypp_in, omega, r, order)

    Yp_initial = (Bref * Yp_out(rout) + Binc * Yp_in(rout)) * (rout^2 - 2 * rout +a^2) / (rout^2 + a^2)
    Y_initial = Bref * Y_out(rout) + Binc * Y_in(rout)

    return Y_initial, Yp_initial
end

function Y_initial_values_inf_p2(m::Int, a, omega, lambda, Binc, Bref, rsout; order::Int=-1)
    #=
    The default order is set to 3 since we have derived the analytical expressions of the them.
    If you want to use a higher order (>3), you can set it by passing the `order` keyword argument.
    =#
    _default_order = _DEFAULT_infinity_expansion_order
    order = (order == -1 ? _default_order : order)

    Ypp_out(ord) = Ypp_outgoing_inf_p2(m, a, omega, lambda, ord)
    Ypp_in(ord) = Ypp_ingoing_inf_p2(m, a, omega, lambda, ord)

    rout = r_from_rstar(a, rsout)
    Yp_out(r) = Yp_outgoing_inf_p2(Ypp_out, omega, r, order)
    Y_out(r) = Y_outgoing_inf_p2(Ypp_out, omega, r, order)
    Yp_in(r) = Yp_ingoing_inf_p2(Ypp_in, r, order)
    Y_in(r) = Y_ingoing_inf_p2(Ypp_in, r, order)

    Yp_initial = (Bref * Yp_out(rout) + Binc * Yp_in(rout)) * (rout^2 - 2 * rout +a^2) / (rout^2 + a^2)
    Y_initial = Bref * Y_out(rout) + Binc * Y_in(rout)

    return Y_initial, Yp_initial
end

function Y_initial_values_hor_m2(m::Int, a, omega, lambda, Cinc, Cref, rsin; order::Int=-1, ifrs = false)
    #=
    The default order is set to 3 since we have derived the analytical expressions of the them.
    Currently, the horizon expansion is truncated up to order of 3.
    =#
    _default_order = _DEFAULT_horizon_expansion_order
    order = (order == -1 ? _default_order : order)

    Ypp_out(ord) = Ypp_outgoing_hor_m2(m, a, omega, lambda, ord)
    Ypp_in(ord) = Ypp_ingoing_hor_m2(m, a, omega, lambda, ord)

    rin = r_from_rstar(a, rsin)
    Yp_out(r) = Yp_outgoing_hor_m2(Ypp_out, a, r, order)
    Y_out(r) = Y_outgoing_hor_m2(Ypp_out, a, r, order)
    Yp_in(r) = Yp_ingoing_hor_m2(Ypp_in, a, omega, m, r, order)
    Y_in(r) = Y_ingoing_hor_m2(Ypp_in, a, omega, m, r, order)

    if ifrs
        return Cref * Y_in(rin) + Cinc * Y_out(rin), Cref * Yp_in(rin) + Cinc * Yp_out(rin)
    else
        return Cref * Y_in(rin) + Cinc * Y_out(rin), (Cref * Yp_in(rin) + Cinc * Yp_out(rin)) * (rin^2 - 2 * rin +a^2) / (rin^2 + a^2)
    end
end

function Y_initial_values_inf_up_m2(m::Int, a, omega, lambda, rsout; order::Int=-1)
    #=
    The default order is set to 3 since we have derived the analytical expressions of the them.
    Currently, the horizon expansion is truncated up to order of 3.
    =#
    _default_order = _DEFAULT_horizon_expansion_order
    order = (order == -1 ? _default_order : order)

    Ypp_out(ord) = Ypp_outgoing_inf_m2(m, a, omega, lambda, ord)

    rout = r_from_rstar(a, rsout)
    Yp_out(r) = Yp_outgoing_inf_m2(Ypp_out, r, order)
    Y_out(r) = Y_outgoing_inf_m2(Ypp_out, r, order)

    return Y_out(rout), Yp_out(rout) * (rout^2 - 2 * rout +a^2) / (rout^2 + a^2)
end

function Y_initial_values_inf_up_p2(m::Int, a, omega, lambda, rsout; order::Int=-1)
    #=
    The default order is set to 3 since we have derived the analytical expressions of the them.
    Currently, the horizon expansion is truncated up to order of 3.
    =#
    _default_order = _DEFAULT_horizon_expansion_order
    order = (order == -1 ? _default_order : order)

    Ypp_out(ord) = Ypp_outgoing_inf_p2(m, a, omega, lambda, ord)

    rout = r_from_rstar(a, rsout)
    Yp_out(r) = Yp_outgoing_inf_p2(Ypp_out, omega, r, order)
    Y_out(r) = Y_outgoing_inf_p2(Ypp_out, omega, r, order)

    return Y_out(rout), Yp_out(rout) * (rout^2 - 2 * rout +a^2) / (rout^2 + a^2)
end

function Y_initial_values_hor_p2(m::Int, a, omega, lambda, Cinc, Cref, rsin; order::Int=-1, ifrs = false)
    #=
    The default order is set to 3 since we have derived the analytical expressions of the them.
    Currently, the horizon expansion is truncated up to order of 3.
    =#
    _default_order = _DEFAULT_horizon_expansion_order
    order = (order == -1 ? _default_order : order)

    Ypp_out(ord) = Ypp_outgoing_hor_p2(m, a, omega, lambda, ord)
    Ypp_in(ord) = Ypp_ingoing_hor_p2(m, a, omega, lambda, ord)

    rin = r_from_rstar(a, rsin)
    Yp_out(r) = Yp_outgoing_hor_p2(Ypp_out, a, omega, m, r, order)
    Y_out(r) = Y_outgoing_hor_p2(Ypp_out, a, omega, m, r, order)
    Yp_in(r) = Yp_ingoing_hor_p2(Ypp_in, a, r, order)
    Y_in(r) = Y_ingoing_hor_p2(Ypp_in, a, r, order)

    if ifrs
        return Cref * Y_in(rin) + Cinc * Y_out(rin), Cref * Yp_in(rin) + Cinc * Yp_out(rin)
    else
        return Cref * Y_in(rin) + Cinc * Y_out(rin), (Cref * Yp_in(rin) + Cinc * Yp_out(rin)) * (rin^2 - 2 * rin +a^2) / (rin^2 + a^2)
    end
end

function Y_and_Yp_ingoing_hor_m2(m::Int, a, omega, lambda; order::Int=-1)
    _default_order = _DEFAULT_horizon_expansion_order
    order = (order == -1 ? _default_order : order)
    Ypp_in(ord) = Ypp_ingoing_hor_m2(m, a, omega, lambda, ord)

    Yp_in(r) = Yp_ingoing_hor_m2(Ypp_in, a, omega, m, r, order)
    Y_in(r) = Y_ingoing_hor_m2(Ypp_in, a, omega, m, r, order)

    return Y_in, Yp_in
end

function Y_and_Yp_outgoing_hor_m2(m::Int, a, omega, lambda; order::Int=-1)
    _default_order = _DEFAULT_horizon_expansion_order
    order = (order == -1 ? _default_order : order)
    Ypp_out(ord) = Ypp_outgoing_hor_m2(m, a, omega, lambda, ord)

    Yp_out(r) = Yp_outgoing_hor_m2(Ypp_out, a, r, order)
    Y_out(r) = Y_outgoing_hor_m2(Ypp_out, a, r, order)

    return Y_out, Yp_out
end

function Y_and_Yp_ingoing_hor_p2(m::Int, a, omega, lambda; order::Int=-1)
    _default_order = _DEFAULT_horizon_expansion_order
    order = (order == -1 ? _default_order : order)
    Ypp_in(ord) = Ypp_ingoing_hor_p2(m, a, omega, lambda, ord)

    Yp_in(r) = Yp_ingoing_hor_p2(Ypp_in, a, r, order)
    Y_in(r) = Y_ingoing_hor_p2(Ypp_in, a, r, order)

    return Y_in, Yp_in
end

function Y_and_Yp_outgoing_hor_p2(m::Int, a, omega, lambda; order::Int=-1)
    _default_order = _DEFAULT_horizon_expansion_order
    order = (order == -1 ? _default_order : order)
    Ypp_out(ord) = Ypp_outgoing_hor_p2(m, a, omega, lambda, ord)

    Yp_out(r) = Yp_outgoing_hor_p2(Ypp_out, a, omega, m, r, order)
    Y_out(r) = Y_outgoing_hor_p2(Ypp_out, a, omega, m, r, order)

    return Y_out, Yp_out
end

function Y_and_Yp_outgoing_inf_m2(m::Int, a, omega, lambda; order::Int=-1)
    _default_order = _DEFAULT_infinity_expansion_order
    order = (order == -1 ? _default_order : order)
    Ypp_out(ord) = Ypp_outgoing_inf_m2(m, a, omega, lambda, ord)

    Yp_out(r) = Yp_outgoing_inf_m2(Ypp_out, r, order)
    Y_out(r) = Y_outgoing_inf_m2(Ypp_out, r, order)

    return Y_out, Yp_out
end

function Y_and_Yp_outgoing_inf_p2(m::Int, a, omega, lambda; order::Int=-1)
    _default_order = _DEFAULT_infinity_expansion_order
    order = (order == -1 ? _default_order : order)
    Ypp_out(ord) = Ypp_outgoing_inf_p2(m, a, omega, lambda, ord)

    Yp_out(r) = Yp_outgoing_inf_p2(Ypp_out, omega, r, order)
    Y_out(r) = Y_outgoing_inf_p2(Ypp_out, omega, r, order)

    return Y_out, Yp_out
end

end
