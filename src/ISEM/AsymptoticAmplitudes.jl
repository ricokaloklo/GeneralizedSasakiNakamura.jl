module AsymptoticAmplitudes

include("TeukolskyTransformation.jl")
using .TeukolskyTransformation

export C_generic, C_inf, C_zero, _f1, _f2, _f3, _f4,
       B_trans, B_inc, B_ref, C_trans, C_ref, C_inc

@inline _scaled(x, factor) = x * factor

@inline
function _iota(x, s, epsilon, tau, kappa)
    ι = exp(2im * kappa * x * epsilon) * (1 - x)^(s - 1im * tau + 1im * epsilon + 1) * (- x)^(1 - s - 1im * tau - 1im * epsilon)
    return ι
end

@inline
function inner_product(P1, P2, x)
    P1x, P1x_prime, _, _ = P1(x)
    P2x, P2x_prime, _, _ = P2(x)
    return P1x * P2x_prime - P2x * P1x_prime
end

@inline
function C_generic(x, P1, P2, P3, P4)
    return inner_product(P1, P2, x) / inner_product(P3, P4, x)
end

@inline
function C_inf(x, P1, P2, s, epsilon, tau, kappa)
    return _iota(x, s, epsilon, tau, kappa) * inner_product(P1, P2, x) / (2im * epsilon * kappa)
end

@inline
function C_inf(rho, P1, P2, s, epsilon, tau, kappa, z)
    return (-z) * _iota(-rho/z, s, epsilon, tau, kappa) * inner_product(P1, P2, rho) / (2im * epsilon * kappa)
end

@inline
function C_zero(x, P1, P2, s, epsilon, tau, kappa)
    return _iota(x, s, epsilon, tau, kappa) * inner_product(P1, P2, x) / (s + 1im * (epsilon + tau))
end

@inline
function C_zero(rho, P1, P2, s, epsilon, tau, kappa, z)
    return (-z) * _iota(-rho/z, s, epsilon, tau, kappa) * inner_product(P1, P2, rho) / (s + 1im * (epsilon + tau))
end

@inline
function _f1(s, epsilon, tau, kappa)
    rp = 1 + kappa
    p = kappa * (epsilon + tau) / (2 * rp)
    return (4.0)^s * kappa^(2im * p + 2 * s) * exp(1im * p * rp)
end

@inline
function _f2(epsilon, tau, kappa)
    rp = 1 + kappa
    p = kappa * (epsilon + tau) / (2 * rp)
    return kappa^(- 2im * p) * exp(- 1im * p * rp)
end

@inline
function _f3(epsilon, kappa)
    ω = epsilon / 2
    rp = 1 + kappa
    return 2.0 * kappa^(2im * ω + 1) * exp(1im * ω * rp)
end

@inline
function _f4(s, epsilon, kappa)
    ω = epsilon / 2
    rp = 1 + kappa
    return (2.0)^(1 + 2 * s) * kappa^(2 * s - 2im * ω + 1) * exp(- 1im * ω * rp)
end

@inline
function B_trans(s, epsilon, tau, kappa)
    return _f1(s, epsilon, tau, kappa)
end

@inline
function B_inc(x, P1, P2, s, epsilon, tau, kappa)
    return C_inf(x, P1, P2, s, epsilon, tau, kappa) * _f3(epsilon, kappa)
end

@inline
function B_ref(x, P1, P2, s, epsilon, tau, kappa)
    return C_inf(x, P1, P2, s, epsilon, tau, kappa) * _f4(s, epsilon, kappa)
end

@inline
function B_inc(CInfIn, epsilon, tau, kappa)
    return _scaled(CInfIn, _f3(epsilon, kappa))
end

@inline
function B_ref(CInfOut, s, epsilon, tau, kappa)
    return _scaled(CInfOut, _f4(s, epsilon, kappa))
end

@inline
function C_trans(s, epsilon, kappa)
    return _f4(s, epsilon, kappa)
end

@inline
function C_ref(x, P1, P2, s, epsilon, tau, kappa)
    return C_zero(x, P1, P2, s, epsilon, tau, kappa) * _f1(s, epsilon, tau, kappa)
end

@inline
function C_inc(x, P1, P2, s, epsilon, tau, kappa)
    return C_zero(x, P1, P2, s, epsilon, tau, kappa) * _f2(epsilon, tau, kappa)
end

@inline
function C_ref(CZeroIn, s, epsilon, tau, kappa)
    return _scaled(CZeroIn, _f1(s, epsilon, tau, kappa))
end

@inline
function C_inc(CZeroOut, epsilon, tau, kappa)
    return _scaled(CZeroOut, _f2(epsilon, tau, kappa))
end

end
