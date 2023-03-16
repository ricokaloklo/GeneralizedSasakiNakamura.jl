module InitialConditions

using ForwardDiff
using ..Kerr
using ..Coordinates
using ..AsymptoticExpansionCoefficients

export Xup_initialconditions, Xin_initialconditions

const I = 1im # Mathematica being Mathematica

function Xup_initialconditions(s::Int, m::Int, a, omega, lambda, rsout; order::Int=-1)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    _default_order = 3
    order = (order == -1 ? _default_order : order)

    coeffs = zeros(ComplexF64, order+1)
    for i in 0:order
        coeffs[i+1] = outgoing_coefficient_at_inf(s, m, a, omega, lambda, i)
    end

    rout = r_from_rstar(a, rsout)
    function fansatz(r)
        ans = 0.0
        for i in 0:order
            ans += coeffs[i+1]/((omega*r)^i)
        end
        return ans
    end
    dfansatz_dr(r) = ForwardDiff.derivative(fansatz, r)
    _fansatz = fansatz(rout)
    _dfansatz_dr = dfansatz_dr(rout)
    phase = exp(1im * omega * rsout)

    return phase*_fansatz, phase*(1im*omega*_fansatz + (Delta(a, rout)/(rout^2 + a^2))*_dfansatz_dr)
end

function Xin_initialconditions(s::Int, m::Int, a, omega, lambda, rsin; order::Int=-1)
    #=
    For Xin, which we obtain by integrating from r_* -> -inf (or r -> r_+),

    Write Xin = \sum_j C^{H}_{-} (r - r_+)^j

    =#
    _default_order = 1
    order = (order == -1 ? _default_order : order)

    coeffs = zeros(ComplexF64, order+1)
    for i in 0:order
        coeffs[i+1] = ingoing_coefficient_at_hor(s, m, a, omega, lambda, i)
    end

    rin = r_from_rstar(a, rsin)
    function gansatz(r)
        ans = 0.0
        for i in 0:order
            ans += coeffs[i+1]*(r-r_plus(a))^i
        end
        return ans
    end
    dgansatz_dr(r) = ForwardDiff.derivative(gansatz, r)
    _gansatz = gansatz(rin)
    _dgansatz_dr = dgansatz_dr(rin)
    p = omega - m*omega_horizon(a)
    phase = exp(-1im * p * rsin)

    return phase*_gansatz, phase*(-1im*p*_gansatz + (Delta(a, rin)/(rin^2 + a^2))*_dgansatz_dr)
end

end