module InitialConditions

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
    _highest_order_implemented = 3
    order = (order == -1 ? _highest_order_implemented : order)
    # The following can be implemented in a more elegant (and general) way
    # But whatever, leave this to future me
    A1 = 0.0
    A2 = 0.0
    A3 = 0.0

    if order >= 1
        A1 = outgoing_coefficient_at_inf(s, m, a, omega, lambda, 1)
    end
    if order >= 2
        A2 = outgoing_coefficient_at_inf(s, m, a, omega, lambda, 2)
    end
    if order >= 3
        A3 = outgoing_coefficient_at_inf(s, m, a, omega, lambda, 3)
    end

    rout = r_from_rstar(a, rsout)
    fansatz = 1 + A1/(omega*rout) + A2/(omega*rout)^2 + A3/(omega*rout)^3
    dfansatz_dr = -A1/(omega*rout^2) - (2*A2)/(omega^2 * rout^3) - (3*A3)/(omega^3 * rout^4)
    phase = exp(1im * omega * rsout)

    return phase*fansatz, phase*(1im*omega*fansatz + (Delta(a, rout)/(rout^2 + a^2))*dfansatz_dr)
end

function Xin_initialconditions(s::Int, m::Int, a, omega, lambda, rsin; order::Int=-1)
    #=
    For Xin, which we obtain by integrating from r_* -> -inf (or r -> r_+),
    we simply start our integration very close to the EH
    i.e. gansatz = 1
    =#
    _highest_order_implemented = 0
    order = (order == -1 ? _highest_order_implemented : order)
    gansatz = 1
    dgansatz_dr = 0
    p = omega - m*omega_horizon(a)
    phase = exp(-1im * p * rsin)

    return phase*gansatz, phase*(-1im*p*gansatz)
end

end