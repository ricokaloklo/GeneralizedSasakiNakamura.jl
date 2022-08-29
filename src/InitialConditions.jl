module InitialConditions

using ..Kerr
using ..Coordinates

export Xinf_initialconditions, Xhor_initialconditions

const I = 1im # Mathematica being Mathematica

function Xinf_initialconditions(s::Int, m::Int, a, omega, lambda, rsout; order::Int=-1)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    _highest_order_implemented = 3
    order = (order == -1 ? _highest_order_implemented : order)
    A1 = 0.0
    A2 = 0.0
    A3 = 0.0

    if s == 0
        # s = 0
        if order >= 1
            A1 = begin
                (1/2)*I*(lambda + 2*a*m*omega)
            end
        end
        if order >= 2
            A2 = begin
                (1/8)*(-lambda^2 + lambda*(2 - 4*a*m*omega) + 4*omega*(I + a*m + a*m*(2*I - a*m)*omega))
            end
        end
        if order >= 3
            A3 = begin
                (1/2)*(1 - I*a*m)*omega - (1/48)*I*((-6 + lambda)*(-2 + lambda)*lambda + 2*lambda*(-18*I + a*m*(-16 + 3*lambda))*omega + 4*a*(3*a*m^2*(-2 + lambda) + 2*a*(-1 + lambda) - 6*I*m*(1 + lambda))*omega^2 + 8*a*m*(-8 - 6*I*a*m + a^2*(2 + m^2))*omega^3)
            end
        end
    elseif s == 1
        # s = +1
        if order >= 1
            A1 = begin
                (1/2)*I*(2 + lambda + 2*a*m*omega)
            end
        end
        if order >= 2
            A2 = begin
                (1/8)*((-lambda)*(2 + lambda) - 4*a*m*(1 + lambda)*omega - 4*a*m*(-2*I + a*m)*omega^2)
            end
        end
        if order >= 3
            A3 = begin
                (-(1/48))*I*((-4 + lambda)*lambda*(2 + lambda) + 2*(-12*I*(2 + lambda) + a*m*(-12 + lambda*(-4 + 3*lambda)))*omega + 4*a*(-6*I*m*(2 + lambda) + a*(6 + (2 + 3*m^2)*lambda))*omega^2 + 8*a*m*(-8 - 6*I*a*m + a^2*(2 + m^2))*omega^3)
            end
        end
    elseif s == -1
        # s = -1
        if order >= 1
            A1 = begin
                (I*(lambda^2 + 2*a*m*(2 + lambda)*omega))/(2*lambda)
            end
        end
        if order >= 2
            A2 = begin
                -(((-2 + lambda)*lambda^2 + 4*a*m*(-2 + lambda + lambda^2)*omega + 4*a*(-2*I*m*lambda + a*(2 + m^2*(4 + lambda)))*omega^2)/(8*lambda))
            end
        end
        if order >= 3
            A3 = begin
                -((I*((-6 + lambda)*(-2 + lambda)*lambda^2 + 48*a*m*omega + 2*lambda*(-12*I*lambda + a*m*(-16 + lambda*(-10 + 3*lambda)))*omega + 4*a*(-6*I*m*(4 + lambda^2) + a*(-12 + 3*m^2*(-2 + lambda)*(4 + lambda) + 2*lambda*(4 + lambda)))*omega^2 + 
                8*a*(-8*m*lambda - 6*I*a*(-2 + m^2*(2 + lambda)) + a^2*m*(m^2*(6 + lambda) + 2*(9 + lambda)))*omega^3))/(48*lambda))
            end
        end
    elseif s == 2
        # s = +2
        if order >= 1
            A1 = begin
                (1/2)*I*(6 + lambda + 2*a*m*omega)
            end
        end
        if order >= 2
            A2 = begin
                (1/8)*(-lambda^2 - 2*lambda*(5 + 2*a*m*omega) - 4*(6 + (3*I + 5*a*m)*omega + a*m*(-2*I + a*m)*omega^2))
            end
        end
        if order >= 3
            A3 = begin
                (-(1/48))*I*(lambda^3 + 2*lambda^2*(5 + 3*a*m*omega) + 4*lambda*(6 + (3*I + 10*a*m)*omega + a*(2*a - 6*I*m + 3*a*m^2)*omega^2) + 
                8*a*omega*(13*a*omega + 6*a*m^2*(1 - I*omega)*omega + a^2*m^3*omega^2 + m*(2 - 9*I*omega + 2*(-4 + a^2)*omega^2)))
            end
        end
    elseif s == -2
        # s = -2
        if order >= 1
            A1 = begin
                (1/6)*I*(-6 + 7*lambda) + I*a*m*omega + (2*I*(-3 + lambda)*lambda*(2 + lambda) + 24*(-3 - 3*I*a*m + lambda)*omega)/(-3*lambda*(2 + lambda) + 36*omega*(I - a*m + a^2*omega))
            end
        end
        if order >= 2
            A2 = begin
                -((lambda^2*(2 + lambda)^2 + 4*a*m*lambda*(16 + lambda*(14 + lambda))*omega + 4*(36 + a*(a*(10 - 11*lambda)*lambda - 2*I*m*(2 + lambda)*(6 + lambda) + a*m^2*(60 + lambda*(30 + lambda))))*omega^2 + 16*a*(-6*m + a*(3*I - 9*I*m^2 + a*m*(-15 + 3*m^2 - 7*lambda)))*omega^3 - 
                48*a^3*(-2*I*m + a*(-4 + m^2))*omega^4)/(8*lambda*(2 + lambda) - 96*omega*(I - a*m + a^2*omega)))
            end
        end
        if order >= 3
            A3 = begin
                -((1/(48*(lambda*(2 + lambda) - 12*omega*(I - a*m + a^2*omega))))*(I*(-4 + lambda)*lambda^2*(2 + lambda)^2 + 2*I*a*m*lambda*(-96 + lambda*(-44 + lambda*(32 + 3*lambda)))*omega + 
                4*I*(36*(-4 + lambda) + a*(-6*I*m*lambda*(2 + lambda)^2 + a*lambda*(-60 + (40 - 13*lambda)*lambda) + 3*a*m^2*(-48 + lambda*(40 + lambda*(24 + lambda)))))*omega^2 + 
                8*I*a*(-6*I*a*(-3*(2 + lambda) + m^2*(1 + lambda)*(18 + lambda)) - 4*m*(-9 + lambda*(13 + 2*lambda)) + a^2*m*(108 - lambda*(44 + 31*lambda) + m^2*(144 + lambda*(44 + lambda))))*omega^3 - 
                48*a*(16*m + 28*I*a*m^2 - 2*a^2*m*(5 + 7*m^2 - 7*lambda) - I*a^3*(2*m^4 + 2*(-9 + 5*lambda) - m^2*(32 + 7*lambda)))*omega^4 - 96*I*a^3*m*(-8 - 6*I*a*m + a^2*(-10 + m^2))*omega^5))
            end
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end

    rout = r_from_rstar(a, rsout)
    fansatz = 1 + A1/(omega*rout) + A2/(omega*rout)^2 + A3/(omega*rout)^3
    dfansatz_dr = -A1/(omega*rout^2) - (2*A2)/(omega^2 * rout^3) - (3*A3)/(omega^3 * rout^4)
    phase = exp(1im * omega * rsout)

    return phase*fansatz, phase*(1im*omega*fansatz + (Delta(a, rout)/(rout^2 + a^2))*dfansatz_dr)
end

function Xhor_initialconditions(s::Int, m::Int, a, omega, lambda, rsin; order::Int=-1)
    #=
    For Xhor, which we obtain by integrating from r_* -> -inf (or r -> r_+),
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