module AsymptoticExpansionCoefficients

using TaylorSeries
using ..Kerr
using ..Coordinates

export outgoing_coefficient_at_inf, ingoing_coefficient_at_inf
export outgoing_coefficient_at_hor, ingoing_coefficient_at_hor

const I = 1im # Mathematica being Mathematica
_DEFAULTDATATYPE = ComplexF64 # Double precision by default

function PminusInf_z(s::Int, m::Int, a, omega, lambda, z)
    if s == 0
        return begin
            -((2*I*(omega + a^2*z^4*(-I + a^2*omega) + z^2*(I + 2*a^2*omega)))/
            ((1 + a^2*z^2)*(1 - 2*z + a^2*z^2)))
        end
    elseif s == 1
        return begin
            ((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - (2*z*(1 - 2*z + a^2*z^2))/
            (1 + a^2*z^2)^2 + (2*a*z^2*(1 - 2*z + a^2*z^2)*((-I)*m*(1 + 3*a^2*z^2) + 
            a*z*(3 + 3*z + 2*lambda + 2*a^2*z^2*(1 + lambda))))/((1 + a^2*z^2)*
            (2 - 2*I*a*m*z - 2*I*a^3*m*z^3 + lambda + a^4*z^4*(1 + lambda) + 
            a^2*z^2*(3 + 2*z + 2*lambda))) - 2*I*omega))/(1 - 2*z + a^2*z^2)
        end
    elseif s == -1
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 + 
            (2*a*z^2*(1 - 2*z + a^2*z^2)*(I*m*(1 + 3*a^2*z^2) + 
            a*z*(-1 + 3*z + 2*a^2*z^2*(-1 + lambda) + 2*lambda)))/
            ((1 + a^2*z^2)*(2*I*a*m*z + 2*I*a^3*m*z^3 + a^4*z^4*(-1 + lambda) + lambda + 
            a^2*z^2*(-1 + 2*z + 2*lambda))) - 2*I*omega))
        end
    elseif s == 2
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 - 2*I*omega + 
            (8*a*z^2*(1 - 2*z + a^2*z^2)*(-6*a*m^2*z + I*m*(-4 + 9*a^2*z^2 - lambda + 
            z*(-6 - 12*I*a^2*omega)) + a*(6*a^2*z^3 + I*(1 + lambda)*omega + z^2*(-9 - 9*I*a^2*omega) + 
            z*(3 + 6*I*omega - 6*a^2*omega^2))))/((1 + a^2*z^2)*(24 + 10*lambda + lambda^2 - 
            4*I*a*m*(6*z^2 + 2*z*(4 + lambda) + 3*I*omega) + 12*I*omega + 24*a^3*m*z^2*(I*z + 2*omega) + 
            12*a^4*z^2*(z^2 - 2*I*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*z*(1 + lambda)*omega + 3*omega^2)))))
        end
    elseif s == -2
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 - 2*I*omega + 
            (8*a*z^2*(1 - 2*z + a^2*z^2)*(I*m*(6*z + lambda) + 3*a^2*m*z*(-3*I*z + 4*omega) - 
            a*(9*z^2 + z*(-3 + 6*m^2 + 6*I*omega) + I*(-3 + lambda)*omega) + 3*a^3*z*(2*z^2 + 3*I*z*omega - 2*omega^2)))/
            ((1 + a^2*z^2)*(2*lambda + lambda^2 - 12*I*omega + 24*a^3*m*z^2*((-I)*z + 2*omega) + 4*a*m*(6*I*z^2 + 2*I*z*lambda + 3*omega) + 
            12*a^4*z^2*(z^2 + 2*I*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 + 6*I*omega) + 2*I*z*(-3 + lambda)*omega + 3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function QminusInf_z(s::Int, m::Int, a, omega, lambda, z)
    if s == 0
        return begin
            -((z^2*(-2*z*(-1 + lambda) + lambda - 4*a^2*z^3*(1 + lambda) - 2*a^4*z^5*(3 + lambda) + 2*a*m*omega + 
            a^6*z^6*(1 - m^2 + lambda + 2*a*m*omega) + a^2*z^4*(8 + a^2*(2 - 2*m^2 + 3*lambda) + 
            6*a^3*m*omega) + z^2*(-4 + a^2*(1 - m^2 + 3*lambda) + 6*a^3*m*omega)))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2))
        end
    elseif s == 1
        return begin
            (z^2*(-4 - lambda^2 + 2*a^12*z^12*(1 + lambda) - 8*a*m*omega - 4*I*a^10*z^11*(1 + lambda)*
            (-3*I + a^2*omega) - 2*lambda*(2 + a*m*omega) - a^4*z^6*(52 + 32*lambda + 2*I*a*m*(31 + 6*lambda) + 
            a^2*(71 + 10*lambda^2 - 3*m^2*(3 + 2*lambda) + 8*lambda*(7 - 6*I*omega)) + 
            4*a^3*m*(17 + 5*lambda)*omega) + 2*a^6*z^9*(8 - 6*I*a*m + 
            a^2*(-3 + 3*m^2 - 4*lambda + lambda^2 + 4*I*omega) + 2*I*a^4*(m^2 - 5*(1 + lambda))*omega + 
            a^3*m*(4*I - I*m^2 + I*lambda + 2*omega)) + 
            2*z*((2 + lambda)^2 + 2*I*a^2*(-1 + m^2 - lambda)*omega + a*m*(2*I + I*lambda + 4*omega)) + 
            2*a^3*z^5*(14*I*m + a*(69 - m^2 + 48*lambda + 6*lambda^2 + 28*I*omega) + 
            4*I*a^3*(3*m^2 - 5*(1 + lambda))*omega + 3*a^2*m*(6*I - I*m^2 + 2*I*lambda + 6*omega)) + 
            2*a^4*z^7*(-4 + 8*I*a*m + 4*a^2*(10 + m^2 + 7*lambda + lambda^2 + 5*I*omega) + 
            4*I*a^4*(2*m^2 - 5*(1 + lambda))*omega + a^3*m*(14*I - 3*I*m^2 + 4*I*lambda + 10*omega)) + 
            2*a^2*z^3*(36 - 2*m^2 + 26*lambda + 4*lambda^2 + 12*I*omega + 2*I*a^2*(4*m^2 - 5*(1 + lambda))*omega + 
            a*m*(10*I - I*m^2 + 4*I*lambda + 14*omega)) - 
            a^8*z^10*(-6*I*a*m - 8*(1 + 2*lambda) + 2*a^3*m*(3 + lambda)*omega + 
            a^2*(1 - m^2*(-1 + lambda) + lambda^2 - 4*I*omega - 8*I*lambda*omega)) - 
            a*z^2*(4*I*m*(2 + lambda) + 2*a^2*m*(19 + 5*lambda)*omega + 
            a*(28 + 24*lambda + 5*lambda^2 - m^2*(4 + lambda) + 4*I*omega - 8*I*lambda*omega)) - 
            a^2*z^4*(24*(2 + lambda) + 2*I*a*m*(23 + 6*lambda) + 4*a^3*m*(18 + 5*lambda)*omega + 
            a^2*(67 + 54*lambda + 10*lambda^2 - m^2*(11 + 4*lambda) + 8*I*omega - 32*I*lambda*omega)) + 
            a^6*z^8*(4 + 8*lambda - 2*I*a*m*(9 + 2*lambda) - 2*a^3*m*(16 + 5*lambda)*omega + 
            a^2*(-31 - 24*lambda - 5*lambda^2 + m^2*(1 + 4*lambda) + 8*I*omega + 32*I*lambda*omega))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(2 - 2*I*a*m*z - 2*I*a^3*m*z^3 + lambda + 
            a^4*z^4*(1 + lambda) + a^2*z^2*(3 + 2*z + 2*lambda)))
        end
    elseif s == -1
        return begin
            (z^2*((-1 + 2*z)*lambda^2 + a^2*z*(m^2*(-4*z^2 + z*(2 + lambda) - 4*I*omega) + 
            lambda*(-24*z^3 + 4*z^2*(5 + 2*lambda) + z*(-4 - 5*lambda + 8*I*omega) - 4*I*omega)) + 
            2*I*a*m*lambda*(-z + 2*z^2 + I*omega) + 2*a^12*z^11*(-1 + lambda)*(z - 2*I*omega) - 
            2*a^11*m*z^10*(-5 + lambda)*omega + 2*a^9*m*z^8*(-3*I*z^2 + I*z*(-2 + m^2 - lambda + 10*I*omega) - 
            5*(-4 + lambda)*omega) + 2*I*a^7*m*z^6*(6*z^3 + z^2*(5 + 2*lambda) + 
            z*(-6 + 3*m^2 - 4*lambda + 30*I*omega) + 10*I*(-3 + lambda)*omega) + 
            2*I*a^5*m*z^4*(-8*z^3 + z^2*(19 + 6*lambda) + 3*z*(-2 + m^2 - 2*lambda + 10*I*omega) + 
            10*I*(-2 + lambda)*omega) + 2*I*a^3*m*z^2*(-14*z^3 + z^2*(11 + 6*lambda) + 
            z*(-2 + m^2 - 4*lambda + 10*I*omega) + 5*I*(-1 + lambda)*omega) - 
            a^10*z^9*(12*z^2*(-1 + lambda) + 4*I*(-4 + m^2 + 5*lambda)*omega + 
            z*(5 - m^2*(-3 + lambda) - 4*lambda + lambda^2 + 16*I*omega - 8*I*lambda*omega)) + 
            a^4*z^3*(-8*z^4 - 4*z^3*(-3 + 8*lambda) - 2*z^2*(3 + m^2 - 24*lambda - 6*lambda^2 - 8*I*omega) - 
            4*I*(-1 + 4*m^2 + 5*lambda)*omega + z*(1 - 14*lambda - 10*lambda^2 + m^2*(3 + 4*lambda) - 16*I*omega + 
            32*I*lambda*omega)) + a^8*z^7*(8*z^3*(-3 + 2*lambda) + 
            2*z^2*(9 + 3*m^2 - 8*lambda + lambda^2 + 8*I*omega) - 8*I*(-3 + 2*m^2 + 5*lambda)*omega + 
            z*(-3 - 4*lambda - 5*lambda^2 + m^2*(-7 + 4*lambda) - 48*I*omega + 32*I*lambda*omega)) + 
            a^6*z^5*(16*z^4 + 4*z^3*(-3 + 2*lambda) + 8*z^2*(m^2 + 3*lambda + lambda^2 + 4*I*omega) - 
            8*I*(-2 + 3*m^2 + 5*lambda)*omega + z*(1 - 16*lambda - 10*lambda^2 + m^2*(-3 + 6*lambda) - 48*I*omega + 
            48*I*lambda*omega))))/((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*
            (2*I*a*m*z + 2*I*a^3*m*z^3 + a^4*z^4*(-1 + lambda) + lambda + a^2*z^2*(-1 + 2*z + 2*lambda)))
        end
    elseif s == 2
        return begin
            (z^2*(24*a^12*z^12 - lambda^3 - 48*I*a^10*z^11*(-3*I + a^2*omega) - 2*lambda^2*(8 + a*m*omega) + 
            4*lambda*(-21 - 3*(I + 4*a*m)*omega + 7*a^2*omega^2) - 
            12*a^8*z^10*(-24 + 6*I*a*m - a^2*(5 + 3*m^2 + 20*I*omega) + 4*a^4*omega^2) + 
            8*I*a^6*z^9*(24*I + 18*a*m + 3*I*a^2*(9 + m^2 + 12*I*omega) + 
            a^3*m*(19 - 3*m^2 + 4*lambda + 18*I*omega) + a^4*(-31 + 15*m^2 - 4*lambda - 24*I*omega)*omega - 
            24*a^5*m*omega^2 + 12*a^6*omega^3) + 8*(-18 - (9*I + 23*a*m)*omega + 
            a*(-3*I*m + a*(11 - 3*m^2))*omega^2 + 3*a^3*m*omega^3) + 
            a^3*z^6*(384*I*m + 8*I*a^2*m*(-177 + 18*m^2 - 53*lambda - 2*lambda^2 + 24*I*omega) + 
            16*a*(-33 + 12*m^2 - 10*lambda - lambda^2 - 36*I*omega) + 
            a^3*(-756 - 48*m^4 - lambda^3 + lambda^2*(-41 + 16*I*omega) + lambda*(-334 + 364*I*omega) + 
            m^2*(564 + 106*lambda + lambda^2 - 756*I*omega) + 780*I*omega) + 
            2*a^4*m*(-506 + 126*m^2 - 132*lambda - lambda^2 + 588*I*omega)*omega + 
            4*a^5*(49 - 93*m^2 + 37*lambda - 144*I*omega)*omega^2 + 168*a^6*m*omega^3) + 
            a*z^4*(96*I*m + 8*I*a^2*m*(-202 + 9*m^2 - 65*lambda - 4*lambda^2) - 
            8*a*(120 + 50*lambda + 5*lambda^2 + 72*I*omega) + a^3*(-1068 - 24*m^4 - 3*lambda^3 + 
            lambda^2*(-74 + 32*I*omega) + 2*m^2*(246 + 58*lambda + lambda^2 - 288*I*omega) + 
            lambda*(-512 + 388*I*omega) + 488*I*omega) + 6*a^4*m*(-224 + 36*m^2 - 60*lambda - lambda^2 + 
            172*I*omega)*omega + 12*a^5*(38 - 34*m^2 + 17*lambda - 48*I*omega)*omega^2 + 216*a^6*m*omega^3) + 
            4*I*a*z^3*(12*m*(3 + lambda) + I*a*(-408 - 194*lambda - 27*lambda^2 - lambda^3 + 16*m^2*(4 + lambda) - 
            168*I*omega) - 12*a^4*m*(17 + lambda)*omega^2 + 96*a^5*omega^3 + 
            2*a^3*omega*(-53 - 26*lambda - 3*lambda^2 + m^2*(79 + 7*lambda) + 60*I*omega + 24*I*lambda*omega) - 
            2*a^2*m*(-98 - 35*lambda - 3*lambda^2 + m^2*(10 + lambda) + 128*I*omega + 32*I*lambda*omega)) + 
            2*I*a^3*z^5*(24*m*(9 + 4*lambda) - 24*a^4*m*(25 + lambda)*omega^2 + 288*a^5*omega^3 + 
            4*a^3*omega*(-89 - 35*lambda - 3*lambda^2 + m^2*(113 + 8*lambda) + 51*I*omega + 27*I*lambda*omega) - 
            4*a^2*m*(-143 - 44*lambda - 3*lambda^2 + m^2*(23 + 2*lambda) + 181*I*omega + 43*I*lambda*omega) + 
            I*a*(-924 - 394*lambda - 47*lambda^2 - lambda^3 + 4*m^2*(91 + 16*lambda) - 300*I*omega + 
            84*I*lambda*omega)) + 4*a^6*z^8*(36 - 24*m^2 + 6*a^3*m*(4*m^2 - 3*(4 + lambda - 6*I*omega))*
            omega - 2*a^4*(7 + 15*m^2 - 5*lambda + 24*I*omega)*omega^2 + 12*a^5*m*omega^3 + 
            2*a*m*(-44*I + 9*I*m^2 - 11*I*lambda + 36*omega) - 
            a^2*(39 + 6*m^4 + 20*lambda + 2*lambda^2 + m^2*(-65 - 8*lambda + 78*I*omega) - 142*I*omega - 
            22*I*lambda*omega + 48*omega^2)) + 8*I*a^4*z^7*(-12*I + 6*a*m*(-5 + lambda) + 
            I*a^2*(m^2*(62 + 8*lambda) - 3*(24 + lambda^2 + 2*lambda*(5 - I*omega) + 10*I*omega)) + 
            a^3*m*(88 + lambda^2 - m^2*(16 + lambda) + lambda*(23 - 18*I*omega) - 66*I*omega) - 
            2*a^5*m*(49 + lambda)*omega^2 + 48*a^6*omega^3 + a^4*omega*(-73 - lambda^2 + 3*m^2*(23 + lambda) - 
            14*I*omega + 10*I*lambda*(2*I + omega))) + 2*z*((9 + lambda)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            8*I*a^3*m*(13 + lambda)*omega^2 + 48*I*a^4*omega^3 - 
            4*a^2*omega*(12*I + I*lambda^2 - 2*I*m^2*(10 + lambda) + 19*omega + 7*lambda*(I + omega)) + 
            4*a*m*(24*I + I*lambda^2 + 31*omega + lambda*(10*I + 7*omega))) + 
            z^2*(-8*I*a*m*(60 + 23*lambda + 2*lambda^2 + 6*I*omega) - 12*(24 + 10*lambda + lambda^2 + 12*I*omega) + 
            6*a^3*m*(-134 + 10*m^2 - 36*lambda - lambda^2 + 44*I*omega)*omega + 
            4*a^4*(85 - 45*m^2 + 31*lambda - 48*I*omega)*omega^2 + 120*a^5*m*omega^3 + 
            a^2*(-3*lambda^3 + lambda^2*(-57 + 16*I*omega) + lambda*(-342 + 100*I*omega) + 
            m^2*(152 + 42*lambda + lambda^2 - 132*I*omega) - 12*(54 + 3*I*omega + 4*omega^2)))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(24 + 10*lambda + lambda^2 - 
            4*I*a*m*(6*z^2 + 2*z*(4 + lambda) + 3*I*omega) + 12*I*omega + 24*a^3*m*z^2*(I*z + 2*omega) + 
            12*a^4*z^2*(z^2 - 2*I*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*z*(1 + lambda)*omega + 3*omega^2)))
        end
    elseif s == -2
        return begin
            (z^2*(-((-1 + 2*z)*(-2 + 6*z - lambda)*(2*lambda + lambda^2 - 12*I*omega)) - 48*a^11*m*z^8*omega*(3*z^2 - omega^2) + 
            24*a^12*z^10*(z^2 - 2*I*z*omega + 4*omega^2) - 4*a^10*z^8*(36*z^3 - 3*z^2*(5 + 3*m^2 + 8*I*omega) + 6*(-17 + 5*m^2 - lambda)*omega^2 + 
            2*z*omega*(39*I + 3*I*m^2 - 4*I*lambda + 24*omega)) + 8*a^9*m*z^6*(9*I*z^4 + z^2*(-72 + 12*m^2 - 7*lambda - 6*I*omega)*omega + 
            2*I*z*(-3 + lambda)*omega^2 + 21*omega^3 + z^3*(-3*I + 3*I*m^2 - 4*I*lambda + 30*omega)) - 
            2*I*a^7*m*z^4*(72*z^5 + 4*z^4*(9*m^2 - 11*lambda + 12*I*omega) + I*z^2*(-426 + 126*m^2 - 92*lambda - lambda^2 - 12*I*omega)*omega - 
            24*z*(-3 + lambda)*omega^2 + 108*I*omega^3 + 4*z^3*(12 + 15*lambda + lambda^2 - m^2*(12 + lambda) + 150*I*omega + 14*I*lambda*omega)) - 
            2*I*a^3*m*(192*z^6 + 24*z^5*(-7 + 4*lambda) + 4*z^4*(-6 + 9*m^2 - 33*lambda - 4*lambda^2 - 72*I*omega) - 8*z*(-3 + lambda)*omega^2 + 12*I*omega^3 - 
            z^2*omega*(90*I - 30*I*m^2 + 52*I*lambda + 3*I*lambda^2 + 60*omega) + 4*z^3*(6 + 11*lambda + 3*lambda^2 - m^2*(6 + lambda) + 84*I*omega + 
            20*I*lambda*omega)) - 2*I*a^5*m*z^2*(24*z^5*(-9 + lambda) + 4*z^4*(3 + 18*m^2 - 37*lambda - 2*lambda^2 - 96*I*omega) + 
            3*I*z^2*(-96 + 36*m^2 - 36*lambda - lambda^2 + 20*I*omega)*omega - 24*z*(-3 + lambda)*omega^2 + 60*I*omega^3 + 
            4*z^3*(15 + 20*lambda + 3*lambda^2 - m^2*(15 + 2*lambda) + 189*I*omega + 31*I*lambda*omega)) + 
            4*a^8*z^6*(-6*(9 + m^2)*z^3 + 72*z^4 + z^2*(9 - 6*m^4 - 2*lambda^2 + lambda*(-4 - 22*I*omega) + m^2*(33 + 8*lambda + 30*I*omega) + 258*I*omega) + 
            3*(55 - 31*m^2 + 7*lambda)*omega^2 - 2*z*omega*(87*I - I*lambda^2 + 3*I*m^2*(3 + lambda) + 114*omega + 6*lambda*(-2*I + omega))) + 
            2*a*m*(-48*I*z^4 - 24*I*z^3*(-1 + lambda) + 4*I*z^2*(7*lambda + 2*lambda^2 + 18*I*omega) - (12 + 8*lambda + lambda^2 - 12*I*omega)*omega + 
            4*z*((-I)*lambda^2 + 15*omega + lambda*(-2*I + 3*omega))) - a^2*(8*z^4*(10*lambda + 5*lambda^2 - 72*I*omega) - 
            4*z^3*((26 - 16*m^2)*lambda + 15*lambda^2 + lambda^3 - 240*I*omega) + 12*(-2 + 2*m^2 - lambda)*omega^2 + 
            z^2*(3*lambda^3 + lambda^2*(21 - m^2 + 16*I*omega) + lambda*(30 - 34*m^2 - 28*I*omega) + 12*I*(-35 + 5*m^2 + 12*I*omega)*omega) + 
            8*z*omega*(6*I - I*lambda^2 + 15*omega + lambda*(I + 2*I*m^2 + 3*omega))) + 
            a^4*z^2*(96*z^5 + 16*z^4*(-9 + 12*m^2 - 2*lambda - lambda^2 + 36*I*omega) + 60*(3 - 3*m^2 + lambda)*omega^2 + 
            2*z^3*(36 + 66*lambda + 35*lambda^2 + lambda^3 - 4*m^2*(27 + 16*lambda) - 1068*I*omega + 84*I*lambda*omega) - 
            8*z*omega*(39*I - 3*I*lambda^2 + I*m^2*(3 + 7*lambda) + 84*omega + 2*lambda*(-I + 6*omega)) - 
            z^2*(24*m^4 + 3*lambda^3 - 2*m^2*(30 + 50*lambda + lambda^2) + lambda^2*(38 + 32*I*omega) + 4*lambda*(16 + 33*I*omega) - 
            12*(-1 + 134*I*omega + 48*omega^2))) - a^6*z^4*(48*(-3 + 2*m^2)*z^4 + 192*z^5 + 
            8*z^3*(m^2*(30 + 8*lambda) - 3*(lambda^2 + lambda*(2 + 2*I*omega) - 54*I*omega)) + 12*(-42 + 34*m^2 - 9*lambda)*omega^2 + 
            8*z*omega*(87*I - 3*I*lambda^2 + I*m^2*(9 + 8*lambda) + 159*omega + lambda*(-11*I + 15*omega)) + 
            z^2*(48*m^4 + lambda^3 + lambda^2*(29 + 16*I*omega) - m^2*(156 + 98*lambda + lambda^2 + 180*I*omega) + lambda*(54 + 236*I*omega) - 
            12*(-1 + 179*I*omega + 48*omega^2)))))/((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*
            (2*lambda + lambda^2 - 12*I*omega + 24*a^3*m*z^2*((-I)*z + 2*omega) + 4*a*m*(6*I*z^2 + 2*I*z*lambda + 3*omega) + 
            12*a^4*z^2*(z^2 + 2*I*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 + 6*I*omega) + 2*I*z*(-3 + lambda)*omega + 3*omega^2)))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

# Cache mechanism for the ingoing coefficients at infinity
# Initialize the cache with a set of fiducial parameters
_cached_ingoing_coefficients_at_inf_params::NamedTuple{(:s, :m, :a, :omega, :lambda), Tuple{Int, Int, _DEFAULTDATATYPE, _DEFAULTDATATYPE, _DEFAULTDATATYPE}} = (s=-2, m=2, a=0, omega=0.5, lambda=1)
_cached_ingoing_coefficients_at_inf::NamedTuple{(:expansion_coeffs, :Pcoeffs, :Qcoeffs), Tuple{Vector{_DEFAULTDATATYPE}, Vector{_DEFAULTDATATYPE}, Vector{_DEFAULTDATATYPE}}} = (
    expansion_coeffs = [_DEFAULTDATATYPE(1.0)], 
    Pcoeffs = [_DEFAULTDATATYPE(0.0)],
    Qcoeffs = [_DEFAULTDATATYPE(0.0)]
)

function ingoing_coefficient_at_inf(s::Int, m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    for up to (1/r)^3, which is probably more than enough
    
    But we have also shown a recurrence relation where
    one can generate as higher an order as one pleases.
    However, the recurrence relation that we have
    actually depends *all* previous terms so this function
    is designed to be evaluated recursively to build
    the full list of coefficients
    =#
    global _cached_ingoing_coefficients_at_inf_params
    global _cached_ingoing_coefficients_at_inf

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    if order == 0
        return 1.0 # This is always 1.0
    elseif order == 1
        if s == 0
            return (-(1/2))*I*(lambda + 2*a*m*omega)
        elseif s == +1
            return begin
                -((I*(4 + lambda^2 + 8*a*m*omega + 2*lambda*(2 + a*m*omega)))/
                (2*(2 + lambda)))
            end
        elseif s == -1
            return (-(1/2))*I*(lambda + 2*a*m*omega)
        elseif s == +2
            return begin
                (I*(-lambda^3 - 2*lambda^2*(8 + a*m*omega) + 
                4*lambda*(-21 - 3*(I + 4*a*m)*omega + 7*a^2*omega^2) + 
                8*(-18 - (9*I + 23*a*m)*omega + 
                a*(11*a - 3*I*m - 3*a*m^2)*omega^2 + 
                3*a^3*m*omega^3)))/(2*(10*lambda + lambda^2 + 
                12*(2 + (I + a*m)*omega - a^2*omega^2)))
            end
        elseif s == -2
            return (-(1/2))*I*(2 + lambda + 2*a*m*omega)
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    elseif order == 2
        if s == 0
            return begin
                (1/8)*(-((-2 + lambda)*lambda) - 4*(I + a*m*(-1 + lambda))*omega - 
                4*a*m*(2*I + a*m)*omega^2)
            end
        elseif s == +1
            return begin
                -((lambda^3 + 4*lambda^2*(1 + a*m*omega) + 
                8*a*omega*(m*(2 + 2*I*omega) - a*omega + 3*a*m^2*omega) + 
                4*lambda*(1 + a*m*(5 + 2*I*omega)*omega + a^2*(-2 + m^2)*omega^2))/
                (8*(2 + lambda)))
            end
        elseif s == -1
            return begin
                (1/8)*(-lambda^2 + lambda*(2 - 4*a*m*omega) + 
                4*a*omega*(m + 2*a*omega - m*(2*I + a*m)*omega))
            end
        elseif s == +2
            return begin
                -((lambda^4 + 4*lambda^3*(5 + a*m*omega) + 
                4*lambda^2*(37 + 2*a*m*(13 + I*omega)*omega + 
                a^2*(-11 + m^2)*omega^2) - 
                8*lambda*(-60 + 8*a*m*(-11 - 2*I*omega)*omega + 
                a^2*(39 - 19*m^2)*omega^2 + 14*a^3*m*omega^3) - 
                16*(a^2*(34 + m^2*(-49 - 9*I*omega) + 3*I*omega)*
                omega^2 + a^3*m*(43 - 3*m^2 + 6*I*omega)*omega^3 + 
                3*a^4*(-4 + m^2)*omega^4 - 9*(4 + omega^2) + 
                2*a*m*omega*(-44 - 15*I*omega + 3*omega^2)))/
                (8*(10*lambda + lambda^2 + 12*(2 + (I + a*m)*omega - 
                a^2*omega^2))))
            end
        elseif s == -2
            return begin
                (1/8)*((-lambda)*(2 + lambda) - 4*(-3*I + a*m*(1 + lambda))*omega - 
                4*a*m*(2*I + a*m)*omega^2)
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    elseif order == 3
        if s == 0
            return begin
                (1/48)*(I*(-6 + lambda)*(-2 + lambda)*lambda + 
                2*(12 - 18*lambda + I*a*m*(12 + lambda*(-16 + 3*lambda)))*omega + 
                4*I*a*(3*a*m^2*(-2 + lambda) + 2*a*(-1 + lambda) + 
                6*I*m*(1 + lambda))*omega^2 + 
                8*I*a*m*(-8 + 6*I*a*m + a^2*(2 + m^2))*omega^3)
            end
        elseif s == +1
            return begin
                (1/(48*(2 + lambda)))*(I*(lambda^4 + 6*a*m*lambda^3*omega + 
                4*lambda^2*(-3 + (6*I + 4*a*m)*omega + 
                a*(-4*a + 6*I*m + 3*a*m^2)*omega^2) + 
                8*lambda*(-2 + (12*I - 5*a*m)*omega + 
                a*(12*I*m + a*(-4 + 9*m^2))*omega^2 + 
                a*m*(-8 + 6*I*a*m + a^2*(-4 + m^2))*omega^3) + 
                16*omega*(6*I + a^3*m*(-1 + 4*m^2)*omega^2 + 
                3*I*a^2*omega*(I - 2*omega + 4*m^2*omega) + 
                a*m*(-3 + 12*I*omega - 8*omega^2))))
            end
        elseif s == -1
            return begin
                (1/48)*I*(lambda^3 + lambda^2*(-8 + 6*a*m*omega) + 
                8*a*omega*(a*omega + m*(2 - 3*a*m*omega + 
                (-8 + 6*I*a*m + a^2*(-4 + m^2))*omega^2)) + 
                4*lambda*(3 + omega*(6*I + a*(-4*a*omega + 
                m*(-8 + 3*(2*I + a*m)*omega)))))
            end
        elseif s == +2
            return begin
                -((I*(-lambda^5 - 2*lambda^4*(10 + 3*a*m*omega) - 
                4*lambda^3*(37 + 2*a*m*(20 + 3*I*omega)*omega + 
                a^2*(-13 + 3*m^2)*omega^2) - 
                8*lambda^2*(60 + 2*a^2*(-29 + 3*m^2*(9 + I*omega))*
                omega^2 + a^3*m*(-31 + m^2)*omega^3 + 
                a*m*omega*(157 + 48*I*omega - 8*omega^2)) + 
                16*lambda*(a^2*(91 + m^2*(-210 - 81*I*omega) + 9*I*omega)*
                omega^2 + 2*a^3*m*(73 - 13*m^2 + 21*I*omega)*
                omega^3 + 3*a^4*(-10 + 7*m^2)*omega^4 - 
                9*(4 + omega^2) + 2*a*m*omega*(-116 - 63*I*omega + 
                29*omega^2)) + 96*a*omega*((-a^3)*m^4*omega^3 + 
                a*omega*(18 + 9*I*omega - 11*a^2*omega^2) + 
                a^2*m^3*omega^2*(-28 - 7*I*omega + a^2*omega^2) + 
                a*m^2*omega*(-70 - 55*I*omega + 2*(7 + 15*a^2)*
                omega^2 + 6*I*a^2*omega^3) + 
                m*(-36 - 36*I*omega + (25 + 47*a^2)*omega^2 + 
                I*(8 + 23*a^2)*omega^3 - 2*a^2*(4 + 5*a^2)*
                omega^4))))/(48*(10*lambda + lambda^2 + 
                12*(2 + (I + a*m)*omega - a^2*omega^2))))
            end
        elseif s == -2
            return begin
                (1/48)*(I*(-4 + lambda)*lambda*(2 + lambda) + 
                2*(6*(-4 + lambda) + I*a*m*(-24 + lambda*(-4 + 3*lambda)))*
                omega + 4*a*(-6*m*(-1 + lambda) + 
                I*a*(-6 + (2 + 3*m^2)*lambda))*omega^2 + 
                8*I*a*m*(-8 + 6*I*a*m + a^2*(2 + m^2))*omega^3)
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    else
        # Evaluate higher order corrections using AD
        # Specifically we use TaylorSeries.jl for a much more performant AD

        _this_params = (s=s, m=m, a=a, omega=omega, lambda=lambda)
        # Check if we can use the cached results
        if _cached_ingoing_coefficients_at_inf_params == _this_params
            expansion_coeffs = _cached_ingoing_coefficients_at_inf.expansion_coeffs
            Pcoeffs = _cached_ingoing_coefficients_at_inf.Pcoeffs
            Qcoeffs = _cached_ingoing_coefficients_at_inf.Qcoeffs
        else
            # Cannot re-use the cached results, re-compute from zero
            expansion_coeffs = [data_type(1.0)] # order 0
            Pcoeffs = [data_type(PminusInf_z(s, m, a, omega, lambda, 0))] # order 0
            Qcoeffs = [data_type(0.0), data_type(0.0)] # the recurrence relation takes Q_{r+1}
        end

        # Compute Pcoeffs to the necessary order
        _P(z) = PminusInf_z(s, m, a, omega, lambda, z)
        _P_taylor = taylor_expand(_P, 0, order=order) # FIXME This is not the most efficient way to do this
        for i in length(Pcoeffs):order
            append!(Pcoeffs, getcoeff(_P_taylor, i))
        end

        # Compute Qcoeffs to the necessary order (to current order + 1)
        _Q(z) = QminusInf_z(s, m, a, omega, lambda, z)
        _Q_taylor = taylor_expand(_Q, 0, order=order+1)
        for i in length(Qcoeffs):order+1
            append!(Qcoeffs, getcoeff(_Q_taylor, i))
        end

        # Note that the expansion coefficients we store is scaled by \omega^{i}
        for i in length(expansion_coeffs):order
            _P0 = Pcoeffs[1] # P0
            sum = 0.0
            for k in 1:i
                sum += (Qcoeffs[k+2] - (i-k)*Pcoeffs[k+1])*(expansion_coeffs[i-k+1]/omega^(i-k))
            end
            append!(expansion_coeffs, omega^(i)*((i*(i-1)*(expansion_coeffs[i]/omega^(i-1)) + sum)/(_P0*i)))
        end

        # Update cache
        _cached_ingoing_coefficients_at_inf_params = _this_params
        _cached_ingoing_coefficients_at_inf = (
            expansion_coeffs = expansion_coeffs,
            Pcoeffs = Pcoeffs,
            Qcoeffs = Qcoeffs
        )

        return expansion_coeffs[order+1]
    end
end

function PplusInf_z(s::Int, m::Int, a, omega, lambda, z)
    if s == 0
        return begin
            (2*I*(omega + a^2*z^4*(I + a^2*omega) + z^2*(-I + 2*a^2*omega)))/
            ((1 + a^2*z^2)*(1 - 2*z + a^2*z^2))
        end
    elseif s == 1
        return begin
            ((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - (2*z*(1 - 2*z + a^2*z^2))/
            (1 + a^2*z^2)^2 + (2*a*z^2*(1 - 2*z + a^2*z^2)*((-I)*m*(1 + 3*a^2*z^2) + 
            a*z*(3 + 3*z + 2*lambda + 2*a^2*z^2*(1 + lambda))))/((1 + a^2*z^2)*
            (2 - 2*I*a*m*z - 2*I*a^3*m*z^3 + lambda + a^4*z^4*(1 + lambda) + 
            a^2*z^2*(3 + 2*z + 2*lambda))) + 2*I*omega))/(1 - 2*z + a^2*z^2)
        end
    elseif s == -1
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 + 
            (2*a*z^2*(1 - 2*z + a^2*z^2)*(I*m*(1 + 3*a^2*z^2) + 
            a*z*(-1 + 3*z + 2*a^2*z^2*(-1 + lambda) + 2*lambda)))/
            ((1 + a^2*z^2)*(2*I*a*m*z + 2*I*a^3*m*z^3 + a^4*z^4*(-1 + lambda) + lambda + 
            a^2*z^2*(-1 + 2*z + 2*lambda))) + 2*I*omega))
        end
    elseif s == 2
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 + 2*I*omega + 
            (8*a*z^2*(1 - 2*z + a^2*z^2)*(-6*a*m^2*z + I*m*(-4 + 9*a^2*z^2 - lambda + 
            z*(-6 - 12*I*a^2*omega)) + a*(6*a^2*z^3 + I*(1 + lambda)*omega + z^2*(-9 - 9*I*a^2*omega) + 
            z*(3 + 6*I*omega - 6*a^2*omega^2))))/((1 + a^2*z^2)*(24 + 10*lambda + lambda^2 - 
            4*I*a*m*(6*z^2 + 2*z*(4 + lambda) + 3*I*omega) + 12*I*omega + 24*a^3*m*z^2*(I*z + 2*omega) + 
            12*a^4*z^2*(z^2 - 2*I*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*z*(1 + lambda)*omega + 3*omega^2)))))
        end
    elseif s == -2
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 + z*(-2 + a^2*z)))/(1 + a^2*z^2)^2 + 2*I*omega + 
            (8*a*z^2*(1 + z*(-2 + a^2*z))*(-3*a*z*(-1 + 2*m^2 + 3*z) + I*m*(6*z + lambda) - 
            I*a*(-3 + 6*z + lambda)*omega + 3*a^2*m*z*(-3*I*z + 4*omega) + 
            3*a^3*z*(2*z^2 + 3*I*z*omega - 2*omega^2)))/((1 + a^2*z^2)*(lambda*(2 + lambda) - 12*I*omega + 
            12*a^4*z^2*(z - (1 - I)*omega)*(z + (1 + I)*omega) + 24*a^3*m*z^2*((-I)*z + 2*omega) + 
            4*a*m*(2*I*z*(3*z + lambda) + 3*omega) + 4*a^2*(-3*z^2*(-1 + 2*m^2 + 2*z) - 
            2*I*z*(-3 + 3*z + lambda)*omega - 3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function QplusInf_z(s::Int, m::Int, a, omega, lambda, z)
    if s == 0
        return begin
            -((z^2*(-2*z*(-1 + lambda) + lambda - 4*a^2*z^3*(1 + lambda) - 2*a^4*z^5*(3 + lambda) + 2*a*m*omega + 
            a^6*z^6*(1 - m^2 + lambda + 2*a*m*omega) + a^2*z^4*(8 + a^2*(2 - 2*m^2 + 3*lambda) + 
            6*a^3*m*omega) + z^2*(-4 + a^2*(1 - m^2 + 3*lambda) + 6*a^3*m*omega)))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2))
        end
    elseif s == 1
        return begin
            (z^2*(2*a^12*z^12*(1 + lambda) + 4*I*a^10*z^11*(1 + lambda)*(3*I + a^2*omega) - 
            (2 + lambda)*(2 + lambda + 2*a*m*omega) - a^6*z^8*(-4 - 8*lambda + 2*I*a*m*(9 + 2*lambda) + 
            a^2*(31 + 5*lambda^2 - m^2*(1 + 4*lambda) + 8*lambda*(3 + 4*I*omega) + 16*I*omega) + 
            10*a^3*m*(-2 + lambda)*omega) - a^4*z^6*(52 + 32*lambda + 2*I*a*m*(31 + 6*lambda) + 
            a^2*(71 + 10*lambda^2 - 3*m^2*(3 + 2*lambda) + 8*lambda*(7 + 6*I*omega) + 48*I*omega) + 
            20*a^3*m*(-1 + lambda)*omega) - a*z^2*(4*I*m*(2 + lambda) - a*m^2*(4 + lambda) + 
            a*(2 + lambda)*(14 + 5*lambda + 8*I*omega) + 10*a^2*m*(1 + lambda)*omega) + 
            2*z*(I*a*m*(2 + lambda) + (2 + lambda)^2 + 2*I*a^2*(2 + m^2 + lambda)*omega) + 
            2*a^6*z^9*(8 - 6*I*a*m + a^2*(-3 + 3*m^2 - 4*lambda + lambda^2 - 8*I*omega) - 
            I*a^3*m*(-4 + m^2 - lambda - 10*I*omega) + 2*I*a^4*(6 + m^2 + 5*lambda)*omega) + 
            2*a^4*z^7*(-4 + 8*I*a*m + a^3*m*(14*I - 3*I*m^2 + 4*I*lambda - 30*omega) + 
            4*a^2*(10 + m^2 + 7*lambda + lambda^2 - 4*I*omega) + 4*I*a^4*(7 + 2*m^2 + 5*lambda)*omega) + 
            2*I*a^3*z^5*(14*m + a^2*(-3*m^3 + 6*m*(3 + lambda + 5*I*omega)) + 
            I*a*(-69 + m^2 - 48*lambda - 6*lambda^2 + 8*I*omega) + 4*a^3*(8 + 3*m^2 + 5*lambda)*omega) + 
            2*a^2*z^3*(36 - 2*m^2 + 26*lambda + 4*lambda^2 - I*a*(m^3 - 2*m*(5 + 2*lambda + 5*I*omega)) + 
            2*I*a^2*(9 + 4*m^2 + 5*lambda)*omega) - a^8*z^10*(-6*I*a*m - 8*(1 + 2*lambda) + 
            2*a^3*m*(-3 + lambda)*omega + a^2*(1 - m^2*(-1 + lambda) + lambda^2 + 8*I*lambda*omega)) - 
            a^2*z^4*(24*(2 + lambda) + 2*I*a*m*(23 + 6*lambda) + 20*a^3*m*lambda*omega + 
            a^2*(67 + 54*lambda + 10*lambda^2 - m^2*(11 + 4*lambda) + 48*I*omega + 32*I*lambda*omega))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(2 - 2*I*a*m*z - 2*I*a^3*m*z^3 + lambda + 
            a^4*z^4*(1 + lambda) + a^2*z^2*(3 + 2*z + 2*lambda)))
        end
    elseif s == -1
        return begin
            (z^2*((-1 + 2*z)*lambda^2 + 2*a^12*z^11*(-1 + lambda)*(z + 2*I*omega) - 2*a^11*m*z^10*(1 + lambda)*omega + 
            2*a^9*m*z^8*(-3*I*z^2 + I*z*(-2 + m^2 - lambda - 2*I*omega) - (6 + 5*lambda)*omega) + 
            2*I*a^7*m*z^6*(3*m^2*z + 6*z^3 + z^2*(5 + 2*lambda) - 2*z*(3 + 2*lambda + 5*I*omega) + 
            2*I*(7 + 5*lambda)*omega) + 2*I*a^5*m*z^4*(-8*z^3 + z^2*(19 + 6*lambda) + 
            3*z*(m^2 - 2*(1 + lambda + 3*I*omega)) + 2*I*(8 + 5*lambda)*omega) + 
            2*I*a^3*m*z^2*(-14*z^3 + z^2*(11 + 6*lambda) + z*(m^2 - 2*(1 + 2*lambda + 7*I*omega)) + 
            I*(9 + 5*lambda)*omega) + 2*a*m*(2*I*z^2*lambda - (2 + lambda)*omega + z*((-I)*lambda + 4*omega)) - 
            a^10*z^9*(12*z^2*(-1 + lambda) + 4*I*(5 + m^2 - 5*lambda)*omega + 
            z*(5 - m^2*(-3 + lambda) - 4*lambda + lambda^2 - 12*I*omega + 8*I*lambda*omega)) + 
            a^2*z*(-24*z^3*lambda + m^2*(-4*z^2 + z*(2 + lambda) - 4*I*omega) + 
            4*z^2*(5*lambda + 2*lambda^2 - 6*I*omega) + 4*I*(-1 + lambda)*omega - 
            z*(4*lambda + 5*lambda^2 - 20*I*omega + 8*I*lambda*omega)) + 
            a^8*z^7*(8*z^3*(-3 + 2*lambda) + 2*z^2*(9 + 3*m^2 - 8*lambda + lambda^2 - 4*I*omega) - 
            8*I*(5 + 2*m^2 - 5*lambda)*omega + z*(-3 - 4*lambda - 5*lambda^2 + m^2*(-7 + 4*lambda) + 56*I*omega - 
            32*I*lambda*omega)) + a^4*z^3*(-8*z^4 - 4*z^3*(-3 + 8*lambda) - 
            2*z^2*(3 + m^2 - 24*lambda - 6*lambda^2 + 28*I*omega) - 4*I*(5 + 4*m^2 - 5*lambda)*omega + 
            z*(1 - 14*lambda - 10*lambda^2 + m^2*(3 + 4*lambda) + 72*I*omega - 32*I*lambda*omega)) + 
            a^6*z^5*(16*z^4 + 4*z^3*(-3 + 2*lambda) + 8*z^2*(m^2 + 3*lambda + lambda^2 - 5*I*omega) - 
            8*I*(5 + 3*m^2 - 5*lambda)*omega + z*(1 - 16*lambda - 10*lambda^2 + m^2*(-3 + 6*lambda) + 96*I*omega - 
            48*I*lambda*omega))))/((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*
            (2*I*a*m*z + 2*I*a^3*m*z^3 + a^4*z^4*(-1 + lambda) + lambda + a^2*z^2*(-1 + 2*z + 2*lambda)))
        end
    elseif s == 2
        return begin
            (z^2*(24*a^12*z^12 - lambda^3 + 48*I*a^10*z^11*(3*I + a^2*omega) - 2*lambda^2*(8 + a*m*omega) + 
            8*I*a^6*z^9*(24*I + 18*a*m + 3*I*a^2*(9 + m^2) + 
            a^3*m*(19 - 3*m^2 + 4*lambda - 30*I*omega) + a^4*(23 + 3*m^2 - 4*lambda + 24*I*omega)*omega) + 
            4*lambda*(-21 - (3*I + 8*a*m)*omega + 3*a^2*omega^2) + 
            12*a^8*z^10*(24 - 6*I*a*m + a^2*(5 + 3*m^2 - 8*I*omega) - 12*a^3*m*omega + 8*a^4*omega^2) - 
            2*I*a^3*z^5*(-24*m*(9 + 4*lambda) + a*(-4*I*m^2*(91 + 16*lambda) + 
            I*(924 + 47*lambda^2 + lambda^3 + lambda*(394 - 84*I*omega) + 732*I*omega)) + 
            4*a^2*m*(-143 - 3*lambda^2 + m^2*(23 + 2*lambda) + lambda*(-44 + 31*I*omega) + 313*I*omega) + 
            4*a^3*(5 + 3*lambda^2 - m^2*(41 + 8*lambda) + 5*lambda*(7 - 3*I*omega) - 219*I*omega)*omega + 
            24*a^4*m*(1 + lambda)*omega^2) + 24*(-6 - (3*I + 5*a*m)*omega - 
            a*(I*m + a*(-3 + m^2))*omega^2 + a^3*m*omega^3) - 
            4*I*a*z^3*(-12*m*(3 + lambda) + I*a*(408 + 194*lambda + 27*lambda^2 + lambda^3 - 16*m^2*(4 + lambda) + 
            240*I*omega) + 12*a^4*m*(1 + lambda)*omega^2 + 2*a^3*omega*(17 + 26*lambda + 3*lambda^2 - 
            m^2*(31 + 7*lambda) - 132*I*omega - 12*I*lambda*omega) + 
            2*a^2*m*(-98 - 35*lambda - 3*lambda^2 + m^2*(10 + lambda) + 164*I*omega + 20*I*lambda*omega)) + 
            4*a^6*z^8*(36 - 24*m^2 + 2*a^3*m*(-100 + 12*m^2 - 7*lambda + 6*I*omega)*omega + 
            6*a^4*(21 - 5*m^2 + lambda)*omega^2 + 12*a^5*m*omega^3 + 
            2*a*m*(-44*I + 9*I*m^2 - 11*I*lambda + 12*omega) - 
            a^2*(39 + 6*m^4 + 20*lambda + 2*lambda^2 + m^2*(-65 - 8*lambda + 30*I*omega) + 170*I*omega - 
            22*I*lambda*omega)) + 2*z*((9 + lambda)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            8*I*a^3*m*(1 + lambda)*omega^2 + 4*a*m*(24*I + 10*I*lambda + I*lambda^2 + 27*omega + 3*lambda*omega) - 
            4*a^2*omega*(6*I + 7*I*lambda + I*lambda^2 - 2*I*m^2*(4 + lambda) + 27*omega + 3*lambda*omega)) + 
            a^3*z^6*(384*I*m + 16*a*(-33 + 12*m^2 - 10*lambda - lambda^2 - 36*I*omega) + 
            8*I*a^2*m*(-177 + 18*m^2 - 53*lambda - 2*lambda^2 + 96*I*omega) + 
            2*a^4*m*(-810 + 126*m^2 - 100*lambda - lambda^2 + 12*I*omega)*omega + 
            12*a^5*(83 - 31*m^2 + 7*lambda)*omega^2 + 168*a^6*m*omega^3 - 
            a^3*(756 + 48*m^4 + lambda^3 + lambda^2*(41 - 16*I*omega) - 
            m^2*(564 + 106*lambda + lambda^2 - 180*I*omega) + lambda*(334 - 364*I*omega) + 948*I*omega - 
            576*omega^2)) + a*z^4*(96*I*m + 8*I*a^2*m*(-202 + 9*m^2 - 65*lambda - 4*lambda^2 + 
            72*I*omega) - 8*a*(120 + 50*lambda + 5*lambda^2 + 72*I*omega) + 
            6*a^4*m*(-256 + 36*m^2 - 44*lambda - lambda^2 - 20*I*omega)*omega - 12*a^5*(-78 + 34*m^2 - 9*lambda)*
            omega^2 + 216*a^6*m*omega^3 + a^3*(-1068 - 24*m^4 - 3*lambda^3 + 
            2*m^2*(246 + 58*lambda + lambda^2) + lambda^2*(-74 + 32*I*omega) + lambda*(-512 + 388*I*omega) - 
            568*I*omega + 576*omega^2)) - 8*I*a^4*z^7*(12*I - 6*a*m*(-5 + lambda) + 
            2*a^5*m*(1 + lambda)*omega^2 + a^4*omega*(-23 + 20*lambda + lambda^2 - 3*m^2*(7 + lambda) - 138*I*omega - 
            6*I*lambda*omega) + a^3*m*(-88 - 23*lambda - lambda^2 + m^2*(16 + lambda) + 206*I*omega + 
            14*I*lambda*omega) + a^2*(-2*I*m^2*(31 + 4*lambda) + 3*I*(24 + 10*lambda + lambda^2 + 46*I*omega - 
            2*I*lambda*omega))) + z^2*(-12*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            8*I*a*m*(60 + 23*lambda + 2*lambda^2 - 18*I*omega) + 2*a^3*m*(-346 + 30*m^2 - 76*lambda - 
            3*lambda^2 - 60*I*omega)*omega - 60*a^4*(-7 + 3*m^2 - lambda)*omega^2 + 120*a^5*m*omega^3 + 
            a^2*(-3*lambda^3 + lambda^2*(-57 + 16*I*omega) + m^2*(152 + 42*lambda + lambda^2 + 60*I*omega) + 
            lambda*(-342 + 100*I*omega) + 12*(-54 - 23*I*omega + 12*omega^2)))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(24 + 10*lambda + lambda^2 - 
            4*I*a*m*(6*z^2 + 2*z*(4 + lambda) + 3*I*omega) + 12*I*omega + 24*a^3*m*z^2*(I*z + 2*omega) + 
            12*a^4*z^2*(z^2 - 2*I*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*z*(1 + lambda)*omega + 3*omega^2)))
        end
    elseif s == -2
        return begin
            (z^2*(-((-1 + 2*z)*(-2 + 6*z - lambda)*(2*lambda + lambda^2 - 12*I*omega)) + 
            48*a^11*m*z^8*omega^2*(4*I*z + omega) + 24*a^12*z^9*(z^3 + 2*I*z^2*omega - 2*z*omega^2 - 
            4*I*omega^3) + 8*I*a^9*m*z^6*(9*z^4 + z^3*(-3 + 3*m^2 - 4*lambda + 18*I*omega) + 
            3*z^2*(-4*I*m^2 + 3*I*lambda - 18*omega)*omega + 2*z*(45 + lambda)*omega^2 - 21*I*omega^3) + 
            a^4*z*(96*z^6 + 16*z^5*(-9 + 12*m^2 - 2*lambda - lambda^2 + 36*I*omega) - 
            z^3*(12 + 24*m^4 + 3*lambda^3 + lambda^2*(38 + 32*I*omega) + 4*lambda*(16 + 33*I*omega) - 
            2*m^2*(30 + 50*lambda + lambda^2 + 288*I*omega) - 552*I*omega) + 
            2*z^4*(36 + 35*lambda^2 + lambda^3 - 4*m^2*(27 + 16*lambda) + lambda*(66 + 84*I*omega) - 636*I*omega) - 
            8*I*z^2*(3 - 3*lambda^2 + m^2*(51 + 7*lambda) + lambda*(-2 - 24*I*omega) + 36*I*omega)*omega + 
            4*z*(-39 - 45*m^2 + 31*lambda + 48*I*omega)*omega^2 - 96*I*omega^3) - 
            2*I*a^7*m*z^4*(72*z^5 + 4*z^3*(12 + lambda^2 - m^2*(12 + lambda) + 3*lambda*(5 + 6*I*omega) - 
            6*I*omega) + 4*z^4*(9*m^2 - 11*lambda + 36*I*omega) + 
            I*z^2*(6 + 126*m^2 - 124*lambda - lambda^2 - 588*I*omega)*omega - 24*z*(21 + lambda)*omega^2 + 
            108*I*omega^3) - a^6*z^3*(48*(-3 + 2*m^2)*z^5 + 192*z^6 + 
            8*z^4*(m^2*(30 + 8*lambda) - 3*(lambda^2 + lambda*(2 + 2*I*omega) - 18*I*omega)) + 
            z^3*(12 + 48*m^4 + lambda^3 + lambda^2*(29 + 16*I*omega) + lambda*(54 + 236*I*omega) - 
            m^2*(156 + 98*lambda + lambda^2 + 756*I*omega) - 420*I*omega) + 
            8*I*z^2*(3 - 3*lambda^2 + m^2*(81 + 8*lambda) + lambda*(-11 - 27*I*omega) + 57*I*omega)*omega + 
            12*z*(30 + 34*m^2 - 17*lambda - 48*I*omega)*omega^2 + 384*I*omega^3) - 
            4*a^10*z^7*(36*z^4 - 3*z^3*(5 + 3*m^2 - 20*I*omega) + 
            2*z*(27 + 15*m^2 - 5*lambda - 24*I*omega)*omega^2 + 96*I*omega^3 - 
            2*z^2*omega*(15*I - 15*I*m^2 + 4*I*lambda + 24*omega)) - 
            2*I*a^3*m*(192*z^6 + 24*z^5*(-7 + 4*lambda) + 4*z^4*(-6 + 9*m^2 - 33*lambda - 4*lambda^2) + 
            3*I*z^2*(-6 + 10*m^2 - 28*lambda - lambda^2 - 44*I*omega)*omega - 8*z*(9 + lambda)*omega^2 + 
            12*I*omega^3 - 4*z^3*(-6 - 11*lambda - 3*lambda^2 + m^2*(6 + lambda) - 32*I*lambda*omega)) - 
            2*I*a^5*m*z^2*(24*z^5*(-9 + lambda) + 4*z^4*(3 + 18*m^2 - 37*lambda - 2*lambda^2 - 24*I*omega) + 
            3*z^2*(36*I*m^2 - I*(52*lambda + lambda^2 + 172*I*omega))*omega - 24*z*(13 + lambda)*omega^2 + 
            60*I*omega^3 + 4*z^3*(15 + 20*lambda + 3*lambda^2 - m^2*(15 + 2*lambda) + 9*I*omega + 43*I*lambda*omega)) - 
            a^2*(8*z^4*(10*lambda + 5*lambda^2 - 72*I*omega) - 4*z^3*((26 - 16*m^2)*lambda + 15*lambda^2 + lambda^3 - 
            168*I*omega) + 4*(6 + 6*m^2 - 7*lambda)*omega^2 + 8*I*z*omega*(lambda - lambda^2 + 2*m^2*(6 + lambda) + 
            9*I*omega - 7*I*lambda*omega) + z^2*(3*lambda^3 + lambda^2*(21 - m^2 + 16*I*omega) + 
            lambda*(30 - 34*m^2 - 28*I*omega) - 12*I*(15 + 11*m^2 + 4*I*omega)*omega)) + 
            4*a^8*z^5*(72*z^5 - 6*z^4*(9 + m^2 - 12*I*omega) + z*(-99 - 93*m^2 + 37*lambda + 144*I*omega)*
            omega^2 - 144*I*omega^3 - 2*I*z^2*omega*(-9 - 12*lambda - lambda^2 + 3*m^2*(19 + lambda) + 54*I*omega - 
            10*I*lambda*omega) - z^3*(-9 + 6*m^4 + 4*lambda + 2*lambda^2 - m^2*(33 + 8*lambda + 78*I*omega) + 
            54*I*omega + 22*I*lambda*omega + 48*omega^2)) + 2*a*m*(-48*I*z^4 - 24*I*z^3*(-1 + lambda) - 
            (12 + 16*lambda + lambda^2 - 12*I*omega)*omega + 4*z^2*(7*I*lambda + 2*I*lambda^2 + 6*omega) + 
            4*z*((-I)*lambda^2 + 3*omega + lambda*(-2*I + 7*omega)))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(2*lambda + lambda^2 - 12*I*omega + 
            24*a^3*m*z^2*((-I)*z + 2*omega) + 4*a*m*(6*I*z^2 + 2*I*z*lambda + 3*omega) + 
            12*a^4*z^2*(z^2 + 2*I*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*z*(-3 + lambda)*omega + 3*omega^2)))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

# Cache mechanism for the outgoing coefficients at infinity
# Initialize the cache with a set of fiducial parameters
_cached_outgoing_coefficients_at_inf_params::NamedTuple{(:s, :m, :a, :omega, :lambda), Tuple{Int, Int, _DEFAULTDATATYPE, _DEFAULTDATATYPE, _DEFAULTDATATYPE}} = (s=-2, m=2, a=0, omega=0.5, lambda=1)
_cached_outgoing_coefficients_at_inf::NamedTuple{(:expansion_coeffs, :Pcoeffs, :Qcoeffs), Tuple{Vector{_DEFAULTDATATYPE}, Vector{_DEFAULTDATATYPE}, Vector{_DEFAULTDATATYPE}}} = (
    expansion_coeffs = [_DEFAULTDATATYPE(1.0)], 
    Pcoeffs = [_DEFAULTDATATYPE(0.0)],
    Qcoeffs = [_DEFAULTDATATYPE(0.0)]
)

function outgoing_coefficient_at_inf(s::Int, m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    for up to (1/r)^3, which is probably more than enough
    
    But we have also shown a recurrence relation where
    one can generate as higher an order as one pleases.
    However, the recurrence relation that we have
    actually depends *all* previous terms so this function
    is designed to be evaluated recursively to build
    the full list of coefficients
    =#
    global _cached_outgoing_coefficients_at_inf_params
    global _cached_outgoing_coefficients_at_inf

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    if order == 0
        return 1.0 # This is always 1.0
    elseif order == 1
        if s == 0
            return (1/2)*I*(lambda + 2*a*m*omega)
        elseif s == +1
            return (1/2)*I*(2 + lambda + 2*a*m*omega)
        elseif s == -1
            return (I*(lambda^2 + 2*a*m*(2 + lambda)*omega))/(2*lambda)
        elseif s == +2
            return (1/2)*I*(6 + lambda + 2*a*m*omega)
        elseif s == -2
            return begin
                (1/6)*I*(-6 + 7*lambda) + I*a*m*omega + 
                (2*I*(-3 + lambda)*lambda*(2 + lambda) + 
                24*(-3 - 3*I*a*m + lambda)*omega)/(-3*lambda*(2 + lambda) + 
                36*omega*(I - a*m + a^2*omega))
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    elseif order == 2
        if s == 0
            return begin
                (1/8)*(-lambda^2 + lambda*(2 - 4*a*m*omega) + 
                4*omega*(I + a*m + a*m*(2*I - a*m)*omega))
            end
        elseif s == +1
            return begin
                (1/8)*(-lambda^2 - 2*lambda*(1 + 2*a*m*omega) - 
                4*a*omega*(m - 2*a*omega - 2*I*m*omega + a*m^2*omega))
            end
        elseif s == -1
            return begin
                -((1/(8*lambda))*((-2 + lambda)*lambda^2 + 4*a*m*(-2 + lambda + lambda^2)*omega + 
                4*a*(-2*I*m*lambda + a*(2 - 2*lambda + m^2*(4 + lambda)))*omega^2))
            end
        elseif s == +2
            return begin
                (1/8)*(-lambda^2 - 2*lambda*(5 + 2*a*m*omega) - 
                4*(6 + (3*I + 5*a*m)*omega + a*m*(-2*I + a*m)*omega^2))
            end
        elseif s == -2
            return begin
                -((1/(8*lambda*(2 + lambda) - 96*omega*(I - a*m + a^2*omega)))*
                (lambda^2*(2 + lambda)^2 + 4*a*m*lambda*(16 + lambda*(14 + lambda))*
                omega + 4*(36 + a*(a*(10 - 11*lambda)*lambda - 
                2*I*m*(2 + lambda)*(6 + lambda) + 
                a*m^2*(60 + lambda*(30 + lambda))))*omega^2 + 
                16*a*(-6*m + a*(3*I - 9*I*m^2 + 
                a*m*(-15 + 3*m^2 - 7*lambda)))*omega^3 - 
                48*a^3*(-2*I*m + a*(-4 + m^2))*omega^4))
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    elseif order == 3
        if s == 0
            return begin
                (1/2)*(1 - I*a*m)*omega - (1/48)*I*((-6 + lambda)*(-2 + lambda)*lambda + 
                2*lambda*(-18*I + a*m*(-16 + 3*lambda))*omega + 
                4*a*(3*a*m^2*(-2 + lambda) + 2*a*(-1 + lambda) - 
                6*I*m*(1 + lambda))*omega^2 + 8*a*m*(-8 - 6*I*a*m + 
                a^2*(2 + m^2))*omega^3)
            end
        elseif s == +1
            return begin
                (-(1/48))*I*(lambda^3 + lambda^2*(-2 + 6*a*m*omega) + 
                4*lambda*(-2 - 2*(3*I + a*m)*omega + 
                a*(-4*a - 6*I*m + 3*a*m^2)*omega^2) + 
                8*omega*(-6*I + a^3*m*(-4 + m^2)*omega^2 + 
                3*a^2*omega*(-1 - 2*I*m^2*omega) - 
                a*m*(3 + 6*I*omega + 8*omega^2)))
            end
        elseif s == -1
            return begin
                -((1/(48*lambda))*(I*((-6 + lambda)*(-2 + lambda)*lambda^2 + 48*a*m*omega + 
                2*lambda*(-12*I*lambda + a*m*(-16 + lambda*(-10 + 3*lambda)))*omega + 
                4*a*(3*a*m^2*(-2 + lambda)*(4 + lambda) - 
                4*a*(3 + (-2 + lambda)*lambda) - 6*I*m*(4 + lambda^2))*omega^2 + 
                8*a*(-8*m*lambda - 6*I*a*(-2 + m^2*(2 + lambda)) + 
                a^2*m*(6 - 4*lambda + m^2*(6 + lambda)))*omega^3)))
            end
        elseif s == +2
            return begin
                (-(1/48))*I*(lambda^3 + 2*lambda^2*(5 + 3*a*m*omega) + 
                4*lambda*(6 + (3*I + 10*a*m)*omega + 
                a*(2*a - 6*I*m + 3*a*m^2)*omega^2) + 
                8*a*omega*(a*omega + 6*a*m^2*(1 - I*omega)*omega + 
                a^2*m^3*omega^2 + m*(2 - 9*I*omega + 2*(-4 + a^2)*
                omega^2)))
            end
        elseif s == -2
            return begin
                -((1/(48*(lambda*(2 + lambda) - 12*omega*(I - a*m + a^2*omega))))*
                (I*(-4 + lambda)*lambda^2*(2 + lambda)^2 + 
                2*I*a*m*lambda*(-96 + lambda*(-44 + lambda*(32 + 3*lambda)))*
                omega + 4*I*(36*(-4 + lambda) + 
                a*(-6*I*m*lambda*(2 + lambda)^2 + 
                a*lambda*(-60 + (40 - 13*lambda)*lambda) + 
                3*a*m^2*(-48 + lambda*(40 + lambda*(24 + lambda)))))*
                omega^2 + 8*I*a*(-6*I*a*(-3*(2 + lambda) + 
                m^2*(1 + lambda)*(18 + lambda)) - 
                4*m*(-9 + lambda*(13 + 2*lambda)) + 
                a^2*m*(108 - lambda*(44 + 31*lambda) + 
                m^2*(144 + lambda*(44 + lambda))))*omega^3 - 
                48*a*(16*m + 28*I*a*m^2 - 2*a^2*m*
                (5 + 7*m^2 - 7*lambda) - 
                I*a^3*(2*m^4 + 2*(-9 + 5*lambda) - 
                m^2*(32 + 7*lambda)))*omega^4 - 
                96*I*a^3*m*(-8 - 6*I*a*m + a^2*(-10 + m^2))*
                omega^5))
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    else
        # Evaluate higher order corrections using AD

        _this_params = (s=s, m=m, a=a, omega=omega, lambda=lambda)
        # Check if we can use the cached results
        if _cached_outgoing_coefficients_at_inf_params == _this_params
            expansion_coeffs = _cached_outgoing_coefficients_at_inf.expansion_coeffs
            Pcoeffs = _cached_outgoing_coefficients_at_inf.Pcoeffs
            Qcoeffs = _cached_outgoing_coefficients_at_inf.Qcoeffs
        else
            # Cannot re-use the cached results, re-compute from zero
            expansion_coeffs = [data_type(1.0)] # order 0
            Pcoeffs = [data_type(PplusInf_z(s, m, a, omega, lambda, 0))] # order 0
            Qcoeffs = [data_type(0.0), data_type(0.0)] # the recurrence relation takes Q_{r+1}
        end

        # Compute Pcoeffs to the necessary order
        _P(z) = PplusInf_z(s, m, a, omega, lambda, z)
        _P_taylor = taylor_expand(_P, 0, order=order)
        for i in length(Pcoeffs):order
            append!(Pcoeffs, getcoeff(_P_taylor, i))
        end

        # Compute Qcoeffs to the necessary order (to current order + 1)
        _Q(z) = QplusInf_z(s, m, a, omega, lambda, z)
        _Q_taylor = taylor_expand(_Q, 0, order=order+1)
        for i in length(Qcoeffs):order+1
            append!(Qcoeffs, getcoeff(_Q_taylor, i))
        end

        # Note that the expansion coefficients we store is scaled by \omega^{i}
        for i in length(expansion_coeffs):order
            _P0 = Pcoeffs[1] # P0
            sum = 0.0
            for k in 1:i
                sum += (Qcoeffs[k+2] - (i-k)*Pcoeffs[k+1])*(expansion_coeffs[i-k+1]/omega^(i-k))
            end
            append!(expansion_coeffs, omega^(i)*((i*(i-1)*(expansion_coeffs[i]/omega^(i-1)) + sum)/(_P0*i)))
        end

        # Update cache
        _cached_outgoing_coefficients_at_inf_params = _this_params
        _cached_outgoing_coefficients_at_inf = (
            expansion_coeffs = expansion_coeffs,
            Pcoeffs = Pcoeffs,
            Qcoeffs = Qcoeffs
        )

        return expansion_coeffs[order+1]
    end
end

function PplusH(s::Int, m::Int, a, omega, lambda, x)
    if s == 0
        return begin
            ((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            2*I*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)))/(2*sqrt(1 - a^2) + x)
        end
    elseif s == 1
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(-3*I*a^2*m*(1 + sqrt(1 - a^2) + x) - 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(-2*I*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))) + 2*I*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)))
        end
    elseif s == -1
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(3*I*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(2*I*a^3*m*(1 + sqrt(1 - a^2) + x) + 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            lambda + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) + 
            2*I*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)))
        end
    elseif s == 2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            2*I*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(I*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I + 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - I*(1 + sqrt(1 - a^2) + x)^2*
            (1 + lambda)*omega) + a^3*(-6 + 9*I*(1 + sqrt(1 - a^2) + x)*omega + 
            6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((-(1 + sqrt(1 - a^2) + x)^4)*
            (24 + 10*lambda + lambda^2 + 12*I*omega) - 24*a^3*m*(1 + sqrt(1 - a^2) + x)*
            (I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*
            (6*I + 2*I*(1 + sqrt(1 - a^2) + x)*(4 + lambda) - 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 + 2*I*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 - 6*I*omega) - 2*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    elseif s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            2*I*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega) - 
            (8*a*x*(2*sqrt(1 - a^2) + x)*((-I)*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + I*(1 + sqrt(1 - a^2) + x)^2*
            (-3 + lambda)*omega) + a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 
            6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^4*
            (2*lambda + lambda^2 - 12*I*omega) + 24*a^3*m*(1 + sqrt(1 - a^2) + x)*
            (-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*
            (6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 + 6*I*omega) + 2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function QplusH(s::Int, m::Int, a, omega, lambda, x)
    if s == 0
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            (x^2*(1 + sqrt(1 - a^2) + x)^2*(2*sqrt(1 - a^2) + x)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(a^4 - 4*a^2*(1 + sqrt(1 - a^2) + x) - 
            (-3 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)^3) + 
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(x*(2*sqrt(1 - a^2) + x)*lambda - 
            (a*m - (a^2 + (1 + sqrt(1 - a^2) + x)^2)*omega)^2))))
        end
    elseif s == 1
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            ((2*I*a*x*(2*sqrt(1 - a^2) + x)*(-3*I*a^2*m*(1 + sqrt(1 - a^2) + x) - 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda)))*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + 
            omega))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-2*I*a^3*m*(1 + sqrt(1 - a^2) + x) - 2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*(2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda)))) - 
            (-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            (2*x*(2*sqrt(1 - a^2) + x)*(a^4 - a^2*(1 + sqrt(1 - a^2) + x) - 
            (-2 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)^3) + 
            ((1 + sqrt(1 - a^2) + x)^2*(-1 + 2*sqrt(1 - a^2) + 2*x) + 
            a^2*(1 + 2*sqrt(1 - a^2) + 2*x))^2 + (2*a*x*(2*sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^2*
            (-1 + 2*sqrt(1 - a^2) + 2*x) + a^2*(1 + 2*sqrt(1 - a^2) + 2*x))*
            (-3*I*a^2*m*(1 + sqrt(1 - a^2) + x) - I*m*(1 + sqrt(1 - a^2) + x)^3 + 
            2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*(3 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*(-2*I*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (2*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(-lambda - I*(1 + sqrt(1 - a^2) + x)*omega) + 
            a^8*(1 + lambda)*(2 + (1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            2*I*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-6 + (1 + sqrt(1 - a^2) + x)*
            (4 + m^2 - lambda - 4*I*omega) - 3*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2) - 2*I*a*m*(1 + sqrt(1 - a^2) + x)^5*
            (3 + (1 + sqrt(1 - a^2) + x)*(-5 + 2*lambda) + (1 + sqrt(1 - a^2) + x)^2*
            (2 - lambda + 2*I*omega) - I*(1 + sqrt(1 - a^2) + x)^3*(3 + lambda)*omega + 
            (1 + sqrt(1 - a^2) + x)^4*omega^2) - 2*I*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (13 + (1 + sqrt(1 - a^2) + x)*(-17 + 2*lambda) + (1 + sqrt(1 - a^2) + x)^2*
            (6 + m^2 - 2*lambda - 2*I*omega) - 3*I*(1 + sqrt(1 - a^2) + x)^3*(2 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^4*omega^2) + (1 + sqrt(1 - a^2) + x)^6*
            (-3*(2 + lambda) - (1 + sqrt(1 - a^2) + x)^2*lambda*(2 + lambda) + 
            2*(1 + sqrt(1 - a^2) + x)*(2 + 3*lambda + lambda^2) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda)*omega^2) + a^6*(1 + sqrt(1 - a^2) + x)*(-16*(1 + lambda) + 
            (1 + sqrt(1 - a^2) + x)*(5 + m^2*(-1 + lambda) + 6*lambda - lambda^2 + 2*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(5 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*
            (I + 2*I*m^2 + omega)) + a^4*(1 + sqrt(1 - a^2) + x)^2*
            (11 + 25*lambda + 2*(1 + sqrt(1 - a^2) + x)*(-5 + 3*m^2 - 13*lambda + lambda^2 - 2*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^2*(5 + 4*lambda - 3*lambda^2 + m^2*(3 + 2*lambda) - 8*I*omega) + 
            3*(1 + sqrt(1 - a^2) + x)^4*(3 + 2*lambda)*omega^2 + 4*(1 + sqrt(1 - a^2) + x)^3*omega*
            (I + 2*I*m^2 + omega)) + a^2*(1 + sqrt(1 - a^2) + x)^3*
            (30 + 5*(1 + sqrt(1 - a^2) + x)*(-7 + 2*lambda) - 4*(1 + sqrt(1 - a^2) + x)^2*
            (-2 + m^2 + lambda - lambda^2 - 3*I*omega) + (1 + sqrt(1 - a^2) + x)^3*
            (2 - 2*lambda - 3*lambda^2 + m^2*(4 + lambda) - 10*I*omega) + (1 + sqrt(1 - a^2) + x)^5*
            (7 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^4*omega*(I + 2*I*m^2 + omega))))/
            ((1 + sqrt(1 - a^2) + x)^2*(-2*I*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))))))
        end
    elseif s == -1
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            ((2*I*a*x*(2*sqrt(1 - a^2) + x)*(3*I*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + 
            omega))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (2*I*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) - 
            (-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            ((a^2 - (1 + sqrt(1 - a^2) + x)^2)^2 + 2*x*(1 + sqrt(1 - a^2) + x)*
            (2*sqrt(1 - a^2) + x)*(-3*a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(a^2 - (1 + sqrt(1 - a^2) + x)^2)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(3*I*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (2*I*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) - 
            ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(2*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*
            (2 - lambda + I*(1 + sqrt(1 - a^2) + x)*omega) + a^8*(-1 + lambda)*
            (2 + (1 + sqrt(1 - a^2) + x)^2*omega^2) + 2*I*a^5*m*(1 + sqrt(1 - a^2) + x)^3*
            (-2 + m^2 - lambda + 4*I*omega - 3*I*(1 + sqrt(1 - a^2) + x)*omega + 
            3*I*(1 + sqrt(1 - a^2) + x)*lambda*omega + 3*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 
            (1 + sqrt(1 - a^2) + x)^6*lambda*(-3 - (1 + sqrt(1 - a^2) + x)^2*lambda + 
            2*(1 + sqrt(1 - a^2) + x)*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^5*(sqrt(1 - a^2) + x + 
            2*(1 + sqrt(1 - a^2) + x)*lambda - (1 + sqrt(1 - a^2) + x)^2*(lambda + 2*I*omega) + 
            I*(1 + sqrt(1 - a^2) + x)^3*(1 + lambda)*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            2*I*a^3*m*(1 + sqrt(1 - a^2) + x)^3*(1 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda) + 
            (1 + sqrt(1 - a^2) + x)^2*(-2 + m^2 - 2*lambda + 2*I*omega) + 
            3*I*(1 + sqrt(1 - a^2) + x)^3*lambda*omega + 3*(1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            a^6*(1 + sqrt(1 - a^2) + x)*(8 - 8*lambda + (1 + sqrt(1 - a^2) + x)*
            (-1 + m^2*(-3 + lambda) - lambda^2 - 2*I*omega) + (1 + sqrt(1 - a^2) + x)^3*(-3 + 4*lambda)*
            omega^2 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(-I - 2*I*m^2 + omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^2*(-11 + 9*lambda + 2*(1 + sqrt(1 - a^2) + x)*
            (1 + 3*m^2 + lambda + lambda^2 + 2*I*omega) + (1 + sqrt(1 - a^2) + x)^2*
            (1 - 6*lambda - 3*lambda^2 + m^2*(-1 + 2*lambda) + 8*I*omega) + 3*(1 + sqrt(1 - a^2) + x)^4*
            (-1 + 2*lambda)*omega^2 + 4*(1 + sqrt(1 - a^2) + x)^3*omega*(-I - 2*I*m^2 + omega)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^3*(6 - 3*(1 + sqrt(1 - a^2) + x)*(1 + 2*lambda) - 
            4*(1 + sqrt(1 - a^2) + x)^2*(m^2 - 3*lambda - lambda^2 + 3*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(-4*lambda - 3*lambda^2 + m^2*(2 + lambda) + 10*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^5*(-1 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^4*omega*
            (-I - 2*I*m^2 + omega))))/((1 + sqrt(1 - a^2) + x)^2*
            (2*I*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))))))
        end
    elseif s == 2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - 
            (8*I*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*
            (I*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*I + 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - 
            I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega) + 
            a^3*(-6 + 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(24 + 10*lambda + lambda^2 + 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(-6*I - 2*I*(1 + sqrt(1 - a^2) + x)*(4 + lambda) + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) - 12*a^4*(-1 + 2*I*(1 + sqrt(1 - a^2) + x)*omega + 
            2*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - (((1 + sqrt(1 - a^2) + x)^2*(-1 + 3*sqrt(1 - a^2) + 3*x) + 
            a^2*(1 + 3*sqrt(1 - a^2) + 3*x))^2 + x*(2*sqrt(1 - a^2) + x)*
            (3*a^4 + (1 + sqrt(1 - a^2) + x)^3*(8 - 3*(1 + sqrt(1 - a^2) + x))) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^2*(-1 + 3*sqrt(1 - a^2) + 3*x) + 
            a^2*(1 + 3*sqrt(1 - a^2) + 3*x))*(I*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I + 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - 
            I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega) + 
            a^3*(-6 + 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*((-(1 + sqrt(1 - a^2) + x)^4)*(24 + 10*lambda + lambda^2 + 
            12*I*omega) - 24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I + 2*(1 + sqrt(1 - a^2) + x)*
            omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*
            (4 + lambda) - 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 + 2*I*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*
            omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 - 6*I*omega) - 2*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(3 + 3*I*(1 + sqrt(1 - a^2) + x)*omega - 
            4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - (1 + sqrt(1 - a^2) + x)^6*
            (24 + 10*lambda + lambda^2 + 12*I*omega)*(-12 - (1 + sqrt(1 - a^2) + x)^2*lambda + 
            2*(1 + sqrt(1 - a^2) + x)*(4 + lambda) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-54*I - 2*I*(1 + sqrt(1 - a^2) + x)*
            (-55 + 3*m^2 - 4*lambda - 18*I*omega) + 2*(1 + sqrt(1 - a^2) + x)^2*
            (4 + 12*m^2 - 8*lambda + 27*I*omega)*omega - 6*I*(1 + sqrt(1 - a^2) + x)^3*(8 + lambda)*
            omega^2 + 45*(1 + sqrt(1 - a^2) + x)^4*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 - 2*I*(1 + sqrt(1 - a^2) + x)^3*
            omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 4*a^6*(1 + sqrt(1 - a^2) + x)*
            (60 - 3*(1 + sqrt(1 - a^2) + x)*(16 + 3*m^2 + 18*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(1 + 36*m^2 - 8*lambda + 18*I*omega)*omega^2 - 
            2*I*(1 + sqrt(1 - a^2) + x)^4*(7 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 - 2*(1 + sqrt(1 - a^2) + x)^2*omega*(-40*I + 9*I*m^2 - 4*I*lambda + 9*omega)) + 
            2*a*m*(1 + sqrt(1 - a^2) + x)^4*(144*I + 8*I*(1 + sqrt(1 - a^2) + x)*
            (-38 + lambda) + (1 + sqrt(1 - a^2) + x)^4*(40 + 20*lambda + lambda^2 + 24*I*omega)*omega + 
            4*I*(1 + sqrt(1 - a^2) + x)^5*(4 + lambda)*omega^2 - 6*(1 + sqrt(1 - a^2) + x)^6*
            omega^3 + 4*(1 + sqrt(1 - a^2) + x)^2*(16*I + 3*I*lambda + 2*I*lambda^2 + 6*omega) - 
            4*I*(1 + sqrt(1 - a^2) + x)^3*(-12 + lambda + lambda^2 - 14*I*omega - 5*I*lambda*omega)) - 
            2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*(312*I + (1 + sqrt(1 - a^2) + x)^3*
            (84 + 30*m^2 - 52*lambda - lambda^2 + 36*I*omega)*omega - 4*I*(1 + sqrt(1 - a^2) + x)^4*
            (19 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3 + 
            12*(1 + sqrt(1 - a^2) + x)*(-63*I + 3*I*m^2 - 3*I*lambda + 32*omega) + 
            4*(1 + sqrt(1 - a^2) + x)^2*(71*I + I*lambda^2 - I*m^2*(10 + lambda) - 104*omega + 
            lambda*(9*I + 16*omega))) + a^2*(1 + sqrt(1 - a^2) + x)^3*
            (576 + 96*(1 + sqrt(1 - a^2) + x)*(-12 + 4*m^2 - 3*I*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^2*(18 - 3*lambda^2 + m^2*(-58 + 8*lambda) + 
            lambda*(-30 - 2*I*omega) + 46*I*omega) - 8*I*(1 + sqrt(1 - a^2) + x)^4*
            (-lambda^2 + 2*m^2*(7 + lambda) + lambda*(2 + 5*I*omega) + 11*I*omega)*omega + 
            2*(1 + sqrt(1 - a^2) + x)^5*(-34 + 24*m^2 - 20*lambda - lambda^2 - 24*I*omega)*omega^2 - 
            8*I*(1 + sqrt(1 - a^2) + x)^6*(1 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^7*
            omega^4 + (1 + sqrt(1 - a^2) + x)^3*(lambda^3 + 36*lambda*(4 + I*omega) + 
            2*lambda^2*(11 - 8*I*omega) - m^2*(-136 + 42*lambda + lambda^2 - 36*I*omega) - 
            8*(-18 + 19*I*omega + 6*omega^2))) + a^4*(1 + sqrt(1 - a^2) + x)^2*
            (-672 + (1 + sqrt(1 - a^2) + x)*(960 - 72*m^2 + 624*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^4*(44 + 180*m^2 - 62*lambda - lambda^2 + 36*I*omega)*omega^2 - 
            8*I*(1 + sqrt(1 - a^2) + x)^5*(5 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^6*
            omega^4 + 8*(1 + sqrt(1 - a^2) + x)^3*omega*(56*I + I*lambda^2 - 3*I*m^2*(9 + lambda) - 
            46*omega + lambda*(6*I + 8*omega)) + 4*(1 + sqrt(1 - a^2) + x)^2*
            (6*m^4 + m^2*(7 - 8*lambda + 54*I*omega) + 2*(-15 + 10*lambda + lambda^2 - 144*I*omega - 9*I*
            lambda*omega + 48*omega^2)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*(4 + lambda) - 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + 12*a^4*(-1 + 2*I*(1 + sqrt(1 - a^2) + x)*
            omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))))/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4))
        end
    elseif s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - 
            (8*I*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*
            ((-I)*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*lambda) + 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 
            24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            ((a^2*(-1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^3)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(-a^4 - 8*a^2*(1 + sqrt(1 - a^2) + x) + 
            (1 + sqrt(1 - a^2) + x)^4) + (8*a*x*(2*sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(a^2*(-1 + sqrt(1 - a^2) + x) + 
            (1 + sqrt(1 - a^2) + x)^3)*((-I)*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 
            24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(-3 + 3*I*(1 + sqrt(1 - a^2) + x)*omega + 
            4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (72*I - 4*I*(1 + sqrt(1 - a^2) + x)*(15 + 9*m^2 - 13*lambda) + 
            4*(1 + sqrt(1 - a^2) + x)^2*(I*m^2*(6 + lambda) - I*(3 + lambda^2 + lambda*
            (9 + 16*I*omega) + 24*I*omega)) + (1 + sqrt(1 - a^2) + x)^3*
            (-36 + 30*m^2 - 44*lambda - lambda^2 - 36*I*omega)*omega + 4*I*(1 + sqrt(1 - a^2) + x)^4*
            (3 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 + 2*I*(1 + sqrt(1 - a^2) + x)^3*
            omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 4*a^6*(1 + sqrt(1 - a^2) + x)*
            (12 - 3*(1 + sqrt(1 - a^2) + x)*(-4 + 3*m^2 + 6*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(-39 + 36*m^2 - 8*lambda - 18*I*omega)*omega^2 + 
            2*I*(1 + sqrt(1 - a^2) + x)^4*(3 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(24*I + 9*I*m^2 - 4*I*lambda + 15*omega)) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-18*I + 2*(1 + sqrt(1 - a^2) + x)^2*
            (-36 + 12*m^2 - 8*lambda - 27*I*omega)*omega + 6*I*(1 + sqrt(1 - a^2) + x)^3*(4 + lambda)*
            omega^2 + 45*(1 + sqrt(1 - a^2) + x)^4*omega^3 + 2*(1 + sqrt(1 - a^2) + x)*
            (9*I + 3*I*m^2 - 4*I*lambda + 30*omega)) - 2*a*m*(1 + sqrt(1 - a^2) + x)^4*
            (-48*I - 24*I*(1 + sqrt(1 - a^2) + x)*(-2 + lambda) + 
            4*I*(1 + sqrt(1 - a^2) + x)^2*(7*lambda + 2*lambda^2 + 6*I*omega) - 
            (1 + sqrt(1 - a^2) + x)^4*(12*lambda + lambda^2 - 24*I*omega)*omega + 
            4*I*(1 + sqrt(1 - a^2) + x)^5*lambda*omega^2 + 6*(1 + sqrt(1 - a^2) + x)^6*omega^3 - 
            4*I*(1 + sqrt(1 - a^2) + x)^3*(lambda + lambda^2 + 6*I*omega + 5*I*lambda*omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^3*(24*m^4*(1 + sqrt(1 - a^2) + x) + 
            48*(-4 + 3*I*omega) + 8*(1 + sqrt(1 - a^2) + x)*(9 + lambda^2 + lambda*(2 + 13*I*omega) - 
            72*I*omega) - (1 + sqrt(1 - a^2) + x)^3*(60 + 54*lambda + lambda^2 + 36*I*omega)*omega^2 + 
            8*I*(1 + sqrt(1 - a^2) + x)^4*(-3 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 + 8*(1 + sqrt(1 - a^2) + x)^2*omega*(24*I - 6*I*lambda - I*lambda^2 + 18*omega + 
            8*lambda*omega) + 4*m^2*(30 - (1 + sqrt(1 - a^2) + x)*(33 + 8*lambda + 54*I*omega) + 
            6*I*(1 + sqrt(1 - a^2) + x)^2*(5 + lambda)*omega + 45*(1 + sqrt(1 - a^2) + x)^3*
            omega^2)) + (1 + sqrt(1 - a^2) + x)^6*((-1 + sqrt(1 - a^2) + x)*
            (1 + sqrt(1 - a^2) + x)*lambda^3 + lambda^2*(12 - 12*(1 + sqrt(1 - a^2) + x) + 
            2*(1 + sqrt(1 - a^2) + x)^2 - (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            12*I*omega*(-12 + 8*(1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            2*lambda*(-12 + 4*(1 + sqrt(1 - a^2) + x)*(2 - 3*I*omega) + 
            6*I*(1 + sqrt(1 - a^2) + x)^2*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^4*(96 + 8*(1 + sqrt(1 - a^2) + x)*
            (m^2*(6 + 8*lambda) - 3*(2 + lambda^2 + lambda*(2 + 2*I*omega) - 22*I*omega)) - 96*I*omega + 
            2*(1 + sqrt(1 - a^2) + x)^4*(6 + 24*m^2 - 12*lambda - lambda^2 + 24*I*omega)*omega^2 + 
            8*I*(1 + sqrt(1 - a^2) + x)^5*(-3 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^6*
            omega^4 + (1 + sqrt(1 - a^2) + x)^2*(lambda^3 + lambda*(24 - 34*m^2 - 4*I*omega) + 
            lambda^2*(14 - m^2 + 16*I*omega) - 12*I*(22 + 3*m^2 - 4*I*omega)*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^3*omega*((-I)*lambda^2 + 2*I*m^2*(3 + lambda) + 3*omega + 
            lambda*(2*I + 5*omega)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I - 2*(1 + sqrt(1 - a^2) + x)*omega) - 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + 12*a^4*(-1 - 2*I*(1 + sqrt(1 - a^2) + x)*
            omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

# Cache mechanism for the outgoing coefficients at horizon
# Initialize the cache with a set of fiducial parameters
_cached_outgoing_coefficients_at_hor_params::NamedTuple{(:s, :m, :a, :omega, :lambda), Tuple{Int, Int, _DEFAULTDATATYPE, _DEFAULTDATATYPE, _DEFAULTDATATYPE}} = (s=-2, m=2, a=0, omega=0.5, lambda=1)
_cached_outgoing_coefficients_at_hor::NamedTuple{(:expansion_coeffs, :Pcoeffs, :Qcoeffs), Tuple{Vector{_DEFAULTDATATYPE}, Vector{_DEFAULTDATATYPE}, Vector{_DEFAULTDATATYPE}}} = (
    expansion_coeffs = [_DEFAULTDATATYPE(1.0)], 
    Pcoeffs = [_DEFAULTDATATYPE(0.0)],
    Qcoeffs = [_DEFAULTDATATYPE(0.0)]
)

function outgoing_coefficient_at_hor(s::Int, m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    global _cached_outgoing_coefficients_at_hor_params
    global _cached_outgoing_coefficients_at_hor

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    _this_params = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    # Check if we can use the cached results
    if _cached_outgoing_coefficients_at_hor_params == _this_params
        expansion_coeffs = _cached_outgoing_coefficients_at_hor.expansion_coeffs
        Pcoeffs = _cached_outgoing_coefficients_at_hor.Pcoeffs
        Qcoeffs = _cached_outgoing_coefficients_at_hor.Qcoeffs
    else
        # Cannot re-use the cached results, re-compute from zero
        expansion_coeffs = [data_type(1.0)] # order 0
        Pcoeffs = [data_type(PplusH(s, m, a, omega, lambda, 0))] # order 0
        Qcoeffs = [data_type(0.0)] # order 0
    end

    if order > 0
        # Compute series expansion coefficients for P and Q
        _P(x) = PplusH(s, m, a, omega, lambda, x)
        _Q(x) = QplusH(s, m, a, omega, lambda, x)
        _P_taylor = taylor_expand(_P, 0, order=order)
        _Q_taylor = taylor_expand(_Q, 0, order=order)

        for i in length(Pcoeffs):order
            append!(Pcoeffs, getcoeff(_P_taylor, i))
            append!(Qcoeffs, getcoeff(_Q_taylor, i))       
        end
    end
    # Define the indicial polynomial
    indicial(nu) = nu*(nu - 1) + Pcoeffs[1]*nu + Qcoeffs[1]

    if order > 0
        # Evaluate the C coefficients
        for i in length(expansion_coeffs):order
            sum = 0.0
            for r in 0:i-1
                sum += expansion_coeffs[r+1]*(r*Pcoeffs[i-r+1] + Qcoeffs[i-r+1])
            end
            append!(expansion_coeffs, -sum/indicial(i))
        end
    end

    # Update cache
    _cached_outgoing_coefficients_at_hor_params = _this_params
    _cached_outgoing_coefficients_at_hor = (
        expansion_coeffs = expansion_coeffs,
        Pcoeffs = Pcoeffs,
        Qcoeffs = Qcoeffs
    )
    return expansion_coeffs[order+1]
end

function PminusH(s::Int, m::Int, a, omega, lambda, x)
    if s == 0
        return begin
            ((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            I*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega)))/(2*sqrt(1 - a^2) + x)
        end
    elseif s == 1
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(-3*I*a^2*m*(1 + sqrt(1 - a^2) + x) - 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(-2*I*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))) + I*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega)))
        end
    elseif s == -1
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(3*I*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(2*I*a^3*m*(1 + sqrt(1 - a^2) + x) + 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            lambda + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) + 
            I*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega)))
        end
    elseif s == 2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            I*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(I*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I + 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - I*(1 + sqrt(1 - a^2) + x)^2*
            (1 + lambda)*omega) + a^3*(-6 + 9*I*(1 + sqrt(1 - a^2) + x)*omega + 
            6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((-(1 + sqrt(1 - a^2) + x)^4)*
            (24 + 10*lambda + lambda^2 + 12*I*omega) - 24*a^3*m*(1 + sqrt(1 - a^2) + x)*
            (I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*
            (6*I + 2*I*(1 + sqrt(1 - a^2) + x)*(4 + lambda) - 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 + 2*I*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 - 6*I*omega) - 2*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    elseif s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            I*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega) - 
            (8*a*x*(2*sqrt(1 - a^2) + x)*((-I)*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + I*(1 + sqrt(1 - a^2) + x)^2*
            (-3 + lambda)*omega) + a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 
            6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^4*
            (2*lambda + lambda^2 - 12*I*omega) + 24*a^3*m*(1 + sqrt(1 - a^2) + x)*
            (-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*
            (6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 + 6*I*omega) + 2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function QminusH(s::Int, m::Int, a, omega, lambda, x)
    if s == 0
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            (x^2*(1 + sqrt(1 - a^2) + x)^2*(2*sqrt(1 - a^2) + x)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(a^4 - 4*a^2*(1 + sqrt(1 - a^2) + x) - 
            (-3 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)^3) + 
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(x*(2*sqrt(1 - a^2) + x)*lambda - 
            (a*m - (a^2 + (1 + sqrt(1 - a^2) + x)^2)*omega)^2))))
        end
    elseif s == 1
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-((2*I*a*x*(2*sqrt(1 - a^2) + x)*(-3*I*a^2*m*(1 + sqrt(1 - a^2) + x) - 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda)))*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + 
            omega))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-2*I*a^3*m*(1 + sqrt(1 - a^2) + x) - 2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*(2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda))))) - 
            (-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            (2*x*(2*sqrt(1 - a^2) + x)*(a^4 - a^2*(1 + sqrt(1 - a^2) + x) - 
            (-2 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)^3) + 
            ((1 + sqrt(1 - a^2) + x)^2*(-1 + 2*sqrt(1 - a^2) + 2*x) + 
            a^2*(1 + 2*sqrt(1 - a^2) + 2*x))^2 + (2*a*x*(2*sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^2*
            (-1 + 2*sqrt(1 - a^2) + 2*x) + a^2*(1 + 2*sqrt(1 - a^2) + 2*x))*
            (-3*I*a^2*m*(1 + sqrt(1 - a^2) + x) - I*m*(1 + sqrt(1 - a^2) + x)^3 + 
            2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*(3 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*(-2*I*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (2*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(-lambda - I*(1 + sqrt(1 - a^2) + x)*omega) + 
            a^8*(1 + lambda)*(2 + (1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            2*I*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-6 + (1 + sqrt(1 - a^2) + x)*
            (4 + m^2 - lambda - 4*I*omega) - 3*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2) - 2*I*a*m*(1 + sqrt(1 - a^2) + x)^5*
            (3 + (1 + sqrt(1 - a^2) + x)*(-5 + 2*lambda) + (1 + sqrt(1 - a^2) + x)^2*
            (2 - lambda + 2*I*omega) - I*(1 + sqrt(1 - a^2) + x)^3*(3 + lambda)*omega + 
            (1 + sqrt(1 - a^2) + x)^4*omega^2) - 2*I*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (13 + (1 + sqrt(1 - a^2) + x)*(-17 + 2*lambda) + (1 + sqrt(1 - a^2) + x)^2*
            (6 + m^2 - 2*lambda - 2*I*omega) - 3*I*(1 + sqrt(1 - a^2) + x)^3*(2 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^4*omega^2) + (1 + sqrt(1 - a^2) + x)^6*
            (-3*(2 + lambda) - (1 + sqrt(1 - a^2) + x)^2*lambda*(2 + lambda) + 
            2*(1 + sqrt(1 - a^2) + x)*(2 + 3*lambda + lambda^2) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda)*omega^2) + a^6*(1 + sqrt(1 - a^2) + x)*(-16*(1 + lambda) + 
            (1 + sqrt(1 - a^2) + x)*(5 + m^2*(-1 + lambda) + 6*lambda - lambda^2 + 2*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(5 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*
            (I + 2*I*m^2 + omega)) + a^4*(1 + sqrt(1 - a^2) + x)^2*
            (11 + 25*lambda + 2*(1 + sqrt(1 - a^2) + x)*(-5 + 3*m^2 - 13*lambda + lambda^2 - 2*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^2*(5 + 4*lambda - 3*lambda^2 + m^2*(3 + 2*lambda) - 8*I*omega) + 
            3*(1 + sqrt(1 - a^2) + x)^4*(3 + 2*lambda)*omega^2 + 4*(1 + sqrt(1 - a^2) + x)^3*omega*
            (I + 2*I*m^2 + omega)) + a^2*(1 + sqrt(1 - a^2) + x)^3*
            (30 + 5*(1 + sqrt(1 - a^2) + x)*(-7 + 2*lambda) - 4*(1 + sqrt(1 - a^2) + x)^2*
            (-2 + m^2 + lambda - lambda^2 - 3*I*omega) + (1 + sqrt(1 - a^2) + x)^3*
            (2 - 2*lambda - 3*lambda^2 + m^2*(4 + lambda) - 10*I*omega) + (1 + sqrt(1 - a^2) + x)^5*
            (7 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^4*omega*(I + 2*I*m^2 + omega))))/
            ((1 + sqrt(1 - a^2) + x)^2*(-2*I*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))))))
        end
    elseif s == -1
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-((2*I*a*x*(2*sqrt(1 - a^2) + x)*(3*I*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + 
            omega))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (2*I*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))) - 
            (-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            ((a^2 - (1 + sqrt(1 - a^2) + x)^2)^2 + 2*x*(1 + sqrt(1 - a^2) + x)*
            (2*sqrt(1 - a^2) + x)*(-3*a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(a^2 - (1 + sqrt(1 - a^2) + x)^2)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(3*I*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            I*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (2*I*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) - 
            ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(2*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*
            (2 - lambda + I*(1 + sqrt(1 - a^2) + x)*omega) + a^8*(-1 + lambda)*
            (2 + (1 + sqrt(1 - a^2) + x)^2*omega^2) + 2*I*a^5*m*(1 + sqrt(1 - a^2) + x)^3*
            (-2 + m^2 - lambda + 4*I*omega - 3*I*(1 + sqrt(1 - a^2) + x)*omega + 
            3*I*(1 + sqrt(1 - a^2) + x)*lambda*omega + 3*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 
            (1 + sqrt(1 - a^2) + x)^6*lambda*(-3 - (1 + sqrt(1 - a^2) + x)^2*lambda + 
            2*(1 + sqrt(1 - a^2) + x)*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            2*I*a*m*(1 + sqrt(1 - a^2) + x)^5*(sqrt(1 - a^2) + x + 
            2*(1 + sqrt(1 - a^2) + x)*lambda - (1 + sqrt(1 - a^2) + x)^2*(lambda + 2*I*omega) + 
            I*(1 + sqrt(1 - a^2) + x)^3*(1 + lambda)*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            2*I*a^3*m*(1 + sqrt(1 - a^2) + x)^3*(1 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda) + 
            (1 + sqrt(1 - a^2) + x)^2*(-2 + m^2 - 2*lambda + 2*I*omega) + 
            3*I*(1 + sqrt(1 - a^2) + x)^3*lambda*omega + 3*(1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            a^6*(1 + sqrt(1 - a^2) + x)*(8 - 8*lambda + (1 + sqrt(1 - a^2) + x)*
            (-1 + m^2*(-3 + lambda) - lambda^2 - 2*I*omega) + (1 + sqrt(1 - a^2) + x)^3*(-3 + 4*lambda)*
            omega^2 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(-I - 2*I*m^2 + omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^2*(-11 + 9*lambda + 2*(1 + sqrt(1 - a^2) + x)*
            (1 + 3*m^2 + lambda + lambda^2 + 2*I*omega) + (1 + sqrt(1 - a^2) + x)^2*
            (1 - 6*lambda - 3*lambda^2 + m^2*(-1 + 2*lambda) + 8*I*omega) + 3*(1 + sqrt(1 - a^2) + x)^4*
            (-1 + 2*lambda)*omega^2 + 4*(1 + sqrt(1 - a^2) + x)^3*omega*(-I - 2*I*m^2 + omega)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^3*(6 - 3*(1 + sqrt(1 - a^2) + x)*(1 + 2*lambda) - 
            4*(1 + sqrt(1 - a^2) + x)^2*(m^2 - 3*lambda - lambda^2 + 3*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(-4*lambda - 3*lambda^2 + m^2*(2 + lambda) + 10*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^5*(-1 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^4*omega*
            (-I - 2*I*m^2 + omega))))/((1 + sqrt(1 - a^2) + x)^2*
            (2*I*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*I*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))))))
        end
    elseif s == 2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 + 
            (8*I*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*
            (I*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*I + 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - 
            I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega) + 
            a^3*(-6 + 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(24 + 10*lambda + lambda^2 + 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(-6*I - 2*I*(1 + sqrt(1 - a^2) + x)*(4 + lambda) + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) - 12*a^4*(-1 + 2*I*(1 + sqrt(1 - a^2) + x)*omega + 
            2*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - (((1 + sqrt(1 - a^2) + x)^2*(-1 + 3*sqrt(1 - a^2) + 3*x) + 
            a^2*(1 + 3*sqrt(1 - a^2) + 3*x))^2 + x*(2*sqrt(1 - a^2) + x)*
            (3*a^4 + (1 + sqrt(1 - a^2) + x)^3*(8 - 3*(1 + sqrt(1 - a^2) + x))) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^2*(-1 + 3*sqrt(1 - a^2) + 3*x) + 
            a^2*(1 + 3*sqrt(1 - a^2) + 3*x))*(I*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I + 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - 
            I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega) + 
            a^3*(-6 + 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*((-(1 + sqrt(1 - a^2) + x)^4)*(24 + 10*lambda + lambda^2 + 
            12*I*omega) - 24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I + 2*(1 + sqrt(1 - a^2) + x)*
            omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*
            (4 + lambda) - 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 + 2*I*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*
            omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 - 6*I*omega) - 2*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(3 + 3*I*(1 + sqrt(1 - a^2) + x)*omega - 
            4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - (1 + sqrt(1 - a^2) + x)^6*
            (24 + 10*lambda + lambda^2 + 12*I*omega)*(-12 - (1 + sqrt(1 - a^2) + x)^2*lambda + 
            2*(1 + sqrt(1 - a^2) + x)*(4 + lambda) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-54*I - 2*I*(1 + sqrt(1 - a^2) + x)*
            (-55 + 3*m^2 - 4*lambda - 18*I*omega) + 2*(1 + sqrt(1 - a^2) + x)^2*
            (4 + 12*m^2 - 8*lambda + 27*I*omega)*omega - 6*I*(1 + sqrt(1 - a^2) + x)^3*(8 + lambda)*
            omega^2 + 45*(1 + sqrt(1 - a^2) + x)^4*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 - 2*I*(1 + sqrt(1 - a^2) + x)^3*
            omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 4*a^6*(1 + sqrt(1 - a^2) + x)*
            (60 - 3*(1 + sqrt(1 - a^2) + x)*(16 + 3*m^2 + 18*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(1 + 36*m^2 - 8*lambda + 18*I*omega)*omega^2 - 
            2*I*(1 + sqrt(1 - a^2) + x)^4*(7 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 - 2*(1 + sqrt(1 - a^2) + x)^2*omega*(-40*I + 9*I*m^2 - 4*I*lambda + 9*omega)) + 
            2*a*m*(1 + sqrt(1 - a^2) + x)^4*(144*I + 8*I*(1 + sqrt(1 - a^2) + x)*
            (-38 + lambda) + (1 + sqrt(1 - a^2) + x)^4*(40 + 20*lambda + lambda^2 + 24*I*omega)*omega + 
            4*I*(1 + sqrt(1 - a^2) + x)^5*(4 + lambda)*omega^2 - 6*(1 + sqrt(1 - a^2) + x)^6*
            omega^3 + 4*(1 + sqrt(1 - a^2) + x)^2*(16*I + 3*I*lambda + 2*I*lambda^2 + 6*omega) - 
            4*I*(1 + sqrt(1 - a^2) + x)^3*(-12 + lambda + lambda^2 - 14*I*omega - 5*I*lambda*omega)) - 
            2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*(312*I + (1 + sqrt(1 - a^2) + x)^3*
            (84 + 30*m^2 - 52*lambda - lambda^2 + 36*I*omega)*omega - 4*I*(1 + sqrt(1 - a^2) + x)^4*
            (19 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3 + 
            12*(1 + sqrt(1 - a^2) + x)*(-63*I + 3*I*m^2 - 3*I*lambda + 32*omega) + 
            4*(1 + sqrt(1 - a^2) + x)^2*(71*I + I*lambda^2 - I*m^2*(10 + lambda) - 104*omega + 
            lambda*(9*I + 16*omega))) + a^2*(1 + sqrt(1 - a^2) + x)^3*
            (576 + 96*(1 + sqrt(1 - a^2) + x)*(-12 + 4*m^2 - 3*I*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^2*(18 - 3*lambda^2 + m^2*(-58 + 8*lambda) + 
            lambda*(-30 - 2*I*omega) + 46*I*omega) - 8*I*(1 + sqrt(1 - a^2) + x)^4*
            (-lambda^2 + 2*m^2*(7 + lambda) + lambda*(2 + 5*I*omega) + 11*I*omega)*omega + 
            2*(1 + sqrt(1 - a^2) + x)^5*(-34 + 24*m^2 - 20*lambda - lambda^2 - 24*I*omega)*omega^2 - 
            8*I*(1 + sqrt(1 - a^2) + x)^6*(1 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^7*
            omega^4 + (1 + sqrt(1 - a^2) + x)^3*(lambda^3 + 36*lambda*(4 + I*omega) + 
            2*lambda^2*(11 - 8*I*omega) - m^2*(-136 + 42*lambda + lambda^2 - 36*I*omega) - 
            8*(-18 + 19*I*omega + 6*omega^2))) + a^4*(1 + sqrt(1 - a^2) + x)^2*
            (-672 + (1 + sqrt(1 - a^2) + x)*(960 - 72*m^2 + 624*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^4*(44 + 180*m^2 - 62*lambda - lambda^2 + 36*I*omega)*omega^2 - 
            8*I*(1 + sqrt(1 - a^2) + x)^5*(5 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^6*
            omega^4 + 8*(1 + sqrt(1 - a^2) + x)^3*omega*(56*I + I*lambda^2 - 3*I*m^2*(9 + lambda) - 
            46*omega + lambda*(6*I + 8*omega)) + 4*(1 + sqrt(1 - a^2) + x)^2*
            (6*m^4 + m^2*(7 - 8*lambda + 54*I*omega) + 2*(-15 + 10*lambda + lambda^2 - 144*I*omega - 9*I*
            lambda*omega + 48*omega^2)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*(4 + lambda) - 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + 12*a^4*(-1 + 2*I*(1 + sqrt(1 - a^2) + x)*
            omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))))/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4))
        end
    elseif s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 + 
            (8*I*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*
            ((-I)*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*lambda) + 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 
            24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            ((a^2*(-1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^3)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(-a^4 - 8*a^2*(1 + sqrt(1 - a^2) + x) + 
            (1 + sqrt(1 - a^2) + x)^4) + (8*a*x*(2*sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(a^2*(-1 + sqrt(1 - a^2) + x) + 
            (1 + sqrt(1 - a^2) + x)^3)*((-I)*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 
            24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(-3 + 3*I*(1 + sqrt(1 - a^2) + x)*omega + 
            4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (72*I - 4*I*(1 + sqrt(1 - a^2) + x)*(15 + 9*m^2 - 13*lambda) + 
            4*(1 + sqrt(1 - a^2) + x)^2*(I*m^2*(6 + lambda) - I*(3 + lambda^2 + lambda*
            (9 + 16*I*omega) + 24*I*omega)) + (1 + sqrt(1 - a^2) + x)^3*
            (-36 + 30*m^2 - 44*lambda - lambda^2 - 36*I*omega)*omega + 4*I*(1 + sqrt(1 - a^2) + x)^4*
            (3 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 + 2*I*(1 + sqrt(1 - a^2) + x)^3*
            omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 4*a^6*(1 + sqrt(1 - a^2) + x)*
            (12 - 3*(1 + sqrt(1 - a^2) + x)*(-4 + 3*m^2 + 6*I*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(-39 + 36*m^2 - 8*lambda - 18*I*omega)*omega^2 + 
            2*I*(1 + sqrt(1 - a^2) + x)^4*(3 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(24*I + 9*I*m^2 - 4*I*lambda + 15*omega)) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-18*I + 2*(1 + sqrt(1 - a^2) + x)^2*
            (-36 + 12*m^2 - 8*lambda - 27*I*omega)*omega + 6*I*(1 + sqrt(1 - a^2) + x)^3*(4 + lambda)*
            omega^2 + 45*(1 + sqrt(1 - a^2) + x)^4*omega^3 + 2*(1 + sqrt(1 - a^2) + x)*
            (9*I + 3*I*m^2 - 4*I*lambda + 30*omega)) - 2*a*m*(1 + sqrt(1 - a^2) + x)^4*
            (-48*I - 24*I*(1 + sqrt(1 - a^2) + x)*(-2 + lambda) + 
            4*I*(1 + sqrt(1 - a^2) + x)^2*(7*lambda + 2*lambda^2 + 6*I*omega) - 
            (1 + sqrt(1 - a^2) + x)^4*(12*lambda + lambda^2 - 24*I*omega)*omega + 
            4*I*(1 + sqrt(1 - a^2) + x)^5*lambda*omega^2 + 6*(1 + sqrt(1 - a^2) + x)^6*omega^3 - 
            4*I*(1 + sqrt(1 - a^2) + x)^3*(lambda + lambda^2 + 6*I*omega + 5*I*lambda*omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^3*(24*m^4*(1 + sqrt(1 - a^2) + x) + 
            48*(-4 + 3*I*omega) + 8*(1 + sqrt(1 - a^2) + x)*(9 + lambda^2 + lambda*(2 + 13*I*omega) - 
            72*I*omega) - (1 + sqrt(1 - a^2) + x)^3*(60 + 54*lambda + lambda^2 + 36*I*omega)*omega^2 + 
            8*I*(1 + sqrt(1 - a^2) + x)^4*(-3 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 + 8*(1 + sqrt(1 - a^2) + x)^2*omega*(24*I - 6*I*lambda - I*lambda^2 + 18*omega + 
            8*lambda*omega) + 4*m^2*(30 - (1 + sqrt(1 - a^2) + x)*(33 + 8*lambda + 54*I*omega) + 
            6*I*(1 + sqrt(1 - a^2) + x)^2*(5 + lambda)*omega + 45*(1 + sqrt(1 - a^2) + x)^3*
            omega^2)) + (1 + sqrt(1 - a^2) + x)^6*((-1 + sqrt(1 - a^2) + x)*
            (1 + sqrt(1 - a^2) + x)*lambda^3 + lambda^2*(12 - 12*(1 + sqrt(1 - a^2) + x) + 
            2*(1 + sqrt(1 - a^2) + x)^2 - (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            12*I*omega*(-12 + 8*(1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            2*lambda*(-12 + 4*(1 + sqrt(1 - a^2) + x)*(2 - 3*I*omega) + 
            6*I*(1 + sqrt(1 - a^2) + x)^2*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^4*(96 + 8*(1 + sqrt(1 - a^2) + x)*
            (m^2*(6 + 8*lambda) - 3*(2 + lambda^2 + lambda*(2 + 2*I*omega) - 22*I*omega)) - 96*I*omega + 
            2*(1 + sqrt(1 - a^2) + x)^4*(6 + 24*m^2 - 12*lambda - lambda^2 + 24*I*omega)*omega^2 + 
            8*I*(1 + sqrt(1 - a^2) + x)^5*(-3 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^6*
            omega^4 + (1 + sqrt(1 - a^2) + x)^2*(lambda^3 + lambda*(24 - 34*m^2 - 4*I*omega) + 
            lambda^2*(14 - m^2 + 16*I*omega) - 12*I*(22 + 3*m^2 - 4*I*omega)*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^3*omega*((-I)*lambda^2 + 2*I*m^2*(3 + lambda) + 3*omega + 
            lambda*(2*I + 5*omega)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I - 2*(1 + sqrt(1 - a^2) + x)*omega) - 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + 12*a^4*(-1 - 2*I*(1 + sqrt(1 - a^2) + x)*
            omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

# Cache mechanism for the ingoing coefficients at horizon
# Initialize the cache with a set of fiducial parameters
_cached_ingoing_coefficients_at_hor_params::NamedTuple{(:s, :m, :a, :omega, :lambda), Tuple{Int, Int, _DEFAULTDATATYPE, _DEFAULTDATATYPE, _DEFAULTDATATYPE}} = (s=-2, m=2, a=0, omega=0.5, lambda=1)
_cached_ingoing_coefficients_at_hor::NamedTuple{(:expansion_coeffs, :Pcoeffs, :Qcoeffs), Tuple{Vector{_DEFAULTDATATYPE}, Vector{_DEFAULTDATATYPE}, Vector{_DEFAULTDATATYPE}}} = (
    expansion_coeffs = [_DEFAULTDATATYPE(1.0)], 
    Pcoeffs = [_DEFAULTDATATYPE(0.0)],
    Qcoeffs = [_DEFAULTDATATYPE(0.0)]
)

function ingoing_coefficient_at_hor(s::Int, m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    global _cached_ingoing_coefficients_at_hor_params
    global _cached_ingoing_coefficients_at_hor

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    _this_params = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    # Check if we can use the cached results
    if _cached_ingoing_coefficients_at_hor_params == _this_params
        expansion_coeffs = _cached_ingoing_coefficients_at_hor.expansion_coeffs
        Pcoeffs = _cached_ingoing_coefficients_at_hor.Pcoeffs
        Qcoeffs = _cached_ingoing_coefficients_at_hor.Qcoeffs
    else
        # Cannot re-use the cached results, re-compute from zero
        expansion_coeffs = [data_type(1.0)] # order 0
        Pcoeffs = [data_type(PminusH(s, m, a, omega, lambda, 0))] # order 0
        Qcoeffs = [data_type(0.0)] # order 0
    end

    if order > 0
        # Compute series expansion coefficients for P and Q
        _P(x) = PminusH(s, m, a, omega, lambda, x)
        _Q(x) = QminusH(s, m, a, omega, lambda, x)
        _P_taylor = taylor_expand(_P, 0, order=order)
        _Q_taylor = taylor_expand(_Q, 0, order=order)

        for i in length(Pcoeffs):order
            append!(Pcoeffs, getcoeff(_P_taylor, i))
            append!(Qcoeffs, getcoeff(_Q_taylor, i))       
        end
    end
    # Define the indicial polynomial
    indicial(nu) = nu*(nu - 1) + Pcoeffs[1]*nu + Qcoeffs[1]

    if order > 0
        # Evaluate the C coefficients
        for i in length(expansion_coeffs):order
            sum = 0.0
            for r in 0:i-1
                sum += expansion_coeffs[r+1]*(r*Pcoeffs[i-r+1] + Qcoeffs[i-r+1])
            end
            append!(expansion_coeffs, -sum/indicial(i))
        end
    end

    # Update cache
    _cached_ingoing_coefficients_at_hor_params = _this_params
    _cached_ingoing_coefficients_at_hor = (
        expansion_coeffs = expansion_coeffs,
        Pcoeffs = Pcoeffs,
        Qcoeffs = Qcoeffs
    )
    return expansion_coeffs[order+1]
end

end