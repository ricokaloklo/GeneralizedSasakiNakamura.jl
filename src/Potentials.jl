module Potentials

using ..Kerr

export sF, sU, VT

const I = 1im # Mathematica being Mathematica

function sF(s::Int, m::Int, a, omega, lambda, r)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    if s == 0
        # s = 0
        return 0.0
    elseif s == 1
        # s = +1
        return begin
            (2*a*(a^2 + (-2 + r)*r)*((-I)*a^4*m + a*(-3 + r)*r^2 + I*m*r^4 + a^3*(1 + r)))/
            ((a^2 + r^2)^2*(-2*I*a^3*m*r - 2*I*a*m*r^3 + a^4*(1 + lambda) + r^4*(2 + lambda) + 
            a^2*r*(2 + r*(3 + 2*lambda))))
        end
    elseif s == -1
        # s = -1
        return begin
            (2*a*(a^2 + (-2 + r)*r)*(I*a^4*m + a*(-3 + r)*r^2 - I*m*r^4 + a^3*(1 + r)))/
            ((-a)*(a^2 + r^2)^2*(a^3 - 2*I*a^2*m*r + a*(-2 + r)*r - 2*I*m*r^3) + (a^2 + r^2)^4*lambda)
        end
    elseif s == 2
        # s = +2
        return begin
            (8*a*(a^2 + (-2 + r)*r)*(I*m*r^3*(6 + r*(4 + lambda)) + 2*a^3*(3 + lambda + r*(3 - 5*m^2 + lambda - 2*I*omega)) - 2*I*a^5*omega + 
            a^4*m*((-I)*lambda + 4*r*omega) + 2*a*r^2*(-3*(3 + lambda) + r*(3 + 3*m^2 + lambda + 6*I*omega) + I*r^2*omega) + 
            2*a^2*m*r*(-I - 6*I*r + 2*r^2*omega)))/((a^2 + r^2)^2*(r^4*(24 + 10*lambda + lambda^2 + 12*I*omega) - 4*a^5*m*omega - 
            4*a^6*omega^2 + 8*a^3*m*(-I - I*r*lambda + r^2*omega) + 4*a*m*r^2*(-6*I - 2*I*r*(4 + lambda) + 3*r^2*omega) - 
            2*a^2*r*(-8*(3 + lambda) + r*(-12 + 12*m^2 - 6*lambda - lambda^2 + 12*I*omega) + 8*I*r^2*omega + 2*r^3*omega^2) + 
            a^4*(8*m^2 + 2*lambda + lambda^2 - 4*omega*(I + 4*I*r + 2*r^2*omega))))
        end
    elseif s == -2
        # s = -2
        return begin
            (8*a*(a^2 + (-2 + r)*r)*((-I)*m*r^2*(6 + r*lambda) + 3*a^2*m*r*(3*I - 4*r*omega) + 
            a*r*(9 + r*(-3 + 6*m^2 + 6*I*omega) + I*r^2*(-3 + lambda)*omega) + a^3*(-6 - 9*I*r*omega + 6*r^2*omega^2)))/
            (r*(a^2 + r^2)*(r^4*(2*lambda + lambda^2 - 12*I*omega) + 24*a^3*m*r*(-I + 2*r*omega) + 
            4*a*m*r^2*(6*I + 2*I*r*lambda + 3*r^2*omega) + a^4*(12 + 24*I*r*omega - 24*r^2*omega^2) - 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*I*omega) + 2*I*r^2*(-3 + lambda)*omega + 3*r^3*omega^2)))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function sU(s::Int, m::Int, a, omega, lambda, r)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    if s == 0
        # s = 0
        return begin
            (r^2*(a^2 - 2*r + r^2)^2)/(a^2 + r^2)^4 + 
            ((a^2 - 2*r + r^2)*((r*(-2 + 2*r))/(a^2 + r^2)^2 - (4*r^2*(a^2 - 2*r + r^2))/(a^2 + r^2)^3 + 
            (a^2 - 2*r + r^2)/(a^2 + r^2)^2))/(a^2 + r^2) + 
            ((a^2 - 2*r + r^2)*(lambda - ((-a)*m + (a^2 + r^2)*omega)^2/(a^2 - 2*r + r^2)))/(a^2 + r^2)^2
        end
    elseif s == 1
        # s = +1
        return begin
            ((-a^12)*(1 + lambda)*omega^2 + 2*a^11*m*omega*(lambda + I*r*omega) + 
            r^9*(2 + lambda)*((-2 + r)*(2 + lambda) - r^3*omega^2) + 
            2*a^9*m*(I*(1 + m^2*r - r*lambda) + r*(4 + r*(3 + 5*lambda))*omega + 5*I*r^3*omega^2) + 
            2*a^7*m*r*(I*(-2 + r*(-3 + r*(-2 + 3*m^2 - 4*lambda) + 2*lambda)) + 2*r^2*(5 + r*(6 + 5*lambda))*omega + 
            10*I*r^4*omega^2) + 2*a^5*m*r^3*(I*(8 + r*(-5 + 6*lambda + 3*r*(m^2 - 2*(1 + lambda)))) + 
            2*r^2*(3 + r*(9 + 5*lambda))*omega + 10*I*r^4*omega^2) + 
            a^10*((1 + lambda)^2 - 2*I*(1 + r)*omega - r*(2 + r*(7 + 6*lambda))*omega^2 - m^2*(-1 + lambda + 4*I*r*omega)) + 
            2*a*m*r^8*(2 + lambda + I*r*omega)*(2*I + r*(-I + r*omega)) + 
            a^8*r*((1 + lambda)*(2 - 2*lambda + r*(7 + 5*lambda)) + 4*I*(1 + r - 2*r^2)*omega - 
            r^2*(8 + 5*r*(4 + 3*lambda))*omega^2 + m^2*(-6 + r*(-1 - 4*lambda - 16*I*r*omega))) + 
            a^2*r^6*(-8*(2 + lambda) + r*((2 + lambda)*(-8*(1 + lambda) + r*(8 + 5*lambda)) - 2*I*(-3 + r)*(-2 + r)*omega - 
            r^2*(2 + r*(11 + 6*lambda))*omega^2 + m^2*(4 - r*(4 + lambda) - 4*I*r^2*omega))) + 
            2*I*a^3*m*r^5*(10 + r*(3 + 6*lambda + r*(-6 + m^2 - 4*lambda - I*(-2 + r*(12 + 5*lambda))*omega + 
            5*r^2*omega^2))) + a^6*r^2*(-4 - 8*lambda + r*(4 + 19*r - 8*lambda + 28*r*lambda - 8*lambda^2 + 
            10*r*lambda^2 - 4*I*(1 + 3*(-2 + r)*r)*omega - 2*r^2*(6 + 5*r*(3 + 2*lambda))*omega^2 + 
            m^2*(-8 + 3*r*(-3 - 2*lambda - 8*I*r*omega)))) + 
            a^4*r^3*(-8 + r*(-4*(5 + 4*lambda) + r*(-2*(3 + 6*lambda*(2 + lambda) + 10*I*omega) + 
            r*(25 + 2*lambda*(16 + 5*lambda) + 28*I*omega) - 5*r^3*(5 + 3*lambda)*omega^2 - 8*r^2*omega*(I + omega) + 
            m^2*(2 + r*(-11 - 4*lambda - 16*I*r*omega))))))/
            ((a^2 + r^2)^4*(-2*I*a^3*m*r - 2*I*a*m*r^3 + a^4*(1 + lambda) + r^4*(2 + lambda) + 
            a^2*r*(2 + r*(3 + 2*lambda))))
        end
    elseif s == -1
        # s = -1
        return begin
            ((-a^12)*(-1 + lambda)*omega^2 + 2*a^11*m*omega*(-2 + lambda - I*r*omega) + 
            2*a*m*r^8*(lambda - I*r*omega)*(I*(-2 + r) + r^2*omega) + r^9*lambda*((-2 + r)*lambda - r^3*omega^2) + 
            a^10*((-1 + lambda)^2 + 2*I*(1 + r)*omega + r*(-2 + 5*r - 6*r*lambda)*omega^2 + m^2*(3 - lambda + 4*I*r*omega)) - 
            2*I*a^9*m*(1 + r*(2 + m^2 - lambda + I*(4 - 7*r + 5*r*lambda)*omega + 5*r^2*omega^2)) - 
            2*I*a^5*m*r^3*(8 + r*(-17 + 3*r*(2 + m^2 - 2*lambda) + 6*lambda + 2*I*r*(3 + r*(-1 + 5*lambda))*omega + 
            10*r^3*omega^2)) + a^2*r^6*(lambda*(-2*(-2 + r)^2 + r*(-8 + 5*r)*lambda) + 
            2*I*(-3 + r)*(-2 + r)*r*omega + r^3*(-2 + r - 6*r*lambda)*omega^2 + 
            m^2*r*(4 + r*(-2 - lambda + 4*I*r*omega))) + a^8*r*((-1 + lambda)*(6 - 3*r - 2*lambda + 5*r*lambda) + 
            4*I*(-1 + r)*(1 + 2*r)*omega + r^2*(-8 + 5*r*(2 - 3*lambda))*omega^2 + 
            m^2*(-6 + r*(7 - 4*lambda + 16*I*r*omega))) - 2*I*a^3*m*r^5*
            (10 + r*(-9 + 6*lambda + r*(2 + m^2 - 4*lambda + I*(-2 + r*(2 + 5*lambda))*omega + 5*r^2*omega^2))) - 
            2*I*a^7*m*r*(-2 + r*(-7 + 2*lambda + r*(6 + 3*m^2 - 4*lambda + 2*I*(5 - 4*r + 5*r*lambda)*omega + 
            10*r^2*omega^2))) + a^6*r^2*(12 - 8*lambda + 
            r*(3*(-4 + r + 8*lambda) + 2*lambda*(-4*lambda + r*(-6 + 5*lambda)) + 4*I*omega + 12*I*(-2 + r)*r*omega - 
            2*r^2*(6 + 5*r*(-1 + 2*lambda))*omega^2 + m^2*(-8 + 3*r*(1 - 2*lambda + 8*I*r*omega)))) + 
            a^4*r^3*(-8 + r*(12 - 16*lambda + r*(-6 + r + 24*lambda - 8*r*lambda - 12*lambda^2 + 10*r*lambda^2 + 
            4*I*(-1 + r)*(-5 + 2*r)*omega + r^2*(-8 + 5*r - 15*r*lambda)*omega^2 + 
            m^2*(2 + r*(-3 - 4*lambda + 16*I*r*omega))))))/
            ((-a)*(a^2 + r^2)^4*(a^3 - 2*I*a^2*m*r + a*(-2 + r)*r - 2*I*m*r^3) + (a^2 + r^2)^6*lambda)
        end
    elseif s == 2
        # s = +2
        return begin
            (4*a^13*m*omega^3 - 4*a^14*omega^4 + a^12*omega^2*(12 + 12*m^2 + 6*lambda + lambda^2 - 4*I*(1 + 4*r)*omega - 24*r^2*omega^2) - 
            2*a^11*m*omega*(18 + 10*m^2 - 4*lambda + lambda^2 - 32*I*r*omega + 4*I*r*lambda*omega - 16*r^2*omega^2) + 
            r^8*(24 + 10*lambda + lambda^2 + 12*I*omega)*(-12 - r^2*(6 + lambda) + 2*r*(9 + lambda) + r^4*omega^2) - 
            2*I*a^9*m*(4*m^2*(1 + r*(-10 + lambda)) + 4*(15 + lambda) + 4*r^3*(-28 + 5*lambda)*omega^2 + 50*I*r^4*omega^3 + 
            r^2*omega*(4*I - 4*I*lambda - 5*I*lambda^2 + 24*omega) - 4*r*(-18 - 5*lambda + lambda^2 - 9*I*omega + 7*I*lambda*omega)) + 
            4*a^7*m*r*(4*I*(33 + lambda) - 2*I*r*(-45 + 25*m^2 - 42*lambda + 2*lambda^2 - 24*I*omega) + 
            r^3*(22 + 30*m^2 - 28*lambda - 5*lambda^2 - 48*I*omega)*omega - 4*I*r^4*(-16 + 5*lambda)*omega^2 + 40*r^5*omega^3 + 
            2*I*r^2*(-42 + m^2*(10 - 3*lambda) - 5*lambda + 4*lambda^2 + 64*I*omega + 16*I*lambda*omega)) + 
            2*a^3*m*r^4*(-528*I - 8*I*r*(69 + 31*lambda) + 4*I*r^2*(51 + 9*m^2 - 2*lambda - 6*lambda^2 - 52*I*omega) + 
            r^4*(-170 + 30*m^2 - 76*lambda - 5*lambda^2 - 96*I*omega)*omega - 4*I*r^5*(8 + 5*lambda)*omega^2 + 32*r^6*omega^3 + 
            4*r^3*(42*I + 25*I*lambda + 4*I*lambda^2 - I*m^2*(10 + lambda) - 24*omega + 8*lambda*omega)) - 
            2*a^6*r^2*(20*m^4*r^2 - 36*(-12 - 2*lambda + lambda^2 - 4*I*omega) - 2*r^4*(84 + 40*lambda + 5*lambda^2 - 28*I*omega)*omega^2 + 
            80*I*r^5*omega^3 + 30*r^6*omega^4 - 4*r*(-108 - 48*lambda - 8*lambda^2 + lambda^3 - 48*I*omega - 4*I*lambda*omega) - 
            24*r^3*omega*(9*I + 4*I*lambda + 3*omega + lambda*omega) + r^2*(-36 + 78*lambda + 35*lambda^2 + 5*lambda^3 - 204*I*omega + 108*I*lambda*omega + 
            104*omega^2) + m^2*(-112 - 16*r*(23 + lambda) - r^2*(-32 + 42*lambda + 3*lambda^2 + 252*I*omega) - 
            24*I*r^3*(-1 + 2*lambda)*omega + 120*r^4*omega^2)) + 
            2*a*m*r^6*(48*I + 24*I*r*(3 + lambda) - 4*I*r^2*(60 + 23*lambda + 2*lambda^2 - 6*I*omega) - 
            r^4*(76 + 20*lambda + lambda^2 + 24*I*omega)*omega - 4*I*r^5*(4 + lambda)*omega^2 + 6*r^6*omega^3 + 
            4*r^3*(24*I + I*lambda^2 + 29*omega + 5*lambda*(2*I + omega))) + 
            a^2*r^6*(8*r*lambda*(24 + 10*lambda + lambda^2 + 12*I*omega) + 48*(24 + 10*lambda + lambda^2 + 16*I*omega) + 
            2*r^4*(64 + 28*lambda + 3*lambda^2 + 12*I*omega)*omega^2 - 16*I*r^5*omega^3 - 4*r^6*omega^4 - 
            r^2*(360 + 65*lambda^2 + 5*lambda^3 + 2*lambda*(135 + 62*I*omega) + 516*I*omega - 240*omega^2) + 
            m^2*r*(-64*(4 + lambda) + r*(152 + 42*lambda + lambda^2 - 36*I*omega) + 16*I*r^2*(7 + lambda)*omega - 48*r^3*omega^2) + 
            8*r^3*omega*(18*I - 11*omega + lambda*(4*I + omega))) + a^10*(8*m^4 - lambda^3 + lambda^2*(-1 + 6*r^2*omega^2) + 
            lambda*(-6 + 4*I*(5 + 8*r)*omega + 8*r*(1 + 5*r)*omega^2) + 4*omega*(3*I + 34*r*omega + 10*r^2*(2 - I*omega)*omega - 
            20*I*r^3*omega^2 - 15*r^4*omega^3) + m^2*(-24 + lambda^2 - 4*I*(-3 + 32*r)*omega + 2*I*lambda*(7*I + 8*r*omega))) + 
            4*a^5*m*r^2*(-144*I - 12*I*r*(-9 + 11*lambda) - 2*I*r^2*(-171 + 15*m^2 - 64*lambda + 6*lambda^2 + 34*I*omega) + 
            r^4*(-36 + 40*m^2 - 52*lambda - 5*lambda^2 - 72*I*omega)*omega - 4*I*r^5*(-4 + 5*lambda)*omega^2 + 35*r^6*omega^3 - 
            2*I*r^3*(m^2*(10 + 3*lambda) - 3*(-2 + 5*lambda + 2*lambda^2 + 42*I*omega + 2*I*lambda*omega))) + 
            a^4*r^3*(576*(3 + lambda) - 4*r^2*(-3*lambda^2 - 3*lambda^3 + 4*m^2*(-31 + 5*lambda) + 10*lambda*(9 - 2*I*omega) + 
            36*(9 + 5*I*omega)) + r*(-928*m^2 + 12*(48 + 38*lambda + 11*lambda^2 - 44*I*omega)) + 
            r^5*(284 - 180*m^2 + 130*lambda + 15*lambda^2 - 28*I*omega)*omega^2 - 80*I*r^6*omega^3 - 24*r^7*omega^4 + 
            16*r^4*omega*(27*I + I*m^2*(13 + 4*lambda) - 8*omega + 2*lambda*(4*I + omega)) - 
            2*r^3*(12*m^4 + 50*lambda^2 + 5*lambda^3 + 4*lambda*(39 + 35*I*omega) - 2*m^2*(52 + 28*lambda + lambda^2 + 48*I*omega) + 
            4*(27 + 54*I*omega - 32*omega^2))) + a^8*r*(-8*m^4*r + 2*(-17*lambda^2 + lambda^3 + lambda*(-30 - 20*I*omega) + 60*I*omega) + 
            r^3*(224 + 110*lambda + 15*lambda^2 - 108*I*omega)*omega^2 - 160*I*r^4*omega^3 - 80*r^5*omega^4 - 
            4*m^2*(4 - 12*lambda - r*(-36 + lambda^2 + 72*I*omega) - 4*I*r^2*(-17 + 4*lambda)*omega + 30*r^3*omega^2) + 
            16*r^2*omega*(9*I + 20*omega + 2*lambda*(4*I + omega)) - r*(20*lambda^2 + 5*lambda^3 + 4*lambda*(9 + 7*I*omega) + 
            8*(-9 - 51*I*omega + 28*omega^2))))/((a^2 + r^2)^4*((-r^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) + 4*a^5*m*omega + 
            4*a^6*omega^2 + 4*a*m*r^2*(6*I + 2*I*r*(4 + lambda) - 3*r^2*omega) + 8*a^3*m*(I + I*r*lambda - r^2*omega) + 
            2*a^2*r*(-8*(3 + lambda) + r*(-12 + 12*m^2 - 6*lambda - lambda^2 + 12*I*omega) + 8*I*r^2*omega + 2*r^3*omega^2) - 
            a^4*(8*m^2 + 2*lambda + lambda^2 - 4*omega*(I + 4*I*r + 2*r^2*omega))))
        end
    elseif s == -2
        # s = -2
        return begin
            (1/(a^2 + r^2)^4)*((a^2*(-2 + r) + r^3)^2 + (a^2 + (-2 + r)*r)*(-a^4 - 8*a^2*r + r^4) + (8*a*(a^2 + (-2 + r)*r)*(a^2 + r^2)*(a^2*(-2 + r) + r^3)*((-I)*m*r^2*(6 + r*lambda) + 
            3*a^2*m*r*(3*I - 4*r*omega) + a*r*(9 + r*(-3 + 6*m^2 + 6*I*omega) + I*r^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*I*r*omega + 6*r^2*omega^2)))/(r*(r^4*(2*lambda + lambda^2 - 12*I*omega) + 24*a^3*m*r*(-I + 2*r*omega) + 
            4*a*m*r^2*(6*I + 2*I*r*lambda + 3*r^2*omega) + a^4*(12 + 24*I*r*omega - 24*r^2*omega^2) - 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*I*omega) + 2*I*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))) - ((a^2 + r^2)^2*(-24*a^7*m*r^2*omega*(-3 + 3*I*r*omega + 4*r^2*omega^2) + 
            12*a^8*(-2 - 3*r^2*omega^2 + 2*I*r^3*omega^3 + 2*r^4*omega^4) + 
            4*a^6*r*(12 - 3*r*(-4 + 3*m^2 + 6*I*omega) + r^3*(-39 + 36*m^2 - 8*lambda - 18*I*omega)*omega^2 + 
            2*I*r^4*(3 + lambda)*omega^3 + 15*r^5*omega^4 + 2*r^2*omega*(24*I + 9*I*m^2 - 4*I*lambda + 15*omega)) - 
            4*a^5*m*r^2*(-18*I + 2*r^2*(-36 + 12*m^2 - 8*lambda - 27*I*omega)*omega + 6*I*r^3*(4 + lambda)*omega^2 + 
            45*r^4*omega^3 + 2*r*(9*I + 3*I*m^2 - 4*I*lambda + 30*omega)) - 
            2*a*m*r^4*(-48*I - 24*I*r*(-2 + lambda) + 4*I*r^2*(7*lambda + 2*lambda^2 + 6*I*omega) - 
            r^4*(12*lambda + lambda^2 - 24*I*omega)*omega + 4*I*r^5*lambda*omega^2 + 6*r^6*omega^3 - 
            4*I*r^3*(lambda + lambda^2 + 6*I*omega + 5*I*lambda*omega)) + a^4*r^3*(24*m^4*r + 48*(-4 + 3*I*omega) - 
            r^3*(60 + 54*lambda + lambda^2 + 36*I*omega)*omega^2 + 8*I*r^4*(-3 + 2*lambda)*omega^3 + 48*r^5*omega^4 + 
            8*r*(9 + 2*lambda + lambda^2 - 72*I*omega + 13*I*lambda*omega) + 8*r^2*omega*(24*I - 6*I*lambda - I*lambda^2 + 18*omega + 
            8*lambda*omega) + 4*m^2*(30 - r*(33 + 8*lambda + 54*I*omega) + 6*I*r^2*(5 + lambda)*omega + 45*r^3*omega^2)) + 
            r^6*((-2 + r)*r*lambda^3 + lambda^2*(12 - 12*r + 2*r^2 - r^4*omega^2) + 12*I*omega*(-12 + 8*r + r^4*omega^2) - 
            2*lambda*(-12 + 4*r*(2 - 3*I*omega) + 6*I*r^2*omega + r^4*omega^2)) + 
            a^2*r^4*(96 - 96*I*omega + 2*r^4*(6 + 24*m^2 - 12*lambda - lambda^2 + 24*I*omega)*omega^2 + 
            8*I*r^5*(-3 + lambda)*omega^3 + 12*r^6*omega^4 + 8*r^3*omega*(2*I*lambda - I*lambda^2 + 2*I*m^2*(3 + lambda) + 3*omega + 
            5*lambda*omega) + r^2*(lambda^3 + lambda*(24 - 34*m^2 - 4*I*omega) + lambda^2*(14 - m^2 + 16*I*omega) - 
            12*I*(22 + 3*m^2 - 4*I*omega)*omega) + 8*r*(m^2*(6 + 8*lambda) - 3*(2 + 2*lambda + lambda^2 - 22*I*omega + 
            2*I*lambda*omega))) - 2*a^3*m*r^3*(72*I - 4*I*r*(15 + 9*m^2 - 13*lambda) + 
            r^3*(-36 + 30*m^2 - 44*lambda - lambda^2 - 36*I*omega)*omega + 4*I*r^4*(3 + 4*lambda)*omega^2 + 48*r^5*omega^3 + 
            4*r^2*(I*m^2*(6 + lambda) - I*(3 + 9*lambda + lambda^2 + 24*I*omega + 16*I*lambda*omega)))))/
            (r^2*((-r^4)*(2*lambda + lambda^2 - 12*I*omega) + 24*a^3*m*r*(I - 2*r*omega) - 
            4*a*m*r^2*(6*I + 2*I*r*lambda + 3*r^2*omega) + 12*a^4*(-1 - 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*I*omega) + 2*I*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function VT(s::Int, m::Int, a, omega, lambda, r)
    _K = K(m, a, omega, r)
    _Delta = Delta(a, r)
    return lambda - 4im*s*omega*r - (_K^2 - 2im*s*(r-1)*_K)/_Delta
end

end