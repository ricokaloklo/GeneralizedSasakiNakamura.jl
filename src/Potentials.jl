module Potentials

using ..Kerr
using ..Transformation

export sF, sU, VT, radial_Teukolsky_equation

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
            -((2*a*(a^2 + (-2 + r)*r)*(-3*I*a^2*m*r - I*m*r^3 + 2*a^3*(1 + lambda) + a*r*(3 + r*(3 + 2*lambda))))/
            (r*(a^2 + r^2)*(-2*I*a^3*m*r - 2*I*a*m*r^3 + a^4*(1 + lambda) + r^4*(2 + lambda) + 
            a^2*r*(2 + r*(3 + 2*lambda)))))
        end
    elseif s == -1
        # s = -1
        return begin
            -((2*a*(a^2 + (-2 + r)*r)*(3*I*a^2*m*r + I*m*r^3 + 2*a^3*(-1 + lambda) + a*r*(3 + r*(-1 + 2*lambda))))/
            (r*(a^2 + r^2)*(2*I*a^3*m*r + 2*I*a*m*r^3 + a^4*(-1 + lambda) + r^4*lambda + a^2*r*(2 + r*(-1 + 2*lambda)))))
        end
    elseif s == 2
        # s = +2
        return begin
            (8*a*(a^2 + (-2 + r)*r)*(I*m*r^2*(6 + r*(4 + lambda)) - 3*a^2*m*r*(3*I + 4*r*omega) + 
            a*r*(9 + r*(-3 + 6*m^2 - 6*I*omega) - I*r^2*(1 + lambda)*omega) + a^3*(-6 + 9*I*r*omega + 6*r^2*omega^2)))/
            (r*(a^2 + r^2)*(r^4*(24 + 10*lambda + lambda^2 + 12*I*omega) + 24*a^3*m*r*(I + 2*r*omega) + 
            4*a*m*r^2*(-6*I - 2*I*r*(4 + lambda) + 3*r^2*omega) - 12*a^4*(-1 + 2*I*r*omega + 2*r^2*omega^2) - 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*I*omega) - 2*I*r^2*(1 + lambda)*omega + 3*r^3*omega^2)))
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
            (1/(a^2 + r^2)^4)*(r^2*(a^2 + (-2 + r)*r)^2 + (a^2 + (-2 + r)*r)*(a^4 - 4*a^2*r - (-4 + r)*r^3) + 
            (a^2 + r^2)^2*((a^2 + (-2 + r)*r)*lambda - (a*m - (a^2 + r^2)*omega)^2))
        end
    elseif s == 1
        # s = +1
        return begin
            (1/(a^2 + r^2)^4)*(2*(a^2 + (-2 + r)*r)*(a^4 - a^2*r - (-3 + r)*r^3) + 
            (r^2*(-3 + 2*r) + a^2*(-1 + 2*r))^2 + (2*a*(a^2 + (-2 + r)*r)*(a^2 + r^2)*
            (r^2*(-3 + 2*r) + a^2*(-1 + 2*r))*(-3*I*a^2*m*r - I*m*r^3 + 2*a^3*(1 + lambda) + 
            a*r*(3 + r*(3 + 2*lambda))))/(r*(-2*I*a^3*m*r - 2*I*a*m*r^3 + a^4*(1 + lambda) + r^4*(2 + lambda) + 
            a^2*r*(2 + r*(3 + 2*lambda)))) - ((a^2 + r^2)^2*(2*a^7*m*r^2*omega*(-lambda - I*r*omega) + 
            a^8*(1 + lambda)*(2 + r^2*omega^2) - 2*I*a^5*m*r^2*(-6 + r*(4 + m^2 - lambda - 4*I*omega) - 
            3*I*r^2*(1 + lambda)*omega + 3*r^3*omega^2) - 2*I*a*m*r^5*(3 + r*(-5 + 2*lambda) + r^2*(2 - lambda + 2*I*omega) - 
            I*r^3*(3 + lambda)*omega + r^4*omega^2) - 2*I*a^3*m*r^3*(13 + r*(-17 + 2*lambda) + 
            r^2*(6 + m^2 - 2*lambda - 2*I*omega) - 3*I*r^3*(2 + lambda)*omega + 3*r^4*omega^2) + 
            r^6*(-3*(2 + lambda) - r^2*lambda*(2 + lambda) + 2*r*(2 + 3*lambda + lambda^2) + r^4*(2 + lambda)*omega^2) + 
            a^6*r*(-16*(1 + lambda) + r*(5 + m^2*(-1 + lambda) + 6*lambda - lambda^2 + 2*I*omega) + r^3*(5 + 4*lambda)*omega^2 + 
            2*r^2*omega*(I + 2*I*m^2 + omega)) + a^4*r^2*(11 + 25*lambda + 
            2*r*(-5 + 3*m^2 - 13*lambda + lambda^2 - 2*I*omega) + r^2*(5 + 4*lambda - 3*lambda^2 + m^2*(3 + 2*lambda) - 
            8*I*omega) + 3*r^4*(3 + 2*lambda)*omega^2 + 4*r^3*omega*(I + 2*I*m^2 + omega)) + 
            a^2*r^3*(30 + 5*r*(-7 + 2*lambda) - 4*r^2*(-2 + m^2 + lambda - lambda^2 - 3*I*omega) + 
            r^3*(2 - 2*lambda - 3*lambda^2 + m^2*(4 + lambda) - 10*I*omega) + r^5*(7 + 4*lambda)*omega^2 + 
            2*r^4*omega*(I + 2*I*m^2 + omega))))/(r^2*(-2*I*a^3*m*r - 2*I*a*m*r^3 + a^4*(1 + lambda) + 
            r^4*(2 + lambda) + a^2*r*(2 + r*(3 + 2*lambda)))))
        end
    elseif s == -1
        # s = -1
        return begin
            (1/(a^2 + r^2)^4)*((a^2 - r^2)^2 + 2*r*(a^2 + (-2 + r)*r)*(-3*a^2 + r^2) + 
            (2*a*(a^2 + (-2 + r)*r)*(a^2 - r^2)*(a^2 + r^2)*(3*I*a^2*m*r + I*m*r^3 + 2*a^3*(-1 + lambda) + 
            a*r*(3 + r*(-1 + 2*lambda))))/(r*(2*I*a^3*m*r + 2*I*a*m*r^3 + a^4*(-1 + lambda) + r^4*lambda + 
            a^2*r*(2 + r*(-1 + 2*lambda)))) - ((a^2 + r^2)^2*(2*a^7*m*r^2*omega*(2 - lambda + I*r*omega) + 
            a^8*(-1 + lambda)*(2 + r^2*omega^2) + 2*I*a^5*m*r^3*(-2 + m^2 - lambda + 4*I*omega - 3*I*r*omega + 
            3*I*r*lambda*omega + 3*r^2*omega^2) + r^6*lambda*(-3 - r^2*lambda + 2*r*(1 + lambda) + r^4*omega^2) + 
            2*I*a*m*r^5*(-1 + r + 2*r*lambda - r^2*(lambda + 2*I*omega) + I*r^3*(1 + lambda)*omega + r^4*omega^2) + 
            2*I*a^3*m*r^3*(1 + r*(3 + 2*lambda) + r^2*(-2 + m^2 - 2*lambda + 2*I*omega) + 3*I*r^3*lambda*omega + 
            3*r^4*omega^2) + a^6*r*(8 - 8*lambda + r*(-1 + m^2*(-3 + lambda) - lambda^2 - 2*I*omega) + 
            r^3*(-3 + 4*lambda)*omega^2 + 2*r^2*omega*(-I - 2*I*m^2 + omega)) + 
            a^4*r^2*(-11 + 9*lambda + 2*r*(1 + 3*m^2 + lambda + lambda^2 + 2*I*omega) + 
            r^2*(1 - 6*lambda - 3*lambda^2 + m^2*(-1 + 2*lambda) + 8*I*omega) + 3*r^4*(-1 + 2*lambda)*omega^2 + 
            4*r^3*omega*(-I - 2*I*m^2 + omega)) + a^2*r^3*(6 - 3*r*(1 + 2*lambda) - 
            4*r^2*(m^2 - 3*lambda - lambda^2 + 3*I*omega) + r^3*(-4*lambda - 3*lambda^2 + m^2*(2 + lambda) + 10*I*omega) + 
            r^5*(-1 + 4*lambda)*omega^2 + 2*r^4*omega*(-I - 2*I*m^2 + omega))))/
            (r^2*(2*I*a^3*m*r + 2*I*a*m*r^3 + a^4*(-1 + lambda) + r^4*lambda + a^2*r*(2 + r*(-1 + 2*lambda)))))
        end
    elseif s == 2
        # s = +2
        return begin
            (1/(a^2 + r^2)^4)*((a^2 + (-2 + r)*r)*(3*a^4 + (8 - 3*r)*r^3) + 
            (r^2*(-4 + 3*r) + a^2*(-2 + 3*r))^2 + (8*a*(a^2 + (-2 + r)*r)*(a^2 + r^2)*
            (r^2*(-4 + 3*r) + a^2*(-2 + 3*r))*(I*m*r^2*(6 + r*(4 + lambda)) - 3*a^2*m*r*(3*I + 4*r*omega) + 
            a*r*(9 + r*(-3 + 6*m^2 - 6*I*omega) - I*r^2*(1 + lambda)*omega) + a^3*(-6 + 9*I*r*omega + 6*r^2*omega^2)))/
            (r*((-r^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 24*a^3*m*r*(I + 2*r*omega) + 
            4*a*m*r^2*(6*I + 2*I*r*(4 + lambda) - 3*r^2*omega) + 12*a^4*(-1 + 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*I*omega) - 2*I*r^2*(1 + lambda)*omega + 3*r^3*omega^2))) - 
            ((a^2 + r^2)^2*(24*a^7*m*r^2*omega*(3 + 3*I*r*omega - 4*r^2*omega^2) - r^6*(24 + 10*lambda + lambda^2 + 12*I*omega)*
            (-12 - r^2*lambda + 2*r*(4 + lambda) + r^4*omega^2) - 4*a^5*m*r^2*
            (-54*I - 2*I*r*(-55 + 3*m^2 - 4*lambda - 18*I*omega) + 2*r^2*(4 + 12*m^2 - 8*lambda + 27*I*omega)*omega - 
            6*I*r^3*(8 + lambda)*omega^2 + 45*r^4*omega^3) + 12*a^8*(-2 - 3*r^2*omega^2 - 2*I*r^3*omega^3 + 
            2*r^4*omega^4) + 4*a^6*r*(60 - 3*r*(16 + 3*m^2 + 18*I*omega) + r^3*(1 + 36*m^2 - 8*lambda + 18*I*omega)*
            omega^2 - 2*I*r^4*(7 + lambda)*omega^3 + 15*r^5*omega^4 - 2*r^2*omega*(-40*I + 9*I*m^2 - 4*I*lambda + 9*omega)) + 
            2*a*m*r^4*(144*I + 8*I*r*(-38 + lambda) + r^4*(40 + 20*lambda + lambda^2 + 24*I*omega)*omega + 
            4*I*r^5*(4 + lambda)*omega^2 - 6*r^6*omega^3 + 4*r^2*(16*I + 3*I*lambda + 2*I*lambda^2 + 6*omega) - 
            4*I*r^3*(-12 + lambda + lambda^2 - 14*I*omega - 5*I*lambda*omega)) - 
            2*a^3*m*r^3*(312*I + r^3*(84 + 30*m^2 - 52*lambda - lambda^2 + 36*I*omega)*omega - 
            4*I*r^4*(19 + 4*lambda)*omega^2 + 48*r^5*omega^3 + 12*r*(-63*I + 3*I*m^2 - 3*I*lambda + 32*omega) + 
            4*r^2*(71*I + 9*I*lambda + I*lambda^2 - I*m^2*(10 + lambda) - 104*omega + 16*lambda*omega)) + 
            a^2*r^3*(576 + 96*r*(-12 + 4*m^2 - 3*I*omega) + 8*r^2*(18 - 3*lambda^2 + m^2*(-58 + 8*lambda) + 
            lambda*(-30 - 2*I*omega) + 46*I*omega) - 8*I*r^4*(-lambda^2 + 2*m^2*(7 + lambda) + lambda*(2 + 5*I*omega) + 
            11*I*omega)*omega + 2*r^5*(-34 + 24*m^2 - 20*lambda - lambda^2 - 24*I*omega)*omega^2 - 8*I*r^6*(1 + lambda)*omega^3 + 
            12*r^7*omega^4 + r^3*(lambda^3 + 36*lambda*(4 + I*omega) + 2*lambda^2*(11 - 8*I*omega) - 
            m^2*(-136 + 42*lambda + lambda^2 - 36*I*omega) - 8*(-18 + 19*I*omega + 6*omega^2))) + 
            a^4*r^2*(-672 + r*(960 - 72*m^2 + 624*I*omega) + r^4*(44 + 180*m^2 - 62*lambda - lambda^2 + 36*I*omega)*
            omega^2 - 8*I*r^5*(5 + 2*lambda)*omega^3 + 48*r^6*omega^4 + 8*r^3*omega*(56*I + 6*I*lambda + I*lambda^2 - 
            3*I*m^2*(9 + lambda) - 46*omega + 8*lambda*omega) + 4*r^2*(6*m^4 + m^2*(7 - 8*lambda + 54*I*omega) + 
            2*(-15 + 10*lambda + lambda^2 - 144*I*omega - 9*I*lambda*omega + 48*omega^2)))))/
            (r^2*((-r^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 24*a^3*m*r*(I + 2*r*omega) + 
            4*a*m*r^2*(6*I + 2*I*r*(4 + lambda) - 3*r^2*omega) + 12*a^4*(-1 + 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*I*omega) - 2*I*r^2*(1 + lambda)*omega + 3*r^3*omega^2))))
        end
    elseif s == -2
        # s = -2
        return begin
            (1/(a^2 + r^2)^4)*((a^2*(-2 + r) + r^3)^2 + (a^2 + (-2 + r)*r)*(-a^4 - 8*a^2*r + r^4) + 
            (8*a*(a^2 + (-2 + r)*r)*(a^2 + r^2)*(a^2*(-2 + r) + r^3)*((-I)*m*r^2*(6 + r*lambda) + 
            3*a^2*m*r*(3*I - 4*r*omega) + a*r*(9 + r*(-3 + 6*m^2 + 6*I*omega) + I*r^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*I*r*omega + 6*r^2*omega^2)))/(r*(r^4*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*r*(-I + 2*r*omega) + 4*a*m*r^2*(6*I + 2*I*r*lambda + 3*r^2*omega) + 
            a^4*(12 + 24*I*r*omega - 24*r^2*omega^2) - 4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))) - 
            ((a^2 + r^2)^2*(-24*a^7*m*r^2*omega*(-3 + 3*I*r*omega + 4*r^2*omega^2) + 
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
            12*I*(22 + 3*m^2 - 4*I*omega)*omega) + 8*r*(m^2*(6 + 8*lambda) - 
            3*(2 + 2*lambda + lambda^2 - 22*I*omega + 2*I*lambda*omega))) - 
            2*a^3*m*r^3*(72*I - 4*I*r*(15 + 9*m^2 - 13*lambda) + r^3*(-36 + 30*m^2 - 44*lambda - lambda^2 - 36*I*omega)*
            omega + 4*I*r^4*(3 + 4*lambda)*omega^2 + 48*r^5*omega^3 + 4*r^2*(I*m^2*(6 + lambda) - 
            I*(3 + 9*lambda + lambda^2 + 24*I*omega + 16*I*lambda*omega)))))/
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

function radial_Teukolsky_equation(s, m, a, omega, lambda, r, R, dRdr, d2Rdr2)
    Delta(a, r)*d2Rdr2(r) + (2*(s+1)*(r-1))*dRdr(r) - VT(s, m, a, omega, lambda, r)*R(r)
end

function VSN(s::Int, m::Int, a, omega, lambda, r)
    _eta(x) = eta(s, m, a, omega, lambda, x)
    _drdrstar(x) = Delta(a, x)/(x^2 + a^2)
    _d2rdrstar2(x) = (2*(x^2 - a^2))/((x^2 + a^2)^2)

    _detadrstar(x) = _drdrstar(x)*eta_prime(s, m, a, omega, lambda, x)
    _d2etadrstar2(x) = _drdrstar(x)*(eta_prime(s, m, a, omega, lambda, x)*_d2rdrstar2(x) + _drdrstar(x)*eta_primeprime(s, m, a, omega, lambda, x))

    _sU(x) = sU(s, m, a, omega, lambda, x)

    return begin
        _sU(r) - ( 0.5/_eta(r) * _d2etadrstar2(r) - 0.75/(_eta(r)^2)*(_detadrstar(r))^2 )
    end
end

end