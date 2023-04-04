module Solutions

using DifferentialEquations
using ..Kerr
using ..Transformation
using ..Coordinates
using ..Potentials
using ..AsymptoticExpansionCoefficients
using ..InitialConditions
using ..ConversionFactors

const I = 1im # Mathematica being Mathematica
_DEFAULTDATATYPE = Float64 # Double precision by default

function PhiRePhiIm_from_XXprime(X, Xprime)
    PhiRe = atan(imag(X), real(X))
    PhiIm = -0.5*log(real(X)^2 + imag(X)^2)

    PhiReprime = imag(Xprime/X)
    PhiImprime = -real(Xprime/X)

    return PhiRe, PhiIm, PhiReprime, PhiImprime
end

function XXprime_from_PhiRePhiImsoln(Xsoln)
    X(rs) = exp(1im*(Xsoln(rs)[1] + 1im*Xsoln(rs)[2]))
    Xprime(rs) = X(rs)*(1im*Xsoln(rs)[3] - Xsoln(rs)[4])

    return X, Xprime
end

function GSN_magn_phase_eqns!(du, u, p, rs)
    r = r_from_rstar(p.a, rs)
    _sF = sF(p.s, p.m, p.a, p.omega, p.lambda, r)
    _sU = sU(p.s, p.m, p.a, p.omega, p.lambda, r)
    
    #=
    We write X = exp(I*(PhiRe + I*PhiIm)) = exp(-PhiIm) exp(I*PhiRe)
    Therefore PhiRe is the phase and exp(-PhiIm) is the magnitude
    Substitute X in this form into the GSN equation will give
    two coupled *non-linear* ODE

    u[1] = PhiRe
    u[2] = PhiIm
    u[3] = PhiRe' = u[1]'
    u[4] = PhiIm' = u[2]'
    =#
    du[1] = u[3]
    du[2] = u[4]
    du[3] = 2*u[3]*u[4] + real(_sF)*u[3] - imag(_sF)*u[4] + imag(_sU)
    du[4] = -u[3]^2 + u[4]^2 + real(_sF)*u[4] + imag(_sF)*u[3] - real(_sU)
end

function solve_Xup(s::Int, m::Int, a, omega, lambda, rsin, rsout; dtype=_DEFAULTDATATYPE, odealgo=Vern9(), reltol=1e-12, abstol=1e-12)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsout, "rsout ($rsout) must be larger than rsin ($rsin)"))
    end
    # Initial conditions at rs = rsout, the outer boundary
    Xup_rsout, Xupprime_rsout = Xup_initialconditions(s, m, a, omega, lambda, rsout)
    # Convert initial conditions for Xup for PhiRe PhiIm
    PhiRe, PhiIm, PhiReprime, PhiImprime = PhiRePhiIm_from_XXprime(Xup_rsout, Xupprime_rsout)
    u0 = [dtype(PhiRe); dtype(PhiIm); dtype(PhiReprime); dtype(PhiImprime)]
    rsspan = (rsout, rsin) # Integrate from rsout to rsin *inward*
    p = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    odeprob = ODEProblem(GSN_magn_phase_eqns!, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
    return odesoln
end

function solve_Xin(s::Int, m::Int, a, omega, lambda, rsin, rsout; dtype=_DEFAULTDATATYPE, odealgo=Vern9(), reltol=1e-12, abstol=1e-12)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsin, "rsin ($rsin) must be smaller than rsout ($rsout)"))
    end
    # Initial conditions at rs = rsin, the inner boundary; this should be very close to EH
    Xin_rsin, Xinprime_rsin = Xin_initialconditions(s, m, a, omega, lambda, rsin)
    # Convert initial conditions for Xin for PhiRe PhiIm
    PhiRe, PhiIm, PhiReprime, PhiImprime = PhiRePhiIm_from_XXprime(Xin_rsin, Xinprime_rsin)
    u0 = [dtype(PhiRe); dtype(PhiIm); dtype(PhiReprime); dtype(PhiImprime)]
    rsspan = (rsin, rsout) # Integrate from rsin to rsout *outward*
    p = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    odeprob = ODEProblem(GSN_magn_phase_eqns!, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
    return odesoln
end

function Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(s, m, a, omega, lambda, r)
    #=
    Here we use explicit form for the conversion matrix to
    faciliate cancellations
    =#
    M11 = M12 = M21 = M22 = 0.0
    if s == 0
        M11 = 1/sqrt(a^2 + r^2)
        M21 = -(r/(a^2 + r^2)^(3/2))
        M22 = sqrt(a^2 + r^2)/(a^2 + (-2 + r)*r)
    elseif s == +1
        M11 = begin
            -((I*r*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            ((-a^3)*m*r - a*m*r^3 + r^5*omega + 
            a^4*(2*I + r*omega) + 2*a^2*r*(-2*I + I*r + 
            r^2*omega)))/(sqrt(a^2 + r^2)*
            sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*
            (-2*I*a^3*m*r - 2*I*a*m*r^3 + a^4*(1 + lambda) + 
            r^4*(2 + lambda) + a^2*r*(2 + r*(3 + 2*lambda)))))
        end
        M12 = begin
            (r^2*(a^2 + r^2)^(3/2)*
            sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r)))/
            (sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*
            (-2*I*a^3*m*r - 2*I*a*m*r^3 + a^4*(1 + lambda) + 
            r^4*(2 + lambda) + a^2*r*(2 + r*(3 + 2*lambda))))
        end
        M21 = begin
            -((r*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            (-2*a^7*m*(I + r*omega) + a^8*omega*(2*I + r*omega) + 
            a^5*m*r*(2*I - I*r - 6*r^2*omega) + 
            a*m*r^5*(-4*I + 3*I*r - 2*r^2*omega) - 
            2*a^3*m*r^3*(I - 2*I*r + 3*r^2*omega) + 
            r^6*(2*(2 + lambda) - r*(2 + lambda) - I*r^2*omega + 
            r^3*omega^2) + a^6*(-4 + r*(1 + m^2 - lambda - 
            2*I*omega) + 5*I*r^2*omega + 4*r^3*omega^2) + 
            a^2*r^3*(4 + 4*r*lambda + r^2*(-3 + m^2 - 3*lambda - 
            2*I*omega) - I*r^3*omega + 4*r^4*omega^2) + 
            a^4*r*(8 + 2*r*(-4 + lambda) + 
            r^2*(2*m^2 - 3*lambda - 4*I*omega) + 3*I*r^3*omega + 
            6*r^4*omega^2)))/(sqrt(a^2 + r^2)*
            ((a^2 + (-2 + r)*r)*(a^2 + r^2))^(3/2)*
            (-2*I*a^3*m*r - 2*I*a*m*r^3 + a^4*(1 + lambda) + 
            r^4*(2 + lambda) + a^2*r*(2 + r*(3 + 2*lambda)))))
        end
        M22 = begin
            -((I*r^2*((a^2 + r^2)/(a^2 + (-2 + r)*r))^(3/2)*
            ((-a^3)*m - a*m*r^2 + a^4*omega + r^3*(-I + r*omega) + 
            a^2*(2*I - I*r + 2*r^2*omega)))/
            (sqrt(a^2 + r^2)*sqrt((a^2 + (-2 + r)*r)*
            (a^2 + r^2))*(-2*I*a^3*m*r - 2*I*a*m*r^3 + 
            a^4*(1 + lambda) + r^4*(2 + lambda) + 
            a^2*r*(2 + r*(3 + 2*lambda)))))
        end
    elseif s == -1
        M11 = begin
            (I*r*sqrt(a^2 + r^2)*((-a^3)*m*r - a*m*r^3 + 
            r^5*omega + a^4*(-2*I + r*omega) + 
            2*a^2*r*(2*I - I*r + r^2*omega)))/
            (sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*
            (2*I*a^3*m*r + 2*I*a*m*r^3 + a^4*(-1 + lambda) + 
            r^4*lambda + a^2*r*(2 + r*(-1 + 2*lambda))))
        end
        M12 = begin
            (r^2*sqrt(a^2 + r^2)*sqrt((a^2 + r^2)/
            (a^2 + (-2 + r)*r))*sqrt((a^2 + (-2 + r)*r)*
            (a^2 + r^2)))/(2*I*a^3*m*r + 2*I*a*m*r^3 + 
            a^4*(-1 + lambda) + r^4*lambda + 
            a^2*r*(2 + r*(-1 + 2*lambda)))
        end
        M21 = begin
            (r*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            (a^8*omega*(2*I - r*omega) + 2*a^7*m*(-I + r*omega) + 
            a*m*r^5*(-2*I + I*r + 2*r^2*omega) + 
            2*a^3*m*r^3*(I + 3*r^2*omega) + 
            a^5*m*r*(4*I - 3*I*r + 6*r^2*omega) + 
            a^6*r*(1 - m^2 + lambda - 4*I*omega + 7*I*r*omega - 
            4*r^2*omega^2) - a^4*r^2*(2*(2 + lambda) + 
            r*(-2 + 2*m^2 - 3*lambda + 10*I*omega) - 9*I*r^2*omega + 
            6*r^3*omega^2) + a^2*r^3*(4 - 4*r*(1 + lambda) - 
            r^2*(-1 + m^2 - 3*lambda + 8*I*omega) + 5*I*r^3*omega - 
            4*r^4*omega^2) + r^6*((-2 + r)*lambda - 
            r*omega*(2*I - I*r + r^2*omega))))/
            ((a^2 + r^2)^(3/2)*sqrt((a^2 + (-2 + r)*r)*
            (a^2 + r^2))*(2*I*a^3*m*r + 2*I*a*m*r^3 + 
            a^4*(-1 + lambda) + r^4*lambda + 
            a^2*r*(2 + r*(-1 + 2*lambda))))
        end
        M22 = begin
            -((r^2*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*
            (-r - (I*(a^2 + r^2)*((-a)*m + (a^2 + r^2)*omega))/
            (a^2 + (-2 + r)*r)))/(sqrt(a^2 + r^2)*
            (2*I*a^3*m*r + 2*I*a*m*r^3 + a^4*(-1 + lambda) + 
            r^4*lambda + a^2*r*(2 + r*(-1 + 2*lambda)))))
        end
    elseif s == +2
        M11 = begin
            (r^2*(a^2 + (-2 + r)*r)*(-4*a^5*m*r*(I + r*omega) - 
            2*a^3*m*r^2*(-3*I + 2*I*r + 4*r^2*omega) + 
            a*m*(2*I*r^4 - 4*r^6*omega) + 
            2*a^6*(-5 + 2*I*r*omega + r^2*omega^2) + 
            a^4*r*(32 + r*(-24 + 2*m^2 - lambda - 6*I*omega) + 
            10*I*r^2*omega + 6*r^3*omega^2) + 
            r^4*(-12 + 2*r*(9 + lambda) - 
            r^2*(6 + lambda + 6*I*omega) + 2*I*r^3*omega + 
            2*r^4*omega^2) + 2*a^2*r^2*(-12 + r*(23 + lambda) + 
            r^2*(-10 + m^2 - lambda - 6*I*omega) + 4*I*r^3*omega + 
            3*r^4*omega^2)))/
            (((a^2 + (-2 + r)*r)^2*(a^2 + r^2))^(3/2)*
            ((-r^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            24*a^3*m*r*(I + 2*r*omega) + 4*a*m*r^2*
            (6*I + 2*I*r*(4 + lambda) - 3*r^2*omega) + 
            12*a^4*(-1 + 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*r^2*(1 + lambda)*omega + 3*r^3*omega^2)))
        end
        M12 = begin
            (2*r^3*(a^2 + (-2 + r)*r)*(a^2 + r^2)^2*
            ((-I)*a*m*r + a^2*(-2 + I*r*omega) + 
            r*(3 - r + I*r^2*omega)))/
            (((a^2 + (-2 + r)*r)^2*(a^2 + r^2))^(3/2)*
            ((-r^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            24*a^3*m*r*(I + 2*r*omega) + 4*a*m*r^2*
            (6*I + 2*I*r*(4 + lambda) - 3*r^2*omega) + 
            12*a^4*(-1 + 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*r^2*(1 + lambda)*omega + 3*r^3*omega^2)))
        end
        M21 = begin
            -((I*r*(-2*a^7*m*r*(-2 + 4*I*r*omega + 3*r^2*omega^2) - 
            2*a^5*m*r^2*(-4 + r*(3 + m^2 - lambda + 2*I*omega) + 
            4*I*r^2*omega + 9*r^3*omega^2) + 
            2*a*m*r^5*(-10 + r*(3 - 2*lambda) + 
            r^2*(2 + lambda + 2*I*omega) + 4*I*r^3*omega - 
            3*r^4*omega^2) - 2*a^3*m*r^3*
            (12 + r^2*(3 + m^2 - 2*lambda) + r*(-17 + 2*lambda) - 
            4*I*r^3*omega + 9*r^4*omega^2) + 
            2*a^8*(-6*I - 2*r*omega + 2*I*r^2*omega^2 + 
            r^3*omega^3) + 2*a^6*r*(16*I + 
            I*r*(-9 + 2*m^2 - 2*lambda + 4*I*omega) + 
            r^2*(-8 + 3*m^2 - lambda + I*omega)*omega + 
            6*I*r^3*omega^2 + 4*r^4*omega^3) + 
            r^6*(-4*I*lambda - 2*r^3*(4 + lambda + 5*I*omega)*omega + 
            2*r^5*omega^3 - 12*r*(I + omega) + 
            r^2*(6*I + I*lambda + 18*omega + 4*lambda*omega)) + 
            a^4*r^2*(-16*I + 2*I*r*(10 + m^2 + 6*lambda - 
            12*I*omega) + 2*r^3*(-14 + 6*m^2 - 3*lambda - 
            3*I*omega)*omega + 12*I*r^4*omega^2 + 12*r^5*omega^3 + 
            r^2*(-2*I - 4*I*m^2 - 7*I*lambda + 22*omega + 
            4*lambda*omega)) + 2*a^2*r^4*(-4*I*lambda + 
            I*r*(-4 + 3*m^2 + 6*lambda + 6*I*omega) + 
            3*r^3*(-4 + m^2 - lambda - 3*I*omega)*omega + 
            2*I*r^4*omega^2 + 4*r^5*omega^3 + 
            r^2*(5*I - 4*I*m^2 - I*lambda + 24*omega + 
            4*lambda*omega))))/
            (((a^2 + (-2 + r)*r)^2*(a^2 + r^2))^(3/2)*
            ((-r^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            24*a^3*m*r*(I + 2*r*omega) + 4*a*m*r^2*
            (6*I + 2*I*r*(4 + lambda) - 3*r^2*omega) + 
            12*a^4*(-1 + 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*r^2*(1 + lambda)*omega + 3*r^3*omega^2))))
        end
        M22 = begin
            (r^2*(a^2 + r^2)^2*(-4*a^3*m*r*(I + r*omega) - 
            2*a*m*r^2*(I - 3*I*r + 2*r^2*omega) + 
            2*a^4*(-3 + 2*I*r*omega + r^2*omega^2) + 
            r^3*(2*lambda - r*(2 + lambda + 10*I*omega) + 2*r^3*omega^2) + 
            a^2*r*(8 + 2*m^2*r - r*lambda + 2*I*r*omega + 
            4*I*r^2*omega + 4*r^3*omega^2)))/
            (((a^2 + (-2 + r)*r)^2*(a^2 + r^2))^(3/2)*
            ((-r^4)*(24 + 10*lambda + lambda^2 + 12*I*omega) - 
            24*a^3*m*r*(I + 2*r*omega) + 4*a*m*r^2*
            (6*I + 2*I*r*(4 + lambda) - 3*r^2*omega) + 
            12*a^4*(-1 + 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*I*omega) - 
            2*I*r^2*(1 + lambda)*omega + 3*r^3*omega^2)))
        end
    elseif s == -2
        M11 = begin
            (r^2*(-4*a^5*m*r*(-I + r*omega) - 
            2*a*m*r^4*(I + 2*r^2*omega) - 2*a^3*m*r^2*
            (3*I - 2*I*r + 4*r^2*omega) + 
            2*a^6*(-5 - 2*I*r*omega + r^2*omega^2) + 
            a^4*r*(32 + r*(-20 + 2*m^2 - lambda + 6*I*omega) - 
            10*I*r^2*omega + 6*r^3*omega^2) + 
            r^4*(-12 + 2*r*(5 + lambda) - 
            r^2*(2 + lambda - 6*I*omega) - 2*I*r^3*omega + 
            2*r^4*omega^2) + 2*a^2*r^2*(-12 + r*(19 + lambda) + 
            r^2*(-6 + m^2 - lambda + 6*I*omega) - 4*I*r^3*omega + 
            3*r^4*omega^2)))/((a^2 + (-2 + r)*r)^3*
            ((a^2 + r^2)/(a^2 + (-2 + r)*r)^2)^(3/2)*
            ((-r^4)*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*r*(I - 2*r*omega) - 4*a*m*r^2*
            (6*I + 2*I*r*lambda + 3*r^2*omega) + 
            12*a^4*(-1 - 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*r^2*(-3 + lambda)*omega + 3*r^3*omega^2)))
        end
        M12 = begin
            (2*r^3*(a^2 + (-2 + r)*r)*
            sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r)^2)*
            (I*a*m*r + a^2*(-2 - I*r*omega) + 
            r*(3 - r - I*r^2*omega)))/
            ((-r^4)*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*r*(I - 2*r*omega) - 
            4*a*m*r^2*(6*I + 2*I*r*lambda + 3*r^2*omega) + 
            12*a^4*(-1 - 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))
        end
        M21 = begin
            (I*r*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r)^2)*
            (2*a^7*m*r*(2 + 4*I*r*omega - 3*r^2*omega^2) - 
            2*a^5*m*r^2*(4 + r*(-1 + m^2 - lambda + 6*I*omega) - 
            12*I*r^2*omega + 9*r^3*omega^2) - 
            2*a^3*m*r^4*(-5 + 2*lambda + r*(3 + m^2 - 2*lambda + 
            16*I*omega) - 12*I*r^2*omega + 9*r^3*omega^2) - 
            2*a*m*r^5*(6 + r*(-7 + 2*lambda) + 
            r^2*(2 - lambda + 10*I*omega) - 4*I*r^3*omega + 
            3*r^4*omega^2) + 2*a^8*(6*I - 2*r*omega - 
            2*I*r^2*omega^2 + r^3*omega^3) + 
            2*a^6*r*(-36*I + r^2*(-12 + 3*m^2 - lambda + 3*I*omega)*
            omega - 10*I*r^3*omega^2 + 4*r^4*omega^3 + 
            r*(21*I - 2*I*m^2 + 2*I*lambda + 4*omega)) + 
            r^5*(-48*I + 12*I*r*(6 + lambda) - 
            12*I*r^2*(3 + lambda - 3*I*omega) - 
            2*r^4*(4 + lambda - 9*I*omega)*omega - 8*I*r^5*omega^2 + 
            2*r^6*omega^3 + r^3*(6*I + 3*I*lambda + 34*omega + 
            4*lambda*omega)) + 2*a^2*r^3*(-48*I + 
            4*I*r*(27 + 2*lambda) + I*r^2*(-72 + m^2 - 14*lambda + 
            30*I*omega) + r^4*(-16 + 3*m^2 - 3*lambda + 
            21*I*omega)*omega - 14*I*r^5*omega^2 + 4*r^6*omega^3 + 
            r^3*(15*I + 5*I*lambda + 48*omega + 4*lambda*omega)) + 
            a^4*r^2*(144*I + 2*I*r*(-90 + 3*m^2 - 8*lambda) + 
            2*r^3*(-22 + 6*m^2 - 3*lambda + 15*I*omega)*omega - 
            36*I*r^4*omega^2 + 12*r^5*omega^3 + 
            r^2*(54*I - 4*I*m^2 + 11*I*lambda + 70*omega + 
            4*lambda*omega))))/((a^2 + r^2)^2*
            ((-r^4)*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*r*(I - 2*r*omega) - 4*a*m*r^2*
            (6*I + 2*I*r*lambda + 3*r^2*omega) + 
            12*a^4*(-1 - 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*r^2*(-3 + lambda)*omega + 3*r^3*omega^2)))
        end
        M22 = begin
            (r^2*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r)^2)*
            (-4*a^3*m*r*(-I + r*omega) - 2*a*m*r^2*
            (3*I - I*r + 2*r^2*omega) + 
            2*a^4*(-3 - 2*I*r*omega + r^2*omega^2) + 
            a^2*r*(24 + r*(-12 + 2*m^2 - lambda + 6*I*omega) - 
            12*I*r^2*omega + 4*r^3*omega^2) + 
            r^2*(-24 + 2*r*(12 + lambda) - 
            r^2*(6 + lambda - 18*I*omega) - 8*I*r^3*omega + 
            2*r^4*omega^2)))/
            ((-r^4)*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*r*(I - 2*r*omega) - 
            4*a*m*r^2*(6*I + 2*I*r*lambda + 3*r^2*omega) + 
            12*a^4*(-1 - 2*I*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*I*omega) + 
            2*I*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))
        end
    else
        # No explicit form available; revert to general-s expression
        return _Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix_general_s(s, m, a, omega, lambda, r)
    end
    
    return [M11 M12; M21 M22]
end

function _Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix_general_s(s, m, a, omega, lambda, r)
    drstar_dr(r) = (r^2 + a^2)/Delta(a, r)
    chi_conversion_factor(r) = (1.0/sqrt((r^2 + a^2) * (Delta(a, r)^s)))
    dchi_conversion_factor_dr(r) = begin
        ((-a^2)*(r - s + r*s) + r^2*(2 + s - r*(1 + s)))/
        ((a^2 + (-2 + r)*r)*sqrt((a^2 + (-2 + r)*r)^s)*(a^2 + r^2)^(3/2))
    end
    _eta(r) = eta(s, m, a, omega, lambda, r)
    _alpha(r) = alpha(s, m, a, omega, lambda, r)
    _alpha_prime(r) = alpha_prime(s, m, a, omega, lambda, r)
    _beta(r) = beta(s, m, a, omega, lambda, r)
    _beta_prime(r) = beta_prime(s, m, a, omega, lambda, r)
    _Delta(r) = Delta(a, r)
    _VT(r) = VT(s, m, a, omega, lambda, r)
    _1_over_eta(r) = 1.0/_eta(r)

    conversion_matrix_from_X_dXdrs_to_X_dXdr(r) = [1 0 ; 0 drstar_dr(r)]
    conversion_matrix_from_X_dXdr_to_chi_dchidr(r) = [chi_conversion_factor(r) 0 ; dchi_conversion_factor_dr(r) chi_conversion_factor(r)]
    conversion_matrix_from_chi_dchidr_to_R_dRdr(r) = _1_over_eta(r) * [_alpha(r)+_beta_prime(r)*(_Delta(r)^(s+1)) -_beta(r)*(_Delta(r)^(s+1)) ; -(_alpha_prime(r)+_beta(r)*_VT(r)*(_Delta(r)^(s))) _alpha(r)]
    # **Left multiplication**
    overall_conversion_matrix = conversion_matrix_from_chi_dchidr_to_R_dRdr(r) * conversion_matrix_from_X_dXdr_to_chi_dchidr(r) * conversion_matrix_from_X_dXdrs_to_X_dXdr(r)
    return overall_conversion_matrix
end

function Teukolsky_radial_function_from_Sasaki_Nakamura_function(Xsoln)
    # Unpack the parameters
    s = Xsoln.prob.p.s
    m = Xsoln.prob.p.m
    a = Xsoln.prob.p.a
    omega = Xsoln.prob.p.omega
    lambda = Xsoln.prob.p.lambda

    #=
    First convert [X(rs), dX/drs(rs)] to [X(r), dX/dr(r)], this is done by
    [X(r), dX/dr(r)]^T = [1, 0; 0, drstar/dr ] * [X(rs), dX/drs(rs)]^T

    Then we convert [X(r), dX/dr(r)] to [chi(r), dchi/dr(r)] by
    [chi(r), dchi/dr(r)]^T = [chi_conversion_factor(r), 0 ; dchi_conversion_factor_dr, chi_conversion_factor] * [X(r), dX/dr(r)]^T

    After that we convert [chi(r), dchi/dr(r)] to [R(r), dR/dr(r)] by
    [R(r), dR/dr(r)]^T = 1/eta(r) * [ alpha + beta_prime*Delta^(s+1), -beta*Delta^(s+1) ; -(alpha_prime + beta*VT*Delta^s), alpha ] *  [chi(r), dchi/dr(r)]^T

    Therefore the overall conversion matrix is 'just' (one matrix for each r)
    the multiplication of each conversion matrix
    =#

    overall_conversion_matrix(r) = Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(s, m, a, omega, lambda, r)
    # Reconstruct X, Xprime from PhiRePhiIm solution
    X, Xprime = XXprime_from_PhiRePhiImsoln(Xsoln)
    Rsoln = (r -> overall_conversion_matrix(r) * [X(rstar_from_r(a, r)); Xprime(rstar_from_r(a, r))])
    return Rsoln
end

function d2Rdr2_from_Rsoln(s::Int, m::Int, a, omega, lambda, Rsoln, r)
    #=
    Using the radial Teukolsky equation we can solve for d2Rdr2 from R and dRdr using

        d2Rdr2 = VT/\Delta R - (2(s+1)(r-M))/\Delta dRdr
    =#
    # NOTE DO NOT USE THE DOT PRODUCT IN LINEAR ALGEBRA
    R, dRdr = Rsoln(r)
    return (VT(s, m, a, omega, lambda, r)/Delta(a, r))*R - ((2*(s+1)*(r-1))/Delta(a,r))*dRdr
end

function scaled_Wronskian(Rin_soln, Rup_soln, r, s, a)
    # The scaled Wronskian is given by W = Delta^{s+1} * det([Rin Rup; Rin' Rup'])
    Rin, Rin_prime = Rin_soln(r)
    Rup, Rup_prime = Rup_soln(r)
    return Delta(a, r)^(s+1) * (Rin*Rup_prime - Rup*Rin_prime)
end

function _extract_asymptotic_amplitude_from_Xsoln(osc_variable, sign, Xsoln, rs_extraction)
    X, Xprime = XXprime_from_PhiRePhiImsoln(Xsoln)
    # This is an internal template function
    return ((exp((-1*sign)*1im*osc_variable*rs_extraction)/(2*1im*osc_variable))*(1im*osc_variable*X(rs_extraction) + sign*Xprime(rs_extraction)))
end

function CrefCinc_SN_from_Xup(Xupsoln, rsin)
    # Extract oscillation variable from Xupsoln
    m = Xupsoln.prob.p.m
    a = Xupsoln.prob.p.a
    omega = Xupsoln.prob.p.omega
    p = omega - m*omega_horizon(a)

    Cinc_SN = _extract_asymptotic_amplitude_from_Xsoln(p, 1, Xupsoln, rsin)
    Cref_SN = _extract_asymptotic_amplitude_from_Xsoln(p, -1, Xupsoln, rsin)
    return Cref_SN, Cinc_SN
end

function _correction_factor(func, omega_times_r, order)
    # This is an internal helper function
    _cf = 0.0
    for j in collect(0:1:order)
        _cf += func(j)/((omega_times_r)^j)
    end
    return _cf
end

function BrefBinc_SN_from_Xin(Xinsoln, rsout; order=3)
    # Unpack the parameters
    s = Xinsoln.prob.p.s
    m = Xinsoln.prob.p.m
    a = Xinsoln.prob.p.a
    omega = Xinsoln.prob.p.omega
    lambda = Xinsoln.prob.p.lambda

    Bref_SN = _extract_asymptotic_amplitude_from_Xsoln(omega, 1, Xinsoln, rsout)
    Binc_SN = _extract_asymptotic_amplitude_from_Xsoln(omega, -1, Xinsoln, rsout)

    # Compute the high-order correction factors
    _r = r_from_rstar(a, rsout)
    _outgoing_coefficient_func(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
    _ingoing_coefficient_func(ord) = ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
    cf_outgoing = _correction_factor(_outgoing_coefficient_func, omega*_r, order)
    cf_ingoing = _correction_factor(_ingoing_coefficient_func, omega*_r, order)

    return Bref_SN/cf_outgoing, Binc_SN/cf_ingoing
end

end