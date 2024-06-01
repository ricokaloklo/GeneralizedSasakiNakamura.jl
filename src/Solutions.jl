module Solutions

using DifferentialEquations
using ForwardDiff
using HypergeometricFunctions
using ..Kerr
using ..Transformation
using ..Coordinates
using ..Potentials
using ..AsymptoticExpansionCoefficients
using ..InitialConditions
using ..ConversionFactors

const I = 1im # Mathematica being Mathematica
_DEFAULTDATATYPE = ComplexF64 # Double precision by default
_DEFAULTSOLVER = Vern9()
_DEFAULTTOLERANCE = 1e-12

# First-order non-linear ODE form of the GSN equation
function GSN_Riccati_eqn!(du, u, p, rs)
    r = r_from_rstar(p.a, rs)
    _sF = sF(p.s, p.m, p.a, p.omega, p.lambda, r)
    _sU = sU(p.s, p.m, p.a, p.omega, p.lambda, r)
    
    #=
    We write X = exp(I*Phi)
    Substitute X in this form into the GSN equation will give
    a Riccati equation, a first-order non-linear equation

    u[1] = Phi
    u[2] = dPhidrs
    =#
    du[1] = u[2]
    du[2] = -1im*_sU + _sF*u[2] - 1im*u[2]*u[2]
end

# Prufer phase function
function GSN_Prufer_eqn!(du, u, p, rs)
    r = r_from_rstar(p.a, rs)
    _sF = sF(p.s, p.m, p.a, p.omega, p.lambda, r)
    _sU = sU(p.s, p.m, p.a, p.omega, p.lambda, r)

    #=
    u[1] = -i log X = Phi = k cot(k*rs + \tilde{P})
    u[2] = \tilde{P} is the Prufer phase function
    =#
    du[1] = -1im * p.omega*cot(p.omega*rs + u[2])
    du[2] = -(p.omega + cot(p.omega*rs + u[2])*_sF + _sU/p.omega)*(sin(p.omega*rs + u[2]))^2
end

# Second-order linear ODE form of the GSN equation
function GSN_linear_eqn!(du, u, p, rs)
    r = r_from_rstar(p.a, rs)
    _sF = sF(p.s, p.m, p.a, p.omega, p.lambda, r)
    _sU = sU(p.s, p.m, p.a, p.omega, p.lambda, r)

    #=
    Using the convention for DifferentialEquations
    u[1] = X(rs)
    u[2] = dX/drs = X'
    therefore
    X'' - sF X' - sU X = 0 => u[2]' - sF u[2] - sU u[1] = 0 => u[2]' = sF u[2] + sU u[1]
    =#
    du[1] = u[2]
    du[2] = _sF*u[2] + _sU*u[1]
end

function PhiPhiprime_from_XXprime(X, Xprime)
    Phi = -1im*log(X)
    Phiprime = -1im*Xprime/X

    return Phi, Phiprime
end

function XXprime_from_PhiPhiprime(Phi, Phiprime)
    X = exp(1im*Phi)
    Xprime = 1im*X*Phiprime

    return X, Xprime
end

function Xsoln_from_Phisoln(Phisoln)
    return rs -> XXprime_from_PhiPhiprime(Phisoln(rs)[1], Phisoln(rs)[2])
end

function solve_Phiup(s::Int, m::Int, a, omega, lambda, rsin, rsout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsout, "rsout ($rsout) must be larger than rsin ($rsin)"))
    end
    # Initial conditions at rs = rsout, the outer boundary
    Xup_rsout, Xupprime_rsout = Xup_initialconditions(s, m, a, omega, lambda, rsout; order=initialconditions_order)
    # Convert initial conditions for Xup for Phi
    Phi, Phiprime = PhiPhiprime_from_XXprime(Xup_rsout, Xupprime_rsout)
    u0 = [dtype(Phi); dtype(Phiprime)]
    rsspan = (rsout, rsin) # Integrate from rsout to rsin *inward*
    p = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    odeprob = ODEProblem(GSN_Riccati_eqn!, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
    return odesoln
end

function solve_Phiin(s::Int, m::Int, a, omega, lambda, rsin, rsout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsin, "rsin ($rsin) must be smaller than rsout ($rsout)"))
    end
    # Initial conditions at rs = rsin, the inner boundary; this should be very close to EH
    Xin_rsin, Xinprime_rsin = Xin_initialconditions(s, m, a, omega, lambda, rsin; order=initialconditions_order)
    # Convert initial conditions for Xin for PhiRe PhiIm
    Phi, Phiprime = PhiPhiprime_from_XXprime(Xin_rsin, Xinprime_rsin)
    u0 = [dtype(Phi); dtype(Phiprime)]
    rsspan = (rsin, rsout) # Integrate from rsin to rsout *outward*
    p = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    odeprob = ODEProblem(GSN_Riccati_eqn!, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
    return odesoln
end

function solve_Xup(s::Int, m::Int, a, omega, lambda, rsin, rsout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsout, "rsout ($rsout) must be larger than rsin ($rsin)"))
    end
    # Initial conditions at rs = rsout, the outer boundary
    Xup_rsout, Xupprime_rsout = Xup_initialconditions(s, m, a, omega, lambda, rsout; order=initialconditions_order)
    u0 = [dtype(Xup_rsout); dtype(Xupprime_rsout)]
    rsspan = (rsout, rsin) # Integrate from rsout to rsin *inward*
    p = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    odeprob = ODEProblem(GSN_linear_eqn!, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
    return odesoln
end

function solve_Xin(s::Int, m::Int, a, omega, lambda, rsin, rsout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsin, "rsin ($rsin) must be smaller than rsout ($rsout)"))
    end
    # Initial conditions at rs = rsin, the inner boundary; this should be very close to EH
    Xin_rsin, Xinprime_rsin = Xin_initialconditions(s, m, a, omega, lambda, rsin; order=initialconditions_order)
    u0 = [dtype(Xin_rsin); dtype(Xinprime_rsin)]
    rsspan = (rsin, rsout) # Integrate from rsin to rsout *outward*
    p = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    odeprob = ODEProblem(GSN_linear_eqn!, u0, rsspan, p)
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
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end 
    return [M11 M12; M21 M22]
end

function Teukolsky_radial_function_from_Sasaki_Nakamura_function(s::Int, m::Int, a, omega, lambda, Xsoln)
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
    X(rs) = Xsoln(rs)[1]
    Xprime(rs) = Xsoln(rs)[2]
    Rsoln = (r -> overall_conversion_matrix(r) * [X(rstar_from_r(a, r)); Xprime(rstar_from_r(a, r))])
    return Rsoln
end

function Sasaki_Nakamura_function_from_Teukolsky_radial_function(s::Int, m::Int, a, omega, lambda, Rsoln)
    #=
    Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix() returns 
    the conversion matrix for going from (X, dX/drs) to (R, dR/dr). Here we need
    the inverse of that to go from (R, dR/dr) to (X, dX/drs).
    =#
    conversion_matrix(r) = inv(Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(s, m, a, omega, lambda, r))

    Xsoln(rs) = begin
        r = r_from_rstar(a, rs)
        conversion_matrix(r) * Rsoln(r)
    end

    return Xsoln
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

function scaled_Wronskian_Teukolsky(Rin_soln, Rup_soln, r, s, a)
    # The scaled Wronskian is given by W = Delta^{s+1} * det([Rin Rup; Rin' Rup'])
    Rin, Rin_prime = Rin_soln(r)
    Rup, Rup_prime = Rup_soln(r)
    return Delta(a, r)^(s+1) * (Rin*Rup_prime - Rup*Rin_prime)
end

function scaled_Wronskian_GSN(Xin_soln, Xup_soln, rs, s, m, a, omega, lambda)
    r = r_from_rstar(a, rs)
    Xin, dXindrs = Xin_soln(rs)
    Xup, dXupdrs = Xup_soln(rs)
    _eta = eta(s, m, a, omega, lambda, r)

    return (Xin*dXupdrs - dXindrs*Xup)/_eta
end

function scaled_Wronskian_from_Phisolns(Phiin_soln, Phiup_soln, rs, s, m, a, omega, lambda)
    r = r_from_rstar(a, rs)
    Xin = exp(1im*Phiin_soln(rs)[1])
    dXindrs = 1im*exp(1im*Phiinsoln(rs)[1])*Phiinsoln(rs)[2]
    Xup = exp(1im*Phiup_soln(rs)[1])
    dXupdrs = 1im*exp(1im*Phiupsoln(rs)[1])*Phiupsoln(rs)[2]
    _eta = eta(s, m, a, omega, lambda, r)

    return (Xin*dXupdrs - dXindrs*Xup)/_eta
end

function residual_from_Xsoln(s::Int, m::Int, a, omega, lambda, Xsoln)
    # Compute the second derivative using autodiff instead of finite diff
    X(rs) = Xsoln(rs)[1]
    first_deriv(rs) = Xsoln(rs)[2]
    second_deriv(rs) = ForwardDiff.derivative(first_deriv, rs)

    _sF(rs) = sF(s, m, a, omega, lambda, r_from_rstar(a, rs))
    _sU(rs) = sU(s, m, a, omega, lambda, r_from_rstar(a, rs))

    return rs -> second_deriv(rs) - _sF(rs)*first_deriv(rs) - _sU(rs)*X(rs)
end

function residual_RiccatiEqn_from_Phisoln(s::Int, m::Int, a, omega, lambda, Phisoln)
    Phi(rs) = Phisoln(rs)[1] # Does not matter actually!
    first_deriv(rs) = Phisoln(rs)[2]
    second_deriv(rs) = ForwardDiff.derivative(first_deriv, rs)

    _sF(rs) = sF(s, m, a, omega, lambda, r_from_rstar(a, rs))
    _sU(rs) = sU(s, m, a, omega, lambda, r_from_rstar(a, rs))

    return rs -> second_deriv(rs) + 1im*_sU(rs) - _sF(rs)*first_deriv(rs) + 1im*(first_deriv(rs))^2
end

function residual_GSNEqn_from_Phisoln(s::Int, m::Int, a, omega, lambda, Phisoln)
    # Compute the second derivative using autodiff instead of finite diff
    X = (rs -> exp(1im*Phisoln(rs)[1]))
    first_deriv = (rs -> 1im*exp(1im*Phisoln(rs)[1])*Phisoln(rs)[2])
    second_deriv(rs) = ForwardDiff.derivative(first_deriv, rs)

    _sF(rs) = sF(s, m, a, omega, lambda, r_from_rstar(a, rs))
    _sU(rs) = sU(s, m, a, omega, lambda, r_from_rstar(a, rs))

    return rs -> second_deriv(rs) - _sF(rs)*first_deriv(rs) - _sU(rs)*X(rs)
end

function CrefCinc_SN_from_Xup(s::Int, m::Int, a, omega, lambda, Xupsoln, rsin; order=0)
    p = omega - m*omega_horizon(a)

    rin = r_from_rstar(a, rsin)
    # Computing A1, A2, A3, A4
    gin(r) = gansatz(
        ord -> ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r;
        order=order
    )
    A1 = gin(rin) * exp(-1im*p*rsin)
    gout(r) = gansatz(
        ord -> outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
        a,
        r;
        order=order
    )
    A2 = gout(rin) * exp(1im*p*rsin)

    _DeltaOverr2pa2 = Delta(a, rin)/(rin^2 + a^2)
    A3 = (_DeltaOverr2pa2 * ForwardDiff.derivative(gin, rin) - 1im*p*gin(rin)) * exp(-1im*p*rsin)
    A4 = (_DeltaOverr2pa2 * ForwardDiff.derivative(gout, rin) + 1im*p*gout(rin)) * exp(1im*p*rsin)

    X(rs) = Xupsoln(rs)[1]
    Xprime(rs) = Xupsoln(rs)[2]
    C1 = X(rsin)
    C2 = Xprime(rsin)

    return -(A4*C1 - A2*C2)/(A2*A3 - A1*A4), -(-A3*C1 + A1*C2)/(A2*A3 - A1*A4)
end

function BrefBinc_SN_from_Xin(s::Int, m::Int, a, omega, lambda, Xinsoln, rsout; order=3)
    rout = r_from_rstar(a, rsout)
    # Computing A1, A2, A3, A4
    fin(r) = fansatz(
        ord -> ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r;
        order=order
    )
    A1 = fin(rout) * exp(-1im*omega*rsout)
    fout(r) = fansatz(
        ord -> outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
        omega,
        r;
        order=order
    )
    A2 = fout(rout) * exp(1im*omega*rsout)

    _DeltaOverr2pa2 = Delta(a, rout)/(rout^2 + a^2)
    A3 = (_DeltaOverr2pa2 * ForwardDiff.derivative(fin, rout) - 1im*omega*fin(rout)) * exp(-1im*omega*rsout)
    A4 = (_DeltaOverr2pa2 * ForwardDiff.derivative(fout, rout) + 1im*omega*fout(rout)) * exp(1im*omega*rsout)

    X(rs) = Xinsoln(rs)[1]
    Xprime(rs) = Xinsoln(rs)[2]
    C1 = X(rsout)
    C2 = Xprime(rsout)

    return -(-A3*C1 + A1*C2)/(A2*A3 - A1*A4), -(A4*C1 - A2*C2)/(A2*A3 - A1*A4)
end

function semianalytical_Xin(s, m, a, omega, lambda, Xinsoln, rsin, rsout, horizon_expansionorder, infinity_expansionorder, rs)
    _r = r_from_rstar(a, rs) # Evaluate at this r

    if rs < rsin
        # Extend the numerical solution to the analytical ansatz from rsin to horizon

        # Construct the analytical ansatz
        p = omega - m*omega_horizon(a)
        gin(r) = gansatz(
            ord -> ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
            a,
            r;
            order=horizon_expansionorder
        )

        _Xin = gin(_r)*exp(-1im*p*rs) # Evaluate at *this* rs

        _DeltaOverr2pa2 = Delta(a, _r)/(_r^2 + a^2)
        _dXindrs = (_DeltaOverr2pa2 * ForwardDiff.derivative(gin, _r) - 1im*p*gin(_r))*exp(-1im*p*rs)

        return (_Xin, _dXindrs)
    elseif rs > rsout
        # Extend the numerical solution to the analytical ansatz from rsout to infinity

        # Obtain the coefficients by imposing continuity in X and dX/drs
        Bref_SN, Binc_SN = BrefBinc_SN_from_Xin(s, m, a, omega, lambda, Xinsoln, rsout; order=infinity_expansionorder)

        # Construct the analytical ansatz
        fin(r) = fansatz(
            ord -> ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
            omega,
            r;
            order=infinity_expansionorder
        )
        fout(r) = fansatz(
            ord -> outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
            omega,
            r;
            order=infinity_expansionorder
        )

        _Xin = Bref_SN*fout(_r)*exp(1im*omega*rs) + Binc_SN*fin(_r)*exp(-1im*omega*rs)

        _DeltaOverr2pa2 = Delta(a, _r)/(_r^2 + a^2)
        _dXindrs = Bref_SN*(_DeltaOverr2pa2 * ForwardDiff.derivative(fout, _r) + 1im*omega*fout(_r))*exp(1im*omega*rs) + Binc_SN*(_DeltaOverr2pa2 * ForwardDiff.derivative(fin, _r) - 1im*omega*fin(_r))*exp(-1im*omega*rs)
        
        return (_Xin, _dXindrs)
    else
        # Requested rs is within the numerical solution
        return Xinsoln(rs)
    end 
end

function semianalytical_Xup(s, m, a, omega, lambda, Xupsoln, rsin, rsout, horizon_expansionorder, infinity_expansionorder, rs)
    _r = r_from_rstar(a, rs) # Evaluate at this r

    if rs < rsin
        # Extend the numerical solution to the analytical ansatz from rsin to horizon

        # Obtain the coefficients by imposing continuity in X and dX/drs
        Cref_SN, Cinc_SN = CrefCinc_SN_from_Xup(s, m, a, omega, lambda, Xupsoln, rsin; order=horizon_expansionorder)
        # with coefficient Cref_SN for the left-going and Cinc_SN for the right-going mode

        # Construct the analytical ansatz
        p = omega - m*omega_horizon(a)
        gin(r) = gansatz(
            ord -> ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
            a,
            r;
            order=horizon_expansionorder
        )
        gout(r) = gansatz(
            ord -> outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord),
            a,
            r;
            order=horizon_expansionorder
        )

        _Xup = Cinc_SN*gout(_r)*exp(1im*p*rs) + Cref_SN*gin(_r)*exp(-1im*p*rs) # Evaluate at *this* rs

        _DeltaOverr2pa2 = Delta(a, _r)/(_r^2 + a^2)
        _dXupdrs = Cinc_SN*(_DeltaOverr2pa2 * ForwardDiff.derivative(gout, _r) + 1im*p*gout(_r))*exp(1im*p*rs) + Cref_SN*(_DeltaOverr2pa2 * ForwardDiff.derivative(gin, _r) - 1im*p*gin(_r))*exp(-1im*p*rs)

        return (_Xup, _dXupdrs)
    elseif rs > rsout
        # Extend the numerical solution to the analytical ansatz from rsout to infinity

        fout(r) = fansatz(
            ord -> outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord),
            omega,
            r;
            order=infinity_expansionorder
        )

        _Xup = fout(_r)*exp(1im*omega*rs)

        _DeltaOverr2pa2 = Delta(a, _r)/(_r^2 + a^2)
        _dXupdrs = (_DeltaOverr2pa2 * ForwardDiff.derivative(fout, _r) + 1im*omega*fout(_r)) * exp(1im*omega*rs)
        
        return (_Xup, _dXupdrs)
    else
        # Requested rs is within the numerical solution
        return Xupsoln(rs)
    end 
end

# Some helper functions for static modes
function x_from_r(a, r)
    #=
    x \equiv (rp - r)/(rp - rm) is a transformed coordinate that
    is used only internally in the code/static mode computation
    =#
    gamma = sqrt(1 - a^2)
    return ((1+gamma)-r)/(2*gamma)
end

function r_from_x(a, x)
    gamma = sqrt(1 - a^2)
    return 1 + gamma - 2*gamma*x
end

function isnegativeinteger(x)
    isinteger(x) && x < 0
end

# Some helper functions for manipulating hypergeometric functions
function pochhammer(x, n)
    # Wrapper for the pochhammer function implemented in HypergeometricFunctions.jl
    HypergeometricFunctions.pochhammer(x, n)
end

#=
Since there is no implementation of a "regularized" Gauss hypergeometric function in
HypergeometricFunctions.jl yet, we use the fact that when c approaches a negative
integer, 2F1(a,b;c,z)/Gamma(c) approaches the value

lim c->-m 2F1(a, b; c, z)/Gamma(c) = (a)_m (b)_m / m! * z^m * 2F1(a+m, b+m; m+1; z)

which is implemented below
=#
function _2F1_over_Gamma_c(a, b, c, z)
    # NOTE This is true only when c is a negative integer
    z *= (1 + 0im) # Make sure that z is a complex number
    m = abs(c) + 1
    ( pochhammer(a, m) * pochhammer(b, m) / factorial(m) ) * z^m * pFq((a+m, b+m), (m+1, ), z)
end

function _d2F1_over_Gamma_c_dz(a, b, c, z)
    # NOTE This is true only when c is a negative integer
    z *= (1 + 0im) # Make sure that z is a complex number
    m = abs(c) + 1
    # Here we explicitly compute the derivative because ForwardDiff.jl does not work very well with HypergeometricFunctions.jl
    ( pochhammer(a, m) * pochhammer(b, m) / factorial(m) ) * (m*z^(m-1) * pFq((a+m, b+m), (m+1, ), z) + z^m * pFq((a+m+1, b+m+1), (m+2, ), z))
end

function solve_static_Rin(s::Int, l::Int, m::Int, a)
    #=
    This is an internal function that returns the static (omega = 0) solution that satisfies
    the "purely left-going" condition at the horizon, namely the R^in solution.

    Use instead GeneralizedSasakiNakamura.Teukolsky_radial() and set omega = 0 and boundary_condition = IN
    =#
    gamma = sqrt(1 - a^2)
    kappa = I*a*m/gamma

    normalization_const = begin
        (-1)^(-kappa/2) * exp(I*a*m/2) * gamma^(I*a*m/(1+gamma)) * (-4*gamma^2)^(-s)
    end

    Rin(r) = begin
        x = x_from_r(a, r)
        if iszero(kappa) && isnegativeinteger(1-s)
            # Regularization
            normalization_const * x^(-s + kappa/2) * (1-x)^(kappa/2) * _2F1_over_Gamma_c(-l+kappa, 1+l+kappa, 1-s, x)
        else
            normalization_const * x^(-s + kappa/2) * (1-x)^(kappa/2) * pFq((-l+kappa, 1+l+kappa), (1-s+kappa,), x)
        end
    end

    dRindr(r) = begin
        x = x_from_r(a, r)
        if iszero(kappa) && isnegativeinteger(1-s)
            # Regularization
            normalization_const * (-1/(2*gamma)) * ( (1/2) * (1-x)^(-1+kappa/2) * x^(-1-s+kappa/2) * (2*s*(x-1) + kappa - 2*kappa*x) * _2F1_over_Gamma_c(-l+kappa, 1+l+kappa, 1-s, x) + x^(-s + kappa/2) * (1-x)^(kappa/2) * _d2F1_over_Gamma_c_dz(-l+kappa, 1+l+kappa, 1-s, x) )
        else
            normalization_const * (-1/(2*gamma)) * ( (1/2) * (1-x)^(-1+kappa/2) * x^(-1-s+kappa/2) * (2*s*(x-1) + kappa - 2*kappa*x) * pFq((-l+kappa, 1+l+kappa), (1-s+kappa,), x) + x^(-s + kappa/2) * (1-x)^(kappa/2) * pochhammer(-l+kappa, 1) * pochhammer(1+l+kappa, 1) / pochhammer(1-s+kappa, 1) * pFq((-l+kappa+1, 1+l+kappa+1), (1-s+kappa+1,), x) )
        end
    end

    return (r -> [Rin(r);  dRindr(r)])
end

function solve_static_Rup(s::Int, l::Int, m::Int, a)
    #=
    This is an internal function that returns the static (omega = 0) solution that satisfies
    the "purely right-going" condition at spatial infinity, namely the R^up solution.

    Use instead GeneralizedSasakiNakamura.Teukolsky_radial() and set omega = 0 and boundary_condition = UP
    =#
    gamma = sqrt(1 - a^2)
    kappa = I*a*m/gamma

    normalization_const = begin
        (-2*gamma)^(-s-l-1)
    end

    Rup(r) = begin
        x = x_from_r(a, r)
        normalization_const * x^(-s-l-1-kappa/2) * (x-1)^(kappa/2) * pFq((1+l+kappa, 1+s+l), (2+2*l,),1/x)
    end
    
    dRupdr(r) = begin
        x = x_from_r(a, r)
        normalization_const * (-1/(2*gamma)) * ( ( (1/2)*(x-1)^(-1+kappa/2) * x^(-2-l-s-kappa/2) * (-2*(1+l+s)*(x-1) + kappa) ) * pFq((1+l+kappa, 1+s+l), (2+2*l,),1/x) + x^(-s-l-1-kappa/2) * (x-1)^(kappa/2) * pochhammer(1+l+kappa, 1) * pochhammer(1+s+l, 1) / pochhammer(2+2*l, 1) * pFq((2+l+kappa, 2+s+l), (3+2*l,),1/x) * (-1/x^2) )
    end

    return (r-> [Rup(r); dRupdr(r)])
end



end