using DifferentialEquations

include("coordinates.jl")
include("potentials.jl")
include("initialconditions.jl")

function generalized_sasaki_nakamura_equation!(du, u, p, rs)
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

function solve_Xinf(s::Int, m::Int, a, omega, lambda, rsin, rsout; reltol=1e-8, abstol=1e-8)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsout, "rsout ($rsout) must be larger than rsin ($rsin)"))
    end
    # Initial conditions at rs = rsout, the outer boundary
    Xinf_rsout, Xinfprime_rsout = Xinf_initialconditions(s, m, a, omega, lambda, rsout)
    u0 = [Xinf_rsout; Xinfprime_rsout]
    rsspan = (rsout, rsin) # Integrate from rsout to rsin *inward*
    p = (s=s, m=m, a=a, omega=a, lambda=lambda)
    odeprob = ODEProblem(generalized_sasaki_nakamura_equation!, u0, rsspan, p)
    #=
    For now, we will stick to the 'simple' Runge-Kutta 4th order algorithm
    In numerical recipes it was described that Bulirsch-Stoer algorithm is
    'superior' to RK4 but it was very much disputed.
    =#
    odealgo = RK4()
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
end

function solve_XH(s::Int, m::Int, a, omega, lambda, rsin, rsout; reltol=1e-8, abstol=1e-8)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsin, "rsin ($rsin) must be smaller than rsout ($rsout)"))
    end
    # Initial conditions at rs = rsin, the inner boundary; this should be very close to EH
    XH_rsin, XHprime_rsin = XH_initialconditions(s, m, a, omega, lambda, rsin)
    u0 = [XH_rsin; XHprime_rsin]
    rsspan = (rsin, rsout) # Integrate from rsin to rsout *outward*
    p = (s=s, m=m, a=a, omega=a, lambda=lambda)
    odeprob = ODEProblem(generalized_sasaki_nakamura_equation!, u0, rsspan, p)
    odealgo = RK4()
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
end
