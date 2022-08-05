using DifferentialEquations

include("kerr.jl")
include("transformation.jl")
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

function Teukolsky_radial_function_from_Sasaki_Nakamura_function(s::Int, m::Int, a, omega, lambda, Xsoln; rspan=nothing)
    local rsgrid, rgrid, X_Xprime_vector
    if isnothing(rspan)
        # Use the grid that was used to perform the integration
        rsgrid = Xsoln.t
        rgrid = (x -> r_from_rstar(a, x)).(rsgrid)
        X_Xprime_vector = Xsoln.u
    else
        # Use the supplied rgrid and perform interpolation/extrapolation if needed
        rgrid = rspan
        rsgrid = (x -> rstar_from_r(a, x)).(rgrid)
        X_Xprime_vector = Xsoln.(rsgrid)
    end
    # Now X and Xprime are function of r, i.e. X = X(r), Xprime = dX/dr = dX/dr_* dr_*/dr_
    X = [X_Xprime_vector[i][1] for i in 1:length(rsgrid)]
    drstar_dr = (x -> (x^2 + a^2)/Delta(a, x)).(rgrid)
    Xprime = [X_Xprime_vector[i][2] for i in 1:length(rsgrid)] .* drstar_dr

    # Now transform from X, Xprime to chi, chiprime
    chi_conversion_factor(r) = 1.0/sqrt((r^2 + a^2) * Delta(a, r)^s)
    dchi_conversion_factor_dr(r) = begin
        ((-a^2)*(r - s + r*s) + r^2*(2 + s - r*(1 + s)))/
        ((a^2 + (-2 + r)*r)*sqrt((a^2 + (-2 + r)*r)^s)*(a^2 + r^2)^(3/2))
    end
    _chi_conv_array = chi_conversion_factor.(rgrid)
    _dchi_conv_dr_array = dchi_conversion_factor_dr.(rgrid)
    chi = X .* _chi_conv_array
    chiprime = Xprime .* _chi_conv_array + Xprime .* _dchi_conv_dr_array

    # And then transform from chi, chiprime to R, Rprime
    _eta = (x -> eta(s, m, a, omega, lambda, x)).(rgrid)
    _alpha = (x -> alpha(s, m, a, omega, lambda, x)).(rgrid)
    _alpha_prime = (x -> alpha_prime(s, m, a, omega, lambda, x)).(rgrid)
    _beta = (x -> beta(s, m, a, omega, lambda, x)).(rgrid)
    _beta_prime = (x -> beta_prime(s, m, a, omega, lambda, x)).(rgrid)
    _Delta = (x -> Delta(a, x)).(rgrid)
    _VT = (x -> VT(s, m, a, omega, lambda, x)).(rgrid)

    _inv_eta = 1.0 ./ _eta

    R = _inv_eta .* ( (_alpha .+ _beta_prime.*((_Delta).^(s+1))).*chi + (-_beta.*((_Delta).^(s+1))).*chiprime )
    Rprime = _inv_eta .* ( -(_alpha_prime .+ _beta.*_VT.*((_Delta).^(s))).*chi + (_alpha).*chiprime )

    return rgrid, rsgrid, R, Rprime
end