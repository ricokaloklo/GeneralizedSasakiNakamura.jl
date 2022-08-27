using OrdinaryDiffEq

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

function solve_Xinf(s::Int, m::Int, a, omega, lambda, rsin, rsout; reltol=1e-10, abstol=1e-10)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsout, "rsout ($rsout) must be larger than rsin ($rsin)"))
    end
    # Initial conditions at rs = rsout, the outer boundary
    Xinf_rsout, Xinfprime_rsout = Xinf_initialconditions(s, m, a, omega, lambda, rsout)
    u0 = [Xinf_rsout; Xinfprime_rsout]
    rsspan = (rsout, rsin) # Integrate from rsout to rsin *inward*
    p = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    odeprob = ODEProblem(generalized_sasaki_nakamura_equation!, u0, rsspan, p)
    odealgo = RK4()
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
end

function solve_Xhor(s::Int, m::Int, a, omega, lambda, rsin, rsout; reltol=1e-10, abstol=1e-10)
    # Sanity check
    if rsin > rsout
        throw(DomainError(rsin, "rsin ($rsin) must be smaller than rsout ($rsout)"))
    end
    # Initial conditions at rs = rsin, the inner boundary; this should be very close to EH
    Xhor_rsin, Xhorprime_rsin = Xhor_initialconditions(s, m, a, omega, lambda, rsin)
    u0 = [Xhor_rsin; Xhorprime_rsin]
    rsspan = (rsin, rsout) # Integrate from rsin to rsout *outward*
    p = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    odeprob = ODEProblem(generalized_sasaki_nakamura_equation!, u0, rsspan, p)
    odealgo = RK4()
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)
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
    [R(r), dR/dr(r)]^T = 1/eta(r) * [ alpha + beta_prime*Delta^(s+1), -beta*Delta^s ; -(alpha_prime + beta*VT*Delta^s), a ] *  [chi(r), dchi/dr(r)]^T

    Therefore the overall conversion matrix is 'just' (one matrix for each r)
    the multiplication of each conversion matrix
    =#

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
    overall_conversion_matrix(r) = conversion_matrix_from_chi_dchidr_to_R_dRdr(r) * conversion_matrix_from_X_dXdr_to_chi_dchidr(r) * conversion_matrix_from_X_dXdrs_to_X_dXdr(r)

    _interpolant_in_r(r) = Xsoln(rstar_from_r(a, r))
    Rsoln = (r -> overall_conversion_matrix(r) * _interpolant_in_r(r))
    return Rsoln
end