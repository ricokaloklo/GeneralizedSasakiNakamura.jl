module ComplexFrequencies

using ..Kerr
using ..Coordinates
using ..Potentials
using ..AsymptoticExpansionCoefficients
using ..InitialConditions

using DifferentialEquations

_DEFAULTDATATYPE = ComplexF64 # Double precision by default
_DEFAULTSOLVER = Vern9()
_DEFAULTTOLERANCE = 1e-12

function drdrho!(du, u, p, rho)
    #=
    u[1] = r
    du[1] = dr/drho = exp(i\beta) \Delta/(r^2 + a^2)

    p should be a named tuple with the following fields: (a, beta)

    NOTE This is an autonomous equation, so the RHS does not depend on rho
    =#
    du[1] = exp(1im*p.beta) * Delta(p.a, u[1]) / (u[1]^2 + p.a^2)
end

function solve_r_from_rho(
    a, beta, rs_mp, rho_end;
    dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER,
    reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE
)
    p = (a=a, beta=beta)
    # Initial condition at rho = 0
    r0 = r_from_rstar(a, rs_mp)
    u0 = [dtype(r0)]

    rhospan = (0, rho_end)
    odeprob = ODEProblem(drdrho!, u0, rhospan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    r_from_rho(rho) = odesoln(rho)[1] # Flatten the output
    return r_from_rho
end

# Power of multiple dispatch
function solve_r_from_rho(
    a, beta_neg, beta_pos,
    rs_mp, rho_neg_end, rho_pos_end;
    dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER,
    reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE
)
    # Obtain r_from_rho for positive rho
    r_from_rhopos = solve_r_from_rho(
        a, beta_pos, rs_mp, rho_pos_end;
        dtype=dtype, odealgo=odealgo,
        reltol=reltol, abstol=abstol
    )

    # Obtain r_from_rho for negative rho
    r_from_rhoneg = solve_r_from_rho(
        a, beta_neg, rs_mp, rho_neg_end;
        dtype=dtype, odealgo=odealgo,
        reltol=reltol, abstol=abstol
    )

    # Stitch together the two solutions
    # This is not very readable but it works!
    r_from_rho(rho) = rho >= 0 ? rho <= rho_pos_end ? r_from_rhopos(rho) : NaN : rho >= rho_neg_end ? r_from_rhoneg(rho) : NaN
    return r_from_rho
end

function GSN_eqn_in_rho!(du, u, p, rho)
    r = p.r_from_rho(rho)
    _sF = exp(1im*p.beta)*sF(p.s, p.m, p.a, p.omega, p.lambda, r)
    _sU = exp(2im*p.beta)*sU(p.s, p.m, p.a, p.omega, p.lambda, r)

    #=
    Using the convention for DifferentialEquations
    u[1] = X(rho)
    u[2] = dX/drho = X'
    therefore
    X'' - sF X' - sU X = 0 => u[2]' - sF u[2] - sU u[1] = 0 => u[2]' = sF u[2] + sU u[1]
    =#
    du[1] = u[2]
    du[2] = _sF*u[2] + _sU*u[1]
end

function Xup_in_rho_initialconditions(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rs_mp, rhoout; order::Int=-1)
    outgoing_coeff_func(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
    fout(r) = fansatz(outgoing_coeff_func, omega, r; order=order)
    dfout_dr(r) = dfansatz_dr(outgoing_coeff_func, omega, r; order=order)

    rout = r_from_rho(rhoout)
    _fansatz = fout(rout)
    _dfansatz_dr = dfout_dr(rout)
    rsout = rs_mp + rhoout * exp(1im*beta)
    phase = exp(1im * omega * rsout)

    Xrho_out = phase*_fansatz
    dXdrho_out = exp(1im*beta)*phase*(1im*omega*_fansatz + (Delta(a, rout)/(rout^2 + a^2))*_dfansatz_dr)

    return Xrho_out, dXdrho_out
end

function solve_X_in_rho(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rhospan, initial_conditions; dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    p = (s=s, m=m, a=a, beta=beta, omega=omega, lambda=lambda, r_from_rho=r_from_rho)

    odeprob = ODEProblem(GSN_eqn_in_rho!, initial_conditions, rhospan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    return odesoln
end

function solve_Xup_in_rho(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin, rhoout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rhoin > rhoout
        throw(DomainError(rhoout, "rhoout ($rhoout) must be larger than rhoin ($rhoin)"))
    end
    if rhoout < 0
        throw(DomainError(rhoout, "rhoout ($rhoout) must be positive"))
    end

    # Initial conditions at rho = rhoout, the outer boundary
    Xup_rhoout, Xupprime_rhoout = Xup_in_rho_initialconditions(s, m, a, beta_pos, omega, lambda, r_from_rho, rs_mp, rhoout; order=initialconditions_order)
    u0 = [dtype(Xup_rhoout); dtype(Xupprime_rhoout)]

    # Solve the ODE to the matching point no matter what
    odesoln_pos = solve_X_in_rho(s, m, a, beta_pos, omega, lambda, r_from_rho, (rhoout, 0), u0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)

    if rhoin < 0
        v0 = [dtype(1); dtype(exp(1im*beta_neg))] .* [dtype(1); dtype(exp(-1im*beta_pos))] .* odesoln_pos(0)
        # Continue the integration
        odesoln_neg = solve_X_in_rho(s, m, a, beta_neg, omega, lambda, r_from_rho, (0, rhoin), v0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    else
        # Should not be used
        odesoln_neg = rho -> [NaN; NaN]
    end

    # Stitch together the two solutions
    # NOTE When rho = 0, the returned solution is [X, dX/drs] and not [X, dX/drho] anymore
    odesoln(rho) = rho >= 0 ? rho != 0 ? odesoln_pos(rho) : [dtype(1); dtype(exp(-1im*beta_pos))] .* odesoln_pos(rho) : odesoln_neg(rho)

    return odesoln, odesoln_pos, odesoln_neg
end

function Xin_in_rho_initialconditions(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rs_mp, rhoin; order::Int=-1)
    ingoing_coeff_func(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord)
    gin(r) = gansatz(ingoing_coeff_func, a, r; order=order)
    dgin_dr(r) = dgansatz_dr(ingoing_coeff_func, a, r; order=order)

    rin = r_from_rho(rhoin)
    _gansatz = gin(rin)
    _dgansatz_dr = dgin_dr(rin)
    rsin = rs_mp + rhoin * exp(1im*beta)
    p = omega - m*omega_horizon(a)
    phase = exp(-1im * p * rsin)

    Xrho_in = phase*_gansatz
    dXdrho_in = exp(1im*beta)*phase*(-1im*p*_gansatz + (Delta(a, rin)/(rin^2 + a^2))*_dgansatz_dr)

    return Xrho_in, dXdrho_in
end

function solve_Xin_in_rho(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin, rhoout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rhoin > rhoout
        throw(DomainError(rhoin, "rhoin ($rhoin) must be smaller than rhoout ($rhoout)"))
    end
    if rhoin > 0
        throw(DomainError(rhoin, "rhoin ($rhoin) must be negative"))
    end

    # Initial conditions at rho = rhoin, the inner boundary
    Xin_rhoin, Xinprime_rhoin = Xin_in_rho_initialconditions(s, m, a, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin; order=initialconditions_order)
    u0 = [dtype(Xin_rhoin); dtype(Xinprime_rhoin)]

    # Solve the ODE to the matching point no matter what
    odesoln_neg = solve_X_in_rho(s, m, a, beta_neg, omega, lambda, r_from_rho, (rhoin, 0), u0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)

    if rhoout > 0
        v0 = [dtype(1); dtype(exp(1im*beta_pos))] .* [dtype(1); dtype(exp(-1im*beta_neg))] .* odesoln_neg(0)
        # Continue the integration
        odesoln_pos = solve_X_in_rho(s, m, a, beta_pos, omega, lambda, r_from_rho, (0, rhoout), v0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    else
        # Should not be used
        odesoln_pos = rho -> [NaN; NaN]
    end

    # Stitch together the two solutions
    # NOTE When rho = 0, the returned solution is [X, dX/drs] and not [X, dX/drho] anymore
    odesoln(rho) = rho <= 0 ? rho != 0 ? odesoln_neg(rho) : [dtype(1); dtype(exp(-1im*beta_neg))] .* odesoln_neg(rho) : odesoln_pos(rho)

    return odesoln, odesoln_pos, odesoln_neg
end

end