module ComplexFrequencies

using ..Kerr
using ..Coordinates
using ..Potentials
using ..AsymptoticExpansionCoefficients
using ..InitialConditions
using ..Solutions

using DifferentialEquations

_DEFAULTDATATYPE = Solutions._DEFAULTDATATYPE
_DEFAULTSOLVER = Solutions._DEFAULTSOLVER
_DEFAULTTOLERANCE = Solutions._DEFAULTTOLERANCE

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
    r_from_rho(rho) = rho >= 0 ? rho <= rho_pos_end ? r_from_rhopos(rho) : NaN : rho >= rho_neg_end ? r_from_rhoneg(rho) : NaN
    return r_from_rho
end

function GSN_linear_eqn!(du, u, p, rho)
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

function GSN_Riccati_eqn!(du, u, p, rho)
    r = p.r_from_rho(rho)
    _sF = exp(1im*p.beta)*sF(p.s, p.m, p.a, p.omega, p.lambda, r)
    _sU = exp(2im*p.beta)*sU(p.s, p.m, p.a, p.omega, p.lambda, r)
    
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

function solve_X_in_rho(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rhospan, initial_conditions; dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    p = (s=s, m=m, a=a, beta=beta, omega=omega, lambda=lambda, r_from_rho=r_from_rho)

    odeprob = ODEProblem(GSN_linear_eqn!, initial_conditions, rhospan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    return odesoln
end

function solve_Phi_in_rho(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rhospan, initial_conditions; dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    p = (s=s, m=m, a=a, beta=beta, omega=omega, lambda=lambda, r_from_rho=r_from_rho)

    odeprob = ODEProblem(GSN_Riccati_eqn!, initial_conditions, rhospan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    return odesoln
end

function Xup_initialconditions(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rs_mp, rhoout; order::Int=-1)
    outgoing_coeff_func(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
    fout(r) = fansatz(outgoing_coeff_func, omega, r; order=order)
    dfout_dr(r) = dfansatz_dr(outgoing_coeff_func, omega, r; order=order)

    rout = r_from_rho(rhoout)
    _fansatz = fout(rout)
    _dfansatz_dr = dfout_dr(rout)
    # rsout = rs_mp + rhoout * exp(1im*beta)
    phase = exp(1im * abs(omega) * rhoout) * exp(1im * omega * rs_mp)

    Xrho_out = phase*_fansatz
    dXdrho_out = exp(1im*beta)*phase*(1im*omega*_fansatz + (Delta(a, rout)/(rout^2 + a^2))*_dfansatz_dr)

    return Xrho_out, dXdrho_out
end

function solve_Xup(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin, rhoout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rhoin > rhoout
        throw(DomainError(rhoout, "rhoout ($rhoout) must be larger than rhoin ($rhoin)"))
    end
    if rhoout < 0
        throw(DomainError(rhoout, "rhoout ($rhoout) must be positive"))
    end

    # Initial conditions at rho = rhoout, the outer boundary
    Xup_rhoout, Xupprime_rhoout = Xup_initialconditions(s, m, a, beta_pos, omega, lambda, r_from_rho, rs_mp, rhoout; order=initialconditions_order)
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
    # NOTE This *always* return X(rho) and dX(rho)/drs
    odesoln(rho) = rho >= 0 ? [dtype(1); dtype(exp(-1im*beta_pos))] .* odesoln_pos(rho) : [dtype(1); dtype(exp(-1im*beta_neg))] .* odesoln_neg(rho)

    return odesoln, odesoln_pos, odesoln_neg
end

function solve_Phiup(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin, rhoout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rhoin > rhoout
        throw(DomainError(rhoout, "rhoout ($rhoout) must be larger than rhoin ($rhoin)"))
    end
    if rhoout < 0
        throw(DomainError(rhoout, "rhoout ($rhoout) must be positive"))
    end

    # Initial conditions at rho = rhoout, the outer boundary
    Xup_rhoout, Xupprime_rhoout = Xup_initialconditions(s, m, a, beta_pos, omega, lambda, r_from_rho, rs_mp, rhoout; order=initialconditions_order)
    # Convert initial conditions for Xup for Phi
    Phi, Phiprime = Solutions.PhiPhiprime_from_XXprime(Xup_rhoout, Xupprime_rhoout)
    u0 = [dtype(Phi); dtype(Phiprime)]

    # Solve the ODE to the matching point no matter what
    odesoln_pos = solve_Phi_in_rho(s, m, a, beta_pos, omega, lambda, r_from_rho, (rhoout, 0), u0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)

    if rhoin < 0
        _Phi, _Phiprime = odesoln_pos(0)
        _X, _Xprime = Solutions.XXprime_from_PhiPhiprime(_Phi, _Phiprime)
        v0_Phi, v0_Phiprime = Solutions.PhiPhiprime_from_XXprime(
            _X,
            dtype(exp(1im*beta_neg)) * dtype(exp(-1im*beta_pos)) * _Xprime
        )
        # Continue the integration
        odesoln_neg = solve_Phi_in_rho(s, m, a, beta_neg, omega, lambda, r_from_rho, (0, rhoin), [v0_Phi; v0_Phiprime]; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    else
        # Should not be used
        odesoln_neg = rho -> [NaN; NaN]
    end

    # Stitch together the two solutions
    # NOTE This *always* return Phi(rho) and dPhi(rho)/drs
    odesoln(rho) = rho >= 0 ? [dtype(1); dtype(exp(-1im*beta_pos))] .* odesoln_pos(rho) : [dtype(1); dtype(exp(-1im*beta_neg))] .* odesoln_neg(rho)

    return odesoln, odesoln_pos, odesoln_neg
end

function Xin_initialconditions(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rs_mp, rhoin; order::Int=-1)
    ingoing_coeff_func(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord)
    gin(r) = gansatz(ingoing_coeff_func, a, r; order=order)
    dgin_dr(r) = dgansatz_dr(ingoing_coeff_func, a, r; order=order)

    rin = r_from_rho(rhoin)
    _gansatz = gin(rin)
    _dgansatz_dr = dgin_dr(rin)
    # rsin = rs_mp + rhoin * exp(1im*beta)
    p = omega - m*omega_horizon(a)
    phase = exp(-1im * abs(p) * rhoin) * exp(-1im * p * rs_mp)

    Xrho_in = phase*_gansatz
    dXdrho_in = exp(1im*beta)*phase*(-1im*p*_gansatz + (Delta(a, rin)/(rin^2 + a^2))*_dgansatz_dr)

    return Xrho_in, dXdrho_in
end

function solve_Xin(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin, rhoout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rhoin > rhoout
        throw(DomainError(rhoin, "rhoin ($rhoin) must be smaller than rhoout ($rhoout)"))
    end
    if rhoin > 0
        throw(DomainError(rhoin, "rhoin ($rhoin) must be negative"))
    end

    # Initial conditions at rho = rhoin, the inner boundary
    Xin_rhoin, Xinprime_rhoin = Xin_initialconditions(s, m, a, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin; order=initialconditions_order)
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
    # NOTE This *always* return X(rho) and dX(rho)/drs
    odesoln(rho) = rho >= 0 ? [dtype(1); dtype(exp(-1im*beta_pos))] .* odesoln_pos(rho) : [dtype(1); dtype(exp(-1im*beta_neg))] .* odesoln_neg(rho)

    return odesoln, odesoln_pos, odesoln_neg
end

function solve_Phiin(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin, rhoout; initialconditions_order=-1, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Sanity check
    if rhoin > rhoout
        throw(DomainError(rhoin, "rhoin ($rhoin) must be smaller than rhoout ($rhoout)"))
    end
    if rhoin > 0
        throw(DomainError(rhoin, "rhoin ($rhoin) must be negative"))
    end

    # Initial conditions at rho = rhoin, the inner boundary
    Xin_rhoin, Xinprime_rhoin = Xin_initialconditions(s, m, a, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin; order=initialconditions_order)
    # Convert initial conditions for Xin for Phi
    Phi, Phiprime = Solutions.PhiPhiprime_from_XXprime(Xin_rhoin, Xinprime_rhoin)
    u0 = [dtype(Phi); dtype(Phiprime)]

    # Solve the ODE to the matching point no matter what
    odesoln_neg = solve_Phi_in_rho(s, m, a, beta_neg, omega, lambda, r_from_rho, (rhoin, 0), u0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    
    if rhoout > 0
        _Phi, _Phiprime = odesoln_neg(0)
        _X, _Xprime = Solutions.XXprime_from_PhiPhiprime(_Phi, _Phiprime)
        v0_Phi, v0_Phiprime = Solutions.PhiPhiprime_from_XXprime(
            _X,
            dtype(exp(1im*beta_pos)) * dtype(exp(-1im*beta_neg)) * _Xprime
        )
        # Continue the integration
        odesoln_pos = solve_Phi_in_rho(s, m, a, beta_pos, omega, lambda, r_from_rho, (0, rhoout), [v0_Phi; v0_Phiprime]; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    else
        # Should not be used
        odesoln_pos = rho -> [NaN; NaN]
    end

    # Stitch together the two solutions
    # NOTE This *always* return Phi(rho) and dPhi(rho)/drs
    odesoln(rho) = rho >= 0 ? [dtype(1); dtype(exp(-1im*beta_pos))] .* odesoln_pos(rho) : [dtype(1); dtype(exp(-1im*beta_neg))] .* odesoln_neg(rho)

    return odesoln, odesoln_pos, odesoln_neg
end

function BrefBinc_SN_from_Xin(s::Int, m::Int, a, beta, omega, lambda, Xinsoln, r_from_rho, rs_mp, rhoout; order=10)
    rout = r_from_rho(rhoout)
    # rsout = rs_mp + rhoout * exp(1im*beta)

    ingoing_coeff_func(ord) = ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
    fin(r) = fansatz(ingoing_coeff_func, omega, r; order=order)
    dfin_dr(r) = dfansatz_dr(ingoing_coeff_func, omega, r; order=order)
    _fin = fin(rout)
    _dfin_dr = dfin_dr(rout)
    _phase_in = exp(-1im * abs(omega) * rhoout) * exp(-1im * omega * rs_mp)

    outgoing_coeff_func(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
    fout(r) = fansatz(outgoing_coeff_func, omega, r; order=order)
    dfout_dr(r) = dfansatz_dr(outgoing_coeff_func, omega, r; order=order)
    _fout = fout(rout)
    _dfout_dr = dfout_dr(rout)
    _phase_out = exp(1im * abs(omega) * rhoout) * exp(1im * omega * rs_mp)

    # Computing A1, A2, A3, A4
    A1 = _fin * _phase_in
    A2 = _fout * _phase_out
    A3 = exp(1im*beta)*_phase_in*(-1im*omega*_fin + (Delta(a, rout)/(rout^2 + a^2))*_dfin_dr)
    A4 = exp(1im*beta)*_phase_out*(1im*omega*_fout + (Delta(a, rout)/(rout^2 + a^2))*_dfout_dr)

    C1 = Xinsoln(rhoout)[1]
    C2 = exp(1im*beta)*Xinsoln(rhoout)[2]

    return -(-A3*C1 + A1*C2)/(A2*A3 - A1*A4), -(A4*C1 - A2*C2)/(A2*A3 - A1*A4)
end

function CrefCinc_SN_from_Xup(s::Int, m::Int, a, beta, omega, lambda, Xupsoln, r_from_rho, rs_mp, rhoin; order=10)
    p = omega - m*omega_horizon(a)
    rin = r_from_rho(rhoin)
    # rsin = rs_mp + rhoin * exp(1im*beta)

    ingoing_coeff_func(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord)
    gin(r) = gansatz(ingoing_coeff_func, a, r; order=order)
    dgin_dr(r) = dgansatz_dr(ingoing_coeff_func, a, r; order=order)
    _gin = gin(rin)
    _dgin_dr = dgin_dr(rin)
    _phase_in = exp(-1im * abs(p) * rhoin) * exp(-1im * p * rs_mp)

    outgoing_coeff_func(ord) = outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord)
    gout(r) = gansatz(outgoing_coeff_func, a, r; order=order)
    dgout_dr(r) = dgansatz_dr(outgoing_coeff_func, a, r; order=order)
    _gout = gout(rin)
    _dgout_dr = dgout_dr(rin)
    _phase_out = exp(1im * abs(p) * rhoin) * exp(1im * p * rs_mp)

    # Computing A1, A2, A3, A4
    A1 = _gin * _phase_in
    A2 = _gout * _phase_out
    A3 = exp(1im*beta)*_phase_in*(-1im*p*_gin + (Delta(a, rin)/(rin^2 + a^2))*_dgin_dr)
    A4 = exp(1im*beta)*_phase_out*(1im*p*_gout + (Delta(a, rin)/(rin^2 + a^2))*_dgout_dr)

    C1 = Xupsoln(rhoin)[1]
    C2 = exp(1im*beta)*Xupsoln(rhoin)[2]

    return -(A4*C1 - A2*C2)/(A2*A3 - A1*A4), -(-A3*C1 + A1*C2)/(A2*A3 - A1*A4)
end

function semianalytical_Xin(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, Xinsoln, r_from_rho, rs_mp, rhoin, rhoout, horizon_expansionorder, infinity_expansionorder, rho)
    if rho < rhoin
        # Extend the numerical solution to the analytical ansatz from rhoin to horizon
        p = omega - m*omega_horizon(a)
        _r = r_from_rho(rho)
        if isnan(_r)
            # Resolve r_from_rho
            r_from_rho = solve_r_from_rho(a, beta_neg, beta_pos, rs_mp, rho, rhoout)
            _r = r_from_rho(rho)
        end
        # _rs = rs_mp + rho * exp(1im*beta_neg)

        # Construct the analytical ansatz
        ingoing_coeff_func_hor(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord)
        gin(r) = gansatz(ingoing_coeff_func_hor, a, r; order=horizon_expansionorder)
        dgin_dr(r) = dgansatz_dr(ingoing_coeff_func_hor, a, r; order=horizon_expansionorder)
        _gin = gin(_r)
        _dgin_dr = dgin_dr(_r)
        _phase_in = exp(-1im * abs(p) * rho) * exp(-1im * p * rs_mp)

        # NOTE There is no exp(1im*beta) since the derivative is wrt rstar
        _Xin = _gin * _phase_in
        _dXin_drs = _phase_in*(-1im*p*_gin + (Delta(a, _r)/(_r^2 + a^2))*_dgin_dr)

        return (_Xin, _dXin_drs)

    elseif rho > rhoout
        # Extend the numerical solution to the analytical ansatz from rhoout to infinity

        Bref_SN, Binc_SN = BrefBinc_SN_from_Xin(s, m, a, beta_pos, omega, lambda, Xinsoln, r_from_rho, rs_mp, rhoout; order=infinity_expansionorder)
        _r = r_from_rho(rho)
        if isnan(_r)
            # Resolve r_from_rho
            r_from_rho = solve_r_from_rho(a, beta_neg, beta_pos, rs_mp, rhoin, rho)
            _r = r_from_rho(rho)
        end
        # _rs = rs_mp + rho * exp(1im*beta_pos)

        # Construct the analytical ansatz
        ingoing_coeff_func(ord) = ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
        fin(r) = fansatz(ingoing_coeff_func, omega, r; order=infinity_expansionorder)
        dfin_dr(r) = dfansatz_dr(ingoing_coeff_func, omega, r; order=infinity_expansionorder)
        _fin = fin(_r)
        _dfin_dr = dfin_dr(_r)
        _phase_in = exp(-1im * abs(omega) * rho) * exp(-1im * omega * rs_mp)

        outgoing_coeff_func(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
        fout(r) = fansatz(outgoing_coeff_func, omega, r; order=infinity_expansionorder)
        dfout_dr(r) = dfansatz_dr(outgoing_coeff_func, omega, r; order=infinity_expansionorder)
        _fout = fout(_r)
        _dfout_dr = dfout_dr(_r)
        _phase_out = exp(1im * abs(omega) * rho) * exp(1im * omega * rs_mp)

        _Xin = Bref_SN*_fout*_phase_out + Binc_SN*_fin*_phase_in
        _dXin_drs = Bref_SN*_phase_out*(1im*omega*_fout + (Delta(a, _r)/(_r^2 + a^2))*_dfout_dr) + Binc_SN*_phase_in*(-1im*omega*_fin + (Delta(a, _r)/(_r^2 + a^2))*_dfin_dr)

        return (_Xin, _dXin_drs)
    else
        # Return the numerical solution
        return Xinsoln(rho)
    end
end

function semianalytical_Xup(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, Xupsoln, r_from_rho, rs_mp, rhoin, rhoout, horizon_expansionorder, infinity_expansionorder, rho)
    if rho < rhoin
        # Extend the numerical solution to the analytical ansatz from rhoin to horizon

        Cref_SN, Cinc_SN = CrefCinc_SN_from_Xup(s, m, a, beta_neg, omega, lambda, Xupsoln, r_from_rho, rs_mp, rhoin; order=horizon_expansionorder)

        p = omega - m*omega_horizon(a)
        _r = r_from_rho(rho)
        if isnan(_r)
            # Resolve r_from_rho
            r_from_rho = solve_r_from_rho(a, beta_neg, beta_pos, rs_mp, rho, rhoout)
            _r = r_from_rho(rho)
        end
        # _rs = rs_mp + rho * exp(1im*beta_neg)

        # Construct the analytical ansatz
        ingoing_coeff_func(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord)
        gin(r) = gansatz(ingoing_coeff_func, a, r; order=horizon_expansionorder)
        dgin_dr(r) = dgansatz_dr(ingoing_coeff_func, a, r; order=horizon_expansionorder)
        _gin = gin(_r)
        _dgin_dr = dgin_dr(_r)
        _phase_in = exp(-1im * abs(p) * rho) * exp(-1im * p * rs_mp)

        outgoing_coeff_func(ord) = outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord)
        gout(r) = gansatz(outgoing_coeff_func, a, r; order=horizon_expansionorder)
        dgout_dr(r) = dgansatz_dr(outgoing_coeff_func, a, r; order=horizon_expansionorder)
        _gout = gout(_r)
        _dgout_dr = dgout_dr(_r)
        _phase_out = exp(1im * abs(p) * rho) * exp(1im * p * rs_mp)

        _Xup = Cref_SN*_gin*_phase_in + Cinc_SN*_gout*_phase_out
        _dXup_drs = Cref_SN*_phase_in*(-1im*p*_gin + (Delta(a, _r)/(_r^2 + a^2))*_dgin_dr) + Cinc_SN*_phase_out*(1im*p*_gout + (Delta(a, _r)/(_r^2 + a^2))*_dgout_dr)

        return (_Xup, _dXup_drs)
    elseif rho > rhoout
        # Extend the numerical solution to the analytical ansatz from rhoout to infinity
        _r = r_from_rho(rho)
        if isnan(_r)
            # Resolve r_from_rho
            r_from_rho = solve_r_from_rho(a, beta_neg, beta_pos, rs_mp, rhoin, rho)
            _r = r_from_rho(rho)
        end

        # Construct the analytical ansatz
        outgoing_coeff_func_inf(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord)
        fout(r) = fansatz(outgoing_coeff_func_inf, omega, r; order=infinity_expansionorder)
        dfout_dr(r) = dfansatz_dr(outgoing_coeff_func_inf, omega, r; order=infinity_expansionorder)
        _fout = fout(_r)
        _dfout_dr = dfout_dr(_r)
        _phase_out = exp(1im * abs(omega) * rho) * exp(1im * omega * rs_mp)

        _Xup = _fout*_phase_out
        _dXup_drs = _phase_out*(1im*omega*_fout + (Delta(a, _r)/(_r^2 + a^2))*_dfout_dr)

        return (_Xup, _dXup_drs)
    else
        # Return the numerical solution
        return Xupsoln(rho)
    end
end

end