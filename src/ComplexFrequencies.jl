module ComplexFrequencies

using ..Kerr
using ..Coordinates
using ..Potentials
using ..AsymptoticExpansionCoefficients
using ..InitialConditions
using ..Solutions

using DifferentialEquations
using Optimization
using OptimizationOptimJL
using StaticArrays

_DEFAULTDATATYPE = Solutions._DEFAULTDATATYPE
_DEFAULTSOLVER = Solutions._DEFAULTSOLVER
_DEFAULTTOLERANCE = Solutions._DEFAULTTOLERANCE

function determine_sign(x)
    if real(x) >= 0
        return sign(1.0)
    else
        return sign(-1.0)
    end
end

function drdrho(u, p, rho)
    #=
    u[1] = r
    du[1] = dr/drho = exp(i\beta) \Delta/(r^2 + a^2)

    p should be a named tuple with the following fields: (a, beta)

    NOTE This is an autonomous equation, so the RHS does not depend on rho
    =#
    return SA[p.sign*exp(1im*p.beta) * Delta(p.a, u[1]) / (u[1]^2 + p.a^2)]
end

function solve_r_from_rho(
    a, beta, rs_mp, rho_end; sign=1,
    dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER,
    reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE
)
    p = (a=a, beta=beta, sign=sign)
    # Initial condition at rho = 0
    r0 = r_from_rstar(a, rs_mp)
    u0 = SA[dtype(r0)]

    rhospan = (0, rho_end)
    odeprob = ODEProblem(drdrho, u0, rhospan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    r_from_rho(rho) = odesoln(rho)[1] # Flatten the output
    return r_from_rho
end

# Power of multiple dispatch
function solve_r_from_rho(
    a, beta_neg, beta_pos,
    rs_mp, rho_neg_end, rho_pos_end;
    sign_neg=1, sign_pos=1,
    dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER,
    reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE
)
    # Obtain r_from_rho for positive rho
    r_from_rhopos = solve_r_from_rho(
        a, beta_pos, rs_mp, rho_pos_end; sign=sign_pos,
        dtype=dtype, odealgo=odealgo,
        reltol=reltol, abstol=abstol
    )

    # Obtain r_from_rho for negative rho
    r_from_rhoneg = solve_r_from_rho(
        a, beta_neg, rs_mp, rho_neg_end; sign=sign_neg,
        dtype=dtype, odealgo=odealgo,
        reltol=reltol, abstol=abstol
    )

    # Stitch together the two solutions
    r_from_rho(rho) = rho >= 0 ? rho <= rho_pos_end ? r_from_rhopos(rho) : NaN : rho >= rho_neg_end ? r_from_rhoneg(rho) : NaN
    return r_from_rho
end

function solve_r_from_rho(
    s, m, a, omega, lambda,
    beta_neg, beta_pos,
    rho_neg_end, rho_pos_end;
    sign_neg=1, sign_pos=1,
    dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER,
    reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE
)
    # Some helper functions
    p = omega - m*omega_horizon(a)

    solve_r_from_rho_rsmp(rsmp) = solve_r_from_rho(
        a, beta_neg, beta_pos, rsmp, rho_neg_end, rho_pos_end;
        sign_neg=sign_neg, sign_pos=sign_pos,
        dtype=dtype, odealgo=odealgo,
        reltol=reltol, abstol=abstol
    )

    # The loss function, to be used in the optimization.
    function asymptotic_behavior_of_sFsU(u, _)
        #=
        This function solves for r(\rho) and then checks if the asymptotic behavior of
        sF and sU in the limit \rho \to \pm \infty are satisfied.

        We want
            sF \to 0 as \rho \to \pm \infty
            sU \to -ω^2 as \rho \to \infty, and sU -> -p^2 as \rho \to -\infty
        =#
        # Solve for r(ρ)
        sol = solve_r_from_rho_rsmp(u[1])

        # Get the asymptotic behavior of sF and sU as \rho \to \infty
        sF_rhoout = sF(s, m, a, omega, lambda, sol(rho_pos_end))
        sU_rhoout = sU(s, m, a, omega, lambda, sol(rho_pos_end)) + omega^2
        
        # Get the asymptotic behavior of sF and sU as \rho \to -\infty
        # NOTE Since r approaches r_+ counterclockwise in the complex r plane,
        # we check the "average" behavior for convergence
        rhos = rho_neg_end:0.1:rho_neg_end + abs(rho_neg_end)/10
        sF_rhoin = sum(sF.(s, m, a, omega, lambda, sol.(rhos)))/length(rhos)
        sU_rhoin = sum(sU.(s, m, a, omega, lambda, sol.(rhos)))/length(rhos) + p^2
    
        return abs(sF_rhoout^2 + sF_rhoin^2 + sU_rhoout^2 + sU_rhoin^2)
    end

    optim_func = OptimizationFunction(asymptotic_behavior_of_sFsU, Optimization.AutoFiniteDiff())
    optim_prob = OptimizationProblem(optim_func, [0.0], (), lb = [rho_neg_end/10], ub = [rho_pos_end/100])
    optim_soln = solve(optim_prob, BFGS(); maxiters=50)

    # Check if the optimization was successful
    if optim_soln.retcode == ReturnCode.Success
        return solve_r_from_rho_rsmp(optim_soln.u[1]), optim_soln.u[1]
    else
        # If the optimization was not successful, use rsmp = 0
        return solve_r_from_rho_rsmp(0), 0
    end
end

function GSN_linear_eqn(u, p, rho)
    r = p.r_from_rho(rho)
    _sF = p.sign * exp(1im*p.beta)*sF(p.s, p.m, p.a, p.omega, p.lambda, r)
    _sU = exp(2im*p.beta)*sU(p.s, p.m, p.a, p.omega, p.lambda, r)

    #=
    Using the convention for DifferentialEquations
    u[1] = X(rho)
    u[2] = dX/drho = X'
    therefore
    X'' - sF X' - sU X = 0 => u[2]' - sF u[2] - sU u[1] = 0 => u[2]' = sF u[2] + sU u[1]
    =#
    return SA[u[2], _sF*u[2] + _sU*u[1]]
end

function GSN_Riccati_eqn(u, p, rho)
    r = p.r_from_rho(rho)
    _sF = p.sign * exp(1im*p.beta)*sF(p.s, p.m, p.a, p.omega, p.lambda, r)
    _sU = exp(2im*p.beta)*sU(p.s, p.m, p.a, p.omega, p.lambda, r)
    
    #=
    We write X = exp(I*Phi)
    Substitute X in this form into the GSN equation will give
    a Riccati equation, a first-order non-linear equation

    u[1] = Phi
    u[2] = dPhidrs
    =#
    return SA[u[2], -1im*_sU + _sF*u[2] - 1im*u[2]*u[2]]
end

function solve_X_in_rho(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, sign, rhospan, initial_conditions; dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    p = (s=s, m=m, a=a, beta=beta, omega=omega, lambda=lambda, sign=sign, r_from_rho=r_from_rho)

    odeprob = ODEProblem(GSN_linear_eqn, initial_conditions, rhospan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    return odesoln
end

function solve_Phi_in_rho(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, sign, rhospan, initial_conditions; dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    p = (s=s, m=m, a=a, beta=beta, omega=omega, lambda=lambda, sign=sign, r_from_rho=r_from_rho)

    odeprob = ODEProblem(GSN_Riccati_eqn, initial_conditions, rhospan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    return odesoln
end

function Xup_initialconditions(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rs_mp, rhoout; order::Int=-1, dtype=_DEFAULTDATATYPE)
    outgoing_coeff_func(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord; data_type=dtype)
    fout(r) = fansatz(outgoing_coeff_func, omega, r; order=order)
    dfout_dr(r) = dfansatz_dr(outgoing_coeff_func, omega, r; order=order)

    rout = r_from_rho(rhoout)
    _fansatz = fout(rout)
    _dfansatz_dr = dfout_dr(rout)

    phase = exp(1im * abs(omega) * determine_sign(omega)*rhoout) * exp(1im * omega * rs_mp)

    Xrho_out = phase*_fansatz
    dXdrho_out = determine_sign(omega)*exp(1im*beta)*phase*(1im*omega*_fansatz + (Delta(a, rout)/(rout^2 + a^2))*_dfansatz_dr)

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
    Xup_rhoout, Xupprime_rhoout = Xup_initialconditions(s, m, a, beta_pos, omega, lambda, r_from_rho, rs_mp, rhoout; order=initialconditions_order, dtype=dtype)
    u0 = SA[dtype(Xup_rhoout); dtype(Xupprime_rhoout)]

    # Solve the ODE to the matching point no matter what
    odesoln_pos = solve_X_in_rho(s, m, a, beta_pos, omega, lambda, r_from_rho, determine_sign(omega), (rhoout, 0), u0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)

    if rhoin < 0
        v0 = SA[dtype(1); dtype(determine_sign(omega - m*omega_horizon(a))*exp(1im*beta_neg))] .* SA[dtype(1); dtype(determine_sign(omega)*exp(-1im*beta_pos))] .* odesoln_pos(0)
        # Continue the integration
        odesoln_neg = solve_X_in_rho(s, m, a, beta_neg, omega, lambda, r_from_rho, determine_sign(omega - m*omega_horizon(a)), (0, rhoin), v0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    else
        # Should not be used
        odesoln_neg = rho -> SA[NaN; NaN]
    end

    # Stitch together the two solutions
    # NOTE This *always* return X(rho) and dX(rho)/drs
    odesoln(rho) = rho >= 0 ? SA[dtype(1); dtype(determine_sign(omega)*exp(-1im*beta_pos))] .* odesoln_pos(rho) : SA[dtype(1); dtype(determine_sign(omega - m*omega_horizon(a))*exp(-1im*beta_neg))] .* odesoln_neg(rho)

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
    Xup_rhoout, Xupprime_rhoout = Xup_initialconditions(s, m, a, beta_pos, omega, lambda, r_from_rho, rs_mp, rhoout; order=initialconditions_order, dtype=dtype)
    # Convert initial conditions for Xup for Phi
    Phi, Phiprime = Solutions.PhiPhiprime_from_XXprime(Xup_rhoout, Xupprime_rhoout)
    u0 = SA[dtype(Phi); dtype(Phiprime)]

    # Solve the ODE to the matching point no matter what
    odesoln_pos = solve_Phi_in_rho(s, m, a, beta_pos, omega, lambda, r_from_rho, determine_sign(omega), (rhoout, 0), u0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)

    if rhoin < 0
        _Phi, _Phiprime = odesoln_pos(0)
        _X, _Xprime = Solutions.XXprime_from_PhiPhiprime(_Phi, _Phiprime)
        v0_Phi, v0_Phiprime = Solutions.PhiPhiprime_from_XXprime(
            _X,
            dtype(determine_sign(omega - m*omega_horizon(a))*exp(1im*beta_neg)) * dtype(determine_sign(omega)*exp(-1im*beta_pos)) * _Xprime
        )
        # Continue the integration
        odesoln_neg = solve_Phi_in_rho(s, m, a, beta_neg, omega, lambda, r_from_rho, determine_sign(omega - m*omega_horizon(a)), (0, rhoin), SA[v0_Phi; v0_Phiprime]; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    else
        # Should not be used
        odesoln_neg = rho -> SA[NaN; NaN]
    end

    # Stitch together the two solutions
    # NOTE This *always* return Phi(rho) and dPhi(rho)/drs
    odesoln(rho) = rho >= 0 ? SA[dtype(1); dtype(determine_sign(omega)*exp(-1im*beta_pos))] .* odesoln_pos(rho) : SA[dtype(1); dtype(determine_sign(omega - m*omega_horizon(a))*exp(-1im*beta_neg))] .* odesoln_neg(rho)

    return odesoln, odesoln_pos, odesoln_neg
end

function Xin_initialconditions(s::Int, m::Int, a, beta, omega, lambda, r_from_rho, rs_mp, rhoin; order::Int=-1, dtype=_DEFAULTDATATYPE)
    ingoing_coeff_func(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord; data_type=dtype)
    gin(r) = gansatz(ingoing_coeff_func, a, r; order=order)
    dgin_dr(r) = dgansatz_dr(ingoing_coeff_func, a, r; order=order)

    rin = r_from_rho(rhoin)
    _gansatz = gin(rin)
    _dgansatz_dr = dgin_dr(rin)
    # rsin = rs_mp + rhoin * exp(1im*beta)
    p = omega - m*omega_horizon(a)
    phase = exp(-1im * abs(p) * determine_sign(p)*rhoin) * exp(-1im * p * rs_mp)

    Xrho_in = phase*_gansatz
    dXdrho_in = determine_sign(p)*exp(1im*beta)*phase*(-1im*p*_gansatz + (Delta(a, rin)/(rin^2 + a^2))*_dgansatz_dr)

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
    Xin_rhoin, Xinprime_rhoin = Xin_initialconditions(s, m, a, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin; order=initialconditions_order, dtype=dtype)
    u0 = SA[dtype(Xin_rhoin); dtype(Xinprime_rhoin)]

    # Solve the ODE to the matching point no matter what
    odesoln_neg = solve_X_in_rho(s, m, a, beta_neg, omega, lambda, r_from_rho, determine_sign(omega - m*omega_horizon(a)), (rhoin, 0), u0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)

    if rhoout > 0
        v0 = SA[dtype(1); dtype(determine_sign(omega)*exp(1im*beta_pos))] .* SA[dtype(1); dtype(determine_sign(omega - m*omega_horizon(a))*exp(-1im*beta_neg))] .* odesoln_neg(0)
        # Continue the integration
        odesoln_pos = solve_X_in_rho(s, m, a, beta_pos, omega, lambda, r_from_rho, determine_sign(omega), (0, rhoout), v0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    else
        # Should not be used
        odesoln_pos = rho -> SA[NaN; NaN]
    end

    # Stitch together the two solutions
    # NOTE This *always* return X(rho) and dX(rho)/drs
    odesoln(rho) = rho > 0 ? SA[dtype(1); dtype(determine_sign(omega)*exp(-1im*beta_pos))] .* odesoln_pos(rho) : SA[dtype(1); dtype(determine_sign(omega - m*omega_horizon(a))*exp(-1im*beta_neg))] .* odesoln_neg(rho)

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
    Xin_rhoin, Xinprime_rhoin = Xin_initialconditions(s, m, a, beta_neg, omega, lambda, r_from_rho, rs_mp, rhoin; order=initialconditions_order, dtype=type)
    # Convert initial conditions for Xin for Phi
    Phi, Phiprime = Solutions.PhiPhiprime_from_XXprime(Xin_rhoin, Xinprime_rhoin)
    u0 = SA[dtype(Phi); dtype(Phiprime)]

    # Solve the ODE to the matching point no matter what
    odesoln_neg = solve_Phi_in_rho(s, m, a, beta_neg, omega, lambda, r_from_rho, determine_sign(omega - m*omega_horizon(a)), (rhoin, 0), u0; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    
    if rhoout > 0
        _Phi, _Phiprime = odesoln_neg(0)
        _X, _Xprime = Solutions.XXprime_from_PhiPhiprime(_Phi, _Phiprime)
        v0_Phi, v0_Phiprime = Solutions.PhiPhiprime_from_XXprime(
            _X,
            dtype(determine_sign(omega)*exp(1im*beta_pos)) * dtype(determine_sign(omega - m*omega_horizon(a))*exp(-1im*beta_neg)) * _Xprime
        )
        # Continue the integration
        odesoln_pos = solve_Phi_in_rho(s, m, a, beta_pos, omega, lambda, r_from_rho, determine_sign(omega), (0, rhoout), SA[v0_Phi; v0_Phiprime]; dtype=dtype, odealgo=odealgo, reltol=reltol, abstol=abstol)
    else
        # Should not be used
        odesoln_pos = rho -> SA[NaN; NaN]
    end

    # Stitch together the two solutions
    # NOTE This *always* return Phi(rho) and dPhi(rho)/drs
    odesoln(rho) = rho > 0 ? SA[dtype(1); dtype(determine_sign(omega)*exp(-1im*beta_pos))] .* odesoln_pos(rho) : SA[dtype(1); dtype(determine_sign(omega - m*omega_horizon(a))*exp(-1im*beta_neg))] .* odesoln_neg(rho)

    return odesoln, odesoln_pos, odesoln_neg
end

function BrefBinc_SN_from_Xin(s::Int, m::Int, a, beta, omega, lambda, Xinsoln, r_from_rho, rs_mp, rhoout; order=10, dtype=_DEFAULTDATATYPE)
    rout = r_from_rho(rhoout)

    ingoing_coeff_func(ord) = ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord; data_type=dtype)
    fin(r) = fansatz(ingoing_coeff_func, omega, r; order=order)
    dfin_dr(r) = dfansatz_dr(ingoing_coeff_func, omega, r; order=order)
    _fin = fin(rout)
    _dfin_dr = dfin_dr(rout)
    _phase_in = exp(-1im * abs(omega) * determine_sign(omega)*rhoout) * exp(-1im * omega * rs_mp)

    outgoing_coeff_func(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord; data_type=dtype)
    fout(r) = fansatz(outgoing_coeff_func, omega, r; order=order)
    dfout_dr(r) = dfansatz_dr(outgoing_coeff_func, omega, r; order=order)
    _fout = fout(rout)
    _dfout_dr = dfout_dr(rout)
    _phase_out = exp(1im * abs(omega) * determine_sign(omega)*rhoout) * exp(1im * omega * rs_mp)

    # Computing A1, A2, A3, A4
    A1 = _fin * _phase_in
    A2 = _fout * _phase_out
    A3 = determine_sign(omega)*exp(1im*beta)*_phase_in*(-1im*omega*_fin + (Delta(a, rout)/(rout^2 + a^2))*_dfin_dr)
    A4 = determine_sign(omega)*exp(1im*beta)*_phase_out*(1im*omega*_fout + (Delta(a, rout)/(rout^2 + a^2))*_dfout_dr)

    C1 = Xinsoln(rhoout)[1]
    C2 = determine_sign(omega)*exp(1im*beta)*Xinsoln(rhoout)[2]

    return -(-A3*C1 + A1*C2)/(A2*A3 - A1*A4), -(A4*C1 - A2*C2)/(A2*A3 - A1*A4)
end

function CrefCinc_SN_from_Xup(s::Int, m::Int, a, beta, omega, lambda, Xupsoln, r_from_rho, rs_mp, rhoin; order=10, dtype=_DEFAULTDATATYPE)
    p = omega - m*omega_horizon(a)
    rin = r_from_rho(rhoin)
    # rsin = rs_mp \pm rhoin * exp(1im*beta)

    ingoing_coeff_func(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord; data_type=dtype)
    gin(r) = gansatz(ingoing_coeff_func, a, r; order=order)
    dgin_dr(r) = dgansatz_dr(ingoing_coeff_func, a, r; order=order)
    _gin = gin(rin)
    _dgin_dr = dgin_dr(rin)
    _phase_in = exp(-1im * abs(p) * determine_sign(p)*rhoin) * exp(-1im * p * rs_mp)

    outgoing_coeff_func(ord) = outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord; data_type=dtype)
    gout(r) = gansatz(outgoing_coeff_func, a, r; order=order)
    dgout_dr(r) = dgansatz_dr(outgoing_coeff_func, a, r; order=order)
    _gout = gout(rin)
    _dgout_dr = dgout_dr(rin)
    _phase_out = exp(1im * abs(p) * determine_sign(p)*rhoin) * exp(1im * p * rs_mp)

    # Computing A1, A2, A3, A4
    A1 = _gin * _phase_in
    A2 = _gout * _phase_out
    A3 = determine_sign(p)*exp(1im*beta)*_phase_in*(-1im*p*_gin + (Delta(a, rin)/(rin^2 + a^2))*_dgin_dr)
    A4 = determine_sign(p)*exp(1im*beta)*_phase_out*(1im*p*_gout + (Delta(a, rin)/(rin^2 + a^2))*_dgout_dr)

    C1 = Xupsoln(rhoin)[1]
    C2 = determine_sign(p)*exp(1im*beta)*Xupsoln(rhoin)[2]

    return -(A4*C1 - A2*C2)/(A2*A3 - A1*A4), -(-A3*C1 + A1*C2)/(A2*A3 - A1*A4)
end

function semianalytical_Xin(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, Xinsoln, r_from_rho, rs_mp, rhoin, rhoout, horizon_expansionorder, infinity_expansionorder, rho; dtype=_DEFAULTDATATYPE)
    if rho < rhoin
        # Extend the numerical solution to the analytical ansatz from rhoin to horizon
        p = omega - m*omega_horizon(a)
        _r = r_from_rho(rho)
        if isnan(_r)
            # Resolve r_from_rho
            r_from_rho = solve_r_from_rho(a, beta_neg, beta_pos, rs_mp, rho, rhoout; sign_neg=determine_sign(p), sign_pos=determine_sign(omega))
            _r = r_from_rho(rho)
        end
        # _rs = rs_mp + rho * exp(1im*beta_neg)

        # Construct the analytical ansatz
        ingoing_coeff_func_hor(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord; data_type=dtype)
        gin(r) = gansatz(ingoing_coeff_func_hor, a, r; order=horizon_expansionorder)
        dgin_dr(r) = dgansatz_dr(ingoing_coeff_func_hor, a, r; order=horizon_expansionorder)
        _gin = gin(_r)
        _dgin_dr = dgin_dr(_r)
        _phase_in = exp(-1im * abs(p) * determine_sign(p)* rho) * exp(-1im * p * rs_mp)

        # NOTE There is no exp(1im*beta) since the derivative is wrt rstar
        _Xin = _gin * _phase_in
        _dXin_drs = _phase_in*(-1im*p*_gin + (Delta(a, _r)/(_r^2 + a^2))*_dgin_dr)

        return (_Xin, _dXin_drs)

    elseif rho > rhoout
        # Extend the numerical solution to the analytical ansatz from rhoout to infinity

        Bref_SN, Binc_SN = BrefBinc_SN_from_Xin(s, m, a, beta_pos, omega, lambda, Xinsoln, r_from_rho, rs_mp, rhoout; order=infinity_expansionorder, dtype=dtype)

        p = omega - m*omega_horizon(a)
        _r = r_from_rho(rho)
        if isnan(_r)
            # Resolve r_from_rho
            r_from_rho = solve_r_from_rho(a, beta_neg, beta_pos, rs_mp, rhoin, rho; sign_neg=determine_sign(p), sign_pos=determine_sign(omega))
            _r = r_from_rho(rho)
        end
        # _rs = rs_mp + rho * exp(1im*beta_pos)

        # Construct the analytical ansatz
        ingoing_coeff_func(ord) = ingoing_coefficient_at_inf(s, m, a, omega, lambda, ord; data_type=dtype)
        fin(r) = fansatz(ingoing_coeff_func, omega, r; order=infinity_expansionorder)
        dfin_dr(r) = dfansatz_dr(ingoing_coeff_func, omega, r; order=infinity_expansionorder)
        _fin = fin(_r)
        _dfin_dr = dfin_dr(_r)
        _phase_in = exp(-1im * abs(omega) * determine_sign(omega)*rho) * exp(-1im * omega * rs_mp)

        outgoing_coeff_func(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord; data_type=dtype)
        fout(r) = fansatz(outgoing_coeff_func, omega, r; order=infinity_expansionorder)
        dfout_dr(r) = dfansatz_dr(outgoing_coeff_func, omega, r; order=infinity_expansionorder)
        _fout = fout(_r)
        _dfout_dr = dfout_dr(_r)
        _phase_out = exp(1im * abs(omega) * determine_sign(omega)*rho) * exp(1im * omega * rs_mp)

        _Xin = Bref_SN*_fout*_phase_out + Binc_SN*_fin*_phase_in
        _dXin_drs = Bref_SN*_phase_out*(1im*omega*_fout + (Delta(a, _r)/(_r^2 + a^2))*_dfout_dr) + Binc_SN*_phase_in*(-1im*omega*_fin + (Delta(a, _r)/(_r^2 + a^2))*_dfin_dr)

        return (_Xin, _dXin_drs)
    else
        # Return the numerical solution
        return Xinsoln(rho)
    end
end

function semianalytical_Xup(s::Int, m::Int, a, beta_pos, beta_neg, omega, lambda, Xupsoln, r_from_rho, rs_mp, rhoin, rhoout, horizon_expansionorder, infinity_expansionorder, rho; dtype=_DEFAULTDATATYPE)
    if rho < rhoin
        # Extend the numerical solution to the analytical ansatz from rhoin to horizon

        Cref_SN, Cinc_SN = CrefCinc_SN_from_Xup(s, m, a, beta_neg, omega, lambda, Xupsoln, r_from_rho, rs_mp, rhoin; order=horizon_expansionorder, dtype=dtype)

        p = omega - m*omega_horizon(a)
        _r = r_from_rho(rho)
        if isnan(_r)
            # Resolve r_from_rho
            r_from_rho = solve_r_from_rho(a, beta_neg, beta_pos, rs_mp, rho, rhoout; sign_neg=determine_sign(p), sign_pos=determine_sign(omega))
            _r = r_from_rho(rho)
        end
        # _rs = rs_mp + rho * exp(1im*beta_neg)

        # Construct the analytical ansatz
        ingoing_coeff_func(ord) = ingoing_coefficient_at_hor(s, m, a, omega, lambda, ord; data_type=dtype)
        gin(r) = gansatz(ingoing_coeff_func, a, r; order=horizon_expansionorder)
        dgin_dr(r) = dgansatz_dr(ingoing_coeff_func, a, r; order=horizon_expansionorder)
        _gin = gin(_r)
        _dgin_dr = dgin_dr(_r)
        _phase_in = exp(-1im * abs(p) * determine_sign(p)* rho) * exp(-1im * p * rs_mp)

        outgoing_coeff_func(ord) = outgoing_coefficient_at_hor(s, m, a, omega, lambda, ord; data_type=dtype)
        gout(r) = gansatz(outgoing_coeff_func, a, r; order=horizon_expansionorder)
        dgout_dr(r) = dgansatz_dr(outgoing_coeff_func, a, r; order=horizon_expansionorder)
        _gout = gout(_r)
        _dgout_dr = dgout_dr(_r)
        _phase_out = exp(1im * abs(p) * determine_sign(p)* rho) * exp(1im * p * rs_mp)

        _Xup = Cref_SN*_gin*_phase_in + Cinc_SN*_gout*_phase_out
        _dXup_drs = Cref_SN*_phase_in*(-1im*p*_gin + (Delta(a, _r)/(_r^2 + a^2))*_dgin_dr) + Cinc_SN*_phase_out*(1im*p*_gout + (Delta(a, _r)/(_r^2 + a^2))*_dgout_dr)

        return (_Xup, _dXup_drs)
    elseif rho > rhoout
        # Extend the numerical solution to the analytical ansatz from rhoout to infinity
        p = omega - m*omega_horizon(a)
        _r = r_from_rho(rho)
        if isnan(_r)
            # Resolve r_from_rho
            r_from_rho = solve_r_from_rho(a, beta_neg, beta_pos, rs_mp, rhoin, rho; sign_neg=determine_sign(p), sign_pos=determine_sign(omega))
            _r = r_from_rho(rho)
        end

        # Construct the analytical ansatz
        outgoing_coeff_func_inf(ord) = outgoing_coefficient_at_inf(s, m, a, omega, lambda, ord; data_type=dtype)
        fout(r) = fansatz(outgoing_coeff_func_inf, omega, r; order=infinity_expansionorder)
        dfout_dr(r) = dfansatz_dr(outgoing_coeff_func_inf, omega, r; order=infinity_expansionorder)
        _fout = fout(_r)
        _dfout_dr = dfout_dr(_r)
        _phase_out = exp(1im * abs(omega) * determine_sign(omega)* rho) * exp(1im * omega * rs_mp)

        _Xup = _fout*_phase_out
        _dXup_drs = _phase_out*(1im*omega*_fout + (Delta(a, _r)/(_r^2 + a^2))*_dfout_dr)

        return (_Xup, _dXup_drs)
    else
        # Return the numerical solution
        return Xupsoln(rho)
    end
end

end