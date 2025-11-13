module SolutionsY

using DifferentialEquations
using StaticArrays

using ..GeneralizedSasakiNakamura: BoundaryCondition, IN, UP
using ..GeneralizedSasakiNakamura: GSN_radial
using ..Coordinates: r_from_rstar
using ..AsymptoticExpansionCoefficientsY
using ..InitialConditionsY

export Y_solution, Ypp_m2, Ypp_p2

_DEFAULTDATATYPE = ComplexF64
_DEFAULTSOLVER = AutoVern9(Rosenbrock23(autodiff=false))
_DEFAULTTOLERANCE = 1e-12
_DEFAULT_rsin = -50
_DEFAULT_rsout = 500

struct YSolutionResult
    basis_type::String
    mode::NamedTuple
    asymptotic::NamedTuple
    solution::NamedTuple
end

function GSN_source_eqn_m2(u, p, rs)
    a = p.a
    ω = p.omega
    m = p.m
    X = p.X
    r = r_from_rstar(a, rs)
    rm = 1 - sqrt(1 - a^2)
    rp = 1 + sqrt(1 - a^2)
    phase = exp(-1im*(ω * rs - a * m * log((r - rm)/(r - rp))/(rm - rp)))

    du1 = u[2]
    du2 = 2*(r^2 - a^2)*u[2] / (r^2 + a^2)^2 +
          X(rs) * (r^2-2*r+a^2)^2 * phase / (r^2 * (r^2 + a^2)^(5/2))

    return SA[du1, du2]
end

function GSN_source_eqn_p2(u, p, rs)
    a = p.a
    ω = p.omega
    m = p.m
    X = p.X
    r = r_from_rstar(a, rs)
    rm = 1 - sqrt(1 - a^2)
    rp = 1 + sqrt(1 - a^2)
    phase = exp(1im*(ω * rs - a * m * log((r - rm)/(r - rp))/(rm - rp)))

    du1 = u[2]
    du2 = 2*(r^2 - a^2)*u[2] / (r^2 + a^2)^2 +
          X(rs) * (r^2-2*r+a^2)^2 * phase / (r^2 * (r^2 + a^2)^(5/2))

    return SA[du1, du2]
end

function Y_solution_in_m2(l::Int, m::Int, a, omega, rsin, rsout, order_inf, order_hor, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Precompute GSN radial solution and initial data
    rin = r_from_rstar(a, rsin)
    rout = r_from_rstar(a, rsout)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(-2, l, m, a, omega, IN, rsin, rsout; horizon_expansion_order=order_hor+2, infinity_expansion_order=order_inf+3, method=method)
    λ = X.mode.lambda
    Binc = X.incidence_amplitude
    Bref = X.reflection_amplitude
    Y0, Yp0 = Y_initial_values_inf_m2(m, a, omega, λ, Binc, Bref, rsout; order=order_inf)
    Y_ingoing_hor_func, Yp_ingoing_hor_func = Y_and_Yp_ingoing_hor_m2(m, a, omega, λ; order=order_hor)
    u0 = SA[dtype(Y0); dtype(Yp0)]
    rsspan = (rsout, rsin)  # Integrate inward
    p = (m=m, a=a, omega=omega, X=X)
    odeprob = ODEProblem(GSN_source_eqn_m2, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    # Get boundary values
    Y_in = odesoln(rsin)[1]
    Yp_in = odesoln(rsin)[2] * (rin^2 + a^2) / (rin^2 - 2 * rin+ a^2)

    # Define final callable function
    function Y_inf(rs)
        r = r_from_rstar(a, rs)
        if rs > rsout
            Y_inf_val, Yp_inf_val = Y_initial_values_inf_m2(m, a, omega, λ, Binc, Bref, rs)
            return Y_inf_val, Yp_inf_val * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        elseif rs < rsin
            Y_hor_val = Y_in - (Y_ingoing_hor_func(r) - Y_ingoing_hor_func(rin)) + Yp_in * (r - rin)
            Yp_hor_val = Yp_in + (Yp_ingoing_hor_func(r) - Yp_ingoing_hor_func(rin))
            return Y_hor_val, Yp_hor_val
        else
            sol = odesoln(rs)
            return sol[1], sol[2] * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        end
    end

    Y_value_hor, Yp_value_hor = Y_inf(-200)

    function Y_hor(rs)
        r = r_from_rstar(a, rs)
        rp = 1 + sqrt(1 - a^2)
        Y = Y_inf(rs)[1]- Yp_value_hor * r + Yp_value_hor * rp - Y_value_hor
        Yp = Y_inf(rs)[2] - Yp_value_hor
        return Y, Yp
    end

    return YSolutionResult(
        "IN",
        (s=-2, l=l, m=m, a=a, omega=omega, lambda=λ),
        (rin=rin, rout=rout, Binc=Binc, Bref=Bref),
        (X=rs -> X(rs), Y_inf=rs -> Y_inf(rs)[1], Yp_inf=rs -> Y_inf(rs)[2], Y_hor=rs -> Y_hor(rs)[1], Yp_hor=rs -> Y_hor(rs)[2])
    )
end

function Y_solution_in_p2(l::Int, m::Int, a, omega, rsin, rsout, order_inf, order_hor, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Precompute GSN radial solution and initial data
    rin = r_from_rstar(a, rsin)
    rout = r_from_rstar(a, rsout)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(2, l, m, a, omega, IN, rsin, rsout; horizon_expansion_order=order_hor+2, infinity_expansion_order=order_inf+3, method=method)
    λ = X.mode.lambda
    Binc = X.incidence_amplitude
    Bref = X.reflection_amplitude
    Y0, Yp0 = Y_initial_values_inf_p2(m, a, omega, λ, Binc, Bref, rsout; order=order_inf)
    Y_ingoing_hor_func, Yp_ingoing_hor_func = Y_and_Yp_ingoing_hor_p2(m, a, omega, λ; order=order_hor)
    u0 = SA[dtype(Y0); dtype(Yp0)]
    rsspan = (rsout, rsin)  # Integrate inward
    p = (m=m, a=a, omega=omega, X=X)
    odeprob = ODEProblem(GSN_source_eqn_p2, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    # Get boundary values
    Y_in = odesoln(rsin)[1]
    Yp_in = odesoln(rsin)[2] * (rin^2 + a^2) / (rin^2 - 2 * rin+ a^2)

    # Define final callable function
    function Y_inf(rs)
        r = r_from_rstar(a, rs)
        if rs > rsout
            Y_inf_val, Yp_inf_val = Y_initial_values_inf_p2(m, a, omega, λ, Binc, Bref, rs)
            return Y_inf_val, Yp_inf_val * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        elseif rs < rsin
            Y_hor_val = Y_in - (Y_ingoing_hor_func(r) - Y_ingoing_hor_func(rin)) + Yp_in * (r - rin)
            Yp_hor_val = Yp_in + (Yp_ingoing_hor_func(r) - Yp_ingoing_hor_func(rin))
            return Y_hor_val, Yp_hor_val
        else
            sol = odesoln(rs)
            return sol[1], sol[2] * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        end
    end

    Y_value_hor, Yp_value_hor = Y_inf(-200)

    function Y_hor(rs)
        r = r_from_rstar(a, rs)
        rp = 1 + sqrt(1 - a^2)
        Y = Y_inf(rs)[1]- Yp_value_hor * r + Yp_value_hor * rp - Y_value_hor
        Yp = Y_inf(rs)[2] - Yp_value_hor
        return Y, Yp
    end

    return YSolutionResult(
        "IN",
        (s=2, l=l, m=m, a=a, omega=omega, lambda=λ),
        (rin=rin, rout=rout, Binc=Binc, Bref=Bref),
        (X=rs -> X(rs), Y_inf=rs -> Y_inf(rs)[1], Yp_inf=rs -> Y_inf(rs)[2], Y_hor=rs -> Y_hor(rs)[1], Yp_hor=rs -> Y_hor(rs)[2])
    )
end

function Y_solution_up_m2(l::Int, m::Int, a, omega, rsin, rsout, order_inf, order_hor, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Precompute GSN radial solution and initial data
    rin = r_from_rstar(a, rsin)
    rout = r_from_rstar(a, rsout)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(-2, l, m, a, omega, UP, rsin, rsout; horizon_expansion_order=order_hor+2, infinity_expansion_order=order_inf+3, method=method)
    λ = X.mode.lambda
    Cinc = X.incidence_amplitude
    Cref = X.reflection_amplitude
    Y0, Yp0 = Y_initial_values_hor_m2(m, a, omega, λ, Cinc, Cref, rsin; order=order_hor)
    Y_outgoing_inf_func, Yp_outgoing_inf_func = Y_and_Yp_outgoing_inf_m2(m, a, omega, λ; order=order_inf)
    u0 = SA[dtype(Y0); dtype(Yp0)]
    rsspan = (rsin, rsout)  # Integrate outward
    p = (m=m, a=a, omega=omega, X=X)
    odeprob = ODEProblem(GSN_source_eqn_m2, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    # Get boundary values
    Y_up = odesoln(rsout)[1] - Y_outgoing_inf_func(rout)
    Yp_up = odesoln(rsout)[2] * (rout^2 + a^2) / (rout^2 - 2 * rout+ a^2) - Yp_outgoing_inf_func(rout)

    # Define final callable function
    function Y_hor(rs)
        r = r_from_rstar(a, rs)
        if rs < rsin
            Y_hor_val, Yp_hor_val = Y_initial_values_hor_m2(m, a, omega, λ, Cinc, Cref, rs; ifrs=true)
            return Y_hor_val, Yp_hor_val
        elseif rs > rsout
            Y_inf_val = Y_up + Y_outgoing_inf_func(r) + Yp_up * (r - rout)
            Yp_inf_val = Yp_up + Yp_outgoing_inf_func(r)
            return Y_inf_val, Yp_inf_val
        else
            sol = odesoln(rs)
            return sol[1], sol[2] * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        end
    end

    function Y_inf(rs)
        r = r_from_rstar(a, rs)
        if rs <= rsout
            Y = Y_hor(rs)[1] + Yp_up * (rout - r) - Y_up
        else
            Y = Y_outgoing_inf_func(r)
        end
        Yp = Y_hor(rs)[2] - Yp_up
        return Y, Yp
    end

    return YSolutionResult(
        "UP",
        (s=-2, l=l, m=m, a=a, omega=omega, lambda=λ),
        (rin=rin, rout=rout, Cinc=Cinc, Cref=Cref),
        (X=rs -> X(rs), Y_inf=rs -> Y_inf(rs)[1], Yp_inf=rs -> Y_inf(rs)[2], Y_hor=rs -> Y_hor(rs)[1], Yp_hor=rs -> Y_hor(rs)[2])
    )
end

function Y_solution_up_p2(l::Int, m::Int, a, omega, rsin, rsout, order_inf, order_hor, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Precompute GSN radial solution and initial data
    rin = r_from_rstar(a, rsin)
    rout = r_from_rstar(a, rsout)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(2, l, m, a, omega, UP, rsin, rsout; horizon_expansion_order=order_hor+2, infinity_expansion_order=order_inf+3, method=method)
    λ = X.mode.lambda
    Cinc = X.incidence_amplitude
    Cref = X.reflection_amplitude
    Y0, Yp0 = Y_initial_values_hor_p2(m, a, omega, λ, Cinc, Cref, rsin; order=order_hor)
    Y_outgoing_inf_func, Yp_outgoing_inf_func = Y_and_Yp_outgoing_inf_p2(m, a, omega, λ; order=order_inf)
    u0 = SA[dtype(Y0); dtype(Yp0)]
    rsspan = (rsin, rsout)  # Integrate outward
    p = (m=m, a=a, omega=omega, X=X)
    odeprob = ODEProblem(GSN_source_eqn_p2, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    # Get boundary values
    Y_up = odesoln(rsout)[1] - Y_outgoing_inf_func(rout)
    Yp_up = odesoln(rsout)[2] * (rout^2 + a^2) / (rout^2 - 2 * rout+ a^2) - Yp_outgoing_inf_func(rout)

    # Define final callable function
    function Y_hor(rs)
        r = r_from_rstar(a, rs)
        if rs < rsin
            Y_hor_val, Yp_hor_val = Y_initial_values_hor_p2(m, a, omega, λ, Cinc, Cref, rs; ifrs=true)
            return Y_hor_val, Yp_hor_val
        elseif rs > rsout
            Y_inf_val = Y_up + Y_outgoing_inf_func(r) + Yp_up * (r - rout)
            Yp_inf_val = Yp_up + Yp_outgoing_inf_func(r)
            return Y_inf_val, Yp_inf_val
        else
            sol = odesoln(rs)
            return sol[1], sol[2] * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        end
    end

    function Y_inf(rs)
        r = r_from_rstar(a, rs)
        if rs <= rsout
            Y = Y_hor(rs)[1] + Yp_up * (rout - r) - Y_up
        else
            Y = Y_outgoing_inf_func(r)
        end
        Yp = Y_hor(rs)[2] - Yp_up
        return Y, Yp
    end

    return YSolutionResult(
        "UP",
        (s=2, l=l, m=m, a=a, omega=omega, lambda=λ),
        (rin=rin, rout=rout, Cinc=Cinc, Cref=Cref),
        (X=rs -> X(rs), Y_inf=rs -> Y_inf(rs)[1], Yp_inf=rs -> Y_inf(rs)[2], Y_hor=rs -> Y_hor(rs)[1], Yp_hor=rs -> Y_hor(rs)[2])
    )
end

# Wrapper for external use
function Y_solution(s::Int, l::Int, m::Int, a, omega, bc::BoundaryCondition; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order_inf = 6, order_hor = 3)
    if s == 2
        if bc == IN
            return Y_solution_in_p2(l, m, a, omega, rsin, rsout, order_inf, order_hor)
        elseif bc == UP
            return Y_solution_up_p2(l, m, a, omega, rsin, rsout, order_inf, order_hor)
        else
            throw(DomainError(bc, "Boundary condition must be IN or UP"))
        end
    elseif s == -2
        if bc == IN
            return Y_solution_in_m2(l, m, a, omega, rsin, rsout, order_inf, order_hor)
        elseif bc == UP
            return Y_solution_up_m2(l, m, a, omega, rsin, rsout, order_inf, order_hor)
        else
            throw(DomainError(bc, "Boundary condition must be IN or UP"))
        end
    else 
        throw(DomainError(s, "Only s = +2/-2 is supported"))
    end
end

function Ypp_m2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=6)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(-2, l, m, a, omega, IN, rsin, rsout; method=method)
    Binc = X.incidence_amplitude
    Bref = X.reflection_amplitude

    function Ypp_num(rs)
        r = r_from_rstar(a, rs)
        rm = 1 - sqrt(1 - a^2)
        rp = 1 + sqrt(1 - a^2)
        phase = exp(-1im * (omega * rs - a * m * log((r - rm)/(r - rp))/(rm - rp)))
        return X(rs) * phase / (r^2 * sqrt(r^2 + a^2))
    end

    function Ypp_ana(rs; order=order)
        r = r_from_rstar(a, rs)
        if order < 0
            throw(DomainError(order, "Only non-negative expansion order is supported"))
        end
        Ypp_out = 0.0+0.0im
        Ypp_in = 0.0+0.0im
        for n in 0:order
            Ypp_out += Ypp_outgoing_inf_m2(m, a, omega, X.mode.lambda, n) / r^n
            Ypp_in += Ypp_ingoing_inf_m2(m, a, omega, X.mode.lambda, n) / r^n
        end
        Ypp_asym = Bref * Ypp_out / r^3 + Binc * Ypp_in * exp(- 2im * omega * r) / r^(3+4im*omega)
        return Ypp_asym
    end

    return (Ypp_num=Ypp_num, Ypp_ana=Ypp_ana)
    
end

function Ypp_p2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=6)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(2, l, m, a, omega, IN, rsin, rsout; method=method)

    function Ypp_num(rs)
        r = r_from_rstar(a, rs)
        rm = 1 - sqrt(1 - a^2)
        rp = 1 + sqrt(1 - a^2)
        phase = exp(1im * (omega * rs - a * m * log((r - rm)/(r - rp))/(rm - rp)))
        return X(rs) * phase / (r^2 * sqrt(r^2 + a^2))
    end

    function Ypp_ana(rs; order=order)
        r = r_from_rstar(a, rs)
        if order < 0
            throw(DomainError(order, "Only non-negative expansion order is supported"))
        end
        Ypp_out = 0.0+0.0im
        for n in 0:order
            Ypp_out += Ypp_outgoing_inf_p2(m, a, omega, X.mode.lambda, n) / r^n
        end
        Ypp_asym = Ypp_out * exp(2im * omega * r) / r^(3-4im*omega)
        return Ypp_asym
    end

    return (Ypp_num=Ypp_num, Ypp_ana=Ypp_ana)
    
end

end