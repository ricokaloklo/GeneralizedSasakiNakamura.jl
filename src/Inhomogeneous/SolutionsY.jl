module SolutionsY

using DifferentialEquations
using StaticArrays

using ..GeneralizedSasakiNakamura: BoundaryCondition, IN, UP
using ..GeneralizedSasakiNakamura: GSN_radial
using ..Coordinates: r_from_rstar
using ..AsymptoticExpansionCoefficientsY
using ..InitialConditionsY

export Y_solution, Ypp_asymptotic_inf_error, Ypp_asymptotic_hor_error

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

    X = GSN_radial(-2, l, m, a, omega, IN, rsin, rsout; horizon_expansion_order=order_hor+2, infinity_expansion_order=order_inf+3, method="linear")
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
    Y0_in = Y_ingoing_hor_func(rin)
    Yp0_in = Yp_ingoing_hor_func(rin)
    y1 = Yp_in - Yp0_in
    y0 = Y_in - Y0_in - (Yp_in - Yp0_in) * rin

    # Define final callable function
    function Y_inf(rs)
        r = r_from_rstar(a, rs)
        if rs > rsout
            Y_inf_func, Yp_inf_func = Y_initial_values_inf_m2(m, a, omega, λ, Binc, Bref, rs)
            return Y_inf_func, Yp_inf_func * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        elseif rs < rsin
            Y_hor_func = Y_ingoing_hor_func(r) + y1 * r + y0
            Yp_hor_func = Yp_ingoing_hor_func(r) + y1
            return Y_hor_func, Yp_hor_func
        else
            sol = odesoln(rs)
            return sol[1], sol[2] * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        end
    end

    function Y_hor(rs)
        r = r_from_rstar(a, rs)
        Y = Y_inf(rs)[1] - y1 * r - y0
        Yp = Y_inf(rs)[2] - y1
        return Y, Yp
    end

    return YSolutionResult(
        "IN",
        (s=-2, l=l, m=m, a=a, omega=omega, lambda=λ),
        (rin=rin, rout=rout, Binc=Binc, Bref=Bref, y0=y0, y1=y1),
        (X=rs -> X(rs), Y_inf=rs -> Y_inf(rs)[1], Yp_inf=rs -> Y_inf(rs)[2], Y_hor=rs -> Y_hor(rs)[1], Yp_hor=rs -> Y_hor(rs)[2])
    )
end

function Y_solution_in_p2(l::Int, m::Int, a, omega, rsin, rsout, order_inf, order_hor, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Precompute GSN radial solution and initial data
    rin = r_from_rstar(a, rsin)
    rout = r_from_rstar(a, rsout)

    X = GSN_radial(2, l, m, a, omega, IN, rsin, rsout; horizon_expansion_order=order_hor+2, infinity_expansion_order=order_inf+3, method="linear")
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
    Y0_in = Y_ingoing_hor_func(rin)
    Yp0_in = Yp_ingoing_hor_func(rin)
    y1 = Yp_in - Yp0_in
    y0 = Y_in - Y0_in - (Yp_in - Yp0_in) * rin

    # Define final callable function
    function Y_inf(rs)
        r = r_from_rstar(a, rs)
        if rs > rsout
            Y_inf_func, Yp_inf_func = Y_initial_values_inf_p2(m, a, omega, λ, Binc, Bref, rs)
            return Y_inf_func, Yp_inf_func * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        elseif rs < rsin
            Y_hor_func = Y_ingoing_hor_func(r) + y1 * r + y0
            Yp_hor_func = Yp_ingoing_hor_func(r) + y1
            return Y_hor_func, Yp_hor_func
        else
            sol = odesoln(rs)
            return sol[1], sol[2] * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        end
    end

    function Y_hor(rs)
        r = r_from_rstar(a, rs)
        Y = Y_inf(rs)[1] - y1 * r - y0
        Yp = Y_inf(rs)[2] - y1
        return Y, Yp
    end

    return YSolutionResult(
        "IN",
        (s=2, l=l, m=m, a=a, omega=omega, lambda=λ),
        (rin=rin, rout=rout, Binc=Binc, Bref=Bref, y0=y0, y1=y1),
        (X=rs -> X(rs), Y_inf=rs -> Y_inf(rs)[1], Yp_inf=rs -> Y_inf(rs)[2], Y_hor=rs -> Y_hor(rs)[1], Yp_hor=rs -> Y_hor(rs)[2])
    )
end

function Y_solution_up_m2(l::Int, m::Int, a, omega, rsin, rsout, order_inf, order_hor, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Precompute GSN radial solution and initial data
    rin = r_from_rstar(a, rsin)
    rout = r_from_rstar(a, rsout)

    X = GSN_radial(-2, l, m, a, omega, UP, rsin, rsout; horizon_expansion_order=order_hor+2, infinity_expansion_order=order_inf+3, method="linear")
    λ = X.mode.lambda
    Cinc = X.incidence_amplitude
    Cref = X.reflection_amplitude
    Y0, Yp0 = Y_initial_values_inf_up_m2(m, a, omega, λ, rsout; order=order_inf)
    Y_outgoing_hor_func, Yp_outgoing_hor_func = Y_and_Yp_outgoing_hor_m2(m, a, omega, λ; order=order_hor)
    Y_ingoing_hor_func, Yp_ingoing_hor_func = Y_and_Yp_ingoing_hor_m2(m, a, omega, λ; order=order_hor)
    u0 = SA[dtype(Y0); dtype(Yp0)]
    rsspan = (rsout, rsin)  # Integrate outward
    p = (m=m, a=a, omega=omega, X=X)
    odeprob = ODEProblem(GSN_source_eqn_m2, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    # Get boundary values
    Y_up = odesoln(rsin)[1]
    Yp_up = odesoln(rsin)[2] * (rin^2 + a^2) / (rin^2 - 2 * rin + a^2)
    Y0_up = Cinc * Y_outgoing_hor_func(rin) + Cref * Y_ingoing_hor_func(rin)
    Yp0_up = Cinc * Yp_outgoing_hor_func(rin) + Cref * Yp_ingoing_hor_func(rin)
    y1 = Yp_up - Yp0_up
    y0 = Y_up - Y0_up - (Yp_up - Yp0_up) * rin

    # Define final callable function
    function Y_inf(rs)
        r = r_from_rstar(a, rs)
        if rs > rsout
            Y_inf_func, Yp_inf_func = Y_initial_values_inf_up_m2(m, a, omega, λ, rs)
            return Y_inf_func, Yp_inf_func * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        elseif rs < rsin
            Y_hor_func = Cinc * Y_outgoing_hor_func(rin) + Cref * Y_ingoing_hor_func(rin) + y1 * r + y0
            Yp_hor_func = Cinc * Yp_outgoing_hor_func(rin) + Cref * Yp_ingoing_hor_func(rin) + y1
            return Y_hor_func, Yp_hor_func
        else
            sol = odesoln(rs)
            return sol[1], sol[2] * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        end
    end

    function Y_hor(rs)
        r = r_from_rstar(a, rs)
        Y = Y_inf(rs)[1] - y1 * r - y0
        Yp = Y_inf(rs)[2] - y1
        return Y, Yp
    end

    return YSolutionResult(
        "UP",
        (s=-2, l=l, m=m, a=a, omega=omega, lambda=λ),
        (rin=rin, rout=rout, Cinc=Cinc, Cref=Cref, y0=y0, y1=y1),
        (X=rs -> X(rs), Y_inf=rs -> Y_inf(rs)[1], Yp_inf=rs -> Y_inf(rs)[2], Y_hor=rs -> Y_hor(rs)[1], Yp_hor=rs -> Y_hor(rs)[2])
    )
end

function Y_solution_up_p2(l::Int, m::Int, a, omega, rsin, rsout, order_inf, order_hor, dtype=_DEFAULTDATATYPE, odealgo=_DEFAULTSOLVER, reltol=_DEFAULTTOLERANCE, abstol=_DEFAULTTOLERANCE)
    # Precompute GSN radial solution and initial data
    rin = r_from_rstar(a, rsin)
    rout = r_from_rstar(a, rsout)

    X = GSN_radial(2, l, m, a, omega, UP, rsin, rsout; horizon_expansion_order=order_hor+2, infinity_expansion_order=order_inf+3, method="linear")
    λ = X.mode.lambda
    Cinc = X.incidence_amplitude
    Cref = X.reflection_amplitude
    Y0, Yp0 = Y_initial_values_inf_up_p2(m, a, omega, λ, rsout; order=order_inf)
    Y_outgoing_hor_func, Yp_outgoing_hor_func = Y_and_Yp_outgoing_hor_p2(m, a, omega, λ; order=order_hor)
    Y_ingoing_hor_func, Yp_ingoing_hor_func = Y_and_Yp_ingoing_hor_p2(m, a, omega, λ; order=order_hor)
    u0 = SA[dtype(Y0); dtype(Yp0)]
    rsspan = (rsout, rsin)  # Integrate outward
    p = (m=m, a=a, omega=omega, X=X)
    odeprob = ODEProblem(GSN_source_eqn_p2, u0, rsspan, p)
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    # Get boundary values
    Y_up = odesoln(rsin)[1]
    Yp_up = odesoln(rsin)[2] * (rin^2 + a^2) / (rin^2 - 2 * rin + a^2)
    Y0_up = Cinc * Y_outgoing_hor_func(rin) + Cref * Y_ingoing_hor_func(rin)
    Yp0_up = Cinc * Yp_outgoing_hor_func(rin) + Cref * Yp_ingoing_hor_func(rin)
    y1 = Yp_up - Yp0_up
    y0 = Y_up - Y0_up - (Yp_up - Yp0_up) * rin

    # Define final callable function
    function Y_inf(rs)
        r = r_from_rstar(a, rs)
        if rs > rsout
            Y_inf_func, Yp_inf_func = Y_initial_values_inf_up_m2(m, a, omega, λ, rs)
            return Y_inf_func, Yp_inf_func * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        elseif rs < rsin
            Y_hor_func = Cinc * Y_outgoing_hor_func(rin) + Cref * Y_ingoing_hor_func(rin) + y1 * r + y0
            Yp_hor_func = Cinc * Yp_outgoing_hor_func(rin) + Cref * Yp_ingoing_hor_func(rin) + y1
            return Y_hor_func, Yp_hor_func
        else
            sol = odesoln(rs)
            return sol[1], sol[2] * (r^2 + a^2) / (r^2 - 2 * r + a^2)
        end
    end

    function Y_hor(rs)
        r = r_from_rstar(a, rs)
        Y = Y_inf(rs)[1] - y1 * r - y0
        Yp = Y_inf(rs)[2] - y1
        return Y, Yp
    end

    return YSolutionResult(
        "UP",
        (s=2, l=l, m=m, a=a, omega=omega, lambda=λ),
        (rin=rin, rout=rout, Cinc=Cinc, Cref=Cref, y0=y0, y1=y1),
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

function Ypp_in_inf_m2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=6)
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

    function relative_error(rs)
        return abs((Ypp_num(rs) - Ypp_ana(rs)) / Ypp_num(rs))
    end

    return relative_error
end

function Ypp_in_hor_m2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=3)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(-2, l, m, a, omega, IN, rsin, rsout; method=method, horizon_expansion_order=5)
    lambda = X.mode.lambda

    function Ypp_num(rs)
        r = r_from_rstar(a, rs)
        rm = 1 - sqrt(1 - a^2)
        rp = 1 + sqrt(1 - a^2)
        phase = exp(-1im * (omega * rs - a * m * log((r - rm)/(r - rp))/(rm - rp)))
        return X(rs) * phase / (r^2 * sqrt(r^2 + a^2))
    end

    function Ypp_ana(rs; order=order)
        rp = 1 + sqrt(1 - a^2)
        q = (a * rp * m + 2 * a^2 * omega - 4 * rp * omega) / (rp * sqrt(1 - a^2))
        r = r_from_rstar(a, rs)
        if order < 0
            throw(DomainError(order, "Only non-negative expansion order is supported"))
        end
        Ypp_in = 0.0+0.0im
        for n in 0:order
            Ypp_in += Ypp_ingoing_hor_m2(m, a, omega, lambda, n) * (r - rp)^n
        end
        Ypp_asym = Ypp_in * (r - rp)^(im*q)
        return Ypp_asym
    end

    function relative_error(rs)
        return abs((Ypp_num(rs) - Ypp_ana(rs)) / Ypp_num(rs))
    end

    return relative_error
    
end

function Ypp_up_inf_m2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=6)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(-2, l, m, a, omega, UP, rsin, rsout; method=method)

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
        for n in 0:order
            Ypp_out += Ypp_outgoing_inf_m2(m, a, omega, X.mode.lambda, n) / r^n
        end
        Ypp_asym = Ypp_out / r^3
        return Ypp_asym
    end

    function relative_error(rs)
        return abs((Ypp_num(rs) - Ypp_ana(rs)) / Ypp_num(rs))
    end

    return relative_error
end

function Ypp_up_hor_m2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=3)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(-2, l, m, a, omega, UP, rsin, rsout; method=method)
    Cinc = X.incidence_amplitude
    Cref = X.reflection_amplitude

    function Ypp_num(rs)
        r = r_from_rstar(a, rs)
        rm = 1 - sqrt(1 - a^2)
        rp = 1 + sqrt(1 - a^2)
        phase = exp(-1im * (omega * rs - a * m * log((r - rm)/(r - rp))/(rm - rp)))
        return X(rs) * phase / (r^2 * sqrt(r^2 + a^2))
    end

    function Ypp_ana(rs; order=order)
        r = r_from_rstar(a, rs)
        rp = 1 + sqrt(1 - a^2)
        q = (a * rp * m + 2 * a^2 * omega - 4 * rp * omega) / (rp * sqrt(1 - a^2))
        if order < 0
            throw(DomainError(order, "Only non-negative expansion order is supported"))
        end
        Ypp_out = 0.0 + 0.0im
        Ypp_in = 0.0 + 0.0im
        for n in 0:order
            Ypp_out += Ypp_outgoing_hor_m2(m, a, omega, X.mode.lambda, n) * (r - rp)^n
            Ypp_in += Ypp_ingoing_hor_m2(m, a, omega, X.mode.lambda, n) * (r - rp)^n
        end
        Ypp_asym = Cinc * Ypp_out + Cref * Ypp_in * (r - rp)^(im * q)
        return Ypp_asym
    end

    function relative_error(rs)
        return abs((Ypp_num(rs) - Ypp_ana(rs)) / Ypp_num(rs))
    end

    return relative_error
end

function Ypp_in_inf_p2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=6)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(2, l, m, a, omega, IN, rsin, rsout; method=method)
    Binc = X.incidence_amplitude
    Bref = X.reflection_amplitude

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
        Ypp_in = 0.0+0.0im
        for n in 0:order
            Ypp_out += Ypp_outgoing_inf_p2(m, a, omega, X.mode.lambda, n) / r^n
            Ypp_in += Ypp_ingoing_inf_p2(m, a, omega, X.mode.lambda, n) / r^n
        end
        Ypp_asym = Bref * Ypp_out * exp(2im * omega * r) / r^(3-4im*omega) + Binc * Ypp_in / r^3
        return Ypp_asym
    end

    function relative_error(rs)
        return abs((Ypp_num(rs) - Ypp_ana(rs)) / Ypp_num(rs))
    end

    return relative_error
end

function Ypp_in_hor_p2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=3)
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
        rp = 1 + sqrt(1 - a^2)
        if order < 0
            throw(DomainError(order, "Only non-negative expansion order is supported"))
        end
        Ypp_in = 0.0+0.0im
        for n in 0:order
            Ypp_in += Ypp_ingoing_hor_p2(m, a, omega, X.mode.lambda, n) * (r - rp)^n
        end
        Ypp_asym = Ypp_in
        return Ypp_asym
    end

    function relative_error(rs)
        return abs((Ypp_num(rs) - Ypp_ana(rs)) / Ypp_num(rs))
    end

    return relative_error
end

function Ypp_up_inf_p2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=6)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(2, l, m, a, omega, UP, rsin, rsout; method=method)

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

    function relative_error(rs)
        return abs((Ypp_num(rs) - Ypp_ana(rs)) / Ypp_num(rs))
    end

    return relative_error
end

function Ypp_up_hor_p2(l::Int, m::Int, a, omega; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout, order=3)
    if abs(omega) > 1.0
        method = "Riccati"
    else
        method = "linear"
    end
    X = GSN_radial(2, l, m, a, omega, UP, rsin, rsout; method=method)
    Cinc = X.incidence_amplitude
    Cref = X.reflection_amplitude

    function Ypp_num(rs)
        r = r_from_rstar(a, rs)
        rm = 1 - sqrt(1 - a^2)
        rp = 1 + sqrt(1 - a^2)
        phase = exp(1im * (omega * rs - a * m * log((r - rm)/(r - rp))/(rm - rp)))
        return X(rs) * phase / (r^2 * sqrt(r^2 + a^2))
    end

    function Ypp_ana(rs; order=order)
        rp = 1 + sqrt(1 - a^2)
        r = r_from_rstar(a, rs)
        q = (a * rp * m + 2 * a^2 * omega - 4 * rp * omega) / (rp * sqrt(1 - a^2))
        if order < 0
            throw(DomainError(order, "Only non-negative expansion order is supported"))
        end
        Ypp_out = 0.0 + 0.0im
        Ypp_in = 0.0 + 0.0im
        for n in 0:order
            Ypp_out += Ypp_outgoing_hor_p2(m, a, omega, X.mode.lambda, n) * (r - rp)^n
            Ypp_in += Ypp_ingoing_hor_p2(m, a, omega, X.mode.lambda, n) * (r - rp)^n
        end
        Ypp_asym = Cinc * Ypp_out * (r - rp)^(-im * q) + Cref * Ypp_in
        return Ypp_asym
    end

    function relative_error(rs)
        return abs((Ypp_num(rs) - Ypp_ana(rs)) / Ypp_num(rs))
    end

    return relative_error
end

function Ypp_asymptotic_hor_error(s::Int, l::Int, m::Int, a, omega, bc::BoundaryCondition; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout)
    if s == 2
        if bc == 1
            return Ypp_in_hor_p2(l, m, a, omega; rsin = rsin, rsout = rsout, order = 3)
        elseif bc == 2
            return Ypp_up_hor_p2(l, m, a, omega; rsin = rsin, rsout = rsout, order = 3)
        else
            throw(DomainError(bc, "Boundary condition should be IN or UP "))
        end
    elseif s == -2
        if bc == 1
            return Ypp_in_hor_m2(l, m, a, omega; rsin = rsin, rsout = rsout, order = 3)
        elseif bc == 2
            return Ypp_up_hor_m2(l, m, a, omega; rsin = rsin, rsout = rsout, order = 3)
        else
            throw(DomainError(bc, "Boundary condition should be IN or UP"))
        end
    else
        throw(DomainError(s, "Only s = +2/-2 are supported"))
    end
end

function Ypp_asymptotic_inf_error(s::Int, l::Int, m::Int, a, omega, bc; rsin=_DEFAULT_rsin, rsout=_DEFAULT_rsout)
    if s == 2
        if bc == 1
            return Ypp_in_inf_p2(l, m, a, omega; rsin = rsin, rsout = rsout, order = 6)
        elseif bc == 2
            return Ypp_up_inf_p2(l, m, a, omega; rsin = rsin, rsout = rsout, order = 6)
        else
            throw(DomainError(bc, "Boundary condition should be IN or UP "))
        end
    elseif s == -2
        if bc == 1
            return Ypp_in_inf_m2(l, m, a, omega; rsin = rsin, rsout = rsout, order = 6)
        elseif bc == 2
            return Ypp_up_inf_m2(l, m, a, omega; rsin = rsin, rsout = rsout, order = 6)
        else
            throw(DomainError(bc, "Boundary condition should be IN or UP"))
        end
    else
        throw(DomainError(s, "Only s = +2/-2 are supported"))
    end
end

end
