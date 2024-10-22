module ComplexFrequencies

using ..Kerr
using ..Coordinates

using DifferentialEquations

const I = 1im # Mathematica being Mathematica
_DEFAULTDATATYPE = ComplexF64 # Double precision by default
_DEFAULTSOLVER = Vern9()
_DEFAULTTOLERANCE = 1e-12

function drdrho!(du, u, p, rho)
    #=
    u[1] = r
    du[1] = dr/drho = exp(i\beta) \Delta/(r^2 + a^2)

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

    r_from_rho(rho) = odesoln(rho)[1]
    return r_from_rho
end

end