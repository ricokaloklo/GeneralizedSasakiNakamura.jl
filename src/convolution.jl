using DifferentialEquations
using FiniteDifferences

include("kerr.jl")

function Levin!(du, u, p, r)
    # Compute q'(r) using finite differencing (5th order)
    qprime = central_fdm(5, 1)(p.q, r)
    # p'(r) = f(r) - A(r)p(r)
    du[1] = p.f(r) - (1im*qprime)*u[1]
end

function convolution_integral(s::Int, a, R_homo, sourceterm_without_phasing, sourceterm_phasing, lower_limit, upper_limit; reltol=1e-10, abstol=1e-10)
    #=
    Calculate the convolution integral

    I = \int_{a}^{b} R(r) * Delta(r)^s |T(r)|*exp(i * phasing(r)) dr

    using Levin method instead of using a traditioanl quadrature rule. This is because
    the integrand is highly oscillatory and therefore a very very small stepsize has to
    be used in order to obtain a reasonable estimate, otherwise imperfect cancellations
    would accumulate to give (very) wrong estimate of the integral.

    In short, Levin method transform the integration of a rapidly oscillatory integrand
    into solving (possible a system of) ordinary differential equation(s) of non-oscillatory
    functions.

    Following the Levin paper [Fast integration of rapidly oscillatory functions (1994)], we write

    I = \int_{a}^{b} f(r) S(r) dr
    
    where f(r) is a non-oscillatory function and S(r) is the oscillatory part.
    We can set f(r) = R(r) \Delta^{s}(r) * |T(r)| and S(r) \equiv exp[i(\omega t(r) - m phi(r)],
    or S \equiv exp[i q(r)] where q(r) \equiv \omega t(r) - m phi(r). We verified that indeed
    f(r) is only slowly varying with r.

    Since S(r) is of the form of an exponential, we can apply the Levin method that

    dS/dr = exp[i q(r)] iq'(r) = iq'(r) S

    with (m=1 in Levin's paper) a scalar function A(r) \equiv iq'(r). Therefore we only need to solve
    a scalar ODE.

    Now we seek for a function p(r) such that

    I = \int_{a}^{b} f(r) S(r) dr 
      = \int_{a}^{b} d/dr ( p(r) * S(r) ) dr
      = p(b)S(b) - p(a)S(a)
    by the Fundamental Theorem of Calculus

    Notice that d/dr ( p(r) * S(r) ) = p'(r) S(r) + p(r) S'(r) by simple product rule
    However by construnction we have S'(r) = A(r) S(r) and so
    d/dr ( p(r) * S(r) ) = p'(r) S(r) + A(r) * S(r) p(r) = [p'(r) + A(r)*p(r)] * S(r) = f(r) * S(r)

    Effectively we transform the integration into a simple ODE problem
    where we find a solution to p(r) satisfying

    p'(r) + A(r)*p(r) = f(r)

    subject to an arbitrary initial condition at either endpoint of the integration domain (i.e. a or b).
    Usually one sets p(a) = 0 such that the original integral I is simply given by

    I = p(b) * S(b)
=#
    f(r) = R_homo(r) * Delta(a, r)^s * sourceterm_without_phasing(r)
    q(r) = sourceterm_phasing(r)

    u0 = [0.0+0.0im] # Complex output
    rspan = (lower_limit, upper_limit)
    p = (f=f, q=q)
    odeprob = ODEProblem(Levin!, u0, rspan, p)
    odealgo = RK4() # Seems to be fine
    odesoln = solve(odeprob, odealgo; reltol=reltol, abstol=abstol)

    return first(odesoln(last(rspan))) * exp(1im*sourceterm_phasing(last(rspan)))
end