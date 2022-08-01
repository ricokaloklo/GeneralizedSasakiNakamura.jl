using Roots
using Interpolations

include("kerr.jl")

function rstar_from_rp(a, r_from_rp)
    rp = r_plus(a)
    rm = r_minus(a)
    return rp + r_from_rp + (2*1*rp)/(rp-rm) * log(r_from_rp/(2*1)) - (2*1*rm)/(rp-rm) * log((r_from_rp+rp-rm)/(2*1))
end

function rstar_from_r(a, r)
    rp = r_plus(a)
    return rstar_from_rp(a, r-rp)
end

function r_from_rstar(a, rstar)
    rp = r_plus(a)
    #=
    To find r' that solves the equation rstar_from_r(r') = rstar,
    we first write r' = rp + h' and instead solve for h', then add back rp
    i.e. we solve the equation rstar_from_rp(h') = rstar
    =#
    f(x) = rstar_from_rp(a, x) - rstar
    if rstar <= 0
        #=
        When rstar <= 0, it is more efficient to use bisection method,
        this is because in this case h' is bounded (weakly),
        h' cannot be smaller than 0 (since r=rp maps to rstar=-Inf),
        and suppose h'_u solves the equation
        rstar_from_rp(h'_u) = 0.0, in which h'_u is a function of |a|
        The maximum of h'_u occurs when |a| -> 1 with value ~ 1.328
        =#
        return rp + find_zero(f, (0, 1.4))
    else
        # Use secant method instead; for large rstar, rstar \approx r
        return rp + find_zero(f, rstar)
    end
end

function build_r_from_rstar_interpolant(a, rsin, rsout, rsstep::Float64=0.01)
    rsgrid = collect(rsin:rsstep:rsout)
    return linear_interpolation(rsgrid, (x -> r_from_rstar(a, x)).(rsgrid))
end