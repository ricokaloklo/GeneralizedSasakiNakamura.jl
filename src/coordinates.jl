using Roots

include("kerr.jl")

function rstar_from_rp(M, a, r_from_rp)
    rp = r_plus(M, a)
    rm = r_minus(M, a)
    return rp + r_from_rp + (2*M*rp)/(rp-rm) * log(r_from_rp/(2*M)) - (2*M*rm)/(rp-rm) * log((r_from_rp+rp-rm)/(2*M))
end

function rstar_from_r(M, a, r)
    rp = r_plus(M, a)
    return rstar_from_rp(M, a, r-rp)
end

function r_from_rstar(M, a, rstar)
    rp = r_plus(M, a)
    #=
    To find r' that solves the equation rstar_from_r(r') = rstar,
    we first write r' = rp + h' and instead solve for h', then add back rp
    i.e. we solve the equation rstar_from_rp(h') = rstar
    =#
    f(x) = rstar_from_rp(M, a, x) - rstar
    if rstar < 0
        # Use bisection method
        return rp + find_zero(f, (0, 2.6))
    else
        # Use secant method instead; for large rstar, rstar \approx r
        return rp + find_zero(f, rstar)
    end
end
