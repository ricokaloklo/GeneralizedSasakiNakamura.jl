module Coordinates

using ..Kerr
using Roots
using Interpolations

export rstar_from_r, r_from_rstar

function rstar_from_rp(a, r_from_rp)
    rp = r_plus(a)
    rm = r_minus(a)

    if abs(a) < 1
        return rp + r_from_rp + (2*1*rp)/(rp-rm) * log(r_from_rp/(2*1)) - (2*1*rm)/(rp-rm) * log((r_from_rp+rp-rm)/(2*1))
    elseif abs(a) == 1
        return rp + r_from_rp - log(4) + 2*(log(r_from_rp) - 1/(r_from_rp))
    else
        throw(ArgumentError("a must be in the range [-1, 1]"))
    end
end

@doc raw"""
    rstar_from_r(a, r)

Convert a Boyer-Lindquist coordinate `r` to the corresponding tortoise coordinate `rstar`.
"""
function rstar_from_r(a, r)
    rp = r_plus(a)
    return rstar_from_rp(a, r-rp)
end

@doc raw"""
    r_from_rstar(a, rstar)

Convert a tortoise coordinate `rstar` to the corresponding Boyer-Lindquist coordiante `r`. 
It uses a bisection method when `rstar <= 0`, and Newton method otherwise.

The function assumes that $r \geq r_{+}$ where $r_{+}$ is the outer event horizon.
"""
function r_from_rstar(a, rstar)
    rp = r_plus(a)
    #=
    To find r' that solves the equation rstar_from_r(r') = rstar,
    we first write r' = rp + h' and instead solve for h', then add back rp
    i.e. we solve the equation rstar_from_rp(h') = rstar
    =#

    # The root-finding algorithm might try a negative x, which is not allowed
    # We rectify this by taking an absolute value, i.e. we solve for distance from rp
    f(x) = rstar_from_rp(a, abs(x)) - rstar
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
        # Use Newton method instead; for large rstar, rstar \approx r
        return rp + abs(find_zero((f, x -> sign(x)*((rp + abs(x))^2 + a^2)/Delta(a, rp+abs(x))), rstar, Roots.Newton()))
    end
end

function build_r_from_rstar_interpolant(a, rsin, rsout; rsstep::Float64=0.01)
    rsgrid = collect(rsin:rsstep:rsout)
    return linear_interpolation(rsgrid, (x -> r_from_rstar(a, x)).(rsgrid))
end

end