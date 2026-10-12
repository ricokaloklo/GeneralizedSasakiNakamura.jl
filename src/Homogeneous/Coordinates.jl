module Coordinates

using ..Kerr
using Roots
using Interpolations

export rstar_from_r, r_from_rstar

_is_exact_extremal_spin(a) = isreal(a) && isone(abs(real(a)))

# Exact extremal spin (|a| = 1, r_+ = 1): closed-form r* as a function of h = r - r_+, its h
# derivative, and the inversion by Newton's method with a bisection fallback.
function _rstar_from_exact_extremal_rp(h)
    h < 0 && throw(DomainError(h, "r-r_+ must be nonnegative"))
    h == 0 && return -Inf
    return 1 + h - log(4) + 2*(log(h) - 1/h)
end

function _drstar_dh_exact_extremal(h)
    h > 0 || throw(DomainError(h, "r-r_+ must be positive"))
    return 1 + 2/h + 2/(h*h)
end

function _exact_extremal_initial_h(rstar)
    if rstar < 0
        return max(2 / max(-rstar, 1), sqrt(eps(Float64)))
    elseif rstar > 10
        return max(rstar - 2 * log(max(rstar, 2)) + log(4) - 1, sqrt(eps(Float64)))
    end
    return 1.3
end

function _exact_extremal_bracket_h(rstar)
    low = 0.0
    high = rstar > 0 ? max(2.0, rstar + 2.0) : 2.0
    fhigh = _rstar_from_exact_extremal_rp(high) - rstar
    expansions = 0
    while (!isfinite(fhigh) || fhigh < 0) && expansions < 80
        high *= 2
        fhigh = _rstar_from_exact_extremal_rp(high) - rstar
        expansions += 1
    end
    (!isfinite(fhigh) || fhigh < 0) &&
        error("failed to bracket exact extremal rstar=$rstar")

    for _ in 1:260
        mid = 0.5 * (low + high)
        fmid = _rstar_from_exact_extremal_rp(mid) - rstar
        !isfinite(fmid) && (low = mid; continue)
        abs(fmid) <= 1e-12 && return mid
        if fmid < 0
            low = mid
        else
            high = mid
        end
        high - low <= 4eps(max(high, 1.0)) && return 0.5 * (low + high)
    end
    return 0.5 * (low + high)
end

function _exact_extremal_newton_h(rstar, h0)
    h = h0
    for _ in 1:90
        h > 0 || error("exact extremal Newton stepped to nonpositive h")
        f = _rstar_from_exact_extremal_rp(h) - rstar
        abs(f) <= 1e-12 && return h
        d = _drstar_dh_exact_extremal(h)
        isfinite(f) && isfinite(d) && d != 0 ||
            error("nonfinite exact extremal Newton state")
        hnew = h - f / d
        isfinite(hnew) && hnew > 0 ||
            error("exact extremal Newton produced invalid h")
        h = hnew
    end
    error("exact extremal Newton did not converge")
end

function _r_from_rstar_exact_extremal(rstar)
    h = try
        hnewton = _exact_extremal_newton_h(rstar, _exact_extremal_initial_h(rstar))
        residual = abs(_rstar_from_exact_extremal_rp(hnewton) - rstar)
        residual <= 1e-10 ? hnewton : _exact_extremal_bracket_h(rstar)
    catch
        _exact_extremal_bracket_h(rstar)
    end
    return 1 + h
end

function rstar_from_rp(a, r_from_rp)
    if _is_exact_extremal_spin(a)
        return _rstar_from_exact_extremal_rp(r_from_rp)
    elseif abs(a) < 1
        rp = r_plus(a)
        rm = r_minus(a)
        return rp + r_from_rp + (2*1*rp)/(rp-rm) * log(r_from_rp/(2*1)) - (2*1*rm)/(rp-rm) * log((r_from_rp+rp-rm)/(2*1))
    else
        throw(ArgumentError("a must be in the range [-1, 1]"))
    end
end

@doc raw"""
    rstar_from_r(a, r)

Convert a Boyer-Lindquist coordinate `r` to the corresponding tortoise coordinate `rstar`.
"""
function rstar_from_r(a, r)
    _is_exact_extremal_spin(a) && return _rstar_from_exact_extremal_rp(r - 1)

    rp = r_plus(a)
    return rstar_from_rp(a, r-rp)
end

@doc raw"""
    r_from_rstar(a, rstar)

Convert a tortoise coordinate `rstar` to the corresponding Boyer-Lindquist coordinate `r`. 
It uses a bisection method when `rstar <= 0`, and Newton method otherwise.

The function assumes that $r \geq r_{+}$ where $r_{+}$ is the outer event horizon.
"""
function r_from_rstar(a, rstar)
    _is_exact_extremal_spin(a) && return _r_from_rstar_exact_extremal(rstar)

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
