module InitialConditions

using ..Kerr
using ..Coordinates
using ..AsymptoticExpansionCoefficients

export Xup_initialconditions, Xin_initialconditions
export fansatz, gansatz, dfansatz_dr, dgansatz_dr

const I = 1im # Mathematica being Mathematica
_DEFAULTDATATYPE = AsymptoticExpansionCoefficients._DEFAULTDATATYPE

# The expansion at infinity is an asymptotic series in 1/(omega r): its terms shrink
# only down to an optimal order of about omega r and grow afterwards, so summing a
# fixed number of terms is wrong whenever that number exceeds the optimal order.
# _truncated_sum adds terms while they keep shrinking and stops at the first growing
# term (from the third term on, to tolerate a small or vanishing low-order
# coefficient), when the term is below round-off relative to the sum, or at `order`,
# which therefore acts as a cap rather than a fixed count.
function _truncated_sum(term, order)
    ans = term(0)
    prev = Inf
    for i in 1:order
        t = term(i)
        at = abs(t)
        if i >= 3 && at > prev
            break
        end
        ans += t
        at <= eps(Float64) * abs(ans) && break
        prev = at
    end
    return ans
end

function fansatz(func, omega, r; order=3)
    # A template function that gives the asymptotic expansion at infinity
    z = omega * r
    return _truncated_sum(i -> func(i) / z^i, order)
end

function dfansatz_dr(func, omega, r; order=3)
    # A template function that gives the derivative of the asymptotic expansion at infinity
    return _truncated_sum(i -> i == 0 ? zero(func(0)) : -i * func(i) / (omega^i * r^(i + 1)), order)
end

function gansatz(func, a, r; order=1)
    # A template function that gives the asymptotic expansion at horizon
    ans = 0.0
    for i in 0:order
        ans += func(i)*(r-r_plus(a))^i
    end
    return ans
end

function dgansatz_dr(func, a, r; order=1)
    # A template function that gives the derivative of the asymptotic expansion at horizon
    ans = 0.0
    for i in 1:order
        ans += i*func(i)*(r-r_plus(a))^(i-1)
    end
    return ans
end

function Xup_initialconditions(s::Int, m::Int, a, omega, lambda, rsout; order::Int=-1, dtype=_DEFAULTDATATYPE)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    _default_order = 3
    order = (order == -1 ? _default_order : order)

    outgoing_coeff_func = coefficient_sequence(outgoing_coefficient_at_inf, s, m, a, omega, lambda; data_type=dtype)
    fout(r) = fansatz(outgoing_coeff_func, omega, r; order=order)
    dfout_dr(r) = dfansatz_dr(outgoing_coeff_func, omega, r; order=order)
    rout = r_from_rstar(a, rsout)

    _fansatz = fout(rout)
    _dfansatz_dr = dfout_dr(rout)
    phase = exp(1im * omega * rsout)

    return phase*_fansatz, phase*(1im*omega*_fansatz + (Delta(a, rout)/(rout^2 + a^2))*_dfansatz_dr)
end

function Xin_initialconditions(s::Int, m::Int, a, omega, lambda, rsin; order::Int=-1, dtype=_DEFAULTDATATYPE)
    #=
    For Xin, which we obtain by integrating from r_* -> -inf (or r -> r_+),

    Write Xin = \sum_j C^{H}_{-} (r - r_+)^j

    =#
    _default_order = 0
    order = (order == -1 ? _default_order : order)

    ingoing_coeff_func = coefficient_sequence(ingoing_coefficient_at_hor, s, m, a, omega, lambda; data_type=dtype)
    gin(r) = gansatz(ingoing_coeff_func, a, r; order=order)
    dgin_dr(r) = dgansatz_dr(ingoing_coeff_func, a, r; order=order)
    rin = r_from_rstar(a, rsin)

    _gansatz = gin(rin)
    _dgansatz_dr = dgin_dr(rin)
    p = omega - m*omega_horizon(a)
    phase = exp(-1im * p * rsin)

    return phase*_gansatz, phase*(-1im*p*_gansatz + (Delta(a, rin)/(rin^2 + a^2))*_dgansatz_dr)
end

end
