module Transformation

using ForwardDiff
using ..Kerr
using ..Potentials

export alpha, alpha_prime, beta, beta_prime, eta, eta_prime

const I = 1im # Mathematica being Mathematica

function alpha(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return
    elseif s == +1
        return
    elseif s == -1
        return
    elseif s == +2
        return
    elseif s == -2
        return begin
            (1/(r^2*(a^2 + (-2 + r)*r)))*(4*a^3*m*r*(-I + r*omega) + 2*a*m*r^2*(3*I - I*r + 2*r^2*omega) + 
            a^4*(6 + 4*I*r*omega - 2*r^2*omega^2) + a^2*r*(-24 + r*(12 - 2*m^2 + lambda - 6*I*omega) + 12*I*r^2*omega - 4*r^3*omega^2) + 
            r^2*(24 - 2*r*(12 + lambda) + r^2*(6 + lambda - 18*I*omega) + 8*I*r^3*omega - 2*r^4*omega^2))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function beta(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return
    elseif s == +1
        return
    elseif s == -1
        return
    elseif s == +2
        return
    elseif s == -2
        return begin
            (2*(a^2 + (-2 + r)*r)*(I*a*m*r + a^2*(-2 - I*r*omega) + r*(3 - r - I*r^2*omega)))/r
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function alpha_prime(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return
    elseif s == +1
        return
    elseif s == -1
        return
    elseif s == +2
        return
    elseif s == -2
        return begin
            -((1/(r^3*(a^2 + (-2 + r)*r)^2))*(2*(-2*I*a^5*m*r + a^6*(6 + 2*I*r*omega) + 
            2*a^4*r*(-15 + r*(6 - 4*I*omega) + r^2*omega^2 + r^3*omega^2) + 
            2*a^2*r^2*(24 + r*(-18 + m^2 + 3*I*omega) + r^2*(3 - m^2 + 6*I*omega) + 2*r^4*omega^2 - r^3*omega*(3*I + 2*omega)) + 
            2*r^3*(-12 + 12*r + r^2*(-3 - 9*I*omega) + 8*I*r^3*omega + r^5*omega^2 - r^4*omega*(2*I + 3*omega)) + 
            a*m*r^3*(-6*I + 6*I*r + r^2*(-I + 4*omega)) + a^3*m*r^2*(8*I - r*(5*I + 4*omega)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function beta_prime(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return
    elseif s == +1
        return
    elseif s == -1
        return
    elseif s == +2
        return
    elseif s == -2
        return begin
            2*(-6 + 2*I*a*m*(-1 + r) + (2*a^4)/r^2 + 10*r + r^2*(-3 + 6*I*omega) - 4*I*r^3*omega + 
            a^2*(-3 - 2*I*(-1 + 2*r)*omega))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function eta_coefficient(s::Int, m::Int, a, omega, lambda, order)
    if s == 0
        return
    elseif s == +1
        return
    elseif s == -1
        return
    elseif s == +2
        return 
    elseif s == -2
        if order == 0
            return 2*lambda + lambda^2 - 12*I*omega + 12*a*m*omega - 12*a^2*omega^2
        elseif order == -1
            return 8*I*a*m*lambda + 24*I*a^2*omega - 8*I*a^2*lambda*omega
        elseif order == -2
            return 12*a^2 + 24*I*a*m - 24*a^2*m^2 - 24*I*a^2*omega + 48*a^3*m*omega - 24*a^4*omega^2
        elseif order == -3
            return -24*a^2 - 24*I*a^3*m + 24*I*a^4*omega
        elseif order == -4
            return 12*a^4
        else
            return 0
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function eta(s::Int, m::Int, a, omega, lambda, r)
    # eta(r) = c0 + c1/r + c2/r^2 + c3/r^3 + c4/r^4
    c0 = eta_coefficient(s, m, a, omega, lambda, 0)
    c1 = eta_coefficient(s, m, a, omega, lambda, -1)
    c2 = eta_coefficient(s, m, a, omega, lambda, -2)
    c3 = eta_coefficient(s, m, a, omega, lambda, -3)
    c4 = eta_coefficient(s, m, a, omega, lambda, -4)
    return c0 + c1/r + c2/r^2 + c3/r^3 + c4/r^4
end

function eta_prime(s::Int, m::Int, a, omega, lambda, r)
    # FIXME it might be a bit quicker to derive this by hand
    return ForwardDiff.derivative(x -> eta(s, m, a, omega, lambda, x), r)
end

end