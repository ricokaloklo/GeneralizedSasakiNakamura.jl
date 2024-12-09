module Transformation

using ForwardDiff
using ..Kerr

export alpha, alpha_prime, beta, beta_prime
export eta_coefficient, eta, eta_prime, eta_primeprime

const I = 1im # Mathematica being Mathematica

function alpha(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return 1
    elseif s == +1
        return begin
            (sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*((-I)*a^3*m - I*a*m*r^2 + I*a^4*omega + 
            r^3*(1 + I*r*omega) + a^2*(-2 + r + 2*I*r^2*omega)))/(r^2*sqrt(a^2 + r^2))
        end
    elseif s == -1
        return begin
            (sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*(-r - (I*(a^2 + r^2)*((-a)*m + (a^2 + r^2)*omega))/
            (a^2 + (-2 + r)*r)))/(r^2*sqrt(a^2 + r^2))
        end
    elseif s == +2
        return begin
            (1/(r^2*(a^2 + (-2 + r)*r)))*(4*a^3*m*r*(I + r*omega) + 2*a*m*r^2*(I - 3*I*r + 2*r^2*omega) - 
            2*a^4*(-3 + 2*I*r*omega + r^2*omega^2) + r^3*(-2*lambda + r*(2 + lambda + 10*I*omega) - 2*r^3*omega^2) - 
            a^2*r*(8 + 2*m^2*r - r*lambda + 2*I*r*omega + 4*I*r^2*omega + 4*r^3*omega^2))
        end
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
        return 0
    elseif s == +1
        return begin
            ((a^2 + r^2)/(a^2 + (-2 + r)*r))^(3/2)/(r^2*sqrt(a^2 + r^2))
        end
    elseif s == -1
        return begin
            (sqrt(a^2 + r^2)*sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2)))/r^2
        end
    elseif s == +2
        return begin
            (-2*I*a*m*r + a^2*(-4 + 2*I*r*omega) + 2*r*(3 - r + I*r^2*omega))/(r*(a^2 + (-2 + r)*r)^3)
        end
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
        return 0
    elseif s == +1
        return begin
            -((I*sqrt(a^2 + r^2)*(-2*a^5*m + a^3*m*(5 - 3*r)*r - a*m*(-1 + r)*r^3 + 2*a^6*omega - 
            r^4*(I - 3*r*omega + r^2*omega) - a^2*r*(10*I - 9*I*r + r^2*(I + 2*omega)) + 
            a^4*(4*I + 3*r^2*omega - r*(I + 5*omega))))/(r^3*(a^2 + (-2 + r)*r)^2*
            sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))))
        end
    elseif s == -1
        return begin
            (I*(a^2 + r^2)^(3/2)*(-2*a^5*m + a^3*m*(5 - 3*r)*r - a*m*(-1 + r)*r^3 + 2*a^6*omega + 
            r^3*(-2*I + I*r + 3*r^2*omega - r^3*omega) + a^4*r*(-I + (-5 + 3*r)*omega) + 
            a^2*r^2*(3*I - r*(I + 2*omega))))/(r^3*((a^2 + (-2 + r)*r)*(a^2 + r^2))^(3/2))
        end
    elseif s == +2
        return begin
            -((1/(r^3*(a^2 + (-2 + r)*r)^2))*(2*(2*I*a^5*m*r + a^3*m*r^2*(-8*I + r*(9*I - 4*omega)) + 
            a^6*(6 - 2*I*r*omega) + 2*r^5*(1 + 5*I*omega + (-3 + r)*r^2*omega^2) + 
            2*a^4*r*(-11 + r*(6 + 4*I*omega) + r^3*omega^2 + r^2*omega*(-2*I + omega)) + 
            2*a^2*r^2*(8 + r*(-6 + m^2 + I*omega) - r^2*(1 + m^2 + 6*I*omega) + 2*r^4*omega^2 - 
            r^3*omega*(I + 2*omega)) + a*m*r^3*(-2*I + 2*I*r + r^2*(-3*I + 4*omega)))))
        end
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
        return 0
    elseif s == +1
        return begin
            -((sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*(2*a^4 + 3*(-1 + r)*r^3 + a^2*r*(-7 + 5*r)))/
            (r^3*(a^2 + (-2 + r)*r)^2*sqrt(a^2 + r^2)))
        end
    elseif s == -1
        return begin
            (-2*a^6 - 3*a^4*(-1 + r)*r + 2*a^2*r^3 + (-1 + r)*r^5)/
            (r^3*sqrt(a^2 + r^2)*sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2)))
        end
    elseif s == +2
        return begin
            (1/(r^2*(a^2 + (-2 + r)*r)^4))*(2*(2*a^4 + 6*I*a*m*(-1 + r)*r^2 + 
            a^2*r*(-16 + r*(13 + 6*I*omega) - 4*I*r^2*omega) + 
            r^2*(18 - 22*r + r^2*(5 + 2*I*omega) - 4*I*r^3*omega)))
        end
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
        if order == 0
            return 1
        elseif order == -1
            return 0
        elseif order == -2
            return 0
        elseif order == -3
            return 0
        elseif order == -4
            return 0
        else
            return 0
        end
    elseif s == +1
        if order == 0
            return -2 - lambda
        elseif order == -1
            return 2*I*a*m
        elseif order == -2
            return -3*a^2 - 2*a^2*lambda
        elseif order == -3
            return -2*a^2 + 2*I*a^3*m
        elseif order == -4
            return -a^4 - a^4*lambda
        else
            return 0
        end
    elseif s == -1
        if order == 0
            return -lambda
        elseif order == -1
            return -2*I*a*m
        elseif order == -2
            return a^2 - 2*a^2*lambda
        elseif order == -3
            return -2*a^2 - 2*I*a^3*m
        elseif order == -4
            return a^4 - a^4*lambda
        else
            return 0
        end
    elseif s == +2
        if order == 0
            return 24 + 10*lambda + lambda^2 + 12*I*omega + 12*a*m*omega - 12*a^2*omega^2
        elseif order == -1
            return -32*I*a*m - 8*I*a*m*lambda + 8*I*a^2*omega + 8*I*a^2*lambda*omega
        elseif order == -2
            return 12*a^2 - 24*I*a*m - 24*a^2*m^2 + 24*I*a^2*omega + 48*a^3*m*omega - 24*a^4*omega^2
        elseif order == -3
            return -24*a^2 + 24*I*a^3*m - 24*I*a^4*omega
        elseif order == -4
            return 12*a^4
        else
            return 0
        end
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
    # eta'(r) = -c1/r^2 - 2c2/r^3 - 3c3/r^4 - 4c4/r^5
    c1 = eta_coefficient(s, m, a, omega, lambda, -1)
    c2 = eta_coefficient(s, m, a, omega, lambda, -2)
    c3 = eta_coefficient(s, m, a, omega, lambda, -3)
    c4 = eta_coefficient(s, m, a, omega, lambda, -4)
    return -c1/r^2 - 2c2/r^3 - 3c3/r^4 - 4c4/r^5
end

function eta_primeprime(s::Int, m::Int, a, omega, lambda, r)
    # eta''(r) = 2c1/r^3 + 6c2/r^4 + 12c3/r^5 + 20c4/r^6
    c1 = eta_coefficient(s, m, a, omega, lambda, -1)
    c2 = eta_coefficient(s, m, a, omega, lambda, -2)
    c3 = eta_coefficient(s, m, a, omega, lambda, -3)
    c4 = eta_coefficient(s, m, a, omega, lambda, -4)
    return 2*c1/r^3 + 6*c2/r^4 + 12*c3/r^5 + 20*c4/r^6
end

end