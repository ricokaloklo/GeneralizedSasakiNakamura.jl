global const I = 1im # Mathematica being Mathematica

function alpha(s::Int, m::Int, a, omega, lambda, r)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    if s == 0
        # s = 0
        return 1.0
    elseif s == 1
        # s = +1
        return begin
            (sqrt(1/(a^2 + (-2 + r)*r))*(-2*a^2 - I*a^3*m + a^2*r - I*a*m*r^2 + r^3 + I*(a^2 + r^2)^2*omega))/(a^2 + r^2)
        end
    elseif s == -1
        # s = -1
        return begin
            (I*a*m + r*(-1 + (2*r)/(a^2 + r^2)) - I*(a^2 + r^2)*omega)/sqrt(a^2 + (-2 + r)*r)
        end
    elseif s == 2
        # s = +2
        return begin
            (4*a^7*m*omega - 2*a^8*omega^2 + 12*a^3*m*r^2*(I - I*r + r^2*omega) + 
            2*a*m*r^4*(I - 3*I*r + 2*r^2*omega) + 2*a^5*m*(5*I - 3*I*r + 6*r^2*omega) + 
            a^6*(2 - 2*m^2 + lambda - 10*I*omega - 8*r^2*omega^2) + 
            r^5*(-2*lambda + r*(2 + lambda + 10*I*omega) - 2*r^3*omega^2) + 
            a^4*(8 - 2*r*lambda + r^2*(6 - 4*m^2 + 3*lambda - 10*I*omega) - 12*r^4*omega^2) + 
            a^2*r^2*(-24 - 4*r*lambda + r^2*(6 - 2*m^2 + 3*lambda + 10*I*omega) - 8*r^4*omega^2))/
            ((a^2 + (-2 + r)*r)*(a^2 + r^2)^2)
        end
    elseif s == -2
        # s = -2
        return begin
            (1/(r^2*(a^2 + (-2 + r)*r)))*(4*a^3*m*r*(-I + r*omega) + 2*a*m*r^2*(3*I - I*r + 2*r^2*omega) + a^4*(6 + 4*I*r*omega - 2*r^2*omega^2) + 
            a^2*r*(-24 + r*(12 - 2*m^2 + lambda - 6*I*omega) + 12*I*r^2*omega - 4*r^3*omega^2) + 
            r^2*(24 - 2*r*(12 + lambda) + r^2*(6 + lambda - 18*I*omega) + 8*I*r^3*omega - 2*r^4*omega^2))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function beta(s::Int, m::Int, a, omega, lambda, r)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    if s == 0
        # s = 0
        return 0.0
    elseif s == 1
        # s = +1
        return (1/(a^2 + (-2 + r)*r))^(3/2)
    elseif s == -1
        # s = -1
        return sqrt(a^2 + (-2 + r)*r)
    elseif s == 2
        # s = +2
        return begin
            (2*I*((-a^3)*m - a*m*r^2 + a^4*omega + r^2*(-3*I + I*r + r^2*omega) + 
            a^2*(I + I*r + 2*r^2*omega)))/((a^2 + (-2 + r)*r)^3*(a^2 + r^2))
        end
    elseif s == -2
        # s = -2
        return begin
            (2*(a^2 + (-2 + r)*r)*(I*a*m*r + a^2*(-2 - I*r*omega) + r*(3 - r - I*r^2*omega)))/r
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end