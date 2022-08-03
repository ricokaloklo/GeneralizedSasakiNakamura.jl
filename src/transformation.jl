const I = 1im # Mathematica being Mathematica

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

function alpha_prime(s::Int, m::Int, a, omega, lambda, r)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    if s == 0
        # s = 0
        return 0.0
    elseif s == 1
        # s = +1
        return begin
            (sqrt(1/(a^2 + (-2 + r)*r))*(a^2 - 2*I*a*m*r + 3*r^2 + 4*I*r*(a^2 + r^2)*omega))/(a^2 + r^2) - 
            (2*r*sqrt(1/(a^2 + (-2 + r)*r))*(-2*a^2 - I*a^3*m + a^2*r - I*a*m*r^2 + r^3 + 
            I*(a^2 + r^2)^2*omega))/(a^2 + r^2)^2 - ((-2 + 2*r)*(1/(a^2 + (-2 + r)*r))^(3/2)*
            (-2*a^2 - I*a^3*m + a^2*r - I*a*m*r^2 + r^3 + I*(a^2 + r^2)^2*omega))/(2*(a^2 + r^2))
        end
    elseif s == -1
        # s = -1
        return begin
            (-1 + (2*r)/(a^2 + r^2) + r*(-((4*r^2)/(a^2 + r^2)^2) + 2/(a^2 + r^2)) - 2*I*r*omega)/
            sqrt(a^2 + (-2 + r)*r) - ((-2 + 2*r)*(I*a*m + r*(-1 + (2*r)/(a^2 + r^2)) - I*(a^2 + r^2)*omega))/
            (2*(a^2 + (-2 + r)*r)^(3/2))
        end
    elseif s == 2
        # s = +2
        return begin
            -((1/((a^2 + (-2 + r)*r)^2*(a^2 + r^2)^3))*
            (2*(38*I*a^5*m*(-1 + r)*r^2 + a^9*m*(3*I - 4*omega) + 
            2*a^10*(1 + r)*omega^2 + 2*r^8*(1 + 5*I*omega + 
            (-3 + r)*r^2*omega^2) + 2*a^2*r^4*(36 - 24*r - 
            r^3*(m^2 + 10*I*omega) + r^2*(2 + m^2 + 20*I*omega) - 
            11*r^4*omega^2 + 5*r^5*omega^2) + 
            2*a^6*(-4 + 24*r - 3*r^3*(m^2 + 10*I*omega) + 
            r^2*(-2 + 3*m^2 + 20*I*omega) - 6*r^4*omega^2 + 
            10*r^5*omega^2) + a^4*(-64*r^2 + 
            6*r^4*(m^2 + 10*I*omega) - 6*r^5*(m^2 + 10*I*omega) - 
            28*r^6*omega^2 + 20*r^7*omega^2) - 
            2*a^8*(1 + m^2*(-1 + r) + 5*I*(-1 + 2*r)*omega - 
            r^2*(1 + 5*r)*omega^2) + a*m*r^6*(-2*I + 2*I*r + 
            r^2*(-3*I + 4*omega)) - 2*a^7*m*(5*I - 9*I*r + 
            r^2*(-3*I + 4*omega)) + 2*a^3*m*r^4*(-15*I + 11*I*r + 
            r^2*(-3*I + 4*omega)))))
        end
    elseif s == -2
        # s = -2
        return begin
            -((1/(r^3*(a^2 + (-2 + r)*r)^2))*(2*(-2*I*a^5*m*r + a^6*(6 + 2*I*r*omega) + 
            2*a^4*r*(-15 + r*(6 - 4*I*omega) + r^2*omega^2 + r^3*omega^2) + 
            2*a^2*r^2*(24 + r*(-18 + m^2 + 3*I*omega) + r^2*(3 - m^2 + 6*I*omega) + 2*r^4*omega^2 - 
            r^3*omega*(3*I + 2*omega)) + 2*r^3*(-12 + 12*r + r^2*(-3 - 9*I*omega) + 8*I*r^3*omega + r^5*omega^2 - 
            r^4*omega*(2*I + 3*omega)) + a*m*r^3*(-6*I + 6*I*r + r^2*(-I + 4*omega)) + 
            a^3*m*r^2*(8*I - r*(5*I + 4*omega)))))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function beta_prime(s::Int, m::Int, a, omega, lambda, r)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    if s == 0
        # s = 0
        return 0.0
    elseif s == 1
        # s = +1
        return (-(3/2))*(-2 + 2*r)*(1/(a^2 + (-2 + r)*r))^(5/2)
    elseif s == -1
        # s = -1
        return (-2 + 2*r)/(2*sqrt(a^2 + (-2 + r)*r))
    elseif s == 2
        # s = +2
        return begin
            -((2*I*(-6*a^5*m*(-1 + r) - 12*a^3*m*(-1 + r)*r^2 - 
            6*a*m*(-1 + r)*r^4 + r^4*(18*I - 22*I*r + 
            r^2*(5*I - 2*omega) + 4*r^3*omega) + 
            a^4*(-6*I + 10*I*r + r^2*(3*I - 14*omega) + 12*r^3*omega) + 
            a^2*r^2*(-4*I - 12*I*r + r^2*(9*I - 10*omega) + 
            12*r^3*omega) + a^6*(-I + (-6 + 4*r)*omega)))/
            ((a^2 + (-2 + r)*r)^4*(a^2 + r^2)^2))
        end
    elseif s == -2
        # s = -2
        return begin
            2*(-6 + 2*I*a*m*(-1 + r) + (2*a^4)/r^2 + 10*r + r^2*(-3 + 6*I*omega) - 4*I*r^3*omega + 
            a^2*(-3 - 2*I*(-1 + 2*r)*omega))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end