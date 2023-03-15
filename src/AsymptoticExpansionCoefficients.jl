module AsymptoticExpansionCoefficients

using ForwardDiff
using ..Kerr
using ..Coordinates

export outgoing_coefficient_at_inf, ingoing_coefficient_at_inf
export outgoing_coefficient_at_hor, ingoing_coefficient_at_hor

const I = 1im # Mathematica being Mathematica

# A helper function to evaluate nth derivative using AD with recursion
function nth_derivative(f, n::Int)
    if n > 1
        _ith_derivative(x) = ForwardDiff.derivative(nth_derivative(f, n-1), x)
    elseif n == 1
        # Stop the recursion, evaluate (1st order) derivative normally
        return x -> ForwardDiff.derivative(f, x)
    else
        throw(DomainError(n, "Does not understand the order n"))
    end
    return _ith_derivative
end

function ingoing_coefficient_at_inf(s::Int, m::Int, a, omega, lambda, order::Int)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    for up to (1/r)^3, which is probably more than enough
    
    But we have also shown a recurrence relation where
    one can generate as higher an order as one pleases.
    However, the recurrence relation that we have
    actually depends *all* previous terms so this function
    is designed to be evaluated recursively to build
    the full list of coefficients
    =#
    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    if order == 0
        return 1.0 # This is always 1.0
    elseif order == 1
        if s == 0
            return (-(1/2))*I*(lambda + 2*a*m*omega)
        elseif s == +1
            return begin
                -((I*(4 + lambda^2 + 8*a*m*omega + 2*lambda*(2 + a*m*omega)))/
                (2*(2 + lambda)))
            end
        elseif s == -1
            return (-(1/2))*I*(lambda + 2*a*m*omega)
        elseif s == +2
            return begin
                (I*(-lambda^3 - 2*lambda^2*(8 + a*m*omega) + 
                4*lambda*(-21 - 3*(I + 4*a*m)*omega + 7*a^2*omega^2) + 
                8*(-18 - (9*I + 23*a*m)*omega + 
                a*(11*a - 3*I*m - 3*a*m^2)*omega^2 + 
                3*a^3*m*omega^3)))/(2*(10*lambda + lambda^2 + 
                12*(2 + (I + a*m)*omega - a^2*omega^2)))
            end
        elseif s == -2
            return (-(1/2))*I*(2 + lambda + 2*a*m*omega)
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    elseif order == 2
        if s == 0
            return begin
                (1/8)*(-((-2 + lambda)*lambda) - 4*(I + a*m*(-1 + lambda))*omega - 
                4*a*m*(2*I + a*m)*omega^2)
            end
        elseif s == +1
            return begin
                -((lambda^3 + 4*lambda^2*(1 + a*m*omega) + 
                8*a*omega*(m*(2 + 2*I*omega) - a*omega + 3*a*m^2*omega) + 
                4*lambda*(1 + a*m*(5 + 2*I*omega)*omega + a^2*(-2 + m^2)*omega^2))/
                (8*(2 + lambda)))
            end
        elseif s == -1
            return begin
                (1/8)*(-lambda^2 + lambda*(2 - 4*a*m*omega) + 
                4*a*omega*(m + 2*a*omega - m*(2*I + a*m)*omega))
            end
        elseif s == +2
            return begin
                -((lambda^4 + 4*lambda^3*(5 + a*m*omega) + 
                4*lambda^2*(37 + 2*a*m*(13 + I*omega)*omega + 
                a^2*(-11 + m^2)*omega^2) - 
                8*lambda*(-60 + 8*a*m*(-11 - 2*I*omega)*omega + 
                a^2*(39 - 19*m^2)*omega^2 + 14*a^3*m*omega^3) - 
                16*(a^2*(34 + m^2*(-49 - 9*I*omega) + 3*I*omega)*
                omega^2 + a^3*m*(43 - 3*m^2 + 6*I*omega)*omega^3 + 
                3*a^4*(-4 + m^2)*omega^4 - 9*(4 + omega^2) + 
                2*a*m*omega*(-44 - 15*I*omega + 3*omega^2)))/
                (8*(10*lambda + lambda^2 + 12*(2 + (I + a*m)*omega - 
                a^2*omega^2))))
            end
        elseif s == -2
            return begin
                (1/8)*((-lambda)*(2 + lambda) - 4*(-3*I + a*m*(1 + lambda))*omega - 
                4*a*m*(2*I + a*m)*omega^2)
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    elseif order == 3
        if s == 0
            return begin
                (1/48)*(I*(-6 + lambda)*(-2 + lambda)*lambda + 
                2*(12 - 18*lambda + I*a*m*(12 + lambda*(-16 + 3*lambda)))*omega + 
                4*I*a*(3*a*m^2*(-2 + lambda) + 2*a*(-1 + lambda) + 
                6*I*m*(1 + lambda))*omega^2 + 
                8*I*a*m*(-8 + 6*I*a*m + a^2*(2 + m^2))*omega^3)
            end
        elseif s == +1
            return begin
                (1/(48*(2 + lambda)))*(I*(lambda^4 + 6*a*m*lambda^3*omega + 
                4*lambda^2*(-3 + (6*I + 4*a*m)*omega + 
                a*(-4*a + 6*I*m + 3*a*m^2)*omega^2) + 
                8*lambda*(-2 + (12*I - 5*a*m)*omega + 
                a*(12*I*m + a*(-4 + 9*m^2))*omega^2 + 
                a*m*(-8 + 6*I*a*m + a^2*(-4 + m^2))*omega^3) + 
                16*omega*(6*I + a^3*m*(-1 + 4*m^2)*omega^2 + 
                3*I*a^2*omega*(I - 2*omega + 4*m^2*omega) + 
                a*m*(-3 + 12*I*omega - 8*omega^2))))
            end
        elseif s == -1
            return begin
                (1/48)*I*(lambda^3 + lambda^2*(-8 + 6*a*m*omega) + 
                8*a*omega*(a*omega + m*(2 - 3*a*m*omega + 
                (-8 + 6*I*a*m + a^2*(-4 + m^2))*omega^2)) + 
                4*lambda*(3 + omega*(6*I + a*(-4*a*omega + 
                m*(-8 + 3*(2*I + a*m)*omega)))))
            end
        elseif s == +2
            return begin
                -((I*(-lambda^5 - 2*lambda^4*(10 + 3*a*m*omega) - 
                4*lambda^3*(37 + 2*a*m*(20 + 3*I*omega)*omega + 
                a^2*(-13 + 3*m^2)*omega^2) - 
                8*lambda^2*(60 + 2*a^2*(-29 + 3*m^2*(9 + I*omega))*
                omega^2 + a^3*m*(-31 + m^2)*omega^3 + 
                a*m*omega*(157 + 48*I*omega - 8*omega^2)) + 
                16*lambda*(a^2*(91 + m^2*(-210 - 81*I*omega) + 9*I*omega)*
                omega^2 + 2*a^3*m*(73 - 13*m^2 + 21*I*omega)*
                omega^3 + 3*a^4*(-10 + 7*m^2)*omega^4 - 
                9*(4 + omega^2) + 2*a*m*omega*(-116 - 63*I*omega + 
                29*omega^2)) + 96*a*omega*((-a^3)*m^4*omega^3 + 
                a*omega*(18 + 9*I*omega - 11*a^2*omega^2) + 
                a^2*m^3*omega^2*(-28 - 7*I*omega + a^2*omega^2) + 
                a*m^2*omega*(-70 - 55*I*omega + 2*(7 + 15*a^2)*
                omega^2 + 6*I*a^2*omega^3) + 
                m*(-36 - 36*I*omega + (25 + 47*a^2)*omega^2 + 
                I*(8 + 23*a^2)*omega^3 - 2*a^2*(4 + 5*a^2)*
                omega^4))))/(48*(10*lambda + lambda^2 + 
                12*(2 + (I + a*m)*omega - a^2*omega^2))))
            end
        elseif s == -2
            return begin
                (1/48)*(I*(-4 + lambda)*lambda*(2 + lambda) + 
                2*(6*(-4 + lambda) + I*a*m*(-24 + lambda*(-4 + 3*lambda)))*
                omega + 4*a*(-6*m*(-1 + lambda) + 
                I*a*(-6 + (2 + 3*m^2)*lambda))*omega^2 + 
                8*I*a*m*(-8 + 6*I*a*m + a^2*(2 + m^2))*omega^3)
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    else
        # In the future, evaluate the coefficient numerically
        # But for now, return 0
        return 0.0
    end
end

function outgoing_coefficient_at_inf(s::Int, m::Int, a, omega, lambda, order::Int)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    for up to (1/r)^3, which is probably more than enough
    
    But we have also shown a recurrence relation where
    one can generate as higher an order as one pleases.
    However, the recurrence relation that we have
    actually depends *all* previous terms so this function
    is designed to be evaluated recursively to build
    the full list of coefficients
    =#
    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    if order == 0
        return 1.0 # This is always 1.0
    elseif order == 1
        if s == 0
            return (1/2)*I*(lambda + 2*a*m*omega)
        elseif s == +1
            return (1/2)*I*(2 + lambda + 2*a*m*omega)
        elseif s == -1
            return (I*(lambda^2 + 2*a*m*(2 + lambda)*omega))/(2*lambda)
        elseif s == +2
            return (1/2)*I*(6 + lambda + 2*a*m*omega)
        elseif s == -2
            return begin
                (1/6)*I*(-6 + 7*lambda) + I*a*m*omega + 
                (2*I*(-3 + lambda)*lambda*(2 + lambda) + 
                24*(-3 - 3*I*a*m + lambda)*omega)/(-3*lambda*(2 + lambda) + 
                36*omega*(I - a*m + a^2*omega))
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    elseif order == 2
        if s == 0
            return begin
                (1/8)*(-lambda^2 + lambda*(2 - 4*a*m*omega) + 
                4*omega*(I + a*m + a*m*(2*I - a*m)*omega))
            end
        elseif s == +1
            return begin
                (1/8)*(-lambda^2 - 2*lambda*(1 + 2*a*m*omega) - 
                4*a*omega*(m - 2*a*omega - 2*I*m*omega + a*m^2*omega))
            end
        elseif s == -1
            return begin
                -((1/(8*lambda))*((-2 + lambda)*lambda^2 + 4*a*m*(-2 + lambda + lambda^2)*omega + 
                4*a*(-2*I*m*lambda + a*(2 - 2*lambda + m^2*(4 + lambda)))*omega^2))
            end
        elseif s == +2
            return begin
                (1/8)*(-lambda^2 - 2*lambda*(5 + 2*a*m*omega) - 
                4*(6 + (3*I + 5*a*m)*omega + a*m*(-2*I + a*m)*omega^2))
            end
        elseif s == -2
            return begin
                -((1/(8*lambda*(2 + lambda) - 96*omega*(I - a*m + a^2*omega)))*
                (lambda^2*(2 + lambda)^2 + 4*a*m*lambda*(16 + lambda*(14 + lambda))*
                omega + 4*(36 + a*(a*(10 - 11*lambda)*lambda - 
                2*I*m*(2 + lambda)*(6 + lambda) + 
                a*m^2*(60 + lambda*(30 + lambda))))*omega^2 + 
                16*a*(-6*m + a*(3*I - 9*I*m^2 + 
                a*m*(-15 + 3*m^2 - 7*lambda)))*omega^3 - 
                48*a^3*(-2*I*m + a*(-4 + m^2))*omega^4))
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    elseif order == 3
        if s == 0
            return begin
                (1/2)*(1 - I*a*m)*omega - (1/48)*I*((-6 + lambda)*(-2 + lambda)*lambda + 
                2*lambda*(-18*I + a*m*(-16 + 3*lambda))*omega + 
                4*a*(3*a*m^2*(-2 + lambda) + 2*a*(-1 + lambda) - 
                6*I*m*(1 + lambda))*omega^2 + 8*a*m*(-8 - 6*I*a*m + 
                a^2*(2 + m^2))*omega^3)
            end
        elseif s == +1
            return begin
                (-(1/48))*I*(lambda^3 + lambda^2*(-2 + 6*a*m*omega) + 
                4*lambda*(-2 - 2*(3*I + a*m)*omega + 
                a*(-4*a - 6*I*m + 3*a*m^2)*omega^2) + 
                8*omega*(-6*I + a^3*m*(-4 + m^2)*omega^2 + 
                3*a^2*omega*(-1 - 2*I*m^2*omega) - 
                a*m*(3 + 6*I*omega + 8*omega^2)))
            end
        elseif s == -1
            return begin
                -((1/(48*lambda))*(I*((-6 + lambda)*(-2 + lambda)*lambda^2 + 48*a*m*omega + 
                2*lambda*(-12*I*lambda + a*m*(-16 + lambda*(-10 + 3*lambda)))*omega + 
                4*a*(3*a*m^2*(-2 + lambda)*(4 + lambda) - 
                4*a*(3 + (-2 + lambda)*lambda) - 6*I*m*(4 + lambda^2))*omega^2 + 
                8*a*(-8*m*lambda - 6*I*a*(-2 + m^2*(2 + lambda)) + 
                a^2*m*(6 - 4*lambda + m^2*(6 + lambda)))*omega^3)))
            end
        elseif s == +2
            return begin
                (-(1/48))*I*(lambda^3 + 2*lambda^2*(5 + 3*a*m*omega) + 
                4*lambda*(6 + (3*I + 10*a*m)*omega + 
                a*(2*a - 6*I*m + 3*a*m^2)*omega^2) + 
                8*a*omega*(a*omega + 6*a*m^2*(1 - I*omega)*omega + 
                a^2*m^3*omega^2 + m*(2 - 9*I*omega + 2*(-4 + a^2)*
                omega^2)))
            end
        elseif s == -2
            return begin
                -((1/(48*(lambda*(2 + lambda) - 12*omega*(I - a*m + a^2*omega))))*
                (I*(-4 + lambda)*lambda^2*(2 + lambda)^2 + 
                2*I*a*m*lambda*(-96 + lambda*(-44 + lambda*(32 + 3*lambda)))*
                omega + 4*I*(36*(-4 + lambda) + 
                a*(-6*I*m*lambda*(2 + lambda)^2 + 
                a*lambda*(-60 + (40 - 13*lambda)*lambda) + 
                3*a*m^2*(-48 + lambda*(40 + lambda*(24 + lambda)))))*
                omega^2 + 8*I*a*(-6*I*a*(-3*(2 + lambda) + 
                m^2*(1 + lambda)*(18 + lambda)) - 
                4*m*(-9 + lambda*(13 + 2*lambda)) + 
                a^2*m*(108 - lambda*(44 + 31*lambda) + 
                m^2*(144 + lambda*(44 + lambda))))*omega^3 - 
                48*a*(16*m + 28*I*a*m^2 - 2*a^2*m*
                (5 + 7*m^2 - 7*lambda) - 
                I*a^3*(2*m^4 + 2*(-9 + 5*lambda) - 
                m^2*(32 + 7*lambda)))*omega^4 - 
                96*I*a^3*m*(-8 - 6*I*a*m + a^2*(-10 + m^2))*
                omega^5))
            end
        else
            throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
        end
    else
        # In the future, evaluate the coefficient numerically
        # But for now, return 0
        return 0.0
    end
end

function outgoing_coefficient_at_hor(s::Int, m::Int, a, omega, lambda, order::Int)
    #=
    We have derived a similar recurrence relation for the expansion
    coefficients about the horizon. However the explicit expressions
    contain so many terms that they might be better off evaluated
    numerically instead.
    =#
    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    if order == 0
        return 1.0 # This is always 1.0
    else
        return 0.0 # We leave the implementation in the future
    end
end

function PminusH_r_minus_rp(s::Int, m::Int, a, omega, lambda, x)
    if s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/(2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + I*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega) - 
            (8*a*x*(2*sqrt(1 - a^2) + x)*((-I)*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*I*omega) + 24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function QminusH_r_minus_rp_sq(s::Int, m::Int, a, omega, lambda, x)
    if s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 + 
            (8*I*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*((-I)*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*lambda) + 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 
            I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2))) - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*((a^2*(-1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^3)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(-a^4 - 8*a^2*(1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^4) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*(a^2*(-1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^3)*
            ((-I)*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*I - 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*I*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*I*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-I + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*I*(1 + sqrt(1 - a^2) + x)*omega - 24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(-24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*
            (-3 + 3*I*(1 + sqrt(1 - a^2) + x)*omega + 4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (72*I - 4*I*(1 + sqrt(1 - a^2) + x)*(15 + 9*m^2 - 13*lambda) + 4*(1 + sqrt(1 - a^2) + x)^2*(I*m^2*(6 + lambda) - I*(3 + lambda^2 + lambda*(9 + 16*I*omega) + 24*I*omega)) + 
            (1 + sqrt(1 - a^2) + x)^3*(-36 + 30*m^2 - 44*lambda - lambda^2 - 36*I*omega)*omega + 4*I*(1 + sqrt(1 - a^2) + x)^4*(3 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 + 2*I*(1 + sqrt(1 - a^2) + x)^3*omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 
            4*a^6*(1 + sqrt(1 - a^2) + x)*(12 - 3*(1 + sqrt(1 - a^2) + x)*(-4 + 3*m^2 + 6*I*omega) + (1 + sqrt(1 - a^2) + x)^3*(-39 + 36*m^2 - 8*lambda - 18*I*omega)*omega^2 + 
            2*I*(1 + sqrt(1 - a^2) + x)^4*(3 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*omega^4 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(24*I + 9*I*m^2 - 4*I*lambda + 15*omega)) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-18*I + 2*(1 + sqrt(1 - a^2) + x)^2*(-36 + 12*m^2 - 8*lambda - 27*I*omega)*omega + 6*I*(1 + sqrt(1 - a^2) + x)^3*(4 + lambda)*omega^2 + 
            45*(1 + sqrt(1 - a^2) + x)^4*omega^3 + 2*(1 + sqrt(1 - a^2) + x)*(9*I + 3*I*m^2 - 4*I*lambda + 30*omega)) - 2*a*m*(1 + sqrt(1 - a^2) + x)^4*
            (-48*I - 24*I*(1 + sqrt(1 - a^2) + x)*(-2 + lambda) + 4*I*(1 + sqrt(1 - a^2) + x)^2*(7*lambda + 2*lambda^2 + 6*I*omega) - (1 + sqrt(1 - a^2) + x)^4*(12*lambda + lambda^2 - 24*I*omega)*
            omega + 4*I*(1 + sqrt(1 - a^2) + x)^5*lambda*omega^2 + 6*(1 + sqrt(1 - a^2) + x)^6*omega^3 - 4*I*(1 + sqrt(1 - a^2) + x)^3*(lambda + lambda^2 + 6*I*omega + 5*I*lambda*omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^3*(24*m^4*(1 + sqrt(1 - a^2) + x) + 48*(-4 + 3*I*omega) + 8*(1 + sqrt(1 - a^2) + x)*(9 + lambda^2 + lambda*(2 + 13*I*omega) - 72*I*omega) - 
            (1 + sqrt(1 - a^2) + x)^3*(60 + 54*lambda + lambda^2 + 36*I*omega)*omega^2 + 8*I*(1 + sqrt(1 - a^2) + x)^4*(-3 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^4 + 
            8*(1 + sqrt(1 - a^2) + x)^2*omega*(24*I - 6*I*lambda - I*lambda^2 + 18*omega + 8*lambda*omega) + 4*m^2*(30 - (1 + sqrt(1 - a^2) + x)*(33 + 8*lambda + 54*I*omega) + 
            6*I*(1 + sqrt(1 - a^2) + x)^2*(5 + lambda)*omega + 45*(1 + sqrt(1 - a^2) + x)^3*omega^2)) + (1 + sqrt(1 - a^2) + x)^6*
            ((-1 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)*lambda^3 + lambda^2*(12 - 12*(1 + sqrt(1 - a^2) + x) + 2*(1 + sqrt(1 - a^2) + x)^2 - 
            (1 + sqrt(1 - a^2) + x)^4*omega^2) + 12*I*omega*(-12 + 8*(1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            2*lambda*(-12 + 4*(1 + sqrt(1 - a^2) + x)*(2 - 3*I*omega) + 6*I*(1 + sqrt(1 - a^2) + x)^2*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^4*(96 + 8*(1 + sqrt(1 - a^2) + x)*(m^2*(6 + 8*lambda) - 3*(2 + lambda^2 + lambda*(2 + 2*I*omega) - 22*I*omega)) - 96*I*omega + 
            2*(1 + sqrt(1 - a^2) + x)^4*(6 + 24*m^2 - 12*lambda - lambda^2 + 24*I*omega)*omega^2 + 8*I*(1 + sqrt(1 - a^2) + x)^5*(-3 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^6*omega^4 + 
            (1 + sqrt(1 - a^2) + x)^2*(lambda^3 + lambda*(24 - 34*m^2 - 4*I*omega) + lambda^2*(14 - m^2 + 16*I*omega) - 12*I*(22 + 3*m^2 - 4*I*omega)*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^3*omega*((-I)*lambda^2 + 2*I*m^2*(3 + lambda) + 3*omega + lambda*(2*I + 5*omega)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(2*lambda + lambda^2 - 12*I*omega) + 24*a^3*m*(1 + sqrt(1 - a^2) + x)*(I - 2*(1 + sqrt(1 - a^2) + x)*omega) - 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*I + 2*I*(1 + sqrt(1 - a^2) + x)*lambda + 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 - 2*I*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*I*omega) + 2*I*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*omega^2))))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

# Cache mechanism for the ingoing coefficients at horizon
# Initialize the cache with a set of fiducial parameters
_cached_ingoing_coefficients_at_hor_params = (s=-2, m=2, a=0, omega=0.5, lambda=1)
_cached_ingoing_coefficients_at_hor = (
    expansion_coeffs = [ComplexF64(1.0)], 
    Pcoeffs = [ComplexF64(PminusH_r_minus_rp(
        _cached_ingoing_coefficients_at_hor_params.s,
        _cached_ingoing_coefficients_at_hor_params.m,
        _cached_ingoing_coefficients_at_hor_params.a,
        _cached_ingoing_coefficients_at_hor_params.omega,
        _cached_ingoing_coefficients_at_hor_params.lambda,
        0
        ))],
    Qcoeffs = [ComplexF64(0.0)]
)

function ingoing_coefficient_at_hor(s::Int, m::Int, a, omega, lambda, order::Int)
    global _cached_ingoing_coefficients_at_hor_params
    global _cached_ingoing_coefficients_at_hor

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    _this_params = (s=s, m=m, a=a, omega=omega, lambda=lambda)
    # Check if we can use the cached results
    if _cached_ingoing_coefficients_at_hor_params == _this_params
        expansion_coeffs = _cached_ingoing_coefficients_at_hor.expansion_coeffs
        Pcoeffs = _cached_ingoing_coefficients_at_hor.Pcoeffs
        Qcoeffs = _cached_ingoing_coefficients_at_hor.Qcoeffs
        _start_order = length(expansion_coeffs)
    else
        # Cannot re-use the cached results, re-compute from zero
        expansion_coeffs = [ComplexF64(1.0)] # order 0
        Pcoeffs = [ComplexF64(PminusH_r_minus_rp(s, m, a, omega, lambda, 0))] # order 0
        Qcoeffs = [ComplexF64(0.0)] # order 0
        _start_order = 1
    end

    if order > 0
        # Compute series expansion coefficients for P and Q
        _P(x) = PminusH_r_minus_rp(s, m, a, omega, lambda, x)
        _Q(x) = QminusH_r_minus_rp_sq(s, m, a, omega, lambda, x)

        for i in _start_order:order
            println("Expanding $_start_order")
            println(Pcoeffs)
            println(Qcoeffs)
            append!(Pcoeffs, nth_derivative(_P, i)(0)/factorial(i))
            append!(Qcoeffs, nth_derivative(_Q, i)(0)/factorial(i))
            println("Expansion done!")
            println(Pcoeffs)
            println(Qcoeffs)           
        end
    end
    # Define the indicial polynomial
    indicial(nu) = nu*(nu - 1) + Pcoeffs[1]*nu + Qcoeffs[1]

    if order > 0
        # Evaluate the C coefficients
        for i in _start_order:order
            sum = 0.0
            for r in 0:i-1
                sum += expansion_coeffs[r+1]*(r*Pcoeffs[i-r+1] + Qcoeffs[i-r+1])
            end
            append!(expansion_coeffs, -sum/indicial(i))
        end
    end

    # Update cache
    _cached_ingoing_coefficients_at_hor_params = _this_params
    _cached_ingoing_coefficients_at_hor = (
        expansion_coeffs = expansion_coeffs,
        Pcoeffs = Pcoeffs,
        Qcoeffs = Qcoeffs
    )
    return expansion_coeffs[order+1]
end

end