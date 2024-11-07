module ConversionFactors

const I = 1im

#=
The conversion factors here are *always* defined as
the ratio of the coefficient in Teukosky equation *over*
the coefficient in the Sasaki-Nakamura formalism
=#

function Ctrans(s::Int, m::Int, a, omega, lambda)
    if s == 0
        return 1
    elseif s == +1
        inv = 2*I*omega
        return 1/inv
    elseif s == -1
        return (2*I*omega)/lambda
    elseif s == +2
        inv = -(4*omega^2)
        return 1/inv
    elseif s == -2
        return begin
            (4*omega^2)/(-lambda*(2+lambda)+12*I*omega-12*a*m*omega+12*(a*omega)^2)
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function Binc(s::Int, m::Int, a, omega, lambda)
    if s == 0
        return 1
    elseif s == +1
        return -((2*I*omega)/(2 + lambda))
    elseif s == -1
        inv = -(2*I*omega)
        return 1/inv
    elseif s == +2
        return begin
            (4*omega^2)/(-24 - 10*lambda - lambda^2 - 12*I*omega - 12*a*m*omega + 12*a^2*omega^2)
        end
    elseif s == -2
        inv = -(4*omega^2)
        return 1/inv
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function Btrans(s::Int, m::Int, a, omega, lambda)
    if s == 0
        inv = sqrt(2)*sqrt(1 + sqrt(1 - a^2))
        return 1/inv
    elseif s == +1
        return begin
            (sqrt(2)*sqrt(1 + sqrt(1 - a^2))*((-a)*m + 2*(1 + sqrt(1 - a^2))*omega))/
            (2*a*m + 2*I*(2 + lambda))
        end
    elseif s == -1
        inv = begin
            (1/(1 + sqrt(1 - a^2))^(5/2))*(4*sqrt(2)*(2 - 2*a^2 + 2*sqrt(1 - a^2) - 
            a^2*sqrt(1 - a^2) + 2*I*a*m - I*a^3*m + 2*I*a*sqrt(1 - a^2)*m - 8*I*omega + 
            6*I*a^2*omega - 8*I*sqrt(1 - a^2)*omega + 2*I*a^2*sqrt(1 - a^2)*omega))
        end
        return 1/inv
    elseif s == +2
        return begin
            (2*sqrt(2)*(1 + sqrt(1 - a^2))^(3/2)*(2*a*(1 + sqrt(1 - a^2))*m*(1 + 8*I*omega) - 
            8*I*(1 + sqrt(1 - a^2))*omega*(-I + 4*omega) + 
            a^3*m*(-2 - sqrt(1 - a^2) - 4*I*(3 + sqrt(1 - a^2))*omega) + 
            I*a^4*(m^2 + 2*(I - 2*omega)*omega) + 2*a^2*((-I)*(1 + sqrt(1 - a^2))*m^2 + 
            omega*(5 + 3*sqrt(1 - a^2) + 8*I*(2 + sqrt(1 - a^2))*omega))))/
            (16*a*(1 + sqrt(1 - a^2))*m*(11 + 2*lambda + 6*I*omega) + 8*I*(1 + sqrt(1 - a^2))*
            (24 + lambda*(10 + lambda) + 12*I*omega) - 36*I*a^5*m*omega + 12*I*a^6*omega^2 - 
            8*a^3*m*(18 + 7*sqrt(1 - a^2) + (3 + sqrt(1 - a^2))*lambda - 6*I*sqrt(1 - a^2)*omega) + 
            a^4*(24*I*m^2 + I*(4 + lambda)*(6 + lambda) + 32*sqrt(1 - a^2)*omega + 8*sqrt(1 - a^2)*lambda*omega + 
            12*(5 + 2*lambda)*omega + 48*I*omega^2) + 4*a^2*(-24*I*(2 + sqrt(1 - a^2)) - 
            12*I*(1 + sqrt(1 - a^2))*m^2 - 10*I*(2 + sqrt(1 - a^2))*lambda - 
            I*(2 + sqrt(1 - a^2))*lambda^2 - 8*(1 + sqrt(1 - a^2))*lambda*omega + 
            4*omega*(1 - 2*sqrt(1 - a^2) - 6*I*(1 + sqrt(1 - a^2))*omega)))
        end
    elseif s == -2
        inv = begin
            -((1/(1 + sqrt(1 - a^2))^(3/2))*(4*(sqrt(2)*(-4 + 6*a^2 - 2*a^4 - 4*sqrt(1 - a^2) + 
            4*a^2*sqrt(1 - a^2) - 6*I*a*m + 6*I*a^3*m - 6*I*a*sqrt(1 - a^2)*m + 
            3*I*a^3*sqrt(1 - a^2)*m + 2*a^2*m^2 - a^4*m^2 + 2*a^2*sqrt(1 - a^2)*m^2 + 
            24*I*omega - 30*I*a^2*omega + 6*I*a^4*omega + 24*I*sqrt(1 - a^2)*omega - 
            18*I*a^2*sqrt(1 - a^2)*omega - 16*a*m*omega + 12*a^3*m*omega - 16*a*sqrt(1 - a^2)*m*omega + 
            4*a^3*sqrt(1 - a^2)*m*omega + 32*omega^2 - 32*a^2*omega^2 + 4*a^4*omega^2 + 
            32*sqrt(1 - a^2)*omega^2 - 16*a^2*sqrt(1 - a^2)*omega^2))))
        end
        return 1/inv
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

function Cinc(s::Int, m::Int, a, omega, lambda)
    if s == 0
        return 1/(sqrt(2)*sqrt(1 + sqrt(1 - a^2)))
    elseif s == +1
        inv = begin
            (4*sqrt(2)*((-I)*a*(1 + sqrt(1 - a^2))*m + a^2*(-1 - 2*I*omega) + 
            (1 + sqrt(1 - a^2))*(1 + 4*I*omega)))/(1 + sqrt(1 - a^2))^(3/2)
        end
        return 1/inv
    elseif s == -1
        return begin
            (-4*a*m + 3*a^3*m - 4*a*sqrt(1 - a^2)*m + a^3*sqrt(1 - a^2)*m + 16*omega - 16*a^2*omega + 
            2*a^4*omega + 16*sqrt(1 - a^2)*omega - 8*a^2*sqrt(1 - a^2)*omega)/
            (sqrt(2)*sqrt(1 + sqrt(1 - a^2))*(2 - a^2 + 2*sqrt(1 - a^2))*(a*m - I*lambda))
        end
    elseif s == +2
        inv = begin
            -((1/(1 + sqrt(1 - a^2))^(3/2))*(4*sqrt(2)*(-2*a*(1 + sqrt(1 - a^2))*m*(-3*I + 8*omega) + 
            a^3*m*(-3*I*(2 + sqrt(1 - a^2)) + 4*(3 + sqrt(1 - a^2))*omega) - 
            a^4*(2 + m^2 + 6*I*omega - 4*omega^2) + 4*(1 + sqrt(1 - a^2))*(-1 - 6*I*omega + 8*omega^2) + 
            2*a^2*(3 + 2*sqrt(1 - a^2) + (1 + sqrt(1 - a^2))*m^2 + 
            omega*(3*I*(5 + 3*sqrt(1 - a^2)) - 8*(2 + sqrt(1 - a^2))*omega)))))
        end
        return 1/inv
    elseif s == -2
        return begin
            -((2*sqrt(2)*(1 + sqrt(1 - a^2))^(3/2)*(-2*I*a*m + 2*I*a^3*m - 2*I*a*sqrt(1 - a^2)*m + 
            I*a^3*sqrt(1 - a^2)*m + 2*a^2*m^2 - a^4*m^2 + 2*a^2*sqrt(1 - a^2)*m^2 + 8*I*omega - 
            10*I*a^2*omega + 2*I*a^4*omega + 8*I*sqrt(1 - a^2)*omega - 6*I*a^2*sqrt(1 - a^2)*omega - 
            16*a*m*omega + 12*a^3*m*omega - 16*a*sqrt(1 - a^2)*m*omega + 4*a^3*sqrt(1 - a^2)*m*omega + 
            32*omega^2 - 32*a^2*omega^2 + 4*a^4*omega^2 + 32*sqrt(1 - a^2)*omega^2 - 
            16*a^2*sqrt(1 - a^2)*omega^2))/((1 + sqrt(1 - a^2))^4*(lambda*(2 + lambda) - 12*I*omega) + 
            24*a^3*(1 + sqrt(1 - a^2))*m*(-I + 2*(1 + sqrt(1 - a^2))*omega) + 
            4*a*(1 + sqrt(1 - a^2))^2*m*(6*I + 2*I*(1 + sqrt(1 - a^2))*lambda + 
            3*(1 + sqrt(1 - a^2))^2*omega) - 4*a^2*(1 + sqrt(1 - a^2))*
            (6 + (1 + sqrt(1 - a^2))*(-3 + 6*m^2) + 2*I*(1 + sqrt(1 - a^2))*
            (3 + (1 + sqrt(1 - a^2))*(-3 + lambda))*omega + 3*(1 + sqrt(1 - a^2))^3*omega^2) + 
            12*a^4*(1 - 2*(1 + sqrt(1 - a^2))*omega*(-I + omega + sqrt(1 - a^2)*omega))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

# These functions are redundant and are defined for convenience
# The conversion factor for Cref is identical to Btrans
function Cref(s::Int, m::Int, a, omega, lambda)
    return Btrans(s, m, a, omega, lambda)
end

# The conversion factor for Bref is identical to Ctrans
function Bref(s::Int, m::Int, a, omega, lambda)
    return Ctrans(s, m, a, omega, lambda)
end

function TeukolskyStarobinsky_abs_Csq(s::Int, m::Int, a, omega, lambda)
    if abs(s) != 2
        throw(AssertionError("|C|^2 is only defined for s = +/-2"))
    end

    # Here we assume that lambda is computed with s = -2
    if s == 2
        # −sAlm = sAlm + 2s
        # Compute what lambda should be for s = -2
        lambda_minus2 = lambda + 4
    else
        lambda_minus2 = lambda
    end

    E = lambda_minus2 + 2*a*m*omega - a^2*omega^2 + 2
    Q = E + a^2*omega^2 - 2*a*m*omega
    return begin
        (Q^2 + 4*a*omega*m - 4*a^2*omega^2)*((Q-2)^2 + 36*a*omega*m - 36*a^2*omega^2) + (2*Q - 1)*(96*a^2*omega^2 - 48a*omega*m) + 144*omega^2(1-a^2)
    end
end

function TeukolskyStarobinsky_Bsq(s::Int, m::Int, a, omega, lambda)
    if abs(s) != 1
        throw(AssertionError("|B|^2 is only defined for s = +/-1"))
    end

    # Here we assume that lambda is computed with s = -1
    if s == 1
        # −sAlm = sAlm + 2s
        # Compute what lambda should be for s = -1
        lambda_minus1 = lambda + 2
    else
        lambda_minus1 = lambda
    end

    E = lambda_minus1 + 2*a*m*omega - a^2*omega^2
    Q = E + a^2*omega^2 - 2*a*m*omega
    return begin
        Q^2 + 4*a*omega*m - 4*a^2*omega^2
    end
end

end