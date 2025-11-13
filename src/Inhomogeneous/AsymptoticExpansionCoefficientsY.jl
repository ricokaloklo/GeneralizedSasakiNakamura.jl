module AsymptoticExpansionCoefficientsY

using SpecialFunctions
using ..AsymptoticExpansionCoefficients: outgoing_coefficient_at_inf, ingoing_coefficient_at_inf, outgoing_coefficient_at_hor, ingoing_coefficient_at_hor 

export Ypp_outgoing_inf_m2, Ypp_ingoing_inf_m2, Ypp_outgoing_hor_m2, Ypp_ingoing_hor_m2, Ypp_outgoing_inf_p2, Ypp_ingoing_inf_p2, Ypp_outgoing_hor_p2, Ypp_ingoing_hor_p2

_DEFAULTDATATYPE = ComplexF64 # Double precision by default
_DEFAULT_infinity_expansion_order = 6

#=  Calaulate the asymptotic coefficients f_plus_inf and f_minus_inf deined in Eq.(F9) of RL2024 where they are denoted by a_+/-.
    Note that the expressions for order=1,2,3 are given explicitly in Appendix G of RL2024. 
    For order >= 4, the coefficients are calculated the recurrence relation given in Eq.(F12) of RL2024.
    One can also check that the coefficients are 0 for order=0.
=#

function a_plus_inf(s::Int, m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)

    C_plus_inf = outgoing_coefficient_at_inf(s, m, a, omega, lambda, order)

    # The relation between C_plus_inf and a_plus_inf is given by Eq.(F13) of RL2024

    aplus_inf = data_type(C_plus_inf / omega^order)

    return aplus_inf
end

function a_minus_inf(s::Int, m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)

    C_minus_inf = ingoing_coefficient_at_inf(s, m, a, omega, lambda, order)

    # The relation between C_minus_inf and a_minus_inf is given by Eq.(F13) of RL2024

    aminus_inf = data_type(C_minus_inf / omega^order)

    return aminus_inf
end

function a_plus_hor(s::Int, m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)

    C_plus_hor = outgoing_coefficient_at_hor(s, m, a, omega, lambda, order)

    # The relation between C_plus_hor and a_plus_hor is given by Eq.(F19) of RL2024

    aplus_hor = data_type(C_plus_hor * omega^order)

    return aplus_hor
end

function a_minus_hor(s::Int, m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)

    C_minus_hor = ingoing_coefficient_at_hor(s, m, a, omega, lambda, order)

    # The relation between C_minus_hor and a_minus_hor is given by Eq.(F19) of RL2024

    aminus_hor = data_type(C_minus_hor * omega^order)

    return aminus_hor
end

# The a_n, b_n, d_n coefficients

# The complete Bell polynomial B_n is needed to compute the coefficients a_n and d_n

function Complete_Bell_Poly(x_vec::Vector{}; data_type=_DEFAULTDATATYPE)
    # Compute the complete Bell polynomial B_n(x) for a vector of coefficients x_vec
    x_vec = ComplexF64.(x_vec)  # Ensure the vector is of type ComplexF64
    n = length(x_vec)
    B_n = data_type(0)
    #=  Use the recurrence relation defined by
        B_{n+1}(x_1,...,x_{n+1}) = sum_{i=0}^{n} binomial(n, i) * x[i+1] * B_{n-i}(x_1,...,x_{n-i})
        which is equivalent to
        B_n(x_1,...,x_n) = sum_{i=0}^{n-1} binomial(n-1, i) * x[i+1] * B_{n-1-i}(x_1,...,x_{n-1-i})
    =#
    if n == 0
        return data_type(1)  # B_0(x) = 1.0 + 0.0im
    else
        for i in 0:n-1
            x_vec_reduced = x_vec[1:n-1-i] # Reduce the vector by removing the last i elements
            # Compute the binomial coefficient
            B_n += binomial(n-1, i) * x_vec[i+1] * Complete_Bell_Poly(x_vec_reduced)
        end
        return B_n
    end
end

function a_coefficient_m2(m, a, order::Int; data_type=_DEFAULTDATATYPE)
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    P_vec = ComplexF64[1im*a*m*(rp^i-rm^i)*gamma(i)/(rm-rp) for i in 1:order]
    a_ord = Complete_Bell_Poly(P_vec)/gamma(order+1)
    return data_type(a_ord)
end

function a_coefficient_p2(m, a, order::Int; data_type=_DEFAULTDATATYPE)
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    P_vec = ComplexF64[1im*a*m*(rp^i-rm^i)*gamma(i)/(rp-rm) for i in 1:order]
    a_ord = Complete_Bell_Poly(P_vec)/gamma(order+1)
    return data_type(a_ord)
end

function b_coefficient(a, order::Int; data_type=_DEFAULTDATATYPE)
    if isodd(order)
        return data_type(0)
    end
    return data_type(binomial(-1/2, Int.(order/2)) * a^order)
end

function d_coefficient_m2(a, omega, order::Int; data_type=_DEFAULTDATATYPE)
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    Q_vec = ComplexF64[4im*omega*gamma(i)*sum(rp^j*rm^(i-j) for j in 0:i) for i in 1:order]
    d_ord = Complete_Bell_Poly(Q_vec)/gamma(order+1)
    return d_ord
end

function d_coefficient_p2(a, omega, order::Int; data_type=_DEFAULTDATATYPE)
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    Q_vec = ComplexF64[-4im*omega*gamma(i)*sum(rp^j*rm^(i-j) for j in 0:i) for i in 1:order]
    d_ord = Complete_Bell_Poly(Q_vec)/gamma(order+1)
    return d_ord
end


# The asymptotic coefficients Ypp_plus_inf and Ypp_minus_inf

function Ypp_outgoing_inf_m2(m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    # Compute the outgoing expansion coefficients of Y'' at infinity

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    if order == 0
        return data_type(1)
    else
        fplus = data_type[a_plus_inf(-2, m, a, omega, lambda, ord) for ord in 0:order]
        a_ord = data_type[a_coefficient_m2(m, a, ord) for ord in 0:order]
        b_ord = data_type[b_coefficient(a, ord) for ord in 0:order]

        Ypp_outgoing = sum(
            sum(
                a_ord[i + 1] * b_ord[j + 1] * fplus[order - i - j + 1]
                for j in 0:order-i
            )
            for i in 0:order
        )
        return Ypp_outgoing
    end
end

function Ypp_ingoing_inf_p2(m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    # Compute the outgoing expansion coefficients of Y'' at infinity

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    if order == 0
        return data_type(1)
    else
        fminus = data_type[a_minus_inf(2, m, a, omega, lambda, ord) for ord in 0:order]
        a_ord = data_type[a_coefficient_p2(m, a, ord) for ord in 0:order]
        b_ord = data_type[b_coefficient(a, ord) for ord in 0:order]

        Ypp_ingoing = sum(
            sum(
                a_ord[i + 1] * b_ord[j + 1] * fminus[order - i - j + 1]
                for j in 0:order-i
            )
            for i in 0:order
        )
        return Ypp_ingoing
    end
end

function Ypp_ingoing_inf_m2(m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    # Compute the ingoing expansion coefficients of Y'' at infinity

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    YIm0 = exp(4im * omega * log(2))
    if order == 0
        return YIm0
    else
        fminus = data_type[a_minus_inf(-2, m, a, omega, lambda, ord) for ord in 0:order]
        a_ord = data_type[a_coefficient_m2(m, a, ord) for ord in 0:order]
        b_ord = data_type[b_coefficient(a, ord) for ord in 0:order]
        d_ord = data_type[d_coefficient_m2(a, omega, ord) for ord in 0:order]

        Ypp_ingoing = sum(
            sum(
                sum(a_ord[i + 1] * b_ord[j + 1] * fminus[k + 1] * d_ord[order - i - j - k + 1]
                    for k in 0:order-i-j)
                for j in 0:order-i)
            for i in 0:order)
        return Ypp_ingoing * YIm0
    end
end

function Ypp_outgoing_inf_p2(m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    # Compute the ingoing expansion coefficients of Y'' at infinity

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    YIm0 = exp(- 4im * omega * log(2))
    if order == 0
        return YIm0
    else
        fplus = data_type[a_plus_inf(2, m, a, omega, lambda, ord) for ord in 0:order]
        a_ord = data_type[a_coefficient_p2(m, a, ord) for ord in 0:order]
        b_ord = data_type[b_coefficient(a, ord) for ord in 0:order]
        d_ord = data_type[d_coefficient_p2(a, omega, ord) for ord in 0:order]

        Ypp_ingoing = sum(
            sum(
                sum(a_ord[i + 1] * b_ord[j + 1] * fplus[k + 1] * d_ord[order - i - j - k + 1]
                    for k in 0:order-i-j)
                for j in 0:order-i)
            for i in 0:order)
        return Ypp_ingoing * YIm0
    end
end

function Ypp_outgoing_hor_m2(m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    # Compute the outgoing expansion coefficients of Y'' at the horizon

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    ap = sqrt(1 - a^2)
    YHp0 = exp(-(im * m * (a^2 - (-1 + ap) * log(1 - a^2))) / (2 * a)) / (sqrt(2) * (1 + ap)^(5/2))

    if order == 0
        return YHp0
    elseif order == 1
        CHp1 = a_plus_hor(-2, m, a, omega, lambda, 1; data_type=data_type)
        return (-(1/2) + (2*(-1 + ap))/a^2 - (im*a*m)/(2*ap) + CHp1) * YHp0
    elseif order == 2
        CHp1 = a_plus_hor(-2, m, a, omega, lambda, 1; data_type=data_type)
        CHp2 = a_plus_hor(-2, m, a, omega, lambda, 2; data_type=data_type)
        return ((1 / (8 * a^4 * (-1 + a^2))) * (48*(-1 + ap) + a^2*(66 - 42*ap - 3*a^2*(7 + 2*ap) 
                - 7*im*a*(-1 + ap)*m - 2*im*a^3*(4 + ap)*m + a^4*(3 + m^2)) - 4*a^2*(4*(-1 + ap) 
                + a^2*(3 + a^2 - 4*ap - im*a*ap*m)) * CHp1 + 8*a^4*(-1 + a^2)*CHp2)) * YHp0
    elseif order == 3
        CHp1 = a_plus_hor(-2, m, a, omega, lambda, 1; data_type=data_type)
        CHp2 = a_plus_hor(-2, m, a, omega, lambda, 2; data_type=data_type)
        CHp3 = a_plus_hor(-2, m, a, omega, lambda, 3; data_type=data_type)
        return ((1 / (48 * a^6 * (-1 + a^2)^2)) * (768*(-1 + ap) + a^2*(96*(21 - 17*ap) + 6*a^2*(-283 + 163*ap) 
                + im*a^5*(97 - 66*ap)*m - 120*im*a*(-1 + ap)*m + 5*im*a^3*(-47 + 35*ap)*m - 3*a^8*(5 + m^2) 
                + im*a^7*m*(9*(2 + ap) + ap*m^2) + a^4*(405 - 132*ap - 9*(-1 + ap)*m^2) + 6*a^6*(10 + 3*ap 
                + (-1 + 2*ap)*m^2)) + 6*a^2*(-1 + a^2)*(48*(-1 + ap) + a^2*(66 - 42*ap - 3*a^2*(7 + 2*ap) 
                - 7*im*a*(-1 + ap)*m - 2*im*a^3*(4 + ap)*m + a^4*(3 + m^2)) * CHp1 - 4*a^2*(4*(-1 + ap) 
                + a^2*(3 + a^2 - 4*ap - im*a*ap*m)) * CHp2 + 8*a^4*(-1 + a^2)*CHp3))) * YHp0
    else
        throw(DomainError(order, "Currently horizontal expansion orders only up to 3 are supported"))
    end
end

function Ypp_outgoing_hor_p2(m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    # Compute the outgoing expansion coefficients of Y'' at the horizon

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    ap = sqrt(1 - a^2)
    YHp0 = (2^(-((-2*im*a*(1 + ap)*m + a^2*(-1 - 4*im*omega) + (1 + ap)*(1 + 8*im*omega)) / (2*(1 - a^2 + ap)))) 
            * (1 - a^2)^(((im*a*(m - 2*a*omega)) / (2*(1 - a^2 + ap)))) * exp(-((im*(-a^3*m + a*(1 + ap)*m 
            - 8*(1 + ap)*omega + 4*a^2*(2 + ap)*omega)) / (2*(1 - a^2 + ap))))) / (1 + ap)^(5/2)

    if order == 0
        return YHp0
    elseif order == 1
        CHp1 = a_plus_hor(2, m, a, omega, lambda, 1; data_type=data_type)
        return ((1 / (2 * ap * (1 - a^2 + ap))) * (-5 - ap + im*a^3*m + 4*im*omega + 4*im*ap*omega 
                + a^2*(5 + ap - 2*im*(3 + 2*ap)*omega) - 2*(-1 + a^2)*(1 + ap)*CHp1)) * YHp0
    elseif order == 2
        CHp1 = a_plus_hor(2, m, a, omega, lambda, 1; data_type=data_type)
        CHp2 = a_plus_hor(2, m, a, omega, lambda, 2; data_type=data_type)
        return ((1 / (8 * a^4 * (-1 + a^2)^2)) * ( -48*(-1 + ap) + a^2*(-114 + 90*ap + a*(15*im*(-1 + ap)*m 
                - 7*im*a^2*(-3 + 2*ap)*m - 2*a^4*m*(im*(3 + ap) + (-2 + 4*ap)*omega) - a*(-87 + 36*ap 
                + 2*im*(-5 + 17*ap)*omega + 8*(1 + ap)*omega^2) + a^5*(3 + m^2 - 8*omega*(im + 2*omega)) 
                + 2*a^3*(-3*(4 + ap) + (-1 + ap)*m^2 + omega*(-im + 18*im*ap + 10*omega + 8*ap*omega)) )) 
                + 4*a^2*(-1 + a^2)*(4 - 4*ap + a^2*(-3 + 4*ap + im*a*(-1 + ap)*m + a^2*(-1 + 4*im*omega) 
                - 2*im*(1 + ap)*omega)) * CHp1 + 8*a^4*(-1 + a^2)^2*CHp2 )) * YHp0
    elseif order == 3
        CHp1 = a_plus_hor(2, m, a, omega, lambda, 1; data_type=data_type)
        CHp2 = a_plus_hor(2, m, a, omega, lambda, 2; data_type=data_type)
        CHp3 = a_plus_hor(2, m, a, omega, lambda, 3; data_type=data_type)
        return ((1 / (48 * a^6 * (-1 + a^2)^3)) * ( -768*(-1 + ap) + a^2*(96*(-29 + 25*ap) 
                + a*(im*a^2*(709 - 577*ap)*m + 264*im*(-1 + ap)*m + a^4*m*(im*(-619 + 358*ap) 
                + 6*(20 - 21*ap)*omega) + 6*a*(619 - 435*ap - 64*im*(-1 + ap)*omega) + a^8*m*(9*im*(1 + ap) 
                + im*(-3 + ap)*m^2 + 84*omega + 24*ap*(1 - 2*im*omega)*omega) + a^6*m*(165*im - 54*im*ap 
                - 4*im*(-1 + ap)*m^2 + 6*omega*(-34 + 19*ap + 2*im*(1 + ap)*omega)) + a^9*(-15 + m^2*(-3 
                + 12*im*omega) + 4*omega*(9*im + 4*(3 - 4*im*omega)*omega)) + a^3*(3*(-701 + 370*ap) 
                + 42*(-1 + ap)*m^2 + 8*im*omega*(-130 + 88*ap + (1 + ap)*omega*(9*im + 4*omega))) + 3*a^7*(25 
                + 6*ap + m^2*(-9 + 2*ap + 6*im*(-1 + ap)*omega) + 2*im*omega*(-58 - 15*ap + 2*omega*(11*im 
                + 20*im*ap + 12*omega + 8*ap*omega))) + a^5*(345 - 150*ap + (72 - 51*ap)*m^2 + 2*omega*(484*im 
                - 115*im*ap + 2*omega*(39 + 75*ap - 2*im*(15 + 13*ap)*omega))))) - 6*a^2*(-1 + a^2)*(48*(-1 + ap) 
                + a^2*(114 - 90*ap - 15*im*a*(-1 + ap)*m + 7*im*a^3*(-3 + 2*ap)*m + 2*a^5*m*(im*(3 + ap) 
                + (-2 + 4*ap)*omega) + a^2*(-87 + 36*ap + 2*im*(-5 + 17*ap)*omega + 8*(1 + ap)*omega^2) 
                - a^6*(3 + m^2 - 8*omega*(im + 2*omega)) + 2*a^4*(3*(4 + ap) - (-1 + ap)*m^2 + im*omega 
                - 2*omega*(9*im*ap + (5 + 4*ap)*omega)))) * CHp1 + 24*a^4*(-1 + a^2)^2*(4 - 4*ap + a^2*(-3 
                + 4*ap + im*a*(-1 + ap)*m + a^2*(-1 + 4*im*omega) - 2*im*(1 + ap)*omega)) * CHp2 
                + 48*a^6*(-1 + a^2)^3*CHp3 )) * YHp0
    else
        throw(DomainError(order, "Currently horizontal expansion orders only up to 3 are supported"))
    end
end

function Ypp_ingoing_hor_m2(m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    # Compute the outgoing expansion coefficients of Y'' at the horizon

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    ap = sqrt(1 - a^2)
    YHm0 = (2^(-((2*im*a*(1 + ap)*m + a^2*(-1 + 4*im*omega) + (1 + ap)*(1 - 8*im*omega)) / (2*(1 - a^2 + ap)))) 
            * (1 - a^2)^(-((im*a*(m - 2*a*omega)) / (2*(1 - a^2 + ap)))) * exp((im*(-a^3*m + a*(1 + ap)*m 
            - 8*(1 + ap)*omega + 4*a^2*(2 + ap)*omega)) / (2*(1 - a^2 + ap)))) / (1 + ap)^(5/2)

    if order == 0
        return YHm0
    elseif order == 1
        CHm1 = a_minus_hor(-2, m, a, omega, lambda, 1; data_type=data_type)
        return ((1 / (2 * ap * (1 - a^2 + ap))) * (-5 - ap - im*a^3*m - 4*im*omega - 4*im*ap*omega 
                + a^2*(5 + ap + 2*im*(3 + 2*ap)*omega) - 2*(-1 + a^2)*(1 + ap)*CHm1)) * YHm0
    elseif order == 2
        CHm1 = a_minus_hor(-2, m, a, omega, lambda, 1; data_type=data_type)
        CHm2 = a_minus_hor(-2, m, a, omega, lambda, 2; data_type=data_type)
        return ((1 / (8 * a^4 * (-1 + a^2)^2)) * (-48 * (-1 + ap) +a^2 * (-114 + 90 * ap + a * (
                -15 * im * (-1 + ap) * m + 7 * im * a^2 * (-3 + 2 * ap) * m + 2 * a^4 * m * (im * (3 + ap) + (2 - 4 * ap) * omega) +
                a^5 * (3 + m^2 + 8 * (im - 2 * omega) * omega) - a * (-87 + 36 * ap + 2 * im * (5 - 17 * ap) * omega + 8 * (1 + ap) * omega^2) +
                2 * a^3 * (-3 * (4 + ap) + (-1 + ap) * m^2 + omega * (im - 18 * im * ap + 10 * omega + 8 * ap * omega)))) 
                + 4 * a^2 * (-1 + a^2) * (4 - 4 * ap + a^2 * (-3 + 4 * ap - im * a * (-1 + ap) * m + a^2 * (-1 - 4 * im * omega) 
                + 2 * im * (1 + ap) * omega)) * CHm1 + 8 * a^4 * (-1 + a^2)^2 * CHm2)) * YHm0
    elseif order == 3
        CHm1 = a_minus_hor(-2, m, a, omega, lambda, 1; data_type=data_type)
        CHm2 = a_minus_hor(-2, m, a, omega, lambda, 2; data_type=data_type)
        CHm3 = a_minus_hor(-2, m, a, omega, lambda, 3; data_type=data_type)
        return ((1 / (48 * a^6 * (-1 + a^2)^3)) * ( -768*(-1 + ap) + a^2*(96*(-29 + 25*ap) + a*(-264*im*(-1 + ap)*m 
                + im*a^2*(-709 + 577*ap)*m + a^4*m*(619*im - 358*im*ap + 6*(20 - 21*ap)*omega) + 6*a*(619 - 435*ap 
                + 64*im*(-1 + ap)*omega) + a^8*m*(-9*im*(1 + ap) - im*(-3 + ap)*m^2 + 84*omega + 24*ap*(1 + 2*im*omega)*omega) 
                + a^9*(-3*(5 + m^2) - 12*im*(3 + m^2)*omega + 48*omega^2 + 64*im*omega^3) + a^6*m*(3*im*(-55 + 18*ap) 
                + 4*im*(-1 + ap)*m^2 + 6*omega*(-34 + 19*ap - 2*im*(1 + ap)*omega)) + a^5*(345 - 150*ap + (72 - 51*ap)*m^2 
                + 2*im*(-484 + 115*ap)*omega + 4*omega^2*(39 + 75*ap + 2*im*(15 + 13*ap)*omega)) + 3*a^7*(25 + 6*ap 
                + m^2*(-9 + 2*ap - 6*im*(-1 + ap)*omega) - 2*im*omega*(-58 - 15*ap - 2*im*(11 + 20*ap)*omega + 8*(3 
                + 2*ap)*omega^2)) + a^3*(3*(-701 + 370*ap) + 42*(-1 + ap)*m^2 - 8*im*omega*(-130 + 88*ap + (1 + ap)
                *omega*(-9*im + 4*omega))) )) + 6*a^2*(-1 + a^2)*(-48*(-1 + ap) + a^2*(-114 + 90*ap - 15*im*a*(-1 + ap)*m 
                + 7*im*a^3*(-3 + 2*ap)*m + 2*a^5*m*(im*(3 + ap) + (2 - 4*ap)*omega) + a^6*(3 + m^2 + 8*(im - 2*omega)*omega) 
                - a^2*(-87 + 36*ap + 2*im*(5 - 17*ap)*omega + 8*(1 + ap)*omega^2) + 2*a^4*(-3*(4 + ap) + (-1 + ap)*m^2 
                + omega*(im - 18*im*ap + 10*omega + 8*ap*omega)) )) * CHm1 + 4*a^2*(-1 + a^2)*(4 - 4*ap + a^2*(-3 + 4*ap 
                - im*a*(-1 + ap)*m + a^2*(-1 - 4*im*omega) + 2*im*(1 + ap)*omega)) * CHm2 + 2*a^2*(-1 + a^2)*CHm3 )) * YHm0
    else
        throw(DomainError(order, "Currently horizontal expansion orders only up to 3 are supported"))
    end
end

function Ypp_ingoing_hor_p2(m::Int, a, omega, lambda, order::Int; data_type=_DEFAULTDATATYPE)
    # Compute the outgoing expansion coefficients of Y'' at the horizon

    if order < 0
        throw(DomainError(order, "Only positive expansion order is supported"))
    end

    ap = sqrt(1 - a^2)
    YHm0 = exp((im * m * (a^2 - (-1 + ap) * log(1 - a^2))) / (2 * a)) / (sqrt(2) * (1 + ap)^(5/2))

    if order == 0
        return YHm0
    elseif order == 1
        CHm1 = a_minus_hor(2, m, a, omega, lambda, 1; data_type=data_type)
        return (-(1/2) + (2*(-1 + ap))/a^2 + (im*a*m)/(2*ap) + CHm1) * YHm0
    elseif order == 2
        CHm1 = a_minus_hor(2, m, a, omega, lambda, 1; data_type=data_type)
        CHm2 = a_minus_hor(2, m, a, omega, lambda, 2; data_type=data_type)
        return ((1 / (8 * a^4 * (-1 + a^2))) * (48*(-1 + ap) + a^2*(66 - 42*ap - 3*a^2*(7 + 2*ap) 
                + 7*im*a*(-1 + ap)*m + 2*im*a^3*(4 + ap)*m + a^4*(3 + m^2)) - 4*a^2*(4*(-1 + ap) 
                + a^2*(3 + a^2 - 4*ap + im*a*ap*m)) * CHm1 + 8*a^4*(-1 + a^2)*CHm2)) * YHm0
    elseif order == 3
        CHm1 = a_minus_hor(2, m, a, omega, lambda, 1; data_type=data_type)
        CHm2 = a_minus_hor(2, m, a, omega, lambda, 2; data_type=data_type)
        CHm3 = a_minus_hor(2, m, a, omega, lambda, 3; data_type=data_type)
        return ((1 / (48 * a^6 * (-1 + a^2)^2)) * (768*(-1 + ap) + a^2*(96*(21 - 17*ap) + 6*a^2*(-283 + 163*ap) 
                + 120*im*a*(-1 + ap)*m - 5*im*a^3*(-47 + 35*ap)*m + im*a^5*(-97 + 66*ap)*m - 3*a^8*(5 + m^2) 
                - im*a^7*m*(9*(2 + ap) + ap*m^2) + a^4*(405 - 132*ap - 9*(-1 + ap)*m^2) + 6*a^6*(10 + 3*ap 
                + (-1 + 2*ap)*m^2)) + 6*a^2*(-1 + a^2)*(48*(-1 + ap) + a^2*(66 - 42*ap - 3*a^2*(7 + 2*ap) 
                + 7*im*a*(-1 + ap)*m + 2*im*a^3*(4 + ap)*m + a^4*(3 + m^2)) * CHm1 - 4*a^2*(4*(-1 + ap) 
                + a^2*(3 + a^2 - 4*ap + im*a*ap*m)) * CHm2 + 8*a^4*(-1 + a^2)*CHm3))) * YHm0
    else
        throw(DomainError(order, "Currently horizontal expansion orders only up to 3 are supported"))
    end
end

end