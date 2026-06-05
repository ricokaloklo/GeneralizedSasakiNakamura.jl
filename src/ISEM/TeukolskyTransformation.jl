module TeukolskyTransformation

using StaticArrays

export r_to_x, x_to_r, Tx, Ex
export Delta, Vr, Vrp, d2R, d3R, d4R, d5R, coe1_funcs, coe2_funcs
export PlusTwo_to_MinusTwo_Inf, PlusTwo_to_MinusTwo_Hor, MinusTwo_to_PlusTwo_Inf, MinusTwo_to_PlusTwo_Hor
export calculate_q1_q2, d2P, calculate_q3_q4, calculate_q5_q6, _relative_residual
export contour_phase, contour_q1, contour_q2, contour_residual
export evaluate_power_series, evaluate_inverse_power_series, evaluate_frobenius_series
export contour_PlusTwo_to_MinusTwo_Inf, contour_PlusTwo_to_MinusTwo_Hor, contour_MinusTwo_to_PlusTwo_Inf, contour_MinusTwo_to_PlusTwo_Hor
export P_to_Y_coeffs_minus2_in, P_to_Y_coeffs_plus2_up
export GSN_to_Y_coeffs_minus2_in, GSN_to_Y_coeffs_plus2_up, GSN_to_Y_solution_from_matrix
export P_to_GSN_coefficients_from_matrix, P_to_GSN_solution_from_matrix
export Teukolsky_to_GSN_solution

function r_to_x(r, r_plus, r_minus)
    return (r_plus - r) / (r_plus - r_minus)
end

function r_to_x(r, kappa)
    return (1 + kappa - r) / (2 * kappa)
end

function x_to_r(x, r_plus, r_minus)
    return r_plus - x * (r_plus - r_minus)
end

function x_to_r(x, kappa)
    return 1 + kappa - 2 * x * kappa
end

function Delta(r, a)
    return r^2 + a^2 - 2 * r
end

function Teukolsky_to_GSN_solution(r_from_rstar, teukolsky_to_gsn_matrix, Rsoln)
    return rs -> begin
        r = r_from_rstar(rs)
        vals = Rsoln(r)
        teuk_vec = SA[vals[1], vals[2]]
        teukolsky_to_gsn_matrix(r) * teuk_vec
    end
end

@inline _is_complex_frequency(omega) = !(omega isa Real)

# Direct Kerr P -> GSN coefficient formulas generated from ISEM-T111.
# They exclude the common Teukolsky prefactor T(x).
@inline _direct_power(x, y) = x^y



function _P_to_GSN_direct_coefficients_kerr(s::Int, m::Int, a, omega, lambda, r)

    af = Float64(a)

    omegaf = Float64(omega)

    lambdaf = Float64(lambda)

    rf = Float64(r)

    if s == 0

        return _P_to_GSN_direct_coeffs_s0(af, omegaf, m, lambdaf, rf)

    elseif s == 1

        return _P_to_GSN_direct_coeffs_s1(af, omegaf, m, lambdaf, rf)

    elseif s == -1

        return _P_to_GSN_direct_coeffs_sminus1(af, omegaf, m, lambdaf, rf)

    elseif s == 2

        return _P_to_GSN_direct_coeffs_s2(af, omegaf, m, lambdaf, rf)

    elseif s == -2

        return _P_to_GSN_direct_coeffs_sminus2(af, omegaf, m, lambdaf, rf)

    else

        error("unsupported spin")

    end

end



function _P_to_GSN_direct_coeffs_s0(a, omega, m, lambda, r)
    A0 = sqrt(_direct_power(a,2) + _direct_power(r,2))
    A1 = 0
    B0 = ((_direct_power(a,2) + (-2 + r)*r)*(2*r - ((_direct_power(a,2) + _direct_power(r,2))*(complex(0,-4)*omega + complex(0,2)*sqrt(1 - _direct_power(a,2))*omega + (complex(0,2)*(-1 + r)*(-(a*m) + 2*omega*r))/(_direct_power(a,2) + (-2 + r)*r)))/sqrt(1 - _direct_power(a,2))))/(2.0*_direct_power(_direct_power(a,2) + _direct_power(r,2),1.5))
    B1 = -0.5*(_direct_power(a,2) + (-2 + r)*r)/sqrt(-((-1 + _direct_power(a,2))*(_direct_power(a,2) + _direct_power(r,2))))
    return A0, A1, B0, B1
end

function _P_to_GSN_direct_coeffs_s1(a, omega, m, lambda, r)
    A0 = (sqrt((_direct_power(a,2) + _direct_power(r,2))/(1 - _direct_power(a,2)))*(_direct_power(a,4)*(1 + complex(0,2)*omega) - complex(0,1)*_direct_power(a,3)*m*(1 + sqrt(1 - _direct_power(a,2)) - r) + complex(0,1)*a*m*_direct_power(r,2)*(-1 - sqrt(1 - _direct_power(a,2)) + r) + (-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(r,2)*(1 + complex(0,2)*omega*r) + _direct_power(a,2)*(-1 - sqrt(1 - _direct_power(a,2)) + _direct_power(r,2) + complex(0,2)*omega*r*(-1 + sqrt(1 - _direct_power(a,2)) + r))))/_direct_power(r,2)
    A1 = -0.5*((_direct_power(a,2) + (-2 + r)*r)*_direct_power(_direct_power(a,2) + _direct_power(r,2),1.5))/(sqrt(1 - _direct_power(a,2))*_direct_power(r,2))
    B0 = (-(_direct_power(a,8)*(complex(0,-1) + 2*omega)*(complex(0,2) + omega*r)) + _direct_power(a,7)*m*(complex(0,2)*(1 + sqrt(1 - _direct_power(a,2))) + (complex(0,-3) + (3 + sqrt(1 - _direct_power(a,2)))*omega)*r - omega*_direct_power(r,2)) + _direct_power(r,6)*(-4*sqrt(1 - _direct_power(a,2)) + sqrt(1 - _direct_power(a,2))*lambda*(-2 + r) + 2*sqrt(1 - _direct_power(a,2))*r + complex(0,1)*(-1 + sqrt(1 - _direct_power(a,2)))*omega*r - 2*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2)*_direct_power(r,2)) + _direct_power(a,3)*m*_direct_power(r,3)*(complex(0,-2) + (complex(0,-2)*(-3 + sqrt(1 - _direct_power(a,2))) + 4*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (complex(0,-3) + (5 + 3*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) - 3*omega*_direct_power(r,3)) + _direct_power(a,5)*m*r*(complex(0,-3)*(1 + sqrt(1 - _direct_power(a,2))) + 2*(complex(0,1)*(4 + sqrt(1 - _direct_power(a,2))) + (-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (complex(0,-6) + (7 + 3*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) - 3*omega*_direct_power(r,3)) + a*m*_direct_power(r,5)*(complex(0,1)*(1 + 3*sqrt(1 - _direct_power(a,2))) - complex(0,2)*sqrt(1 - _direct_power(a,2))*r + 2*(-1 + sqrt(1 - _direct_power(a,2)))*omega*r + (1 + sqrt(1 - _direct_power(a,2)))*omega*_direct_power(r,2) - omega*_direct_power(r,3)) + _direct_power(a,6)*(2*(1 + sqrt(1 - _direct_power(a,2))) + (4 + sqrt(1 - _direct_power(a,2)) + sqrt(1 - _direct_power(a,2))*lambda - (1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + complex(0,11)*omega - complex(0,5)*sqrt(1 - _direct_power(a,2))*omega)*r + (_direct_power(m,2) - 2*(2 + complex(0,4)*omega + (-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2)))*_direct_power(r,2) + 3*(complex(0,1) - 2*omega)*omega*_direct_power(r,3)) - _direct_power(a,2)*_direct_power(r,3)*(4 + (-2 + 6*sqrt(1 - _direct_power(a,2)) + 4*sqrt(1 - _direct_power(a,2))*lambda - complex(0,8)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (-5*sqrt(1 - _direct_power(a,2)) - 3*sqrt(1 - _direct_power(a,2))*lambda + (1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) - complex(0,1)*omega + complex(0,3)*sqrt(1 - _direct_power(a,2))*omega)*_direct_power(r,2) - (_direct_power(m,2) - 6*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) + omega*(complex(0,-1) + 2*omega)*_direct_power(r,4)) - _direct_power(a,4)*r*(4*(1 + sqrt(1 - _direct_power(a,2))) + 2*(-2 + sqrt(1 - _direct_power(a,2))*lambda - complex(0,4)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (-4 - 4*sqrt(1 - _direct_power(a,2)) - 3*sqrt(1 - _direct_power(a,2))*lambda + 2*(1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) - complex(0,13)*omega + complex(0,9)*sqrt(1 - _direct_power(a,2))*omega)*_direct_power(r,2) + (2 - 2*_direct_power(m,2) + complex(0,4)*omega + 6*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) + 3*omega*(complex(0,-1) + 2*omega)*_direct_power(r,4)))/(sqrt(1 - _direct_power(a,2))*_direct_power(r,3)*_direct_power(_direct_power(a,2) + _direct_power(r,2),1.5))
    B1 = (complex(0,-0.5)*(_direct_power(a,2) + (-2 + r)*r)*(-(_direct_power(a,3)*m*r) - a*m*_direct_power(r,3) + omega*_direct_power(r,5) + _direct_power(a,4)*(complex(0,2) + omega*r) + 2*_direct_power(a,2)*r*(complex(0,-2) + r*(complex(0,1) + omega*r))))/(_direct_power(r,3)*sqrt(-((-1 + _direct_power(a,2))*(_direct_power(a,2) + _direct_power(r,2)))))
    return A0, A1, B0, B1
end

function _P_to_GSN_direct_coeffs_sminus1(a, omega, m, lambda, r)
    A0 = (sqrt((_direct_power(a,2) + _direct_power(r,2))/(1 - _direct_power(a,2)))*(_direct_power(a,4)*(-1 - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*omega) + complex(0,1)*_direct_power(a,3)*m*(-1 + sqrt(1 - _direct_power(a,2)) + r) + complex(0,1)*a*m*_direct_power(r,2)*(-1 + sqrt(1 - _direct_power(a,2)) + r) + _direct_power(r,2)*(1 + sqrt(1 - _direct_power(a,2)) - complex(0,2)*omega*r*(1 - sqrt(1 - _direct_power(a,2)) + sqrt(1 - _direct_power(a,2))*r)) + _direct_power(a,2)*(1 - sqrt(1 - _direct_power(a,2)) - _direct_power(r,2) - complex(0,2)*omega*r*(1 - sqrt(1 - _direct_power(a,2)) - r + 2*sqrt(1 - _direct_power(a,2))*r))))/(_direct_power(r,2)*(_direct_power(a,2) + (-2 + r)*r))
    A1 = -0.5*_direct_power(_direct_power(a,2) + _direct_power(r,2),1.5)/(sqrt(1 - _direct_power(a,2))*_direct_power(r,2))
    B0 = (-(_direct_power(a,8)*(complex(0,-1) + 2*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*(complex(0,-2) + omega*r)) + _direct_power(a,7)*m*(complex(0,-2)*(-1 + sqrt(1 - _direct_power(a,2))) + 3*(complex(0,-1) + (-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + omega*_direct_power(r,2)) + a*m*_direct_power(r,5)*(complex(0,1) - complex(0,3)*sqrt(1 - _direct_power(a,2)) + complex(0,2)*sqrt(1 - _direct_power(a,2))*r - 2*(-1 + sqrt(1 - _direct_power(a,2)))*omega*r + (-1 + 3*sqrt(1 - _direct_power(a,2)))*omega*_direct_power(r,2) + omega*_direct_power(r,3)) + _direct_power(a,5)*m*r*(complex(0,3)*(-1 + sqrt(1 - _direct_power(a,2))) + (complex(0,-2)*(-4 + sqrt(1 - _direct_power(a,2))) - 2*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (complex(0,-6) + (-7 + 9*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) + 3*omega*_direct_power(r,3)) + _direct_power(a,3)*m*_direct_power(r,3)*(complex(0,-2) + (complex(0,2)*(3 + sqrt(1 - _direct_power(a,2))) - 4*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (complex(0,-3) + (-5 + 9*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) + 3*omega*_direct_power(r,3)) + _direct_power(a,6)*(2*(-1 + sqrt(1 - _direct_power(a,2))) + (-4 - sqrt(1 - _direct_power(a,2)) + sqrt(1 - _direct_power(a,2))*lambda - (-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + complex(0,11)*omega - complex(0,11)*sqrt(1 - _direct_power(a,2))*omega)*r - (-4 + _direct_power(m,2) - complex(0,4)*(-2 + 3*sqrt(1 - _direct_power(a,2)))*omega - 2*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,2) + omega*(complex(0,3) + (6 - 8*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,3)) + _direct_power(a,2)*_direct_power(r,3)*(4 + (-4*sqrt(1 - _direct_power(a,2))*lambda + 2*(-1 + sqrt(1 - _direct_power(a,2)))*(1 + complex(0,4)*omega))*r + (-sqrt(1 - _direct_power(a,2)) + 3*sqrt(1 - _direct_power(a,2))*lambda + _direct_power(m,2) - sqrt(1 - _direct_power(a,2))*_direct_power(m,2) + complex(0,1)*omega - complex(0,13)*sqrt(1 - _direct_power(a,2))*omega)*_direct_power(r,2) - (_direct_power(m,2) - 2*omega*(complex(0,2)*sqrt(1 - _direct_power(a,2)) + 3*(-1 + sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,3) + omega*(complex(0,1) + (2 - 8*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,4)) + _direct_power(a,4)*r*(4 - 4*sqrt(1 - _direct_power(a,2)) - 2*sqrt(1 - _direct_power(a,2))*lambda*r + 4*(-1 + sqrt(1 - _direct_power(a,2)))*(1 + complex(0,2)*omega)*r + (-4 - 2*sqrt(1 - _direct_power(a,2)) + 3*sqrt(1 - _direct_power(a,2))*lambda - 2*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + complex(0,13)*omega - complex(0,23)*sqrt(1 - _direct_power(a,2))*omega)*_direct_power(r,2) - 2*(-1 + _direct_power(m,2) - complex(0,2)*(-1 + 3*sqrt(1 - _direct_power(a,2)))*omega - 3*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) + 3*omega*(complex(0,1) + (2 - 4*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,4)) + _direct_power(r,6)*(sqrt(1 - _direct_power(a,2))*lambda*(-2 + r) - omega*r*(complex(0,1)*(1 + sqrt(1 - _direct_power(a,2))) + 2*omega*r*(1 - sqrt(1 - _direct_power(a,2)) + sqrt(1 - _direct_power(a,2))*r))))/(_direct_power(r,3)*(_direct_power(a,2) + (-2 + r)*r)*sqrt(-((-1 + _direct_power(a,2))*_direct_power(_direct_power(a,2) + _direct_power(r,2),3))))
    B1 = (complex(0,0.5)*(-(_direct_power(a,3)*m*r) - a*m*_direct_power(r,3) + omega*_direct_power(r,5) + _direct_power(a,4)*(complex(0,-2) + omega*r) + 2*_direct_power(a,2)*r*(complex(0,2) + r*(complex(0,-1) + omega*r))))/(_direct_power(r,3)*sqrt(-((-1 + _direct_power(a,2))*(_direct_power(a,2) + _direct_power(r,2)))))
    return A0, A1, B0, B1
end

function _P_to_GSN_direct_coeffs_s2(a, omega, m, lambda, r)
    A0 = (sqrt(-((-1 + _direct_power(a,2))*(_direct_power(a,2) + _direct_power(r,2))))*(_direct_power(a,2)*r*(-8 + 16*sqrt(1 - _direct_power(a,2)) - (sqrt(1 - _direct_power(a,2))*(lambda - 2*(_direct_power(m,2) + complex(0,6)*omega)) + 2*(6 + 4*sqrt(1 - _direct_power(a,2)) - _direct_power(m,2) + complex(0,8)*omega))*r + (4 - 2*_direct_power(m,2) + complex(0,2)*(2 + sqrt(1 - _direct_power(a,2)))*omega + 4*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,2) + 4*omega*(complex(0,-1) + omega)*_direct_power(r,3)) + 2*_direct_power(a,4)*(-3*sqrt(1 - _direct_power(a,2)) + 2*(complex(0,-1) + omega)*r*(complex(0,2) + omega*r)) + 2*a*m*_direct_power(r,2)*(complex(0,1)*(1 + sqrt(1 - _direct_power(a,2))) + complex(0,1)*r*(-4 + sqrt(1 - _direct_power(a,2)) + r) + omega*r*(2 - 2*sqrt(1 - _direct_power(a,2)) - (1 + sqrt(1 - _direct_power(a,2)))*r + _direct_power(r,2))) + 2*_direct_power(a,3)*m*r*(complex(0,-2)*(1 + sqrt(1 - _direct_power(a,2))) + r*(complex(0,4) + omega*(-3 - sqrt(1 - _direct_power(a,2)) + r))) + _direct_power(r,2)*(-12*(-1 + sqrt(1 - _direct_power(a,2))) + r*(4*(-1 + 4*sqrt(1 - _direct_power(a,2)) + complex(0,3)*omega) - sqrt(1 - _direct_power(a,2))*lambda*(-2 + r) - 6*sqrt(1 - _direct_power(a,2))*r + 4*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2)*_direct_power(r,2) + complex(0,2)*sqrt(1 - _direct_power(a,2))*omega*(-6 + (-2 + r)*r)))))/((-1 + _direct_power(a,2))*_direct_power(r,2))
    A1 = ((_direct_power(a,2) + (-2 + r)*r)*sqrt((_direct_power(a,2) + _direct_power(r,2))/(1 - _direct_power(a,2)))*(2*_direct_power(a,2) - 3*r + complex(0,1)*a*(m - a*omega)*r + _direct_power(r,2) - complex(0,1)*omega*_direct_power(r,3)))/r
    B0 = (2*_direct_power(a,7)*m*r*(complex(0,-1)*(5 + 2*sqrt(1 - _direct_power(a,2))) + (complex(0,9) - 2*(3 + sqrt(1 - _direct_power(a,2)))*omega)*r + omega*(6 + complex(0,1)*(5 + sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) - complex(0,1)*_direct_power(omega,2)*_direct_power(r,3)) + 2*_direct_power(a,8)*(-6*sqrt(1 - _direct_power(a,2)) + (10 - complex(0,1)*(-10 + 3*sqrt(1 - _direct_power(a,2)))*omega)*r + 4*omega*(complex(0,-1) + omega)*_direct_power(r,2) - complex(0,2)*_direct_power(omega,2)*(complex(0,-1) + omega)*_direct_power(r,3)) + _direct_power(a,5)*m*_direct_power(r,2)*(complex(0,24) + (complex(0,-68) - complex(0,2)*sqrt(1 - _direct_power(a,2)) - complex(0,1)*(1 + 2*sqrt(1 - _direct_power(a,2)))*lambda + complex(0,2)*(1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + 18*omega - 10*sqrt(1 - _direct_power(a,2))*omega)*r + complex(0,1)*(32 + lambda - 2*_direct_power(m,2) + complex(0,8)*(3 + sqrt(1 - _direct_power(a,2)))*omega + 8*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,2) + 2*omega*(13 + complex(0,1)*(11 + 3*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,3) - complex(0,6)*_direct_power(omega,2)*_direct_power(r,4)) + _direct_power(a,6)*r*(-20 + 52*sqrt(1 - _direct_power(a,2)) + (-64 - 38*sqrt(1 - _direct_power(a,2)) - 4*sqrt(1 - _direct_power(a,2))*lambda + 4*(1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) - complex(0,76)*omega + complex(0,52)*sqrt(1 - _direct_power(a,2))*omega)*r + (48 + 2*lambda - complex(0,10)*(-6 + sqrt(1 - _direct_power(a,2)))*omega + complex(0,1)*(2 + sqrt(1 - _direct_power(a,2)))*lambda*omega + 4*(-4 + 3*sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2) - complex(0,4)*_direct_power(m,2)*(complex(0,-2) + (2 + sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,2) + 2*omega*(complex(0,-10) + complex(0,2)*_direct_power(m,2) + (10 + sqrt(1 - _direct_power(a,2)))*omega - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) - complex(0,12)*_direct_power(omega,2)*(complex(0,-1) + omega)*_direct_power(r,4)) + a*m*_direct_power(r,5)*(complex(0,8)*(-1 + 2*sqrt(1 - _direct_power(a,2))) + (complex(0,-2)*(-15 + sqrt(1 - _direct_power(a,2))) + complex(0,2)*(1 + 2*sqrt(1 - _direct_power(a,2)))*lambda + 4*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (complex(0,-4)*(6 + sqrt(1 - _direct_power(a,2))) - complex(0,1)*(3 + 2*sqrt(1 - _direct_power(a,2)))*lambda + 2*(-1 + 5*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) + complex(0,1)*(6 + lambda + complex(0,8)*omega + 8*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) + 2*omega*(1 + complex(0,1)*(1 + sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,4) - complex(0,2)*_direct_power(omega,2)*_direct_power(r,5)) + 2*_direct_power(a,3)*m*_direct_power(r,3)*(complex(0,6)*(-1 + sqrt(1 - _direct_power(a,2))) + (complex(0,31) - complex(0,7)*sqrt(1 - _direct_power(a,2)) + complex(0,1)*(1 + 2*sqrt(1 - _direct_power(a,2)))*lambda + 6*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (complex(0,-35) - complex(0,1)*sqrt(1 - _direct_power(a,2)) - complex(0,2)*(1 + sqrt(1 - _direct_power(a,2)))*lambda + complex(0,1)*(1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + 4*omega)*_direct_power(r,2) + complex(0,1)*(10 + lambda - _direct_power(m,2) + complex(0,10)*omega + complex(0,2)*sqrt(1 - _direct_power(a,2))*omega - 8*_direct_power(omega,2) + 8*sqrt(1 - _direct_power(a,2))*_direct_power(omega,2))*_direct_power(r,3) + omega*(8 + complex(0,1)*(7 + 3*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,4) - complex(0,3)*_direct_power(omega,2)*_direct_power(r,5)) + _direct_power(a,4)*_direct_power(r,2)*(64 - 80*sqrt(1 - _direct_power(a,2)) + 2*(66*sqrt(1 - _direct_power(a,2)) + (-1 + 7*sqrt(1 - _direct_power(a,2)))*lambda - (1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + complex(0,50)*omega - complex(0,50)*sqrt(1 - _direct_power(a,2))*omega)*r - (92 + 50*sqrt(1 - _direct_power(a,2)) + lambda*(4 + 9*sqrt(1 - _direct_power(a,2)) + complex(0,6)*omega) - complex(0,24)*(-5 + 3*sqrt(1 - _direct_power(a,2)))*omega + 12*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2) + complex(0,2)*_direct_power(m,2)*(complex(0,5) + 2*(-1 + sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,2) + (40 + 4*lambda + complex(0,64)*omega + complex(0,4)*sqrt(1 - _direct_power(a,2))*omega + complex(0,1)*(4 + 3*sqrt(1 - _direct_power(a,2)))*lambda*omega - 32*_direct_power(omega,2) + 20*sqrt(1 - _direct_power(a,2))*_direct_power(omega,2) - complex(0,4)*_direct_power(m,2)*(complex(0,-2) + (3 + 2*sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,3) + 2*omega*(complex(0,-8) + complex(0,4)*_direct_power(m,2) + (8 + 3*sqrt(1 - _direct_power(a,2)))*omega - complex(0,6)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,4) - complex(0,12)*_direct_power(omega,2)*(complex(0,-1) + omega)*_direct_power(r,5)) + _direct_power(r,5)*(24*(-1 + sqrt(1 - _direct_power(a,2))) - 4*(-9 + 15*sqrt(1 - _direct_power(a,2)) + (-1 + 2*sqrt(1 - _direct_power(a,2)))*lambda - complex(0,6)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + 2*(6 + lambda)*(-1 + 3*sqrt(1 - _direct_power(a,2)) - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) - (6*sqrt(1 - _direct_power(a,2)) + lambda*(sqrt(1 - _direct_power(a,2)) + complex(0,2)*omega) + complex(0,4)*(2 + sqrt(1 - _direct_power(a,2)))*omega + 12*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) + complex(0,1)*sqrt(1 - _direct_power(a,2))*(6 + lambda + complex(0,4)*omega)*omega*_direct_power(r,4) + 2*_direct_power(omega,2)*(sqrt(1 - _direct_power(a,2)) - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,5)) + _direct_power(a,2)*_direct_power(r,3)*(48*(-1 + sqrt(1 - _direct_power(a,2))) - 4*(-23 + 35*sqrt(1 - _direct_power(a,2)) + (-1 + 3*sqrt(1 - _direct_power(a,2)))*lambda - complex(0,12)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + 2*(-8 + 62*sqrt(1 - _direct_power(a,2)) + (1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + complex(0,46)*omega - complex(0,46)*sqrt(1 - _direct_power(a,2))*omega + 2*lambda*(-1 + 5*sqrt(1 - _direct_power(a,2)) - complex(0,1)*(-1 + sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,2) - 2*(lambda*(2 + 3*sqrt(1 - _direct_power(a,2)) + complex(0,4)*omega) + _direct_power(m,2)*(-1 + 2*sqrt(1 - _direct_power(a,2)) + complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*omega) + (complex(0,-3) + 2*omega)*(complex(0,1)*(6 + 5*sqrt(1 - _direct_power(a,2))) + 6*(-1 + sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,3) + (12 - complex(0,2)*(-12 - 7*sqrt(1 - _direct_power(a,2)) + 2*(1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2))*omega + 4*(-4 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2) + lambda*(2 + complex(0,1)*(2 + 3*sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,4) + 2*omega*(complex(0,-2) + complex(0,2)*_direct_power(m,2) + (2 + 3*sqrt(1 - _direct_power(a,2)))*omega - complex(0,6)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,5) - complex(0,4)*_direct_power(omega,2)*(complex(0,-1) + omega)*_direct_power(r,6)))/(sqrt(1 - _direct_power(a,2))*_direct_power(r,3)*_direct_power(_direct_power(a,2) + _direct_power(r,2),1.5))
    B1 = ((_direct_power(a,2) + (-2 + r)*r)*(-10*_direct_power(a,6) + 4*_direct_power(a,4)*(8 + complex(0,1)*a*(-m + a*omega))*r + _direct_power(a,2)*(-24 + a*(complex(0,6)*m + a*(-24 - lambda - complex(0,6)*omega + 2*_direct_power(m - a*omega,2))))*_direct_power(r,2) + 2*_direct_power(a,2)*(23 + lambda + complex(0,1)*a*(-2*m + 5*a*omega))*_direct_power(r,3) + 2*(-6 + a*(complex(0,1)*m + a*(-10 - lambda - complex(0,6)*omega + (m - 3*a*omega)*(m - a*omega))))*_direct_power(r,4) + 2*(9 + lambda + complex(0,4)*_direct_power(a,2)*omega)*_direct_power(r,5) - (6 + lambda + 2*omega*(complex(0,3) + 2*a*m - 3*_direct_power(a,2)*omega))*_direct_power(r,6) + complex(0,2)*omega*_direct_power(r,7) + 2*_direct_power(omega,2)*_direct_power(r,8)))/(2.0*sqrt(1 - _direct_power(a,2))*_direct_power(r,2)*_direct_power(_direct_power(a,2) + _direct_power(r,2),1.5))
    return A0, A1, B0, B1
end

function _P_to_GSN_direct_coeffs_sminus2(a, omega, m, lambda, r)
    A0 = (sqrt((_direct_power(a,2) + _direct_power(r,2))/(1 - _direct_power(a,2)))*(_direct_power(a,2)*r*(-8*(1 + 2*sqrt(1 - _direct_power(a,2))) + (sqrt(1 - _direct_power(a,2))*(lambda - 2*(_direct_power(m,2) + complex(0,8)*omega)) + 2*(-6 + 2*sqrt(1 - _direct_power(a,2)) + _direct_power(m,2) + complex(0,8)*omega))*r - 2*(-2 + _direct_power(m,2) + omega*(complex(0,2) - complex(0,7)*sqrt(1 - _direct_power(a,2)) + 2*omega - 2*sqrt(1 - _direct_power(a,2))*omega))*_direct_power(r,2) + 4*omega*(complex(0,1) + omega - 2*sqrt(1 - _direct_power(a,2))*omega)*_direct_power(r,3)) + 2*_direct_power(a,4)*(3*sqrt(1 - _direct_power(a,2)) - 2*(complex(0,-1) + (-1 + sqrt(1 - _direct_power(a,2)))*omega)*r*(complex(0,-2) + omega*r)) + 2*_direct_power(a,3)*m*r*(complex(0,-2)*(-1 + sqrt(1 - _direct_power(a,2))) + r*(complex(0,-4) + omega*(-3 + 3*sqrt(1 - _direct_power(a,2)) + r))) + 2*a*m*_direct_power(r,2)*(complex(0,1)*(-1 + sqrt(1 - _direct_power(a,2))) + complex(0,1)*(4 + sqrt(1 - _direct_power(a,2)) - r)*r + omega*r*(2 - 2*sqrt(1 - _direct_power(a,2)) + r*(-1 + 3*sqrt(1 - _direct_power(a,2)) + r))) + _direct_power(r,2)*(12*(1 + sqrt(1 - _direct_power(a,2))) + r*(-4*(1 + 2*sqrt(1 - _direct_power(a,2)) + complex(0,3)*omega) + sqrt(1 - _direct_power(a,2))*lambda*(-2 + r) + 2*sqrt(1 - _direct_power(a,2))*r + 4*_direct_power(omega,2)*_direct_power(r,2)*(-1 + sqrt(1 - _direct_power(a,2)) - sqrt(1 - _direct_power(a,2))*r) + complex(0,6)*sqrt(1 - _direct_power(a,2))*omega*(2 + (-4 + r)*r)))))/(_direct_power(r,2)*_direct_power(_direct_power(a,2) + (-2 + r)*r,2))
    A1 = (sqrt((_direct_power(a,2) + _direct_power(r,2))/(1 - _direct_power(a,2)))*(complex(0,-1)*a*m*r + _direct_power(a,2)*(2 + complex(0,1)*omega*r) + r*(-3 + r + complex(0,1)*omega*_direct_power(r,2))))/(r*(_direct_power(a,2) + (-2 + r)*r))
    B0 = (2*_direct_power(a,7)*m*r*(complex(0,1)*(-5 + 2*sqrt(1 - _direct_power(a,2))) + (complex(0,9) - 6*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + omega*(-6 - complex(0,5)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) - complex(0,1)*_direct_power(omega,2)*_direct_power(r,3)) + 2*_direct_power(a,8)*(-6*sqrt(1 - _direct_power(a,2)) + (-10 - complex(0,1)*(-10 + 7*sqrt(1 - _direct_power(a,2)))*omega)*r + 4*omega*(complex(0,-1) + (-1 + sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,2) + 2*_direct_power(omega,2)*(1 + complex(0,1)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,3)) + _direct_power(a,5)*m*_direct_power(r,2)*(complex(0,24) + (complex(0,-64) - complex(0,6)*sqrt(1 - _direct_power(a,2)) + complex(0,1)*(-1 + 2*sqrt(1 - _direct_power(a,2)))*lambda - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) - 18*omega + 18*sqrt(1 - _direct_power(a,2))*omega)*r + (complex(0,28) + complex(0,1)*lambda - complex(0,2)*_direct_power(m,2) - 24*(-1 + sqrt(1 - _direct_power(a,2)))*omega + complex(0,8)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,2) + 2*omega*(-13 - complex(0,1)*(-11 + 15*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,3) - complex(0,6)*_direct_power(omega,2)*_direct_power(r,4)) + _direct_power(a,6)*r*(20 + 52*sqrt(1 - _direct_power(a,2)) + (64 - 22*sqrt(1 - _direct_power(a,2)) - 4*sqrt(1 - _direct_power(a,2))*lambda + 4*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) - complex(0,76)*omega + complex(0,52)*sqrt(1 - _direct_power(a,2))*omega)*r + (-40 - 2*lambda - complex(0,2)*(-26 + 23*sqrt(1 - _direct_power(a,2)))*omega - complex(0,1)*(-2 + 3*sqrt(1 - _direct_power(a,2)))*lambda*omega - 16*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2) + 2*_direct_power(m,2)*(4 + complex(0,4)*(-1 + sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,2) + 2*omega*(complex(0,-10) + complex(0,2)*_direct_power(m,2) + 5*(-2 + 3*sqrt(1 - _direct_power(a,2)))*omega - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) + 4*_direct_power(omega,2)*(3 + complex(0,1)*(-3 + 4*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,4)) + a*m*_direct_power(r,5)*(complex(0,-8)*(1 + 2*sqrt(1 - _direct_power(a,2))) + 2*(complex(0,11) + complex(0,9)*sqrt(1 - _direct_power(a,2)) - complex(0,1)*(-1 + 2*sqrt(1 - _direct_power(a,2)))*lambda + 2*omega - 2*sqrt(1 - _direct_power(a,2))*omega)*r + (complex(0,1)*(-3 + 2*sqrt(1 - _direct_power(a,2)))*lambda + 2*(complex(0,-2)*(3 + sqrt(1 - _direct_power(a,2))) + omega + 7*sqrt(1 - _direct_power(a,2))*omega))*_direct_power(r,2) + (complex(0,2) + complex(0,1)*lambda + 8*omega + complex(0,8)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) + 2*omega*(-1 - complex(0,1)*(-1 + 5*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,4) - complex(0,2)*_direct_power(omega,2)*_direct_power(r,5)) + 2*_direct_power(a,3)*m*_direct_power(r,3)*(complex(0,-6)*(1 + sqrt(1 - _direct_power(a,2))) + (complex(0,3)*(9 + 5*sqrt(1 - _direct_power(a,2))) - complex(0,1)*(-1 + 2*sqrt(1 - _direct_power(a,2)))*lambda - 6*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + (complex(0,-27) - complex(0,7)*sqrt(1 - _direct_power(a,2)) + complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*lambda - complex(0,1)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) - 4*omega + 16*sqrt(1 - _direct_power(a,2))*omega)*_direct_power(r,2) + (complex(0,6) + complex(0,1)*lambda - complex(0,1)*_direct_power(m,2) + 10*omega - 6*sqrt(1 - _direct_power(a,2))*omega - complex(0,8)*_direct_power(omega,2) + complex(0,8)*sqrt(1 - _direct_power(a,2))*_direct_power(omega,2))*_direct_power(r,3) + omega*(-8 - complex(0,1)*(-7 + 15*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,4) - complex(0,3)*_direct_power(omega,2)*_direct_power(r,5)) + _direct_power(a,4)*_direct_power(r,2)*(-16*(4 + 5*sqrt(1 - _direct_power(a,2))) + 2*(-4 + 38*sqrt(1 - _direct_power(a,2)) + lambda + 7*sqrt(1 - _direct_power(a,2))*lambda - (-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + complex(0,50)*omega - complex(0,38)*sqrt(1 - _direct_power(a,2))*omega)*r + (lambda*(4 - 9*sqrt(1 - _direct_power(a,2)) + complex(0,2)*(-3 + 4*sqrt(1 - _direct_power(a,2)))*omega) + 2*(38 - 7*sqrt(1 - _direct_power(a,2)) + complex(0,2)*(-24 + 37*sqrt(1 - _direct_power(a,2)))*omega + 6*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2) + _direct_power(m,2)*(-5 - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)))*_direct_power(r,2) + (4*_direct_power(m,2)*(2 + complex(0,1)*(-3 + 4*sqrt(1 - _direct_power(a,2)))*omega) + lambda*(-4 - complex(0,1)*(-4 + 9*sqrt(1 - _direct_power(a,2)))*omega) - 8*(3 + complex(0,1)*(-6 + 7*sqrt(1 - _direct_power(a,2)))*omega + (-4 + 7*sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2)))*_direct_power(r,3) + 2*omega*(complex(0,-8) + complex(0,4)*_direct_power(m,2) + (-8 + 21*sqrt(1 - _direct_power(a,2)))*omega - complex(0,6)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,4) + 12*_direct_power(omega,2)*(1 + complex(0,1)*(-1 + 2*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,5)) + _direct_power(r,5)*(24*(1 + sqrt(1 - _direct_power(a,2))) - 4*(5 + 7*sqrt(1 - _direct_power(a,2)) + lambda + 2*sqrt(1 - _direct_power(a,2))*lambda - complex(0,6)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + 2*(2 + 6*sqrt(1 - _direct_power(a,2)) + lambda + 3*sqrt(1 - _direct_power(a,2))*lambda - complex(0,4)*(-1 + 7*sqrt(1 - _direct_power(a,2)))*omega - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*lambda*omega)*_direct_power(r,2) + (-2*sqrt(1 - _direct_power(a,2)) - sqrt(1 - _direct_power(a,2))*lambda + complex(0,32)*sqrt(1 - _direct_power(a,2))*omega + complex(0,2)*(-1 + 4*sqrt(1 - _direct_power(a,2)))*lambda*omega + 12*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,3) - 3*sqrt(1 - _direct_power(a,2))*omega*(complex(0,2) + complex(0,1)*lambda + 8*omega)*_direct_power(r,4) + 2*_direct_power(omega,2)*(3*sqrt(1 - _direct_power(a,2)) - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,5) + complex(0,4)*sqrt(1 - _direct_power(a,2))*_direct_power(omega,3)*_direct_power(r,6)) + _direct_power(a,2)*_direct_power(r,3)*(48*(1 + sqrt(1 - _direct_power(a,2))) - 4*(19 + 23*sqrt(1 - _direct_power(a,2)) + lambda + 3*sqrt(1 - _direct_power(a,2))*lambda - complex(0,12)*(-1 + sqrt(1 - _direct_power(a,2)))*omega)*r + 2*(22*sqrt(1 - _direct_power(a,2)) + (-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(m,2) + complex(0,38)*omega - complex(0,74)*sqrt(1 - _direct_power(a,2))*omega + 2*lambda*(1 + 5*sqrt(1 - _direct_power(a,2)) - complex(0,1)*(-1 + sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,2) + 2*(10 - 3*sqrt(1 - _direct_power(a,2)) - complex(0,14)*omega + complex(0,64)*sqrt(1 - _direct_power(a,2))*omega - 12*_direct_power(omega,2) + 12*sqrt(1 - _direct_power(a,2))*_direct_power(omega,2) + _direct_power(m,2)*(-1 - 2*sqrt(1 - _direct_power(a,2)) - complex(0,2)*(-1 + sqrt(1 - _direct_power(a,2)))*omega) + lambda*(2 - 3*sqrt(1 - _direct_power(a,2)) + complex(0,4)*(-1 + 2*sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,3) + (-4 + complex(0,2)*(8 - 15*sqrt(1 - _direct_power(a,2)) + (-2 + 4*sqrt(1 - _direct_power(a,2)))*_direct_power(m,2))*omega - 16*(-1 + 4*sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2) + lambda*(-2 - complex(0,1)*(-2 + 9*sqrt(1 - _direct_power(a,2)))*omega))*_direct_power(r,4) + 2*omega*(complex(0,-2) + complex(0,2)*_direct_power(m,2) + (-2 + 13*sqrt(1 - _direct_power(a,2)))*omega - complex(0,6)*(-1 + sqrt(1 - _direct_power(a,2)))*_direct_power(omega,2))*_direct_power(r,5) + 4*_direct_power(omega,2)*(1 + complex(0,1)*(-1 + 4*sqrt(1 - _direct_power(a,2)))*omega)*_direct_power(r,6)))/(sqrt(1 - _direct_power(a,2))*_direct_power(r,3)*_direct_power(_direct_power(a,2) + (-2 + r)*r,2)*_direct_power(_direct_power(a,2) + _direct_power(r,2),1.5))
    B1 = (-10*_direct_power(a,6) + 4*_direct_power(a,4)*(8 + complex(0,1)*a*(m - a*omega))*r + _direct_power(a,2)*(-24 + a*(complex(0,-6)*m + a*(-20 - lambda + complex(0,6)*omega + 2*_direct_power(m - a*omega,2))))*_direct_power(r,2) + 2*_direct_power(a,2)*(19 + lambda + complex(0,1)*a*(2*m - 5*a*omega))*_direct_power(r,3) + 2*(-6 + a*(complex(0,-1)*m + a*(-6 - lambda + complex(0,6)*omega + (m - 3*a*omega)*(m - a*omega))))*_direct_power(r,4) + 2*(5 + lambda - complex(0,4)*_direct_power(a,2)*omega)*_direct_power(r,5) - (2 + lambda - 2*omega*(complex(0,3) - 2*a*m + 3*_direct_power(a,2)*omega))*_direct_power(r,6) - complex(0,2)*omega*_direct_power(r,7) + 2*_direct_power(omega,2)*_direct_power(r,8))/(2.0*sqrt(1 - _direct_power(a,2))*_direct_power(r,2)*(_direct_power(a,2) + (-2 + r)*r)*_direct_power(_direct_power(a,2) + _direct_power(r,2),1.5))
    return A0, A1, B0, B1
end

function P_to_GSN_coefficients_from_matrix(teukolsky_from_gsn_matrix, s, m, a, omega, lambda; coefficient_precision = nothing)
    if !_is_complex_frequency(omega)
        kappa = sqrt(1 - a^2)
        epsilon = 2 * omega
        tau = (epsilon - m * a) / kappa
        Tfun, _, _ = Tx(s, epsilon, kappa, tau)
        x_of_r = r -> r_to_x(r, kappa)

        return r -> begin
            Tval = Tfun(x_of_r(r))
            A0, A1, B0, B1 = _P_to_GSN_direct_coefficients_kerr(s, m, a, omega, lambda, r)
            return Tval * A0, Tval * A1, Tval * B0, Tval * B1
        end
    end

    if coefficient_precision !== nothing
        return r -> setprecision(coefficient_precision) do
            rb = big(r)
            ab = big(a)
            omegab = big(omega)
            lambdab = big(lambda)
            kappa = sqrt(1 - ab^2)
            epsilon = 2 * omegab
            tau = (epsilon - m * ab) / kappa
            Tfun, _, TpTfun = Tx(s, epsilon, kappa, tau)
            x = r_to_x(rb, kappa)
            Tval = Tfun(x)
            Xi = TpTfun(x)
            M = teukolsky_from_gsn_matrix(rb)
            M11, M12, M21, M22 = M[1, 1], M[1, 2], M[2, 1], M[2, 2]
            detM = M11 * M22 - M12 * M21
            half_kappa_inv = inv(2 * kappa)

            A0 = Tval * (M22 + M12 * Xi * half_kappa_inv) / detM
            A1 = Tval * M12 * half_kappa_inv / detM
            B0 = Tval * (-M21 - M11 * Xi * half_kappa_inv) / detM
            B1 = -Tval * M11 * half_kappa_inv / detM
            return ComplexF64(A0), ComplexF64(A1), ComplexF64(B0), ComplexF64(B1)
        end
    end

    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    Tfun, _, TpTfun = Tx(s, epsilon, kappa, tau)
    x_of_r = r -> r_to_x(r, kappa)

    return r -> begin
        x = x_of_r(r)
        Tval = Tfun(x)
        Xi = TpTfun(x)
        M = teukolsky_from_gsn_matrix(r)
        M11, M12, M21, M22 = M[1, 1], M[1, 2], M[2, 1], M[2, 2]
        detM = M11 * M22 - M12 * M21
        half_kappa_inv = inv(2 * kappa)

        A0 = Tval * (M22 + M12 * Xi * half_kappa_inv) / detM
        A1 = Tval * M12 * half_kappa_inv / detM
        B0 = Tval * (-M21 - M11 * Xi * half_kappa_inv) / detM
        B1 = -Tval * M11 * half_kappa_inv / detM
        return A0, A1, B0, B1
    end
end

function P_to_GSN_solution_from_matrix(r_from_rstar, teukolsky_from_gsn_matrix, Psoln, s, m, a, omega, lambda; coefficient_precision = nothing)
    coeffs = P_to_GSN_coefficients_from_matrix(teukolsky_from_gsn_matrix, s, m, a, omega, lambda; coefficient_precision = coefficient_precision)
    kappa = sqrt(1 - a^2)
    x_of_r = r -> r_to_x(r, kappa)

    return rs -> begin
        r = r_from_rstar(rs)
        P, Px, _, error = Psoln(x_of_r(r))
        A0, A1, B0, B1 = coeffs(r)
        X = A0 * P + A1 * Px
        Xp = B0 * P + B1 * Px
        return (X, Xp, error)
    end
end

function _GSN_to_Y_coefficients_from_maps(p_to_y_coeffs, p_to_gsn_coeffs, gsn_unit_scale)
    f1, f2, f3, _, _ = p_to_y_coeffs

    return r -> begin
        A0, A1, B0, B1 = p_to_gsn_coeffs(r)
        D = A0 * B1 - A1 * B0
        SD = gsn_unit_scale * D
        F1 = f1(r)
        F2 = f2(r)
        F3 = f3(r)
        return (
            F1 * B1 / SD,
            -F1 * A1 / SD,
            (F2 * B1 - F3 * B0) / SD,
            (-F2 * A1 + F3 * A0) / SD,
        )
    end
end

function GSN_to_Y_coeffs_minus2_in(teukolsky_from_gsn_matrix, m, a, omega, lambda, Btrans; coefficient_precision = nothing)
    p_to_y_coeffs = P_to_Y_coeffs_minus2_in(m, a, omega, lambda, Btrans)
    p_to_gsn_coeffs = P_to_GSN_coefficients_from_matrix(teukolsky_from_gsn_matrix, -2, m, a, omega, lambda; coefficient_precision = coefficient_precision)
    return _GSN_to_Y_coefficients_from_maps(p_to_y_coeffs, p_to_gsn_coeffs, Btrans)
end

function GSN_to_Y_coeffs_plus2_up(teukolsky_from_gsn_matrix, m, a, omega, lambda, Ctrans; coefficient_precision = nothing)
    p_to_y_coeffs = P_to_Y_coeffs_plus2_up(m, a, omega, lambda, Ctrans)
    p_to_gsn_coeffs = P_to_GSN_coefficients_from_matrix(teukolsky_from_gsn_matrix, 2, m, a, omega, lambda; coefficient_precision = coefficient_precision)
    return _GSN_to_Y_coefficients_from_maps(p_to_y_coeffs, p_to_gsn_coeffs, Ctrans)
end

function GSN_to_Y_solution_from_matrix(rstar_from_r, teukolsky_from_gsn_matrix, gsn_solution, s, m, a, omega, lambda, branch_transmission_amplitude; coefficient_precision = nothing)
    if s == -2
        coeffs = GSN_to_Y_coeffs_minus2_in(teukolsky_from_gsn_matrix, m, a, omega, lambda, branch_transmission_amplitude; coefficient_precision = coefficient_precision)
    elseif s == 2
        coeffs = GSN_to_Y_coeffs_plus2_up(teukolsky_from_gsn_matrix, m, a, omega, lambda, branch_transmission_amplitude; coefficient_precision = coefficient_precision)
    else
        error("GSN_to_Y_solution_from_matrix supports only s = -2 and s = +2.")
    end

    return r -> begin
        Xvals = gsn_solution(rstar_from_r(r))
        X = Xvals[1]
        Xp = Xvals[2]
        C11, C12, C21, C22 = coeffs(r)
        Y = C11 * X + C12 * Xp
        Yp = C21 * X + C22 * Xp
        return (Y, Yp, X, zero(abs(Y)))
    end
end

function Vr(r, s, a, omega, m, lambda)
    K = (r^2 + a^2) * omega - m * a
    V_r = lambda - 4im * s * omega * r - (K^2 - 2im * s * (r - 1.0) * K) / Delta(r, a)
    return V_r
end

function Vrp(r, s, a, omega, m)
    return begin
        (1/((a^2 + (-2+r)*r)^2))*(2*a*m*(a*m*(-1+r) - im*(-2+a^2-(-2+r)*r)*s) + 4*a*m*(a-r)*
        (a+r)*omega - 2*im*(a^4+2*a^2*(1+(-4+r)*r)+r^2*(6+(-4+r)*r))*s*omega - 2*(a-im*r)*
        (a+im*r)*((-3+r)*r^2+a^2*(1+r))*omega^2)
    end
end

function Tx(s, epsilon, kappa, tau)
    T = x -> exp(1im * epsilon * kappa * x) * (-x)^(-s - 1im * (epsilon + tau) / 2.0) * (1.0 - x)^(1im * (epsilon - tau) / 2.0)
    Tp = x -> T(x) * (-2 * s + 1im * (epsilon - 2 * x * epsilon * kappa + 2 * x^2 * epsilon * kappa + tau - 2 * x * tau) / (x - 1)) / (2 * x)
    TpT = x -> (-2 * s + 1im * (epsilon - 2 * x * epsilon * kappa + 2 * x^2 * epsilon * kappa + tau - 2 * x * tau) / (x - 1)) / (2 * x)
    return T, Tp, TpT
end

function d2R(R, dR, s, a, omega, m, lambda)
    return r -> Vr(r, s, a, omega, m, lambda) * R(r) / Delta(r, a) - 2 * (r - 1) * (s + 1) * dR(r) / Delta(r, a)
end

function d2R(r, R, dR, s, a, omega, m, lambda)
    return Vr(r, s, a, omega, m, lambda) * R / Delta(r, a) - 2 * (r - 1) * (s + 1) * dR / Delta(r, a)
end

function coe1_funcs(s, a, omega, m, lambda)
    V = r -> Vr(r, s, a, omega, m, lambda)
    Vp = r -> Vrp(r, s, a, omega, m)
    coe0 = r -> begin
        -4*a^7*m*omega^3+a^8*omega^4+2*a^6*omega^2*(-13+3*m^2+4*s+2*omega*(-1im+r*(-6*1im+r*omega))-3*V(r))+
        4*a^5*m*omega*(10-m^2+omega*(3*1im+r*(11*1im-3*r*omega))+3*V(r))+2*a^2*(-36*1im*r^5*omega^3+2*r^6*
        omega^4+m^2*(8+6*s+r*(-2*1im*r*(-11+4*r)*omega+3*r^3*omega^2+(-2+r)*(1+2*s-3*V(r))))+r^4*omega^2*
        (-107+16*s+50*1im*omega-9*V(r))+2*(-3+s)*s*(-2*1im*omega+V(r))+4*(2*1im*omega+V(r))+4*1im*r^3*omega*
        (24-5*s+s^2+1im*omega*(-53+14*s-3*V(r))+8*V(r)-2*s*V(r))+r^2*(24-4*1im*(51+s*(-22+5*s))*omega+36*
        (-2+s)*omega^2+2*(11-5*s+s^2+8*1im*(-4+s)*omega)*V(r)+V(r)^2)+4*1im*r*omega*(16+5*s^2+5*V(r)-2*s*
        (11+V(r)))-2*r*(24+V(r)*(22+2*(-5+s)*s+V(r))))+a^4*(24+m^4-4*1im*r^2*(-11+18*r)*omega^3+6*r^4*omega^4+
        18*V(r)-4*s*V(r)+V(r)^2-2*m^2*(7+4*s+2*omega*(3*1im+r*(4*1im-3*r*omega))+3*V(r))+4*1im*omega*(r*s*
        (-5+s-2*V(r))-5*(3+V(r))+8*r*(3+V(r))+s*(7-s+2*V(r)))+2*omega^2*(8+6*s+r*(54-32*s+r*(-73+14*s-9*V(r))+
        6*V(r))))+4*a^3*m*(22*1im*(-1+r)*r^2*omega^2-3*r^4*omega^3-1im*(-1+r)*(m^2+5*(2+V(r))-s*(-2+s+2*V(r)))-
        omega*(8+6*s+r*(26-18*s+r*(-28+m^2+4*s-6*V(r))+6*V(r))))-4*a*m*(1im*(25-11*r)*r^4*omega^2+r^6*omega^3+
        1im*(-1+r)*(-4+2*(-3+s)*s-2*r*(6+s^2+5*V(r)-2*s*(2+V(r)))+r^2*(6+s^2+5*V(r)-2*s*(2+V(r))))+r^2*omega*
        (18*(-2+s)+r*(58-14*s+4*r*s+6*V(r)-3*r*(6+V(r)))))+r*(-4*1im*r^5*(-13+6*r)*omega^3+r^7*omega^4+2*r^3*
        omega^2*(30*(-4+s)+r*(158-24*s+r*(-47+6*s-3*V(r))+6*V(r)))+(-2+r)*(4*(-2+s)*(-1+s)*V(r)-2*r*(24+V(r)*
        (26+4*(-4+s)*s+V(r)))+r^2*(24+V(r)*(26+4*(-4+s)*s+V(r))))+4*1im*r*omega*(-6*(-4+s)*(-1+s)+r^3*(s*(-5+
        s-2*V(r))+8*(3+V(r)))+2*r*(52+5*s^2+11*V(r)-2*s*(12+V(r)))+r^2*(-95-27*V(r)+s*(25-5*s+6*V(r)))))
    end
    coe1 = r -> begin
        4*(a^2+(-2+r)*r)*(3*1im*a^5*m*omega^2+r^4*(18-10*r-3*s+3*r*s)*omega^2-1im*a^6*omega^3-1im*r^6*
        omega^3+1im*a^3*m*(-8+m^2+s+2*1im*(2+3*r*(-2+s)-3*s)*omega+6*r^2*omega^2-V(r))+a^4*omega*
        (-3*1im*m^2+(2-10*r-3*s+3*r*s)*omega-3*1im*r^2*omega^2+1im*(10-3*s+V(r)))-(-1+r)*(-2+s)*
        (-2+2*s^2-2*r*(3+s*(-3+2*s)+V(r))+r^2*(3+s*(-3+2*s)+V(r)))+1im*r^2*omega*(26+4*(-6+s)*s-2*r*
        (25+s*(-19+4*s)+V(r))+r^2*(18+s*(-15+4*s)+V(r)))+a*m*(2*r^2*(-10-3*r*(-2+s)+3*s)*omega+3*1im*
        r^4*omega^2-1im*(-2+4*(-2+s)*s-2*r*(6+s*(-9+4*s)+V(r))+r^2*(6+s*(-9+4*s)+V(r))))+a^2*(2*(-1+r)*
        r^2*(-10+3*s)*omega^2-3*1im*r^4*omega^3+(-1+r)*(m^2*(-2+3*s)+(-2+s)*(-5+3*s-V(r)))-1im*omega*
        (2-4*(-2+s)*s+2*r*(15+s*(-15+4*s)+V(r))+r^2*(3*m^2-2*(14+s*(-9+2*s)+V(r))))))
    end
    coe2 = r -> begin
        96*(-2+r)*(-1+r)*r-4*1im*a^3*m*(10+m^2-(-2+s)*s)+4*a^2*(-1+r)*(24+m^2*(1+2*s))-4*1im*a*m*
        (8+2*s*(-7+2*s)-6*r*(6+(-4+s)*s)+3*r^2*(6+(-4+s)*s))-4*1im*(a^4*(-24+4*m^2-(-5+s)*s)+2*a^2*
        (-16+r*(102-72*r+m^2*(-11+6*r))+22*s+r*(-44+15*r)*s+(-5+(10-3*r)*r)*s^2)-2*1im*a^3*m*(13-
        9*s+r*(-28+m^2+4*s))-2*1im*a*m*r*(18*(-2+s)+r*(87-21*s+4*r*(-9+2*s)))+r*(12*(-4+s)*(-1+s)+
        20*r^2*(19+(-5+s)*s)-5*r^3*(24+(-5+s)*s)-6*r*(52+s*(-24+5*s))))*omega+4*(11*1im*a^5*m+22*
        1im*a^3*m*r*(-2+3*r)+5*1im*a*m*r^3*(-20+11*r)+a^4*(27-16*s+r*(-73+6*m^2+14*s))+r^3*(60*
        (-4+s)+r*(395-60*s+3*r*(-47+6*s)))+2*a^2*r*(18*(-2+s)+r*(159-42*s+r*(-107+3*m^2+16*s))))*
        omega^2-8*1im*(a^2+r^2)*(3*a^4-3*1im*a^3*m*r-3*1im*a*m*r^3+3*r^3*(-13+7*r)+a^2*r*(-11+24*r))*
        omega^3+8*r*(a^2+r^2)^3*omega^4+4*(a^2+(-2+r)*r)*(-3-1im*a*m+3*r+s-r*s+1im*(a^2+r^2)*omega)*
        V(r)^2-2*(a^2+(-2+r)*r)*(-4+26*r+6*s+4*r*(-4+s)*s-2*s^2-2*1im*a*m*(-1+r)*(-5+2*s)+r^2*(-13-
        2*s*(-4+s+2*1im*omega)+22*1im*omega)-6*a^3*m*omega-6*a*m*r^2*omega+4*1im*r^3*(-4+s)*omega+
        3*a^4*omega^2+3*r^4*omega^2+a^2*(-9+3*m^2+2*s+2*1im*(5+2*r*(-4+s)-2*s)*omega+6*r^2*omega^2))*
        Vp(r)+2*V(r)*(2*((-1+r)*(32*(-2+r)*r+(-4-25*(-2+r)*r)*s+(6+11*(-2+r)*r)*s^2-2*(-1+r)^2*s^3)+
        1im*r^2*(92+2*r*(-79+29*r)+(-36+(62-25*r)*r)*s+4*(-1+r)^2*s^2)*omega+3*1im*a^5*m*omega^2+
        r^4*(33-19*r+3*(-1+r)*s)*omega^2-1im*a^6*omega^3-1im*r^6*omega^3+a^3*m*(1im*(-13+m^2+3*s)+
        2*(-5-3*r*(-4+s)+3*s)*omega+6*1im*r^2*omega^2)+a*m*(-1im*(8+21*(-2+r)*r-12*s-15*(-2+r)*r*s+
        4*(-1+r)^2*s^2)+2*r^2*(-19-3*r*(-4+s)+3*s)*omega+3*1im*r^4*omega^2)+a^2*((-1+r)*(32+m^2*
        (-5+3*s)+s*(-21+5*s))-1im*(-8+r*(94+(-76+3*m^2)*r)+12*s+2*r*(-23+15*r)*s-4*(-1+r)^2*s^2)*
        omega+2*(-1+r)*r^2*(-19+3*s)*omega^2-3*1im*r^4*omega^3)+a^4*omega*(18*1im-3*1im*m^2+s*
        (-5*1im+3*(-1+r)*omega)+omega*(5+r*(-19-3*1im*r*omega))))+(a^2+(-2+r)*r)^2*Vp(r))
    end
    coe3 = r -> begin
        4*(8*(-2+r)*r*(1+3*(-2+r)*r)+(8-(-2+r)*r*(20+39*(-2+r)*r))*s+(-4+(-2+r)*r*(32+39*(-2+r)*
        r))*s^2-4*(-1+r)^2*(2+5*(-2+r)*r)*s^3+4*(-1+r)^4*s^4)-4*1im*r^2*(128+r*(-456+(389-96*r)*r)-
        178*s+r*(476+r*(-395+101*r))*s-(-1+r)*(70+r*(-120+47*r))*s^2+8*(-1+r)^3*s^3)*omega-2*r^4*
        (408+r*(-502+147*r)-150*s+4*(55-19*r)*r*s+12*(-1+r)^2*s^2)*omega^2-4*a^7*m*omega^3+4*1im*
        r^6*(25+2*r*(-6+s)-2*s)*omega^3+a^8*omega^4+r^8*omega^4+2*a^6*omega^2*(-33+3*m^2+10*s+2*
        1im*(-1+2*r*(-6+s)-2*s)*omega+2*r^2*omega^2)-4*a^5*m*omega*(-22+m^2+6*s+1im*(-3-23*r+6*
        (-1+r)*s)*omega+3*r^2*omega^2)+4*a^3*m*(-1im*(-1+r)*(22+m^2*(1+2*s)+s*(-32+9*s))-(8+r*(90+
        (-76+m^2)*r)+14*s+2*r*(-37+26*r)*s-12*(-1+r)^2*s^2)*omega-2*1im*(-1+r)*r^2*(-23+6*s)*
        omega^2-3*r^4*omega^3)+4*a*m*(1im*(-1+r)*(4-18*(-2+r)*r+2*(1+17*(-2+r)*r)*s-9*(2+3*(-2+
        r)*r)*s^2+8*(-1+r)^2*s^3)+2*r^2*(58+r*(-85+27*r)-41*s+(63-23*r)*r*s+6*(-1+r)^2*s^2)*omega-
        1im*r^4*(49-23*r+6*(-1+r)*s)*omega^2-r^6*omega^3)+a^4*(64-22*m^2+m^4-44*s+4*m^2*s+12*s^2+
        4*1im*(-45+80*r+57*s-61*r*s+15*(-1+r)*s^2+m^2*(-3-10*r+6*(-1+r)*s))*omega+2*(8+3*r*(58+
        (-71+2*m^2)*r)+14*s+4*r*(-29+24*r)*s-12*(-1+r)^2*s^2)*omega^2+4*1im*r^2*(23+6*r*(-6+s)-
        6*s)*omega^3+6*r^4*omega^4)+2*a^2*(4*(4+20*(-2+r)*r+(-18-29*(-2+r)*r)*s+(20+23*(-2+r)*r)*
        s^2-6*(-1+r)^2*s^3)-4*1im*(-2+r*(-72+(197-88*r)*r)-s+r*(104+r*(-188+81*r))*s-(-1+r)*(9+
        r*(-46+31*r))*s^2+4*(-1+r)^3*s^3)*omega+r^2*(-232+676*r-327*r^2+2*(82+3*r*(-56+27*r))*
        s-24*(-1+r)^2*s^2)*omega^2+2*1im*r^4*(49+6*r*(-6+s)-6*s)*omega^3+2*r^6*omega^4+m^2*(8+
        2*(7-6*s)*s+r*(6+8*s*(-4+3*s))+r^2*(-3-4*s*(-4+3*s+3*1im*omega)+46*1im*omega)+4*1im*
        r^3*(-5+3*s)*omega+3*r^4*omega^2))-(a^2+(-2+r)*r)*(2*(-6*(-2+s)*(-1+s)-2*1im*a*m*(-1+r)*
        (-7+4*s)+2*r*(25-22*s+6*s^2)+r^2*(-25-2*s*(-11+3*s+4*1im*omega)+34*1im*omega)-6*a^3*m*
        omega-6*a*m*r^2*omega+8*1im*r^3*(-3+s)*omega+3*a^4*omega^2+3*r^4*omega^2+a^2*(-13+3*m^2+
        4*s+2*1im*(7+4*r*(-3+s)-4*s)*omega+6*r^2*omega^2))*V(r)-(a^2+(-2+r)*r)*V(r)^2+4*(a^2+
        (-2+r)*r)*((-1+r)*(-2+s)+1im*(a*m-(a^2+r^2)*omega))*Vp(r))
    end
    return (coe0, coe1, coe2, coe3)
end

function coe2_funcs(s, a, omega, m, lambda)
    V = r -> Vr(r, s, a, omega, m, lambda)
    Vp = r -> Vrp(r, s, a, omega, m)
    coe0 = r -> begin
        (1/((a^2+(-2+r)*r)^4))*(-4*a^7*m*omega^3+a^8*omega^4+4*a^3*m*(1im*(-1+r)*(6+5*m^2+s*(6+s))+
        (32+10*s-r*(42+22*s+r*(m^2-4*(1+s))))*omega+2*1im*(-1+r)*r^2*omega^2-3*r^4*omega^3)+4*a*m*
        (-1im*(-1+r)*(2*(2+s)*(3+s)-2*r*(6+s*(4+s))+r^2*(6+s*(4+s)))+2*r^2*(-2+7*s+r*(3-5*s+r*(3+
        2*s)))*omega+1im*r^4*(13+r)*omega^2-r^6*omega^3)+r^2*omega*(4*1im*(-6*s*(3+s)+r^3*s*(3+s)+
        2*r*s*(16+5*s)+r^2*(3-5*s*(3+s)))-2*r^2*((26-5*r)*r+2*(9+(-4+r)*r)*s)*omega+4*1im*r^4*(-9+
        2*r)*omega^2+r^6*omega^3)+a^4*(m^2*(10+m^2+8*s)-4*1im*(-3+m^2*(-15+8*r)+s*(1+s-r*(3+s)))*
        omega+2*(-2*(16+5*s)+r*(30+r*(3+6*m^2-10*s)+32*s))*omega^2+4*1im*r^2*(1+6*r)*omega^3+6*r^4*
        omega^4)+2*a^2*(-m^2*(32+10*s+3*(-2+r)*r*(9+2*s))-2*1im*(r*(-24+r*(18+m^2*(-1+8*r)))-6*r*
        (6+(-6+r)*r)*s-2*(-1+r)*(1+(-4+r)*r)*s^2+2*(6+5*s))*omega+r^2*(8-28*s+r*(4+r*(9+3*m^2-8*s)+
        40*s))*omega^2+2*1im*r^4*(-13+6*r)*omega^3+2*r^6*omega^4)+2*a^6*omega^2*(-1+3*m^2-4*s+2*
        omega*(5*1im+r*(2*1im+r*omega)))-4*a^5*m*omega*(2+m^2+omega*(15*1im+r*(-1im+3*r*omega)))-
        (a^2+(-2+r)*r)*V(r)*(2*(-2*(2+s)*(3+s)+2*1im*a*m*(-1+r)*(7+2*s)+2*r*(9+2*s*(4+s))+r^2*(-
        9-2*s*(4+s-2*1im*omega)+2*1im*omega)-6*a^3*m*omega-6*a*m*r^2*omega-4*1im*r^3*(2+s)*omega+
        3*a^4*omega^2+3*r^4*omega^2+a^2*(3+3*m^2+2*s-2*1im*(-7+4*r+2*(-1+r)*s)*omega+6*r^2*
        omega^2))-(a^2+(-2+r)*r)*V(r)))
    end
    coe1 = r -> begin
        (1/((a^2+(-2+r)*r)^3))*4*(-((-1+r)*(1+s)*(2+s)*(2*(3+s)+(-2+r)*r*(3+2*s)))-1im*r^2*
        (2+4*s*(3+s)+r^2*(2+s)*(3+4*s)-2*r*(1+s*(9+4*s)))*omega-3*1im*a^5*m*omega^2+r^4*(6-
        3*s+r*(2+3*s))*omega^2+1im*a^6*omega^3+1im*r^6*omega^3+a^4*omega*(1im*(2+3*m^2+s)+2*
        (-5+r)*omega+3*(-1+r)*s*omega+3*1im*r^2*omega^2)-1im*a^3*m*(4+m^2+3*s-2*1im*(-10-3*
        s+3*r*(2+s))*omega+6*r^2*omega^2)+a*m*(1im*(22+4*s*(5+s)-2*r*(2+s)*(9+4*s)+r^2*(2+s)*
        (9+4*s))+2*r^2*(2+3*s-3*r*(2+s))*omega-3*1im*r^4*omega^2)+a^2*((-1+r)*(3*(1+s)*(2+s)+
        m^2*(10+3*s))+1im*(-22-4*s*(5+s)+2*r*(3+s)*(5+4*s)+r^2*(3*m^2-2*(2+s)*(1+2*s)))*omega+
        2*(-1+r)*r^2*(2+3*s)*omega^2+3*1im*r^4*omega^3)+(a^2+(-2+r)*r)*(-((-1+r)*(2+s))+1im*
        (a*m-(a^2+r^2)*omega))*V(r))
    end
    coe2 = r -> begin
        (1/((a^2+(-2+r)*r)^5))*(4*(2*a^8*omega^3*(1im+omega)+a^7*m*omega^2*(1im+2*(-4+r)*omega)+
        a^4*(-m^2*(-1+r)*(47+2*m^2+22*s)+1im*(2*m^2*(60+r*(-83+16*r))-r^2*s*(3+s)+2*(2+s)*(12+s)-
        2*r*(6+s)*(5+3*s))*omega+(-8*(16+5*s)+r*(226+108*s+r*(-93-12*m^2*(-1+r)-72*s+r*(9+14*
        s))))*omega^2+2*1im*(2-23*r)*r^2*omega^3)+a^5*m*(1im*(6+5*m^2+s*(6+s))+2*(-29+m^2*(-4+
        3*r)-11*s+4*r*(3+s))*omega-1im*(120+(-122+r)*r)*omega^2+6*(-2+r)*r^2*omega^3)+a^3*m*
        (-1im*(72+60*(-2+r)*r+5*m^2*(8+7*(-2+r)*r)+6*(11+9*(-2+r)*r)*s+2*(6+5*(-2+r)*r)*s^2)+
        2*(8*(16+5*s)+r*(-258-92*s+r*(164+m^2*(-2+3*r)+(70-4*r)*s)))*omega-1im*r^2*(8+r*(-68+
        5*r))*omega^2+6*r^5*omega^3)+r^2*omega*(-1im*(24*s*(3+s)+3*r^4*s*(3+s)+10*r^2*s*(16+
        5*s)-4*r*s*(43+14*s)-6*r^3*(-2+3*s*(3+s)))+r^3*(26+28*s+r*(29-5*r+2*(-4+r)*s))*omega-
        2*1im*(-6+r)*r^4*(3+r)*omega^2-2*r^6*omega^3)+a*m*(1im*(96+30*(-2+r)*r*(4+(-2+r)*r)+
        2*(40+(-2+r)*r*(47+10*(-2+r)*r))*s+(16+5*(-2+r)*r*(4+(-2+r)*r))*s^2)-2*r^2*(8+3*r*
        (2+r)*(-3+4*r)-28*s+r*(52+r*(-25+8*r))*s)*omega-3*1im*r^5*(18+r)*omega^2+2*r^6*(2+r)*
        omega^3)+a^2*(m^2*(-1+r)*(9*(-2+r)*r*(9+2*s)+8*(16+5*s))+1im*(-16*(2+s)*(3+s)+2*r^3*
        (60-11*m^2+84*s+22*s^2)+5*r^4*(8*m^2-s*(3+s))+4*r^2*(-60+m^2-5*s*(15+4*s))+4*r*(60+s*
        (65+16*s)))*omega-r^2*(8*(-2+7*s)+r*(20-88*s+r*(75+r*(3+6*m^2-10*s)+80*s)))*omega^2-
        2*1im*r^5*(7+2*r)*omega^3-4*r^6*omega^4)+a^6*omega*(1im*s*(3+s)+11*omega-2*m^2*(4*
        1im+3*(-2+r)*omega)+omega*(40*1im*omega+r*(7+6*s-26*1im*omega+4*r*omega*(1im+omega)))))+
        2*(a^2+(-2+r)*r)*(2*(a^2+(-2+r)*r)*(-((-1+r)*(3+s))+1im*(a*m-(a^2+r^2)*omega))*V(r)^2-
        (a^2+(-2+r)*r)*(-2*(2+s)*(3+s)+2*1im*a*m*(-1+r)*(7+2*s)+2*r*(9+2*s*(4+s))+r^2*(-9-2*
        s*(4+s-2*1im*omega)+2*1im*omega)-6*a^3*m*omega-6*a*m*r^2*omega-4*1im*r^3*(2+s)*omega+
        3*a^4*omega^2+3*r^4*omega^2+a^2*(3+3*m^2+2*s-2*1im*(-7+4*r+2*(-1+r)*s)*omega+6*r^2*
        omega^2))*Vp(r)+V(r)*(2*(-((-1+r)*(3+s)*(16+8*(-2+r)*r+12*s+7*(-2+r)*r*s+2*(-1+r)^2*
        s^2))-1im*r^2*(4+4*s*(4+s)-2*r*(3+s)*(1+4*s)+r^2*(2+s)*(9+4*s))*omega-3*1im*a^5*m*
        omega^2+r^4*(9+5*r+3*(-1+r)*s)*omega^2+1im*a^6*omega^3+1im*r^6*omega^3+a^4*omega*(3*
        1im*(2+m^2+s)+(-19+5*r+3*(-1+r)*s)*omega+3*1im*r^2*omega^2)-1im*a^3*m*(11+m^2+5*s-
        2*1im*(-19-3*s+3*r*(4+s))*omega+6*r^2*omega^2)+a*m*(1im*(53*(-2+r)*r+27*(-2+r)*r*
        s+4*(-1+r)^2*s^2+32*(2+s))+2*r^2*(5+3*s-3*r*(4+s))*omega-3*1im*r^4*omega^2)+a^2*
        ((-1+r)*(m^2*(19+3*s)+(3+s)*(8+5*s))+1im*(-4*(4+s)^2+r*(86+50*s+8*s^2)+r^2*(3*m^2-
        2*(2+s)*(3+2*s)))*omega+2*(-1+r)*r^2*(5+3*s)*omega^2+3*1im*r^4*omega^3))+(a^2+(-2+
        r)*r)^2*Vp(r))))
    end
    coe3 = r -> begin
        (1/((a^2+(-2+r)*r)^4))*(4*(1+s)*(2+s)*(48+15*(-2+r)*r*(4+(-2+r)*r)+4*(-1+r)^2*(7+4*
        (-2+r)*r)*s+4*(-1+r)^4*s^2)+4*1im*r^2*(-8+r*(16+r*(-7+24*r))-70*s+r*(148+r*(-121+
        59*r))*s+(-1+r)*(46+r*(-64+39*r))*s^2+8*(-1+r)^3*s^3)*omega-2*r^4*(r*(82+7*r)-6*s+
        4*r*(-1+7*r)*s+12*(-1+r)^2*s^2)*omega^2-4*a^7*m*omega^3-4*1im*r^6*(13+2*(-1+r)*s)*
        omega^3+a^8*omega^4+r^8*omega^4-4*a^5*m*omega*(14+m^2+6*s-1im*(13*(-3+r)+6*(-1+r)*
        s)*omega+3*r^2*omega^2)+4*a^3*m*(1im*(-1+r)*(74+m^2*(13+2*s)+3*s*(24+5*s))+(192-r*
        (266+(-52+m^2)*r)+98*s+2*r*(-79+26*r)*s+12*(-1+r)^2*s^2)*omega+2*1im*(-1+r)*r^2*(13+
        6*s)*omega^2-3*r^4*omega^3)+4*a*m*(-1im*(-1+r)*(188+114*(-2+r)*r+214*s+142*(-2+r)*
        r*s+(74+59*(-2+r)*r)*s^2+8*(-1+r)^2*s^3)+2*r^2*(6+r^2*(3+s)*(11+6*s)+s*(23+6*s)-3*
        r*(7+s*(15+4*s)))*omega+1im*r^4*(13*(1+r)+6*(-1+r)*s)*omega^2-r^6*omega^3)+2*a^6*
        omega^2*(3+3*m^2+2*s-4*1im*(-1+r)*s*omega+2*omega*(13*1im+r^2*omega))+2*a^2*(-12*
        (1+s)*(2+s)*(2*(3+s)+(-2+r)*r*(5+2*s))+4*1im*(-94+(188-119*r)*r-107*s+r*(240+r*
        (-164+15*r))*s+(-1+r)*(37+r*(-58+15*r))*s^2+4*(-1+r)^3*s^3)*omega-r^2*(24+r*(-116+
        11*r)+92*s+18*r*(-8+3*r)*s+24*(-1+r)^2*s^2)*omega^2-2*1im*r^4*(13+6*(-1+r)*s)*omega^3+
        2*r^6*omega^4+m^2*(2*r*(167+4*s*(22+3*s))-2*(96+s*(49+6*s))+r^2*(-167-4*s*(22+3*s-
        3*1im*omega)+26*1im*omega)-4*1im*r^3*(13+3*s)*omega+3*r^4*omega^2))+a^4*(m^4+12*
        (1+s)*(2+s)-4*1im*(-49-9*s*(5+s)+r*(24+s*(29+9*s)))*omega-2*(192+(-198+r)*r+98*s+
        4*r*(-35+6*r)*s+12*(-1+r)^2*s^2)*omega^2-4*1im*r^2*(-13+6*(-1+r)*s)*omega^3+6*r^4*
        omega^4+2*m^2*(25+2*s*(5-6*1im*(-1+r)*omega)+2*omega*(39*1im+r*(-26*1im+3*r*omega))))-
        (a^2+(-2+r)*r)*(2*(-6*(2+s)*(3+s)+2*1im*a*m*(-1+r)*(13+4*s)+2*r*(29+26*s+6*s^2)-6*
        a^3*m*omega-6*a*m*r^2*omega-8*1im*r^3*(2+s)*omega+3*a^4*omega^2+3*r^4*omega^2+r^2*
        (-29-26*s-6*s^2+2*1im*(3+4*s)*omega)+a^2*(7+3*m^2+4*s-2*1im*(-13+8*r+4*(-1+r)*s)*
        omega+6*r^2*omega^2))*V(r)-(a^2+(-2+r)*r)*V(r)^2+4*(a^2+(-2+r)*r)*((-1+r)*(2+s)+1im*
        (-a*m+(a^2+r^2)*omega))*Vp(r)))
    end
    return (coe0, coe1, coe2, coe3)
end

function PlusTwo_to_MinusTwo_Inf(Pplus, m, a, omega, lambda)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    (coe0, coe1, coe2, coe3) = coe1_funcs(2, a, omega, m, lambda)
    Tplus, Tpplus, _ = Tx(2, epsilon, kappa, tau)
    Tminus, Tpminus, TpTminus = Tx(-2, epsilon, kappa, tau)
    f4p = 32 * kappa^(5 - 2im * omega) * exp(- 1im * omega * (1 + kappa))
    f4m = kappa^(- 3 - 2im * omega) * exp(- 1im * omega * (1 + kappa)) / 8
    function Rplus(r)
        x = r_to_x(r, kappa)
        P, Pp, Ppp, error = Pplus(x)
        R = Tplus(x) * P / f4p
        Rp = - (Tpplus(x) * P + Tplus(x) * Pp) / (2 * kappa) / f4p
        Rpp = d2R(r, R, Rp, 2, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    function Pminus(x)
        r = x_to_r(x, kappa)
        R, Rp, Rpp, error = Rplus(r)
        Rm = (coe0(r) * R + coe1(r) * Rp) / epsilon^4
        Rmp = (coe2(r) * R + coe3(r) * Rp) / epsilon^4
        P = (Rm / Tminus(x)) * f4m
        Pp = (- 2 * kappa * Rmp / Tminus(x) - TpTminus(x) * Rm / Tminus(x)) * f4m
        Ppp = d2P(x, P, Pp, -2, epsilon, kappa, tau, lambda + 4)
        return (P, Pp, Ppp, error)
    end
    return Rplus, Pminus
end

function PlusTwo_to_MinusTwo_Hor(Pplus, m, a, omega, lambda)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    rp = 1 + kappa
    p = kappa * (epsilon + tau) / (2 * rp)
    err_scale = kappa / (4 * rp)
    (coe0, coe1, coe2, coe3) = coe1_funcs(2, a, omega, m, lambda)
    Tplus, Tpplus, _ = Tx(2, epsilon, kappa, tau)
    Tminus, Tpminus, TpTminus = Tx(-2, epsilon, kappa, tau)
    f1p = 16 * kappa^(2im * p + 4) * exp(1im * p * rp)
    f1m = kappa^(2im * p - 4) * exp(1im * p * rp) / 16
    Q = lambda + 6
    Csq = (Q^2 + 4*a*m*omega-4*a^2*omega^2)*((Q-2)^2+
            36*a*m*omega-36*a^2*omega^2)+(2*Q-1)*(96*
            a^2*omega^2-48*a*omega*m)+144*omega^2*(1-a^2)
    horizon_factor = 256 * rp^4 * (p + 4im * err_scale) * p * (p^2 + 4 * err_scale^2)
    function Rplus(r)
        x = r_to_x(r, kappa)
        P, Pp, Ppp, error = Pplus(x)
        R = Tplus(x) * P / f1p
        Rp = - (Tpplus(x) * P + Tplus(x) * Pp) / (2 * kappa) / f1p
        Rpp = d2R(r, R, Rp, 2, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    function Pminus(x)
        r = x_to_r(x, kappa)
        R, Rp, Rpp, error = Rplus(r)
        Rm = (coe0(r) * R + coe1(r) * Rp) * horizon_factor / Csq
        Rmp = (coe2(r) * R + coe3(r) * Rp) * horizon_factor / Csq
        P = (Rm / Tminus(x)) * f1m
        Pp = (- 2 * kappa * Rmp / Tminus(x) - TpTminus(x) * Rm / Tminus(x)) * f1m
        Ppp = d2P(x, P, Pp, -2, epsilon, kappa, tau, lambda + 4)
        return (P, Pp, Ppp, error)
    end
    return Rplus, Pminus
end

function MinusTwo_to_PlusTwo_Inf(Pminus, m, a, omega, lambda)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    (coe0, coe1, coe2, coe3) = coe2_funcs(-2, a, omega, m, lambda)
    Tminus, Tpminus, TpTminus = Tx(-2, epsilon, kappa, tau)
    Tplus, Tpplus, TpTplus = Tx(2, epsilon, kappa, tau)
    f3 = 2.0 * kappa^(2im * omega + 1) * exp(1im * omega * (1 + kappa))
    function Rminus(r)
        x = r_to_x(r, kappa)
        P, Pp, Ppp, error = Pminus(x)
        R = Tminus(x) * P / f3
        Rp = - (Tpminus(x) * P + Tminus(x) * Pp) / (2 * kappa) / f3
        Rpp = d2R(r, R, Rp, -2, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    function Pplus(x)
        r = x_to_r(x, kappa)
        R, Rp, Rpp, error = Rminus(r)
        Rpl = (coe0(r) * R + coe1(r) * Rp) / epsilon^4
        Rplp = (coe2(r) * R + coe3(r) * Rp) / epsilon^4
        P = (Rpl / Tplus(x)) * f3
        Pp = (- 2 * kappa * Rplp / Tplus(x) - TpTplus(x) * Rpl / Tplus(x)) * f3
        Ppp = d2P(x, P, Pp, 2, epsilon, kappa, tau, lambda - 4)
        return (P, Pp, Ppp, error)
    end
    return Rminus, Pplus
end

function MinusTwo_to_PlusTwo_Hor(Pminus, m, a, omega, lambda)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    rp = 1 + kappa
    p = kappa * (epsilon + tau) / (2 * rp)
    err_scale = kappa / (4 * rp)
    (coe0, coe1, coe2, coe3) = coe2_funcs(-2, a, omega, m, lambda)
    Tminus, Tpminus, TpTminus = Tx(-2, epsilon, kappa, tau)
    Tplus, Tpplus, TpTplus = Tx(2, epsilon, kappa, tau)
    f2 = kappa^(- 2im * p) * exp(- 1im * p * rp)
    Q = lambda + 2
    Csq = (Q^2 + 4*a*m*omega-4*a^2*omega^2)*((Q-2)^2+
            36*a*m*omega-36*a^2*omega^2)+(2*Q-1)*(96*
            a^2*omega^2-48*a*omega*m)+144*omega^2*(1-a^2)
    horizon_factor = 256 * rp^4 * (p - 4im * err_scale) * p * (p^2 + 4 * err_scale^2)
    function Rminus(r)
        x = r_to_x(r, kappa)
        P, Pp, Ppp, error = Pminus(x)
        R = Tminus(x) * P / f2
        Rp = - (Tpminus(x) * P + Tminus(x) * Pp) / (2 * kappa) / f2
        Rpp = d2R(r, R, Rp, -2, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end
    function Pplus(x)
        r = x_to_r(x, kappa)
        R, Rp, Rpp, error = Rminus(r)
        Rpl = (coe0(r) * R + coe1(r) * Rp) * horizon_factor / Csq
        Rplp = (coe2(r) * R + coe3(r) * Rp) * horizon_factor / Csq
        P = (Rpl / Tplus(x)) * f2
        Pp = (- 2 * kappa * Rplp / Tplus(x) - TpTplus(x) * Rpl / Tplus(x)) * f2
        Ppp = d2P(x, P, Pp, 2, epsilon, kappa, tau, lambda - 4)
        return (P, Pp, Ppp, error)
    end
    return Rminus, Pplus
end

function Ex(x, s)
    return (1 - x)^(-(s + 1) / 2) * (- x)^(-(s + 1) / 2)
end

@inline function _rstar_from_r(a, r)
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    return r + (2 * rp) / (rp - rm) * log((r - rp) / 2) - (2 * rm) / (rp - rm) * log((r - rm) / 2)
end

@inline function _y_phase_argument(r, a, m, omega)
    kappa = sqrt(1 - a^2)
    rp = 1 + kappa
    rm = 1 - kappa
    return omega * _rstar_from_r(a, r) + a * m * log((r - rm) / (r - rp)) / (rp - rm)
end

@inline function _y_phase_argument_prime(r, a, m, omega)
    Δ = Delta(r, a)
    return ((r^2 + a^2) * omega - a * m) / Δ
end

@inline function _y_phase_argument_second(r, a, m, omega)
    Δ = Delta(r, a)
    Δp = 2 * (r - 1)
    K = (r^2 + a^2) * omega - a * m
    Kp = 2 * r * omega
    return (Kp * Δ - K * Δp) / (Δ^2)
end

@inline function _tp_over_t_prime(s, epsilon, kappa, tau, x)
    num = epsilon - 2 * x * epsilon * kappa + 2 * x^2 * epsilon * kappa + tau - 2 * x * tau
    den = x - 1
    A = -2 * s + 1im * num / den
    Aprime = 1im * (((-2 * epsilon * kappa + 4 * x * epsilon * kappa - 2 * tau) * den) - num) / (den^2)
    return Aprime / (2 * x) - A / (2 * x^2)
end

function P_to_Y_coeffs_minus2_in(m, a, omega, lambda, Btrans)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    Tfun, _, TpTfun = Tx(-2, epsilon, kappa, tau)
    x_of_r = r -> r_to_x(r, kappa)
    phase = r -> _y_phase_argument(r, a, m, omega)
    phase_p = r -> _y_phase_argument_prime(r, a, m, omega)
    phase_pp = r -> _y_phase_argument_second(r, a, m, omega)
    xspace_data = r -> begin
        x = x_of_r(r)
        q1, q2 = calculate_q1_q2(x, -2, epsilon, kappa, tau, lambda)
        Tval = Tfun(x)
        TpTval = TpTfun(x)
        phase_r = phase(r)
        phase_pr = phase_p(r)
        phase_ppr = phase_pp(r)
        U = TpTval - 1im * (-2 * kappa * phase_pr) + 4 * kappa / r
        Up = _tp_over_t_prime(-2, epsilon, kappa, tau, x) - 1im * (4 * kappa^2 * phase_ppr) + 8 * kappa^2 / r^2
        pref = Btrans * sqrt(r^2 + a^2) * Tval / (4 * kappa^2)
        return (; x, q1, q2, Tval, TpTval, phase_r, phase_pr, phase_ppr, U, Up, pref)
    end
    u = r -> begin
        data = xspace_data(r)
        return -1im * data.phase_pr - 2 / r - data.TpTval / (2 * kappa)
    end
    G = r -> begin
        data = xspace_data(r)
        return Btrans * exp(-1im * data.phase_r) * data.Tval / r^2
    end
    func1 = r -> G(r)
    func2 = r -> G(r) * u(r)
    func3 = r -> -G(r) / (2 * kappa)
    func4 = r -> begin
        data = xspace_data(r)
        return data.pref * (data.Up + data.U^2 - data.q2 / (data.x * (1 - data.x)))
    end
    func5 = r -> begin
        data = xspace_data(r)
        return data.pref * (2 * data.U - data.q1 / (data.x * (1 - data.x)))
    end
    return func1, func2, func3, func4, func5
end

function P_to_Y_coeffs_plus2_up(m, a, omega, lambda, Ctrans)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    Tfun, _, TpTfun = Tx(2, epsilon, kappa, tau)
    x_of_r = r -> r_to_x(r, kappa)
    phase = r -> _y_phase_argument(r, a, m, omega)
    phase_p = r -> _y_phase_argument_prime(r, a, m, omega)
    phase_pp = r -> _y_phase_argument_second(r, a, m, omega)
    xspace_data = r -> begin
        x = x_of_r(r)
        q1, q2 = calculate_q1_q2(x, 2, epsilon, kappa, tau, lambda)
        phase_r = phase(r)
        phase_pr = phase_p(r)
        phase_ppr = phase_pp(r)
        Tval = Tfun(x)
        TpTval = TpTfun(x)
        u = 1im * phase_pr + 2 * (2 * r - 2) / Delta(r, a) - 2 / r - TpTval / (2 * kappa)
        up = 1im * phase_ppr + 2 * (2 * Delta(r, a) - (2 * r - 2)^2) / (Delta(r, a)^2) + 2 / r^2 + _tp_over_t_prime(2, epsilon, kappa, tau, x) / (4 * kappa^2)
        ypref = Ctrans * exp(1im * phase_r) * Delta(r, a)^2 * Tval / r^2
        xpref = Ctrans * sqrt(r^2 + a^2) * Delta(r, a)^2 * Tval / (4 * kappa^2)
        return (; x, q1, q2, u, up, ypref, xpref)
    end
    func1 = r -> begin
        data = xspace_data(r)
        return data.ypref
    end
    func2 = r -> begin
        data = xspace_data(r)
        return data.ypref * data.u
    end
    func3 = r -> begin
        data = xspace_data(r)
        return -data.ypref / (2 * kappa)
    end
    func4 = r -> begin
        data = xspace_data(r)
        return data.xpref * (4 * kappa^2 * (data.u^2 + data.up) - data.q2 / (data.x * (1 - data.x)))
    end
    func5 = r -> begin
        data = xspace_data(r)
        return data.xpref * (-4 * kappa * data.u - data.q1 / (data.x * (1 - data.x)))
    end
    return func1, func2, func3, func4, func5
end

function calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    q1_x = (1.0 - s - 1im * epsilon - 1im * tau) -
           2.0 * (1.0 - 1im * epsilon * kappa - 1im * tau) * x -
           2im * epsilon * kappa * x^2
    q2_x = lambda - epsilon^2 + 1im * epsilon * kappa * (1.0 - 2 * s) +
           s * (s + 1.0) + tau * (1im + tau) +
           2.0 * epsilon * kappa * (-1im + 1im * s + epsilon - tau) * x

    return q1_x, q2_x
end

function d2P(x, P, dP, s, epsilon, kappa, tau, lambda)
    q1, q2 = calculate_q1_q2(x, s, epsilon, kappa, tau, lambda)
    return - (q1 * dP + q2 * P) / (x) * (1 - x)
end

function calculate_q3_q4(x, s, epsilon, kappa, tau, lambda)
    q3_x = - 2 * s * x + s + 1im * tau - 1im * epsilon *
            (2 * kappa * x * (x - 1) - 2 * x + 3) - 1
    q4_x = lambda - (2 * epsilon - 1im) * (epsilon + 1im * (s + 1im * tau - 1)) / x -
            1im * epsilon * (kappa - 2 * s + 1) + 2 * s +
            2 * kappa * epsilon^2 - 2 * kappa * tau * epsilon

    return q3_x, q4_x
end

function calculate_q5_q6(x, s, epsilon, kappa, tau, lambda)
    q5_x = s * (2 * x - 3) + 1im * (epsilon + tau) + 2im * epsilon * kappa * x^2 -
            2im * (kappa + 1) * x * epsilon - 1

    q6_x = (2 * s + 1) * (s + 1 - 1im * (tau + epsilon)) / x +
            lambda + 1im * (kappa + 1) * (2 * s + 1) * epsilon

    return q5_x, q6_x
end

function _relative_residual(term1, term2, term3)
    denom = max(abs(term1), abs(term2), abs(term3))
    denom == 0 && return 0.0
    return abs(term1 + term2 + term3) / denom
end

function evaluate_power_series(x, coeffs::AbstractVector)
    n = length(coeffs)
    n == 0 && return (0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)

    val = coeffs[end]
    dval = zero(val)
    ddval = zero(val)

    @inbounds for i in (n - 1):-1:1
        ddval = muladd(ddval, x, 2 * dval)
        dval = muladd(dval, x, val)
        val = muladd(val, x, coeffs[i])
    end

    return (val, dval, ddval)
end

function evaluate_inverse_power_series(x, coeffs::AbstractVector)
    n = length(coeffs)
    n == 0 && return (0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)

    u = inv(x)
    g = coeffs[n]
    gu = zero(g)
    guu = zero(g)

    @inbounds for i in (n - 1):-1:1
        guu = muladd(guu, u, 2 * gu)
        gu = muladd(gu, u, g)
        g = muladd(g, u, coeffs[i])
    end

    up = -u^2
    upp = 2 * u^3
    fx = g
    fpx = gu * up
    fppx = guu * up^2 + gu * upp
    return (fx, fpx, fppx)
end

function evaluate_frobenius_series(rho, coeffs::AbstractVector)
    return evaluate_frobenius_series(rho, 0, coeffs)
end

function evaluate_frobenius_series(rho, sigma, coeffs::AbstractVector)
    n = length(coeffs)
    n == 0 && return (0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)

    g = coeffs[end]
    gu = zero(g)
    guu = zero(g)

    @inbounds for i in (n - 1):-1:1
        guu = muladd(guu, rho, 2 * gu)
        gu = muladd(gu, rho, g)
        g = muladd(g, rho, coeffs[i])
    end

    pref = rho^sigma
    prefp = sigma == 0 ? zero(pref) : sigma * pref / rho
    prefpp = sigma == 0 ? zero(pref) : (sigma - 1) * prefp / rho
    y = pref * g
    yp = prefp * g + pref * gu
    ypp = prefpp * g + 2 * prefp * gu + pref * guu
    return (y, yp, ypp)
end

function contour_phase(omega)
    ω = complex(omega)
    iszero(ω) && error("omega must be nonzero for the rotated contour.")
    iszero(imag(ω)) && return 1.0 + 0.0im
    return ω / abs(ω)
end

function contour_q1(rho, s, epsilon, tau, kappa, z)
    return z * (1 - 1im * epsilon - s - 1im * tau) +
           2 * (1 - 1im * epsilon * kappa - 1im * tau) * rho -
           2im * epsilon * kappa * rho^2 / z
end

function contour_q2(rho, s, epsilon, tau, kappa, lambda, z)
    return epsilon^2 - lambda - 1im * epsilon * kappa * (1 - 2 * s) -
           s * (s + 1) - tau * (1im + tau) +
           2 * epsilon * kappa * (-1im + 1im * s + epsilon - tau) * rho / z
end

function contour_residual(rho, Y, Yp, Ypp, s, epsilon, tau, kappa, lambda, z)
    A = rho * (z + rho)
    B = contour_q1(rho, s, epsilon, tau, kappa, z)
    C = contour_q2(rho, s, epsilon, tau, kappa, lambda, z)
    return abs(A * Ypp + B * Yp + C * Y)/max(abs(A * Ypp), abs(B * Yp), abs(C * Y))
end

function contour_PlusTwo_to_MinusTwo_Inf(Pplus, m, a, omega, lambda, z)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    (coe0, coe1, coe2, coe3) = coe1_funcs(2, a, omega, m, lambda)
    Tplus, Tpplus, _ = Tx(2, epsilon, kappa, tau)
    Tminus, Tpminus, TpTminus = Tx(-2, epsilon, kappa, tau)
    f4p = 32 * kappa^(5 - 2im * omega) * exp(-1im * omega * (1 + kappa))
    f4m = kappa^(-3 - 2im * omega) * exp(-1im * omega * (1 + kappa)) / 8

    function Rplus(rho)
        x = - rho / z
        r = x_to_r(x, kappa)
        P, Pρ, Pρρ, error = Pplus(rho)
        Px = - z * Pρ
        R = Tplus(x) * P / f4p
        Rp = - (Tpplus(x) * P + Tplus(x) * Px) / (2 * kappa) / f4p
        Rpp = d2R(r, R, Rp, 2, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end

    function Pminus(rho)
        x = - rho / z
        r = x_to_r(x, kappa)
        R, Rp, Rpp, error = Rplus(rho)
        Rm = (coe0(r) * R + coe1(r) * Rp) / epsilon^4
        Rmp = (coe2(r) * R + coe3(r) * Rp) / epsilon^4
        P = (Rm / Tminus(x)) * f4m
        Pp = (-2 * kappa * Rmp / Tminus(x) - TpTminus(x) * Rm / Tminus(x)) * f4m
        Pρ = - z * Pp
        Ppp = d2P(x, P, Pp, -2, epsilon, kappa, tau, lambda + 4)
        Pρρ = Ppp * (-z)^2
        return (P, Pρ, Pρρ, error)
    end

    return Rplus, Pminus
end

function contour_PlusTwo_to_MinusTwo_Hor(Pplus, m, a, omega, lambda, z)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    rp = 1 + kappa
    p = kappa * (epsilon + tau) / (2 * rp)
    err_scale = kappa / (4 * rp)
    (coe0, coe1, coe2, coe3) = coe1_funcs(2, a, omega, m, lambda)
    Tplus, Tpplus, _ = Tx(2, epsilon, kappa, tau)
    Tminus, Tpminus, TpTminus = Tx(-2, epsilon, kappa, tau)
    f1p = 16 * kappa^(2im * p + 4) * exp(1im * p * rp)
    f1m = kappa^(2im * p - 4) * exp(1im * p * rp) / 16
    Q = lambda + 6
    Csq = (Q^2 + 4*a*m*omega-4*a^2*omega^2)*((Q-2)^2+
            36*a*m*omega-36*a^2*omega^2)+(2*Q-1)*(96*
            a^2*omega^2-48*a*omega*m)+144*omega^2*(1-a^2)
    horizon_factor = 256 * rp^4 * (p + 4im * err_scale) * p * (p^2 + 4 * err_scale^2)

    function Rplus(rho)
        x = - rho / z
        r = x_to_r(x, kappa)
        P, Pρ, Pρρ, error = Pplus(rho)
        Px = - z * Pρ
        R = Tplus(x) * P / f1p
        Rp = - (Tpplus(x) * P + Tplus(x) * Px) / (2 * kappa) / f1p
        Rpp = d2R(r, R, Rp, 2, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end

    function Pminus(rho)
        x = - rho / z
        r = x_to_r(x, kappa)
        R, Rp, Rpp, error = Rplus(rho)
        Rm = (coe0(r) * R + coe1(r) * Rp) * horizon_factor / Csq
        Rmp = (coe2(r) * R + coe3(r) * Rp) * horizon_factor / Csq
        P = (Rm / Tminus(x)) * f1m
        Pp = (- 2 * kappa * Rmp / Tminus(x) - TpTminus(x) * Rm / Tminus(x)) * f1m
        Pρ = - z * Pp
        Ppp = d2P(x, P, Pp, -2, epsilon, kappa, tau, lambda + 4)
        Pρρ = Ppp * (-z)^2
        return (P, Pρ, Pρρ, error)
    end

    return Rplus, Pminus
end

function contour_MinusTwo_to_PlusTwo_Inf(Pminus, m, a, omega, lambda, z)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    (coe0, coe1, coe2, coe3) = coe2_funcs(-2, a, omega, m, lambda)
    Tminus, Tpminus, TpTminus = Tx(-2, epsilon, kappa, tau)
    Tplus, Tpplus, TpTplus = Tx(2, epsilon, kappa, tau)
    f3 = 2.0 * kappa^(2im * omega + 1) * exp(1im * omega * (1 + kappa))

    function Rminus(rho)
        x = - rho / z
        r = x_to_r(x, kappa)
        P, Pρ, Pρρ, error = Pminus(rho)
        Px = - z * P
        R = Tminus(x) * P / f3
        Rp = - (Tpminus(x) * P + Tminus(x) * Px) / (2 * kappa) / f3
        Rpp = d2R(r, R, Rp, -2, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end

    function Pplus(rho)
        x = - rho / z
        r = x_to_r(x, kappa)
        R, Rp, Rpp, error = Rminus(rho)
        Rpl = (coe0(r) * R + coe1(r) * Rp) / epsilon^4
        Rplp = (coe2(r) * R + coe3(r) * Rp) / epsilon^4
        P = (Rpl / Tplus(x)) * f3
        Pp = (-2 * kappa * Rplp / Tplus(x) - TpTplus(x) * Rpl / Tplus(x)) * f3
        Pρ = - z * Pp
        Ppp = d2P(x, P, Pp, 2, epsilon, kappa, tau, lambda - 4)
        Pρρ = Ppp * (-z)^2
        return (P, Pρ, Pρρ, error)
    end

    return Rminus, Pplus
end

function contour_MinusTwo_to_PlusTwo_Hor(Pminus, m, a, omega, lambda, z)
    kappa = sqrt(1 - a^2)
    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    rp = 1 + kappa
    p = kappa * (epsilon + tau) / (2 * rp)
    err_scale = kappa / (4 * rp)
    (coe0, coe1, coe2, coe3) = coe2_funcs(-2, a, omega, m, lambda)
    Tminus, Tpminus, TpTminus = Tx(-2, epsilon, kappa, tau)
    Tplus, Tpplus, TpTplus = Tx(2, epsilon, kappa, tau)
    f2 = kappa^(- 2im * p) * exp(- 1im * p * rp)
    Q = lambda + 2
    Csq = (Q^2 + 4*a*m*omega-4*a^2*omega^2)*((Q-2)^2+
            36*a*m*omega-36*a^2*omega^2)+(2*Q-1)*(96*
            a^2*omega^2-48*a*omega*m)+144*omega^2*(1-a^2)
    horizon_factor = 256 * rp^4 * (p - 4im * err_scale) * p * (p^2 + 4 * err_scale^2)

    function Rminus(rho)
        x = - rho / z
        r = x_to_r(x, kappa)
        P, Pρ, Pρρ, error = Pminus(rho)
        Px = - z * P
        R = Tminus(x) * P / f2
        Rp = - (Tpminus(x) * P + Tminus(x) * Px) / (2 * kappa) / f2
        Rpp = d2R(r, R, Rp, -2, a, omega, m, lambda)
        return (R, Rp, Rpp, error)
    end

    function Pplus(rho)
        x = - rho / z
        r = x_to_r(x, kappa)
        R, Rp, Rpp, error = Rminus(rho)
        Rpl = (coe0(r) * R + coe1(r) * Rp) * horizon_factor / Csq
        Rplp = (coe2(r) * R + coe3(r) * Rp) * horizon_factor / Csq
        P = (Rpl / Tplus(x)) * f2
        Pp = (- 2 * kappa * Rplp / Tplus(x) - TpTplus(x) * Rpl / Tplus(x)) * f2
        Pρ = - z * Pp
        Ppp = d2P(x, P, Pp, 2, epsilon, kappa, tau, lambda - 4)
        Pρρ = Ppp * (-z)^2
        return (P, Pρ, Pρρ, error)
    end

    return Rminus, Pplus
end

end
