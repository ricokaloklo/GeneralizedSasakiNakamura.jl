function Delta(a, r)
    r^2 - 2*r + a^2
end

function r_plus(a)
    1 + sqrt(1 - a^2)
end

function r_minus(a)
    1 - sqrt(1 - a^2)
end

function omega_horizon(a)
    rp = r_plus(a)
    return a/(2*1*rp)
end

function K(m::Int, a, omega, r)
    (r^2 + a^2)*omega - m*a
end

function radial_Teukolsky_equation(s, m, a, omega, lambda, r, R, dRdr, d2Rdr2)
    ((1 + s)*Delta(a ,r)^s*dRdr(r)*(2*(r-1)) + Delta(a, r)^(1 + s)*d2Rdr2(r))*(Delta(a,r)^(-s)) - VT(s, m, a, omega, lambda, r)*R(r)
end