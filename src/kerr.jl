function Delta(M, a, r)
    r^2 - 2*M*r + a^2
end

function r_plus(M, a)
    M + sqrt(M^2 - a^2)
end

function r_minus(M, a)
    M - sqrt(M^2 - a^2)
end

function omega_horizon(M, a)
    rp = r_plus(M, a)
    return a/(2*M*rp)
end

function K(m::Int, a, omega)
    (r^2 + a^2)*omega - m*a
end