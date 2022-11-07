module Transformation

using ..Kerr
using ..Potentials

export alpha, alpha_prime, beta, beta_prime, eta, eta_prime

const I = 1im # Mathematica being Mathematica

function alpha(s::Int, m::Int, a, omega, lambda, r)
    throw(DomainError(s, "Not implemented yet"))
end

function beta(s::Int, m::Int, a, omega, lambda, r)
    throw(DomainError(s, "Not implemented yet"))
end

function alpha_prime(s::Int, m::Int, a, omega, lambda, r)
    throw(DomainError(s, "Not implemented yet"))
end

function beta_prime(s::Int, m::Int, a, omega, lambda, r)
    throw(DomainError(s, "Not implemented yet"))
end

function eta(s::Int, m::Int, a, omega, lambda, r)
    throw(DomainError(s, "Not implemented yet"))
end

function eta_prime(s::Int, m::Int, a, omega, lambda, r)
    throw(DomainError(s, "Not implemented yet"))
end

end