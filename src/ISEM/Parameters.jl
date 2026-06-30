module Parameters

using SpinWeightedSpheroidalHarmonics

export TeukolskyParameters, MSEParameters, isem_parameters
export RhoTeukolskyParameters, RhoISEMParameters, rho_parameters

struct TeukolskyParameters
    s::Int
    l::Int
    m::Int
    a
    omega
end

const MSEParameters = TeukolskyParameters
const RhoTeukolskyParameters = TeukolskyParameters

function Base.show(io::IO, ::MIME"text/plain", p::TeukolskyParameters)
    println(io, "TeukolskyParameters:")
    println(io, "  s      = ", p.s,)
    println(io, "  l      = ", p.l,)
    println(io, "  m      = ", p.m,)
    println(io, "  a      = ", p.a,)
    println(io, "  omega  = ", p.omega,)
end

struct ISEMParameters
    s::Int
    kappa
    epsilon
    tau
    lambda
end

struct RhoISEMParameters
    s
    z
    kappa
    epsilon
    tau
    lambda
end

function Base.show(io::IO, ::MIME"text/plain", para::ISEMParameters)
    println(io, "ISEMParameters:")
    println(io, "  s       = ", para.s,)
    println(io, "  kappa   = ", para.kappa,)
    println(io, "  epsilon = ", para.epsilon,)
    println(io, "  tau     = ", para.tau,)
    println(io, "  lambda  = ", para.lambda,)
end

@inline function _core_parameters(s, l, m, a, omega)
    lambda = SpinWeightedSpheroidalHarmonics.spin_weighted_spheroidal_eigenvalue(s, l, m, a * omega)
    kappa = sqrt(1 - a^2)
    return lambda, kappa
end

function isem_parameters(params::TeukolskyParameters)
    s = params.s
    l = params.l
    m = params.m
    a = params.a
    omega = params.omega

    lambda, kappa = _core_parameters(s, l, m, a, omega)
    if isapprox(kappa, 0.0, atol=1e-10)
        error("Parameter 'a' is too close to 1, leading to kappa ≈ 0. Special handling for extreme Kerr black hole (a=1) is required for tau calculation.")
    end

    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa

    return ISEMParameters(s, kappa, epsilon, tau, lambda)
end

function rho_parameters(params::TeukolskyParameters)
    s = params.s
    l = params.l
    m = params.m
    a = params.a
    omega = params.omega

    lambda, kappa = _core_parameters(s, l, m, a, omega)
    isapprox(kappa, 0.0, atol=1e-10) && error("Parameter 'a' is too close to 1, leading to kappa ≈ 0.")

    epsilon = 2 * omega
    tau = (epsilon - m * a) / kappa
    z = iszero(imag(complex(omega))) ? 1.0 + 0.0im : complex(omega) / abs(complex(omega))

    return RhoISEMParameters(s, z, kappa, epsilon, tau, lambda)
end

end
