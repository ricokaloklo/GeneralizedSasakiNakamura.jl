module ConvolutionIntegrals

using LinearAlgebra
using SpinWeightedSpheroidalHarmonics
using KerrGeodesics

using ..GeneralizedSasakiNakamura: BoundaryCondition, IN, UP
using ..ISEM: Y_radial
using ..Coordinates: rstar_from_r
using ..SolutionsY: Y_solution
using ..GridSampling
import ..GridSampling: _isem_gsn_incidence_amplitude

export convolution_integral_trapezoidal, convolution_integral_levin
export convolution_integral_trapezoidal_isem, convolution_integral_levin_isem, convolution_integral_circular_equatorial_isem
export EccentricFluxCache, InclinedFluxCache

function horizon_factor(ω, a, m)
    rp = 1 + sqrt(1 - a^2)
    κ = ω - a * m / (2 * rp)
    ϵ0 = sqrt(1 - a^2) / (4 * rp)
    return ω / (κ * (2 * rp)^3 * (κ^2 + 4 * ϵ0^2) * 64pi)
end

const _INHOMOGENEOUS_RESONANCE_TOL = 1e-12

@inline _omega_horizon(a, m) = m * a / (2 * (1 + sqrt(1 - a^2)))
@inline _skip_infinity_mode(ω; tol = _INHOMOGENEOUS_RESONANCE_TOL) = abs(ω) < tol
@inline _skip_horizon_mode(a, m, ω; tol = _INHOMOGENEOUS_RESONANCE_TOL) = abs(ω - _omega_horizon(a, m)) < tol
@inline _skip_radiative_mode(s, a, m, ω; tol = _INHOMOGENEOUS_RESONANCE_TOL) =
    (s == -2 && _skip_infinity_mode(ω; tol = tol)) || (s == 2 && _skip_horizon_mode(a, m, ω; tol = tol))

function _zero_radiative_mode(ω, trajectory; reason = "")
    return Dict(
        "Amplitude" => 0.0 + 0.0im,
        "omega" => ω,
        "EnergyFlux" => 0.0,
        "AngularMomentumFlux" => 0.0,
        "CarterConstantFlux" => 0.0,
        "Trajectory" => trajectory,
        "YSolution" => nothing,
        "SWSH" => nothing,
        "SkippedMode" => true,
        "SkippedReason" => reason,
    )
end

const _generic_trapezoidal_cache = Dict{Tuple, Any}()
const _eccentric_trapezoidal_cache = Dict{Tuple, Any}()
const _inclined_trapezoidal_cache = Dict{Tuple, Any}()
const _generic_trapezoidal_last_key = Ref{Any}(nothing)
const _generic_trapezoidal_last_result = Ref{Any}(nothing)
const _generic_levin_last_key = Ref{Any}(nothing)
const _generic_levin_last_result = Ref{Any}(nothing)
const _eccentric_levin_last_key = Ref{Any}(nothing)
const _eccentric_levin_last_result = Ref{Any}(nothing)
const _inclined_levin_last_key = Ref{Any}(nothing)
const _inclined_levin_last_result = Ref{Any}(nothing)

mutable struct EccentricFluxCache
    samples::Dict{Int, Dict}
end

mutable struct InclinedFluxCache
    samples::Dict{Int, Dict}
end

EccentricFluxCache() = EccentricFluxCache(Dict{Int, Dict}())
InclinedFluxCache() = InclinedFluxCache(Dict{Int, Dict}())

function _cached_eccentric_sample!(cache::EccentricFluxCache, KG_master::Dict, N::Int)
    return get!(cache.samples, N) do
        subsample_eccentric_sample(KG_master, N)
    end
end

function _cached_inclined_sample!(cache::InclinedFluxCache, KG_master::Dict, K::Int)
    return get!(cache.samples, K) do
        subsample_inclined_sample(KG_master, K)
    end
end

function _generic_trapezoidal_context(a, p, e, x, s, l, m, n, k, N_sample, K_sample)
    key = (a, p, e, x, s, l, m, n, k, N_sample, K_sample)
    return get!(_generic_trapezoidal_cache, key) do
        KG = Kerr_Geodesics(a, p, e, x)
        Frequencies = KG["Frequencies"]
        Γ = Frequencies["ϒt"]
        ϒr = Frequencies["ϒr"]
        ϒθ = Frequencies["ϒθ"]
        ϒφ = Frequencies["ϒϕ"]
        omega = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
        KG_samp = kerr_geo_generic_sample(KG, N_sample, K_sample)
        SH = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Ysol = _isem_y_solution(s, l, m, a, omega)
        Ysamp = s == 2 ? y_sample_p2_isem(Ysol, KG_samp) : y_sample_m2_isem(Ysol, KG_samp)
        SHsamp = swsh_sample(SH, KG_samp)
        carter_samp = carter_ingredients_sample(KG_samp, a, m, omega)
        carter_factor = trapezoidal_1d_integral(carter_samp)
        return (KG = KG, KG_samp = KG_samp, omega = omega, Γ = Γ, ϒθ = ϒθ,
                Ysol = Ysol, Ysamp = Ysamp, SH = SH, SHsamp = SHsamp,
                carter_factor = carter_factor)
    end
end

function _eccentric_trapezoidal_context(a, p, e, s, l, m, n, N_sample)
    key = (a, p, e, s, l, m, n, N_sample)
    return get!(_eccentric_trapezoidal_cache, key) do
        KG = Kerr_Geodesics(a, p, e, 1.0)
        Frequencies = KG["Frequencies"]
        Γ = Frequencies["ϒt"]
        ϒθ = Frequencies["ϒθ"]
        ϒφ = Frequencies["ϒϕ"]
        ϒr = Frequencies["ϒr"]
        omega = (m * ϒφ + n * ϒr) / Γ
        KG_samp = kerr_geo_eccentric_sample(KG, N_sample)
        SH = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Ysol = _isem_y_solution(s, l, m, a, omega)
        Ysamp = s == 2 ? y_sample_p2_isem(Ysol, KG_samp) : y_sample_m2_isem(Ysol, KG_samp)
        SHsamp = swsh_sample(SH, KG_samp)
        return (KG = KG, KG_samp = KG_samp, omega = omega, Γ = Γ, ϒθ = ϒθ, Ysol = Ysol, Ysamp = Ysamp, SH = SH, SHsamp = SHsamp)
    end
end

function _inclined_trapezoidal_context(KG_sample::Dict, s::Int, l::Int, m::Int, k::Int, K_interval::Int64, cache = nothing)
    key = (KG_sample["a"], KG_sample["p"], KG_sample["x"], s, l, m, k, K_interval)
    return get!(_inclined_trapezoidal_cache, key) do
        Frequencies = KG_sample["Frequencies"]
        Γ = Frequencies["ϒt"]
        ϒθ = Frequencies["ϒθ"]
        ϒφ = Frequencies["ϒϕ"]
        a = KG_sample["a"]
        ω = (m * ϒφ + k * ϒθ) / Γ
        KG_samp = cache isa InclinedFluxCache ? _cached_inclined_sample!(cache, KG_sample, K_interval) : subsample_inclined_sample(KG_sample, K_interval)
        Ysol = _isem_y_solution(s, l, m, a, ω)
        SH = spin_weighted_spheroidal_harmonic(s, l, m, a * ω; method = "jacobi")
        Y, Yp, X, _ = Ysol.Y_solution(KG_sample["p"])
        Ydic = Dict(
            "params" => (s = s, l = l, m = m, a = a, omega = ω, lambda = SH.lambda),
            (s == 2 ? "Cinc" : "Binc") => _isem_gsn_incidence_amplitude(Ysol),
            "Y" => Y,
            "Yp" => Yp,
            "X" => X,
        )
        SHsamp = swsh_sample(SH, KG_samp)
        return (KG_samp = KG_samp, Ysol = Ysol, Ydic = Ydic, SH = SH, SHsamp = SHsamp, omega = ω, a = a, m = m, Trajectory = KG_sample["Trajectory"])
    end
end

function _inclined_flux_from_sample(KG_samp::Dict, Ysol, Ydic::Dict, SH_samp::Dict, s::Int, a, ω, m, k::Int)
    carter_factor = trapezoidal_1d_integral(carter_ingredients_sample(KG_samp, a, m, ω))
    if s == 2
        inte = integrand_inclined_sample_p2(KG_samp, Ydic, SH_samp, k)
        integral = trapezoidal_1d_integral(inte)
        hf = horizon_factor(ω, a, m)
        return Dict(
            "Amplitude" => integral,
            "omega" => ω,
            "EnergyFlux" => hf * abs2(integral),
            "AngularMomentumFlux" => hf * m * abs2(integral) / ω,
            "CarterConstantFlux" => 2 * hf * abs2(integral) * (carter_factor + k * KG_samp["Frequencies"]["ϒθ"]) / ω,
            "Trajectory" => KG_samp,
            "YSolution" => Ysol,
            "SWSH" => SH_samp,
        )
    elseif s == -2
        inte = integrand_inclined_sample_m2(KG_samp, Ydic, SH_samp, k)
        integral = trapezoidal_1d_integral(inte)
        return Dict(
            "Amplitude" => integral,
            "omega" => ω,
            "EnergyFlux" => abs2(integral) / (4.0pi * ω^2),
            "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * ω^3),
            "CarterConstantFlux" => abs2(integral) * (carter_factor + k * KG_samp["Frequencies"]["ϒθ"]) / (2.0pi * ω^3),
            "Trajectory" => KG_samp,
            "YSolution" => Ysol,
            "SWSH" => SH_samp,
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_generic_trapezoidal_isem(a, p, e, x, s, l, m, n, k, N_sample, K_sample)
    key = (a, p, e, x, s, l, m, n, k, N_sample, K_sample)
    if _generic_trapezoidal_last_key[] === key
        return _generic_trapezoidal_last_result[]
    end
    KG = Kerr_Geodesics(a, p, e, x)
    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    omega = (m * Frequencies["ϒϕ"] + n * Frequencies["ϒr"] + k * Frequencies["ϒθ"]) / Γ
    if _skip_radiative_mode(s, a, m, omega)
        result = _zero_radiative_mode(omega, KG; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        _generic_trapezoidal_last_key[] = key
        _generic_trapezoidal_last_result[] = result
        return result
    end
    ctx = _generic_trapezoidal_context(a, p, e, x, s, l, m, n, k, N_sample, K_sample)

    result = if s == 2
        inte_in_samp_p2 = integrand_generic_sample_p2(ctx.KG_samp, ctx.Ysamp, ctx.SHsamp, n, k)
        integral_p2 = trapezoidal_2d_integral(inte_in_samp_p2)
        hf = horizon_factor(ctx.omega, a, m)
        Dict(
            "Amplitude" => integral_p2,
            "omega" => ctx.omega,
            "EnergyFlux" => hf * abs2(integral_p2),
            "AngularMomentumFlux" => hf * m * abs2(integral_p2) / ctx.omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_p2) * (ctx.carter_factor + k * ctx.ϒθ) / ctx.omega,
            "Trajectory" => ctx.KG,
            "YSolution" => ctx.Ysol,
            "SWSH" => ctx.SH
        )
    elseif s == -2
        inte_in_samp_m2 = integrand_generic_sample_m2(ctx.KG_samp, ctx.Ysamp, ctx.SHsamp, n, k)
        integral_m2 = trapezoidal_2d_integral(inte_in_samp_m2)
        Dict(
            "Amplitude" => integral_m2,
            "omega" => ctx.omega,
            "EnergyFlux" => abs2(integral_m2)/(4.0pi*ctx.omega^2),
            "AngularMomentumFlux" => m * abs2(integral_m2) / (4.0pi * ctx.omega^3),
            "CarterConstantFlux" => abs2(integral_m2) * (ctx.carter_factor + k * ctx.ϒθ) / (2.0pi*ctx.omega^3),
            "Trajectory" => ctx.KG,
            "YSolution" => ctx.Ysol,
            "SWSH" => ctx.SH
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
    _generic_trapezoidal_last_key[] = key
    _generic_trapezoidal_last_result[] = result
    return result
end

function convolution_integral_generic_levin_isem(a, p, e, x, s, l, m, n, k, N_sample, K_sample)
    key = (a, p, e, x, s, l, m, n, k, N_sample, K_sample)
    if _generic_levin_last_key[] === key
        return _generic_levin_last_result[]
    end
    KG = Kerr_Geodesics(a, p, e, x)
    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
    if _skip_radiative_mode(s, a, m, omega)
        result = _zero_radiative_mode(omega, KG; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        _generic_levin_last_key[] = key
        _generic_levin_last_result[] = result
        return result
    end
    KG_samp = kerr_geo_generic_sample_cheby(KG, N_sample, K_sample)
    KG_trap = kerr_geo_generic_sample(KG, N_sample, K_sample)
    carter_samp = carter_ingredients_sample(KG_trap, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    result = if s == 2
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        Yup_soln = _isem_y_solution(s, l, m, a, omega)
        Yup_samp = y_sample_p2_isem(Yup_soln, KG_samp)
        Jpp_grid_up, Jpm_grid_up, Jmp_grid_up, Jmm_grid_up, drphase_up, dθphase_up, rphaseL_up, rphaseR_up, θphaseL_up, θphaseR_up, prefactor_up = integrand_generic_sample_cheby_p2(KG_samp, Yup_samp, SH_p2_samp, n, k)
        integralpp_up = levin_2d_integral(Jpp_grid_up, drphase_up, dθphase_up, rphaseL_up, rphaseR_up, θphaseL_up, θphaseR_up)
        integralpm_up = levin_2d_integral(Jpm_grid_up, drphase_up, - dθphase_up, rphaseL_up, rphaseR_up, - θphaseL_up, - θphaseR_up)
        integralmp_up = levin_2d_integral(Jmp_grid_up, - drphase_up, dθphase_up, - rphaseL_up, - rphaseR_up, θphaseL_up, θphaseR_up)
        integralmm_up = levin_2d_integral(Jmm_grid_up, - drphase_up, - dθphase_up, - rphaseL_up, - rphaseR_up, - θphaseL_up, - θphaseR_up)
        integral_up = (integralpp_up + integralpm_up + integralmp_up + integralmm_up) * prefactor_up
        hf = horizon_factor(omega, a, m)
        Dict(
            "Amplitude" => integral_up,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_up),
            "AngularMomentumFlux" => hf * m * abs2(integral_up) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_up) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yin_soln = _isem_y_solution(s, l, m, a, omega)
        Yin_samp = y_sample_m2_isem(Yin_soln, KG_samp)
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        Jpp_grid_in, Jpm_grid_in, Jmp_grid_in, Jmm_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_generic_sample_cheby_m2(KG_samp, Yin_samp, SH_m2_samp, n, k)
        integralpp_in = levin_2d_integral(Jpp_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in)
        integralpm_in = levin_2d_integral(Jpm_grid_in, drphase_in, - dθphase_in, rphaseL_in, rphaseR_in, - θphaseL_in, - θphaseR_in)
        integralmp_in = levin_2d_integral(Jmp_grid_in, - drphase_in, dθphase_in, - rphaseL_in, - rphaseR_in, θphaseL_in, θphaseR_in)
        integralmm_in = levin_2d_integral(Jmm_grid_in, - drphase_in, - dθphase_in, - rphaseL_in, - rphaseR_in, - θphaseL_in, - θphaseR_in)
        integral_in = (integralpp_in + integralpm_in + integralmp_in + integralmm_in) * prefactor_in
        Dict(
            "Amplitude" => integral_in,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_in)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_in) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_in) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
    _generic_levin_last_key[] = key
    _generic_levin_last_result[] = result
    return result
end

@inline function _isem_y_solution(s, l, m, a, omega)
    if s == -2
        return Y_radial(s, l, m, a, omega, IN)
    elseif s == 2
        return Y_radial(s, l, m, a, omega, UP)
    else
        error("ISEM Y interfaces currently support only s = -2 and s = 2.")
    end
end

@inline function _generic_isem_y_solution(s, l, m, a, omega)
    return _isem_y_solution(s, l, m, a, omega)
end

function _generic_flux_from_sample(KG_samp::Dict, Y_samp::Dict, SH_samp::Dict, s::Int, l::Int, m::Int, n::Int, k::Int, a::Float64, ω::Float64, ϒθ::Float64)
    carter_factor = trapezoidal_1d_integral(carter_ingredients_sample(KG_samp, a, m, ω))
    if s == 2
        inte = GridSampling._generic_integrand_sample_p2_fast(KG_samp, Y_samp, SH_samp, n, k)
        integral = trapezoidal_2d_integral(inte)
        hf = horizon_factor(ω, a, m)
        return Dict(
            "Amplitude" => integral,
            "omega" => ω,
            "EnergyFlux" => hf * abs2(integral),
            "AngularMomentumFlux" => hf * m * abs2(integral) / ω,
            "CarterConstantFlux" => 2 * hf * abs2(integral) * (carter_factor + k * ϒθ) / ω,
            "Trajectory" => KG_samp,
            "YSolution" => Y_samp,
            "SWSH" => SH_samp,
        )
    elseif s == -2
        inte = GridSampling._generic_integrand_sample_m2_fast(KG_samp, Y_samp, SH_samp, n, k)
        integral = trapezoidal_2d_integral(inte)
        return Dict(
            "Amplitude" => integral,
            "omega" => ω,
            "EnergyFlux" => abs2(integral) / (4.0pi * ω^2),
            "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * ω^3),
            "CarterConstantFlux" => abs2(integral) * (carter_factor + k * ϒθ) / (2.0pi * ω^3),
            "Trajectory" => KG_samp,
            "YSolution" => Y_samp,
            "SWSH" => SH_samp,
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

mutable struct GenericM2FluxCache
    samples::Dict{Tuple{Int, Int}, Dict}
    mcaches::Dict{Tuple{Int, Int, Int}, Any}
end

GenericM2FluxCache() = GenericM2FluxCache(Dict{Tuple{Int, Int}, Dict}(), Dict{Tuple{Int, Int, Int}, Any}())
GenericFluxCache() = GenericM2FluxCache()

struct GenericM2MCache
    r::Vector{Float64}
    rs::Vector{Float64}
    theta::Vector{Float64}
    dtr::Vector{Float64}
    dttheta::Vector{Float64}
    dphir::Vector{Float64}
    dphitheta::Vector{Float64}
    cross::Vector{Function}
    initial_phases::NTuple{4, Float64}
    wr::Vector{Float64}
    wtheta::Vector{Float64}
    qr::Vector{Float64}
    qtheta::Vector{Float64}
    m_phase_r::Vector{ComplexF64}
    m_phase_theta::Vector{ComplexF64}
    m_log_phase::Vector{ComplexF64}
    Np::Vector{Float64}
    Nm::Vector{Float64}
    Np2::Vector{Float64}
    Nm2::Vector{Float64}
    Lp::Vector{Float64}
    Lm::Vector{Float64}
    Lp2::Vector{Float64}
    Lm2::Vector{Float64}
    st::Vector{Float64}
    ct::Vector{Float64}
    invst::Vector{Float64}
    invst2::Vector{Float64}
    termM::Vector{ComplexF64}
    uthetap::Vector{Float64}
    uthetam::Vector{Float64}
    rho::Matrix{ComplexF64}
    rhobar::Matrix{ComplexF64}
    invrho::Matrix{ComplexF64}
    rho_minus::Matrix{ComplexF64}
    rho_plus::Matrix{ComplexF64}
    radial_phase::Vector{ComplexF64}
    radial_tortoise::Vector{ComplexF64}
    theta_phase::Vector{ComplexF64}
    partial::Vector{ComplexF64}
    N::Int
    K::Int
    Gamma::Float64
    E::Float64
    Lz::Float64
    a::Float64
    m::Int
end

function _trap_weights_1d(N::Int, dq::Float64)
    w = fill(dq, N)
    w[1] *= 0.5
    w[end] *= 0.5
    return w
end

function _cached_generic_sample!(cache::GenericM2FluxCache, KG_master::Dict, N::Int, K::Int)
    key = (N, K)
    return get!(cache.samples, key) do
        GridSampling.subsample_generic_sample(KG_master, N, K)
    end
end

function _generic_m2_m_cache(KG_samp::Dict, m::Int)
    r = KG_samp["r"]::Vector{Float64}
    rs = KG_samp["rs"]::Vector{Float64}
    theta = KG_samp["θ"]::Vector{Float64}
    N = KG_samp["N_sample"]::Int
    K = KG_samp["K_sample"]::Int
    a = KG_samp["a"]::Float64
    E = KG_samp["E"]::Float64
    Lz = KG_samp["Lz"]::Float64
    dqr = pi / (N - 1)
    dqtheta = pi / (K - 1)
    wr = _trap_weights_1d(N, dqr)
    wtheta = _trap_weights_1d(K, dqtheta)
    qr = [(i - 1) * dqr for i in 1:N]
    qtheta = [(j - 1) * dqtheta for j in 1:K]
    rp = 1.0 + sqrt(1.0 - a^2)
    rm = 1.0 - sqrt(1.0 - a^2)
    inv2kappa = 1.0 / (2.0 * sqrt(1.0 - a^2))
    dphir = KG_samp["Δφr"]::Vector{Float64}
    dphitheta = KG_samp["Δφθ"]::Vector{Float64}
    m_phase_r = Vector{ComplexF64}(undef, N)
    m_log_phase = Vector{ComplexF64}(undef, N)
    Np = Vector{Float64}(undef, N)
    Nm = Vector{Float64}(undef, N)
    Np2 = Vector{Float64}(undef, N)
    Nm2 = Vector{Float64}(undef, N)
    Lp = Vector{Float64}(undef, N)
    Lm = Vector{Float64}(undef, N)
    Lp2 = Vector{Float64}(undef, N)
    Lm2 = Vector{Float64}(undef, N)
    urp = KG_samp["ur_fwd"]::Vector{Float64}
    urm = KG_samp["ur_rev"]::Vector{Float64}
    @inbounds for i in 1:N
        ri = r[i]
        r2 = ri * ri
        Delta = r2 - 2.0 * ri + a^2
        m_phase_r[i] = exp(-im * m * dphir[i])
        m_log_phase[i] = exp(-im * a * m * log((ri - rp) / (ri - rm)) * inv2kappa)
        numer = E * (r2 + a^2) - a * Lz
        Np[i] = (numer + urp[i]) / Delta
        Nm[i] = (numer + urm[i]) / Delta
        Np2[i] = Np[i] * Np[i]
        Nm2[i] = Nm[i] * Nm[i]
        Lp[i] = (numer - urp[i]) / Delta
        Lm[i] = (numer - urm[i]) / Delta
        Lp2[i] = Lp[i] * Lp[i]
        Lm2[i] = Lm[i] * Lm[i]
    end
    m_phase_theta = Vector{ComplexF64}(undef, K)
    st = Vector{Float64}(undef, K)
    ct = Vector{Float64}(undef, K)
    invst = Vector{Float64}(undef, K)
    invst2 = Vector{Float64}(undef, K)
    termM = Vector{ComplexF64}(undef, K)
    uthetap = KG_samp["uθ_fwd"]::Vector{Float64}
    uthetam = KG_samp["uθ_rev"]::Vector{Float64}
    @inbounds for j in 1:K
        stj = sin(theta[j])
        ctj = cos(theta[j])
        invstj = 1.0 / stj
        invst2j = invstj * invstj
        m_phase_theta[j] = exp(-im * m * dphitheta[j])
        st[j] = stj
        ct[j] = ctj
        invst[j] = invstj
        invst2[j] = invst2j
        termM[j] = im * stj * (a * E - Lz * invst2j)
    end
    rho = Matrix{ComplexF64}(undef, N, K)
    rhobar = Matrix{ComplexF64}(undef, N, K)
    invrho = Matrix{ComplexF64}(undef, N, K)
    rho_minus = Matrix{ComplexF64}(undef, N, K)
    rho_plus = Matrix{ComplexF64}(undef, N, K)
    @inbounds for j in 1:K
        ctj = ct[j]
        for i in 1:N
            rhoij = -1.0 / (r[i] - im * a * ctj)
            rhobarij = -1.0 / (r[i] + im * a * ctj)
            rho[i, j] = rhoij
            rhobar[i, j] = rhobarij
            invrho[i, j] = 1.0 / rhoij
            rho_minus[i, j] = rhoij - rhobarij
            rho_plus[i, j] = rhoij + rhobarij
        end
    end
    radial_phase = Vector{ComplexF64}(undef, N)
    radial_tortoise = Vector{ComplexF64}(undef, N)
    theta_phase = Vector{ComplexF64}(undef, K)
    partial = Vector{ComplexF64}(undef, Threads.maxthreadid())
    return GenericM2MCache(r, rs, theta, KG_samp["Δtr"], KG_samp["Δtθ"], dphir, dphitheta,
        KG_samp["CrossFunction"], KG_samp["initialPhases"], wr, wtheta, qr, qtheta,
        m_phase_r, m_phase_theta, m_log_phase, Np, Nm, Np2, Nm2, Lp, Lm, Lp2, Lm2, st, ct, invst, invst2,
        termM, uthetap, uthetam, rho, rhobar, invrho, rho_minus, rho_plus,
        radial_phase, radial_tortoise, theta_phase, partial, N, K, KG_samp["Γ"], E, Lz, a, m)
end

function _cached_m2_mcache!(cache::GenericM2FluxCache, KG_samp::Dict, m::Int)
    key = (m, KG_samp["N_sample"], KG_samp["K_sample"])
    return get!(cache.mcaches, key) do
        _generic_m2_m_cache(KG_samp, m)
    end
end

@inline function _relative_energy_change(new_energy, old_energy)
    return abs(new_energy - old_energy) / max(min(abs(new_energy), abs(old_energy)), eps(Float64))
end

@inline function _low_flux_cutoff(max_flux::Float64, tol::Float64)
    return min(0.1 * max(max_flux, eps(Float64)) * tol, 1e-16)
end

function _phase_vector_exp!(out::Vector{ComplexF64}, base_omega::Float64, freq_vec::Vector{Float64}, mode_index::Int, q_vec::Vector{Float64}, fixed_phase::Vector{ComplexF64})
    @inbounds for i in eachindex(q_vec)
        out[i] = exp(im * (base_omega * freq_vec[i] + mode_index * q_vec[i])) * fixed_phase[i]
    end
    return out
end

function _phase_vector_recurrence!(out::Vector{ComplexF64}, base_omega::Float64, freq_vec::Vector{Float64}, mode_index::Int, step::Float64, fixed_phase::Vector{ComplexF64})
    qphase = 1.0 + 0.0im
    qstep = exp(im * mode_index * step)
    @inbounds for i in eachindex(freq_vec)
        out[i] = exp(im * base_omega * freq_vec[i]) * qphase * fixed_phase[i]
        qphase *= qstep
    end
    return out
end

@inline function _generic_carter_factor(cache::GenericM2MCache, omega::Float64)
    total = 0.0
    a2ωE = cache.a * cache.a * omega * cache.E
    mLz = cache.m * cache.Lz
    @inbounds for j in 1:cache.K
        ct2 = cache.ct[j] * cache.ct[j]
        total += cache.wtheta[j] * (mLz * ct2 * cache.invst2[j] - a2ωE * ct2) / pi
    end
    return total
end

function _generic_m2_integral_threaded(cache::GenericM2MCache, Y_samp::Dict, SH_samp::Dict, omega::Float64, n::Int, k::Int; recurrence_phase::Bool = true)
    N = cache.N
    K = cache.K
    a = cache.a
    m = cache.m
    dqr = pi / (N - 1)
    dqtheta = pi / (K - 1)
    radial_phase = recurrence_phase ?
        _phase_vector_recurrence!(cache.radial_phase, omega, cache.dtr, n, dqr, cache.m_phase_r) :
        _phase_vector_exp!(cache.radial_phase, omega, cache.dtr, n, cache.qr, cache.m_phase_r)
    radial_tortoise = cache.radial_tortoise
    @inbounds for i in 1:N
        radial_tortoise[i] = exp(im * omega * cache.rs[i]) * cache.m_log_phase[i]
    end
    theta_phase = recurrence_phase ?
        _phase_vector_recurrence!(cache.theta_phase, omega, cache.dttheta, k, dqtheta, cache.m_phase_theta) :
        _phase_vector_exp!(cache.theta_phase, omega, cache.dttheta, k, cache.qtheta, cache.m_phase_theta)
    Binc = Y_samp["Binc"]::ComplexF64
    Y = Y_samp["Y"]::Vector{ComplexF64}
    Yp = Y_samp["Yp"]::Vector{ComplexF64}
    X = Y_samp["X"]::Vector{ComplexF64}
    S0 = SH_samp["S0"]::Vector{ComplexF64}
    S1 = SH_samp["S1"]::Vector{ComplexF64}
    S2 = SH_samp["S2"]::Vector{ComplexF64}
    dt_r, dt_theta, dphi_r, dphi_theta = cache.cross
    qt0, qr0, qtheta0, qphi0 = cache.initial_phases
    Xi0 = Float64(m * (dphi_r(qr0) + dphi_theta(qtheta0) - qphi0) -
                  omega * (dt_r(qr0) + dt_theta(qtheta0) - qt0) -
                  k * qtheta0 - n * qr0)
    prefactor = 4im * pi * omega * exp(im * Xi0) / (cache.Gamma * Binc)
    inv4pi2 = 1.0 / (4.0 * pi^2)
    total = 0.0 + 0.0im
    @inbounds for j in 1:K
        total_j = 0.0 + 0.0im
        stj = cache.st[j]
        invstj = cache.invst[j]
        invst2j = cache.invst2[j]
        ctj = cache.ct[j]
        L1j = -m * invstj + a * omega * stj + ctj * invstj
        L2j = -m * invstj + a * omega * stj + 2.0 * ctj * invstj
        L2pj = m * ctj * invst2j + a * omega * ctj - 2.0 * invst2j
        S0j = S0[j]
        S1j = S1[j]
        S2j = S2[j]
        L2S = S1j + L2j * S0j
        L1Sp = S2j + L1j * S1j
        L1L2S = L1Sp + L2pj * S0j + L2j * S1j + L1j * L2j * S0j
        theta_const = 3im * a * ctj * S0j + 2im * a * stj * S1j - 1im * a * stj * L2j * S0j
        theta_L1 = 3im * a * stj * S0j
        exp_theta_p = theta_phase[j]
        exp_theta_m = conj(theta_phase[j])
        Mpj = cache.termM[j] + cache.uthetap[j]
        Mmj = cache.termM[j] + cache.uthetam[j]
        Mp2 = Mpj * Mpj
        Mm2 = Mmj * Mmj
        @inbounds for i in 1:N
            r = cache.r[i]
            r2 = r * r
            rhoi = cache.rho[i, j]
            rhobari = cache.rhobar[i, j]
            Yi = Y[i]
            Ypi = Yp[i]
            Xi = X[i]
            phase_r = radial_tortoise[i]
            L1pL2pS = cache.invrho[i, j] * L1L2S + theta_const + theta_L1 * L1j
            Wnn = rhobari * L1pL2pS * r2 * Yi * phase_r / 2.0
            termA = L2S + im * a * stj * cache.rho_minus[i, j] * S0j
            termB = 2.0 * Yi + r * Ypi
            termC = L2S * cache.rho_plus[i, j] * r * Yi
            Wnm = (-r * rhobari * (termA * termB + termC) / rhoi) * phase_r
            termD = Xi / (2.0 * sqrt(r2 + a^2))
            termE = (Yi + 2.0 * r * Ypi) * phase_r
            termF = rhoi * r * (2.0 * Yi + r * Ypi) * phase_r
            Wmm = S0j * rhobari * (termD + termE + termF) / rhoi
            Np = cache.Np[i]
            Nm = cache.Nm[i]
            exp_rp = radial_phase[i]
            exp_rm = conj(radial_phase[i])
            Jpp = (Wnn * cache.Np2[i] + Wnm * Np * Mpj + Wmm * Mp2) * exp_rp * exp_theta_p * inv4pi2
            Jpm = (Wnn * cache.Np2[i] + Wnm * Np * Mmj + Wmm * Mm2) * exp_rp * exp_theta_m * inv4pi2
            Jmp = (Wnn * cache.Nm2[i] + Wnm * Nm * Mpj + Wmm * Mp2) * exp_rm * exp_theta_p * inv4pi2
            Jmm = (Wnn * cache.Nm2[i] + Wnm * Nm * Mmj + Wmm * Mm2) * exp_rm * exp_theta_m * inv4pi2
            total_j += cache.wr[i] * (Jpp + Jpm + Jmp + Jmm)
        end
        total += cache.wtheta[j] * total_j
    end
    return prefactor * total
end

function _generic_p2_integral_threaded(cache::GenericM2MCache, Y_samp::Dict, SH_samp::Dict, omega::Float64, lambda::Float64, n::Int, k::Int; recurrence_phase::Bool = true)
    N = cache.N
    K = cache.K
    a = cache.a
    m = cache.m
    dqr = pi / (N - 1)
    dqtheta = pi / (K - 1)
    radial_phase = recurrence_phase ?
        _phase_vector_recurrence!(cache.radial_phase, omega, cache.dtr, n, dqr, cache.m_phase_r) :
        _phase_vector_exp!(cache.radial_phase, omega, cache.dtr, n, cache.qr, cache.m_phase_r)
    radial_tortoise = cache.radial_tortoise
    @inbounds for i in 1:N
        radial_tortoise[i] = exp(-im * omega * cache.rs[i]) * conj(cache.m_log_phase[i])
    end
    theta_phase = recurrence_phase ?
        _phase_vector_recurrence!(cache.theta_phase, omega, cache.dttheta, k, dqtheta, cache.m_phase_theta) :
        _phase_vector_exp!(cache.theta_phase, omega, cache.dttheta, k, cache.qtheta, cache.m_phase_theta)
    Cinc = Y_samp["Cinc"]::ComplexF64
    Y = Y_samp["Y"]::Vector{ComplexF64}
    Yp = Y_samp["Yp"]::Vector{ComplexF64}
    X = Y_samp["X"]::Vector{ComplexF64}
    S0 = SH_samp["S0"]::Vector{ComplexF64}
    S1 = SH_samp["S1"]::Vector{ComplexF64}
    S2 = SH_samp["S2"]::Vector{ComplexF64}
    dt_r, dt_theta, dphi_r, dphi_theta = cache.cross
    qt0, qr0, qtheta0, qphi0 = cache.initial_phases
    rp = 1.0 + sqrt(1.0 - a^2)
    kappa = omega - a * m / (2.0 * rp)
    c0 = 24 + 12im * omega + lambda * (10 + lambda) - 12 * a * omega * (a * omega - m)
    c1 = -32im * a * m - 8im * a * m * lambda + 8im * a^2 * omega * (1 + lambda)
    c2 = 12 * a^2 - 24im * a * m - 24 * a^2 * m^2 + 24im * a^2 * omega + 48 * a^3 * m * omega - 24 * a^4 * omega^2
    c3 = -24im * a^3 * (a * omega - m) - 24 * a^2
    c4 = 12 * a^4
    eta = c0 + c1 / rp + c2 / rp^2 + c3 / rp^3 + c4 / rp^4
    factor = 2 * sqrt(2) * rp^(3 / 2) * (((4 * omega) * (im - 4 * omega) - a * m * (im - 8 * omega) - a^2 * (m^2 + 2im * omega - 4 * omega^2)) * rp^2 +
        a^2 * (im - 4 * omega) * (a * m - 2 * omega) * rp) / (2 * rp^3 * (24 + 10 * lambda + lambda^2 + 12im * omega) -
        rp^2 * (8im * a * m * (11 + 2 * lambda + 6im * omega) + a^2 * (24 + 24 * m^2 + 10 * lambda + lambda^2 - 28im * omega - 16im * lambda * omega + 48 * omega^2)) +
        8im * a^3 * rp * (m * (7 + lambda - 6im * omega) - a * omega * (4 + lambda)) + 12 * a^5 * omega * (a * omega - 3 * m))
    Xi0 = Float64(m * (dphi_r(qr0) + dphi_theta(qtheta0) - qphi0) -
                  omega * (dt_r(qr0) + dt_theta(qtheta0) - qt0) -
                  k * qtheta0 - n * qr0)
    prefactor = -im * pi * factor * eta * exp(im * Xi0) / (cache.Gamma * kappa * Cinc)
    inv4pi2 = 1.0 / (4.0 * pi^2)
    total = 0.0 + 0.0im
    @inbounds for j in 1:K
        total_j = 0.0 + 0.0im
        stj = cache.st[j]
        invstj = cache.invst[j]
        invst2j = cache.invst2[j]
        ctj = cache.ct[j]
        L1j = m * invstj - a * omega * stj + ctj * invstj
        L2j = m * invstj - a * omega * stj + 2.0 * ctj * invstj
        L2pj = -m * ctj * invst2j - a * omega * ctj - 2.0 * invst2j
        S0j = S0[j]
        S1j = S1[j]
        S2j = S2[j]
        L2S = S1j + L2j * S0j
        L1Sp = S2j + L1j * S1j
        L1L2S = L1Sp + L2pj * S0j + L2j * S1j + L1j * L2j * S0j
        theta_const = 3im * a * ctj * S0j + 2im * a * stj * S1j - 1im * a * stj * L2j * S0j
        theta_L1 = 3im * a * stj * S0j
        exp_theta_p = theta_phase[j]
        exp_theta_m = conj(theta_phase[j])
        Mpj = -cache.termM[j] + cache.uthetap[j]
        Mmj = -cache.termM[j] + cache.uthetam[j]
        Mp2 = Mpj * Mpj
        Mm2 = Mmj * Mmj
        @inbounds for i in 1:N
            r = cache.r[i]
            r2 = r * r
            rhoi = cache.rho[i, j]
            rhobari = cache.rhobar[i, j]
            Yi = Y[i]
            Ypi = Yp[i]
            Xi = X[i]
            phase_r = radial_tortoise[i]
            L1pL2pS = cache.invrho[i, j] * L1L2S + theta_const + theta_L1 * L1j
            Wll = 2.0 * rhobari * L1pL2pS * r2 * Yi * phase_r
            termA = L2S + im * a * stj * cache.rho_minus[i, j] * S0j
            termB = 2.0 * Yi + r * Ypi
            termC = L2S * cache.rho_plus[i, j] * r * Yi
            Wlm = 4.0 * (r * rhobari * (termA * termB + termC) / rhoi) * phase_r
            termD = Xi / (2.0 * sqrt(r2 + a^2))
            termE = (Yi + 2.0 * r * Ypi) * phase_r
            termF = rhoi * r * (2.0 * Yi + r * Ypi) * phase_r
            Wmm = 4.0 * S0j * rhobari * (termD + termE + termF) / rhoi
            exp_rp = radial_phase[i]
            exp_rm = conj(radial_phase[i])
            Jpp = (Wll * cache.Lp2[i] + Wlm * cache.Lp[i] * Mpj + Wmm * Mp2) * exp_rp * exp_theta_p * inv4pi2
            Jpm = (Wll * cache.Lp2[i] + Wlm * cache.Lp[i] * Mmj + Wmm * Mm2) * exp_rp * exp_theta_m * inv4pi2
            Jmp = (Wll * cache.Lm2[i] + Wlm * cache.Lm[i] * Mpj + Wmm * Mp2) * exp_rm * exp_theta_p * inv4pi2
            Jmm = (Wll * cache.Lm2[i] + Wlm * cache.Lm[i] * Mmj + Wmm * Mm2) * exp_rm * exp_theta_m * inv4pi2
            total_j += cache.wr[i] * (Jpp + Jpm + Jmp + Jmm)
        end
        total += cache.wtheta[j] * total_j
    end
    return prefactor * total
end

function _generic_m2_integral_cached!(cache::GenericM2FluxCache, KG_master::Dict, Y_samp::Dict, SH_samp::Dict, omega::Float64, m::Int, n::Int, k::Int, N::Int, K::Int)
    KG_samp = _cached_generic_sample!(cache, KG_master, N, K)
    mcache = _cached_m2_mcache!(cache, KG_samp, m)
    return _generic_m2_integral_threaded(mcache, Y_samp, SH_samp, omega, n, k)
end

function _generic_p2_integral_cached!(cache::GenericM2FluxCache, KG_master::Dict, Y_samp::Dict, SH_samp::Dict, omega::Float64, lambda::Float64, m::Int, n::Int, k::Int, N::Int, K::Int)
    KG_samp = _cached_generic_sample!(cache, KG_master, N, K)
    mcache = _cached_m2_mcache!(cache, KG_samp, m)
    return _generic_p2_integral_threaded(mcache, Y_samp, SH_samp, omega, lambda, n, k)
end

function generic_mode_flux_from_master_cached!(cache::GenericM2FluxCache, KG_master::Dict, s::Int, l::Int, m::Int, n::Int, k::Int; N0::Int = 64, K0::Int = 16, Nmax::Int = 2^14, Kmax::Int = 2^12, sample_tol::Float64 = 1e-3, tol::Float64 = 1e-8, max_flux::Float64 = 1.0, mode_abs_floor::Float64 = 0.0, zero_low_flux::Bool = false, threaded_sampling::Bool = false)
    if !(s == -2 || s == 2)
        error("Spin weight s must be either 2 or -2.")
    end
    ispow2(N0) || throw(ArgumentError("N0 must be a power of 2"))
    ispow2(K0) || throw(ArgumentError("K0 must be a power of 2"))
    ispow2(Nmax) || throw(ArgumentError("Nmax must be a power of 2"))
    ispow2(Kmax) || throw(ArgumentError("Kmax must be a power of 2"))
    N0 <= Nmax || throw(ArgumentError("N0 must not exceed Nmax"))
    K0 <= Kmax || throw(ArgumentError("K0 must not exceed Kmax"))

    Frequencies = KG_master["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    a = KG_master["a"]
    ω = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
    if _skip_radiative_mode(s, a, m, ω)
        res = _zero_radiative_mode(ω, KG_master; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        res["N_sample_requested"] = N0
        res["K_sample_requested"] = K0
        res["N_sample"] = N0
        res["K_sample"] = K0
        return res
    end
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * ω; method = "jacobi")
    Ysol = _generic_isem_y_solution(s, l, m, a, ω)

    N = N0
    K = K0
    KG_samp = _cached_generic_sample!(cache, KG_master, N, K)
    Ysamp = if s == -2
        threaded_sampling ? GridSampling.y_sample_m2_isem_threaded(Ysol, KG_samp) : GridSampling.y_sample_m2_isem(Ysol, KG_samp)
    else
        threaded_sampling ? GridSampling.y_sample_p2_isem_threaded(Ysol, KG_samp) : GridSampling.y_sample_p2_isem(Ysol, KG_samp)
    end
    SHsamp = threaded_sampling ? GridSampling.swsh_sample_threaded(SH, KG_samp) : GridSampling.swsh_sample(SH, KG_samp)
    lambda = (Ysamp["params"]).lambda
    integral = s == -2 ? _generic_m2_integral_cached!(cache, KG_master, Ysamp, SHsamp, ω, m, n, k, N, K) : _generic_p2_integral_cached!(cache, KG_master, Ysamp, SHsamp, ω, lambda, m, n, k, N, K)
    hf = s == 2 ? horizon_factor(ω, a, m) : 0.0
    energy = s == -2 ? abs2(integral) / (4.0pi * ω^2) : hf * abs2(integral)
    flux_scale = max(max_flux, eps(Float64))
    effective_sample_tol = min(sample_tol, 3.0 * sqrt(tol))
    min_stop_N = min(max(4N0, 256), Nmax)
    min_stop_K = min(max(4K0, 64), Kmax)

    while N < Nmax || K < Kmax
        N2 = min(2N, Nmax)
        K2 = min(2K, Kmax)
        next_sample = _cached_generic_sample!(cache, KG_master, N2, K2)
        next_Ysamp = threaded_sampling ? GridSampling.refine_generic_y_sample_threaded(Ysol, Ysamp, next_sample) : GridSampling.refine_generic_y_sample(Ysol, Ysamp, next_sample)
        next_SHsamp = threaded_sampling ? GridSampling.refine_generic_swsh_sample_threaded(SH, SHsamp, next_sample) : GridSampling.refine_generic_swsh_sample(SH, SHsamp, next_sample)
        lambda = (next_Ysamp["params"]).lambda
        integral2 = s == -2 ? _generic_m2_integral_cached!(cache, KG_master, next_Ysamp, next_SHsamp, ω, m, n, k, N2, K2) : _generic_p2_integral_cached!(cache, KG_master, next_Ysamp, next_SHsamp, ω, lambda, m, n, k, N2, K2)
        energy2 = s == -2 ? abs2(integral2) / (4.0pi * ω^2) : hf * abs2(integral2)
        N = N2
        K = K2
        excess = abs(energy2) / flux_scale
        factor = excess <= 1.0 ? 1.0 : min(sqrt(excess), 50.0)
        low_flux_cutoff = _low_flux_cutoff(max_flux, tol)
        relE = _relative_energy_change(energy2, energy)
        KG_samp = next_sample
        Ysamp = next_Ysamp
        SHsamp = next_SHsamp
        integral = integral2
        energy = energy2
        significant_for_alias = abs(energy2) > low_flux_cutoff
        allow_stop = (N2 >= min_stop_N && K2 >= min_stop_K) || !significant_for_alias
        has_flux_reference = max_flux > 100 * eps(Float64)
        low_flux_done = mode_abs_floor > 0.0 && abs(energy2) < mode_abs_floor
        suspect_min_N = min_stop_N
        suspect_min_K = min_stop_K
        if has_flux_reference && !low_flux_done
            if excess > 1e8
                suspect_min_N = min(max(suspect_min_N, 4096), Nmax)
                suspect_min_K = min(max(suspect_min_K, 512), Kmax)
            elseif excess > 1e4
                suspect_min_N = min(max(suspect_min_N, 2048), Nmax)
                suspect_min_K = min(max(suspect_min_K, 512), Kmax)
            elseif excess > 100.0
                suspect_min_N = min(max(suspect_min_N, 1024), Nmax)
                suspect_min_K = min(max(suspect_min_K, 256), Kmax)
            elseif excess > 10.0
                suspect_min_N = min(max(suspect_min_N, 2 * min_stop_N), Nmax)
                suspect_min_K = min(max(suspect_min_K, 2 * min_stop_K), Kmax)
            end
        end
        suspect_needs_extra = has_flux_reference && !low_flux_done && (N2 < suspect_min_N || K2 < suspect_min_K)
        if relE == 0.0
            if allow_stop && !suspect_needs_extra
                break
            end
            continue
        end
        if allow_stop && !suspect_needs_extra && (factor * relE <= effective_sample_tol || abs(energy2) < low_flux_cutoff || low_flux_done)
                if zero_low_flux && mode_abs_floor > 0.0 && abs(energy) < mode_abs_floor
                    integral = 0.0 + 0.0im
                    energy = 0.0
                end
                break
        end
    end

    carter_factor = _generic_carter_factor(_cached_m2_mcache!(cache, KG_samp, m), ω)
    res = Dict(
        "Amplitude" => integral,
        "omega" => ω,
        "EnergyFlux" => energy,
        "AngularMomentumFlux" => s == -2 ? m * abs2(integral) / (4.0pi * ω^3) : hf * m * abs2(integral) / ω,
        "CarterConstantFlux" => s == -2 ? abs2(integral) * (carter_factor + k * ϒθ) / (2.0pi * ω^3) : 2 * hf * abs2(integral) * (carter_factor + k * ϒθ) / ω,
        "Trajectory" => KG_samp,
        "YSolution" => Ysamp,
        "SWSH" => SHsamp,
        "N_sample_requested" => N0,
        "K_sample_requested" => K0,
        "N_sample" => N,
        "K_sample" => K,
    )
    if zero_low_flux && mode_abs_floor > 0.0 && energy == 0.0
        res["LowFluxZeroed"] = true
    end
    return res
end

function generic_mode_flux_from_master(KG_master::Dict, s::Int, l::Int, m::Int, n::Int, k::Int; N0::Int = 64, K0::Int = 16, Nmax::Int = 2^14, Kmax::Int = 2^12, sample_tol::Float64 = 1e-3, tol::Float64 = 1e-8, max_flux::Float64 = 1.0, mode_abs_floor::Float64 = 0.0, zero_low_flux::Bool = false, threaded_sampling::Bool = false)
    ispow2(N0) || throw(ArgumentError("N0 must be a power of 2"))
    ispow2(K0) || throw(ArgumentError("K0 must be a power of 2"))
    ispow2(Nmax) || throw(ArgumentError("Nmax must be a power of 2"))
    ispow2(Kmax) || throw(ArgumentError("Kmax must be a power of 2"))
    N0 <= Nmax || throw(ArgumentError("N0 must not exceed Nmax"))
    K0 <= Kmax || throw(ArgumentError("K0 must not exceed Kmax"))

    Frequencies = KG_master["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    a = KG_master["a"]
    ω = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
    if _skip_radiative_mode(s, a, m, ω)
        res = _zero_radiative_mode(ω, KG_master; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        res["N_sample_requested"] = N0
        res["K_sample_requested"] = K0
        res["N_sample"] = N0
        res["K_sample"] = K0
        return res
    end
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * ω; method = "jacobi")
    Ysol = _generic_isem_y_solution(s, l, m, a, ω)

    N = N0
    K = K0
    KG_samp = GridSampling.subsample_generic_sample(KG_master, N, K)
    Ysamp = if s == 2
        threaded_sampling ? GridSampling.y_sample_p2_isem_threaded(Ysol, KG_samp) : GridSampling.y_sample_p2_isem(Ysol, KG_samp)
    else
        threaded_sampling ? GridSampling.y_sample_m2_isem_threaded(Ysol, KG_samp) : GridSampling.y_sample_m2_isem(Ysol, KG_samp)
    end
    SHsamp = threaded_sampling ? GridSampling.swsh_sample_threaded(SH, KG_samp) : GridSampling.swsh_sample(SH, KG_samp)
    res = _generic_flux_from_sample(KG_samp, Ysamp, SHsamp, s, l, m, n, k, a, ω, ϒθ)
    flux_scale = max(max_flux, eps(Float64))
    effective_sample_tol = min(sample_tol, 3.0 * sqrt(tol))
    min_stop_N = min(max(4N0, 256), Nmax)
    min_stop_K = min(max(4K0, 64), Kmax)
    while N < Nmax || K < Kmax
        N2 = min(2N, Nmax)
        K2 = min(2K, Kmax)
        next_sample = GridSampling.subsample_generic_sample(KG_master, N2, K2)
        next_Ysamp = threaded_sampling ? GridSampling.refine_generic_y_sample_threaded(Ysol, Ysamp, next_sample) : GridSampling.refine_generic_y_sample(Ysol, Ysamp, next_sample)
        next_SHsamp = threaded_sampling ? GridSampling.refine_generic_swsh_sample_threaded(SH, SHsamp, next_sample) : GridSampling.refine_generic_swsh_sample(SH, SHsamp, next_sample)
        res2 = _generic_flux_from_sample(next_sample, next_Ysamp, next_SHsamp, s, l, m, n, k, a, ω, ϒθ)
        N = N2
        K = K2
        factor = max(abs(res2["EnergyFlux"]) / flux_scale, 1.0)
        low_flux_cutoff = _low_flux_cutoff(max_flux, tol)
        relE = _relative_energy_change(res2["EnergyFlux"], res["EnergyFlux"])
        if relE == 0.0
            KG_samp = next_sample
            Ysamp = next_Ysamp
            SHsamp = next_SHsamp
            res = res2
            if N >= min_stop_N && K >= min_stop_K && abs(res2["EnergyFlux"]) < flux_scale && abs(res2["EnergyFlux"]) < low_flux_cutoff
                break
            end
            continue
        end
        KG_samp = next_sample
        Ysamp = next_Ysamp
        SHsamp = next_SHsamp
        res = res2
        low_flux_done = mode_abs_floor > 0.0 && abs(res2["EnergyFlux"]) < mode_abs_floor
        significant_for_alias = abs(res2["EnergyFlux"]) > low_flux_cutoff
        allow_stop = (N >= min_stop_N && K >= min_stop_K) || !significant_for_alias
        if allow_stop && (factor * relE <= effective_sample_tol || abs(res2["EnergyFlux"]) < low_flux_cutoff || low_flux_done) && abs(res2["EnergyFlux"]) < flux_scale
            if zero_low_flux && low_flux_done
                res["Amplitude"] = 0.0 + 0.0im
                res["EnergyFlux"] = 0.0
                res["AngularMomentumFlux"] = 0.0
                res["CarterConstantFlux"] = 0.0
                res["LowFluxZeroed"] = true
            end
            break
        end
    end

    res["N_sample_requested"] = N0
    res["K_sample_requested"] = K0
    res["N_sample"] = N
    res["K_sample"] = K
    return res
end

function generic_mode_flux(a, p, e, x, s::Int, l::Int, m::Int, n::Int, k::Int; N0::Int = 64, K0::Int = 16, Nmax::Int = 2^14, Kmax::Int = 2^12, sample_tol::Float64 = 1e-3, tol::Float64 = 1e-8, max_flux::Float64 = 1.0, mode_abs_floor::Float64 = 0.0, zero_low_flux::Bool = false, threaded_sampling::Bool = false)
    KG = Kerr_Geodesics(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end
    KG_master = GridSampling.kerr_geo_generic_sample_dense(KG, Nmax, Kmax)
    return generic_mode_flux_from_master(KG_master, s, l, m, n, k; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, sample_tol = sample_tol, tol = tol, max_flux = max_flux, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
end

function trapezoidal_1d_integral(f::Vector{<:Number})
    # Get grid size: M = number of sampling points (1D)
    M = length(f)
    M < 2 && return 0.0  # Need at least 2 points for integration
    
    # Step size: integration range [0, π] for 1D phase (qr or qθ)
    Δq = π / (M - 1)  # Uniform step between adjacent points
    
    # --------------------------
    # Assign weights based on position
    # --------------------------
    # Initialize weight vector with all 1.0 (internal points default)
    weight = ones(Float64, M)
    
    # Edge points (first and last) have weight 0.5 (trapezoidal rule)
    weight[1] = 0.5    # First point
    weight[end] = 0.5  # Last point
    
    # --------------------------
    # Compute 1D trapezoidal integral
    # --------------------------
    # Element-wise product of function values and weights
    weighted_f = f .* weight
    # Sum all weighted values and multiply by step size
    integral_value = sum(weighted_f) * Δq
    
    return integral_value
end

function trapezoidal_2d_integral(f::Matrix{<:Number})
    # Get grid size: N×K (N=radial points, K=polar points)
    N, K = size(f)
    
    # Calculate step sizes: integration range [0, π] for both qr and qθ
    Δqr = π / (N - 1)  # Step in radial direction (qr)
    Δqθ = π / (K - 1)  # Step in polar direction (qθ)
    
    total = sum(f)
    edge = sum(@view f[1, :]) + sum(@view f[N, :]) + sum(@view f[2:N-1, 1]) + sum(@view f[2:N-1, K])
    corner = f[1, 1] + f[1, K] + f[N, 1] + f[N, K]
    return Δqr * Δqθ * (total - 0.5 * edge - 0.25 * corner)
end

function cheb_nodes_and_D(n::Int)
    # Chebyshev extremal nodes on reference interval [-1, 1]
    x_ref = [cos(π * k / (n - 1)) for k in 0:n-1]
    
    # Construct Chebyshev differentiation matrix
    D_ref = zeros(Float64, n, n)
    for i in 1:n, j in 1:n
        if i != j
            # Weight factors (2 for endpoints, 1 otherwise)
            c_i = (i == 1 || i == n) ? 2.0 : 1.0
            c_j = (j == 1 || j == n) ? 2.0 : 1.0
            # Off-diagonal elements (standard Chebyshev differentiation formula)
            D_ref[i,j] = (-1)^(i+j) * c_i / (c_j * (x_ref[i] - x_ref[j]))
        end
    end
    # Diagonal elements enforce row sum = 0 (property of differentiation matrices)
    for i in 1:n
        D_ref[i,i] = -sum(D_ref[i,:])
    end
    return x_ref, D_ref
end

function levin_1d_integral(f_vals::Vector{<:Number}, gprime_vals::Vector{<:Number}, g0::Number, gπ::Number)
    k = length(f_vals)
    @assert length(gprime_vals) == k "gprime length must match f ($(length(gprime_vals)) ≠ $k)"
    
    # Fixed integration interval [0, π]
    a, b = 0.0, π
    scale = 2.0 / (b - a)  # Scaling factor for differentiation matrix
    
    # Generate Chebyshev nodes (reference [-1,1] and physical [0,π]) and differentiation matrix
    x_ref, D_ref = cheb_nodes_and_D(k)
    x_phys = (b - a) * (x_ref .+ 1.0) ./ 2.0 .+ a  # Map [-1,1] → [0,π]
    D = scale * D_ref  # Physical-space differentiation matrix
    
    # Prepare complex-valued inputs for numerical stability
    rhs = convert(Vector{ComplexF64}, f_vals)
    gprime_hat = convert(Vector{ComplexF64}, gprime_vals)
    
    # Construct system matrix: A = D + i·diag(g') (core of Levin method)
    A = D .+ 1im .* Diagonal(gprime_hat)
    
    # Solve linear system using pivoted QR decomposition (numerically robust)
    F = qr(A, Val(true))  # Pivoted QR factorization
    R = F.R               # Upper triangular factor
    piv = F.p             # Pivot indices
    diagR = abs.(diag(R)) # Magnitudes of diagonal elements (singularity measure)
    Anorm = maximum(diagR)# Estimate of matrix norm
    tol_sing = Anorm * 1e-12  # Singularity tolerance
    
    # Compute solution p(x) at Chebyshev nodes
    pvals = zeros(ComplexF64, k)
    if !all(diagR .< tol_sing)
        l = count(≥(tol_sing), diagR)  # Number of significant singular values
        y = F.Q[:, 1:l]' * rhs          # Project RHS onto significant subspace
        ytop = R[1:l, 1:l] \ y          # Solve reduced system
        colp = zeros(ComplexF64, k)
        colp[piv[1:l]] = ytop           # Map solution back to original indices
        pvals .= colp
    end
    
    # Identify boundary indices using physical coordinates (most reliable method)
    idx_x0 = findmin(x_phys)[2]  # Index of x=0 in physical nodes
    idx_xπ = findmax(x_phys)[2]  # Index of x=π in physical nodes
    
    # Core integral result from Levin's method: p(b)e^{ig(b)} - p(a)e^{ig(a)}
    term_b = pvals[idx_xπ] * exp(1im * gπ)
    term_a = pvals[idx_x0] * exp(1im * g0)
    return term_b - term_a
end

function levin_2d_integral(f_vals::Matrix{<:Number}, g1_prime_vals::Vector{<:Number}, g2_prime_vals::Vector{<:Number}, g1_0::Number, g1_π::Number, g2_0::Number, g2_π::Number)
    # Get grid dimensions (n×k, non-square allowed)
    n, k = size(f_vals)
    
    # Validate input dimensions
    @assert length(g1_prime_vals) == n "g1_prime_vals length must match x-dimension n ($(length(g1_prime_vals)) ≠ $n)"
    @assert length(g2_prime_vals) == k "g2_prime_vals length must match y-dimension k ($(length(g2_prime_vals)) ≠ $k)"
    
    # ----------------------
    # Step 1: Solve y-direction Levin equation for each x (n rows)
    # ----------------------
    # Generate y-direction Chebyshev nodes and differentiation matrix
    y_ref, D_ref_y = cheb_nodes_and_D(k)
    y_phys = π .* (y_ref .+ 1.0) ./ 2.0  # y ∈ [0,π]
    scale_y = 2.0 / π  # Scale factor for y-direction
    D_y = scale_y * D_ref_y  # y-differentiation matrix
    
    # Storage for y-solutions: p(x,y) has size n×k
    p_vals = zeros(ComplexF64, n, k)
    
    for i in 1:n  # Iterate over x-nodes (each row of f_vals)
        # Extract f(xᵢ, y) for fixed xᵢ (i-th row)
        f_row = f_vals[i, :]
        rhs_y = convert(Vector{ComplexF64}, f_row)
        
        # Construct y-matrix: Aᵧ = Dᵧ + i·diag(∂g₂/∂y)
        A_y = D_y .+ 1im .* Diagonal(convert(Vector{ComplexF64}, g2_prime_vals))
        
        # Solve Aᵧ·p = f using pivoted QR
        F_y = qr(A_y, Val(true))
        R_y = F_y.R
        piv_y = F_y.p
        diagR_y = abs.(diag(R_y))
        Anorm_y = maximum(diagR_y)
        tol_sing_y = Anorm_y * 1e-16
        
        p_row = zeros(ComplexF64, k)
        if !all(diagR_y .< tol_sing_y)
            l_y = count(≥(tol_sing_y), diagR_y)
            y_proj = F_y.Q[:, 1:l_y]' * rhs_y
            y_sol = R_y[1:l_y, 1:l_y] \ y_proj
            colp = zeros(ComplexF64, k)
            colp[piv_y[1:l_y]] = y_sol
            p_row .= colp
        end
        p_vals[i, :] = p_row
    end
    
    # ----------------------
    # Step 2: Reduce to x-direction problem using y-boundaries
    # ----------------------
    # Find y-boundary indices (y=0 and y=π)
    idx_y0 = findmin(y_phys)[2]  # Index of y=0 in y_phys
    idx_yπ = findmax(y_phys)[2]  # Index of y=π in y_phys
    
    # Compute h(x) = p(x, π)e^{ig₂(π)} - p(x, 0)e^{ig₂(0)} (length n)
    h_vals = [p_vals[i, idx_yπ] * exp(1im * g2_π) - p_vals[i, idx_y0] * exp(1im * g2_0) for i in 1:n]
    
    # ----------------------
    # Step 3: Solve x-direction Levin equation
    # ----------------------
    # Generate x-direction Chebyshev nodes and differentiation matrix
    x_ref, D_ref_x = cheb_nodes_and_D(n)
    x_phys = π .* (x_ref .+ 1.0) ./ 2.0  # x ∈ [0,π]
    scale_x = 2.0 / π  # Scale factor for x-direction
    D_x = scale_x * D_ref_x  # x-differentiation matrix
    
    # Construct x-matrix: Aₓ = Dₓ + i·diag(∂g₁/∂x)
    A_x = D_x .+ 1im .* Diagonal(convert(Vector{ComplexF64}, g1_prime_vals))
    rhs_x = convert(Vector{ComplexF64}, h_vals)
    
    # Solve Aₓ·q = h using pivoted QR
    F_x = qr(A_x, Val(true))
    R_x = F_x.R
    piv_x = F_x.p
    diagR_x = abs.(diag(R_x))
    Anorm_x = maximum(diagR_x)
    tol_sing_x = Anorm_x * 1e-16
    
    q_vals = zeros(ComplexF64, n)
    if !all(diagR_x .< tol_sing_x)
        l_x = count(≥(tol_sing_x), diagR_x)
        x_proj = F_x.Q[:, 1:l_x]' * rhs_x
        x_sol = R_x[1:l_x, 1:l_x] \ x_proj
        colp_x = zeros(ComplexF64, n)
        colp_x[piv_x[1:l_x]] = x_sol
        q_vals .= colp_x
    end
    
    # ----------------------
    # Step 4: Compute final 2D integral from x-boundaries
    # ----------------------
    idx_x0 = findmin(x_phys)[2]  # Index of x=0 in x_phys
    idx_xπ = findmax(x_phys)[2]  # Index of x=π in x_phys
    integral_val = q_vals[idx_xπ] * exp(1im * g1_π) - q_vals[idx_x0] * exp(1im * g1_0)
    
    return integral_val
end

function convolution_integral_generic_trapezoidal(a, p, e, x, s, l, m, n, k, N_sample, K_sample)

    KG = Kerr_Geodesics(a, p, e, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(500.0, 10pi / abs(omega))
    KG_samp = kerr_geo_generic_sample(KG, N_sample, K_sample)
    carter_samp = carter_ingredients_sample(KG_samp, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    if s == 2
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        Yup_samp = y_sample_p2(Yup_soln, KG_samp)
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        inte_in_samp_p2 = integrand_generic_sample_p2(KG_samp, Yup_samp, SH_p2_samp, n, k)
        integral_p2 = trapezoidal_2d_integral(inte_in_samp_p2)
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_p2,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_p2),
            "AngularMomentumFlux" => hf * m * abs2(integral_p2) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_p2) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        Yin_samp = y_sample_m2(Yin_soln, KG_samp)
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        inte_in_samp_m2 = integrand_generic_sample_m2(KG_samp, Yin_samp, SH_m2_samp, n, k)
        integral_m2 = trapezoidal_2d_integral(inte_in_samp_m2)
        return Dict(
            "Amplitude" => integral_m2,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_m2)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_m2) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_m2) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_generic_levin(a, p, e, x, s, l, m, n, k, N_sample, K_sample)

    KG = Kerr_Geodesics(a, p, e, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(500.0, 10pi / abs(omega))
    KG_samp = kerr_geo_generic_sample_cheby(KG, N_sample, K_sample)
    KG_trap = kerr_geo_generic_sample(KG, N_sample, K_sample)
    carter_samp = carter_ingredients_sample(KG_trap, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    if s == 2
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        Yup_samp = y_sample_p2(Yup_soln, KG_samp)
        Jpp_grid_up, Jpm_grid_up, Jmp_grid_up, Jmm_grid_up, drphase_up, dθphase_up, rphaseL_up, rphaseR_up, θphaseL_up, θphaseR_up, prefactor_up = integrand_generic_sample_cheby_p2(KG_samp, Yup_samp, SH_p2_samp, n, k)
        integralpp_up = levin_2d_integral(Jpp_grid_up, drphase_up, dθphase_up, rphaseL_up, rphaseR_up, θphaseL_up, θphaseR_up)
        integralpm_up = levin_2d_integral(Jpm_grid_up, drphase_up, - dθphase_up, rphaseL_up, rphaseR_up, - θphaseL_up, - θphaseR_up)
        integralmp_up = levin_2d_integral(Jmp_grid_up, - drphase_up, dθphase_up, - rphaseL_up, - rphaseR_up, θphaseL_up, θphaseR_up)
        integralmm_up = levin_2d_integral(Jmm_grid_up, - drphase_up, - dθphase_up, - rphaseL_up, - rphaseR_up, - θphaseL_up, - θphaseR_up)
        integral_up = (integralpp_up + integralpm_up + integralmp_up + integralmm_up) * prefactor_up
        
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_up,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_up),
            "AngularMomentumFlux" => hf * m * abs2(integral_up) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_up) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        Yin_samp = y_sample_m2(Yin_soln, KG_samp)
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        Jpp_grid_in, Jpm_grid_in, Jmp_grid_in, Jmm_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_generic_sample_cheby_m2(KG_samp, Yin_samp, SH_m2_samp, n, k)
        integralpp_in = levin_2d_integral(Jpp_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in)
        integralpm_in = levin_2d_integral(Jpm_grid_in, drphase_in, - dθphase_in, rphaseL_in, rphaseR_in, - θphaseL_in, - θphaseR_in)
        integralmp_in = levin_2d_integral(Jmp_grid_in, - drphase_in, dθphase_in, - rphaseL_in, - rphaseR_in, θphaseL_in, θphaseR_in)
        integralmm_in = levin_2d_integral(Jmm_grid_in, - drphase_in, - dθphase_in, - rphaseL_in, - rphaseR_in, - θphaseL_in, - θphaseR_in)
        integral_in = (integralpp_in + integralpm_in + integralmp_in + integralmm_in) * prefactor_in
        return Dict(
            "Amplitude" => integral_in,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_in)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_in) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_in) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_eccentric_trapezoidal(a, p, e, s, l, m, n, N_sample)
    KG = Kerr_Geodesics(a, p, e, 1.0)  # Assume this handles equatorial case
    
    # Extract frequencies and compute omega (radial-only mode)
    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr) / Γ  # No θ contribution (k omitted)
    
    # Radial solution parameters
    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(500.0, 10pi / abs(omega))
    KG_samp = kerr_geo_eccentric_sample(KG, N_sample)  # Equatorial sampling (θ fixed)
    
    if s == 2
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yup_samp = y_sample_p2(Yup_soln, KG_samp)
        S0_p2 = SH_p2(π/2, 0.0)      
        S1_p2 = SH_p2(π/2, 0.0; theta_derivative=1)
        S2_p2 = (m^2 + (a * omega)^2 - 2 - 2 * a * omega * m - SH_p2.lambda) * S0_p2
        inte_samp_p2 = integrand_eccentric_sample_p2(KG_samp, Yup_samp, (S0_p2, S1_p2, S2_p2), n)
        integral_p2 = trapezoidal_1d_integral(inte_samp_p2)
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_p2,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_p2),
            "AngularMomentumFlux" => hf * m * abs2(integral_p2) / omega,
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")  # θ-fixed harmonics
        Yin_samp = y_sample_m2(Yin_soln, KG_samp) 
        S0_m2 = SH_m2(π/2, 0.0)      
        S1_m2 = SH_m2(π/2, 0.0; theta_derivative=1)
        S2_m2 = (m^2 + (a * omega)^2 + 2 - 2 * a * omega * m - SH_m2.lambda) * S0_m2
        inte_samp_m2 = integrand_eccentric_sample_m2(KG_samp, Yin_samp, (S0_m2, S1_m2, S2_m2), n)
        integral_m2 = trapezoidal_1d_integral(inte_samp_m2)
        return Dict(
            "Amplitude" => integral_m2,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_m2)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_m2) / (4.0pi * omega^3),
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_eccentric_levin(a, p, e, s, l, m, n, N_sample)

    KG = Kerr_Geodesics(a, p, e, 1.0)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr) / Γ
    if _skip_radiative_mode(s, a, m, omega)
        result = _zero_radiative_mode(omega, KG; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        _eccentric_levin_last_key[] = key
        _eccentric_levin_last_result[] = result
        return result
    end
    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(500.0, 10pi / abs(omega))
    KG_samp = kerr_geo_eccentric_sample_cheby(KG, N_sample)

    if s == 2
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yup_samp = y_sample_p2(Yup_soln, KG_samp)
        S0_p2 = SH_p2(π/2, 0.0)      
        S1_p2 = SH_p2(π/2, 0.0; theta_derivative=1)
        S2_p2 = (m^2 + (a * omega)^2 - 2 - 2 * a * omega * m - SH_p2.lambda) * S0_p2
        Jp_up, Jm_up, drphase_up, rphaseL_up, rphaseR_up, prefactor_up = integrand_eccentric_sample_cheby_p2(KG_samp, Yup_samp, (S0_p2, S1_p2, S2_p2), n)
        integralp_p2 = levin_1d_integral(Jp_up, drphase_up, rphaseL_up, rphaseR_up)
        integralm_p2 = levin_1d_integral(Jm_up, - drphase_up, - rphaseL_up, - rphaseR_up)
        integral_p2 = (integralp_p2 + integralm_p2) * prefactor_up
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_p2,
            "omega" => omega,
            "EnergyFlux_hor" => hf * abs2(integral_p2),
            "AngularMomentumFlux_hor" => hf * m * abs2(integral_p2) / omega,
            "CarterConstantFlux_hor" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yin_samp = y_sample_m2(Yin_soln, KG_samp)   
        S0_m2 = SH_m2(π/2, 0.0)      
        S1_m2 = SH_m2(π/2, 0.0; theta_derivative=1)
        S2_m2 = (m^2 + (a * omega)^2 + 2 - 2 * a * omega * m - SH_m2.lambda) * S0_m2
        Jp_in, Jm_in, drphase_in, rphaseL_in, rphaseR_in, prefactor_in = integrand_eccentric_sample_cheby_m2(KG_samp, Yin_samp, (S0_m2, S1_m2, S2_m2), n)
        integralp_m2 = levin_1d_integral(Jp_in, drphase_in, rphaseL_in, rphaseR_in)
        integralm_m2 = levin_1d_integral(Jm_in, - drphase_in, - rphaseL_in, - rphaseR_in)
        integral_m2 = (integralp_m2 + integralm_m2) * prefactor_in
        return Dict(
            "Amplitude" => integral_m2,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_m2)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_m2) / (4.0pi * omega^3),
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_inclined_trapezoidal(a, p, x, s, l, m, k, K_sample)
    KG = Kerr_Geodesics(a, p, 0.0, x)  # Assume this handles radial-fixed case
    
    # Extract frequencies and compute omega (polar-only mode)
    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + k * ϒθ) / Γ  # No radial contribution (n omitted)
    
    # Radial solution parameters (r fixed, but retain for consistency)
    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(200.0, 10pi / abs(omega))
    KG_samp = kerr_geo_inclined_sample(KG, K_sample)
    carter_samp = carter_ingredients_sample(KG_samp, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    if s == 2
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yup = Dict("params" => (s=2, l=l, m=m, a=a, omega=omega, lambda=SH_p2.lambda),
            "Cinc" => Yup_soln.asymptotic.Cinc,
            "Y" => Yup_soln.solution.Y_inf(rstar_from_r(a, p)),
            "Yp" => Yup_soln.solution.Yp_inf(rstar_from_r(a, p)),
            "X" => Yup_soln.solution.X(rstar_from_r(a, p)))
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        inte_samp_up = integrand_inclined_sample_p2(KG_samp, Yup, SH_p2_samp, k)
        integral_up = trapezoidal_1d_integral(inte_samp_up)
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_up,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_up),
            "AngularMomentumFlux" => hf * m * abs2(integral_up) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_up) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yin = Dict("params" => (s=-2, l=l, m=m, a=a, omega=omega, lambda=SH_m2.lambda),
            "Binc" => Yin_soln.asymptotic.Binc,
            "Y" => Yin_soln.solution.Y_hor(rstar_from_r(a, p)),
            "Yp" => Yin_soln.solution.Yp_hor(rstar_from_r(a, p)),
            "X" => Yin_soln.solution.X(rstar_from_r(a, p)))  
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        inte_samp_in = integrand_inclined_sample_m2(KG_samp, Yin, SH_m2_samp, k)
        integral_in = trapezoidal_1d_integral(inte_samp_in)
        return Dict(
            "Amplitude" => integral_in,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_in)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_in) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_in) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_inclined_levin(a, p, x, s, l, m, k, K_sample)

    KG = Kerr_Geodesics(a, p, 0.0, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + k * ϒθ) / Γ
    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(200.0, 10pi / abs(omega))
    KG_samp = kerr_geo_inclined_sample_cheby(KG, K_sample)
    KG_trap = kerr_geo_inclined_sample(KG, K_sample)
    carter_samp = carter_ingredients_sample(KG_trap, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    if s == 2
        Yup_soln = Y_solution(s, l, m, a, omega, UP; rsin = rsin, rsout = rsout)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yup = Dict("params" => (s=2, l=l, m=m, a=a, omega=omega, lambda=SH_p2.lambda),
            "Cinc" => Yup_soln.asymptotic.Cinc,
            "Y" => Yup_soln.solution.Y_inf(rstar_from_r(a, p)),
            "Yp" => Yup_soln.solution.Yp_inf(rstar_from_r(a, p)),
            "X" => Yup_soln.solution.X(rstar_from_r(a, p)))
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        Jp_up, Jm_up, dθphase_up, θphaseL_up, θphaseR_up, prefactor_up = integrand_inclined_sample_cheby_p2(KG_samp, Yup, SH_p2_samp, k)
        integralp_up = levin_1d_integral(Jp_up, dθphase_up, θphaseL_up, θphaseR_up)
        integralm_up = levin_1d_integral(Jm_up, - dθphase_up, - θphaseL_up, - θphaseR_up)
        integral_up = (integralp_up + integralm_up) * prefactor_up
        hf = horizon_factor(omega, a, m)
        return Dict(
            "Amplitude" => integral_up,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_up),
            "AngularMomentumFlux" => hf * m * abs2(integral_up) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_up) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = Y_solution(s, l, m, a, omega, IN; rsin = rsin, rsout = rsout)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yin = Dict("params" => (s=-2, l=l, m=m, a=a, omega=omega, lambda=SH_m2.lambda),
            "Binc" => Yin_soln.asymptotic.Binc,
            "Y" => Yin_soln.solution.Y_hor(rstar_from_r(a, p)),
            "Yp" => Yin_soln.solution.Yp_hor(rstar_from_r(a, p)),
            "X" => Yin_soln.solution.X(rstar_from_r(a, p)))
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        Jp_in, Jm_in, dθphase_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_inclined_sample_cheby_m2(KG_samp, Yin, SH_m2_samp, k)
        integralp_in = levin_1d_integral(Jp_in, dθphase_in, θphaseL_in, θphaseR_in)
        integralm_in = levin_1d_integral(Jm_in, - dθphase_in, - θphaseL_in, - θphaseR_in)
        integral_in = (integralp_in + integralm_in) * prefactor_in
        return Dict(
            "Amplitude" => integral_in,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_in)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_in) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_in) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_circular_equatorial_m2(a, p, l, m)

    KG = Kerr_Geodesics(a, p, 0.0, 1.0)
    E = KG["Energy"]
    Lz = KG["AngularMomentum"]
    Γ = KG["Frequencies"]["ϒt"]
    ϒφ = KG["Frequencies"]["ϒϕ"]
    ω = m * ϒφ / Γ 

    r = p
    θ = π/2
    r2 = r^2
    Δ = r2 - 2r + a^2 
    rp = 1.0 + sqrt(1.0 - a^2)
    rm = 1.0 - sqrt(1.0 - a^2)

    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(500.0, 10π / abs(ω))
    Y_soln = Y_solution(-2, l, m, a, ω, IN; rsin=rsin, rsout=rsout)
    Y = Y_soln.solution.Y_hor(rstar_from_r(a, r))
    Yp = Y_soln.solution.Yp_hor(rstar_from_r(a, r))
    X = Y_soln.solution.X(rstar_from_r(a, r))
    Binc = Y_soln.asymptotic.Binc 

    c = a * ω
    SH = spin_weighted_spheroidal_harmonic(-2, l, m, c)
    S0 = SH(θ, 0.0)
    S1 = SH(θ, 0.0; theta_derivative=1)
    λ = SH.lambda
    S2 = (m^2 + c^2 + 2 - 2c*m - λ) * S0

    ρ = -1.0 / r
    ρbar = -1.0 / r

    sinθ = sin(θ)
    L1 = -m + c
    L2 = -m + c
    L2p = -2.0

    L2S = S1 + L2 * S0
    L1Sp = S2 + L1 * S1
    L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0
    L1pL2pS = -r * L1L2S + 3im * a * L1 * S0 + 2im * a * S1 -1im * a * L2 * S0

    rs = rstar_from_r(a, r)
    log_term = log((r - rp)/(r - rm))
    phase = exp(im*ω*rs - im*a*m*log_term/(2*sqrt(1 - a^2)))

    Wnn = ρbar * L1pL2pS * r2 * Y * phase / 2.0

    term_A = L2S + im*a*sinθ*(ρ - ρbar)*S0 
    term_B = 2Y + r*Yp
    term_C = L2S*(ρ + ρbar)*r*Y
    Wnmbar = (-r * ρbar * (term_A * term_B + term_C) / ρ) * phase
    
    term_D = X / (2*sqrt(r2 + a^2))
    term_E = (Y + 2r*Yp) * phase
    term_F = ρ*r*(2Y + r*Yp) * phase
    Wmbarmbar = S0 * ρbar * (term_D + term_E + term_F) / ρ

    numerator_N = E * (r2 + a^2) - a * Lz 
    N = numerator_N / Δ 
    term_M = im * sinθ * (a*E - Lz / (sinθ^2))
    Mbar = term_M 

    J = Wnn*N^2 + Wnmbar*N*Mbar + Wmbarmbar*Mbar^2
    integral = 4im * π * ω * J / (Binc * Γ)

    return Dict(
        "Amplitude" => integral,
        "omega" => ω,
        "EnergyFlux" => abs2(integral)/(4.0pi*ω^2),
        "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * ω^3),
        "CarterConstantFlux" => 0.0,
        "Trajectory" => KG,
        "YSolution" => Y_soln,
        "SWSH" => SH
    )
end

function convolution_integral_circular_equatorial_p2(a, p, l, m)

    KG = Kerr_Geodesics(a, p, 0.0, 1.0)
    E = KG["Energy"]
    Lz = KG["AngularMomentum"]
    Γ = KG["Frequencies"]["ϒt"]
    ϒφ = KG["Frequencies"]["ϒϕ"]
    ω = m * ϒφ / Γ 

    r = p
    θ = π/2
    r2 = r^2
    Δ = r2 - 2*r + a^2 
    rp = 1.0 + sqrt(1.0 - a^2)
    rm = 1.0 - sqrt(1.0 - a^2)
    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(500.0, 10π / abs(ω))

    Y_soln = Y_solution(2, l, m, a, ω, UP; rsin=rsin, rsout=rsout)
    Y = Y_soln.solution.Y_inf(rstar_from_r(a, r))
    Yp = Y_soln.solution.Yp_inf(rstar_from_r(a, r))
    X = Y_soln.solution.X(rstar_from_r(a, r))
    Cinc = Y_soln.asymptotic.Cinc

    c = a * ω
    SH = spin_weighted_spheroidal_harmonic(2, l, m, c)
    S0 = SH(θ, 0.0)
    S1 = SH(θ, 0.0; theta_derivative=1)
    λ = SH.lambda
    S2 = (m^2 + c^2 - 2 - 2c*m - λ) * S0

    ρ = -1.0 / r
    ρbar = -1.0 / r

    sinθ = sin(θ)
    L1 = m - c
    L2 = m - c
    L2p = -2.0

    L2S = S1 + L2 * S0
    L1Sp = S2 + L1 * S1
    L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0
    L1pL2pS = -r * L1L2S + 3im * a * L1 * S0 + 2im * a * S1 -1im * a * L2 * S0

    rs = rstar_from_r(a, r)
    log_term = log((r - rp)/(r - rm))
    phase = exp(-(im*ω*rs - im*a*m*log_term/(2*sqrt(1 - a^2))))

    Wll = 2.0 * ρbar * L1pL2pS * r2 * Y * phase
    term_A = L2S + im*a*sinθ*(ρ - ρbar)*S0 
    term_B = 2.0 * Y + r*Yp
    term_C = L2S*(ρ + ρbar)*r*Y
    Wlm = 4.0 * (r * ρbar * (term_A * term_B + term_C) / ρ) * phase
    Wmm = 2.0 * S0 * ρbar * (X / sqrt(r^2 + a^2) + 2 * phase * (Y + 2 * r * Y * ρ + r * Yp * (2 + r * ρ))) / ρ

    L = (E * (r2 + a^2) - a * Lz) / Δ 
    M = - im * sinθ * (a*E - Lz / (sinθ^2))

    J = Wll*L^2 + Wlm*L*M + Wmm*M^2

    κ = ω - a * m / (2 * rp)
    c0 = 24 + 12im * ω + λ * (10 + λ) - 12 * a * ω * (a * ω - m)
    c1 = -32im * a * m - 8im * a * m * λ + 8im * a^2 * ω * (1 + λ)
    c2 = 12 * a^2 - 24im * a * m - 24 * a^2 * m^2 + 24im * a^2 * ω + 48 * a^3 * m * ω - 24 * a^4 * ω^2
    c3 = - 24im * a^3 * (a * ω - m) - 24 * a^2
    c4 = 12 * a^4
    η = c0 + c1 / rp + c2 / rp^2 + c3 / rp^3 + c4 / rp^4
    factor = 2*sqrt(2)*rp^(3/2)*(((4*ω)*(im-4*ω)-a*m*(im-8*ω)-a^2*(m^2+2im*ω-4*ω^2))*rp^2
            +a^2*(im-4*ω)*(a*m-2*ω)*rp)/(2*rp^3*(24+10*λ+λ^2+12im*ω)-rp^2*(8im*a*m*(11
            +2*λ+6im*ω)+a^2*(24+24*m^2+10*λ+λ^2-28im*ω-16im*λ*ω+48*ω^2))+8im*a^3*rp*(m
            *(7+λ-6im*ω)-a*ω*(4+λ))+12*a^5*ω*(a*ω-3*m))
    ϵ0 = sqrt(1 - a^2) / (4 * rp)
    energy_factor = ω / (κ * (2 * rp)^3 * (κ^2 + 4 * ϵ0^2) * 64pi)

    integral = - im * π * factor * J * η / (κ * Cinc * Γ)

    return Dict(
        "Amplitude" => integral,
        "omega" => ω,
        "EnergyFlux" => abs2(integral) * energy_factor,
        "AngularMomentumFlux" => abs2(integral) * energy_factor * m / ω,
        "CarterConstantFlux" => 0.0,
        "Trajectory" => KG,
        "YSolution" => Y_soln,
        "SWSH" => SH
    )
end

function convolution_integral_circular_equatorial(a, p, s, l, m)
    if s == 2
        return convolution_integral_circular_equatorial_p2(a, p, l, m)
    elseif s == -2
        return convolution_integral_circular_equatorial_m2(a, p, l, m)
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_trapezoidal(a, p, e, x, s, l, m, n, k; N = 256, K = 64)
    KG = Kerr_Geodesics(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end
    if m == 0 && n == 0 && k == 0
        return Dict("Amplitude" => 0.0 + 0.0im, "omega" => 0.0, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
    end
    if isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0 || k != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_circular_equatorial(a, p, s, l, m)
    elseif isapprox(e, 0.0; atol=1e-12) && !isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_inclined_trapezoidal(a, p, x, s, l, m, k, K)
    elseif !isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if k != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_eccentric_trapezoidal(a, p, e, s, l, m, n, N)
    else
        return convolution_integral_generic_trapezoidal(a, p, e, x, s, l, m, n, k, N, K)
    end
end

function convolution_integral_levin(a, p, e, x, s, l, m, n, k; N = 256, K = 32)
    KG = Kerr_Geodesics(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end
    if m == 0 && n == 0 && k == 0
        return Dict("Amplitude" => 0.0 + 0.0im, "omega" => 0.0, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
    end
    if isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0 || k != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_circular_equatorial(a, p, s, l, m)
    elseif isapprox(e, 0.0; atol=1e-12) && !isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_inclined_levin(a, p, x, s, l, m, k, K)
    elseif !isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if k != 0
            return Dict(
                    "Amplitude" => 0.0 + 0.0im,
                    "omega" => nothing,
                    "EnergyFlux" => 0.0,
                    "AngularMomentumFlux" => 0.0,
                    "CarterConstantFlux" => 0.0,
                    "Trajectory" => KG,
                    "YSolution" => nothing,
                    "SWSH" => nothing
                    )
        end
        return convolution_integral_eccentric_levin(a, p, e, s, l, m, n, N)
    else
        return convolution_integral_generic_levin(a, p, e, x, s, l, m, n, k, N, K)
    end
end

function convolution_integral_circular_equatorial_isem(a, p, s, l, m)
    KG = Kerr_Geodesics(a, p, 0.0, 1.0)
    E = KG["Energy"]
    Lz = KG["AngularMomentum"]
    Γ = KG["Frequencies"]["ϒt"]
    ϒφ = KG["Frequencies"]["ϒϕ"]
    ω = m * ϒφ / Γ
    if _skip_radiative_mode(s, a, m, ω)
        return _zero_radiative_mode(ω, KG; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
    end

    r = p
    θ = π / 2
    r2 = r^2
    Δ = r2 - 2r + a^2
    rp = 1.0 + sqrt(1.0 - a^2)
    rm = 1.0 - sqrt(1.0 - a^2)

    Y_soln = _isem_y_solution(s, l, m, a, ω)
    Y, Yp, X, _ = Y_soln.Y_solution(r)
    A = Y_soln.mode.lambda

    c = a * ω
    SH = spin_weighted_spheroidal_harmonic(s, l, m, c)
    S0 = SH(θ, 0.0)
    S1 = SH(θ, 0.0; theta_derivative = 1)
    S2 = s == -2 ? (m^2 + c^2 + 2 - 2c*m - A) * S0 : (m^2 + c^2 - 2 - 2c*m - A) * S0

    ρ = -1.0 / r
    ρbar = -1.0 / r
    sinθ = sin(θ)

    if s == -2
        L1 = -m + c
        L2 = -m + c
        L2p = -2.0
        L2S = S1 + L2 * S0
        L1Sp = S2 + L1 * S1
        L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0
        L1pL2pS = -r * L1L2S + 3im * a * L1 * S0 + 2im * a * S1 - 1im * a * L2 * S0
        rs = rstar_from_r(a, r)
        log_term = log((r - rp)/(r - rm))
        phase = exp(im*ω*rs - im*a*m*log_term/(2*sqrt(1 - a^2)))

        Wnn = ρbar * L1pL2pS * r2 * Y * phase / 2.0
        term_A = L2S + im*a*sinθ*(ρ - ρbar)*S0
        term_B = 2Y + r*Yp
        term_C = L2S*(ρ + ρbar)*r*Y
        Wnmbar = (-r * ρbar * (term_A * term_B + term_C) / ρ) * phase
        term_D = X / (2*sqrt(r2 + a^2))
        term_E = (Y + 2r*Yp) * phase
        term_F = ρ*r*(2Y + r*Yp) * phase
        Wmbarmbar = S0 * ρbar * (term_D + term_E + term_F) / ρ
        N = (E * (r2 + a^2) - a * Lz) / Δ
        Mbar = im * sinθ * (a*E - Lz / (sinθ^2))
        J = Wnn*N^2 + Wnmbar*N*Mbar + Wmbarmbar*Mbar^2
        integral = 4im * π * ω * J / (_isem_gsn_incidence_amplitude(Y_soln) * Γ)
        return Dict(
            "Amplitude" => integral,
            "omega" => ω,
            "EnergyFlux" => abs2(integral)/(4.0pi*ω^2),
            "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * ω^3),
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Y_soln,
            "SWSH" => SH
        )
    elseif s == 2
        L1 = m - c
        L2 = m - c
        L2p = -2.0
        L2S = S1 + L2 * S0
        L1Sp = S2 + L1 * S1
        L1L2S = L1Sp + L2p * S0 + L2 * S1 + L1 * L2 * S0
        L1pL2pS = -r * L1L2S + 3im * a * L1 * S0 + 2im * a * S1 - 1im * a * L2 * S0
        rs = rstar_from_r(a, r)
        log_term = log((r - rp)/(r - rm))
        phase = exp(-(im*ω*rs - im*a*m*log_term/(2*sqrt(1 - a^2))))

        Wll = 2.0 * ρbar * L1pL2pS * r2 * Y * phase
        term_A = L2S + im*a*sinθ*(ρ - ρbar)*S0
        term_B = 2.0 * Y + r*Yp
        term_C = L2S*(ρ + ρbar)*r*Y
        Wlm = 4.0 * (r * ρbar * (term_A * term_B + term_C) / ρ) * phase
        Wmm = 2.0 * S0 * ρbar * (X / sqrt(r^2 + a^2) + 2 * phase * (Y + 2 * r * Y * ρ + r * Yp * (2 + r * ρ))) / ρ
        L = (E * (r2 + a^2) - a * Lz) / Δ
        M = - im * sinθ * (a*E - Lz / (sinθ^2))
        J = Wll*L^2 + Wlm*L*M + Wmm*M^2
        κ = ω - a * m / (2 * rp)
        c0 = 24 + 12im * ω + A * (10 + A) - 12 * a * ω * (a * ω - m)
        c1 = -32im * a * m - 8im * a * m * A + 8im * a^2 * ω * (1 + A)
        c2 = 12 * a^2 - 24im * a * m - 24 * a^2 * m^2 + 24im * a^2 * ω + 48 * a^3 * m * ω - 24 * a^4 * ω^2
        c3 = -24im * a^3 * (a * ω - m) - 24 * a^2
        c4 = 12 * a^4
        η = c0 + c1 / rp + c2 / rp^2 + c3 / rp^3 + c4 / rp^4
        factor = 2*sqrt(2)*rp^(3/2)*(((4*ω)*(im-4*ω)-a*m*(im-8*ω)-a^2*(m^2+2im*ω-4*ω^2))*rp^2 +
            a^2*(im-4*ω)*(a*m-2*ω)*rp)/(2*rp^3*(24+10*A+A^2+12im*ω)-rp^2*(8im*a*m*(11 +
            2*A+6im*ω)+a^2*(24+24*m^2+10*A+A^2-28im*ω-16im*A*ω+48*ω^2))+8im*a^3*rp*(m *
            (7+A-6im*ω)-a*ω*(4+A))+12*a^5*ω*(a*ω-3*m))
        ϵ0 = sqrt(1 - a^2) / (4 * rp)
        energy_factor = ω / (κ * (2 * rp)^3 * (κ^2 + 4 * ϵ0^2) * 64pi)
        integral = - im * π * factor * J * η / (κ * _isem_gsn_incidence_amplitude(Y_soln) * Γ)
        return Dict(
            "Amplitude" => integral,
            "omega" => ω,
            "EnergyFlux" => abs2(integral) * energy_factor,
            "AngularMomentumFlux" => abs2(integral) * energy_factor * m / ω,
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Y_soln,
            "SWSH" => SH
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_eccentric_trapezoidal_isem(a, p, e, s, l, m, n, N_sample; Nmax::Int = 2^14, kwargs...)
    KG = Kerr_Geodesics(a, p, e, 1.0)
    KG_sample = kerr_geo_eccentric_sample_dense(KG, Nmax)
    return convolution_integral_eccentric_trapezoidal_isem(KG_sample, s, l, m, n, N_sample; Nmax = Nmax, kwargs...)
end

@inline function _eccentric_flux_relerr(flux1, flux2, flux_max, tol)
    if flux1["EnergyFlux"] == flux2["EnergyFlux"]
        return 1.0
    end
    factor = max(abs(flux2["EnergyFlux"]) / flux_max / tol, 10.0)
    return factor * abs(flux1["EnergyFlux"] - flux2["EnergyFlux"]) / min(abs(flux1["EnergyFlux"]), abs(flux2["EnergyFlux"]))
end

function _eccentric_flux_from_sample(KG_samp::Dict, Ysamp::Dict, SH, s, a, ω, m, n)
    if s == 2
        S0 = SH(π/2, 0.0)
        S1 = SH(π/2, 0.0; theta_derivative = 1)
        S2 = (m^2 + (a * ω)^2 - 2 - 2 * a * ω * m - SH.lambda) * S0
        inte_samp = integrand_eccentric_sample_p2(KG_samp, Ysamp, (S0, S1, S2), n)
        integral = trapezoidal_1d_integral(inte_samp)
        hf = horizon_factor(ω, a, m)
        return Dict(
            "Amplitude" => integral,
            "omega" => ω,
            "EnergyFlux" => hf * abs2(integral),
            "AngularMomentumFlux" => hf * m * abs2(integral) / ω,
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG_samp,
            "YSolution" => Ysamp,
            "SWSH" => SH
        )
    elseif s == -2
        S0 = SH(π/2, 0.0)
        S1 = SH(π/2, 0.0; theta_derivative = 1)
        S2 = (m^2 + (a * ω)^2 + 2 - 2 * a * ω * m - SH.lambda) * S0
        inte_samp = integrand_eccentric_sample_m2(KG_samp, Ysamp, (S0, S1, S2), n)
        integral = trapezoidal_1d_integral(inte_samp)
        return Dict(
            "Amplitude" => integral,
            "omega" => ω,
            "EnergyFlux" => abs2(integral) / (4.0pi * ω^2),
            "AngularMomentumFlux" => m * abs2(integral) / (4.0pi * ω^3),
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG_samp,
            "YSolution" => Ysamp,
            "SWSH" => SH
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_eccentric_trapezoidal_isem(KG_sample::Dict, s, l, m, n, N_sample::Int64; Nmax::Int = 2^14, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 0.0, zero_low_flux::Bool = false, threaded_sampling::Bool = false, cache = nothing)
    ispow2(N_sample) || throw(ArgumentError("N_sample must be a power of 2"))
    ispow2(Nmax) || throw(ArgumentError("Nmax must be a power of 2"))
    N_sample <= Nmax || throw(ArgumentError("N_sample must not exceed Nmax"))
    Frequencies = KG_sample["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒφ = Frequencies["ϒϕ"]
    a = KG_sample["a"]
    ω = (m * ϒφ + n * ϒr) / Γ
    if _skip_radiative_mode(s, a, m, ω)
        res = _zero_radiative_mode(ω, KG_sample; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        res["N_sample_requested"] = N_sample
        res["N_sample"] = N_sample
        return res
    end
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * ω; method = "jacobi")
    Ysol = _isem_y_solution(s, l, m, a, ω)

    N = N_sample
    KG_samp = cache isa EccentricFluxCache ? _cached_eccentric_sample!(cache, KG_sample, N) : subsample_eccentric_sample(KG_sample, N)
    Ysamp = if s == 2
        threaded_sampling ? GridSampling.y_sample_p2_isem_threaded(Ysol, KG_samp) : y_sample_p2_isem(Ysol, KG_samp)
    else
        threaded_sampling ? GridSampling.y_sample_m2_isem_threaded(Ysol, KG_samp) : y_sample_m2_isem(Ysol, KG_samp)
    end
    res = _eccentric_flux_from_sample(KG_samp, Ysamp, SH, s, a, ω, m, n)
    flux_scale = max(Float64(max_flux), eps(Float64))
    effective_sample_tol = min(sample_tol, 3.0 * sqrt(Float64(tol)))
    min_stop_N = min(max(4N_sample, 256), Nmax)

    while N < Nmax
        N2 = min(2N, Nmax)
        KG_samp2 = cache isa EccentricFluxCache ? _cached_eccentric_sample!(cache, KG_sample, N2) : subsample_eccentric_sample(KG_sample, N2)
        Ysamp2 = threaded_sampling ? GridSampling.refine_generic_y_sample_threaded(Ysol, Ysamp, KG_samp2) : refine_eccentric_y_sample(Ysol, Ysamp, KG_sample, N)
        res2 = _eccentric_flux_from_sample(KG_samp2, Ysamp2, SH, s, a, ω, m, n)
        low_flux_cutoff = _low_flux_cutoff(Float64(max_flux), Float64(tol))
        relE = _relative_energy_change(res2["EnergyFlux"], res["EnergyFlux"])
        significant_for_alias = abs(res2["EnergyFlux"]) > low_flux_cutoff
        allow_stop = (N2 >= min_stop_N) || !significant_for_alias
        low_flux_done = mode_abs_floor > 0.0 && abs(res2["EnergyFlux"]) < mode_abs_floor
        if relE == 0.0
            N = N2
            KG_samp = KG_samp2
            Ysamp = Ysamp2
            res = res2
            if allow_stop && abs(res2["EnergyFlux"]) < flux_scale && abs(res2["EnergyFlux"]) < low_flux_cutoff
                break
            end
            continue
        end
        factor = max(abs(res2["EnergyFlux"]) / flux_scale, 1.0)
        if allow_stop && (factor * relE <= effective_sample_tol || abs(res2["EnergyFlux"]) < low_flux_cutoff || low_flux_done) && abs(res2["EnergyFlux"]) < flux_scale
            if zero_low_flux && low_flux_done
                res2["Amplitude"] = 0.0 + 0.0im
                res2["EnergyFlux"] = 0.0
                res2["AngularMomentumFlux"] = 0.0
                res2["CarterConstantFlux"] = 0.0
                res2["LowFluxZeroed"] = true
            end
            N = N2
            res = res2
            break
        end
        N = N2
        KG_samp = KG_samp2
        Ysamp = Ysamp2
        res = res2
    end

    res["N_sample_requested"] = N_sample
    res["N_sample"] = N
    return res
end

function convolution_integral_eccentric_levin_isem(a, p, e, s, l, m, n, N_sample)
    key = (a, p, e, s, l, m, n, N_sample)
    if _eccentric_levin_last_key[] === key
        return _eccentric_levin_last_result[]
    end

    KG = Kerr_Geodesics(a, p, e, 1.0)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + n * ϒr) / Γ
    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(500.0, 10pi / abs(omega))
    KG_samp = kerr_geo_eccentric_sample_cheby(KG, N_sample)

    result = if s == 2
        Yup_soln = _isem_y_solution(s, l, m, a, omega)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yup_samp = y_sample_p2_isem(Yup_soln, KG_samp)
        S0_p2 = SH_p2(π/2, 0.0)
        S1_p2 = SH_p2(π/2, 0.0; theta_derivative=1)
        S2_p2 = (m^2 + (a * omega)^2 - 2 - 2 * a * omega * m - SH_p2.lambda) * S0_p2
        Jp_up, Jm_up, drphase_up, rphaseL_up, rphaseR_up, prefactor_up = integrand_eccentric_sample_cheby_p2(KG_samp, Yup_samp, (S0_p2, S1_p2, S2_p2), n)
        integralp_p2 = levin_1d_integral(Jp_up, drphase_up, rphaseL_up, rphaseR_up)
        integralm_p2 = levin_1d_integral(Jm_up, - drphase_up, - rphaseL_up, - rphaseR_up)
        integral_p2 = (integralp_p2 + integralm_p2) * prefactor_up
        hf = horizon_factor(omega, a, m)
        Dict(
            "Amplitude" => integral_p2,
            "omega" => omega,
            "EnergyFlux_hor" => hf * abs2(integral_p2),
            "AngularMomentumFlux_hor" => hf * m * abs2(integral_p2) / omega,
            "CarterConstantFlux_hor" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = _isem_y_solution(s, l, m, a, omega)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yin_samp = y_sample_m2_isem(Yin_soln, KG_samp)
        S0_m2 = SH_m2(π/2, 0.0)
        S1_m2 = SH_m2(π/2, 0.0; theta_derivative=1)
        S2_m2 = (m^2 + (a * omega)^2 + 2 - 2 * a * omega * m - SH_m2.lambda) * S0_m2
        Jp_in, Jm_in, drphase_in, rphaseL_in, rphaseR_in, prefactor_in = integrand_eccentric_sample_cheby_m2(KG_samp, Yin_samp, (S0_m2, S1_m2, S2_m2), n)
        integralp_m2 = levin_1d_integral(Jp_in, drphase_in, rphaseL_in, rphaseR_in)
        integralm_m2 = levin_1d_integral(Jm_in, - drphase_in, - rphaseL_in, - rphaseR_in)
        integral_m2 = (integralp_m2 + integralm_m2) * prefactor_in
        Dict(
            "Amplitude" => integral_m2,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_m2)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_m2) / (4.0pi * omega^3),
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
    _eccentric_levin_last_key[] = key
    _eccentric_levin_last_result[] = result
    return result
end

function convolution_integral_inclined_trapezoidal_isem(KG_sample::Dict, s, l, m, k, K_sample::Int64; Kmax::Int = 2^12, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 0.0, zero_low_flux::Bool = false, threaded_sampling::Bool = false, cache = nothing)
    ispow2(K_sample) || throw(ArgumentError("K_sample must be a power of 2"))
    ispow2(Kmax) || throw(ArgumentError("Kmax must be a power of 2"))
    K_sample <= Kmax || throw(ArgumentError("K_sample must not exceed Kmax"))
    K_interval = K_sample
    Frequencies = KG_sample["Frequencies"]
    Γ = Frequencies["ϒt"]
    a = KG_sample["a"]
    ω = (m * Frequencies["ϒϕ"] + k * Frequencies["ϒθ"]) / Γ
    if _skip_radiative_mode(s, a, m, ω)
        res = _zero_radiative_mode(ω, KG_sample; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        res["K_sample_requested"] = K_sample
        res["K_sample"] = K_sample
        return res
    end
    ctx = _inclined_trapezoidal_context(KG_sample, s, l, m, k, K_interval, cache)
    SHsamp = threaded_sampling ? GridSampling.swsh_sample_threaded(ctx.SH, ctx.KG_samp) : ctx.SHsamp
    res = _inclined_flux_from_sample(ctx.KG_samp, ctx.Ysol, ctx.Ydic, SHsamp, s, ctx.a, ctx.omega, m, k)
    flux_scale = max(Float64(max_flux), eps(Float64))
    effective_sample_tol = min(sample_tol, 3.0 * sqrt(Float64(tol)))
    min_stop_K = min(max(4K_sample, 64), Kmax)

    while K_interval < Kmax
        K2 = min(2 * K_interval, Kmax)
        next_sample = cache isa InclinedFluxCache ? _cached_inclined_sample!(cache, KG_sample, K2) : subsample_inclined_sample(KG_sample, K2)
        next_SHsamp = threaded_sampling ? GridSampling.refine_generic_swsh_sample_threaded(ctx.SH, SHsamp, next_sample) : refine_inclined_swsh_sample(ctx.SH, SHsamp, next_sample)
        res2 = _inclined_flux_from_sample(next_sample, ctx.Ysol, ctx.Ydic, next_SHsamp, s, ctx.a, ctx.omega, m, k)
        low_flux_cutoff = _low_flux_cutoff(Float64(max_flux), Float64(tol))
        relE = _relative_energy_change(res2["EnergyFlux"], res["EnergyFlux"])
        relL = _relative_energy_change(res2["AngularMomentumFlux"], res["AngularMomentumFlux"])
        relQ = _relative_energy_change(res2["CarterConstantFlux"], res["CarterConstantFlux"])
        rel = max(relE, relL, relQ)
        significant_for_alias = abs(res2["EnergyFlux"]) > low_flux_cutoff
        allow_stop = (K2 >= min_stop_K) || !significant_for_alias
        low_flux_done = mode_abs_floor > 0.0 && abs(res2["EnergyFlux"]) < mode_abs_floor
        if rel == 0.0
            K_interval = K2
            SHsamp = next_SHsamp
            ctx = (KG_samp = next_sample, Ysol = ctx.Ysol, Ydic = ctx.Ydic, SH = ctx.SH, SHsamp = SHsamp, omega = ctx.omega, a = ctx.a, m = ctx.m, Trajectory = ctx.Trajectory)
            res = res2
            if allow_stop && abs(res2["EnergyFlux"]) < flux_scale && abs(res2["EnergyFlux"]) < low_flux_cutoff
                break
            end
            continue
        end
        factor = max(abs(res2["EnergyFlux"]) / flux_scale, 1.0)
        if allow_stop && (factor * rel <= effective_sample_tol || abs(res2["EnergyFlux"]) < low_flux_cutoff || low_flux_done) && abs(res2["EnergyFlux"]) < flux_scale
            if zero_low_flux && low_flux_done
                res2["Amplitude"] = 0.0 + 0.0im
                res2["EnergyFlux"] = 0.0
                res2["AngularMomentumFlux"] = 0.0
                res2["CarterConstantFlux"] = 0.0
                res2["LowFluxZeroed"] = true
            end
            K_interval = K2
            res = res2
            break
        end
        K_interval = K2
        SHsamp = next_SHsamp
        ctx = (KG_samp = next_sample, Ysol = ctx.Ysol, Ydic = ctx.Ydic, SH = ctx.SH, SHsamp = SHsamp, omega = ctx.omega, a = ctx.a, m = ctx.m, Trajectory = ctx.Trajectory)
        res = res2
    end

    res["K_sample_requested"] = K_sample
    res["K_sample"] = K_interval
    return res
end

function convolution_integral_inclined_trapezoidal_isem(a, p, x, s, l, m, k, K_sample; Kmax::Int = 2^12, kwargs...)
    KG = Kerr_Geodesics(a, p, 0.0, x)
    KG_sample = kerr_geo_inclined_sample_dense(KG, x, Kmax)
    return convolution_integral_inclined_trapezoidal_isem(KG_sample, s, l, m, k, K_sample; Kmax = Kmax, kwargs...)
end

function convolution_integral_inclined_levin_isem(a, p, x, s, l, m, k, K_sample)
    key = (a, p, x, s, l, m, k, K_sample)
    if _inclined_levin_last_key[] === key
        return _inclined_levin_last_result[]
    end

    KG = Kerr_Geodesics(a, p, 0.0, x)

    Frequencies = KG["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    omega = (m * ϒφ + k * ϒθ) / Γ
    if _skip_radiative_mode(s, a, m, omega)
        result = _zero_radiative_mode(omega, KG; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        _inclined_levin_last_key[] = key
        _inclined_levin_last_result[] = result
        return result
    end

    rsin = rstar_from_r(a, 1+sqrt(1-a^2)+1e-4)
    rsout = max(200.0, 10pi / abs(omega))
    KG_samp = kerr_geo_inclined_sample_cheby(KG, K_sample)
    KG_trap = kerr_geo_inclined_sample(KG, K_sample)
    carter_samp = carter_ingredients_sample(KG_trap, a, m, omega)
    carter_factor = trapezoidal_1d_integral(carter_samp)

    result = if s == 2
        Yup_soln = _isem_y_solution(s, l, m, a, omega)
        SH_p2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Y, Yp, X, _ = Yup_soln.Y_solution(p)
        Yup = Dict("params" => (s=2, l=l, m=m, a=a, omega=omega, lambda=SH_p2.lambda),
            "Cinc" => _isem_gsn_incidence_amplitude(Yup_soln),
            "Y" => Y,
            "Yp" => Yp,
            "X" => X)
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        Jp_up, Jm_up, dθphase_up, θphaseL_up, θphaseR_up, prefactor_up = integrand_inclined_sample_cheby_p2(KG_samp, Yup, SH_p2_samp, k)
        integralp_up = levin_1d_integral(Jp_up, dθphase_up, θphaseL_up, θphaseR_up)
        integralm_up = levin_1d_integral(Jm_up, - dθphase_up, - θphaseL_up, - θphaseR_up)
        integral_up = (integralp_up + integralm_up) * prefactor_up
        hf = horizon_factor(omega, a, m)
        Dict(
            "Amplitude" => integral_up,
            "omega" => omega,
            "EnergyFlux" => hf * abs2(integral_up),
            "AngularMomentumFlux" => hf * m * abs2(integral_up) / omega,
            "CarterConstantFlux" => 2 * hf * abs2(integral_up) * (carter_factor + k * ϒθ) / omega,
            "Trajectory" => KG,
            "YSolution" => Yup_soln,
            "SWSH" => SH_p2
        )
    elseif s == -2
        Yin_soln = _isem_y_solution(s, l, m, a, omega)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Y, Yp, X, _ = Yin_soln.Y_solution(p)
        Yin = Dict("params" => (s=-2, l=l, m=m, a=a, omega=omega, lambda=SH_m2.lambda),
            "Binc" => _isem_gsn_incidence_amplitude(Yin_soln),
            "Y" => Y,
            "Yp" => Yp,
            "X" => X)
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        Jp_in, Jm_in, dθphase_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_inclined_sample_cheby_m2(KG_samp, Yin, SH_m2_samp, k)
        integralp_in = levin_1d_integral(Jp_in, dθphase_in, θphaseL_in, θphaseR_in)
        integralm_in = levin_1d_integral(Jm_in, - dθphase_in, - θphaseL_in, - θphaseR_in)
        integral_in = (integralp_in + integralm_in) * prefactor_in
        Dict(
            "Amplitude" => integral_in,
            "omega" => omega,
            "EnergyFlux" => abs2(integral_in)/(4.0pi*omega^2),
            "AngularMomentumFlux" => m * abs2(integral_in) / (4.0pi * omega^3),
            "CarterConstantFlux" => abs2(integral_in) * (carter_factor + k * ϒθ) / (2.0pi*omega^3),
            "Trajectory" => KG,
            "YSolution" => Yin_soln,
            "SWSH" => SH_m2
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
    _inclined_levin_last_key[] = key
    _inclined_levin_last_result[] = result
    return result
end

function convolution_integral_trapezoidal_isem(a, p, e, x, s, l, m, n, k; N = 256, K = 64, Nmax::Int = 2^14, Kmax::Int = 2^12, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 0.0, zero_low_flux::Bool = false, threaded_sampling::Bool = false)
    KG = Kerr_Geodesics(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end
    if m == 0 && n == 0 && k == 0
        return Dict("Amplitude" => 0.0 + 0.0im, "omega" => 0.0, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
    end
    if isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0 || k != 0
            return Dict("Amplitude" => 0.0 + 0.0im, "omega" => nothing, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
        end
        return convolution_integral_circular_equatorial_isem(a, p, s, l, m)
    elseif isapprox(e, 0.0; atol=1e-12) && !isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0
            return Dict("Amplitude" => 0.0 + 0.0im, "omega" => nothing, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
        end
        return convolution_integral_inclined_trapezoidal_isem(a, p, x, s, l, m, k, K; Kmax = Kmax, tol = tol, sample_tol = sample_tol, max_flux = max_flux, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
    elseif !isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if k != 0
            return Dict("Amplitude" => 0.0 + 0.0im, "omega" => nothing, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
        end
        return convolution_integral_eccentric_trapezoidal_isem(a, p, e, s, l, m, n, N; Nmax = Nmax, tol = tol, sample_tol = sample_tol, max_flux = max_flux, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
    else
        return generic_mode_flux(a, p, e, x, s, l, m, n, k; N0 = N, K0 = K, Nmax = Nmax, Kmax = Kmax, tol = tol, sample_tol = sample_tol, max_flux = max_flux, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
    end
end

function convolution_integral_levin_isem(a, p, e, x, s, l, m, n, k; N = 256, K = 32)
    KG = Kerr_Geodesics(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end
    if m == 0 && n == 0 && k == 0
        return Dict("Amplitude" => 0.0 + 0.0im, "omega" => 0.0, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
    end
    if isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0 || k != 0
            return Dict("Amplitude" => 0.0 + 0.0im, "omega" => nothing, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
        end
        return convolution_integral_circular_equatorial_isem(a, p, s, l, m)
    elseif isapprox(e, 0.0; atol=1e-12) && !isapprox(abs(x), 1.0; atol=1e-12)
        if n != 0
            return Dict("Amplitude" => 0.0 + 0.0im, "omega" => nothing, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
        end
        return convolution_integral_inclined_levin_isem(a, p, x, s, l, m, k, K)
    elseif !isapprox(e, 0.0; atol=1e-12) && isapprox(abs(x), 1.0; atol=1e-12)
        if k != 0
            return Dict("Amplitude" => 0.0 + 0.0im, "omega" => nothing, "EnergyFlux" => 0.0, "AngularMomentumFlux" => 0.0, "CarterConstantFlux" => 0.0, "Trajectory" => KG, "YSolution" => nothing, "SWSH" => nothing)
        end
        return convolution_integral_eccentric_levin_isem(a, p, e, s, l, m, n, N)
    else
        return convolution_integral_generic_levin_isem(a, p, e, x, s, l, m, n, k, N, K)
    end
end

end
