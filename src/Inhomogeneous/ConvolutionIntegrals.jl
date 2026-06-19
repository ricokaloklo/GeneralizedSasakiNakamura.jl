module ConvolutionIntegrals

using LinearAlgebra
using SpinWeightedSpheroidalHarmonics
using KerrGeodesics

using ..ISEM
using ..Coordinates
using ..SolutionsY
using ..GridSampling

export convolution_integral_trapezoidal, convolution_integral_levin
export convolution_integral_trapezoidal_isem, convolution_integral_levin_isem
export convolution_integral_circular_equatorial_isem
export convolution_integral_eccentric_adaptive_levin_isem
export convolution_integral_generic_adaptive_levin_isem

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
    levin_samples::Dict{Int, Dict}
    adaptive_levin_segments::Dict{Tuple{Int, Int, Int}, Any}
    phase_vectors::Dict{Tuple{Int, Int, Int, UInt64}, Any}
    levin_phase_factors::Dict{Tuple{Int, Int, Int, Int, UInt64}, Any}
    levin_factored_phase_factors::Dict{Tuple{Int, Int, Int, Int, UInt64}, Any}
end

mutable struct InclinedFluxCache
    samples::Dict{Int, Dict}
end

EccentricFluxCache() = EccentricFluxCache(Dict{Int, Dict}(), Dict{Int, Dict}(), Dict{Tuple{Int, Int, Int}, Any}(), Dict{Tuple{Int, Int, Int, UInt64}, Any}(), Dict{Tuple{Int, Int, Int, Int, UInt64}, Any}(), Dict{Tuple{Int, Int, Int, Int, UInt64}, Any}())
EccentricFluxCache(samples::Dict{Int, Dict}) = EccentricFluxCache(samples, Dict{Int, Dict}(), Dict{Tuple{Int, Int, Int}, Any}(), Dict{Tuple{Int, Int, Int, UInt64}, Any}(), Dict{Tuple{Int, Int, Int, Int, UInt64}, Any}(), Dict{Tuple{Int, Int, Int, Int, UInt64}, Any}())
InclinedFluxCache() = InclinedFluxCache(Dict{Int, Dict}())

function _cached_eccentric_sample!(cache::EccentricFluxCache, KG_master::Dict, N::Int)
    return get!(cache.samples, N) do
        subsample_eccentric_sample(KG_master, N)
    end
end

function _kg_from_presampled_master(KG_sample::Dict)
    cross = KG_sample["CrossFunction"]
    cross4 = length(cross) == 2 ? (cross[1], nothing, cross[2], nothing) : cross
    dcross = KG_sample["DerivativesCrossFunction"]
    dcross4 = length(dcross) == 2 ? (dcross[1], nothing, dcross[2], nothing) : dcross
    kg = Dict{String, Any}(
        "a" => KG_sample["a"],
        "p" => KG_sample["p"],
        "Energy" => KG_sample["E"],
        "AngularMomentum" => KG_sample["Lz"],
        "Frequencies" => KG_sample["Frequencies"],
        "Trajectory" => KG_sample["Trajectory"],
        "FourVelocity" => KG_sample["FourVelocity"],
        "CrossFunction" => cross4,
        "DerivativesCrossFunction" => dcross4,
        "InitialPhases" => KG_sample["InitialPhases"],
    )
    haskey(KG_sample, "e") && (kg["e"] = KG_sample["e"])
    haskey(KG_sample, "x") && (kg["x"] = KG_sample["x"])
    return kg
end

function _cached_eccentric_cheby_sample!(cache::EccentricFluxCache, KG_master::Dict, N::Int)
    return get!(cache.levin_samples, N) do
        kerr_geo_eccentric_sample_cheby(_kg_from_presampled_master(KG_master), N)
    end
end

function _pow2_refinement_levels(start::Int, stop::Int)
    start <= stop || return Int[]
    levels = Int[]
    n = start
    while true
        push!(levels, n)
        n >= stop && break
        n = min(2n, stop)
    end
    return levels
end

function prewarm_eccentric_levin_samples!(cache::EccentricFluxCache, KG_master::Dict; N0::Int, Nmax::Int)
    for N in _pow2_refinement_levels(N0, Nmax)
        _cached_eccentric_cheby_sample!(cache, KG_master, N)
        _levin_1d_local_plan(N, DEFAULT_LEVIN_LOCAL_ORDER)
    end
    return cache
end

function _cached_inclined_sample!(cache::InclinedFluxCache, KG_master::Dict, K::Int)
    return get!(cache.samples, K) do
        subsample_inclined_sample(KG_master, K)
    end
end

function _generic_trapezoidal_context(a, p, e, x, s, l, m, n, k, N_sample, K_sample)
    key = (a, p, e, x, s, l, m, n, k, N_sample, K_sample)
    return get!(_generic_trapezoidal_cache, key) do
        KG = kerr_geo_orbit(a, p, e, x)
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
        KG = kerr_geo_orbit(a, p, e, 1.0)
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
            (s == 2 ? "Cinc" : "Binc") => GridSampling._isem_gsn_incidence_amplitude(Ysol),
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
    KG = kerr_geo_orbit(a, p, e, x)
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
    KG = kerr_geo_orbit(a, p, e, x)
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
        rphase = _radial_phase_vector(KG_samp, omega, m, n)
        θphase = _polar_phase_vector(KG_samp, omega, m, k)
        Jpp_grid_up, Jpm_grid_up, Jmp_grid_up, Jmm_grid_up, drphase_up, dθphase_up, rphaseL_up, rphaseR_up, θphaseL_up, θphaseR_up, prefactor_up = integrand_generic_sample_cheby_p2(KG_samp, Yup_samp, SH_p2_samp, n, k)
        integralpp_up = levin_2d_integral_local_phase(Jpp_grid_up, drphase_up, dθphase_up, rphase, θphase)
        integralpm_up = levin_2d_integral_local_phase(Jpm_grid_up, drphase_up, - dθphase_up, rphase, -θphase)
        integralmp_up = levin_2d_integral_local_phase(Jmp_grid_up, -drphase_up, dθphase_up, -rphase, θphase)
        integralmm_up = levin_2d_integral_local_phase(Jmm_grid_up, -drphase_up, -dθphase_up, -rphase, -θphase)
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
        rphase = _radial_phase_vector(KG_samp, omega, m, n)
        θphase = _polar_phase_vector(KG_samp, omega, m, k)
        Jpp_grid_in, Jpm_grid_in, Jmp_grid_in, Jmm_grid_in, drphase_in, dθphase_in, rphaseL_in, rphaseR_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_generic_sample_cheby_m2(KG_samp, Yin_samp, SH_m2_samp, n, k)
        integralpp_in = levin_2d_integral_local_phase(Jpp_grid_in, drphase_in, dθphase_in, rphase, θphase)
        integralpm_in = levin_2d_integral_local_phase(Jpm_grid_in, drphase_in, -dθphase_in, rphase, -θphase)
        integralmp_in = levin_2d_integral_local_phase(Jmp_grid_in, -drphase_in, dθphase_in, -rphase, θphase)
        integralmm_in = levin_2d_integral_local_phase(Jmm_grid_in, -drphase_in, -dθphase_in, -rphase, -θphase)
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

function _generic_flux_from_cheby_sample(KG_samp::Dict, Y_samp::Dict, SH_samp::Dict, s::Int, l::Int, m::Int, n::Int, k::Int, a::Float64, ω::Float64, ϒθ::Float64, carter_factor::Float64, cache = nothing)
    N_sample = KG_samp["N_sample"]::Int
    K_sample = KG_samp["K_sample"]::Int
    phase_key = (s, N_sample, K_sample, m, n, k)
    if s == 2
        Jpp, Jpm, Jmp, Jmm, drphase, dθphase, _, _, _, _, prefactor = integrand_generic_sample_cheby_p2(KG_samp, Y_samp, SH_samp, n, k)
        basis_plus = _glevin_basis!(cache, KG_samp, s, Float64(a), 1)
        basis_minus = _glevin_basis!(cache, KG_samp, s, Float64(a), -1)
        radial_factor_vec = _cached_generic_levin_phase_factor!(cache, (:internal_radial_factor, phase_key)) do
            exp.(1im .* _internal_radial_phase(KG_samp, s, a, ω, m).phase_values)
        end
        radial_factor = reshape(radial_factor_vec, :, 1)
        θfactor = _cached_generic_levin_phase_factor!(cache, (:theta, phase_key)) do
            θphase_data = _glevin_theta(basis_plus, m, n, k)
            _levin_1d_local_phase_factor(θphase_data.derivative, θphase_data.phase)
        end
        θfactor_neg = _cached_generic_levin_phase_factor!(cache, (:theta_neg, phase_key)) do
            _conjugate_levin_1d_local_phase_factor(θfactor)
        end
        rfactor_plus = _cached_generic_levin_phase_factor!(cache, (:radial_plus, phase_key)) do
            rphase_plus = _glevin_radial(basis_plus, m, n, k)
            _levin_1d_local_phase_factor_from_values(rphase_plus)
        end
        rfactor_minus = _cached_generic_levin_phase_factor!(cache, (:radial_minus, phase_key)) do
            rphase_minus = _glevin_radial(basis_minus, m, n, k)
            _levin_1d_local_phase_factor_from_values(rphase_minus)
        end
        integralpp, integralmp = levin_2d_integral_local_phase_pair_batched_radial_values_factored(Jpp ./ radial_factor, Jmp ./ radial_factor, θfactor, rfactor_plus, rfactor_minus)
        integralpm, integralmm = levin_2d_integral_local_phase_pair_batched_radial_values_factored(Jpm ./ radial_factor, Jmm ./ radial_factor, θfactor_neg, rfactor_plus, rfactor_minus)
        integral = (integralpp + integralpm + integralmp + integralmm) * prefactor
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
        Jpp, Jpm, Jmp, Jmm, drphase, dθphase, _, _, _, _, prefactor = integrand_generic_sample_cheby_m2(KG_samp, Y_samp, SH_samp, n, k)
        basis_plus = _glevin_basis!(cache, KG_samp, s, Float64(a), 1)
        basis_minus = _glevin_basis!(cache, KG_samp, s, Float64(a), -1)
        radial_factor_vec = _cached_generic_levin_phase_factor!(cache, (:internal_radial_factor, phase_key)) do
            exp.(1im .* _internal_radial_phase(KG_samp, s, a, ω, m).phase_values)
        end
        radial_factor = reshape(radial_factor_vec, :, 1)
        θfactor = _cached_generic_levin_phase_factor!(cache, (:theta, phase_key)) do
            θphase_data = _glevin_theta(basis_plus, m, n, k)
            _levin_1d_local_phase_factor(θphase_data.derivative, θphase_data.phase)
        end
        θfactor_neg = _cached_generic_levin_phase_factor!(cache, (:theta_neg, phase_key)) do
            _conjugate_levin_1d_local_phase_factor(θfactor)
        end
        rfactor_plus = _cached_generic_levin_phase_factor!(cache, (:radial_plus, phase_key)) do
            rphase_plus = _glevin_radial(basis_plus, m, n, k)
            _levin_1d_local_phase_factor_from_values(rphase_plus)
        end
        rfactor_minus = _cached_generic_levin_phase_factor!(cache, (:radial_minus, phase_key)) do
            rphase_minus = _glevin_radial(basis_minus, m, n, k)
            _levin_1d_local_phase_factor_from_values(rphase_minus)
        end
        integralpp, integralmp = levin_2d_integral_local_phase_pair_batched_radial_values_factored(Jpp ./ radial_factor, Jmp ./ radial_factor, θfactor, rfactor_plus, rfactor_minus)
        integralpm, integralmm = levin_2d_integral_local_phase_pair_batched_radial_values_factored(Jpm ./ radial_factor, Jmm ./ radial_factor, θfactor_neg, rfactor_plus, rfactor_minus)
        integral = (integralpp + integralpm + integralmp + integralmm) * prefactor
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
    levin_samples::Dict{Tuple{Int, Int}, Dict}
    levin_phase_factors::Dict{Tuple, Any}
    levin_phase_bases::Dict{Tuple{Int, Int, Int, Int}, Any}
    adaptive_levin_segments::Dict{Tuple, Dict}
    adaptive_levin_phase_bases::Dict{Tuple, Any}
    adaptive_levin_weights::Dict{Tuple, Vector{Float64}}
end

GenericM2FluxCache() = GenericM2FluxCache(
    Dict{Tuple{Int, Int}, Dict}(),
    Dict{Tuple{Int, Int, Int}, Any}(),
    Dict{Tuple{Int, Int}, Dict}(),
    Dict{Tuple, Any}(),
    Dict{Tuple{Int, Int, Int, Int}, Any}(),
    Dict{Tuple, Dict}(),
    Dict{Tuple, Any}(),
    Dict{Tuple, Vector{Float64}}(),
)
GenericM2FluxCache(samples::Dict{Tuple{Int, Int}, Dict}, mcaches::Dict{Tuple{Int, Int, Int}, Any}) = GenericM2FluxCache(
    samples,
    mcaches,
    Dict{Tuple{Int, Int}, Dict}(),
    Dict{Tuple, Any}(),
    Dict{Tuple{Int, Int, Int, Int}, Any}(),
    Dict{Tuple, Dict}(),
    Dict{Tuple, Any}(),
    Dict{Tuple, Vector{Float64}}(),
)
GenericFluxCache() = GenericM2FluxCache()

struct GenericLevinPhaseBasis
    theta_m::Vector{Float64}
    theta_n::Vector{Float64}
    theta_k::Vector{Float64}
    dtheta_m::Vector{Float64}
    dtheta_n::Vector{Float64}
    dtheta_k::Vector{Float64}
    radial_m::Vector{Float64}
    radial_n::Vector{Float64}
    radial_k::Vector{Float64}
end

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

function _cached_generic_cheby_sample!(cache::GenericM2FluxCache, KG_master::Dict, N::Int, K::Int)
    key = (N, K)
    return get!(cache.levin_samples, key) do
        kerr_geo_generic_sample_cheby(_kg_from_presampled_master(KG_master), N, K)
    end
end

function _cached_generic_levin_phase_factor!(build::Function, cache, key::Tuple)
    cache isa GenericM2FluxCache || return build()
    return get!(cache.levin_phase_factors, key) do
        build()
    end
end

function _glevin_freqs(KG_samp::Dict)
    freqs = KG_samp["Frequencies"]
    Γ = freqs["ϒt"]
    return (Ωφ = freqs["ϒϕ"] / Γ, Ωr = freqs["ϒr"] / Γ, Ωθ = freqs["ϒθ"] / Γ)
end

function _glevin_basis_build(KG_samp::Dict, s::Int, a::Float64, sign::Int)
    Ω = _glevin_freqs(KG_samp)
    qr = _radial_phase_nodes(KG_samp)
    qθ = _cheby_phase_nodes_pi(length(KG_samp["Δtθ"]))
    r = KG_samp["r"]
    rs = KG_samp["rs"]
    κ = sqrt(1.0 - a^2)
    rp = 1.0 + κ
    rm = 1.0 - κ
    log_term = log.((r .- rp) ./ (r .- rm))
    σrs = s == -2 ? 1.0 : -1.0
    σlog = s == -2 ? -1.0 : 1.0
    radial_base = sign .* KG_samp["Δtr"] .+ σrs .* rs
    return GenericLevinPhaseBasis(
        Ω.Ωφ .* KG_samp["Δtθ"] .- KG_samp["Δφθ"],
        Ω.Ωr .* KG_samp["Δtθ"],
        Ω.Ωθ .* KG_samp["Δtθ"] .+ qθ,
        Ω.Ωφ .* KG_samp["dtθ"] .- KG_samp["dφθ"],
        Ω.Ωr .* KG_samp["dtθ"],
        Ω.Ωθ .* KG_samp["dtθ"] .+ 1.0,
        Ω.Ωφ .* radial_base .- sign .* KG_samp["Δφr"] .+ σlog .* a .* log_term ./ (2.0κ),
        Ω.Ωr .* radial_base .+ sign .* qr,
        Ω.Ωθ .* radial_base,
    )
end

function _glevin_basis!(cache, KG_samp::Dict, s::Int, a::Float64, sign::Int)
    cache isa GenericM2FluxCache || return _glevin_basis_build(KG_samp, s, a, sign)
    key = (KG_samp["N_sample"]::Int, KG_samp["K_sample"]::Int, s, sign)
    return get!(cache.levin_phase_bases, key) do
        _glevin_basis_build(KG_samp, s, a, sign)
    end
end

function _glevin_theta(basis::GenericLevinPhaseBasis, m::Int, n::Int, k::Int)
    return (
        phase = m .* basis.theta_m .+ n .* basis.theta_n .+ k .* basis.theta_k,
        derivative = m .* basis.dtheta_m .+ n .* basis.dtheta_n .+ k .* basis.dtheta_k,
    )
end

function _glevin_radial(basis::GenericLevinPhaseBasis, m::Int, n::Int, k::Int)
    return m .* basis.radial_m .+ n .* basis.radial_n .+ k .* basis.radial_k
end

function _paired_pow2_refinement_levels(N0::Int, K0::Int, Nmax::Int, Kmax::Int)
    N0 <= Nmax || return Tuple{Int, Int}[]
    K0 <= Kmax || return Tuple{Int, Int}[]
    levels = Tuple{Int, Int}[]
    N = N0
    K = K0
    while true
        push!(levels, (N, K))
        (N >= Nmax && K >= Kmax) && break
        N = min(2N, Nmax)
        K = min(2K, Kmax)
    end
    return levels
end

function prewarm_generic_levin_samples!(cache::GenericM2FluxCache, KG_master::Dict; N0::Int, K0::Int, Nmax::Int, Kmax::Int)
    for (N, K) in _paired_pow2_refinement_levels(N0, K0, Nmax, Kmax)
        _cached_generic_cheby_sample!(cache, KG_master, N, K)
        _levin_1d_local_plan(N, DEFAULT_LEVIN_LOCAL_ORDER)
        _levin_1d_local_plan(K, DEFAULT_LEVIN_LOCAL_ORDER)
    end
    return cache
end

struct GenericAdaptiveLevin2DKey
    level_r::Int
    level_theta::Int
    ib::Int
    jb::Int
end

struct GenericAdaptiveLevin2DLeaf
    key::GenericAdaptiveLevin2DKey
    amp::ComplexF64
    active::Bool
end

@inline _generic_adaptive_radial_children(key::GenericAdaptiveLevin2DKey) = (
    GenericAdaptiveLevin2DKey(key.level_r + 1, key.level_theta, 2 * key.ib, key.jb),
    GenericAdaptiveLevin2DKey(key.level_r + 1, key.level_theta, 2 * key.ib + 1, key.jb),
)

@inline function _generic_adaptive_segment_bounds(key::GenericAdaptiveLevin2DKey)
    sr = 2.0^key.level_r
    st = 2.0^key.level_theta
    return (
        π * key.ib / sr,
        π * (key.ib + 1) / sr,
        π * key.jb / st,
        π * (key.jb + 1) / st,
    )
end

function _generic_adaptive_cheby_nodes(n::Int, lo::Float64, hi::Float64)
    n >= 2 || throw(ArgumentError("local Chebyshev node count must be at least 2"))
    xref = [cos(π * k / (n - 1)) for k in 0:(n - 1)]
    return (hi - lo) .* (xref .+ 1.0) ./ 2.0 .+ lo
end

function _generic_adaptive_cc_weights(n::Int, lo::Float64, hi::Float64)
    n >= 2 || return zeros(Float64, n)
    N = n - 1
    θ = [π * j / N for j in 0:N]
    weights = zeros(Float64, n)
    if N == 1
        weights .= 1.0
    else
        v = ones(Float64, N - 1)
        if iseven(N)
            weights[1] = 1.0 / (N^2 - 1.0)
            weights[end] = weights[1]
            for k in 1:(N ÷ 2 - 1)
                @inbounds for j in 1:(N - 1)
                    v[j] -= 2.0 * cos(2.0 * k * θ[j + 1]) / (4.0 * k^2 - 1.0)
                end
            end
            @inbounds for j in 1:(N - 1)
                v[j] -= cos(N * θ[j + 1]) / (N^2 - 1.0)
            end
        else
            weights[1] = 1.0 / N^2
            weights[end] = weights[1]
            for k in 1:((N - 1) ÷ 2)
                @inbounds for j in 1:(N - 1)
                    v[j] -= 2.0 * cos(2.0 * k * θ[j + 1]) / (4.0 * k^2 - 1.0)
                end
            end
        end
        @inbounds for j in 1:(N - 1)
            weights[j + 1] = 2.0 * v[j] / N
        end
    end
    weights .*= 0.5 * (hi - lo)
    return weights
end

function _generic_adaptive_segment_sample_cheby(KG_in::Dict, key::GenericAdaptiveLevin2DKey, nr::Int, nt::Int)
    KG = haskey(KG_in, "Energy") ? KG_in : _kg_from_presampled_master(KG_in)
    qr_lo, qr_hi, qθ_lo, qθ_hi = _generic_adaptive_segment_bounds(key)
    a = Float64(KG["a"])
    _, r, θ, _ = KG["Trajectory"]
    _, ur, uθ, _ = KG["FourVelocity"]
    Δtr, Δtθ, Δφr, Δφθ = KG["CrossFunction"]
    dtr, dtθ, dφr, dφθ = KG["DerivativesCrossFunction"]
    qt0, qr0, qθ0, qφ0 = KG["InitialPhases"]
    freqs = KG["Frequencies"]
    ϒr = freqs["ϒr"]
    ϒθ = freqs["ϒθ"]
    Γ = freqs["ϒt"]

    qr_vals = _generic_adaptive_cheby_nodes(nr, qr_lo, qr_hi)
    qθ_vals = _generic_adaptive_cheby_nodes(nt, qθ_lo, qθ_hi)
    rq(qr) = r((qr - qr0) / ϒr)
    θq(qθ) = θ((qθ - qθ0) / ϒθ)
    urq(qr) = (r((qr - qr0) / ϒr)^2 + a^2 * cos(θ((qr - qr0) / ϒr))^2) * ur((qr - qr0) / ϒr)
    uθq(qθ) = (r((qθ - qθ0) / ϒθ)^2 + a^2 * cos(θ((qθ - qθ0) / ϒθ))^2) * uθ((qθ - qθ0) / ϒθ)
    theta_odd(f, q) = q >= π / 2 ? f(q) : -f(π - q)
    theta_odd_derivative(f, q) = q >= π / 2 ? f(q) : f(π - q)
    theta_even(f, q) = q >= π / 2 ? f(q) : f(π - q)

    r_vals = [rq(q) for q in qr_vals]
    θ_vals = [θq(q) for q in qθ_vals]
    ur_fwd = [urq(q) for q in qr_vals]
    uθ_fwd = [theta_even(uθq, q) for q in qθ_vals]

    return Dict{String, Any}(
        "qr" => qr_vals,
        "qθ" => qθ_vals,
        "r" => r_vals,
        "rs" => rstar_from_r.(a, r_vals),
        "θ" => θ_vals,
        "Δtr" => [Δtr(q) for q in qr_vals],
        "Δtθ" => [theta_odd(Δtθ, q) for q in qθ_vals],
        "Δφr" => [Δφr(q) for q in qr_vals],
        "Δφθ" => [theta_odd(Δφθ, q) for q in qθ_vals],
        "dtr" => [dtr(q) for q in qr_vals],
        "dtθ" => [theta_odd_derivative(dtθ, q) for q in qθ_vals],
        "dφr" => [dφr(q) for q in qr_vals],
        "dφθ" => [theta_odd_derivative(dφθ, q) for q in qθ_vals],
        "ur_fwd" => ur_fwd,
        "ur_rev" => -ur_fwd,
        "uθ_fwd" => uθ_fwd,
        "uθ_rev" => -uθ_fwd,
        "N_sample" => nr,
        "K_sample" => nt,
        "a" => a,
        "p" => get(KG, "p", NaN),
        "e" => get(KG, "e", NaN),
        "x" => get(KG, "x", NaN),
        "E" => KG["Energy"],
        "Lz" => KG["AngularMomentum"],
        "Γ" => Γ,
        "Frequencies" => freqs,
        "Trajectory" => KG["Trajectory"],
        "FourVelocity" => KG["FourVelocity"],
        "CrossFunction" => KG["CrossFunction"],
        "DerivativesCrossFunction" => KG["DerivativesCrossFunction"],
        "InitialPhases" => (qt0, qr0, qθ0, qφ0),
        "initialPhases" => (qt0, qr0, qθ0, qφ0),
    )
end

function _generic_adaptive_segment_sample!(cache::GenericM2FluxCache, KG::Dict, key::GenericAdaptiveLevin2DKey, nr::Int, nt::Int)
    return get!(cache.adaptive_levin_segments, (key.level_r, key.level_theta, key.ib, key.jb, nr, nt)) do
        _generic_adaptive_segment_sample_cheby(KG, key, nr, nt)
    end
end

function _generic_adaptive_phase_basis_build(KG_samp::Dict, s::Int, a::Float64, sign::Int)
    Ω = _glevin_freqs(KG_samp)
    qr = Float64.(KG_samp["qr"])
    qθ = Float64.(KG_samp["qθ"])
    r = Float64.(KG_samp["r"])
    rs = Float64.(KG_samp["rs"])
    κ = sqrt(1.0 - a^2)
    rp = 1.0 + κ
    rm = 1.0 - κ
    log_term = log.((r .- rp) ./ (r .- rm))
    σrs = s == -2 ? 1.0 : -1.0
    σlog = s == -2 ? -1.0 : 1.0
    radial_base = sign .* KG_samp["Δtr"] .+ σrs .* rs
    return GenericLevinPhaseBasis(
        Ω.Ωφ .* KG_samp["Δtθ"] .- KG_samp["Δφθ"],
        Ω.Ωr .* KG_samp["Δtθ"],
        Ω.Ωθ .* KG_samp["Δtθ"] .+ qθ,
        Ω.Ωφ .* KG_samp["dtθ"] .- KG_samp["dφθ"],
        Ω.Ωr .* KG_samp["dtθ"],
        Ω.Ωθ .* KG_samp["dtθ"] .+ 1.0,
        Ω.Ωφ .* radial_base .- sign .* KG_samp["Δφr"] .+ σlog .* a .* log_term ./ (2.0κ),
        Ω.Ωr .* radial_base .+ sign .* qr,
        Ω.Ωθ .* radial_base,
    )
end

function _generic_adaptive_phase_basis!(cache::GenericM2FluxCache, KG_samp::Dict, key::GenericAdaptiveLevin2DKey, s::Int, a::Float64, sign::Int)
    nr = KG_samp["N_sample"]::Int
    nt = KG_samp["K_sample"]::Int
    return get!(cache.adaptive_levin_phase_bases, (key.level_r, key.level_theta, key.ib, key.jb, nr, nt, s, sign)) do
        _generic_adaptive_phase_basis_build(KG_samp, s, a, sign)
    end
end

@inline function _generic_adaptive_effective_theta_intervals(local_theta_intervals::Int, k::Int)
    required = 2 ^ ceil(Int, log2(max(1, 2 * abs(k) + 1)))
    return max(local_theta_intervals, required)
end

function _generic_adaptive_theta_weights!(cache::GenericM2FluxCache, key::GenericAdaptiveLevin2DKey, nt::Int)
    _, _, qθ_lo, qθ_hi = _generic_adaptive_segment_bounds(key)
    return get!(cache.adaptive_levin_weights, (:theta_cc, key.level_theta, key.jb, nt)) do
        _generic_adaptive_cc_weights(nt, qθ_lo, qθ_hi)
    end
end

function _generic_adaptive_radial_levin_theta_cc(
        f1::AbstractMatrix{ComplexF64},
        f2::AbstractMatrix{ComplexF64},
        theta_weights::AbstractVector{Float64},
        theta_phase::AbstractVector{<:Number},
        radial_factor1,
        radial_factor2)
    nr, nt = size(f1)
    size(f2) == (nr, nt) || throw(ArgumentError("paired radial Levin inputs must have matching size"))
    length(theta_weights) == nt || throw(ArgumentError("theta weights length mismatch"))
    length(theta_phase) == nt || throw(ArgumentError("theta phase length mismatch"))
    out1 = zero(ComplexF64)
    out2 = zero(ComplexF64)
    @inbounds for j in 1:nt
        phasej = exp(1im * theta_phase[j])
        out1 += theta_weights[j] * phasej * _levin_1d_integral_local_phase_factored(collect(@view f1[:, j]), radial_factor1)
        out2 += theta_weights[j] * phasej * _levin_1d_integral_local_phase_factored(collect(@view f2[:, j]), radial_factor2)
    end
    return out1, out2
end

@inline function _generic_adaptive_energy_from_amp(amp::ComplexF64, s::Int, a::Float64, m::Int, ω::Float64)
    if s == -2
        return abs2(amp) / (4.0pi * ω^2)
    elseif s == 2
        return horizon_factor(ω, a, m) * abs2(amp)
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function _generic_adaptive_segment_amplitude!(cache::GenericM2FluxCache, KG::Dict, Ysol, SH, s::Int, a::Float64, ω::Float64, m::Int, n::Int, k::Int, key::GenericAdaptiveLevin2DKey, nr::Int, nt::Int; threaded_sampling::Bool = false)
    KG_samp = _generic_adaptive_segment_sample!(cache, KG, key, nr, nt)
    Ysamp = if s == -2
        threaded_sampling ? GridSampling.y_sample_m2_isem_threaded(Ysol, KG_samp) : GridSampling.y_sample_m2_isem(Ysol, KG_samp)
    else
        threaded_sampling ? GridSampling.y_sample_p2_isem_threaded(Ysol, KG_samp) : GridSampling.y_sample_p2_isem(Ysol, KG_samp)
    end
    SHsamp = threaded_sampling ? GridSampling.swsh_sample_threaded(SH, KG_samp) : GridSampling.swsh_sample(SH, KG_samp)
    Jpp, Jpm, Jmp, Jmm, _, _, _, _, _, _, prefactor = if s == -2
        integrand_generic_sample_cheby_m2(KG_samp, Ysamp, SHsamp, n, k)
    else
        integrand_generic_sample_cheby_p2(KG_samp, Ysamp, SHsamp, n, k)
    end

    basis_plus = _generic_adaptive_phase_basis!(cache, KG_samp, key, s, a, 1)
    basis_minus = _generic_adaptive_phase_basis!(cache, KG_samp, key, s, a, -1)
    radial_internal = reshape(exp.(1im .* _internal_radial_phase(KG_samp, s, a, ω, m).phase_values), :, 1)
    rfactor_plus = _cached_generic_levin_phase_factor!(cache, (:adaptive_radial, key.level_r, key.ib, s, m, n, k, 1, nr)) do
        _levin_1d_local_phase_factor_from_values(_glevin_radial(basis_plus, m, n, k); local_order = nr)
    end
    rfactor_minus = _cached_generic_levin_phase_factor!(cache, (:adaptive_radial, key.level_r, key.ib, s, m, n, k, -1, nr)) do
        _levin_1d_local_phase_factor_from_values(_glevin_radial(basis_minus, m, n, k); local_order = nr)
    end
    theta_weights = _generic_adaptive_theta_weights!(cache, key, nt)
    theta_phase = _glevin_theta(basis_plus, m, n, k).phase
    Fa = ComplexF64.(Jpp ./ radial_internal)
    Fb = ComplexF64.(Jmp ./ radial_internal)
    integralpp, integralmp = _generic_adaptive_radial_levin_theta_cc(Fa, Fb, theta_weights, theta_phase, rfactor_plus, rfactor_minus)
    Fa = ComplexF64.(Jpm ./ radial_internal)
    Fb = ComplexF64.(Jmm ./ radial_internal)
    integralpm, integralmm = _generic_adaptive_radial_levin_theta_cc(Fa, Fb, theta_weights, -theta_phase, rfactor_plus, rfactor_minus)
    qr_lo, qr_hi, _, _ = _generic_adaptive_segment_bounds(key)
    return ComplexF64((integralpp + integralpm + integralmp + integralmm) * prefactor * ((qr_hi - qr_lo) / π))
end

function prewarm_generic_adaptive_levin_segments!(cache::GenericM2FluxCache, KG_master::Dict;
        local_r_intervals::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N,
        local_theta_intervals::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N,
        max_depth::Int = 8)
    nr = local_r_intervals + 1
    nt = local_theta_intervals + 1
    for level in 0:max_depth
        for ib in 0:(2^level - 1)
            key = GenericAdaptiveLevin2DKey(level, 0, ib, 0)
            _generic_adaptive_segment_sample!(cache, KG_master, key, nr, nt)
            _generic_adaptive_theta_weights!(cache, key, nt)
        end
    end
    return cache
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

@inline function _levin_low_flux_budget(max_flux::Float64, tol::Float64, mode_abs_floor::Float64)
    return max(mode_abs_floor, max(max_flux, eps(Float64)) * tol)
end

@inline function _levin_low_flux_stable(current::Float64, previous::Float64, floor::Float64)
    ac = abs(current)
    ap = abs(previous)
    max(ac, ap) <= floor && return true
    minv = min(ac, ap)
    minv <= floor && return false
    return max(ac, ap) / minv <= 2.0
end

@inline function _same_order_or_smaller(current::Float64, previous::Float64, floor::Float64)
    ac = abs(current)
    ap = abs(previous)
    max(ac, ap) <= floor && return true
    minv = min(ac, ap)
    minv <= floor && return false
    return ac <= 10.0 * ap
end

@inline function _suspect_min_sample(min_stop::Int, max_stop::Int, excess::Float64, level100::Int, level1e4::Int, level1e8::Int)
    target = min_stop
    if excess > 1e8
        target = max(target, level1e8)
    elseif excess > 1e4
        target = max(target, level1e4)
    elseif excess > 100.0
        target = max(target, level100)
    elseif excess > 10.0
        target = max(target, 2 * min_stop)
    end

    return min(target, max_stop)
end

@inline function _eccentric_fourier_min_sample(min_stop::Int, max_stop::Int, n::Integer)
    harmonic = abs(n)
    harmonic == 0 && return min_stop
    target = nextpow(2, max(min_stop, 4 * harmonic))
    return min(target, max_stop)
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

@inline function _generic_carter_factor_sample(KG_samp::Dict, m::Int, omega::Float64)
    theta = KG_samp["θ"]::Vector{Float64}
    K = length(theta)
    K <= 1 && return 0.0
    a = KG_samp["a"]::Float64
    E = KG_samp["E"]::Float64
    Lz = KG_samp["Lz"]::Float64
    dθ = pi / (K - 1)
    a2ωE = a * a * omega * E
    mLz = m * Lz
    total = 0.0
    @inbounds for j in 1:K
        st = sin(theta[j])
        ct = cos(theta[j])
        weight = (j == 1 || j == K) ? 0.5 * dθ : dθ
        total += weight * (mLz * ct * ct / (st * st) - a2ωE * ct * ct) / pi
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

function generic_mode_flux_from_master_cached!(cache::GenericM2FluxCache, KG_master::Dict, s::Int, l::Int, m::Int, n::Int, k::Int; N0::Int = 64, K0::Int = 16, Nmax::Int = 2^14, Kmax::Int = 2^12, sample_tol::Float64 = 1e-3, tol::Float64 = 1e-8, max_flux::Float64 = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false)
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
        suspect_min_N = has_flux_reference && !low_flux_done ? _suspect_min_sample(min_stop_N, Nmax, excess, 1024, 2048, 4096) : min_stop_N
        suspect_min_K = has_flux_reference && !low_flux_done ? _suspect_min_sample(min_stop_K, Kmax, excess, 256, 512, 512) : min_stop_K
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

function generic_mode_flux_from_master_cached_levin!(cache::GenericM2FluxCache, KG_master::Dict, s::Int, l::Int, m::Int, n::Int, k::Int; N0::Int = 256, K0::Int = 64, Nmax::Int = DEFAULT_LEVIN_NMAX, Kmax::Int = DEFAULT_LEVIN_KMAX, sample_tol::Float64 = 1e-3, tol::Float64 = 1e-8, max_flux::Float64 = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false, confirm_low_flux::Bool = true)
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
        res["Quadrature"] = "levin"
        return res
    end

    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * ω; method = "jacobi")
    Ysol = _generic_isem_y_solution(s, l, m, a, ω)
    flux_scale = max(max_flux, eps(Float64))
    effective_sample_tol = min(sample_tol, 3.0 * sqrt(tol))
    low_flux_budget = _levin_low_flux_budget(Float64(max_flux), Float64(tol), mode_abs_floor)

    N = N0
    K = K0
    res = nothing
    pending_low_flux_check = false
    while true
        KG_cheby = _cached_generic_cheby_sample!(cache, KG_master, N, K)
        Ysamp = if s == -2
            threaded_sampling ? GridSampling.y_sample_m2_isem_threaded(Ysol, KG_cheby) : GridSampling.y_sample_m2_isem(Ysol, KG_cheby)
        else
            threaded_sampling ? GridSampling.y_sample_p2_isem_threaded(Ysol, KG_cheby) : GridSampling.y_sample_p2_isem(Ysol, KG_cheby)
        end
        SHsamp = threaded_sampling ? GridSampling.swsh_sample_threaded(SH, KG_cheby) : GridSampling.swsh_sample(SH, KG_cheby)
        KG_trap = _cached_generic_sample!(cache, KG_master, min(N, Nmax), min(K, Kmax))
        carter_factor = _generic_carter_factor_sample(KG_trap, m, ω)
        res2 = _generic_flux_from_cheby_sample(KG_cheby, Ysamp, SHsamp, s, l, m, n, k, a, ω, ϒθ, carter_factor, cache)
        absE = abs(res2["EnergyFlux"])
        if res === nothing
            res = res2
            if mode_abs_floor > 0.0 && absE < mode_abs_floor
                break
            end
            if !confirm_low_flux && absE < low_flux_budget
                break
            end
            pending_low_flux_check = absE < low_flux_budget
        else
            relE = _relative_energy_change(res2["EnergyFlux"], res["EnergyFlux"])
            excess = abs(res2["EnergyFlux"]) / flux_scale
            factor = excess <= 1.0 ? 1.0 : min(sqrt(excess), 50.0)
            low_flux_confirmed = absE < low_flux_budget && (!confirm_low_flux || (pending_low_flux_check && _same_order_or_smaller(res2["EnergyFlux"], res["EnergyFlux"], mode_abs_floor)))
            res = res2
            if (mode_abs_floor > 0.0 && absE < mode_abs_floor) || low_flux_confirmed || factor * relE <= effective_sample_tol
                break
            end
            pending_low_flux_check = absE < low_flux_budget
        end
        (N >= Nmax && K >= Kmax) && break
        N = min(2N, Nmax)
        K = min(2K, Kmax)
    end

    if zero_low_flux && mode_abs_floor > 0.0 && abs(res["EnergyFlux"]) < mode_abs_floor
        res["Amplitude"] = 0.0 + 0.0im
        res["EnergyFlux"] = 0.0
        res["AngularMomentumFlux"] = 0.0
        res["CarterConstantFlux"] = 0.0
        res["LowFluxZeroed"] = true
    end
    res["N_sample_requested"] = N0
    res["K_sample_requested"] = K0
    res["N_sample"] = N
    res["K_sample"] = K
    res["Quadrature"] = "levin"
    return res
end

function generic_mode_flux_from_master_cached_adaptive_levin!(cache::GenericM2FluxCache, KG_master::Dict, s::Int, l::Int, m::Int, n::Int, k::Int;
        sample_tol::Float64 = DEFAULT_ADAPTIVE_LEVIN_TOL0,
        tol::Float64 = 1e-8,
        max_flux::Float64 = 1.0,
        mode_abs_floor::Float64 = 1e-16,
        zero_low_flux::Bool = false,
        threaded_sampling::Bool = false,
        confirm_low_flux::Bool = false,
        local_r_intervals::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N,
        local_theta_intervals::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N,
        min_depth::Int = DEFAULT_ADAPTIVE_LEVIN_MIN_DEPTH,
        max_depth::Int = 8,
        depth_tol_max::Float64 = DEFAULT_ADAPTIVE_LEVIN_TOL_MAX)
    if !(s == -2 || s == 2)
        error("Spin weight s must be either 2 or -2.")
    end
    local_r_intervals >= 1 || throw(ArgumentError("local_r_intervals must be positive"))
    local_theta_intervals >= 1 || throw(ArgumentError("local_theta_intervals must be positive"))
    max_depth >= 0 || throw(ArgumentError("max_depth must be nonnegative"))

    Frequencies = KG_master["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒφ = Frequencies["ϒϕ"]
    a = KG_master["a"]
    ω = (m * ϒφ + n * ϒr + k * ϒθ) / Γ
    if _skip_radiative_mode(s, a, m, ω)
        res = _zero_radiative_mode(ω, KG_master; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        res["N_sample_requested"] = local_r_intervals + 1
        effective_local_theta_intervals = _generic_adaptive_effective_theta_intervals(local_theta_intervals, k)
        res["K_sample_requested"] = effective_local_theta_intervals + 1
        res["N_sample"] = local_r_intervals + 1
        res["K_sample"] = effective_local_theta_intervals + 1
        res["Quadrature"] = "adaptive_levin_radial_theta_cc"
        return res
    end

    nr = local_r_intervals + 1
    effective_local_theta_intervals = _generic_adaptive_effective_theta_intervals(local_theta_intervals, k)
    nt = effective_local_theta_intervals + 1
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * ω; method = "jacobi")
    Ysol = _generic_isem_y_solution(s, l, m, a, ω)
    root_key = GenericAdaptiveLevin2DKey(0, 0, 0, 0)
    root = _generic_adaptive_segment_amplitude!(cache, KG_master, Ysol, SH, s, a, ω, m, n, k, root_key, nr, nt; threaded_sampling = threaded_sampling)
    root_energy = abs(_generic_adaptive_energy_from_amp(root, s, a, m, ω))
    segment_evaluations = 1
    if mode_abs_floor > 0.0 && root_energy < mode_abs_floor
        amp = zero_low_flux ? 0.0 + 0.0im : root
        carter_factor = _generic_carter_factor_sample(KG_master, m, ω)
        res = _generic_adaptive_levin_result(amp, KG_master, Ysol, SH, s, a, m, ω, ϒθ, k, carter_factor)
        res["Quadrature"] = "adaptive_levin_radial_theta_cc"
        res["AdaptiveLevin"] = true
        res["AdaptiveLevinStopReason"] = "low_flux_root"
        res["AdaptiveLevinSegmentEvaluations"] = segment_evaluations
        res["N_sample_requested"] = nr
        res["K_sample_requested"] = nt
        res["N_sample"] = nr
        res["K_sample"] = nt
        res["AdaptiveLevinMaxDepth"] = max_depth
        res["AdaptiveLevinMaxLeafDepth"] = 0
        res["AdaptiveLevinEffectiveRIntervals"] = local_r_intervals
        res["AdaptiveLevinEffectiveThetaIntervals"] = effective_local_theta_intervals
        res["AdaptiveLevinLocalThetaIntervals"] = effective_local_theta_intervals
        res["AdaptiveLevinRequestedLocalThetaIntervals"] = local_theta_intervals
        zero_low_flux && (res["LowFluxZeroed"] = true)
        return res
    end

    leaves = GenericAdaptiveLevin2DLeaf[GenericAdaptiveLevin2DLeaf(root_key, root, true)]
    previous_total = root
    split_count = 0
    stop_reason = "max_depth"
    last_relerr = Inf
    while true
        new_leaves = GenericAdaptiveLevin2DLeaf[]
        changed = false
        for leaf in leaves
            if !leaf.active || leaf.key.level_r >= max_depth
                push!(new_leaves, leaf)
                continue
            end
            child_amp = zero(ComplexF64)
            for child in _generic_adaptive_radial_children(leaf.key)
                amp = _generic_adaptive_segment_amplitude!(cache, KG_master, Ysol, SH, s, a, ω, m, n, k, child, nr, nt; threaded_sampling = threaded_sampling)
                segment_evaluations += 1
                child_amp += amp
                child_energy = abs(_generic_adaptive_energy_from_amp(amp, s, a, m, ω))
                push!(new_leaves, GenericAdaptiveLevin2DLeaf(child, amp, !(mode_abs_floor > 0.0 && child_energy < mode_abs_floor)))
            end
            changed = true
            split_count += 1
        end

        current_total = sum(leaf.amp for leaf in new_leaves; init = 0.0 + 0.0im)
        current_energy = abs(_generic_adaptive_energy_from_amp(current_total, s, a, m, ω))
        previous_energy = abs(_generic_adaptive_energy_from_amp(previous_total, s, a, m, ω))
        leaves = new_leaves
        if mode_abs_floor > 0.0 && current_energy < mode_abs_floor
            stop_reason = "low_flux"
            break
        elseif !changed
            stop_reason = "inactive_or_max_depth"
            break
        else
            scale = max(min(current_energy, previous_energy), mode_abs_floor, eps(Float64))
            last_relerr = abs(current_energy - previous_energy) / scale
            current_depth = maximum((leaf.key.level_r for leaf in leaves); init = 0)
            current_tol = _adaptive_levin_depth_tol(current_depth, sample_tol, min_depth, max_depth, depth_tol_max)
            if current_depth >= min_depth && last_relerr < current_tol
                stop_reason = "global_converged"
                break
            end
        end
        if count(leaf -> leaf.active, leaves) == 0
            stop_reason = "all_low_flux"
            break
        end
        previous_total = current_total
    end

    amp = sum(leaf.amp for leaf in leaves; init = 0.0 + 0.0im)
    if zero_low_flux && mode_abs_floor > 0.0 && abs(_generic_adaptive_energy_from_amp(amp, s, a, m, ω)) < mode_abs_floor
        amp = 0.0 + 0.0im
    end
    max_leaf_depth = maximum((leaf.key.level_r for leaf in leaves); init = 0)
    carter_factor = _generic_carter_factor_sample(KG_master, m, ω)
    res = _generic_adaptive_levin_result(amp, KG_master, Ysol, SH, s, a, m, ω, ϒθ, k, carter_factor)
    res["N_sample_requested"] = nr
    res["K_sample_requested"] = nt
    res["N_sample"] = local_r_intervals * 2^max_leaf_depth + 1
    res["K_sample"] = nt
    res["Quadrature"] = "adaptive_levin_radial_theta_cc"
    res["AdaptiveLevin"] = true
    res["AdaptiveLevinLocalRIntervals"] = local_r_intervals
    res["AdaptiveLevinLocalThetaIntervals"] = effective_local_theta_intervals
    res["AdaptiveLevinRequestedLocalThetaIntervals"] = local_theta_intervals
    res["AdaptiveLevinTol0"] = sample_tol
    res["AdaptiveLevinMinDepth"] = min_depth
    res["AdaptiveLevinAlpha"] = _adaptive_levin_depth_alpha(sample_tol, min_depth, max_depth, depth_tol_max)
    res["AdaptiveLevinAlphaMode"] = "derived"
    res["AdaptiveLevinTolMax"] = depth_tol_max
    res["AdaptiveLevinMaxDepth"] = max_depth
    res["AdaptiveLevinMaxLeafDepth"] = max_leaf_depth
    res["AdaptiveLevinEffectiveRIntervals"] = local_r_intervals * 2^max_leaf_depth
    res["AdaptiveLevinEffectiveThetaIntervals"] = effective_local_theta_intervals
    res["AdaptiveLevinAcceptedSegments"] = length(leaves)
    res["AdaptiveLevinSegmentEvaluations"] = segment_evaluations
    res["AdaptiveLevinSplitCount"] = split_count
    res["AdaptiveLevinRelErr"] = last_relerr
    res["AdaptiveLevinStopReason"] = stop_reason
    if zero_low_flux && amp == 0.0 + 0.0im
        res["LowFluxZeroed"] = true
    end
    return res
end

function _generic_adaptive_levin_result(amp::ComplexF64, KG_master::Dict, Ysol, SH, s::Int, a, m::Int, ω, ϒθ, k::Int, carter_factor::Float64)
    if s == -2
        return Dict(
            "Amplitude" => amp,
            "omega" => ω,
            "EnergyFlux" => abs2(amp) / (4.0pi * ω^2),
            "AngularMomentumFlux" => m * abs2(amp) / (4.0pi * ω^3),
            "CarterConstantFlux" => abs2(amp) * (carter_factor + k * ϒθ) / (2.0pi * ω^3),
            "Trajectory" => KG_master,
            "YSolution" => Ysol,
            "SWSH" => SH,
        )
    elseif s == 2
        hf = horizon_factor(ω, a, m)
        return Dict(
            "Amplitude" => amp,
            "omega" => ω,
            "EnergyFlux" => hf * abs2(amp),
            "AngularMomentumFlux" => hf * m * abs2(amp) / ω,
            "CarterConstantFlux" => 2 * hf * abs2(amp) * (carter_factor + k * ϒθ) / ω,
            "Trajectory" => KG_master,
            "YSolution" => Ysol,
            "SWSH" => SH,
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_generic_levin_isem(KG_master::Dict, s, l, m, n, k, N_sample::Int64, K_sample::Int64; Nmax::Int = DEFAULT_LEVIN_NMAX, Kmax::Int = DEFAULT_LEVIN_KMAX, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false, cache = nothing, confirm_low_flux::Bool = true, adaptive::Bool = true, adaptive_local_r_intervals::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, adaptive_local_theta_intervals::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, adaptive_min_depth::Int = DEFAULT_ADAPTIVE_LEVIN_MIN_DEPTH, adaptive_max_depth::Int = 8, adaptive_tol_max::Float64 = DEFAULT_ADAPTIVE_LEVIN_TOL_MAX)
    levin_cache = cache isa GenericM2FluxCache ? cache : GenericFluxCache()
    if adaptive
        return generic_mode_flux_from_master_cached_adaptive_levin!(
            levin_cache,
            KG_master,
            s,
            l,
            m,
            n,
            k;
            sample_tol = sample_tol,
            tol = tol,
            max_flux = max_flux,
            mode_abs_floor = mode_abs_floor,
            zero_low_flux = zero_low_flux,
            threaded_sampling = threaded_sampling,
            confirm_low_flux = confirm_low_flux,
            local_r_intervals = adaptive_local_r_intervals,
            local_theta_intervals = adaptive_local_theta_intervals,
            min_depth = adaptive_min_depth,
            max_depth = adaptive_max_depth,
            depth_tol_max = adaptive_tol_max,
        )
    end
    Nmax = min(Nmax, DEFAULT_LEVIN_NMAX)
    Kmax = min(Kmax, DEFAULT_LEVIN_KMAX)
    N0 = min(max(N_sample, 128), Nmax)
    K0 = min(max(K_sample, 32), Kmax)
    prewarm_generic_levin_samples!(levin_cache, KG_master; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax)
    return generic_mode_flux_from_master_cached_levin!(
        levin_cache,
        KG_master,
        s,
        l,
        m,
        n,
        k;
        N0 = N0,
        K0 = K0,
        Nmax = Nmax,
        Kmax = Kmax,
        sample_tol = sample_tol,
        tol = tol,
        max_flux = max_flux,
        mode_abs_floor = mode_abs_floor,
        zero_low_flux = zero_low_flux,
        threaded_sampling = threaded_sampling,
        confirm_low_flux = confirm_low_flux,
    )
end

function convolution_integral_generic_adaptive_levin_isem(KG_master::Dict, s, l, m, n, k, N_sample::Int64 = 128, K_sample::Int64 = 32; kwargs...)
    return convolution_integral_generic_levin_isem(KG_master, s, l, m, n, k, N_sample, K_sample; adaptive = true, kwargs...)
end

function convolution_integral_generic_levin_isem(a, p, e, x, s, l, m, n, k, N_sample::Int64, K_sample::Int64; Nmax::Int = DEFAULT_LEVIN_NMAX, Kmax::Int = DEFAULT_LEVIN_KMAX, kwargs...)
    KG = kerr_geo_orbit(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end
    KG_master = GridSampling.kerr_geo_generic_sample_dense(KG, Nmax, Kmax)
    return convolution_integral_generic_levin_isem(KG_master, s, l, m, n, k, N_sample, K_sample; Nmax = Nmax, Kmax = Kmax, kwargs...)
end

function convolution_integral_generic_adaptive_levin_isem(a, p, e, x, s, l, m, n, k, N_sample::Int64 = 128, K_sample::Int64 = 32; kwargs...)
    return convolution_integral_generic_levin_isem(a, p, e, x, s, l, m, n, k, N_sample, K_sample; adaptive = true, kwargs...)
end

function generic_mode_flux_from_master(KG_master::Dict, s::Int, l::Int, m::Int, n::Int, k::Int; N0::Int = 64, K0::Int = 16, Nmax::Int = 2^14, Kmax::Int = 2^12, sample_tol::Float64 = 1e-3, tol::Float64 = 1e-8, max_flux::Float64 = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false)
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
        excess = abs(res2["EnergyFlux"]) / flux_scale
        factor = excess <= 1.0 ? 1.0 : min(sqrt(excess), 50.0)
        low_flux_cutoff = _low_flux_cutoff(max_flux, tol)
        relE = _relative_energy_change(res2["EnergyFlux"], res["EnergyFlux"])
        low_flux_done = mode_abs_floor > 0.0 && abs(res2["EnergyFlux"]) < mode_abs_floor
        significant_for_alias = abs(res2["EnergyFlux"]) > low_flux_cutoff
        allow_stop = (N >= min_stop_N && K >= min_stop_K) || !significant_for_alias
        has_flux_reference = max_flux > 100 * eps(Float64)
        suspect_min_N = has_flux_reference && !low_flux_done ? _suspect_min_sample(min_stop_N, Nmax, excess, 1024, 2048, 4096) : min_stop_N
        suspect_min_K = has_flux_reference && !low_flux_done ? _suspect_min_sample(min_stop_K, Kmax, excess, 256, 512, 512) : min_stop_K
        suspect_needs_extra = has_flux_reference && !low_flux_done && (N < suspect_min_N || K < suspect_min_K)
        if relE == 0.0
            KG_samp = next_sample
            Ysamp = next_Ysamp
            SHsamp = next_SHsamp
            res = res2
            if allow_stop && !suspect_needs_extra
                break
            end
            continue
        end
        KG_samp = next_sample
        Ysamp = next_Ysamp
        SHsamp = next_SHsamp
        res = res2
        if allow_stop && !suspect_needs_extra && (factor * relE <= effective_sample_tol || abs(res2["EnergyFlux"]) < low_flux_cutoff || low_flux_done)
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

function generic_mode_flux(a, p, e, x, s::Int, l::Int, m::Int, n::Int, k::Int; N0::Int = 64, K0::Int = 16, Nmax::Int = 2^14, Kmax::Int = 2^12, sample_tol::Float64 = 1e-3, tol::Float64 = 1e-8, max_flux::Float64 = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false)
    KG = kerr_geo_orbit(a, p, e, x)
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

function _diff_matrix_arbitrary_nodes(x::AbstractVector{<:Real})
    n = length(x)
    D = zeros(Float64, n, n)
    w = ones(Float64, n)
    @inbounds for j in 1:n
        prod = 1.0
        xj = x[j]
        for k in 1:n
            k == j && continue
            prod *= xj - x[k]
        end
        w[j] = 1.0 / prod
    end
    @inbounds for i in 1:n
        xi = x[i]
        rowsum = 0.0
        for j in 1:n
            i == j && continue
            Dij = w[j] / (w[i] * (xi - x[j]))
            D[i, j] = Dij
            rowsum += Dij
        end
        D[i, i] = -rowsum
    end
    return D
end

function _phase_values_from_derivative(x::Vector{Float64}, gprime_vals::Vector{<:Number}, g0::Number, gπ::Number)
    n = length(x)
    gp = convert(Vector{ComplexF64}, gprime_vals)
    g = Vector{ComplexF64}(undef, n)
    g[1] = ComplexF64(g0)
    @inbounds for i in 2:n
        g[i] = g[i - 1] + 0.5 * (gp[i - 1] + gp[i]) * (x[i] - x[i - 1])
    end
    drift = g[end] - ComplexF64(gπ)
    if drift != 0
        scale = x[end] == x[1] ? 0.0 : inv(x[end] - x[1])
        @inbounds for i in 1:n
            g[i] -= drift * (x[i] - x[1]) * scale
        end
    end
    return g
end

struct Levin1DLocalSegmentPlan
    idx::UnitRange{Int}
    D::Matrix{Float64}
end

struct Levin1DLocalPlan
    perm::Vector{Int}
    x::Vector{Float64}
    segments::Vector{Levin1DLocalSegmentPlan}
end

struct Levin1DLocalSegmentFactor
    idx::UnitRange{Int}
    qadj::Matrix{ComplexF64}
    r::Matrix{ComplexF64}
    piv::Vector{Int}
    rank::Int
    n::Int
    exp_left::ComplexF64
    exp_right::ComplexF64
end

struct Levin1DLocalPhaseFactor
    perm::Vector{Int}
    segments::Vector{Levin1DLocalSegmentFactor}
end

const DEFAULT_LEVIN_LOCAL_ORDER = 8
const DEFAULT_LEVIN_NMAX = 2^12
const DEFAULT_LEVIN_KMAX = 2^9
const DEFAULT_ADAPTIVE_LEVIN_LOCAL_N = 16
const DEFAULT_ADAPTIVE_LEVIN_TOL0 = 1e-6
const DEFAULT_ADAPTIVE_LEVIN_MIN_DEPTH = 2
const DEFAULT_ADAPTIVE_LEVIN_TOL_MAX = 1e-2
const _levin_1d_local_plan_cache = Dict{Tuple{Int, Int}, Levin1DLocalPlan}()

function _levin_1d_local_plan(n::Int, local_order::Int)
    p = min(max(local_order, 4), n)
    return get!(_levin_1d_local_plan_cache, (n, p)) do
        x_ref = [cos(π * k / (n - 1)) for k in 0:n-1]
        x_phys = π .* (x_ref .+ 1.0) ./ 2.0
        perm = sortperm(x_phys)
        x = Float64.(x_phys[perm])
        segments = Levin1DLocalSegmentPlan[]
        start = 1
        while start < n
            stop = min(start + p - 1, n)
            if 0 < n - stop < 3
                stop = n
            end
            idx = start:stop
            push!(segments, Levin1DLocalSegmentPlan(idx, _diff_matrix_arbitrary_nodes(collect(x[idx]))))
            stop == n && break
            start = stop
        end
        Levin1DLocalPlan(perm, x, segments)
    end
end

function _levin_1d_local_segment(f_vals::Vector{<:Number}, gprime_vals::Vector{<:Number}, g_left::Number, g_right::Number, x::Vector{Float64}, D::Matrix{Float64} = _diff_matrix_arbitrary_nodes(x))
    n = length(f_vals)
    n <= 1 && return zero(ComplexF64)
    A = D .+ 1im .* Diagonal(convert(Vector{ComplexF64}, gprime_vals))
    rhs = convert(Vector{ComplexF64}, f_vals)
    F = qr(A, Val(true))
    R = F.R
    piv = F.p
    diagR = abs.(diag(R))
    Anorm = maximum(diagR)
    tol_sing = Anorm * 1e-12
    pvals = zeros(ComplexF64, n)
    if !all(diagR .< tol_sing)
        rank = count(≥(tol_sing), diagR)
        y = F.Q[:, 1:rank]' * rhs
        sol = R[1:rank, 1:rank] \ y
        pvals[piv[1:rank]] .= sol
    end
    return pvals[end] * exp(1im * g_right) - pvals[1] * exp(1im * g_left)
end

function _levin_1d_local_segment_batch(rhs::Matrix{ComplexF64}, gprime_vals::AbstractVector{ComplexF64}, g_left::Number, g_right::Number, D::Matrix{Float64})
    n, nrhs = size(rhs)
    n <= 1 && return zeros(ComplexF64, nrhs)
    A = D .+ 1im .* Diagonal(gprime_vals)
    F = qr(A, Val(true))
    R = F.R
    piv = F.p
    diagR = abs.(diag(R))
    Anorm = maximum(diagR)
    tol_sing = Anorm * 1e-12
    out = zeros(ComplexF64, nrhs)
    if !all(diagR .< tol_sing)
        rank = count(≥(tol_sing), diagR)
        y = F.Q[:, 1:rank]' * rhs
        sol = R[1:rank, 1:rank] \ y
        exp_left = exp(1im * g_left)
        exp_right = exp(1im * g_right)
        @inbounds for row in 1:rank
            pidx = piv[row]
            if pidx == 1
                for j in 1:nrhs
                    out[j] -= sol[row, j] * exp_left
                end
            elseif pidx == n
                for j in 1:nrhs
                    out[j] += sol[row, j] * exp_right
                end
            end
        end
    end
    return out
end

function _levin_1d_local_segment_factor_batch(rhs::Matrix{ComplexF64}, segment::Levin1DLocalSegmentFactor)
    _, nrhs = size(rhs)
    out = zeros(ComplexF64, nrhs)
    segment.rank == 0 && return out
    y = segment.qadj * rhs
    sol = segment.r \ y
    @inbounds for row in 1:segment.rank
        pidx = segment.piv[row]
        if pidx == 1
            for j in 1:nrhs
                out[j] -= sol[row, j] * segment.exp_left
            end
        elseif pidx == segment.n
            for j in 1:nrhs
                out[j] += sol[row, j] * segment.exp_right
            end
        end
    end
    return out
end

function levin_1d_integral_local(f_vals::Vector{<:Number}, gprime_vals::Vector{<:Number}, g0::Number, gπ::Number; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n = length(f_vals)
    @assert length(gprime_vals) == n "gprime length must match f ($(length(gprime_vals)) ≠ $n)"
    n <= 2 && return levin_1d_integral(f_vals, gprime_vals, g0, gπ)
    plan = _levin_1d_local_plan(n, local_order)
    x = plan.x
    f = collect(f_vals[plan.perm])
    gp = collect(gprime_vals[plan.perm])
    g = _phase_values_from_derivative(x, gp, g0, gπ)
    total = zero(ComplexF64)
    for segment in plan.segments
        idx = segment.idx
        total += _levin_1d_local_segment(collect(f[idx]), collect(gp[idx]), g[first(idx)], g[last(idx)], collect(x[idx]), segment.D)
    end
    return total
end

function _cheby_phase_nodes_pi(n::Int)
    n <= 1 && return zeros(Float64, n)
    return [π * (cos(π * k / (n - 1)) + 1.0) / 2.0 for k in 0:n-1]
end

@inline function _radial_phase_nodes(KG_samp::Dict)
    return haskey(KG_samp, "qr") ? KG_samp["qr"] : _cheby_phase_nodes_pi(length(KG_samp["Δtr"]))
end

function _radial_phase_vector(KG_samp::Dict, omega, m, n)
    qr = _radial_phase_nodes(KG_samp)
    return omega .* KG_samp["Δtr"] .- m .* KG_samp["Δφr"] .+ n .* qr
end

function _radial_phase_derivative_vector(KG_samp::Dict, omega, m, n)
    return omega .* KG_samp["dtr"] .- m .* KG_samp["dφr"] .+ n
end

function _cached_eccentric_radial_phase!(cache::EccentricFluxCache, KG_samp::Dict, m::Int, n::Int, omega)
    key = (KG_samp["N_sample"]::Int, m, n, _omega_key(omega))
    return get!(cache.phase_vectors, key) do
        (rphase = _radial_phase_vector(KG_samp, omega, m, n),
         drphase = _radial_phase_derivative_vector(KG_samp, omega, m, n))
    end
end

function _eccentric_radial_phase(KG_samp::Dict, cache, m::Int, n::Int, omega)
    return cache isa EccentricFluxCache ? _cached_eccentric_radial_phase!(cache, KG_samp, m, n, omega) :
        (rphase = _radial_phase_vector(KG_samp, omega, m, n),
         drphase = _radial_phase_derivative_vector(KG_samp, omega, m, n))
end

function _polar_phase_vector(KG_samp::Dict, omega, m, k)
    qθ = _cheby_phase_nodes_pi(length(KG_samp["Δtθ"]))
    return omega .* KG_samp["Δtθ"] .- m .* KG_samp["Δφθ"] .+ k .* qθ
end

function levin_1d_integral_local_phase(f_vals::Vector{<:Number}, gprime_vals::Vector{<:Number}, g_vals::Vector{<:Number}; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n = length(f_vals)
    @assert length(gprime_vals) == n "gprime length must match f ($(length(gprime_vals)) ≠ $n)"
    @assert length(g_vals) == n "phase length must match f ($(length(g_vals)) ≠ $n)"
    n <= 2 && return levin_1d_integral(f_vals, gprime_vals, g_vals[end], g_vals[1])
    plan = _levin_1d_local_plan(n, local_order)
    x = plan.x
    f = collect(f_vals[plan.perm])
    gp = collect(gprime_vals[plan.perm])
    g = ComplexF64.(g_vals[plan.perm])
    total = zero(ComplexF64)
    for segment in plan.segments
        idx = segment.idx
        total += _levin_1d_local_segment(collect(f[idx]), collect(gp[idx]), g[first(idx)], g[last(idx)], collect(x[idx]), segment.D)
    end
    return total
end

function _levin_1d_local_phase_factor(gprime_vals::Vector{<:Number}, g_vals::Vector{<:Number}; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n = length(gprime_vals)
    @assert length(g_vals) == n "phase length must match gprime"
    plan = _levin_1d_local_plan(n, local_order)
    gp = Vector{ComplexF64}(undef, n)
    g = Vector{ComplexF64}(undef, n)
    @inbounds for i in 1:n
        src = plan.perm[i]
        gp[i] = gprime_vals[src]
        g[i] = g_vals[src]
    end
    factors = Levin1DLocalSegmentFactor[]
    sizehint!(factors, length(plan.segments))
    for segment in plan.segments
        idx = segment.idx
        A = segment.D .+ 1im .* Diagonal(@view gp[idx])
        F = qr(A, Val(true))
        R = F.R
        piv = F.p
        diagR = abs.(diag(R))
        Anorm = maximum(diagR)
        tol_sing = Anorm * 1e-12
        rank = all(diagR .< tol_sing) ? 0 : count(≥(tol_sing), diagR)
        qadj = rank == 0 ? zeros(ComplexF64, 0, length(idx)) : Matrix(F.Q[:, 1:rank]')
        r = rank == 0 ? zeros(ComplexF64, 0, 0) : Matrix(R[1:rank, 1:rank])
        push!(factors, Levin1DLocalSegmentFactor(idx, qadj, r, collect(piv), rank, length(idx), exp(1im * g[first(idx)]), exp(1im * g[last(idx)])))
    end
    return Levin1DLocalPhaseFactor(plan.perm, factors)
end

function _levin_1d_local_phase_factor_from_values(g_vals::Vector{<:Number}; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n = length(g_vals)
    plan = _levin_1d_local_plan(n, local_order)
    g = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        g[i] = real(g_vals[plan.perm[i]])
    end
    factors = Levin1DLocalSegmentFactor[]
    sizehint!(factors, length(plan.segments))
    for segment in plan.segments
        idx = segment.idx
        gp = segment.D * collect(@view g[idx])
        A = segment.D .+ 1im .* Diagonal(ComplexF64.(gp))
        F = qr(A, Val(true))
        R = F.R
        piv = F.p
        diagR = abs.(diag(R))
        Anorm = maximum(diagR)
        tol_sing = Anorm * 1e-12
        rank = all(diagR .< tol_sing) ? 0 : count(≥(tol_sing), diagR)
        qadj = rank == 0 ? zeros(ComplexF64, 0, length(idx)) : Matrix(F.Q[:, 1:rank]')
        r = rank == 0 ? zeros(ComplexF64, 0, 0) : Matrix(R[1:rank, 1:rank])
        push!(factors, Levin1DLocalSegmentFactor(idx, qadj, r, collect(piv), rank, length(idx), exp(1im * g[first(idx)]), exp(1im * g[last(idx)])))
    end
    return Levin1DLocalPhaseFactor(plan.perm, factors)
end

function _conjugate_levin_1d_local_phase_factor(factor::Levin1DLocalPhaseFactor)
    factors = Levin1DLocalSegmentFactor[]
    sizehint!(factors, length(factor.segments))
    for segment in factor.segments
        push!(factors, Levin1DLocalSegmentFactor(
            segment.idx,
            conj.(segment.qadj),
            conj.(segment.r),
            segment.piv,
            segment.rank,
            segment.n,
            conj(segment.exp_left),
            conj(segment.exp_right),
        ))
    end
    return Levin1DLocalPhaseFactor(factor.perm, factors)
end

function _levin_1d_integral_local_phase_factored(f_vals::Vector{<:Number}, factor::Levin1DLocalPhaseFactor)
    n = length(factor.perm)
    f = Vector{ComplexF64}(undef, n)
    @inbounds for i in 1:n
        f[i] = f_vals[factor.perm[i]]
    end
    total = zero(ComplexF64)
    for segment in factor.segments
        segment.rank == 0 && continue
        y = segment.qadj * @view f[segment.idx]
        sol = segment.r \ y
        @inbounds for row in 1:segment.rank
            pidx = segment.piv[row]
            if pidx == 1
                total -= sol[row] * segment.exp_left
            elseif pidx == segment.n
                total += sol[row] * segment.exp_right
            end
        end
    end
    return total
end

function levin_1d_integral_local_phase_values(f_vals::Vector{<:Number}, g_vals::Vector{<:Number}; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    factor = _levin_1d_local_phase_factor_from_values(g_vals; local_order = local_order)
    return _levin_1d_integral_local_phase_factored(f_vals, factor)
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

function levin_2d_integral_local(f_vals::Matrix{<:Number}, g1_prime_vals::Vector{<:Number}, g2_prime_vals::Vector{<:Number}, g1_0::Number, g1_π::Number, g2_0::Number, g2_π::Number; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n, k = size(f_vals)
    @assert length(g1_prime_vals) == n "g1_prime_vals length must match x-dimension n ($(length(g1_prime_vals)) ≠ $n)"
    @assert length(g2_prime_vals) == k "g2_prime_vals length must match y-dimension k ($(length(g2_prime_vals)) ≠ $k)"
    h_vals = Vector{ComplexF64}(undef, n)
    @inbounds for i in 1:n
        h_vals[i] = levin_1d_integral_local(collect(@view f_vals[i, :]), collect(g2_prime_vals), g2_0, g2_π; local_order = local_order)
    end
    return levin_1d_integral_local(h_vals, collect(g1_prime_vals), g1_0, g1_π; local_order = local_order)
end

function levin_2d_integral_local_phase(f_vals::Matrix{<:Number}, g1_prime_vals::Vector{<:Number}, g2_prime_vals::Vector{<:Number}, g1_vals::Vector{<:Number}, g2_vals::Vector{<:Number}; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n, k = size(f_vals)
    @assert length(g1_prime_vals) == n "g1_prime_vals length must match x-dimension n ($(length(g1_prime_vals)) ≠ $n)"
    @assert length(g2_prime_vals) == k "g2_prime_vals length must match y-dimension k ($(length(g2_prime_vals)) ≠ $k)"
    @assert length(g1_vals) == n "g1_vals length must match x-dimension n ($(length(g1_vals)) ≠ $n)"
    @assert length(g2_vals) == k "g2_vals length must match y-dimension k ($(length(g2_vals)) ≠ $k)"
    if n >= 8 && k >= 8
        return levin_2d_integral_local_phase_batched(f_vals, g1_prime_vals, g2_prime_vals, g1_vals, g2_vals; local_order = local_order)
    end
    h_vals = Vector{ComplexF64}(undef, n)
    @inbounds for i in 1:n
        h_vals[i] = levin_1d_integral_local_phase(collect(@view f_vals[i, :]), collect(g2_prime_vals), collect(g2_vals); local_order = local_order)
    end
    return levin_1d_integral_local_phase(h_vals, collect(g1_prime_vals), collect(g1_vals); local_order = local_order)
end

function levin_2d_integral_local_phase_batched(f_vals::Matrix{<:Number}, g1_prime_vals::Vector{<:Number}, g2_prime_vals::Vector{<:Number}, g1_vals::Vector{<:Number}, g2_vals::Vector{<:Number}; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n, k = size(f_vals)
    plan2 = _levin_1d_local_plan(k, local_order)
    gp2 = Vector{ComplexF64}(undef, k)
    g2 = Vector{ComplexF64}(undef, k)
    @inbounds for i in 1:k
        src = plan2.perm[i]
        gp2[i] = g2_prime_vals[src]
        g2[i] = g2_vals[src]
    end
    h_vals = zeros(ComplexF64, n)
    for segment in plan2.segments
        idx = segment.idx
        rhs = Matrix{ComplexF64}(undef, length(idx), n)
        @inbounds for (ii, jj) in enumerate(idx)
            srcj = plan2.perm[jj]
            for row in 1:n
                rhs[ii, row] = f_vals[row, srcj]
            end
        end
        h_vals .+= _levin_1d_local_segment_batch(rhs, @view(gp2[idx]), g2[first(idx)], g2[last(idx)], segment.D)
    end
    return levin_1d_integral_local_phase(h_vals, collect(g1_prime_vals), collect(g1_vals); local_order = local_order)
end

function levin_2d_integral_local_phase_pair_batched(f1::Matrix{<:Number}, f2::Matrix{<:Number}, g1_prime_vals1::Vector{<:Number}, g1_prime_vals2::Vector{<:Number}, g2_prime_vals::Vector{<:Number}, g1_vals1::Vector{<:Number}, g1_vals2::Vector{<:Number}, g2_vals::Vector{<:Number}; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n, k = size(f1)
    size(f2) == (n, k) || throw(ArgumentError("paired 2D Levin inputs must have the same size"))
    plan2 = _levin_1d_local_plan(k, local_order)
    gp2 = Vector{ComplexF64}(undef, k)
    g2 = Vector{ComplexF64}(undef, k)
    @inbounds for i in 1:k
        src = plan2.perm[i]
        gp2[i] = g2_prime_vals[src]
        g2[i] = g2_vals[src]
    end
    h1 = zeros(ComplexF64, n)
    h2 = zeros(ComplexF64, n)
    for segment in plan2.segments
        idx = segment.idx
        rhs = Matrix{ComplexF64}(undef, length(idx), 2n)
        @inbounds for (ii, jj) in enumerate(idx)
            srcj = plan2.perm[jj]
            for row in 1:n
                rhs[ii, row] = f1[row, srcj]
                rhs[ii, n + row] = f2[row, srcj]
            end
        end
        vals = _levin_1d_local_segment_batch(rhs, @view(gp2[idx]), g2[first(idx)], g2[last(idx)], segment.D)
        @inbounds for row in 1:n
            h1[row] += vals[row]
            h2[row] += vals[n + row]
        end
    end
    return (
        levin_1d_integral_local_phase(h1, collect(g1_prime_vals1), collect(g1_vals1); local_order = local_order),
        levin_1d_integral_local_phase(h2, collect(g1_prime_vals2), collect(g1_vals2); local_order = local_order),
    )
end

function levin_2d_integral_local_phase_pair_batched_radial_values(f1::Matrix{<:Number}, f2::Matrix{<:Number}, g2_prime_vals::Vector{<:Number}, g1_vals1::Vector{<:Number}, g1_vals2::Vector{<:Number}, g2_vals::Vector{<:Number}; local_order::Int = DEFAULT_LEVIN_LOCAL_ORDER)
    n, k = size(f1)
    size(f2) == (n, k) || throw(ArgumentError("paired 2D Levin inputs must have the same size"))
    plan2 = _levin_1d_local_plan(k, local_order)
    gp2 = Vector{ComplexF64}(undef, k)
    g2 = Vector{ComplexF64}(undef, k)
    @inbounds for i in 1:k
        src = plan2.perm[i]
        gp2[i] = g2_prime_vals[src]
        g2[i] = g2_vals[src]
    end
    h1 = zeros(ComplexF64, n)
    h2 = zeros(ComplexF64, n)
    for segment in plan2.segments
        idx = segment.idx
        rhs = Matrix{ComplexF64}(undef, length(idx), 2n)
        @inbounds for (ii, jj) in enumerate(idx)
            srcj = plan2.perm[jj]
            for row in 1:n
                rhs[ii, row] = f1[row, srcj]
                rhs[ii, n + row] = f2[row, srcj]
            end
        end
        vals = _levin_1d_local_segment_batch(rhs, @view(gp2[idx]), g2[first(idx)], g2[last(idx)], segment.D)
        @inbounds for row in 1:n
            h1[row] += vals[row]
            h2[row] += vals[n + row]
        end
    end
    return (
        levin_1d_integral_local_phase_values(h1, collect(g1_vals1); local_order = local_order),
        levin_1d_integral_local_phase_values(h2, collect(g1_vals2); local_order = local_order),
    )
end

function levin_2d_integral_local_phase_pair_batched_radial_values_factored(f1::Matrix{<:Number}, f2::Matrix{<:Number}, theta_factor::Levin1DLocalPhaseFactor, radial_factor1::Levin1DLocalPhaseFactor, radial_factor2::Levin1DLocalPhaseFactor)
    n, k = size(f1)
    size(f2) == (n, k) || throw(ArgumentError("paired 2D Levin inputs must have the same size"))
    length(theta_factor.perm) == k || throw(ArgumentError("theta factor length must match y-dimension"))
    length(radial_factor1.perm) == n || throw(ArgumentError("first radial factor length must match x-dimension"))
    length(radial_factor2.perm) == n || throw(ArgumentError("second radial factor length must match x-dimension"))
    h1 = zeros(ComplexF64, n)
    h2 = zeros(ComplexF64, n)
    for segment in theta_factor.segments
        idx = segment.idx
        rhs = Matrix{ComplexF64}(undef, length(idx), 2n)
        @inbounds for (ii, jj) in enumerate(idx)
            srcj = theta_factor.perm[jj]
            for row in 1:n
                rhs[ii, row] = f1[row, srcj]
                rhs[ii, n + row] = f2[row, srcj]
            end
        end
        vals = _levin_1d_local_segment_factor_batch(rhs, segment)
        @inbounds for row in 1:n
            h1[row] += vals[row]
            h2[row] += vals[n + row]
        end
    end
    return (
        _levin_1d_integral_local_phase_factored(h1, radial_factor1),
        _levin_1d_integral_local_phase_factored(h2, radial_factor2),
    )
end

function convolution_integral_generic_trapezoidal(a, p, e, x, s, l, m, n, k, N_sample, K_sample)

    KG = kerr_geo_orbit(a, p, e, x)

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

    KG = kerr_geo_orbit(a, p, e, x)

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
    KG = kerr_geo_orbit(a, p, e, 1.0)  # Assume this handles equatorial case
    
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

    KG = kerr_geo_orbit(a, p, e, 1.0)

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
    KG = kerr_geo_orbit(a, p, 0.0, x)  # Assume this handles radial-fixed case
    
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

    KG = kerr_geo_orbit(a, p, 0.0, x)

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

    KG = kerr_geo_orbit(a, p, 0.0, 1.0)
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

    KG = kerr_geo_orbit(a, p, 0.0, 1.0)
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
    KG = kerr_geo_orbit(a, p, e, x)
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
    KG = kerr_geo_orbit(a, p, e, x)
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
    KG = kerr_geo_orbit(a, p, 0.0, 1.0)
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
        integral = 4im * π * ω * J / (GridSampling._isem_gsn_incidence_amplitude(Y_soln) * Γ)
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
        integral = - im * π * factor * J * η / (κ * GridSampling._isem_gsn_incidence_amplitude(Y_soln) * Γ)
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
    KG = kerr_geo_orbit(a, p, e, 1.0)
    KG_sample = kerr_geo_eccentric_sample_dense(KG, Nmax)
    return convolution_integral_eccentric_trapezoidal_isem(KG_sample, s, l, m, n, N_sample; Nmax = Nmax, kwargs...)
end

function convolution_integral_eccentric_levin_isem(KG_sample::Dict, s, l, m, n, N_sample::Int64; Nmax::Int = DEFAULT_LEVIN_NMAX, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false, cache = nothing)
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
        res["Quadrature"] = "levin"
        return res
    end
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * ω; method = "jacobi")
    Ysol = _isem_y_solution(s, l, m, a, ω)
    levin_cache = cache isa EccentricFluxCache ? cache : EccentricFluxCache()
    flux_scale = max(Float64(max_flux), eps(Float64))
    effective_sample_tol = min(sample_tol, 3.0 * sqrt(Float64(tol)))
    low_flux_budget = _levin_low_flux_budget(Float64(max_flux), Float64(tol), mode_abs_floor)

    N = N_sample
    res = nothing
    pending_low_flux_check = false
    while true
        KG_cheby = _cached_eccentric_cheby_sample!(levin_cache, KG_sample, N)
        Ysamp = if s == 2
            threaded_sampling ? GridSampling.y_sample_p2_isem_threaded(Ysol, KG_cheby) : y_sample_p2_isem(Ysol, KG_cheby)
        else
            threaded_sampling ? GridSampling.y_sample_m2_isem_threaded(Ysol, KG_cheby) : y_sample_m2_isem(Ysol, KG_cheby)
        end
        res2 = _eccentric_flux_from_cheby_sample(KG_cheby, Ysamp, SH, s, a, ω, m, n; cache = levin_cache)
        absE = abs(res2["EnergyFlux"])
        if res === nothing
            res = res2
            if mode_abs_floor > 0.0 && absE < mode_abs_floor
                break
            end
            pending_low_flux_check = absE < low_flux_budget
        else
            relE = _relative_energy_change(res2["EnergyFlux"], res["EnergyFlux"])
            excess = abs(res2["EnergyFlux"]) / flux_scale
            factor = excess <= 1.0 ? 1.0 : min(sqrt(excess), 50.0)
            low_flux_confirmed = pending_low_flux_check && absE < low_flux_budget && _levin_low_flux_stable(res2["EnergyFlux"], res["EnergyFlux"], mode_abs_floor)
            res = res2
            if (mode_abs_floor > 0.0 && absE < mode_abs_floor) || low_flux_confirmed || factor * relE <= effective_sample_tol
                break
            end
            pending_low_flux_check = absE < low_flux_budget
        end
        N >= Nmax && break
        N = min(2N, Nmax)
    end

    if zero_low_flux && mode_abs_floor > 0.0 && abs(res["EnergyFlux"]) < mode_abs_floor
        res["Amplitude"] = 0.0 + 0.0im
        res["EnergyFlux"] = 0.0
        res["AngularMomentumFlux"] = 0.0
        res["CarterConstantFlux"] = 0.0
        res["LowFluxZeroed"] = true
    end
    res["N_sample_requested"] = N_sample
    res["N_sample"] = N
    res["Quadrature"] = "levin"
    return res
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

@inline function _omega_key(ω)
    return reinterpret(UInt64, Float64(ω))
end

function _cached_eccentric_levin_factor!(cache::EccentricFluxCache, N::Int, m::Int, n::Int, ω, sign::Int, gprime_vals::Vector{<:Number}, g_vals::Vector{<:Number})
    key_pos = (N, m, n, 1, _omega_key(ω))
    positive_factor = get!(cache.levin_phase_factors, key_pos) do
        _levin_1d_local_phase_factor(gprime_vals, g_vals)
    end
    sign == 1 && return positive_factor

    key_neg = (N, m, n, -1, _omega_key(ω))
    return get!(cache.levin_phase_factors, key_neg) do
        _conjugate_levin_1d_local_phase_factor(positive_factor)
    end
end

function _eccentric_levin_integral(f_vals::Vector{<:Number}, gprime_vals::Vector{<:Number}, g_vals::Vector{<:Number}, cache, N::Int, m::Int, n::Int, ω, sign::Int)
    if cache isa EccentricFluxCache
        factor = _cached_eccentric_levin_factor!(cache, N, m, n, ω, sign, gprime_vals, g_vals)
        return _levin_1d_integral_local_phase_factored(f_vals, factor)
    end
    return sign == 1 ? levin_1d_integral_local_phase(f_vals, gprime_vals, g_vals) : levin_1d_integral_local_phase(f_vals, -gprime_vals, -g_vals)
end

function _internal_radial_phase(KG_samp::Dict, s::Int, a, ω, m)
    r = KG_samp["r"]
    rs = KG_samp["rs"]
    κ = sqrt(1.0 - a^2)
    rp = 1.0 + κ
    rm = 1.0 - κ
    log_term = log.((r .- rp) ./ (r .- rm))
    ψ = s == -2 ? ω .* rs .- a .* m .* log_term ./ (2.0κ) :
                  .-ω .* rs .+ a .* m .* log_term ./ (2.0κ)
    return (phase_values = ψ, phase_factor = exp.(1im .* ψ))
end

function _cached_eccentric_factored_levin_factor!(cache::EccentricFluxCache, N::Int, m::Int, n::Int, ω, sign::Int, phase_values::Vector{<:Number})
    key = (N, m, n, sign, _omega_key(ω))
    return get!(cache.levin_factored_phase_factors, key) do
        _levin_1d_local_phase_factor_from_values(phase_values)
    end
end

function _eccentric_levin_integral_factored(f_vals::Vector{<:Number}, phase_factor::Vector{<:Number}, phase_values::Vector{<:Number}, cache, N::Int, m::Int, n::Int, ω, sign::Int)
    smooth_vals = f_vals ./ phase_factor
    if cache isa EccentricFluxCache
        factor = _cached_eccentric_factored_levin_factor!(cache, N, m, n, ω, sign, phase_values)
        return _levin_1d_integral_local_phase_factored(smooth_vals, factor)
    end
    factor = _levin_1d_local_phase_factor_from_values(phase_values)
    return _levin_1d_integral_local_phase_factored(smooth_vals, factor)
end

function _eccentric_flux_from_cheby_sample(KG_samp::Dict, Ysamp::Dict, SH, s, a, ω, m, n; cache = nothing)
    N = KG_samp["N_sample"]::Int
    if s == 2
        S0 = SH(π/2, 0.0)
        S1 = SH(π/2, 0.0; theta_derivative = 1)
        S2 = (m^2 + (a * ω)^2 - 2 - 2 * a * ω * m - SH.lambda) * S0
        Jp, Jm, _, _, _, prefactor = integrand_eccentric_sample_cheby_p2(KG_samp, Ysamp, (S0, S1, S2), n)
        phase = _eccentric_radial_phase(KG_samp, cache, m, n, ω)
        rphase = phase.rphase
        internal_phase = _internal_radial_phase(KG_samp, s, a, ω, m)
        integral = (_eccentric_levin_integral_factored(Jp, internal_phase.phase_factor, rphase .+ internal_phase.phase_values, cache, N, m, n, ω, 1) +
                    _eccentric_levin_integral_factored(Jm, internal_phase.phase_factor, .-rphase .+ internal_phase.phase_values, cache, N, m, n, ω, -1)) * prefactor
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
        Jp, Jm, _, _, _, prefactor = integrand_eccentric_sample_cheby_m2(KG_samp, Ysamp, (S0, S1, S2), n)
        phase = _eccentric_radial_phase(KG_samp, cache, m, n, ω)
        rphase = phase.rphase
        internal_phase = _internal_radial_phase(KG_samp, s, a, ω, m)
        integral = (_eccentric_levin_integral_factored(Jp, internal_phase.phase_factor, rphase .+ internal_phase.phase_values, cache, N, m, n, ω, 1) +
                    _eccentric_levin_integral_factored(Jm, internal_phase.phase_factor, .-rphase .+ internal_phase.phase_values, cache, N, m, n, ω, -1)) * prefactor
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

struct AdaptiveEccentricLevinLeaf
    level::Int
    bin::Int
    amp::ComplexF64
    active::Bool
end

@inline function _adaptive_levin_depth_alpha(tol0::Float64, min_depth::Int, max_depth::Int, tol_max::Float64)
    if !(tol0 > 0.0) || !(tol_max > tol0) || max_depth <= min_depth
        return 0.0
    end
    return log2(tol_max / tol0) / max(1, max_depth - min_depth)
end

@inline function _adaptive_levin_depth_tol(depth::Int, tol0::Float64, min_depth::Int, max_depth::Int, tol_max::Float64)
    alpha = _adaptive_levin_depth_alpha(tol0, min_depth, max_depth, tol_max)
    d = max(depth - min_depth, 0)
    return min(tol_max, max(tol0, tol0 * 2.0^(alpha * d)))
end

function _cached_adaptive_eccentric_segment_bundle!(cache, KG::Dict, local_n::Int, level::Int, bin::Int)
    if cache isa EccentricFluxCache
        return get!(cache.adaptive_levin_segments, (local_n, level, bin)) do
            nbin = 2^level
            qlo = π * bin / nbin
            qhi = π * (bin + 1) / nbin
            eccentric_segment_sample_bundle_cheby(KG, qlo, qhi, local_n)
        end
    end
    nbin = 2^level
    qlo = π * bin / nbin
    qhi = π * (bin + 1) / nbin
    return eccentric_segment_sample_bundle_cheby(KG, qlo, qhi, local_n)
end

function prewarm_eccentric_adaptive_levin_segments!(cache::EccentricFluxCache, KG_master::Dict;
        local_n::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N,
        max_depth::Int = 8)
    KG = haskey(KG_master, "Energy") ? KG_master : _kg_from_presampled_master(KG_master)
    max_depth >= 0 || throw(ArgumentError("max_depth must be nonnegative"))
    depth_limit = max_depth
    @inbounds for level in 0:depth_limit
        for bin in 0:(2^level - 1)
            _cached_adaptive_eccentric_segment_bundle!(cache, KG, local_n, level, bin)
        end
    end
    return cache
end

function _adaptive_eccentric_factored_segment(bundle::Dict, Ysamp::Dict, SH, s::Int, a, ω, m::Int, n::Int)
    KG_samp = bundle["sample"]
    if s == -2
        S0 = SH(π/2, 0.0)
        S1 = SH(π/2, 0.0; theta_derivative = 1)
        S2 = (m^2 + (a * ω)^2 + 2 - 2 * a * ω * m - SH.lambda) * S0
        return integrand_eccentric_factored_segment_cheby_m2(KG_samp, bundle["geometry"], bundle["phase_basis"], Ysamp, (S0, S1, S2), m, n, a, ω)
    elseif s == 2
        S0 = SH(π/2, 0.0)
        S1 = SH(π/2, 0.0; theta_derivative = 1)
        S2 = (m^2 + (a * ω)^2 - 2 - 2 * a * ω * m - SH.lambda) * S0
        Jp, Jm, _, _, _, prefactor = integrand_eccentric_sample_cheby_p2(KG_samp, Ysamp, (S0, S1, S2), n)
        phase = _eccentric_radial_phase(KG_samp, nothing, m, n, ω)
        internal = _internal_radial_phase(KG_samp, s, a, ω, m)
        return (
            smooth_p = ComplexF64.(Jp ./ internal.phase_factor),
            smooth_m = ComplexF64.(Jm ./ internal.phase_factor),
            phase_p = ComplexF64.(phase.rphase .+ internal.phase_values),
            phase_m = ComplexF64.(.-phase.rphase .+ internal.phase_values),
            prefactor = ComplexF64(prefactor),
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function _adaptive_eccentric_segment_amplitude!(cache, KG::Dict, Ysol, SH, s::Int, a, ω, m::Int, n::Int, local_n::Int, level::Int, bin::Int; threaded_sampling::Bool = false)
    bundle = _cached_adaptive_eccentric_segment_bundle!(cache, KG, local_n, level, bin)
    KG_samp = bundle["sample"]
    Ysamp = if s == 2
        threaded_sampling ? GridSampling.y_sample_p2_isem_threaded(Ysol, KG_samp) : y_sample_p2_isem(Ysol, KG_samp)
    else
        threaded_sampling ? GridSampling.y_sample_m2_isem_threaded(Ysol, KG_samp) : y_sample_m2_isem(Ysol, KG_samp)
    end
    factored = _adaptive_eccentric_factored_segment(bundle, Ysamp, SH, s, a, ω, m, n)
    jac = (π / 2^level) / π
    factor_p = _levin_1d_local_phase_factor_from_values(factored.phase_p; local_order = local_n)
    factor_m = _levin_1d_local_phase_factor_from_values(factored.phase_m; local_order = local_n)
    return ComplexF64((
        _levin_1d_integral_local_phase_factored(jac .* factored.smooth_p, factor_p) +
        _levin_1d_integral_local_phase_factored(jac .* factored.smooth_m, factor_m)
    ) * factored.prefactor)
end

function _adaptive_eccentric_levin_result(amp::ComplexF64, KG_sample::Dict, Ysol, SH, s::Int, a, ω, m::Int)
    if s == 2
        hf = horizon_factor(ω, a, m)
        return Dict(
            "Amplitude" => amp,
            "omega" => ω,
            "EnergyFlux" => hf * abs2(amp),
            "AngularMomentumFlux" => hf * m * abs2(amp) / ω,
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG_sample,
            "YSolution" => Ysol,
            "SWSH" => SH,
        )
    elseif s == -2
        return Dict(
            "Amplitude" => amp,
            "omega" => ω,
            "EnergyFlux" => abs2(amp) / (4.0pi * ω^2),
            "AngularMomentumFlux" => m * abs2(amp) / (4.0pi * ω^3),
            "CarterConstantFlux" => 0.0,
            "Trajectory" => KG_sample,
            "YSolution" => Ysol,
            "SWSH" => SH,
        )
    else
        error("Spin weight s must be either 2 or -2.")
    end
end

function convolution_integral_eccentric_adaptive_levin_isem(KG_sample::Dict, s, l, m, n;
        tol = 1e-8,
        tol0::Float64 = DEFAULT_ADAPTIVE_LEVIN_TOL0,
        sample_tol = nothing,
        local_n::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N,
        min_depth::Int = DEFAULT_ADAPTIVE_LEVIN_MIN_DEPTH,
        depth_tol_max::Float64 = DEFAULT_ADAPTIVE_LEVIN_TOL_MAX,
        max_depth::Int = 8,
        mode_abs_floor::Float64 = 1e-16,
        zero_low_flux::Bool = false,
        threaded_sampling::Bool = false,
        cache = nothing)
    local_n >= 2 || throw(ArgumentError("local_n must be at least 2"))
    max_depth >= 0 || throw(ArgumentError("max_depth must be nonnegative"))
    adaptive_tol0 = sample_tol === nothing ? tol0 : Float64(sample_tol)
    Frequencies = KG_sample["Frequencies"]
    Γ = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒφ = Frequencies["ϒϕ"]
    a = KG_sample["a"]
    ω = (m * ϒφ + n * ϒr) / Γ
    if _skip_radiative_mode(s, a, m, ω)
        res = _zero_radiative_mode(ω, KG_sample; reason = s == -2 ? "infinity_static_frequency" : "horizon_static_frequency")
        res["N_sample_requested"] = local_n
        res["N_sample"] = local_n
        res["Quadrature"] = "adaptive_levin"
        res["AdaptiveLevin"] = true
        res["AdaptiveLevinMaxDepth"] = max_depth
        res["AdaptiveLevinMaxLeafDepth"] = 0
        res["AdaptiveLevinEffectiveIntervals"] = local_n
        return res
    end

    KG = haskey(KG_sample, "Energy") ? KG_sample : _kg_from_presampled_master(KG_sample)
    levin_cache = cache isa EccentricFluxCache ? cache : EccentricFluxCache()
    SH = spin_weighted_spheroidal_harmonic(s, l, m, a * ω; method = "jacobi")
    Ysol = _isem_y_solution(s, l, m, a, ω)

    segment_evaluations = 0
    root = _adaptive_eccentric_segment_amplitude!(levin_cache, KG, Ysol, SH, s, a, ω, m, n, local_n, 0, 0; threaded_sampling = threaded_sampling)
    segment_evaluations += 1
    root_energy = abs(_adaptive_eccentric_levin_result(root, KG_sample, Ysol, SH, s, a, ω, m)["EnergyFlux"])
    if mode_abs_floor > 0.0 && root_energy < mode_abs_floor
        res = _adaptive_eccentric_levin_result(root, KG_sample, Ysol, SH, s, a, ω, m)
        res["N_sample_requested"] = local_n
        res["N_sample"] = local_n
        res["Quadrature"] = "adaptive_levin"
        res["AdaptiveLevin"] = true
        res["AdaptiveLevinSegmentEvaluations"] = segment_evaluations
        res["AdaptiveLevinAcceptedSegments"] = 1
        res["AdaptiveLevinMaxDepth"] = max_depth
        res["AdaptiveLevinMaxLeafDepth"] = 0
        res["AdaptiveLevinEffectiveIntervals"] = local_n
        return res
    end

    leaves = AdaptiveEccentricLevinLeaf[AdaptiveEccentricLevinLeaf(0, 0, root, true)]
    previous_total = root
    split_count = 0
    stop_reason = "max_depth"
    while true
        new_leaves = AdaptiveEccentricLevinLeaf[]
        changed = false
        for leaf in leaves
            if !leaf.active || leaf.level >= max_depth
                push!(new_leaves, leaf)
                continue
            end
            child_level = leaf.level + 1
            left_bin = 2 * leaf.bin
            right_bin = left_bin + 1
            left = _adaptive_eccentric_segment_amplitude!(levin_cache, KG, Ysol, SH, s, a, ω, m, n, local_n, child_level, left_bin; threaded_sampling = threaded_sampling)
            right = _adaptive_eccentric_segment_amplitude!(levin_cache, KG, Ysol, SH, s, a, ω, m, n, local_n, child_level, right_bin; threaded_sampling = threaded_sampling)
            segment_evaluations += 2
            child_amp = left + right
            child_energy = abs(_adaptive_eccentric_levin_result(child_amp, KG_sample, Ysol, SH, s, a, ω, m)["EnergyFlux"])
            low_flux = mode_abs_floor > 0.0 && child_energy < mode_abs_floor
            push!(new_leaves, AdaptiveEccentricLevinLeaf(child_level, left_bin, left, !low_flux))
            push!(new_leaves, AdaptiveEccentricLevinLeaf(child_level, right_bin, right, !low_flux))
            changed = true
            split_count += 1
        end

        current_total = sum(leaf.amp for leaf in new_leaves; init = 0.0 + 0.0im)
        current_energy = abs(_adaptive_eccentric_levin_result(current_total, KG_sample, Ysol, SH, s, a, ω, m)["EnergyFlux"])
        previous_energy = abs(_adaptive_eccentric_levin_result(previous_total, KG_sample, Ysol, SH, s, a, ω, m)["EnergyFlux"])
        leaves = new_leaves

        if mode_abs_floor > 0.0 && current_energy < mode_abs_floor
            stop_reason = "low_flux"
            break
        elseif !changed
            stop_reason = "inactive_or_max_depth"
            break
        else
            scale = max(min(current_energy, previous_energy), mode_abs_floor, eps(Float64))
            global_relerr = abs(current_energy - previous_energy) / scale
            current_depth = maximum((leaf.level for leaf in leaves); init = 0)
            current_tol = _adaptive_levin_depth_tol(current_depth, adaptive_tol0, min_depth, max_depth, depth_tol_max)
            if current_depth >= min_depth && global_relerr < current_tol
                stop_reason = "global_converged"
                break
            end
        end

        if count(leaf -> leaf.active, leaves) == 0
            stop_reason = "all_low_flux"
            break
        end
        previous_total = current_total
    end

    amp = sum(leaf.amp for leaf in leaves; init = 0.0 + 0.0im)
    res = _adaptive_eccentric_levin_result(amp, KG_sample, Ysol, SH, s, a, ω, m)
    if zero_low_flux && mode_abs_floor > 0.0 && abs(res["EnergyFlux"]) < mode_abs_floor
        res["Amplitude"] = 0.0 + 0.0im
        res["EnergyFlux"] = 0.0
        res["AngularMomentumFlux"] = 0.0
        res["CarterConstantFlux"] = 0.0
        res["LowFluxZeroed"] = true
    end
    res["N_sample_requested"] = local_n
    max_leaf_depth = maximum((leaf.level for leaf in leaves); init = 0)
    res["N_sample"] = local_n * 2^max_leaf_depth
    res["Quadrature"] = "adaptive_levin"
    res["AdaptiveLevin"] = true
    res["AdaptiveLevinLocalN"] = local_n
    res["AdaptiveLevinTol0"] = adaptive_tol0
    res["AdaptiveLevinMinDepth"] = min_depth
    res["AdaptiveLevinAlpha"] = _adaptive_levin_depth_alpha(adaptive_tol0, min_depth, max_depth, depth_tol_max)
    res["AdaptiveLevinAlphaMode"] = "derived"
    res["AdaptiveLevinTolMax"] = depth_tol_max
    res["AdaptiveLevinMaxDepth"] = max_depth
    res["AdaptiveLevinMaxLeafDepth"] = max_leaf_depth
    res["AdaptiveLevinEffectiveIntervals"] = local_n * 2^max_leaf_depth
    res["AdaptiveLevinAcceptedSegments"] = length(leaves)
    res["AdaptiveLevinSegmentEvaluations"] = segment_evaluations
    res["AdaptiveLevinSplitCount"] = split_count
    res["AdaptiveLevinStopReason"] = stop_reason
    return res
end

function convolution_integral_eccentric_trapezoidal_isem(KG_sample::Dict, s, l, m, n, N_sample::Int64; Nmax::Int = 2^14, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false, cache = nothing)
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
    min_stop_N = _eccentric_fourier_min_sample(min(max(4N_sample, 256), Nmax), Nmax, n)

    while N < Nmax
        N2 = min(2N, Nmax)
        KG_samp2 = cache isa EccentricFluxCache ? _cached_eccentric_sample!(cache, KG_sample, N2) : subsample_eccentric_sample(KG_sample, N2)
        Ysamp2 = threaded_sampling ? GridSampling.refine_generic_y_sample_threaded(Ysol, Ysamp, KG_samp2) : refine_eccentric_y_sample(Ysol, Ysamp, KG_sample, N)
        res2 = _eccentric_flux_from_sample(KG_samp2, Ysamp2, SH, s, a, ω, m, n)
        low_flux_cutoff = _low_flux_cutoff(Float64(max_flux), Float64(tol))
        relE = _relative_energy_change(res2["EnergyFlux"], res["EnergyFlux"])
        excess = abs(res2["EnergyFlux"]) / flux_scale
        factor = excess <= 1.0 ? 1.0 : min(sqrt(excess), 50.0)
        significant_for_alias = abs(res2["EnergyFlux"]) > low_flux_cutoff
        allow_stop = (N2 >= min_stop_N) || !significant_for_alias
        low_flux_done = mode_abs_floor > 0.0 && abs(res2["EnergyFlux"]) < mode_abs_floor
        has_flux_reference = Float64(max_flux) > 100 * eps(Float64)
        suspect_min_N = has_flux_reference && !low_flux_done ? _suspect_min_sample(min_stop_N, Nmax, excess, 1024, 2048, 4096) : min_stop_N
        suspect_needs_extra = has_flux_reference && !low_flux_done && N2 < suspect_min_N
        if relE == 0.0
            N = N2
            KG_samp = KG_samp2
            Ysamp = Ysamp2
            res = res2
            if allow_stop && !suspect_needs_extra
                break
            end
            continue
        end
        if allow_stop && !suspect_needs_extra && (factor * relE <= effective_sample_tol || abs(res2["EnergyFlux"]) < low_flux_cutoff || low_flux_done)
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

    KG = kerr_geo_orbit(a, p, e, 1.0)

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
        rphase = _radial_phase_vector(KG_samp, omega, m, n)
        Jp_up, Jm_up, drphase_up, rphaseL_up, rphaseR_up, prefactor_up = integrand_eccentric_sample_cheby_p2(KG_samp, Yup_samp, (S0_p2, S1_p2, S2_p2), n)
        integralp_p2 = levin_1d_integral_local_phase(Jp_up, drphase_up, rphase)
        integralm_p2 = levin_1d_integral_local_phase(Jm_up, -drphase_up, -rphase)
        integral_p2 = (integralp_p2 + integralm_p2) * prefactor_up
        hf = horizon_factor(omega, a, m)
        Dict(
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
        Yin_soln = _isem_y_solution(s, l, m, a, omega)
        SH_m2 = spin_weighted_spheroidal_harmonic(s, l, m, a * omega; method = "jacobi")
        Yin_samp = y_sample_m2_isem(Yin_soln, KG_samp)
        S0_m2 = SH_m2(π/2, 0.0)
        S1_m2 = SH_m2(π/2, 0.0; theta_derivative=1)
        S2_m2 = (m^2 + (a * omega)^2 + 2 - 2 * a * omega * m - SH_m2.lambda) * S0_m2
        rphase = _radial_phase_vector(KG_samp, omega, m, n)
        Jp_in, Jm_in, drphase_in, rphaseL_in, rphaseR_in, prefactor_in = integrand_eccentric_sample_cheby_m2(KG_samp, Yin_samp, (S0_m2, S1_m2, S2_m2), n)
        integralp_m2 = levin_1d_integral_local_phase(Jp_in, drphase_in, rphase)
        integralm_m2 = levin_1d_integral_local_phase(Jm_in, -drphase_in, -rphase)
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

function convolution_integral_inclined_trapezoidal_isem(KG_sample::Dict, s, l, m, k, K_sample::Int64; Kmax::Int = 2^12, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false, cache = nothing)
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
        excess = abs(res2["EnergyFlux"]) / flux_scale
        factor = excess <= 1.0 ? 1.0 : min(sqrt(excess), 50.0)
        significant_for_alias = abs(res2["EnergyFlux"]) > low_flux_cutoff
        allow_stop = (K2 >= min_stop_K) || !significant_for_alias
        low_flux_done = mode_abs_floor > 0.0 && abs(res2["EnergyFlux"]) < mode_abs_floor
        has_flux_reference = Float64(max_flux) > 100 * eps(Float64)
        suspect_min_K = has_flux_reference && !low_flux_done ? _suspect_min_sample(min_stop_K, Kmax, excess, 256, 512, 512) : min_stop_K
        suspect_needs_extra = has_flux_reference && !low_flux_done && K2 < suspect_min_K
        if rel == 0.0
            K_interval = K2
            SHsamp = next_SHsamp
            ctx = (KG_samp = next_sample, Ysol = ctx.Ysol, Ydic = ctx.Ydic, SH = ctx.SH, SHsamp = SHsamp, omega = ctx.omega, a = ctx.a, m = ctx.m, Trajectory = ctx.Trajectory)
            res = res2
            if allow_stop && !suspect_needs_extra
                break
            end
            continue
        end
        if allow_stop && !suspect_needs_extra && (factor * rel <= effective_sample_tol || abs(res2["EnergyFlux"]) < low_flux_cutoff || low_flux_done)
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
    KG = kerr_geo_orbit(a, p, 0.0, x)
    KG_sample = kerr_geo_inclined_sample_dense(KG, x, Kmax)
    return convolution_integral_inclined_trapezoidal_isem(KG_sample, s, l, m, k, K_sample; Kmax = Kmax, kwargs...)
end

function convolution_integral_inclined_levin_isem(a, p, x, s, l, m, k, K_sample)
    key = (a, p, x, s, l, m, k, K_sample)
    if _inclined_levin_last_key[] === key
        return _inclined_levin_last_result[]
    end

    KG = kerr_geo_orbit(a, p, 0.0, x)

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
            "Cinc" => GridSampling._isem_gsn_incidence_amplitude(Yup_soln),
            "Y" => Y,
            "Yp" => Yp,
            "X" => X)
        SH_p2_samp = swsh_sample(SH_p2, KG_samp)
        θphase = _polar_phase_vector(KG_samp, omega, m, k)
        Jp_up, Jm_up, dθphase_up, θphaseL_up, θphaseR_up, prefactor_up = integrand_inclined_sample_cheby_p2(KG_samp, Yup, SH_p2_samp, k)
        integralp_up = levin_1d_integral_local_phase(Jp_up, dθphase_up, θphase)
        integralm_up = levin_1d_integral_local_phase(Jm_up, -dθphase_up, -θphase)
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
            "Binc" => GridSampling._isem_gsn_incidence_amplitude(Yin_soln),
            "Y" => Y,
            "Yp" => Yp,
            "X" => X)
        SH_m2_samp = swsh_sample(SH_m2, KG_samp)
        θphase = _polar_phase_vector(KG_samp, omega, m, k)
        Jp_in, Jm_in, dθphase_in, θphaseL_in, θphaseR_in, prefactor_in = integrand_inclined_sample_cheby_m2(KG_samp, Yin, SH_m2_samp, k)
        integralp_in = levin_1d_integral_local_phase(Jp_in, dθphase_in, θphase)
        integralm_in = levin_1d_integral_local_phase(Jm_in, -dθphase_in, -θphase)
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

function convolution_integral_trapezoidal_isem(a, p, e, x, s, l, m, n, k; N = 256, K = 64, Nmax::Int = 2^14, Kmax::Int = 2^12, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false)
    KG = kerr_geo_orbit(a, p, e, x)
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

function convolution_integral_levin_isem(a, p, e, x, s, l, m, n, k;
        N = 256,
        K = 32,
        Nmax::Int = DEFAULT_LEVIN_NMAX,
        Kmax::Int = DEFAULT_LEVIN_KMAX,
        tol = 1e-8,
        sample_tol::Float64 = 1e-3,
        max_flux = 1.0,
        mode_abs_floor::Float64 = 1e-16,
        zero_low_flux::Bool = false,
        threaded_sampling::Bool = false,
        adaptive::Bool = true,
        adaptive_local_n::Int = DEFAULT_ADAPTIVE_LEVIN_LOCAL_N,
        adaptive_min_depth::Int = DEFAULT_ADAPTIVE_LEVIN_MIN_DEPTH,
        adaptive_max_depth::Int = 8,
        adaptive_tol_max::Float64 = DEFAULT_ADAPTIVE_LEVIN_TOL_MAX)
    KG = kerr_geo_orbit(a, p, e, x)
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
        if adaptive
            return convolution_integral_eccentric_adaptive_levin_isem(
                KG, s, l, m, n;
                tol = tol,
                tol0 = sample_tol,
                local_n = adaptive_local_n,
                min_depth = adaptive_min_depth,
                max_depth = adaptive_max_depth,
                depth_tol_max = adaptive_tol_max,
                mode_abs_floor = mode_abs_floor,
                zero_low_flux = zero_low_flux,
                threaded_sampling = threaded_sampling,
            )
        end
        return convolution_integral_eccentric_levin_isem(
            KG_sample, s, l, m, n, N;
            Nmax = Nmax,
            tol = tol,
            sample_tol = sample_tol,
            max_flux = max_flux,
            mode_abs_floor = mode_abs_floor,
            zero_low_flux = zero_low_flux,
            threaded_sampling = threaded_sampling,
        )
    else
        return convolution_integral_generic_levin_isem(a, p, e, x, s, l, m, n, k, N, K; Nmax = Nmax, Kmax = Kmax, tol = tol, sample_tol = sample_tol, max_flux = max_flux, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, adaptive = adaptive, adaptive_local_r_intervals = adaptive_local_n, adaptive_local_theta_intervals = adaptive_local_n, adaptive_min_depth = adaptive_min_depth, adaptive_max_depth = adaptive_max_depth, adaptive_tol_max = adaptive_tol_max)
    end
end

end
