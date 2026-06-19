module ModeSummation

using KerrGeodesics
using LsqFit
using HDF5
using Printf
using ..GridSampling
using ..ConvolutionIntegrals

export recommended_mode_grid_sizes
export circular_mode_flux
export circular_mode_summation
export eccentric_mode_flux
export eccentric_mode_summation
export inclined_mode_flux
export inclined_mode_summation
export generic_mode_summation

_relerr(a, b) = abs(a - b) / max(abs(b), eps(Float64))

_next_pow2_at_least(n::Integer) = max(1, 2 ^ (ceil(Int, log2(max(1, n))) + 1))

function _fit_shell_estimate(model, idxs, vals, next_index, p0, fallback)
    isempty(vals) && return fallback
    w = min(length(vals), 6)
    xdata = idxs[end - w + 1:end]
    ydata = vals[end - w + 1:end]
    try
        fit = curve_fit(model, xdata, ydata, p0)
        est = abs(model(next_index, fit.param))
        return isfinite(est) ? est : fallback
    catch
        return fallback
    end
end

function _median_sorted(vals::Vector{Float64})
    n = length(vals)
    n == 0 && return 0.0
    sorted = sort(vals)
    mid = n ÷ 2
    return isodd(n) ? sorted[mid + 1] : 0.5 * (sorted[mid] + sorted[mid + 1])
end

function _use_levin_tail(index::Int, shell_values::Vector{Float64}; nmin::Int = 50, tol::Real = 1e-8, trend_window::Int = 8, consecutive::Int = 3)
    index >= nmin || return false
    tol > 0 || return false
    trend_window > 0 || return false
    consecutive > 0 || return false
    length(shell_values) >= max(2 * trend_window, consecutive) || return false
    vals = abs.(shell_values)
    all(isfinite, vals) || return false
    total = abs(sum(shell_values))
    total = max(total, maximum(vals), eps(Float64))
    threshold = Float64(tol)^(1 / 4) * total
    recent = vals[end - consecutive + 1:end]
    all(v -> v <= threshold, recent) || return false
    w = min(trend_window, length(vals) ÷ 2)
    prev = vals[end - 2 * w + 1:end - w]
    curr = vals[end - w + 1:end]
    return _median_sorted(curr) < _median_sorted(prev)
end

@inline _mode_cutoff_threshold(base::Float64, mode_abs_floor::Float64) = max(base, mode_abs_floor)
@inline function _generic_k_shell_threshold(total_branch_energy::Float64, reference_n_shell::Float64, tol::Real, mode_abs_floor::Float64; scale::Real = 1.0)
    total_scale = max(abs(total_branch_energy), eps(Float64))
    shell_scale = max(abs(reference_n_shell), eps(Float64))
    # Use a branch-level scale for k tails. The current n-shell can be tiny in the tail;
    # multiplying it by tol forces pointless scans into numerical high-k noise.
    energy_scale = max(total_scale, sqrt(total_scale * shell_scale))
    return _mode_cutoff_threshold(max(Float64(scale) * Float64(tol) * energy_scale, eps(Float64)), mode_abs_floor)
end
@inline function _tail_levin_mode_floor(mode_abs_floor::Real, levin_mode_abs_floor::Real)
    levin_mode_abs_floor > 0 && return Float64(levin_mode_abs_floor)
    return Float64(mode_abs_floor)
end

@inline function _clear_generic_mode_dependent_caches!(cache)
    # These factors depend on (m,n,k) and should not accumulate across a full
    # mode summation. The m-grid matrices are large enough that retaining them
    # across every n-shell can dominate memory; keep only orbit-level samples,
    # segment geometry, phase bases, and weights.
    empty!(cache.mcaches)
    empty!(cache.levin_phase_factors)
    return nothing
end

@inline function _positive_m_tail_below(shell_m::Float64, shell_k::Float64, threshold::Float64, tol::Real, m::Int, max_m::Int, mode_abs_floor::Float64, tail_cutoff::Bool)
    abs(shell_m) <= threshold || return false
    m > max_m + 5 || return false
    tail_cutoff && return true
    threshold <= 10.0 * max(mode_abs_floor, eps(Float64)) && return true
    return abs(shell_m) < abs(shell_k) * tol
end

@inline function _negative_m_tail_below(shell_m::Float64, threshold::Float64, neg_branch_scale::Real, m::Int, max_m::Int, tail_cutoff::Bool)
    abs(shell_m) <= neg_branch_scale * threshold || return false
    return tail_cutoff || abs(m) >= max_m
end

function _record_adaptive_levin_metadata!(bgrp, mode)
    for key in (
        "AdaptiveLevin",
        "AdaptiveLevinMaxDepth",
        "AdaptiveLevinMaxLeafDepth",
        "AdaptiveLevinEffectiveIntervals",
        "AdaptiveLevinEffectiveRIntervals",
        "AdaptiveLevinEffectiveThetaIntervals",
        "AdaptiveLevinAcceptedSegments",
        "AdaptiveLevinSegmentEvaluations",
        "AdaptiveLevinSplitCount",
        "AdaptiveLevinStopReason",
    )
        haskey(mode, key) && (bgrp[key] = mode[key])
    end
    return nothing
end

function _record_current_mode!(h5, record::Bool, boundary::AbstractString, n::Int, m::Int, l::Int, mode_index::Base.RefValue{Int}, mode, N_sample::Int)
    record || return nothing
    mode_index[] += 1
    nkey = @sprintf("n_%06d", n)
    mkey = @sprintf("m_%+06d", m)
    lkey = @sprintf("l_%04d", l)
    ngrp = haskey(h5, nkey) ? h5[nkey] : create_group(h5, nkey)
    mgrp = haskey(ngrp, mkey) ? ngrp[mkey] : create_group(ngrp, mkey)
    lgrp = haskey(mgrp, lkey) ? mgrp[lkey] : create_group(mgrp, lkey)
    bgrp = haskey(lgrp, boundary) ? lgrp[boundary] : create_group(lgrp, boundary)
    bgrp["source"] = "julia"
    bgrp["mode_index"] = mode_index[]
    bgrp["boundary"] = boundary
    bgrp["l"] = l
    bgrp["m"] = m
    bgrp["n"] = n
    bgrp["N_sample"] = get(mode, "N_sample", N_sample)
    bgrp["Amplitude_re"] = real(get(mode, "Amplitude", 0.0 + 0.0im))
    bgrp["Amplitude_im"] = imag(get(mode, "Amplitude", 0.0 + 0.0im))
    bgrp["EnergyFlux"] = get(mode, "EnergyFlux", NaN)
    bgrp["AngularMomentumFlux"] = get(mode, "AngularMomentumFlux", NaN)
    bgrp["CarterConstantFlux"] = get(mode, "CarterConstantFlux", NaN)
    _record_adaptive_levin_metadata!(bgrp, mode)
    return nothing
end

function recommended_mode_grid_sizes(m::Integer, n::Integer, k::Integer; minN::Int = 256, minK::Int = 64, scale::Int = 10)
    targetN = max(minN - 1, scale * max(abs(m), abs(n)))
    targetK = max(minK - 1, scale * max(abs(m), abs(k)))
    return (
        N = _next_pow2_at_least(targetN),
        K = _next_pow2_at_least(targetK),
        targetN = targetN,
        targetK = targetK,
    )
end

function circular_mode_flux(a, p, l, m)
    flux_inf = ConvolutionIntegrals.convolution_integral_trapezoidal_isem(a, p, 0.0, 1.0, -2, l, m, 0, 0)
    flux_hor = ConvolutionIntegrals.convolution_integral_trapezoidal_isem(a, p, 0.0, 1.0, 2, l, m, 0, 0)
    return (
        l = l,
        m = m,
        infinity = (
            energy_flux = flux_inf["EnergyFlux"],
            angular_momentum_flux = flux_inf["AngularMomentumFlux"],
            carter_constant_flux = flux_inf["CarterConstantFlux"],
        ),
        horizon = (
            energy_flux = flux_hor["EnergyFlux"],
            angular_momentum_flux = flux_hor["AngularMomentumFlux"],
            carter_constant_flux = flux_hor["CarterConstantFlux"],
        ),
    )
end

function eccentric_mode_flux(a, p, e, l, m, n)
    N = recommended_mode_grid_sizes(m, n, 0).N
    if l == abs(m)
        N *= 4
    end
    if l == abs(m) + 1
        N *= 2
    end
    flux_inf = ConvolutionIntegrals.convolution_integral_trapezoidal_isem(a, p, e, 1.0, -2, l, m, n, 0; N = N)
    flux_hor = ConvolutionIntegrals.convolution_integral_trapezoidal_isem(a, p, e, 1.0, 2, l, m, n, 0; N = N)
    return (
        l = l,
        m = m,
        n = n,
        k = 0,
        infinity = (
            energy_flux = flux_inf["EnergyFlux"],
            angular_momentum_flux = flux_inf["AngularMomentumFlux"],
            carter_constant_flux = flux_inf["CarterConstantFlux"],
        ),
        horizon = (
            energy_flux = flux_hor["EnergyFlux"],
            angular_momentum_flux = flux_hor["AngularMomentumFlux"],
            carter_constant_flux = flux_hor["CarterConstantFlux"],
        ),
    )
end

function circular_mode_summation(a, p; tol = 1e-8, lmax = 30, min_consecutive = 5)
    lmin = 2
    min_consecutive < 1 && throw(ArgumentError("min_consecutive must be positive"))

    total_infinity_energy = 0.0
    total_infinity_angular = 0.0
    total_infinity_carter = 0.0
    total_horizon_energy = 0.0
    total_horizon_angular = 0.0
    total_horizon_carter = 0.0
    total_modes = 0
    l_reached = nothing
    below_count = 0

    for l in lmin:lmax
        shell_infinity_energy = 0.0
        shell_infinity_angular = 0.0
        shell_infinity_carter = 0.0
        shell_horizon_energy = 0.0
        shell_horizon_angular = 0.0
        shell_horizon_carter = 0.0
        for m in 1:l
            mode = circular_mode_flux(a, p, l, m)
            shell_infinity_energy += 2 * mode.infinity.energy_flux
            shell_infinity_angular += 2 * mode.infinity.angular_momentum_flux
            shell_infinity_carter += 2 * mode.infinity.carter_constant_flux
            shell_horizon_energy += 2 * mode.horizon.energy_flux
            shell_horizon_angular += 2 * mode.horizon.angular_momentum_flux
            shell_horizon_carter += 2 * mode.horizon.carter_constant_flux
            total_modes += 2
        end
        l_reached = l
        total_infinity_energy += shell_infinity_energy
        total_infinity_angular += shell_infinity_angular
        total_infinity_carter += shell_infinity_carter
        total_horizon_energy += shell_horizon_energy
        total_horizon_angular += shell_horizon_angular
        total_horizon_carter += shell_horizon_carter

        threshold = tol * max(abs(total_infinity_energy), eps(Float64))
        below_count = abs(shell_infinity_energy) <= threshold ? below_count + 1 : 0
        if below_count >= min_consecutive
            break
        end
    end

    return (
        infinity_energy_flux = total_infinity_energy,
        infinity_angular_momentum_flux = total_infinity_angular,
        infinity_carter_constant_flux = total_infinity_carter,
        horizon_energy_flux = total_horizon_energy,
        horizon_angular_momentum_flux = total_horizon_angular,
        horizon_carter_constant_flux = total_horizon_carter,
        total_modes = total_modes,
        tolerance = tol,
        lmin = lmin,
        lmax = lmax,
        min_consecutive = min_consecutive,
        l_reached = l_reached,
    )
end

_record_tag(x) = replace(replace(@sprintf("%.16g", x), "." => "p"), "-" => "m")

function eccentric_mode_summation(a, p, e; N = 64, N0 = N, Nmax = 2^14, tol = 1e-8, lmax = 30, nmax = 500, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, tail_levin = true, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_max_depth::Int = 8)
    if e == 0.0
        return circular_mode_summation(a, p; tol = tol, lmax = lmax, min_consecutive = minimum_consecutive)
    end
    lmax < 2 && throw(ArgumentError("lmax must be at least 2"))
    nmax < 0 && throw(ArgumentError("nmax must be nonnegative"))
    minimum_consecutive < 1 && throw(ArgumentError("minimum_consecutive must be positive"))
    levin_max_depth < 0 && throw(ArgumentError("levin_max_depth must be nonnegative"))
    tol > 0 || throw(ArgumentError("tol must be positive"))
    ispow2(N0) || throw(ArgumentError("N0 must be a power of 2"))
    ispow2(Nmax) || throw(ArgumentError("Nmax must be a power of 2"))
    N0 <= Nmax || throw(ArgumentError("N0 must not exceed Nmax"))
    N = N0
    record_path = record_path === nothing ? "eccentric_mode_data_a_$(_record_tag(a))_p_$(_record_tag(p))_e_$(_record_tag(e)).h5" : record_path

    KG = kerr_geo_orbit(a, p, e, 1.0)
    KG_sample = GridSampling.kerr_geo_eccentric_sample_dense(KG, Nmax)
    Energy_flux_inf = Vector{Float64}()
    Energy_flux_hor = Vector{Float64}()
    n_list_inf = Vector{Int64}()
    n_list_hor = Vector{Int64}()
    model(x, p) = p[1] .* (x .+ 1) .^p[2] .* exp.(-p[3] .* x)
    p0 = [1e-3, 0.0, 0.1]

    total_infinity_energy = 0.0
    total_infinity_angular = 0.0
    total_infinity_carter = 0.0
    total_horizon_energy = 0.0
    total_horizon_angular = 0.0
    total_horizon_carter = 0.0
    total_modes = 0
    n_reached_inf = nothing
    n_reached_hor = nothing
    below_count_n = 0

    shell_n_infinity_energy = 0.0
    shell_n_infinity_energy_last = 0.0
    shell_n_infinity_angular = 0.0
    shell_n_infinity_carter = 0.0
    shell_n_horizon_energy = 0.0
    shell_n_horizon_energy_last = 0.0
    shell_n_horizon_angular = 0.0
    shell_n_horizon_carter = 0.0
    below_count_m = 0
    mode_index = Ref(0)
    record_h5 = record ? h5open(record_path, "w") : nothing
    mode_cache = fast ? ConvolutionIntegrals.EccentricFluxCache() : nothing
    adaptive_local_n = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N
    use_tail_levin_infinity = tail_levin_infinity === nothing ? tail_levin : tail_levin_infinity
    use_tail_levin_horizon = tail_levin_horizon === nothing ? tail_levin : tail_levin_horizon
    if (use_tail_levin_infinity || use_tail_levin_horizon) && fast && nmax >= levin_nmin
        ConvolutionIntegrals.prewarm_eccentric_adaptive_levin_segments!(mode_cache, KG_sample; local_n = adaptive_local_n, max_depth = levin_max_depth)
    end
    ci_kwargs = fast ? (Nmax = Nmax, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, cache = mode_cache) : (Nmax = Nmax, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = false)
    tail_levin_infinity_active = Ref(false)
    tail_levin_horizon_active = Ref(false)
    use_eccentric_tail_levin = (s, n, tail_values) -> begin
        branch_tail_levin = s == -2 ? use_tail_levin_infinity : use_tail_levin_horizon
        branch_tail_levin && fast || return false
        active = s == -2 ? tail_levin_infinity_active : tail_levin_horizon_active
        if active[] || _use_levin_tail(n, tail_values; nmin = levin_nmin, tol = tol)
            active[] = true
            return true
        end
        return false
    end
    ecc_eval = (s, l, m, n, max_flux, tail_values) -> begin
        use_levin = use_eccentric_tail_levin(s, n, tail_values)
        if use_levin
            ConvolutionIntegrals.convolution_integral_eccentric_adaptive_levin_isem(KG_sample, s, l, m, n; local_n = adaptive_local_n, max_depth = levin_max_depth, tol = tol, mode_abs_floor = _tail_levin_mode_floor(mode_abs_floor, levin_mode_abs_floor), zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, cache = mode_cache)
        else
            ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, s, l, m, n, N; tol = tol, sample_tol = sample_tol, max_flux = max_flux, ci_kwargs...)
        end
    end

    try
        if record
            record_h5["meta/a"] = a
            record_h5["meta/p"] = p
            record_h5["meta/e"] = e
            record_h5["meta/N_requested"] = N
            record_h5["meta/Nmax"] = Nmax
            record_h5["meta/tol"] = tol
            record_h5["meta/lmax"] = lmax
            record_h5["meta/nmax"] = nmax
            record_h5["meta/minimum_consecutive"] = minimum_consecutive
            record_h5["meta/sample_tol"] = sample_tol
            record_h5["meta/levin_max_depth"] = levin_max_depth
        end

    for m in 1:lmax
        lmin = max(2, abs(m))
        shell_m_infinity_energy = 0.0
        shell_m_infinity_angular = 0.0
        shell_m_infinity_carter = 0.0
        below_count_l = 0

        for l in lmin:lmax
            mode = ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, -2, l, m, 0, N; tol = tol, sample_tol = sample_tol, ci_kwargs...)
            _record_current_mode!(record_h5, record, "I", 0, m, l, mode_index, mode, N)
            shell_m_infinity_energy += 2 * mode["EnergyFlux"]
            shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
            shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
            total_modes += 2
            layer_threshold = tol * max(abs(shell_m_infinity_energy), eps(Float64))
            below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
            if below_count_l >= minimum_consecutive
                break
            end
        end

        shell_n_infinity_energy += shell_m_infinity_energy
        shell_n_infinity_angular += shell_m_infinity_angular
        shell_n_infinity_carter += shell_m_infinity_carter

        layer_threshold = tol * max(abs(shell_n_infinity_energy), eps(Float64))
        below_count_m = abs(shell_m_infinity_energy) <= layer_threshold ? below_count_m + 1 : 0
        if below_count_m >= minimum_consecutive
            shell_n_infinity_energy_last = shell_n_infinity_energy
            break
        end
    end

    below_count_m = 0
    for m in 1:lmax
        lmin = max(2, abs(m))
        shell_m_horizon_energy = 0.0
        shell_m_horizon_angular = 0.0
        shell_m_horizon_carter = 0.0
        below_count_l = 0

        for l in lmin:lmax
            mode = ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, 2, l, m, 0, N; tol = tol, sample_tol = sample_tol, ci_kwargs...)
            _record_current_mode!(record_h5, record, "H", 0, m, l, mode_index, mode, N)
            shell_m_horizon_energy += 2 * mode["EnergyFlux"]
            shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
            shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
            total_modes += 2
            layer_threshold = tol * max(abs(shell_m_horizon_energy), eps(Float64))
            below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
            if below_count_l >= minimum_consecutive
                break
            end
        end

        shell_n_horizon_energy += shell_m_horizon_energy
        shell_n_horizon_angular += shell_m_horizon_angular
        shell_n_horizon_carter += shell_m_horizon_carter

        layer_threshold = tol * max(abs(shell_n_horizon_energy), eps(Float64))
        below_count_m = abs(shell_m_horizon_energy) <= layer_threshold ? below_count_m + 1 : 0
        if below_count_m >= minimum_consecutive
            shell_n_horizon_energy_last = shell_n_horizon_energy
            break
        end
    end

    push!(Energy_flux_inf, shell_n_infinity_energy)
    push!(Energy_flux_hor, shell_n_horizon_energy)
    push!(n_list_inf, 0)
    push!(n_list_hor, 0)
    total_infinity_energy += shell_n_infinity_energy
    total_infinity_angular += shell_n_infinity_angular
    total_infinity_carter += shell_n_infinity_carter
    total_horizon_energy += shell_n_horizon_energy
    total_horizon_angular += shell_n_horizon_angular
    total_horizon_carter += shell_n_horizon_carter

    for n in 1:nmax
        E_estimate_inf = _fit_shell_estimate(model, n_list_inf, Energy_flux_inf, n, p0, max(abs(shell_n_infinity_energy_last), eps(Float64)))
        shell_n_infinity_energy = 0.0
        shell_n_infinity_angular = 0.0
        shell_n_infinity_carter = 0.0
        below_count_m = 0

        shell_m_infinity_energy = 0.0
        shell_m_infinity_angular = 0.0
        shell_m_infinity_carter = 0.0
        below_count_l = 0

        for l in 2:lmax
            mode = ecc_eval(-2, l, 0, n, 2 * E_estimate_inf, Energy_flux_inf)
            _record_current_mode!(record_h5, record, "I", n, 0, l, mode_index, mode, N)
            shell_m_infinity_energy += 2 * mode["EnergyFlux"]
            shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
            shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
            total_modes += 2
            layer_threshold = tol * max(abs(shell_m_infinity_energy), eps(Float64))
            below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
            if below_count_l >= minimum_consecutive
                break
            end
        end

        shell_n_infinity_energy += shell_m_infinity_energy
        shell_n_infinity_angular += shell_m_infinity_angular
        shell_n_infinity_carter += shell_m_infinity_carter

        Max_flux = 0.0
        Max_m = 1
        for m in 1:lmax
            lmin = max(2, abs(m))
            shell_m_infinity_energy = 0.0
            shell_m_infinity_angular = 0.0
            shell_m_infinity_carter = 0.0
            below_count_l = 0

            for l in lmin:lmax
                mode = ecc_eval(-2, l, m, n, 2 * E_estimate_inf, Energy_flux_inf)
                _record_current_mode!(record_h5, record, "I", n, m, l, mode_index, mode, N)
                shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                if abs(mode["EnergyFlux"]) > Max_flux
                    Max_flux = abs(mode["EnergyFlux"])
                    Max_m = m
                end
                layer_threshold = max(tol * abs(shell_n_infinity_energy_last), eps(Float64))
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_n_infinity_energy += shell_m_infinity_energy
            shell_n_infinity_angular += shell_m_infinity_angular
            shell_n_infinity_carter += shell_m_infinity_carter

            layer_threshold = max(tol * abs(shell_n_infinity_energy_last), eps(Float64))
            below_count_m = (abs(shell_m_infinity_energy) <= layer_threshold) && (shell_n_infinity_energy > min(0.1, exp(1 - 1 / e)) * shell_n_infinity_energy_last) ? below_count_m + 1 : 0
            if below_count_m >= minimum_consecutive
                break
            end
        end

        below_count_m = 0
        for m in -1:-1:-lmax
            lmin = max(2, abs(m))
            shell_m_infinity_energy = 0.0
            shell_m_infinity_angular = 0.0
            shell_m_infinity_carter = 0.0
            below_count_l = 0

            for l in lmin:lmax
                mode = ecc_eval(-2, l, m, n, 2 * E_estimate_inf, Energy_flux_inf)
                _record_current_mode!(record_h5, record, "I", n, m, l, mode_index, mode, N)
                shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                layer_threshold = max(tol * abs(shell_n_infinity_energy_last), eps(Float64))
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_n_infinity_energy += shell_m_infinity_energy
            shell_n_infinity_angular += shell_m_infinity_angular
            shell_n_infinity_carter += shell_m_infinity_carter

            layer_threshold = max(tol * abs(shell_n_infinity_energy_last), eps(Float64))
            below_count_m = (abs(shell_m_infinity_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
            if below_count_m >= minimum_consecutive
                shell_n_infinity_energy_last = shell_n_infinity_energy
                break
            end
        end

        push!(Energy_flux_inf, shell_n_infinity_energy)
        push!(n_list_inf, n)
        total_infinity_energy += shell_n_infinity_energy
        total_infinity_angular += shell_n_infinity_angular
        total_infinity_carter += shell_n_infinity_carter
        n_threshold = max(tol * abs(total_infinity_energy), eps(Float64))
        below_count_n = abs(shell_n_infinity_energy) <= n_threshold ? below_count_n + 1 : 0
        if mode_cache isa ConvolutionIntegrals.EccentricFluxCache
            empty!(mode_cache.phase_vectors)
            empty!(mode_cache.levin_phase_factors)
            empty!(mode_cache.levin_factored_phase_factors)
        end
        if below_count_n >= minimum_consecutive + 1
            n_reached_inf = n
            break
        end
    end

    for n in 1:nmax
        E_estimate_hor = _fit_shell_estimate(model, n_list_hor, Energy_flux_hor, n, p0, max(abs(shell_n_horizon_energy_last), eps(Float64)))
        shell_n_horizon_energy = 0.0
        shell_n_horizon_angular = 0.0
        shell_n_horizon_carter = 0.0
        below_count_m = 0

        shell_m_horizon_energy = 0.0
        shell_m_horizon_angular = 0.0
        shell_m_horizon_carter = 0.0
        below_count_l = 0

        for l in 2:lmax
            mode = ecc_eval(2, l, 0, n, 2 * E_estimate_hor, Energy_flux_hor)
            _record_current_mode!(record_h5, record, "H", n, 0, l, mode_index, mode, N)
            shell_m_horizon_energy += 2 * mode["EnergyFlux"]
            shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
            shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
            total_modes += 2
            layer_threshold = tol * max(abs(shell_m_horizon_energy), eps(Float64))
            below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
            if below_count_l >= minimum_consecutive
                break
            end
        end

        shell_n_horizon_energy += shell_m_horizon_energy
        shell_n_horizon_angular += shell_m_horizon_angular
        shell_n_horizon_carter += shell_m_horizon_carter

        Max_flux = 0.0
        Max_m = 1
        for m in 1:lmax
            lmin = max(2, abs(m))
            shell_m_horizon_energy = 0.0
            shell_m_horizon_angular = 0.0
            shell_m_horizon_carter = 0.0
            below_count_l = 0

            for l in lmin:lmax
                mode = ecc_eval(2, l, m, n, 2 * E_estimate_hor, Energy_flux_hor)
                _record_current_mode!(record_h5, record, "H", n, m, l, mode_index, mode, N)
                shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                if abs(mode["EnergyFlux"]) > Max_flux
                    Max_flux = abs(mode["EnergyFlux"])
                    Max_m = m
                end
                layer_threshold = max(tol * abs(shell_n_horizon_energy_last), eps(Float64))
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_n_horizon_energy += shell_m_horizon_energy
            shell_n_horizon_angular += shell_m_horizon_angular
            shell_n_horizon_carter += shell_m_horizon_carter

            layer_threshold = max(tol * abs(shell_n_horizon_energy_last), eps(Float64))
            below_count_m = (abs(shell_m_horizon_energy) <= layer_threshold) && (shell_n_horizon_energy > min(0.1, exp(1 - 1 / e)) * shell_n_horizon_energy_last) ? below_count_m + 1 : 0
            if below_count_m >= minimum_consecutive

                break
            end
        end

        below_count_m = 0
        for m in -1:-1:-lmax
            lmin = max(2, abs(m))
            shell_m_horizon_energy = 0.0
            shell_m_horizon_angular = 0.0
            shell_m_horizon_carter = 0.0
            below_count_l = 0

            for l in lmin:lmax
                mode = ecc_eval(2, l, m, n, 2 * E_estimate_hor, Energy_flux_hor)
                _record_current_mode!(record_h5, record, "H", n, m, l, mode_index, mode, N)
                shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                layer_threshold = max(tol * abs(shell_n_horizon_energy_last), eps(Float64))
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_n_horizon_energy += shell_m_horizon_energy
            shell_n_horizon_angular += shell_m_horizon_angular
            shell_n_horizon_carter += shell_m_horizon_carter

            layer_threshold = max(tol * abs(shell_n_horizon_energy_last), eps(Float64))
            below_count_m = (abs(shell_m_horizon_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
            if below_count_m >= minimum_consecutive
                shell_n_horizon_energy_last = shell_n_horizon_energy
                break
            end
        end

        push!(Energy_flux_hor, shell_n_horizon_energy)
        push!(n_list_hor, n)
        total_horizon_energy += shell_n_horizon_energy
        total_horizon_angular += shell_n_horizon_angular
        total_horizon_carter += shell_n_horizon_carter
        n_threshold = max(tol * abs(total_horizon_energy), eps(Float64))
        below_count_n = abs(shell_n_horizon_energy) <= n_threshold ? below_count_n + 1 : 0
        if mode_cache isa ConvolutionIntegrals.EccentricFluxCache
            empty!(mode_cache.phase_vectors)
            empty!(mode_cache.levin_phase_factors)
            empty!(mode_cache.levin_factored_phase_factors)
        end
        if below_count_n >= minimum_consecutive + 1
            n_reached_hor = n
            break
        end
    end

    if record
        record_h5["meta/infinity_energy_flux"] = total_infinity_energy
        record_h5["meta/infinity_angular_momentum_flux"] = total_infinity_angular
        record_h5["meta/infinity_carter_constant_flux"] = total_infinity_carter
        record_h5["meta/horizon_energy_flux"] = total_horizon_energy
        record_h5["meta/horizon_angular_momentum_flux"] = total_horizon_angular
        record_h5["meta/horizon_carter_constant_flux"] = total_horizon_carter
        record_h5["meta/n_reached_inf"] = something(n_reached_inf, -1)
        record_h5["meta/n_reached_hor"] = something(n_reached_hor, -1)
        record_h5["meta/total_modes"] = total_modes
        record_h5["meta/infinity_energy_flux_list"] = Energy_flux_inf
        record_h5["meta/horizon_energy_flux_list"] = Energy_flux_hor
    end

    return (
        infinity_energy_flux = total_infinity_energy,
        infinity_angular_momentum_flux = total_infinity_angular,
        infinity_carter_constant_flux = total_infinity_carter,
        horizon_energy_flux = total_horizon_energy,
        horizon_angular_momentum_flux = total_horizon_angular,
        horizon_carter_constant_flux = total_horizon_carter,
        total_modes = total_modes,
        tolerance = tol,
        minimum_consecutive = minimum_consecutive,
        l_max = lmax,
        n_max = nmax,
        N_sample = N,
        N_max = Nmax,
        n_reached_inf = n_reached_inf,
        n_reached_hor = n_reached_hor,
        infinity_energy_flux_list = Energy_flux_inf,
        horizon_energy_flux_list = Energy_flux_hor,
        record_file = record ? record_path : nothing,
    )
    finally
        record && close(record_h5)
    end
end

function _record_current_mode_inclined!(h5, record::Bool, boundary::AbstractString, k::Int, m::Int, l::Int, mode_index::Base.RefValue{Int}, mode, K_sample::Int)
    record || return nothing
    mode_index[] += 1
    kkey = @sprintf("k_%06d", k)
    mkey = @sprintf("m_%+06d", m)
    lkey = @sprintf("l_%04d", l)
    kgrp = haskey(h5, kkey) ? h5[kkey] : create_group(h5, kkey)
    mgrp = haskey(kgrp, mkey) ? kgrp[mkey] : create_group(kgrp, mkey)
    lgrp = haskey(mgrp, lkey) ? mgrp[lkey] : create_group(mgrp, lkey)
    bgrp = haskey(lgrp, boundary) ? lgrp[boundary] : create_group(lgrp, boundary)
    bgrp["source"] = "julia"
    bgrp["mode_index"] = mode_index[]
    bgrp["boundary"] = boundary
    bgrp["s"] = boundary == "I" ? -2 : 2
    bgrp["l"] = l
    bgrp["m"] = m
    bgrp["k"] = k
    bgrp["K_sample"] = get(mode, "K_sample", K_sample)
    bgrp["Amplitude_re"] = real(get(mode, "Amplitude", 0.0 + 0.0im))
    bgrp["Amplitude_im"] = imag(get(mode, "Amplitude", 0.0 + 0.0im))
    bgrp["EnergyFlux"] = get(mode, "EnergyFlux", NaN)
    bgrp["AngularMomentumFlux"] = get(mode, "AngularMomentumFlux", NaN)
    bgrp["CarterConstantFlux"] = get(mode, "CarterConstantFlux", NaN)
    _record_adaptive_levin_metadata!(bgrp, mode)
    return nothing
end

function _record_current_mode_generic!(h5, record::Bool, boundary::AbstractString, n::Int, k::Int, m::Int, l::Int, mode_index::Base.RefValue{Int}, mode, N_sample::Int, K_sample::Int)
    mode_index[] += 1
    record || return nothing
    nkey = @sprintf("n_%06d", n)
    kkey = @sprintf("k_%+06d", k)
    mkey = @sprintf("m_%+06d", m)
    lkey = @sprintf("l_%04d", l)
    ngrp = haskey(h5, nkey) ? h5[nkey] : create_group(h5, nkey)
    kgrp = haskey(ngrp, kkey) ? ngrp[kkey] : create_group(ngrp, kkey)
    mgrp = haskey(kgrp, mkey) ? kgrp[mkey] : create_group(kgrp, mkey)
    lgrp = haskey(mgrp, lkey) ? mgrp[lkey] : create_group(mgrp, lkey)
    bgrp = haskey(lgrp, boundary) ? lgrp[boundary] : create_group(lgrp, boundary)
    bgrp["source"] = "julia"
    bgrp["mode_index"] = mode_index[]
    bgrp["boundary"] = boundary
    bgrp["s"] = boundary == "I" ? -2 : 2
    bgrp["l"] = l
    bgrp["m"] = m
    bgrp["n"] = n
    bgrp["k"] = k
    bgrp["N_sample"] = get(mode, "N_sample", N_sample)
    bgrp["K_sample"] = get(mode, "K_sample", K_sample)
    bgrp["Amplitude_re"] = real(get(mode, "Amplitude", 0.0 + 0.0im))
    bgrp["Amplitude_im"] = imag(get(mode, "Amplitude", 0.0 + 0.0im))
    bgrp["EnergyFlux"] = get(mode, "EnergyFlux", NaN)
    bgrp["AngularMomentumFlux"] = get(mode, "AngularMomentumFlux", NaN)
    bgrp["CarterConstantFlux"] = get(mode, "CarterConstantFlux", NaN)
    return nothing
end

function inclined_mode_flux(a, p, x, l, m, k; K = 16, Kmax = 2^12, sample_tol = 1e-3)
    if x == 1.0 || x == -1.0
        if k != 0
            zero_flux = (energy_flux = 0.0, angular_momentum_flux = 0.0, carter_constant_flux = 0.0)
            return (l = l, m = m, k = k, infinity = zero_flux, horizon = zero_flux)
        end
        mode = circular_mode_flux(x == -1.0 ? -a : a, p, l, m)
        return (l = l, m = m, k = k, infinity = mode.infinity, horizon = mode.horizon)
    end

    KG = kerr_geo_orbit(a, p, 0.0, x)
    KG_sample = GridSampling.kerr_geo_inclined_sample_dense(KG, x, Kmax)
    flux_inf = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; Kmax = Kmax, sample_tol = sample_tol)
    flux_hor = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; Kmax = Kmax, sample_tol = sample_tol)
    return (
        l = l,
        m = m,
        k = k,
        infinity = (
            energy_flux = flux_inf["EnergyFlux"],
            angular_momentum_flux = flux_inf["AngularMomentumFlux"],
            carter_constant_flux = flux_inf["CarterConstantFlux"],
        ),
        horizon = (
            energy_flux = flux_hor["EnergyFlux"],
            angular_momentum_flux = flux_hor["AngularMomentumFlux"],
            carter_constant_flux = flux_hor["CarterConstantFlux"],
        ),
    )
end

function inclined_mode_summation(a, p, x; K = 16, K0 = K, Kmax = 2^12, tol = 1e-8, lmax = 30, kmax = 20, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false)
    lmax < 2 && throw(ArgumentError("lmax must be at least 2"))
    minimum_consecutive < 1 && throw(ArgumentError("minimum_consecutive must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))

    if x == 1.0 || x == -1.0
        return circular_mode_summation(x == -1.0 ? -a : a, p; tol = tol, lmax = lmax, min_consecutive = minimum_consecutive)
    end

    kmax < 0 && throw(ArgumentError("kmax must be nonnegative"))
    ispow2(K0) || throw(ArgumentError("K0 must be a power of 2"))
    ispow2(Kmax) || throw(ArgumentError("Kmax must be a power of 2"))
    K0 <= Kmax || throw(ArgumentError("K0 must not exceed Kmax"))
    K = K0
    record_path = record_path === nothing ? "inclined_mode_data_a_$(_record_tag(a))_p_$(_record_tag(p))_x_$(_record_tag(x)).h5" : record_path

    KG = kerr_geo_orbit(a, p, 0.0, x)
    KG_sample = GridSampling.kerr_geo_inclined_sample_dense(KG, x, Kmax)
    Energy_flux_inf = Float64[]
    Energy_flux_hor = Float64[]
    k_list_inf = Int64[]
    k_list_hor = Int64[]
    model(x, p) = p[1] .* (x .+ 1) .^ p[2] .* exp.(-p[3] .* x)
    p0 = [1e-3, 0.0, 0.1]

    total_infinity_energy = 0.0
    total_infinity_angular = 0.0
    total_infinity_carter = 0.0
    total_horizon_energy = 0.0
    total_horizon_angular = 0.0
    total_horizon_carter = 0.0
    total_modes = 0
    k_reached_inf = nothing
    k_reached_hor = nothing
    below_count_n = 0

    mode_index = Ref(0)
    record_h5 = record ? h5open(record_path, "w") : nothing
    mode_cache = fast ? ConvolutionIntegrals.InclinedFluxCache() : nothing
    ci_kwargs = fast ? (Kmax = Kmax, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, cache = mode_cache) : (Kmax = Kmax, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = false)

    try
        if record
            record_h5["meta/a"] = a
            record_h5["meta/p"] = p
            record_h5["meta/x"] = x
            record_h5["meta/K_requested"] = K
            record_h5["meta/Kmax"] = Kmax
            record_h5["meta/tol"] = tol
            record_h5["meta/lmax"] = lmax
            record_h5["meta/kmax"] = kmax
            record_h5["meta/minimum_consecutive"] = minimum_consecutive
            record_h5["meta/sample_tol"] = sample_tol
        end

        shell_k_infinity_energy = 0.0
        shell_k_infinity_angular = 0.0
        shell_k_infinity_carter = 0.0
        shell_k_horizon_energy = 0.0
        shell_k_horizon_angular = 0.0
        shell_k_horizon_carter = 0.0

        for m in 1:lmax
            lmin = max(2, abs(m))
            shell_m_infinity_energy = 0.0
            shell_m_infinity_angular = 0.0
            shell_m_infinity_carter = 0.0
            shell_m_horizon_energy = 0.0
            shell_m_horizon_angular = 0.0
            shell_m_horizon_carter = 0.0
            below_count_l = 0

            for l in lmin:lmax
                mode_inf = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, 0, K; tol = tol, sample_tol = sample_tol, ci_kwargs...)
                _record_current_mode_inclined!(record_h5, record, "I", 0, m, l, mode_index, mode_inf, K)
                shell_m_infinity_energy += 2 * mode_inf["EnergyFlux"]
                shell_m_infinity_angular += 2 * mode_inf["AngularMomentumFlux"]
                shell_m_infinity_carter += 2 * mode_inf["CarterConstantFlux"]
                total_modes += 2
                layer_threshold = tol * max(abs(shell_m_infinity_energy), eps(Float64))
                below_count_l = abs(mode_inf["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_k_infinity_energy += shell_m_infinity_energy
            shell_k_infinity_angular += shell_m_infinity_angular
            shell_k_infinity_carter += shell_m_infinity_carter

            below_count_l = 0
            for l in lmin:lmax
                mode_hor = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, 0, K; tol = tol, sample_tol = sample_tol, ci_kwargs...)
                _record_current_mode_inclined!(record_h5, record, "H", 0, m, l, mode_index, mode_hor, K)
                shell_m_horizon_energy += 2 * mode_hor["EnergyFlux"]
                shell_m_horizon_angular += 2 * mode_hor["AngularMomentumFlux"]
                shell_m_horizon_carter += 2 * mode_hor["CarterConstantFlux"]
                total_modes += 2
                layer_threshold = tol * max(abs(shell_m_horizon_energy), eps(Float64))
                below_count_l = abs(mode_hor["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_k_horizon_energy += shell_m_horizon_energy
            shell_k_horizon_angular += shell_m_horizon_angular
            shell_k_horizon_carter += shell_m_horizon_carter
        end

        push!(Energy_flux_inf, shell_k_infinity_energy)
        push!(Energy_flux_hor, shell_k_horizon_energy)
        push!(k_list_inf, 0)
        push!(k_list_hor, 0)
        total_infinity_energy += shell_k_infinity_energy
        total_infinity_angular += shell_k_infinity_angular
        total_infinity_carter += shell_k_infinity_carter
        total_horizon_energy += shell_k_horizon_energy
        total_horizon_angular += shell_k_horizon_angular
        total_horizon_carter += shell_k_horizon_carter

        shell_k_infinity_energy_last = shell_k_infinity_energy
        shell_k_horizon_energy_last = shell_k_horizon_energy

        for k in 1:kmax
            E_estimate_inf = _fit_shell_estimate(model, k_list_inf, Energy_flux_inf, k, p0, max(abs(shell_k_infinity_energy_last), eps(Float64)))
            shell_k_infinity_energy = 0.0
            shell_k_infinity_angular = 0.0
            shell_k_infinity_carter = 0.0
            below_count_m = 0

            shell_m_infinity_energy = 0.0
            shell_m_infinity_angular = 0.0
            shell_m_infinity_carter = 0.0
            below_count_l = 0

            for l in 2:lmax
                mode = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, 0, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, ci_kwargs...)
                _record_current_mode_inclined!(record_h5, record, "I", k, 0, l, mode_index, mode, K)
                shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                layer_threshold = max(tol * abs(shell_k_infinity_energy_last), eps(Float64))
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_k_infinity_energy += shell_m_infinity_energy
            shell_k_infinity_angular += shell_m_infinity_angular
            shell_k_infinity_carter += shell_m_infinity_carter

            Max_flux = 0.0
            Max_m = 1
            for m in 1:lmax
                lmin = max(2, abs(m))
                shell_m_infinity_energy = 0.0
                shell_m_infinity_angular = 0.0
                shell_m_infinity_carter = 0.0
                below_count_l = 0

                for l in lmin:lmax
                    mode = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, ci_kwargs...)
                    _record_current_mode_inclined!(record_h5, record, "I", k, m, l, mode_index, mode, K)
                    shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                    shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    if abs(mode["EnergyFlux"]) > Max_flux
                        Max_flux = abs(mode["EnergyFlux"])
                        Max_m = m
                    end
                    layer_threshold = max(tol * abs(shell_k_infinity_energy_last), eps(Float64))
                    below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end

                shell_k_infinity_energy += shell_m_infinity_energy
                shell_k_infinity_angular += shell_m_infinity_angular
                shell_k_infinity_carter += shell_m_infinity_carter

                layer_threshold = max(tol * abs(shell_k_infinity_energy_last), eps(Float64))
                below_count_m = (abs(shell_m_infinity_energy) <= layer_threshold) && (shell_k_infinity_energy > 0.1 * shell_k_infinity_energy_last) ? below_count_m + 1 : 0
                if below_count_m >= minimum_consecutive
                    break
                end
            end

            below_count_m = 0
            for m in -1:-1:-lmax
                lmin = max(2, abs(m))
                shell_m_infinity_energy = 0.0
                shell_m_infinity_angular = 0.0
                shell_m_infinity_carter = 0.0
                below_count_l = 0

                for l in lmin:lmax
                    mode = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, ci_kwargs...)
                    _record_current_mode_inclined!(record_h5, record, "I", k, m, l, mode_index, mode, K)
                    shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                    shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    layer_threshold = max(tol * abs(shell_k_infinity_energy_last), eps(Float64))
                    below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end

                shell_k_infinity_energy += shell_m_infinity_energy
                shell_k_infinity_angular += shell_m_infinity_angular
                shell_k_infinity_carter += shell_m_infinity_carter

                layer_threshold = max(tol * abs(shell_k_infinity_energy_last), eps(Float64))
                below_count_m = (abs(shell_m_infinity_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
                if below_count_m >= minimum_consecutive
                    shell_k_infinity_energy_last = shell_k_infinity_energy
                    break
                end
            end

            push!(Energy_flux_inf, shell_k_infinity_energy)
            push!(k_list_inf, k)
            total_infinity_energy += shell_k_infinity_energy
            total_infinity_angular += shell_k_infinity_angular
            total_infinity_carter += shell_k_infinity_carter
            k_threshold = max(tol * abs(total_infinity_energy), eps(Float64))
            below_count_n = abs(shell_k_infinity_energy) <= k_threshold ? below_count_n + 1 : 0
            if below_count_n >= minimum_consecutive + 1
                k_reached_inf = k
                break
            end
        end

        below_count_n = 0
        for k in 1:kmax
            E_estimate_hor = _fit_shell_estimate(model, k_list_hor, Energy_flux_hor, k, p0, max(abs(shell_k_horizon_energy_last), eps(Float64)))
            shell_k_horizon_energy = 0.0
            shell_k_horizon_angular = 0.0
            shell_k_horizon_carter = 0.0
            below_count_m = 0

            shell_m_horizon_energy = 0.0
            shell_m_horizon_angular = 0.0
            shell_m_horizon_carter = 0.0
            below_count_l = 0

            for l in 2:lmax
                mode = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, 0, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, ci_kwargs...)
                _record_current_mode_inclined!(record_h5, record, "H", k, 0, l, mode_index, mode, K)
                shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                layer_threshold = max(tol * abs(shell_k_horizon_energy_last), eps(Float64))
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_k_horizon_energy += shell_m_horizon_energy
            shell_k_horizon_angular += shell_m_horizon_angular
            shell_k_horizon_carter += shell_m_horizon_carter

            Max_flux = 0.0
            Max_m = 1
            for m in 1:lmax
                lmin = max(2, abs(m))
                shell_m_horizon_energy = 0.0
                shell_m_horizon_angular = 0.0
                shell_m_horizon_carter = 0.0
                below_count_l = 0

                for l in lmin:lmax
                    mode = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, ci_kwargs...)
                    _record_current_mode_inclined!(record_h5, record, "H", k, m, l, mode_index, mode, K)
                    shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                    shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    if abs(mode["EnergyFlux"]) > Max_flux
                        Max_flux = abs(mode["EnergyFlux"])
                        Max_m = m
                    end
                    layer_threshold = max(tol * abs(shell_k_horizon_energy_last), eps(Float64))
                    below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end

                shell_k_horizon_energy += shell_m_horizon_energy
                shell_k_horizon_angular += shell_m_horizon_angular
                shell_k_horizon_carter += shell_m_horizon_carter

                layer_threshold = max(tol * abs(shell_k_horizon_energy_last), eps(Float64))
                below_count_m = (abs(shell_m_horizon_energy) <= layer_threshold) && (shell_k_horizon_energy > 0.1 * shell_k_horizon_energy_last) ? below_count_m + 1 : 0
                if below_count_m >= minimum_consecutive
                    break
                end
            end

            below_count_m = 0
            for m in -1:-1:-lmax
                lmin = max(2, abs(m))
                shell_m_horizon_energy = 0.0
                shell_m_horizon_angular = 0.0
                shell_m_horizon_carter = 0.0
                below_count_l = 0

                for l in lmin:lmax
                    mode = ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, ci_kwargs...)
                    _record_current_mode_inclined!(record_h5, record, "H", k, m, l, mode_index, mode, K)
                    shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                    shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    layer_threshold = max(tol * abs(shell_k_horizon_energy_last), eps(Float64))
                    below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end

                shell_k_horizon_energy += shell_m_horizon_energy
                shell_k_horizon_angular += shell_m_horizon_angular
                shell_k_horizon_carter += shell_m_horizon_carter

                layer_threshold = max(tol * abs(shell_k_horizon_energy_last), eps(Float64))
                below_count_m = (abs(shell_m_horizon_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
                if below_count_m >= minimum_consecutive
                    shell_k_horizon_energy_last = shell_k_horizon_energy
                    break
                end
            end

            push!(Energy_flux_hor, shell_k_horizon_energy)
            push!(k_list_hor, k)
            total_horizon_energy += shell_k_horizon_energy
            total_horizon_angular += shell_k_horizon_angular
            total_horizon_carter += shell_k_horizon_carter
            k_threshold = max(tol * abs(total_horizon_energy), eps(Float64))
            below_count_n = abs(shell_k_horizon_energy) <= k_threshold ? below_count_n + 1 : 0
            if below_count_n >= minimum_consecutive + 1
                k_reached_hor = k
                break
            end
        end

        if record
            record_h5["meta/infinity_energy_flux"] = total_infinity_energy
            record_h5["meta/infinity_angular_momentum_flux"] = total_infinity_angular
            record_h5["meta/infinity_carter_constant_flux"] = total_infinity_carter
            record_h5["meta/horizon_energy_flux"] = total_horizon_energy
            record_h5["meta/horizon_angular_momentum_flux"] = total_horizon_angular
            record_h5["meta/horizon_carter_constant_flux"] = total_horizon_carter
            record_h5["meta/k_reached_inf"] = something(k_reached_inf, -1)
            record_h5["meta/k_reached_hor"] = something(k_reached_hor, -1)
            record_h5["meta/total_modes"] = total_modes
            record_h5["meta/infinity_energy_flux_list"] = Energy_flux_inf
            record_h5["meta/horizon_energy_flux_list"] = Energy_flux_hor
        end

        return (
            infinity_energy_flux = total_infinity_energy,
            infinity_angular_momentum_flux = total_infinity_angular,
            infinity_carter_constant_flux = total_infinity_carter,
            horizon_energy_flux = total_horizon_energy,
            horizon_angular_momentum_flux = total_horizon_angular,
            horizon_carter_constant_flux = total_horizon_carter,
            total_modes = total_modes,
            tolerance = tol,
            minimum_consecutive = minimum_consecutive,
            l_max = lmax,
            k_max = kmax,
            K_sample = K,
            K_max = Kmax,
            k_reached_inf = k_reached_inf,
            k_reached_hor = k_reached_hor,
            infinity_energy_flux_list = Energy_flux_inf,
            horizon_energy_flux_list = Energy_flux_hor,
            record_file = record ? record_path : nothing,
        )
    finally
        record && close(record_h5)
    end
end

function generic_mode_summation(a, p, e, x; N0 = 64, K0 = 16, Nmax = 2^14, Kmax = 2^12, tol = 1e-8, lmax = 30, kmax = 20, nmax = 500, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, neg_branch_scale = 0.1, tail_levin = true, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_max_depth::Int = 8, progress_interval::Int = 0, progress_path = nothing)
    lmax < 2 && throw(ArgumentError("lmax must be at least 2"))
    kmax < 0 && throw(ArgumentError("kmax must be nonnegative"))
    nmax < 0 && throw(ArgumentError("nmax must be nonnegative"))
    minimum_consecutive < 1 && throw(ArgumentError("minimum_consecutive must be positive"))
    levin_max_depth < 0 && throw(ArgumentError("levin_max_depth must be nonnegative"))
    tol > 0 || throw(ArgumentError("tol must be positive"))

    if x == 1.0 || x == -1.0
        return eccentric_mode_summation(x == -1.0 ? -a : a, p, e; N = N0, N0 = N0, Nmax = Nmax, tol = tol, lmax = lmax, nmax = nmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, tail_levin = tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_max_depth = levin_max_depth)
    end
    if e == 0.0
        return inclined_mode_summation(a, p, x; K = K0, K0 = K0, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
    end

    record_path = record_path === nothing ? "generic_mode_data_a_$(_record_tag(a))_p_$(_record_tag(p))_e_$(_record_tag(e))_x_$(_record_tag(x)).h5" : record_path

    KG = kerr_geo_orbit(a, p, e, x)
    if typeof(KG) == Vector{String}
        return KG
    end

    KG_master = GridSampling.kerr_geo_generic_sample_dense(KG, Nmax, Kmax)
    mode_cache = Ref(ConvolutionIntegrals.GenericFluxCache())
    use_tail_levin_infinity = tail_levin_infinity === nothing ? tail_levin : tail_levin_infinity
    use_tail_levin_horizon = tail_levin_horizon === nothing ? tail_levin : tail_levin_horizon
    if (use_tail_levin_infinity || use_tail_levin_horizon) && fast && nmax >= levin_nmin
        ConvolutionIntegrals.prewarm_generic_adaptive_levin_segments!(mode_cache[], KG_master; max_depth = levin_max_depth)
    end
    model(x, p) = p[1] .* (x .+ 1) .^ p[2] .* exp.(-p[3] .* x)
    p0 = [1e-3, 0.0, 0.1]
    mode_index = Ref(0)
    record_h5 = record ? h5open(record_path, "w") : nothing
    progress_start = time()
    progress_last_time = Ref(progress_start)
    progress_last_index = Ref(0)
    progress_io = progress_path === nothing ? nothing : open(progress_path, "w")

    fit_estimate = (idxs, vals, n, fallback) -> begin
        length(vals) < 3 && return fallback
        w = min(length(vals), 6)
        xdata = idxs[end - w + 1:end]
        ydata = vals[end - w + 1:end]
        try
            fit = curve_fit(model, xdata, ydata, p0)
            est = abs(model(n, fit.param))
            isfinite(est) ? est : fallback
        catch
            fallback
        end
    end

    try
        if record
            record_h5["meta/a"] = a
            record_h5["meta/p"] = p
            record_h5["meta/e"] = e
            record_h5["meta/x"] = x
            record_h5["meta/N_requested"] = N0
            record_h5["meta/K_requested"] = K0
            record_h5["meta/Nmax"] = Nmax
            record_h5["meta/Kmax"] = Kmax
            record_h5["meta/tol"] = tol
            record_h5["meta/lmax"] = lmax
            record_h5["meta/kmax"] = kmax
            record_h5["meta/nmax"] = nmax
            record_h5["meta/minimum_consecutive"] = minimum_consecutive
            record_h5["meta/sample_tol"] = sample_tol
            record_h5["meta/neg_branch_scale"] = neg_branch_scale
            record_h5["meta/levin_max_depth"] = levin_max_depth
        end
        if progress_io !== nothing
            println(progress_io, join((
                "mode_index", "batch_modes", "batch_seconds", "batch_modes_per_s", "elapsed_seconds",
                "boundary", "s", "l", "m", "n", "k", "energy_flux", "N_sample", "K_sample",
                "adaptive_max_depth", "adaptive_max_leaf_depth", "adaptive_accepted_segments",
                "quadrature", "record_bytes"
            ), '\t'))
            flush(progress_io)
        end

        current_shell_list = Ref{Vector{Float64}}(Float64[])
        tail_levin_infinity_active = Ref(false)
        tail_levin_horizon_active = Ref(false)
        use_generic_tail_levin = (s, n) -> begin
            branch_tail_levin = s == -2 ? use_tail_levin_infinity : use_tail_levin_horizon
            branch_tail_levin && fast || return false
            active = s == -2 ? tail_levin_infinity_active : tail_levin_horizon_active
            if active[] || _use_levin_tail(n, current_shell_list[]; nmin = levin_nmin, tol = tol)
                active[] = true
                return true
            end
            return false
        end
        eval_mode = (s, l, m, n, k, max_flux, use_levin=false) -> begin
            auto_levin = use_generic_tail_levin(s, n)
            use_levin = use_levin || auto_levin
            mode = if fast
                if use_levin
                    levin_flux_scale = max(abs(sum(current_shell_list[])), abs(max_flux), eps(Float64))
                    levin_floor = _tail_levin_mode_floor(mode_abs_floor, levin_mode_abs_floor)
                    ConvolutionIntegrals.generic_mode_flux_from_master_cached_adaptive_levin!(mode_cache[], KG_master, s, l, m, n, k; sample_tol = sample_tol, tol = tol, max_flux = levin_flux_scale, mode_abs_floor = levin_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, confirm_low_flux = true, max_depth = levin_max_depth)
                else
                    ConvolutionIntegrals.generic_mode_flux_from_master_cached!(mode_cache[], KG_master, s, l, m, n, k; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, sample_tol = sample_tol, tol = tol, max_flux = max_flux, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
                end
            else
                ConvolutionIntegrals.generic_mode_flux_from_master(KG_master, s, l, m, n, k; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, sample_tol = sample_tol, tol = tol, max_flux = max_flux, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = false)
            end
            mode
        end

        record_mode = (s, l, m, n, k, mode) -> begin
            _record_current_mode_generic!(record_h5, record, s == -2 ? "I" : "H", n, k, m, l, mode_index, mode, N0, K0)
            if progress_interval > 0 && mode_index[] % progress_interval == 0
                now_t = time()
                batch_modes = mode_index[] - progress_last_index[]
                batch_seconds = now_t - progress_last_time[]
                batch_rate = batch_seconds > 0 ? batch_modes / batch_seconds : Inf
                record_bytes = record && record_path !== nothing && isfile(record_path) ? filesize(record_path) : 0
                line = join((
                    mode_index[],
                    batch_modes,
                    batch_seconds,
                    batch_rate,
                    now_t - progress_start,
                    s == -2 ? "I" : "H",
                    s,
                    l,
                    m,
                    n,
                    k,
                    get(mode, "EnergyFlux", NaN),
                    get(mode, "N_sample", N0),
                    get(mode, "K_sample", K0),
                    get(mode, "AdaptiveLevinMaxDepth", ""),
                    get(mode, "AdaptiveLevinMaxLeafDepth", ""),
                    get(mode, "AdaptiveLevinAcceptedSegments", ""),
                    get(mode, "Quadrature", ""),
                    record_bytes,
                ), '\t')
                println(line)
                flush(stdout)
                if progress_io !== nothing
                    println(progress_io, line)
                    flush(progress_io)
                end
                progress_last_index[] = mode_index[]
                progress_last_time[] = now_t
            end
            return nothing
        end

        function branch_summation(s::Int)
            shell_list = Float64[]
            current_shell_list[] = shell_list
            n_list = Int[]
            total_energy = 0.0
            total_angular = 0.0
            total_carter = 0.0
            total_modes = 0
            n_reached = nothing
            shell_last = 0.0

            shell_n = 0.0
            shell_n_angular = 0.0
            shell_n_carter = 0.0
            below_count_m = 0
            for m in 1:lmax
                lmin = max(2, abs(m))
                shell_m = 0.0
                shell_m_angular = 0.0
                shell_m_carter = 0.0
                below_count_l = 0
                max_flux = 1e-2
                for l in lmin:lmax
                    mode = eval_mode(s, l, m, 0, 0, max_flux)
                    record_mode(s, l, m, 0, 0, mode)
                    shell_m += 2 * mode["EnergyFlux"]
                    shell_m_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_carter += 2 * mode["CarterConstantFlux"]
                    shell_n += 2 * mode["EnergyFlux"]
                    shell_n_angular += 2 * mode["AngularMomentumFlux"]
                    shell_n_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    below_count_l = abs(mode["EnergyFlux"]) <= tol * max(abs(shell_n), eps(Float64)) ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end
                below_count_m = abs(shell_m) <= tol * max(abs(shell_n), eps(Float64)) ? below_count_m + 1 : 0
                if below_count_m >= minimum_consecutive
                    break
                end
            end

            initial_mode_threshold = _mode_cutoff_threshold(max(tol * max(abs(total_energy + shell_n), abs(shell_n), eps(Float64)), eps(Float64)), mode_abs_floor)
            initial_k_shell_threshold = _generic_k_shell_threshold(total_energy + shell_n, shell_n, tol, mode_abs_floor)
            below_count_k = 0
            for k in 1:kmax
                shell_k = 0.0
                shell_k_angular = 0.0
                shell_k_carter = 0.0
                max_flux = abs(shell_n)

                shell_m = 0.0
                shell_m_angular = 0.0
                shell_m_carter = 0.0
                below_count_l = 0
                for l in 2:lmax
                    mode = eval_mode(s, l, 0, 0, k, max_flux)
                    record_mode(s, l, 0, 0, k, mode)
                    shell_m += 2 * mode["EnergyFlux"]
                    shell_m_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_carter += 2 * mode["CarterConstantFlux"]
                    shell_n += 2 * mode["EnergyFlux"]
                    shell_n_angular += 2 * mode["AngularMomentumFlux"]
                    shell_n_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    below_count_l = abs(mode["EnergyFlux"]) <= initial_mode_threshold ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end
                shell_k += shell_m
                shell_k_angular += shell_m_angular
                shell_k_carter += shell_m_carter

                below_count_m = 0
                Max_flux = 0.0
                Max_m = 1
                for m in 1:lmax
                    lmin = max(2, abs(m))
                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    for l in lmin:lmax
                        mode = eval_mode(s, l, m, 0, k, max_flux)
                        record_mode(s, l, m, 0, k, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        below_count_l = abs(mode["EnergyFlux"]) <= initial_mode_threshold ? below_count_l + 1 : 0
                        if below_count_l >= minimum_consecutive
                            break
                        end
                    end
                    if abs(shell_m) > Max_flux
                        Max_flux = abs(shell_m)
                        Max_m = m
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter
                    below_count_m = _positive_m_tail_below(shell_m, shell_k, initial_mode_threshold, tol, m, Max_m, mode_abs_floor, false) ? below_count_m + 1 : 0
                    if below_count_m >= minimum_consecutive
                        break
                    end
                end

                below_count_m = 0
                for m in -1:-1:-lmax
                    lmin = max(2, abs(m))
                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    for l in lmin:lmax
                        mode = eval_mode(s, l, m, 0, k, max_flux)
                        record_mode(s, l, m, 0, k, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        below_count_l = abs(mode["EnergyFlux"]) <= initial_mode_threshold ? below_count_l + 1 : 0
                        if below_count_l >= minimum_consecutive
                            break
                        end
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter
                    below_count_m = _negative_m_tail_below(shell_m, initial_mode_threshold, neg_branch_scale, m, Max_m, false) ? below_count_m + 1 : 0
                    if below_count_m >= minimum_consecutive
                        break
                    end
                end

                below_count_k = abs(shell_k) <= initial_k_shell_threshold ? below_count_k + 1 : 0
                if below_count_k >= minimum_consecutive
                    break
                end
            end

            push!(shell_list, shell_n)
            push!(n_list, 0)
            total_energy += shell_n
            total_angular += shell_n_angular
            total_carter += shell_n_carter
            shell_last = shell_n
            _clear_generic_mode_dependent_caches!(mode_cache[])

            below_count_n = 0
            for n in 1:nmax
                max_flux = 2 * abs(fit_estimate(n_list, shell_list, n, max(abs(shell_last), eps(Float64))))
                shell_n = 0.0
                shell_n_angular = 0.0
                shell_n_carter = 0.0
                threshold = _mode_cutoff_threshold(max(tol * abs(shell_last), eps(Float64)), mode_abs_floor)
                k_shell_scale = sqrt(max(abs(total_energy), eps(Float64)) * max(abs(shell_last), eps(Float64)))
                k_shell_threshold = _mode_cutoff_threshold(max(tol * k_shell_scale, eps(Float64)), mode_abs_floor)
                negative_k_shell_threshold = _mode_cutoff_threshold(max(neg_branch_scale * tol * abs(total_energy), eps(Float64)), mode_abs_floor)
                tail_cutoff = use_generic_tail_levin(s, n)

                shell_k = 0.0
                shell_k_angular = 0.0
                shell_k_carter = 0.0
                shell_m = 0.0
                shell_m_angular = 0.0
                shell_m_carter = 0.0
                below_count_l = 0
                for l in 2:lmax
                    mode = eval_mode(s, l, 0, n, 0, max_flux)
                    record_mode(s, l, 0, n, 0, mode)
                    shell_m += 2 * mode["EnergyFlux"]
                    shell_m_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_carter += 2 * mode["CarterConstantFlux"]
                    shell_n += 2 * mode["EnergyFlux"]
                    shell_n_angular += 2 * mode["AngularMomentumFlux"]
                    shell_n_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end

                shell_k += shell_m
                shell_k_angular += shell_m_angular
                shell_k_carter += shell_m_carter

                below_count_m = 0
                Max_flux = 0.0
                Max_m = 1
                for m in 1:lmax
                    lmin = max(2, abs(m))
                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    for l in lmin:lmax
                        mode = eval_mode(s, l, m, n, 0, max_flux)
                        record_mode(s, l, m, n, 0, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                        if below_count_l >= minimum_consecutive
                            break
                        end
                    end
                    if abs(shell_m) > Max_flux
                        Max_flux = abs(shell_m)
                        Max_m = m
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter
                    below_count_m = _positive_m_tail_below(shell_m, shell_k, threshold, tol, m, Max_m, mode_abs_floor, tail_cutoff) ? below_count_m + 1 : 0
                    if below_count_m >= minimum_consecutive
                        break
                    end
                end

                below_count_m = 0
                for m in -1:-1:-lmax
                    lmin = max(2, abs(m))
                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    for l in lmin:lmax
                        mode = eval_mode(s, l, m, n, 0, max_flux)
                        record_mode(s, l, m, n, 0, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                        if below_count_l >= minimum_consecutive
                            break
                        end
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter
                    below_count_m = _negative_m_tail_below(shell_m, threshold, neg_branch_scale, m, Max_m, tail_cutoff) ? below_count_m + 1 : 0
                    if below_count_m >= minimum_consecutive
                        break
                    end
                end

                below_count_k = 0
                for k in 1:kmax
                    shell_k = 0.0
                    shell_k_angular = 0.0
                    shell_k_carter = 0.0

                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    for l in 2:lmax
                        mode = eval_mode(s, l, 0, n, k, max_flux)
                        record_mode(s, l, 0, n, k, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                        if below_count_l >= minimum_consecutive
                            break
                        end
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter

                    below_count_m = 0
                    Max_flux = 0.0
                    Max_m = 1
                    for m in 1:lmax
                        lmin = max(2, abs(m))
                        shell_m = 0.0
                        shell_m_angular = 0.0
                        shell_m_carter = 0.0
                        below_count_l = 0
                        for l in lmin:lmax
                            mode = eval_mode(s, l, m, n, k, max_flux)
                            record_mode(s, l, m, n, k, mode)
                            shell_m += 2 * mode["EnergyFlux"]
                            shell_m_angular += 2 * mode["AngularMomentumFlux"]
                            shell_m_carter += 2 * mode["CarterConstantFlux"]
                            shell_n += 2 * mode["EnergyFlux"]
                            shell_n_angular += 2 * mode["AngularMomentumFlux"]
                            shell_n_carter += 2 * mode["CarterConstantFlux"]
                            total_modes += 2
                            below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                            if below_count_l >= minimum_consecutive
                                break
                            end
                        end
                        if abs(shell_m) > Max_flux
                            Max_flux = abs(shell_m)
                            Max_m = m
                        end
                        shell_k += shell_m
                        shell_k_angular += shell_m_angular
                        shell_k_carter += shell_m_carter
                        below_count_m = _positive_m_tail_below(shell_m, shell_k, threshold, tol, m, Max_m, mode_abs_floor, tail_cutoff) ? below_count_m + 1 : 0
                        if below_count_m >= minimum_consecutive
                            break
                        end
                    end

                    below_count_m = 0
                    for m in -1:-1:-lmax
                        lmin = max(2, abs(m))
                        shell_m = 0.0
                        shell_m_angular = 0.0
                        shell_m_carter = 0.0
                        below_count_l = 0
                        for l in lmin:lmax
                            mode = eval_mode(s, l, m, n, k, max_flux)
                            record_mode(s, l, m, n, k, mode)
                            shell_m += 2 * mode["EnergyFlux"]
                            shell_m_angular += 2 * mode["AngularMomentumFlux"]
                            shell_m_carter += 2 * mode["CarterConstantFlux"]
                            shell_n += 2 * mode["EnergyFlux"]
                            shell_n_angular += 2 * mode["AngularMomentumFlux"]
                            shell_n_carter += 2 * mode["CarterConstantFlux"]
                            total_modes += 2
                            below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                            if below_count_l >= minimum_consecutive
                                break
                            end
                        end
                        shell_k += shell_m
                        shell_k_angular += shell_m_angular
                        shell_k_carter += shell_m_carter
                        below_count_m = _negative_m_tail_below(shell_m, threshold, neg_branch_scale, m, Max_m, tail_cutoff) ? below_count_m + 1 : 0
                        if below_count_m >= minimum_consecutive
                            break
                        end
                    end

                    below_count_k = abs(shell_k) <= k_shell_threshold ? below_count_k + 1 : 0
                    if below_count_k >= minimum_consecutive
                        break
                    end
                end

                below_count_k = 0
                for k in -1:-1:-kmax
                    shell_k = 0.0
                    shell_k_angular = 0.0
                    shell_k_carter = 0.0

                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    for l in 2:lmax
                        mode = eval_mode(s, l, 0, n, k, max_flux)
                        record_mode(s, l, 0, n, k, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                        if below_count_l >= minimum_consecutive
                            break
                        end
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter

                    below_count_m = 0
                    Max_flux = 0.0
                    Max_m = 1
                    for m in 1:lmax
                        lmin = max(2, abs(m))
                        shell_m = 0.0
                        shell_m_angular = 0.0
                        shell_m_carter = 0.0
                        below_count_l = 0
                        for l in lmin:lmax
                            mode = eval_mode(s, l, m, n, k, max_flux)
                            record_mode(s, l, m, n, k, mode)
                            shell_m += 2 * mode["EnergyFlux"]
                            shell_m_angular += 2 * mode["AngularMomentumFlux"]
                            shell_m_carter += 2 * mode["CarterConstantFlux"]
                            shell_n += 2 * mode["EnergyFlux"]
                            shell_n_angular += 2 * mode["AngularMomentumFlux"]
                            shell_n_carter += 2 * mode["CarterConstantFlux"]
                            total_modes += 2
                            below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                            if below_count_l >= minimum_consecutive
                                break
                            end
                        end
                        if abs(shell_m) > Max_flux
                            Max_flux = abs(shell_m)
                            Max_m = m
                        end
                        shell_k += shell_m
                        shell_k_angular += shell_m_angular
                        shell_k_carter += shell_m_carter
                        below_count_m = _positive_m_tail_below(shell_m, shell_k, threshold, tol, m, Max_m, mode_abs_floor, tail_cutoff) ? below_count_m + 1 : 0
                        if below_count_m >= minimum_consecutive
                            break
                        end
                    end

                    below_count_m = 0
                    for m in -1:-1:-lmax
                        lmin = max(2, abs(m))
                        shell_m = 0.0
                        shell_m_angular = 0.0
                        shell_m_carter = 0.0
                        below_count_l = 0
                        for l in lmin:lmax
                            mode = eval_mode(s, l, m, n, k, max_flux)
                            record_mode(s, l, m, n, k, mode)
                            shell_m += 2 * mode["EnergyFlux"]
                            shell_m_angular += 2 * mode["AngularMomentumFlux"]
                            shell_m_carter += 2 * mode["CarterConstantFlux"]
                            shell_n += 2 * mode["EnergyFlux"]
                            shell_n_angular += 2 * mode["AngularMomentumFlux"]
                            shell_n_carter += 2 * mode["CarterConstantFlux"]
                            total_modes += 2
                            below_count_l = abs(mode["EnergyFlux"]) <= threshold ? below_count_l + 1 : 0
                            if below_count_l >= minimum_consecutive
                                break
                            end
                        end
                        shell_k += shell_m
                        shell_k_angular += shell_m_angular
                        shell_k_carter += shell_m_carter
                        below_count_m = _negative_m_tail_below(shell_m, threshold, neg_branch_scale, m, Max_m, tail_cutoff) ? below_count_m + 1 : 0
                        if below_count_m >= minimum_consecutive
                            break
                        end
                    end

                    below_count_k = abs(shell_k) <= negative_k_shell_threshold ? below_count_k + 1 : 0
                    if below_count_k >= minimum_consecutive
                        break
                    end
                end

                push!(shell_list, shell_n)
                push!(n_list, n)
                total_energy += shell_n
                total_angular += shell_n_angular
                total_carter += shell_n_carter
                n_threshold = _mode_cutoff_threshold(max(tol * abs(total_energy), eps(Float64)), mode_abs_floor)
                below_count_n = abs(shell_n) <= n_threshold ? below_count_n + 1 : 0
                shell_last = shell_n
                _clear_generic_mode_dependent_caches!(mode_cache[])
                if below_count_n >= minimum_consecutive + 1
                    n_reached = n
                    break
                end
            end

            return (
                energy_flux = total_energy,
                angular_momentum_flux = total_angular,
                carter_constant_flux = total_carter,
                total_modes = total_modes,
                n_reached = n_reached,
                shell_list = shell_list,
            )
        end

        inf = branch_summation(-2)
        hor = branch_summation(2)
        total_modes = inf.total_modes + hor.total_modes

        if record
            record_h5["meta/infinity_energy_flux"] = inf.energy_flux
            record_h5["meta/infinity_angular_momentum_flux"] = inf.angular_momentum_flux
            record_h5["meta/infinity_carter_constant_flux"] = inf.carter_constant_flux
            record_h5["meta/horizon_energy_flux"] = hor.energy_flux
            record_h5["meta/horizon_angular_momentum_flux"] = hor.angular_momentum_flux
            record_h5["meta/horizon_carter_constant_flux"] = hor.carter_constant_flux
            record_h5["meta/n_reached_inf"] = something(inf.n_reached, -1)
            record_h5["meta/n_reached_hor"] = something(hor.n_reached, -1)
            record_h5["meta/total_modes"] = total_modes
            record_h5["meta/infinity_energy_flux_list"] = inf.shell_list
            record_h5["meta/horizon_energy_flux_list"] = hor.shell_list
        end

        return (
            infinity_energy_flux = inf.energy_flux,
            infinity_angular_momentum_flux = inf.angular_momentum_flux,
            infinity_carter_constant_flux = inf.carter_constant_flux,
            horizon_energy_flux = hor.energy_flux,
            horizon_angular_momentum_flux = hor.angular_momentum_flux,
            horizon_carter_constant_flux = hor.carter_constant_flux,
            total_modes = total_modes,
            tolerance = tol,
            minimum_consecutive = minimum_consecutive,
            l_max = lmax,
            k_max = kmax,
            n_max = nmax,
            N_sample = N0,
            K_sample = K0,
            N_max = Nmax,
            K_max = Kmax,
            n_reached_inf = inf.n_reached,
            n_reached_hor = hor.n_reached,
            infinity_energy_flux_list = inf.shell_list,
            horizon_energy_flux_list = hor.shell_list,
            record_file = record ? record_path : nothing,
        )
    finally
        record && close(record_h5)
        progress_io !== nothing && close(progress_io)
    end
end

end
