module ModeSummation

using KerrGeodesics
using LsqFit
using HDF5
using Printf
using ..GridSampling
using ..ConvolutionIntegrals

include("SubmissionScope.jl")

export recommended_mode_grid_sizes
export circular_mode_flux
export circular_mode_summation
export eccentric_mode_summation
export inclined_mode_flux
export inclined_mode_summation
export generic_mode_summation

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

_use_levin_tail(index::Int; nmin::Int = 50) = index >= nmin

function _tail_policy(value)
    value === true && return true
    value === false && return false
    value === :auto && return :auto
    throw(ArgumentError("tail_levin must be true, false, or :auto"))
end

@inline _tail_enabled(policy) = policy !== false

@inline function _uniform_radial_n(N0::Int, Nmax::Int, n::Int)
    base = min(max(4N0, 256), Nmax)
    return ConvolutionIntegrals._eccentric_fourier_min_sample(base, Nmax, n)
end

# Cost heuristic behind tail_levin = :auto: a mode with |n| >= nmin switches to adaptive Levin only when the
# uniform trapezoid grid it would need has at least 4096 radial points and is at least 8x (eccentric) or, in N*K,
# 16x (generic) the adaptive-Levin work.
@inline function _ecc_levin_auto(n::Int, N0::Int, Nmax::Int,
        local_n::Int, max_depth::Int, nmin::Int)
    abs(n) >= nmin || return false
    uniform_n = _uniform_radial_n(N0, Nmax, n)
    adaptive_n = local_n * 2^min(max_depth, 5)
    return uniform_n >= max(4096, 8adaptive_n)
end

@inline function _gen_levin_auto(n::Int, N0::Int, K0::Int, Nmax::Int,
        Kmax::Int, local_n::Int, max_depth::Int, nmin::Int)
    abs(n) >= nmin || return false
    uniform_n = _uniform_radial_n(N0, Nmax, n)
    uniform_k = min(max(4K0, 64), Kmax)
    adaptive_n = local_n * 2^min(max_depth, 5)
    adaptive_k = local_n + 1
    return uniform_n >= 4096 &&
        uniform_n * uniform_k >= 16adaptive_n * adaptive_k
end

@inline _mode_cutoff_threshold(base::Float64, mode_abs_floor::Float64) = max(base, mode_abs_floor)
@inline function _adaptive_branch_mode_floor(total_branch_energy::Float64, tol::Real, base_floor::Real)
    base = Float64(base_floor)
    base > 0 || return base
    scaled = abs(total_branch_energy) * Float64(tol)
    scaled > 0 || return base
    return min(base, scaled)
end
@inline function _generic_k_shell_threshold(total_branch_energy::Float64, reference_n_shell::Float64, tol::Real, mode_abs_floor::Float64; scale::Real = 1.0)
    total_scale = abs(total_branch_energy)
    shell_scale = abs(reference_n_shell)
    # Use a branch-level scale for k tails. The current n-shell can be tiny in the tail;
    # multiplying it by tol forces pointless scans into numerical high-k noise.
    energy_scale = max(total_scale, sqrt(total_scale * shell_scale))
    return _mode_cutoff_threshold(Float64(scale) * Float64(tol) * energy_scale, mode_abs_floor)
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

# At a = 0 an inclined orbit is a rotated equatorial one: mode (l, m, k) carries E_eq(l, j) |d^l_{m,j}|^2 with
# |j| = |m + sign(x) k|, so the l spectrum of a k != 0 branch vanishes below and peaks at l = |j|.
@inline _branch_peak_l(lmin::Int, m::Int, k::Int, x) = max(lmin, abs(m + (x < 0 ? -k : k)))

# For fixed k the same rotation gives |j| = |m + sign(x) k|, which falls to 0 at m* = -sign(x) k (where omega ~ 0 and
# the shells are tiny) and rises again past it. On the side of m* an m branch counts toward the m-shell streak only
# once it is past the mirror image of m = 0 about m*, |m| >= 2|k|.
@inline _before_m_peak(m::Int, k::Int, x) = m * (x < 0 ? k : -k) > 0 && abs(m) < 2 * abs(k)

# Eccentric shells and the k = 0 shells of generic orbits: omega = (m Ups_phi + n Ups_r) / Ups_t vanishes at
# m* = -n Ups_r / Ups_phi, where the m shells are tiny, and the modes past it (the lower sideband) carry flux again:
# on the side of m* an m shell counts toward the m streak only once |m| >= |m*|.
@inline _before_omega_zero(m::Int, m_omega_zero::Float64) = m * m_omega_zero > 0 && abs(m) < abs(m_omega_zero)

# Defaults of the generic-orbit stopping parameters `beta` (tail budget) and `deep` (deep stop).
const DEFAULT_GENERIC_BETA = 1e-4
const DEFAULT_GENERIC_DEEP = 1e-5

# Walk ls, calling eval_l(l) (which returns the mode's energy flux), until minimum_consecutive consecutive modes have
# |flux| <= threshold, or (deep > 0) one mode has |flux| <= deep * threshold. Returns the last l walked.
function _l_walk(eval_l, ls, threshold, minimum_consecutive::Int; deep = 0.0)
    below_count = 0
    last_l = first(ls) - step(ls)
    for l in ls
        last_l = l
        value = abs(eval_l(l))
        below_count = value <= threshold ? below_count + 1 : 0
        (below_count >= minimum_consecutive || (deep > 0 && value <= deep * threshold)) && break
    end
    return last_l
end

# l walk of one k != 0 branch; each mode is evaluated at most once, and every segment stops after minimum_consecutive
# modes below threshold. Infinity (from_bottom = false): up from the peak, then down from the peak towards lmin.
# Walks that search for the peak use search_threshold (the plain shell cutoff) without the deep stop, since a
# larger tail budget cannot bound an unvisited interior peak between the stopped bottom and descent walks.
# Horizon (from_bottom = true): Kerr horizon modes below the peak are not suppressed and often peak at lmin, so first
# walk up from lmin as for k = 0; if that stops below the peak, walk down from the peak to where it stopped, then up
# from the peak.
function _branch_l_walk(eval_l, lmin::Int, lpeak::Int, lmax::Int, threshold, minimum_consecutive::Int, from_bottom::Bool; deep = 0.0, search_threshold = threshold)
    top = min(lpeak, lmax + 1) - 1
    if from_bottom
        l_bottom = _l_walk(eval_l, lmin:lmax, search_threshold, minimum_consecutive)
        lpeak > l_bottom || return nothing
        _l_walk(eval_l, top:-1:(l_bottom + 1), search_threshold, minimum_consecutive)
        _l_walk(eval_l, lpeak:lmax, threshold, minimum_consecutive; deep=deep)
    else
        _l_walk(eval_l, lpeak:lmax, threshold, minimum_consecutive; deep=deep)
        _l_walk(eval_l, top:-1:lmin, search_threshold, minimum_consecutive)
    end
    return nothing
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

# Power-of-two trapezoidal grid sizes (N, K) for mode (m, n, k): about `scale` points per unit of the largest
# harmonic index, at least minN and minK. Also returns the unrounded targets.
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

# Classify the orbit as the public ISEM route does (1e-12 tolerance), before allocating a mode cache.
function _mode_flux_cache_family(e, x)
    equatorial = isapprox(abs(x), 1.0; atol = 1e-12)
    circular = isapprox(e, 0.0; atol = 1e-12)
    return circular ? (equatorial ? :circular : :inclined) : (equatorial ? :eccentric : :generic)
end

# Snapshot the parent under its existing lock; invoke user code outside the lock.
function _mode_flux_parent_stop_requested(scope)
    cancelled, callback = lock(scope.lock) do
        (scope.cancel_requested, scope.stop_requested)
    end
    return cancelled || (callback !== nothing && callback() === true)
end

function _circular_mode_flux_scoped(scope, a, p, l, m; x=1.0)
    family = _mode_flux_cache_family(0.0, x)
    share_geometry = isequal(scope.parameters[3], 0.0)
    flux_inf = with_submission_mode(scope; family = family, share_geometry = share_geometry, mode = (-2, l, m, 0, 0)) do mode_cache
        ConvolutionIntegrals.convolution_integral_trapezoidal_isem(a, p, 0.0, x, -2, l, m, 0, 0; cache = mode_cache, Nmax = scope.Nmax, Kmax = scope.Kmax, trajectory = share_geometry ? submission_orbit!(scope) : nothing)
    end
    flux_hor = with_submission_mode(scope; family = family, share_geometry = share_geometry, mode = (2, l, m, 0, 0)) do mode_cache
        ConvolutionIntegrals.convolution_integral_trapezoidal_isem(a, p, 0.0, x, 2, l, m, 0, 0; cache = mode_cache, Nmax = scope.Nmax, Kmax = scope.Kmax, trajectory = share_geometry ? submission_orbit!(scope) : nothing)
    end
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

# Infinity and horizon fluxes of the circular-orbit mode (l, m) (s = -2 and s = +2), in a submission scope.
function circular_mode_flux(a, p, l, m; submission = nothing, x=1.0, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing, Nmax = submission === nothing ? 2^14 : submission.Nmax, Kmax = submission === nothing ? 2^12 : submission.Kmax)
    scope, owned = ensure_public_submission(submission, a, p, 0.0, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        with_submission_operation(scope) do
            _circular_mode_flux_scoped(scope, a, p, l, m; x = x)
        end
    end
end

function _circular_mode_summation_scoped(scope, a, p; x=1.0, tol = 1e-8, lmax = 30, min_consecutive = 5)
    lmin = 2
    min_consecutive < 1 && throw(ArgumentError("min_consecutive must be positive"))

    total_infinity_energy = 0.0
    total_infinity_angular = 0.0
    total_infinity_carter = 0.0
    total_horizon_energy = 0.0
    total_horizon_angular = 0.0
    total_horizon_carter = 0.0
    total_modes = 0
    l_reached_inf = nothing
    l_reached_hor = nothing
    mode_cache = Dict{Tuple{Int, Int}, Any}()
    get_mode(l, m) = get!(mode_cache, (l, m)) do
        circular_mode_flux(a, p, l, m; x, submission = scope)
    end

    below_count = 0
    for l in lmin:lmax
        shell_infinity_energy = 0.0
        shell_infinity_angular = 0.0
        shell_infinity_carter = 0.0
        for m in 1:l
            mode = get_mode(l, m)
            shell_infinity_energy += 2 * mode.infinity.energy_flux
            shell_infinity_angular += 2 * mode.infinity.angular_momentum_flux
            shell_infinity_carter += 2 * mode.infinity.carter_constant_flux
        end
        l_reached_inf = l
        total_infinity_energy += shell_infinity_energy
        total_infinity_angular += shell_infinity_angular
        total_infinity_carter += shell_infinity_carter

        threshold = tol * abs(total_infinity_energy)
        below_count = abs(shell_infinity_energy) <= threshold ? below_count + 1 : 0
        if below_count >= min_consecutive
            break
        end
    end

    below_count = 0
    for l in lmin:lmax
        shell_horizon_energy = 0.0
        shell_horizon_angular = 0.0
        shell_horizon_carter = 0.0
        for m in 1:l
            mode = get_mode(l, m)
            shell_horizon_energy += 2 * mode.horizon.energy_flux
            shell_horizon_angular += 2 * mode.horizon.angular_momentum_flux
            shell_horizon_carter += 2 * mode.horizon.carter_constant_flux
        end
        l_reached_hor = l
        total_horizon_energy += shell_horizon_energy
        total_horizon_angular += shell_horizon_angular
        total_horizon_carter += shell_horizon_carter

        threshold = tol * abs(total_horizon_energy)
        below_count = abs(shell_horizon_energy) <= threshold ? below_count + 1 : 0
        if below_count >= min_consecutive
            break
        end
    end

    total_modes = 2 * sum(lmin:max(something(l_reached_inf, lmin - 1), something(l_reached_hor, lmin - 1)))

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
        l_reached = max(something(l_reached_inf, lmin - 1), something(l_reached_hor, lmin - 1)),
        l_reached_inf = l_reached_inf,
        l_reached_hor = l_reached_hor,
    )
end

# Total circular-orbit fluxes: sum m = 1..l (doubled for -m) over l shells until min_consecutive shells are
# below tol times the running total, separately at infinity and at the horizon.
function circular_mode_summation(a, p; submission = nothing, x=1.0, tol = 1e-8, lmax = 30, min_consecutive = 5, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing, Nmax = submission === nothing ? 2^14 : submission.Nmax, Kmax = submission === nothing ? 2^12 : submission.Kmax)
    scope, owned = ensure_public_submission(submission, a, p, 0.0, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        with_submission_operation(scope) do
            _circular_mode_summation_scoped(scope, a, p; x = x, tol = tol, lmax = lmax, min_consecutive = min_consecutive)
        end
    end
end

_record_tag(x) = replace(replace(@sprintf("%.16g", x), "." => "p"), "-" => "m")

function _eccentric_mode_summation_scoped(scope, a, p, e; x=1.0, N = 64, N0 = N, Nmax = 2^14, tol = 1e-8, lmax = 30, nmax = 500, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, tail_levin = :auto, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_local_n::Int = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, levin_max_depth::Int = 8)
    if e == 0.0
        return circular_mode_summation(a, p; x, tol = tol, lmax = lmax, min_consecutive = minimum_consecutive, submission = scope)
    end
    lmax < 2 && throw(ArgumentError("lmax must be at least 2"))
    nmax < 0 && throw(ArgumentError("nmax must be nonnegative"))
    minimum_consecutive < 1 && throw(ArgumentError("minimum_consecutive must be positive"))
    levin_local_n < 1 && throw(ArgumentError("levin_local_n must be positive"))
    levin_max_depth < 0 && throw(ArgumentError("levin_max_depth must be nonnegative"))
    tol > 0 || throw(ArgumentError("tol must be positive"))
    ispow2(N0) || throw(ArgumentError("N0 must be a power of 2"))
    ispow2(Nmax) || throw(ArgumentError("Nmax must be a power of 2"))
    N0 <= Nmax || throw(ArgumentError("N0 must not exceed Nmax"))
    N = N0
    record_path = record_path === nothing ? "eccentric_mode_data_a_$(_record_tag(a))_p_$(_record_tag(p))_e_$(_record_tag(e))$(x == 1 ? "" : "_x_$(_record_tag(x))").h5" : record_path

    KG = submission_orbit!(scope)
    radial_over_azimuthal = KG["Frequencies"]["ϒr"] / KG["Frequencies"]["ϒϕ"]
    KG_sample = submission_master!(scope)
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
    record_h5 = record ? _public_track_resource!(scope, h5open(record_path, "w")) : nothing
    adaptive_local_n = levin_local_n
    tail_policy = _tail_policy(tail_levin)
    tail_policy_inf = _tail_policy(tail_levin_infinity === nothing ? tail_policy : tail_levin_infinity)
    tail_policy_hor = _tail_policy(tail_levin_horizon === nothing ? tail_policy : tail_levin_horizon)
    branch_mode_floor_inf = Ref(Float64(mode_abs_floor))
    branch_mode_floor_hor = Ref(Float64(mode_abs_floor))
    levin_floor_tracks_mode_floor = levin_mode_abs_floor == mode_abs_floor
    branch_mode_floor = s -> s == -2 ? branch_mode_floor_inf[] : branch_mode_floor_hor[]
    branch_levin_floor = s -> levin_floor_tracks_mode_floor ? branch_mode_floor(s) : Float64(levin_mode_abs_floor)
    ci_kwargs = fast ? (Nmax = Nmax, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling) : (Nmax = Nmax, zero_low_flux = zero_low_flux, threaded_sampling = false)
    tail_levin_infinity_active = Ref(false)
    tail_levin_horizon_active = Ref(false)
    tail_levin_infinity_start_n = Ref{Union{Nothing, Int}}(nothing)
    tail_levin_horizon_start_n = Ref{Union{Nothing, Int}}(nothing)
    use_eccentric_tail_levin = (s, n) -> begin
        branch_policy = s == -2 ? tail_policy_inf : tail_policy_hor
        _tail_enabled(branch_policy) && fast || return false
        active = s == -2 ? tail_levin_infinity_active : tail_levin_horizon_active
        start_n = s == -2 ? tail_levin_infinity_start_n : tail_levin_horizon_start_n
        if active[]
            return true
        elseif _use_levin_tail(n; nmin = levin_nmin) &&
                (branch_policy === true || _ecc_levin_auto(
                    n, N, Nmax, adaptive_local_n, levin_max_depth,
                    levin_nmin))
            active[] = true
            start_n[] = n
            return true
        end
        return false
    end
    ecc_eval = (s, l, m, n, max_flux, tail_values) -> begin
        use_levin = use_eccentric_tail_levin(s, n)
        floor = branch_mode_floor(s)
        if use_levin
            with_submission_mode(scope; mode = (s, l, m, n, 0)) do mode_cache
                ConvolutionIntegrals.convolution_integral_eccentric_adaptive_levin_isem(KG_sample, s, l, m, n; local_n = adaptive_local_n, max_depth = levin_max_depth, tol = tol, mode_abs_floor = _tail_levin_mode_floor(floor, branch_levin_floor(s)), zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, cache = mode_cache)
            end
        else
            if fast
                with_submission_mode(scope; mode = (s, l, m, n, 0)) do mode_cache
                    ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, s, l, m, n, N; tol = tol, sample_tol = sample_tol, max_flux = max_flux, Nmax = Nmax, mode_abs_floor = floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, cache = mode_cache)
                end
            else
                with_submission_mode(scope; mode = (s, l, m, n, 0)) do mode_cache
                    ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, s, l, m, n, N; tol = tol, sample_tol = sample_tol, max_flux = max_flux, Nmax = Nmax, mode_abs_floor = floor, zero_low_flux = zero_low_flux, threaded_sampling = false, cache = mode_cache)
                end
            end
        end
    end

    try
        if record
            record_h5["meta/a"] = a
            record_h5["meta/p"] = p
            record_h5["meta/e"] = e
            record_h5["meta/x"] = x
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
            mode = with_submission_mode(scope; mode = (-2, l, m, 0, 0)) do mode_cache
                ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, -2, l, m, 0, N; tol = tol, sample_tol = sample_tol, mode_abs_floor = _adaptive_branch_mode_floor(shell_n_infinity_energy + shell_m_infinity_energy, tol, mode_abs_floor), ci_kwargs..., cache = mode_cache)
            end
            _record_current_mode!(record_h5, record, "I", 0, m, l, mode_index, mode, N)
            shell_m_infinity_energy += 2 * mode["EnergyFlux"]
            shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
            shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
            total_modes += 2
            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_m_infinity_energy), _adaptive_branch_mode_floor(shell_n_infinity_energy + shell_m_infinity_energy, tol, mode_abs_floor))
            below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
            if below_count_l >= minimum_consecutive
                break
            end
        end

        shell_n_infinity_energy += shell_m_infinity_energy
        shell_n_infinity_angular += shell_m_infinity_angular
        shell_n_infinity_carter += shell_m_infinity_carter

        layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_infinity_energy), _adaptive_branch_mode_floor(shell_n_infinity_energy, tol, mode_abs_floor))
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
            mode = with_submission_mode(scope; mode = (2, l, m, 0, 0)) do mode_cache
                ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, 2, l, m, 0, N; tol = tol, sample_tol = sample_tol, mode_abs_floor = _adaptive_branch_mode_floor(shell_n_horizon_energy + shell_m_horizon_energy, tol, mode_abs_floor), ci_kwargs..., cache = mode_cache)
            end
            _record_current_mode!(record_h5, record, "H", 0, m, l, mode_index, mode, N)
            shell_m_horizon_energy += 2 * mode["EnergyFlux"]
            shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
            shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
            total_modes += 2
            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_m_horizon_energy), _adaptive_branch_mode_floor(shell_n_horizon_energy + shell_m_horizon_energy, tol, mode_abs_floor))
            below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
            if below_count_l >= minimum_consecutive
                break
            end
        end

        shell_n_horizon_energy += shell_m_horizon_energy
        shell_n_horizon_angular += shell_m_horizon_angular
        shell_n_horizon_carter += shell_m_horizon_carter

        layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_horizon_energy), _adaptive_branch_mode_floor(shell_n_horizon_energy, tol, mode_abs_floor))
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
    branch_mode_floor_inf[] = _adaptive_branch_mode_floor(total_infinity_energy, tol, mode_abs_floor)
    branch_mode_floor_hor[] = _adaptive_branch_mode_floor(total_horizon_energy, tol, mode_abs_floor)

    for n in 1:nmax
        m_omega_zero = -n * radial_over_azimuthal
        E_estimate_inf = _fit_shell_estimate(model, n_list_inf, Energy_flux_inf, n, p0, abs(shell_n_infinity_energy_last))
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
            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_m_infinity_energy), branch_mode_floor_inf[])
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
                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_infinity_energy_last), branch_mode_floor_inf[])
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_n_infinity_energy += shell_m_infinity_energy
            shell_n_infinity_angular += shell_m_infinity_angular
            shell_n_infinity_carter += shell_m_infinity_carter

            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_infinity_energy_last), branch_mode_floor_inf[])
            below_count_m = !_before_omega_zero(m, m_omega_zero) && (abs(shell_m_infinity_energy) <= layer_threshold) && (shell_n_infinity_energy > min(0.1, exp(1 - 1 / e)) * shell_n_infinity_energy_last) ? below_count_m + 1 : 0
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
                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_infinity_energy_last), branch_mode_floor_inf[])
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_n_infinity_energy += shell_m_infinity_energy
            shell_n_infinity_angular += shell_m_infinity_angular
            shell_n_infinity_carter += shell_m_infinity_carter

            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_infinity_energy_last), branch_mode_floor_inf[])
            below_count_m = !_before_omega_zero(m, m_omega_zero) && (abs(shell_m_infinity_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
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
        branch_mode_floor_inf[] = _adaptive_branch_mode_floor(total_infinity_energy, tol, mode_abs_floor)
        n_threshold = _mode_cutoff_threshold(tol * abs(total_infinity_energy), branch_mode_floor_inf[])
        below_count_n = abs(shell_n_infinity_energy) <= n_threshold ? below_count_n + 1 : 0
        if below_count_n >= minimum_consecutive + 1
            n_reached_inf = n
            break
        end
    end

    below_count_n = 0
    for n in 1:nmax
        m_omega_zero = -n * radial_over_azimuthal
        E_estimate_hor = _fit_shell_estimate(model, n_list_hor, Energy_flux_hor, n, p0, abs(shell_n_horizon_energy_last))
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
            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_m_horizon_energy), branch_mode_floor_hor[])
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
                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_horizon_energy_last), branch_mode_floor_hor[])
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_n_horizon_energy += shell_m_horizon_energy
            shell_n_horizon_angular += shell_m_horizon_angular
            shell_n_horizon_carter += shell_m_horizon_carter

            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_horizon_energy_last), branch_mode_floor_hor[])
            below_count_m = !_before_omega_zero(m, m_omega_zero) && (abs(shell_m_horizon_energy) <= layer_threshold) && (shell_n_horizon_energy > min(0.1, exp(1 - 1 / e)) * shell_n_horizon_energy_last) ? below_count_m + 1 : 0
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
                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_horizon_energy_last), branch_mode_floor_hor[])
                below_count_l = abs(mode["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                if below_count_l >= minimum_consecutive
                    break
                end
            end

            shell_n_horizon_energy += shell_m_horizon_energy
            shell_n_horizon_angular += shell_m_horizon_angular
            shell_n_horizon_carter += shell_m_horizon_carter

            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_horizon_energy_last), branch_mode_floor_hor[])
            below_count_m = !_before_omega_zero(m, m_omega_zero) && (abs(shell_m_horizon_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
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
        branch_mode_floor_hor[] = _adaptive_branch_mode_floor(total_horizon_energy, tol, mode_abs_floor)
        n_threshold = _mode_cutoff_threshold(tol * abs(total_horizon_energy), branch_mode_floor_hor[])
        below_count_n = abs(shell_n_horizon_energy) <= n_threshold ? below_count_n + 1 : 0
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
        record_h5["meta/tail_levin_infinity_enabled"] = _tail_enabled(tail_policy_inf) && fast
        record_h5["meta/tail_levin_horizon_enabled"] = _tail_enabled(tail_policy_hor) && fast
        record_h5["meta/tail_levin_infinity_policy"] = string(tail_policy_inf)
        record_h5["meta/tail_levin_horizon_policy"] = string(tail_policy_hor)
        record_h5["meta/tail_levin_start_inf"] = something(tail_levin_infinity_start_n[], -1)
        record_h5["meta/tail_levin_start_hor"] = something(tail_levin_horizon_start_n[], -1)
        record_h5["meta/tail_levin_nmin"] = levin_nmin
        record_h5["meta/tail_levin_local_n"] = levin_local_n
        record_h5["meta/infinity_truncation_floor"] = branch_mode_floor_inf[]
        record_h5["meta/horizon_truncation_floor"] = branch_mode_floor_hor[]
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
        tail_levin_infinity_enabled = _tail_enabled(tail_policy_inf) && fast,
        tail_levin_horizon_enabled = _tail_enabled(tail_policy_hor) && fast,
        tail_levin_infinity_policy = tail_policy_inf,
        tail_levin_horizon_policy = tail_policy_hor,
        tail_levin_start_inf = tail_levin_infinity_start_n[],
        tail_levin_start_hor = tail_levin_horizon_start_n[],
        tail_levin_nmin = levin_nmin,
        tail_levin_local_n = levin_local_n,
        tail_levin_max_depth = levin_max_depth,
        infinity_truncation_floor = branch_mode_floor_inf[],
        horizon_truncation_floor = branch_mode_floor_hor[],
        infinity_energy_flux_list = Energy_flux_inf,
        horizon_energy_flux_list = Energy_flux_hor,
        record_file = record ? record_path : nothing,
    )
    finally
        record && _public_close_resource!(scope, record_h5)
    end
end

# Total eccentric equatorial fluxes: n = 0 shell, then n = 1, 2, ... shells (each over m >= 1, m <= -1 and l),
# until minimum_consecutive + 1 n shells are below tol times the branch total. See Teukolsky_pointparticle_flux.
function eccentric_mode_summation(a, p, e; submission = nothing, x=1.0, N = 64, N0 = N, Nmax = submission === nothing ? 2^14 : submission.Nmax, tol = 1e-8, lmax = 30, nmax = 500, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, tail_levin = :auto, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_local_n::Int = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, levin_max_depth::Int = 8, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing, Kmax = submission === nothing ? 2^12 : submission.Kmax)
    scope, owned = ensure_public_submission(submission, a, p, e, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        with_submission_operation(scope) do
            _eccentric_mode_summation_scoped(scope, a, p, e; x = x, N = N, N0 = N0, Nmax = Nmax, tol = tol, lmax = lmax, nmax = nmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, tail_levin = tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_local_n = levin_local_n, levin_max_depth = levin_max_depth)
        end
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
    for key in ("YInterpolationPoints", "YInterpolationDoublings",
                "YInterpolationFallback", "YInterpolationSelfCheckError")
        haskey(mode, key) && (bgrp[key] = mode[key])
    end
    return nothing
end

function _inclined_mode_flux_scoped(scope, a, p, x, l, m, k; K = 16, Kmax = 2^12, sample_tol = 1e-3)
    if x == 1.0 || x == -1.0
        if k != 0
            zero_flux = (energy_flux = 0.0, angular_momentum_flux = 0.0, carter_constant_flux = 0.0)
            return (l = l, m = m, k = k, infinity = zero_flux, horizon = zero_flux)
        end
        mode = circular_mode_flux(x == -1.0 ? -a : a, p, l, m; submission = x == -1.0 ? nothing : scope, stop_requested = x == -1.0 ? (() -> _mode_flux_parent_stop_requested(scope)) : nothing, lifecycle_observer = scope.lifecycle_observer, geometry_observer = scope.geometry_observer)
        return (l = l, m = m, k = k, infinity = mode.infinity, horizon = mode.horizon)
    end

    KG = submission_orbit!(scope)
    KG_sample = submission_master!(scope)
    flux_inf = with_submission_mode(scope; mode = (-2, l, m, 0, k)) do mode_cache
        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; Kmax = Kmax, sample_tol = sample_tol, cache = mode_cache)
    end
    flux_hor = with_submission_mode(scope; mode = (2, l, m, 0, k)) do mode_cache
        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; Kmax = Kmax, sample_tol = sample_tol, cache = mode_cache)
    end
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

# Infinity and horizon fluxes of the inclined circular mode (l, m, k), in a submission scope.
function inclined_mode_flux(a, p, x, l, m, k; submission = nothing, K = 16, Kmax = submission === nothing ? 2^12 : submission.Kmax, sample_tol = 1e-3, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing, Nmax = submission === nothing ? 2^14 : submission.Nmax)
    scope, owned = ensure_public_submission(submission, a, p, 0.0, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        with_submission_operation(scope) do
            _inclined_mode_flux_scoped(scope, a, p, x, l, m, k; K = K, Kmax = Kmax, sample_tol = sample_tol)
        end
    end
end

function _inclined_mode_summation_scoped(scope, a, p, x; K = 16, K0 = K, Kmax = 2^12, tol = 1e-8, lmax = 30, kmax = 20, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false)
    lmax < 2 && throw(ArgumentError("lmax must be at least 2"))
    minimum_consecutive < 1 && throw(ArgumentError("minimum_consecutive must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))

    if x == 1.0 || x == -1.0
        return circular_mode_summation(x == -1.0 ? -a : a, p; tol = tol, lmax = lmax, min_consecutive = minimum_consecutive, submission = x == -1.0 ? nothing : scope, stop_requested = x == -1.0 ? (() -> _mode_flux_parent_stop_requested(scope)) : nothing, lifecycle_observer = scope.lifecycle_observer, geometry_observer = scope.geometry_observer)
    end

    kmax < 0 && throw(ArgumentError("kmax must be nonnegative"))
    ispow2(K0) || throw(ArgumentError("K0 must be a power of 2"))
    ispow2(Kmax) || throw(ArgumentError("Kmax must be a power of 2"))
    K0 <= Kmax || throw(ArgumentError("K0 must not exceed Kmax"))
    K = K0
    record_path = record_path === nothing ? "inclined_mode_data_a_$(_record_tag(a))_p_$(_record_tag(p))_x_$(_record_tag(x)).h5" : record_path

    KG = submission_orbit!(scope)
    KG_sample = submission_master!(scope)
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
    record_h5 = record ? _public_track_resource!(scope, h5open(record_path, "w")) : nothing
    ci_kwargs = fast ? (Kmax = Kmax, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling) : (Kmax = Kmax, zero_low_flux = zero_low_flux, threaded_sampling = false)
    branch_mode_floor_inf = Ref(Float64(mode_abs_floor))
    branch_mode_floor_hor = Ref(Float64(mode_abs_floor))

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
                mode_inf = with_submission_mode(scope; mode = (-2, l, m, 0, 0)) do mode_cache
                    ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, 0, K; tol = tol, sample_tol = sample_tol, mode_abs_floor = _adaptive_branch_mode_floor(shell_k_infinity_energy + shell_m_infinity_energy, tol, mode_abs_floor), ci_kwargs..., cache = mode_cache)
                end
                _record_current_mode_inclined!(record_h5, record, "I", 0, m, l, mode_index, mode_inf, K)
                shell_m_infinity_energy += 2 * mode_inf["EnergyFlux"]
                shell_m_infinity_angular += 2 * mode_inf["AngularMomentumFlux"]
                shell_m_infinity_carter += 2 * mode_inf["CarterConstantFlux"]
                total_modes += 2
                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_m_infinity_energy), _adaptive_branch_mode_floor(shell_k_infinity_energy + shell_m_infinity_energy, tol, mode_abs_floor))
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
                mode_hor = with_submission_mode(scope; mode = (2, l, m, 0, 0)) do mode_cache
                    ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, 0, K; tol = tol, sample_tol = sample_tol, mode_abs_floor = _adaptive_branch_mode_floor(shell_k_horizon_energy + shell_m_horizon_energy, tol, mode_abs_floor), ci_kwargs..., cache = mode_cache)
                end
                _record_current_mode_inclined!(record_h5, record, "H", 0, m, l, mode_index, mode_hor, K)
                shell_m_horizon_energy += 2 * mode_hor["EnergyFlux"]
                shell_m_horizon_angular += 2 * mode_hor["AngularMomentumFlux"]
                shell_m_horizon_carter += 2 * mode_hor["CarterConstantFlux"]
                total_modes += 2
                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_m_horizon_energy), _adaptive_branch_mode_floor(shell_k_horizon_energy + shell_m_horizon_energy, tol, mode_abs_floor))
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
        branch_mode_floor_inf[] = _adaptive_branch_mode_floor(total_infinity_energy, tol, mode_abs_floor)
        branch_mode_floor_hor[] = _adaptive_branch_mode_floor(total_horizon_energy, tol, mode_abs_floor)

        for k in 1:kmax
            E_estimate_inf = _fit_shell_estimate(model, k_list_inf, Energy_flux_inf, k, p0, abs(shell_k_infinity_energy_last))
            shell_k_infinity_energy = 0.0
            shell_k_infinity_angular = 0.0
            shell_k_infinity_carter = 0.0
            below_count_m = 0

            shell_m_infinity_energy = 0.0
            shell_m_infinity_angular = 0.0
            shell_m_infinity_carter = 0.0
            below_count_l = 0

            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
            _branch_l_walk(2, _branch_peak_l(2, 0, k, x), lmax, layer_threshold, minimum_consecutive, false) do l
                mode = with_submission_mode(scope; mode = (-2, l, 0, 0, k)) do mode_cache
                    ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, 0, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, mode_abs_floor = branch_mode_floor_inf[], ci_kwargs..., cache = mode_cache)
                end
                _record_current_mode_inclined!(record_h5, record, "I", k, 0, l, mode_index, mode, K)
                shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                mode["EnergyFlux"]
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

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
                _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, layer_threshold, minimum_consecutive, false) do l
                    mode = with_submission_mode(scope; mode = (-2, l, m, 0, k)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, mode_abs_floor = branch_mode_floor_inf[], ci_kwargs..., cache = mode_cache)
                    end
                    _record_current_mode_inclined!(record_h5, record, "I", k, m, l, mode_index, mode, K)
                    shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                    shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    if abs(mode["EnergyFlux"]) > Max_flux
                        Max_flux = abs(mode["EnergyFlux"])
                        Max_m = m
                    end
                    mode["EnergyFlux"]
                end

                shell_k_infinity_energy += shell_m_infinity_energy
                shell_k_infinity_angular += shell_m_infinity_angular
                shell_k_infinity_carter += shell_m_infinity_carter

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
                below_count_m = !_before_m_peak(m, k, x) && (abs(shell_m_infinity_energy) <= layer_threshold) && (shell_k_infinity_energy > 0.1 * shell_k_infinity_energy_last) ? below_count_m + 1 : 0
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

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
                _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, layer_threshold, minimum_consecutive, false) do l
                    mode = with_submission_mode(scope; mode = (-2, l, m, 0, k)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, mode_abs_floor = branch_mode_floor_inf[], ci_kwargs..., cache = mode_cache)
                    end
                    _record_current_mode_inclined!(record_h5, record, "I", k, m, l, mode_index, mode, K)
                    shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                    shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    mode["EnergyFlux"]
                end

                shell_k_infinity_energy += shell_m_infinity_energy
                shell_k_infinity_angular += shell_m_infinity_angular
                shell_k_infinity_carter += shell_m_infinity_carter

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
                below_count_m = !_before_m_peak(m, k, x) && (abs(shell_m_infinity_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
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
            branch_mode_floor_inf[] = _adaptive_branch_mode_floor(total_infinity_energy, tol, mode_abs_floor)
            k_threshold = _mode_cutoff_threshold(tol * abs(total_infinity_energy), branch_mode_floor_inf[])
            below_count_n = abs(shell_k_infinity_energy) <= k_threshold ? below_count_n + 1 : 0
            if below_count_n >= minimum_consecutive + 1
                k_reached_inf = k
                break
            end
        end

        below_count_n = 0
        for k in 1:kmax
            E_estimate_hor = _fit_shell_estimate(model, k_list_hor, Energy_flux_hor, k, p0, abs(shell_k_horizon_energy_last))
            shell_k_horizon_energy = 0.0
            shell_k_horizon_angular = 0.0
            shell_k_horizon_carter = 0.0
            below_count_m = 0

            shell_m_horizon_energy = 0.0
            shell_m_horizon_angular = 0.0
            shell_m_horizon_carter = 0.0
            below_count_l = 0

            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
            _branch_l_walk(2, _branch_peak_l(2, 0, k, x), lmax, layer_threshold, minimum_consecutive, true) do l
                mode = with_submission_mode(scope; mode = (2, l, 0, 0, k)) do mode_cache
                    ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, 0, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, mode_abs_floor = branch_mode_floor_hor[], ci_kwargs..., cache = mode_cache)
                end
                _record_current_mode_inclined!(record_h5, record, "H", k, 0, l, mode_index, mode, K)
                shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                mode["EnergyFlux"]
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

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
                _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, layer_threshold, minimum_consecutive, true) do l
                    mode = with_submission_mode(scope; mode = (2, l, m, 0, k)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, mode_abs_floor = branch_mode_floor_hor[], ci_kwargs..., cache = mode_cache)
                    end
                    _record_current_mode_inclined!(record_h5, record, "H", k, m, l, mode_index, mode, K)
                    shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                    shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    if abs(mode["EnergyFlux"]) > Max_flux
                        Max_flux = abs(mode["EnergyFlux"])
                        Max_m = m
                    end
                    mode["EnergyFlux"]
                end

                shell_k_horizon_energy += shell_m_horizon_energy
                shell_k_horizon_angular += shell_m_horizon_angular
                shell_k_horizon_carter += shell_m_horizon_carter

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
                below_count_m = !_before_m_peak(m, k, x) && (abs(shell_m_horizon_energy) <= layer_threshold) && (shell_k_horizon_energy > 0.1 * shell_k_horizon_energy_last) ? below_count_m + 1 : 0
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

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
                _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, layer_threshold, minimum_consecutive, true) do l
                    mode = with_submission_mode(scope; mode = (2, l, m, 0, k)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, mode_abs_floor = branch_mode_floor_hor[], ci_kwargs..., cache = mode_cache)
                    end
                    _record_current_mode_inclined!(record_h5, record, "H", k, m, l, mode_index, mode, K)
                    shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                    shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    mode["EnergyFlux"]
                end

                shell_k_horizon_energy += shell_m_horizon_energy
                shell_k_horizon_angular += shell_m_horizon_angular
                shell_k_horizon_carter += shell_m_horizon_carter

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
                below_count_m = !_before_m_peak(m, k, x) && (abs(shell_m_horizon_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
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
            branch_mode_floor_hor[] = _adaptive_branch_mode_floor(total_horizon_energy, tol, mode_abs_floor)
            k_threshold = _mode_cutoff_threshold(tol * abs(total_horizon_energy), branch_mode_floor_hor[])
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
        record && _public_close_resource!(scope, record_h5)
    end
end

# Total inclined circular fluxes, summed over k shells as the eccentric summation sums n shells.
function inclined_mode_summation(a, p, x; submission = nothing, K = 16, K0 = K, Kmax = submission === nothing ? 2^12 : submission.Kmax, tol = 1e-8, lmax = 30, kmax = 20, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing, Nmax = submission === nothing ? 2^14 : submission.Nmax)
    scope, owned = ensure_public_submission(submission, a, p, 0.0, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        with_submission_operation(scope) do
            _inclined_mode_summation_scoped(scope, a, p, x; K = K, K0 = K0, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
        end
    end
end

function _generic_mode_summation_scoped(scope, a, p, e, x; N0 = 64, K0 = 16, Nmax = 2^14, Kmax = 2^12, tol = 1e-8, lmax = 30, kmax = 20, nmax = 500, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, neg_branch_scale = 0.1, beta = DEFAULT_GENERIC_BETA, deep = DEFAULT_GENERIC_DEEP, tail_levin = :auto, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_local_n::Int = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, levin_max_depth::Int = 8, progress_interval::Int = 0, progress_path = nothing)
    lmax < 2 && throw(ArgumentError("lmax must be at least 2"))
    kmax < 0 && throw(ArgumentError("kmax must be nonnegative"))
    nmax < 0 && throw(ArgumentError("nmax must be nonnegative"))
    minimum_consecutive < 1 && throw(ArgumentError("minimum_consecutive must be positive"))
    levin_local_n < 1 && throw(ArgumentError("levin_local_n must be positive"))
    levin_max_depth < 0 && throw(ArgumentError("levin_max_depth must be nonnegative"))
    tol > 0 || throw(ArgumentError("tol must be positive"))

    if x == 1.0 || x == -1.0
        return eccentric_mode_summation(x == -1.0 ? -a : a, p, e; N = N0, N0 = N0, Nmax = Nmax, tol = tol, lmax = lmax, nmax = nmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, tail_levin = tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_local_n = levin_local_n, levin_max_depth = levin_max_depth, submission = x == -1.0 ? nothing : scope, stop_requested = x == -1.0 ? (() -> _mode_flux_parent_stop_requested(scope)) : nothing, lifecycle_observer = scope.lifecycle_observer, geometry_observer = scope.geometry_observer)
    end
    if e == 0.0
        return inclined_mode_summation(a, p, x; K = K0, K0 = K0, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, submission = isequal(e, 0.0) ? scope : nothing, stop_requested = isequal(e, 0.0) ? nothing : (() -> _mode_flux_parent_stop_requested(scope)), lifecycle_observer = scope.lifecycle_observer, geometry_observer = scope.geometry_observer)
    end

    record_path = record_path === nothing ? "generic_mode_data_a_$(_record_tag(a))_p_$(_record_tag(p))_e_$(_record_tag(e))_x_$(_record_tag(x)).h5" : record_path

    KG = submission_orbit!(scope)
    radial_over_azimuthal = KG["Frequencies"]["ϒr"] / KG["Frequencies"]["ϒϕ"]

    KG_master = submission_master!(scope)
    tail_policy = _tail_policy(tail_levin)
    tail_policy_inf = _tail_policy(tail_levin_infinity === nothing ? tail_policy : tail_levin_infinity)
    tail_policy_hor = _tail_policy(tail_levin_horizon === nothing ? tail_policy : tail_levin_horizon)
    model(x, p) = p[1] .* (x .+ 1) .^ p[2] .* exp.(-p[3] .* x)
    p0 = [1e-3, 0.0, 0.1]
    mode_index = Ref(0)
    record_h5 = record ? _public_track_resource!(scope, h5open(record_path, "w")) : nothing
    progress_start = time()
    progress_last_time = Ref(progress_start)
    progress_last_index = Ref(0)
    progress_io = progress_path === nothing ? nothing : _public_track_resource!(scope, open(progress_path, "w"))

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
        branch_mode_floor_inf = Ref(Float64(mode_abs_floor))
        branch_mode_floor_hor = Ref(Float64(mode_abs_floor))
        levin_floor_tracks_mode_floor = levin_mode_abs_floor == mode_abs_floor
        branch_mode_floor = s -> s == -2 ? branch_mode_floor_inf[] : branch_mode_floor_hor[]
        branch_levin_floor = s -> levin_floor_tracks_mode_floor ? branch_mode_floor(s) : Float64(levin_mode_abs_floor)
        tail_levin_infinity_active = Ref(false)
        tail_levin_horizon_active = Ref(false)
        tail_levin_infinity_start_n = Ref{Union{Nothing, Int}}(nothing)
        tail_levin_horizon_start_n = Ref{Union{Nothing, Int}}(nothing)
        use_generic_tail_levin = (s, n) -> begin
            branch_policy = s == -2 ? tail_policy_inf : tail_policy_hor
            _tail_enabled(branch_policy) && fast || return false
            active = s == -2 ? tail_levin_infinity_active : tail_levin_horizon_active
            start_n = s == -2 ? tail_levin_infinity_start_n : tail_levin_horizon_start_n
            if active[]
                return true
            elseif _use_levin_tail(n; nmin = levin_nmin) &&
                    (branch_policy === true || _gen_levin_auto(
                        n, N0, K0, Nmax, Kmax,
                        levin_local_n,
                        levin_max_depth, levin_nmin))
                active[] = true
                start_n[] = n
                return true
            end
            return false
        end
        eval_mode = (s, l, m, n, k, max_flux, use_levin=false) -> begin
            auto_levin = use_generic_tail_levin(s, n)
            use_levin = use_levin || auto_levin
            floor = branch_mode_floor(s)
            mode = if fast
                if use_levin
                    levin_flux_scale = max(abs(sum(current_shell_list[])), abs(max_flux))
                    levin_floor = _tail_levin_mode_floor(floor, branch_levin_floor(s))
                    with_submission_mode(scope; mode = (s, l, m, n, k)) do mode_cache
                        ConvolutionIntegrals.generic_mode_flux_from_master_cached_adaptive_levin!(mode_cache, KG_master, s, l, m, n, k; sample_tol = sample_tol, tol = tol, max_flux = levin_flux_scale, mode_abs_floor = levin_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, confirm_low_flux = true, local_r_intervals = levin_local_n, local_theta_intervals = levin_local_n, max_depth = levin_max_depth)
                    end
                else
                    with_submission_mode(scope; mode = (s, l, m, n, k)) do mode_cache
                        ConvolutionIntegrals.generic_mode_flux_from_master_cached!(mode_cache, KG_master, s, l, m, n, k; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, sample_tol = sample_tol, tol = tol, max_flux = max_flux, mode_abs_floor = floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
                    end
                end
            else
                with_submission_mode(scope; mode = (s, l, m, n, k)) do mode_cache
                    ConvolutionIntegrals.generic_mode_flux_from_master(KG_master, s, l, m, n, k; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, sample_tol = sample_tol, tol = tol, max_flux = max_flux, mode_abs_floor = floor, zero_low_flux = zero_low_flux, threaded_sampling = false, geometry_owner = submission_geometry!(scope))
                end
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
                    below_count_l = abs(mode["EnergyFlux"]) <= tol * abs(shell_n) ? below_count_l + 1 : 0
                    if s == -2
                        branch_mode_floor_inf[] = _adaptive_branch_mode_floor(shell_n, tol, mode_abs_floor)
                    else
                        branch_mode_floor_hor[] = _adaptive_branch_mode_floor(shell_n, tol, mode_abs_floor)
                    end
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end
                below_count_m = abs(shell_m) <= tol * abs(shell_n) ? below_count_m + 1 : 0
                if below_count_m >= minimum_consecutive
                    break
                end
            end

            initial_floor = _adaptive_branch_mode_floor(total_energy + shell_n, tol, mode_abs_floor)
            initial_mode_threshold = _mode_cutoff_threshold(tol * max(abs(total_energy + shell_n), abs(shell_n)), initial_floor)
            initial_k_shell_threshold = _generic_k_shell_threshold(total_energy + shell_n, shell_n, tol, initial_floor)
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
                _branch_l_walk(2, _branch_peak_l(2, 0, k, x), lmax, initial_mode_threshold, minimum_consecutive, s == 2; deep=deep) do l
                    mode = eval_mode(s, l, 0, 0, k, max_flux)
                    record_mode(s, l, 0, 0, k, mode)
                    shell_m += 2 * mode["EnergyFlux"]
                    shell_m_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_carter += 2 * mode["CarterConstantFlux"]
                    shell_n += 2 * mode["EnergyFlux"]
                    shell_n_angular += 2 * mode["AngularMomentumFlux"]
                    shell_n_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    mode["EnergyFlux"]
                end
                shell_k += shell_m
                shell_k_angular += shell_m_angular
                shell_k_carter += shell_m_carter

                below_count_m = 0
                previous_shell_m = Inf
                Max_flux = 0.0
                Max_m = 1
                for m in 1:lmax
                    lmin = max(2, abs(m))
                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, initial_mode_threshold, minimum_consecutive, s == 2; deep=deep) do l
                        mode = eval_mode(s, l, m, 0, k, max_flux)
                        record_mode(s, l, m, 0, k, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        mode["EnergyFlux"]
                    end
                    if abs(shell_m) > Max_flux
                        Max_flux = abs(shell_m)
                        Max_m = m
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter
                    # an m shell that is still growing does not count toward the stop (retrograde orbits rise toward m < 0)
                    growing_m = abs(shell_m) > abs(previous_shell_m)
                    previous_shell_m = shell_m
                    below_count_m = !growing_m && !_before_m_peak(m, k, x) && _positive_m_tail_below(shell_m, shell_k, initial_mode_threshold, tol, m, Max_m, mode_abs_floor, false) ? below_count_m + 1 : 0
                    if below_count_m >= minimum_consecutive
                        break
                    end
                end

                below_count_m = 0
                previous_shell_m = Inf
                for m in -1:-1:-lmax
                    lmin = max(2, abs(m))
                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, initial_mode_threshold, minimum_consecutive, s == 2; deep=deep) do l
                        mode = eval_mode(s, l, m, 0, k, max_flux)
                        record_mode(s, l, m, 0, k, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        mode["EnergyFlux"]
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter
                    # an m shell that is still growing does not count toward the stop (retrograde orbits rise toward m < 0)
                    growing_m = abs(shell_m) > abs(previous_shell_m)
                    previous_shell_m = shell_m
                    below_count_m = !growing_m && !_before_m_peak(m, k, x) && _negative_m_tail_below(shell_m, initial_mode_threshold, neg_branch_scale, m, Max_m, false) ? below_count_m + 1 : 0
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
            if s == -2
                branch_mode_floor_inf[] = _adaptive_branch_mode_floor(total_energy, tol, mode_abs_floor)
            else
                branch_mode_floor_hor[] = _adaptive_branch_mode_floor(total_energy, tol, mode_abs_floor)
            end
            shell_last = shell_n


            below_count_n = 0
            for n in 1:nmax
                m_omega_zero = -n * radial_over_azimuthal
                max_flux = 2 * abs(fit_estimate(n_list, shell_list, n, abs(shell_last)))
                shell_n = 0.0
                shell_n_angular = 0.0
                shell_n_carter = 0.0
                floor = branch_mode_floor(s)
                search_threshold = _mode_cutoff_threshold(tol * abs(shell_last), floor)
                threshold = max(search_threshold, beta * tol * abs(total_energy))
                k_shell_scale = sqrt(abs(total_energy) * abs(shell_last))
                # k-shell stops are relative to the branch total; the absolute per-mode floor would cut whole
                # k shells when the branch total is small (horizon totals of 1e-8..1e-7).
                k_shell_threshold = tol * k_shell_scale
                negative_k_shell_threshold = neg_branch_scale * tol * abs(total_energy)
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
                    if below_count_l >= minimum_consecutive || (deep > 0 && abs(mode["EnergyFlux"]) <= deep * threshold)
                        break
                    end
                end

                shell_k += shell_m
                shell_k_angular += shell_m_angular
                shell_k_carter += shell_m_carter

                below_count_m = 0
                previous_shell_m = Inf
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
                        if below_count_l >= minimum_consecutive || (deep > 0 && abs(mode["EnergyFlux"]) <= deep * threshold)
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
                    # an m shell that is still growing does not count toward the stop (retrograde orbits rise toward m < 0)
                    growing_m = abs(shell_m) > abs(previous_shell_m)
                    previous_shell_m = shell_m
                    below_count_m = !growing_m && !_before_omega_zero(m, m_omega_zero) && _positive_m_tail_below(shell_m, shell_k, threshold, tol, m, Max_m, mode_abs_floor, tail_cutoff) ? below_count_m + 1 : 0
                    if below_count_m >= minimum_consecutive
                        break
                    end
                end

                below_count_m = 0
                previous_shell_m = Inf
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
                        if below_count_l >= minimum_consecutive || (deep > 0 && abs(mode["EnergyFlux"]) <= deep * threshold)
                            break
                        end
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter
                    # an m shell that is still growing does not count toward the stop (retrograde orbits rise toward m < 0)
                    growing_m = abs(shell_m) > abs(previous_shell_m)
                    previous_shell_m = shell_m
                    below_count_m = !growing_m && !_before_omega_zero(m, m_omega_zero) && _negative_m_tail_below(shell_m, threshold, neg_branch_scale, m, Max_m, tail_cutoff) ? below_count_m + 1 : 0
                    if below_count_m >= minimum_consecutive
                        break
                    end
                end

                below_count_k = 0
                previous_shell_k = Inf
                for k in 1:kmax
                    shell_k = 0.0
                    shell_k_angular = 0.0
                    shell_k_carter = 0.0

                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    _branch_l_walk(2, _branch_peak_l(2, 0, k, x), lmax, threshold, minimum_consecutive, s == 2; deep=deep, search_threshold=search_threshold) do l
                        mode = eval_mode(s, l, 0, n, k, max_flux)
                        record_mode(s, l, 0, n, k, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        mode["EnergyFlux"]
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter

                    below_count_m = 0
                    previous_shell_m = Inf
                    Max_flux = 0.0
                    Max_m = 1
                    for m in 1:lmax
                        lmin = max(2, abs(m))
                        shell_m = 0.0
                        shell_m_angular = 0.0
                        shell_m_carter = 0.0
                        below_count_l = 0
                        _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, threshold, minimum_consecutive, s == 2; deep=deep, search_threshold=search_threshold) do l
                            mode = eval_mode(s, l, m, n, k, max_flux)
                            record_mode(s, l, m, n, k, mode)
                            shell_m += 2 * mode["EnergyFlux"]
                            shell_m_angular += 2 * mode["AngularMomentumFlux"]
                            shell_m_carter += 2 * mode["CarterConstantFlux"]
                            shell_n += 2 * mode["EnergyFlux"]
                            shell_n_angular += 2 * mode["AngularMomentumFlux"]
                            shell_n_carter += 2 * mode["CarterConstantFlux"]
                            total_modes += 2
                            mode["EnergyFlux"]
                        end
                        if abs(shell_m) > Max_flux
                            Max_flux = abs(shell_m)
                            Max_m = m
                        end
                        shell_k += shell_m
                        shell_k_angular += shell_m_angular
                        shell_k_carter += shell_m_carter
                        # an m shell that is still growing does not count toward the stop (retrograde orbits rise toward m < 0)
                        growing_m = abs(shell_m) > abs(previous_shell_m)
                        previous_shell_m = shell_m
                        below_count_m = !growing_m && !_before_m_peak(m, k, x) && _positive_m_tail_below(shell_m, shell_k, threshold, tol, m, Max_m, mode_abs_floor, tail_cutoff) ? below_count_m + 1 : 0
                        if below_count_m >= minimum_consecutive
                            break
                        end
                    end

                    below_count_m = 0
                    previous_shell_m = Inf
                    for m in -1:-1:-lmax
                        lmin = max(2, abs(m))
                        shell_m = 0.0
                        shell_m_angular = 0.0
                        shell_m_carter = 0.0
                        below_count_l = 0
                        _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, threshold, minimum_consecutive, s == 2; deep=deep, search_threshold=search_threshold) do l
                            mode = eval_mode(s, l, m, n, k, max_flux)
                            record_mode(s, l, m, n, k, mode)
                            shell_m += 2 * mode["EnergyFlux"]
                            shell_m_angular += 2 * mode["AngularMomentumFlux"]
                            shell_m_carter += 2 * mode["CarterConstantFlux"]
                            shell_n += 2 * mode["EnergyFlux"]
                            shell_n_angular += 2 * mode["AngularMomentumFlux"]
                            shell_n_carter += 2 * mode["CarterConstantFlux"]
                            total_modes += 2
                            mode["EnergyFlux"]
                        end
                        shell_k += shell_m
                        shell_k_angular += shell_m_angular
                        shell_k_carter += shell_m_carter
                        # an m shell that is still growing does not count toward the stop (retrograde orbits rise toward m < 0)
                        growing_m = abs(shell_m) > abs(previous_shell_m)
                        previous_shell_m = shell_m
                        below_count_m = !growing_m && !_before_m_peak(m, k, x) && _negative_m_tail_below(shell_m, threshold, neg_branch_scale, m, Max_m, tail_cutoff) ? below_count_m + 1 : 0
                        if below_count_m >= minimum_consecutive
                            break
                        end
                    end

                    # a k shell that is still growing does not count toward the stop: shells can dip, then rise again
                    growing = abs(shell_k) > abs(previous_shell_k)
                    previous_shell_k = shell_k
                    below_count_k = !growing && abs(shell_k) <= k_shell_threshold ? below_count_k + 1 : 0
                    if below_count_k >= minimum_consecutive
                        break
                    end
                end

                below_count_k = 0
                previous_shell_k = Inf
                for k in -1:-1:-kmax
                    shell_k = 0.0
                    shell_k_angular = 0.0
                    shell_k_carter = 0.0

                    shell_m = 0.0
                    shell_m_angular = 0.0
                    shell_m_carter = 0.0
                    below_count_l = 0
                    _branch_l_walk(2, _branch_peak_l(2, 0, k, x), lmax, threshold, minimum_consecutive, s == 2; deep=deep, search_threshold=search_threshold) do l
                        mode = eval_mode(s, l, 0, n, k, max_flux)
                        record_mode(s, l, 0, n, k, mode)
                        shell_m += 2 * mode["EnergyFlux"]
                        shell_m_angular += 2 * mode["AngularMomentumFlux"]
                        shell_m_carter += 2 * mode["CarterConstantFlux"]
                        shell_n += 2 * mode["EnergyFlux"]
                        shell_n_angular += 2 * mode["AngularMomentumFlux"]
                        shell_n_carter += 2 * mode["CarterConstantFlux"]
                        total_modes += 2
                        mode["EnergyFlux"]
                    end
                    shell_k += shell_m
                    shell_k_angular += shell_m_angular
                    shell_k_carter += shell_m_carter

                    below_count_m = 0
                    previous_shell_m = Inf
                    Max_flux = 0.0
                    Max_m = 1
                    for m in 1:lmax
                        lmin = max(2, abs(m))
                        shell_m = 0.0
                        shell_m_angular = 0.0
                        shell_m_carter = 0.0
                        below_count_l = 0
                        _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, threshold, minimum_consecutive, s == 2; deep=deep, search_threshold=search_threshold) do l
                            mode = eval_mode(s, l, m, n, k, max_flux)
                            record_mode(s, l, m, n, k, mode)
                            shell_m += 2 * mode["EnergyFlux"]
                            shell_m_angular += 2 * mode["AngularMomentumFlux"]
                            shell_m_carter += 2 * mode["CarterConstantFlux"]
                            shell_n += 2 * mode["EnergyFlux"]
                            shell_n_angular += 2 * mode["AngularMomentumFlux"]
                            shell_n_carter += 2 * mode["CarterConstantFlux"]
                            total_modes += 2
                            mode["EnergyFlux"]
                        end
                        if abs(shell_m) > Max_flux
                            Max_flux = abs(shell_m)
                            Max_m = m
                        end
                        shell_k += shell_m
                        shell_k_angular += shell_m_angular
                        shell_k_carter += shell_m_carter
                        # an m shell that is still growing does not count toward the stop (retrograde orbits rise toward m < 0)
                        growing_m = abs(shell_m) > abs(previous_shell_m)
                        previous_shell_m = shell_m
                        below_count_m = !growing_m && !_before_m_peak(m, k, x) && _positive_m_tail_below(shell_m, shell_k, threshold, tol, m, Max_m, mode_abs_floor, tail_cutoff) ? below_count_m + 1 : 0
                        if below_count_m >= minimum_consecutive
                            break
                        end
                    end

                    below_count_m = 0
                    previous_shell_m = Inf
                    for m in -1:-1:-lmax
                        lmin = max(2, abs(m))
                        shell_m = 0.0
                        shell_m_angular = 0.0
                        shell_m_carter = 0.0
                        below_count_l = 0
                        _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, threshold, minimum_consecutive, s == 2; deep=deep, search_threshold=search_threshold) do l
                            mode = eval_mode(s, l, m, n, k, max_flux)
                            record_mode(s, l, m, n, k, mode)
                            shell_m += 2 * mode["EnergyFlux"]
                            shell_m_angular += 2 * mode["AngularMomentumFlux"]
                            shell_m_carter += 2 * mode["CarterConstantFlux"]
                            shell_n += 2 * mode["EnergyFlux"]
                            shell_n_angular += 2 * mode["AngularMomentumFlux"]
                            shell_n_carter += 2 * mode["CarterConstantFlux"]
                            total_modes += 2
                            mode["EnergyFlux"]
                        end
                        shell_k += shell_m
                        shell_k_angular += shell_m_angular
                        shell_k_carter += shell_m_carter
                        # an m shell that is still growing does not count toward the stop (retrograde orbits rise toward m < 0)
                        growing_m = abs(shell_m) > abs(previous_shell_m)
                        previous_shell_m = shell_m
                        below_count_m = !growing_m && !_before_m_peak(m, k, x) && _negative_m_tail_below(shell_m, threshold, neg_branch_scale, m, Max_m, tail_cutoff) ? below_count_m + 1 : 0
                        if below_count_m >= minimum_consecutive
                            break
                        end
                    end

                    # a k shell that is still growing does not count toward the stop: shells can dip, then rise again
                    growing = abs(shell_k) > abs(previous_shell_k)
                    previous_shell_k = shell_k
                    below_count_k = !growing && abs(shell_k) <= negative_k_shell_threshold ? below_count_k + 1 : 0
                    if below_count_k >= minimum_consecutive
                        break
                    end
                end

                push!(shell_list, shell_n)
                push!(n_list, n)
                total_energy += shell_n
                total_angular += shell_n_angular
                total_carter += shell_n_carter
                if s == -2
                    branch_mode_floor_inf[] = _adaptive_branch_mode_floor(total_energy, tol, mode_abs_floor)
                else
                    branch_mode_floor_hor[] = _adaptive_branch_mode_floor(total_energy, tol, mode_abs_floor)
                end
                n_threshold = _mode_cutoff_threshold(tol * abs(total_energy), branch_mode_floor(s))
                below_count_n = abs(shell_n) <= n_threshold ? below_count_n + 1 : 0
                shell_last = shell_n

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
                truncation_floor = branch_mode_floor(s),
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
            record_h5["meta/tail_levin_infinity_enabled"] = _tail_enabled(tail_policy_inf) && fast
            record_h5["meta/tail_levin_horizon_enabled"] = _tail_enabled(tail_policy_hor) && fast
            record_h5["meta/tail_levin_infinity_policy"] = string(tail_policy_inf)
            record_h5["meta/tail_levin_horizon_policy"] = string(tail_policy_hor)
            record_h5["meta/tail_levin_start_inf"] = something(tail_levin_infinity_start_n[], -1)
            record_h5["meta/tail_levin_start_hor"] = something(tail_levin_horizon_start_n[], -1)
            record_h5["meta/tail_levin_nmin"] = levin_nmin
            record_h5["meta/tail_levin_local_n"] = levin_local_n
            record_h5["meta/infinity_truncation_floor"] = inf.truncation_floor
            record_h5["meta/horizon_truncation_floor"] = hor.truncation_floor
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
            tail_levin_infinity_enabled = _tail_enabled(tail_policy_inf) && fast,
            tail_levin_horizon_enabled = _tail_enabled(tail_policy_hor) && fast,
            tail_levin_infinity_policy = tail_policy_inf,
            tail_levin_horizon_policy = tail_policy_hor,
            tail_levin_start_inf = tail_levin_infinity_start_n[],
            tail_levin_start_hor = tail_levin_horizon_start_n[],
            tail_levin_nmin = levin_nmin,
            tail_levin_local_n = levin_local_n,
            tail_levin_max_depth = levin_max_depth,
            infinity_truncation_floor = inf.truncation_floor,
            horizon_truncation_floor = hor.truncation_floor,
            infinity_energy_flux_list = inf.shell_list,
            horizon_energy_flux_list = hor.shell_list,
            record_file = record ? record_path : nothing,
        )
    finally
        record && _public_close_resource!(scope, record_h5)
        progress_io !== nothing && _public_close_resource!(scope, progress_io)
    end
end

# Total generic-orbit fluxes, summed over n, k, m and l shells with the stopping rules described in
# Teukolsky_pointparticle_flux. Equatorial and e = 0 orbits are passed to the eccentric and inclined summations.
function generic_mode_summation(a, p, e, x; submission = nothing, N0 = 64, K0 = 16, Nmax = submission === nothing ? 2^14 : submission.Nmax, Kmax = submission === nothing ? 2^12 : submission.Kmax, tol = 1e-8, lmax = 30, kmax = 20, nmax = 500, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, neg_branch_scale = 0.1, beta = DEFAULT_GENERIC_BETA, deep = DEFAULT_GENERIC_DEEP, tail_levin = :auto, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_local_n::Int = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, levin_max_depth::Int = 8, progress_interval::Int = 0, progress_path = nothing, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing)
    scope, owned = ensure_public_submission(submission, a, p, e, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        with_submission_operation(scope) do
            _generic_mode_summation_scoped(scope, a, p, e, x; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, nmax = nmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, neg_branch_scale = neg_branch_scale, beta=beta, deep=deep, tail_levin = tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_local_n = levin_local_n, levin_max_depth = levin_max_depth, progress_interval = progress_interval, progress_path = progress_path)
        end
    end
end

end
