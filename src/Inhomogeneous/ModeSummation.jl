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
# 16x (generic) the adaptive-Levin work. The eccentric rule takes that grid per mode from the radial phase rate
# (uniform_n from _mode_radial_trapezoid_n); the generic n-shell rule keeps 4|n|.
@inline function _ecc_levin_auto(uniform_n::Int, n::Int, local_n::Int, max_depth::Int, nmin::Int)
    abs(n) >= nmin || return false
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

# Radial trapezoid intervals a mode needs: its radial phase rate (ConvolutionIntegrals._radial_phase_min_sample),
# not 4|n|. At high eccentricity the rate is far above 4|n| (about 14 n at e = 0.95 and 30 n at e = 0.99, a = 0.9)
# and grows with the counter-rotating |m|, so the tail switch below is decided per mode from it.
function _mode_radial_trapezoid_n(KG::Dict, m::Int, n::Int, k::Int, N0::Int, Nmax::Int)
    F = KG["Frequencies"]
    ω = (m * F["ϒϕ"] + n * F["ϒr"] + k * get(F, "ϒθ", 0.0)) / F["ϒt"]
    return ConvolutionIntegrals._radial_phase_min_sample(KG, ω, m, n, _uniform_radial_n(N0, Nmax, n), typemax(Int))
end

# Per-mode form of _gen_levin_auto with the radial grid from the phase rate.
@inline function _gen_levin_auto_mode(uniform_n::Int, n::Int, K0::Int, Kmax::Int,
        local_n::Int, max_depth::Int, nmin::Int)
    abs(n) >= nmin || return false
    uniform_k = min(max(4K0, 64), Kmax)
    adaptive_n = local_n * 2^min(max_depth, 5)
    adaptive_k = local_n + 1
    return uniform_n >= 4096 && uniform_n * uniform_k >= 16adaptive_n * adaptive_k
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
@inline _branch_peak_l(lmin::Int, m::Int, k::Int, x) = max(lmin, abs(m + (signbit(x) ? -k : k)))

# For fixed k the same rotation gives |j| = |m + sign(x) k|, which falls to 0 at m* = -sign(x) k (where omega ~ 0 and
# the shells are tiny) and rises again past it. On the side of m* an m branch counts toward the m-shell streak only
# once it is past the mirror image of m = 0 about m*, |m| >= 2|k|.
@inline _before_m_peak(m::Int, k::Int, x) = m * (signbit(x) ? k : -k) > 0 && abs(m) < 2 * abs(k)

# Eccentric shells and the k = 0 shells of generic orbits: omega = (m Ups_phi + n Ups_r) / Ups_t vanishes at
# m* = -n Ups_r / Ups_phi, where the m shells are tiny, and the modes past it (the lower sideband) carry flux again:
# on the side of m* an m shell counts toward the m streak only once |m| >= |m*|.
@inline _before_omega_zero(m::Int, m_omega_zero::Float64) = m * m_omega_zero > 0 && abs(m) < abs(m_omega_zero)

# High-l low-frequency modes with negligible flux are skipped on analytic grounds: see
# ConvolutionIntegrals._negligible_low_frequency_mode.
const _negligible_low_frequency_mode = ConvolutionIntegrals._negligible_low_frequency_mode

# Orbit classes for the sum rules, chosen once from (a, p, e, x) before a sum. Rules that decide which modes are
# computed and when a walk stops may depend on the class; the accuracy of a single mode never does.
#
# high_e (e >= HIGH_E_CLASS_MIN, classes E3 and E4): the m walks of the eccentric sum
#   - take the omega = 0 guard (_before_omega_zero) only when the crossing |m*| = n Ups_r / Ups_phi is below lmax;
#     past it the guard pinned the negative walk at m = -lmax for every n >= 40;
#   - start counting a shell's stop at the last m whose shell was above the layer threshold, plus
#     minimum_consecutive, instead of the previous walk's end (the reach ratchet, kept for every other class,
#     where it was added for nearly circular orbits near the ISSO). The ratchet pinned both walks at |m| = lmax
#     from n ~ 30 on: at (0.9, 8, 0.95, 1) and n <= 400 the two changes cut the modes by 63% and moved E_inf by
#     8.6e-10.
const HIGH_E_CLASS_MIN = 0.7
orbit_rules(a, p, e, x) = (high_e = e >= HIGH_E_CLASS_MIN,)

# An eccentric m walk counts a small m shell toward its stop only once |m| is past this fraction of the |m| the previous
# n shell's walk reached. The m shells of a shell can dip and rise again: for the nearly circular (a, p, e) =
# (0.991, 1.59, 0.029) the n = 1 shell is tiny at |m| = 25..29 and carries 1.5e-8 (2e-7 of E_inf) at l = |m| = 30..38,
# which the walk, stopping in the dip, dropped. The fraction lets the reach shrink from shell to shell.
const MWALK_REACH_FRACTION = 1.0

# A mode may stop refining once two refinement levels are below its floor; the branch floor min(1e-16, tol * branch
# total) also sets the walk thresholds, but as a refinement floor it let thousands of horizon modes near 1e-18 keep
# coarse-grid values whose errors add up: at FEWv2 569994 (a = 0.44, p = 22.7, e = 0.84) E_hor moved by 4.5e-7
# (45 tol) between tol 1e-8 and 1e-10. The refinement floor is this fraction of the branch floor (E_hor at tol 1e-8 then
# within 6e-8 of tol 1e-10, same modes); the walk thresholds are unchanged. Spherical and generic modes use it too.
const MODE_REFINE_FLOOR_FRACTION = 0.01

# Length of the below-threshold run that stops an eccentric n walk: minimum_consecutive + 1, or one period
# Ups_phi / Ups_r of the n shells when that is longer (zoom-whirl, see _eccentric_mode_summation_scoped).
_n_stop_run(minimum_consecutive::Int, radial_over_azimuthal::Real) =
    max(minimum_consecutive + 1, ceil(Int, abs(1 / radial_over_azimuthal)))

# Geometric extrapolation of the eccentric n tail. At high eccentricity the n shells fall off slowly and geometrically,
# and the n walk stops (minimum_consecutive + 1 shells below tol times the total) with several tol of flux still in
# the shells it did not evaluate: at e ~ 0.8, tol = 1e-8 that tail is ~ 1e-7 of E_inf. When the walk stopped by
# itself at n >= N_TAIL_MIN_N for e >= N_TAIL_MIN_E, the ratio r of the last N_TAIL_FIT shells (a least-squares fit
# of log|shell| against n) adds shell_last * r / (1 - r), provided those shells share one sign and r < N_TAIL_MAX_RATIO.
const N_TAIL_MIN_E = 0.7
const N_TAIL_MIN_N = 100
const N_TAIL_FIT = 4
const N_TAIL_MAX_RATIO = 0.98

# Returns (tail factor r / (1 - r), r), (0.0, 0.0) when no tail is added; the tail of a flux is its last shell times the factor.
# The shells must also fall steadily over the last period (period = the n stop run): where they oscillate (zoom-whirl)
# a fit over a few shells says nothing about the next peak, and no tail is added.
function _geometric_n_tail_factor(shells::AbstractVector{<:Real}, e::Real, n_reached, period::Int = N_TAIL_FIT)
    (n_reached === nothing || e < N_TAIL_MIN_E || n_reached < N_TAIL_MIN_N) && return (0.0, 0.0)
    return _geometric_shell_tail(shells, period)
end

# Geometric tail factor r / (1 - r) of a sum whose last shells fall steadily (one sign, strictly decreasing in
# magnitude over max(N_TAIL_FIT, period + 1) shells), from a log-linear fit of the last N_TAIL_FIT shells; (0.0, 0.0)
# otherwise or when r >= N_TAIL_MAX_RATIO.
function _geometric_shell_tail(shells::AbstractVector{<:Real}, period::Int = N_TAIL_FIT)
    window = max(N_TAIL_FIT, period + 1)
    length(shells) < window && return (0.0, 0.0)
    recent = Float64.(shells[(end - window + 1):end])
    sgn = sign(recent[end])
    (sgn != 0 && all(x -> sign(x) == sgn, recent) && all(abs(recent[j + 1]) < abs(recent[j]) for j in 1:(window - 1))) || return (0.0, 0.0)
    last_shells = recent[(end - N_TAIL_FIT + 1):end]
    y = log.(abs.(last_shells))
    k = collect(0.0:(N_TAIL_FIT - 1))
    slope = sum((k .- sum(k) / N_TAIL_FIT) .* (y .- sum(y) / N_TAIL_FIT)) / sum(abs2, k .- sum(k) / N_TAIL_FIT)
    r = exp(slope)
    (0 < r < N_TAIL_MAX_RATIO) || return (0.0, 0.0)
    return (r / (1 - r), r)
end

# Defaults of the generic-orbit stopping parameters `beta` (tail budget) and `deep` (deep stop).
const DEFAULT_GENERIC_BETA = 1e-4
const DEFAULT_GENERIC_DEEP = 1e-5

# Walk ls, calling eval_l(l) (which returns the mode's energy flux), until minimum_consecutive consecutive modes have
# |flux| <= threshold, or (deep > 0) one mode has |flux| <= deep * threshold. Returns the last l walked, whether a
# stop rule ended the walk (false when it ran off the end of ls), the below-threshold streak at the end, and the
# first mode's |flux| (0.0 when ls is empty).
function _l_walk(eval_l, ls, threshold, minimum_consecutive::Int; deep = 0.0)
    below_count = 0
    walk_max = 0.0
    last_l = first(ls) - step(ls)
    stopped = false
    first_value = 0.0
    n_walked = 0
    trail = Float64[]
    for l in ls
        last_l = l
        value = abs(eval_l(l))
        n_walked += 1
        n_walked == 1 && (first_value = value)
        _push_trail!(trail, value)
        below_count = value <= threshold ? below_count + 1 : 0
        if below_count >= minimum_consecutive || _deep_stop(value, deep, threshold, walk_max)
            stopped = true
            break
        end
        walk_max = max(walk_max, value)
    end
    return (last_l = last_l, stopped = stopped, below_count = below_count, first_value = first_value, trail = trail)
end

# Keeps the last three |values| of a walk, oldest first.
_push_trail!(trail, value) = (push!(trail, abs(value)); length(trail) > 3 && popfirst!(trail); trail)

# Tail beyond the last of the trailing |values|: v_last r / (1 - r), with r the mean ratio over the trail. When the
# values do not decrease (rounding-level wiggles far below the cutoff, or a sum still rising at the cap), the tail
# is bounded by `cap` more terms of the largest trailing value; a single value counts as one more such term.
function _geometric_tail(trail, cap)
    n = length(trail)
    n == 0 && return 0.0
    last_value = trail[end]
    n == 1 && return last_value
    first_value = trail[1]
    r = first_value > 0 ? (last_value / first_value)^(1 / (n - 1)) : Inf
    return r < 1 ? last_value * r / (1 - r) : cap * maximum(trail)
end

# A sum that reached its cap is truncated when its extrapolated tail exceeds tol times the running boundary total.
_tail_hit(tail, tol, total) = tail > tol * abs(total)

# Tail of one mode (both signs of m) beyond an l walk that ran off the end of ls without its stop rule; 0 otherwise.
@inline _walk_tail(walk, ls) = !walk.stopped && !isempty(ls) && walk.last_l == last(ls) ? 2 * _geometric_tail(walk.trail, last(ls)) : 0.0

# Modes within two steps of a cap, for the tail-based cap-hit test of the eccentric and generic routes, which walk l,
# m and k in many nested loops: (n, m, k) => (l, E) of an l walk, (n, k, m) => m-shell E, (n, k) => k-shell E.
struct CapTracker
    lmax::Int
    kmax::Int
    l_walks::Dict{NTuple{3, Int}, Vector{Tuple{Int, Float64}}}
    m_shells::Dict{NTuple{3, Int}, Float64}
    k_shells::Dict{NTuple{2, Int}, Float64}
end
CapTracker(lmax::Int, kmax::Int) = CapTracker(lmax, kmax, Dict{NTuple{3, Int}, Vector{Tuple{Int, Float64}}}(),
    Dict{NTuple{3, Int}, Float64}(), Dict{NTuple{2, Int}, Float64}())

function _note_cap_mode!(t::CapTracker, l::Int, m::Int, n::Int, k::Int, E)
    isfinite(E) || return t
    l >= t.lmax - 2 && push!(get!(Vector{Tuple{Int, Float64}}, t.l_walks, (n, m, k)), (l, Float64(E)))
    abs(m) >= t.lmax - 2 && (t.m_shells[(n, k, m)] = get(t.m_shells, (n, k, m), 0.0) + 2E)
    t.kmax > 0 && abs(k) >= t.kmax - 2 && (t.k_shells[(n, k)] = get(t.k_shells, (n, k), 0.0) + 2E)
    return t
end

# Tail-based hit of the walks that reached a cap, sorted by their walk index (last three values, oldest first).
function _capped_walk_hit(groups, cap, tol, total, factor)
    for v in values(groups)
        any(e -> e[1] == cap, v) || continue
        sort!(v)
        trail = [factor * abs(e[2]) for e in v[max(1, end - 2):end]]
        _tail_hit(_geometric_tail(trail, cap), tol, total) && return true
    end
    return false
end

# (lmax_hit, kmax_hit) of one boundary from its tracker and its total.
function _cap_hits(t::CapTracker, tol, total)
    lhit = _capped_walk_hit(t.l_walks, t.lmax, tol, total, 2.0)
    if !lhit
        groups = Dict{NTuple{3, Int}, Vector{Tuple{Int, Float64}}}()
        for ((n, k, m), E) in t.m_shells
            push!(get!(Vector{Tuple{Int, Float64}}, groups, (n, k, sign(m))), (abs(m), E))
        end
        lhit = _capped_walk_hit(groups, t.lmax, tol, total, 1.0)
    end
    khit = false
    if t.kmax > 0
        groups = Dict{NTuple{3, Int}, Vector{Tuple{Int, Float64}}}()
        for ((n, k), E) in t.k_shells
            push!(get!(Vector{Tuple{Int, Float64}}, groups, (n, sign(k), 0)), (abs(k), E))
        end
        khit = _capped_walk_hit(groups, t.kmax, tol, total, 1.0)
    end
    return lhit, khit
end

# An n sum that reached nmax without its stop rule is truncated when the geometric tail of its last three shells matters.
_nmax_hit(n_reached, shells, nmax, tol, total) =
    n_reached === nothing && !isempty(shells) && _tail_hit(_geometric_tail(abs.(shells[max(1, end - 2):end]), nmax), tol, total)

# Modes more than this factor below the largest mode of their l walk are taken as zeros by symmetry
# (e.g. every other l at a = 0, |flux| ~ 1e-35 next to neighbours of ~1e-7), not as a decaying tail.
const _STRUCTURAL_ZERO_RATIO = 1e-20

# Deep stop of an l walk: one mode with |flux| <= deep * threshold ends the walk, unless it is a zero by
# symmetry relative to the largest earlier mode of the walk (`walk_max`; none yet means no deep stop), in
# which case only the consecutive-mode rule can stop the walk.
@inline _deep_stop(value, deep, threshold, walk_max) =
    deep > 0 && walk_max > 0 && value <= deep * threshold && value > _STRUCTURAL_ZERO_RATIO * walk_max

# l walk of one k != 0 branch; each mode is evaluated at most once, and every segment stops after minimum_consecutive
# modes below threshold. Infinity (from_bottom = false): up from the peak, then down from the peak towards lmin.
# Walks that search for the peak use search_threshold (the plain shell cutoff) without the deep stop, since a
# larger tail budget cannot bound an unvisited interior peak between the stopped bottom and descent walks.
# Horizon (from_bottom = true): Kerr horizon modes below the peak are not suppressed and often peak at lmin, so first
# walk up from lmin as for k = 0; if that stops below the peak, walk down from the peak to where it stopped, then up
# from the peak.
# Returns the extrapolated l tail (both signs of m) beyond lmax: the geometric tail of an ascending segment that ran
# off the end of lmax, lmax times the mode at lmax when the peak lies beyond lmax and that mode is above threshold,
# 0 otherwise.
function _branch_l_walk(eval_l, lmin::Int, lpeak::Int, lmax::Int, threshold, minimum_consecutive::Int, from_bottom::Bool; deep = 0.0, search_threshold = threshold)
    top = min(lpeak, lmax + 1) - 1
    if from_bottom
        bottom = _l_walk(eval_l, lmin:lmax, search_threshold, minimum_consecutive)
        bottom.stopped || return _walk_tail(bottom, lmin:lmax)
        lpeak > bottom.last_l || return 0.0
        down = _l_walk(eval_l, top:-1:(bottom.last_l + 1), search_threshold, minimum_consecutive)
        up = _l_walk(eval_l, lpeak:lmax, threshold, minimum_consecutive; deep=deep)
        return lpeak > lmax && down.first_value > search_threshold ? 2 * lmax * down.first_value : _walk_tail(up, lpeak:lmax)
    else
        up = _l_walk(eval_l, lpeak:lmax, threshold, minimum_consecutive; deep=deep)
        down = _l_walk(eval_l, top:-1:lmin, search_threshold, minimum_consecutive)
        return lpeak > lmax && down.first_value > search_threshold ? 2 * lmax * down.first_value : _walk_tail(up, lpeak:lmax)
    end
end

# Adaptive Levin modes of the eccentric sum must also change by less than LEVIN_ABS_FRACTION tol times the
# branch total between refinement levels. Without it the relative test (1e-6 relaxed to 1e-2 with depth) stopped
# the dominant mode of a shell early: (l, m, n) = (2, 2, 400) at (0.9, 8, 0.95, 1) was off by 1.7e-5.
const LEVIN_ABS_FRACTION = 0.01

# Automatic l and k caps (lmax = :auto, kmax = :auto).
const AUTO_CAP_MIN = 30
const AUTO_CAP_MAX = 150
const AUTO_CAP_MARGIN = 3
const AUTO_CAP_MAX_RAISES = 2

# Energy flux |F_inf| + |F_hor| of the (l, m = l) modes, l = 2..10, of the equatorial circular orbit at radius r, or
# nothing when that orbit cannot be evaluated.
function _auto_cap_probe_fluxes(a, r, x_probe)
    try
        scope, owned = ensure_public_submission(nothing, a, r, 0.0, x_probe)
        return _with_public_scope(scope, owned) do
            with_submission_operation(scope) do
                map(2:10) do l
                    mode = _circular_mode_flux_scoped(scope, a, r, l, l; x = x_probe)
                    abs(mode.infinity.energy_flux) + abs(mode.horizon.energy_flux)
                end
            end
        end
    catch
        return nothing
    end
end

# Upper estimate of the l needed for tol, from the l falloff of the equatorial circular orbit at the same a and
# radius r = p / (1 + e) (prograde for x >= 0, retrograde otherwise, prograde if that fails). A linear fit of
# log F_l over l = 6..10 gives the slope s < 0; the tail beyond L is F_10 exp(s (L - 10)) / (1 - e^s). The cap is
# the smallest L with that tail below tol * sum(F_2..F_10), plus AUTO_CAP_MARGIN, clamped to
# [AUTO_CAP_MIN, AUTO_CAP_MAX]. Returns (cap, status).
function auto_lmax_estimate(a, p, e, x; tol = 1e-8)
    r = p / (1 + e)
    probes = x >= 0 ? (1.0,) : (-1.0, 1.0)
    F = nothing
    for x_probe in probes
        F = _auto_cap_probe_fluxes(a, r, x_probe)
        F !== nothing && all(f -> isfinite(f) && f > 0, F) && break
        F = nothing
    end
    F === nothing && return (cap = AUTO_CAP_MIN, status = :probe_failed)
    ls = 6:10
    ys = log.(F[ls .- 1])
    lbar = sum(ls) / length(ls)
    ybar = sum(ys) / length(ys)
    slope = sum((ls .- lbar) .* (ys .- ybar)) / sum((ls .- lbar) .^ 2)
    (isfinite(slope) && slope < 0) || return (cap = AUTO_CAP_MAX, status = :no_falloff)
    budget = tol * sum(F)
    tail_scale = F[end] / (1 - exp(slope))
    L = tail_scale <= budget ? 10 : 10 + ceil(Int, log(budget / tail_scale) / slope)
    return (cap = clamp(L + AUTO_CAP_MARGIN, AUTO_CAP_MIN, AUTO_CAP_MAX), status = :ok)
end

# Initial caps: integers pass through; :auto takes the estimate (kmax = lmax estimate).
const AUTO_NMAX_MIN = 100
const AUTO_NMAX_MAX = 20000
const AUTO_NMAX_COEFF = 0.75

# Radial harmonic content: the azimuthal rate at periapsis r_p = p / (1 + e) (at the equatorial crossing for generic
# orbits) over the radial frequency. It grows without bound towards the separatrix (Omega_r -> 0).
function _radial_harmonic_scale(a, p, e, x)
    KG = kerr_geo_orbit(a, p, e, x)
    E = KG["Energy"]
    Lz = KG["AngularMomentum"]
    F = KG["Frequencies"]
    r = p / (1 + e)
    Δ = r^2 - 2r + a^2
    P = E * (r^2 + a^2) - a * Lz
    dφ = a * P / Δ - a * E + Lz
    dt = (r^2 + a^2) * P / Δ - a * (a * E - Lz)
    return abs(dφ / dt) / abs(F["ϒr"] / F["ϒt"])
end

# nmax = c_n (lmax + log(1/tol)) w_p / Omega_r, clamped to [AUTO_NMAX_MIN, AUTO_NMAX_MAX]. c_n = 0.75 gives at least
# 1.5 times the n reached by (0.9, 8, 0.99, 1) at tol 1e-6 (8436), (0.9, 8, 0.95, 1) (1005), (0.9, 6, 0.7, 1) at
# 1e-8 (124) and (0.5, p_sep + 0.05, 0.8, 1) (388).
function auto_nmax_estimate(a, p, e, x, tol, lmax)
    scale = try
        _radial_harmonic_scale(a, p, e, x)
    catch
        return (cap = AUTO_NMAX_MAX, status = :scale_failed)
    end
    isfinite(scale) || return (cap = AUTO_NMAX_MAX, status = :scale_failed)
    return (cap = clamp(ceil(Int, AUTO_NMAX_COEFF * (lmax + log(1 / tol)) * scale), AUTO_NMAX_MIN, AUTO_NMAX_MAX), status = :ok)
end

# Initial caps: integers pass through; :auto takes the estimates (kmax = lmax estimate; nmax from the radial
# harmonic scale with the lmax used). nmax = nothing for routes without an n sum.
function _auto_caps_start(a, p, e, x, tol, lmax, kmax, nmax = nothing)
    auto_l = lmax === :auto
    auto_k = kmax === :auto
    auto_n = nmax === :auto
    (auto_l || auto_k || auto_n) && !(tol > 0) && throw(ArgumentError("tol must be positive"))
    est = (cap = nothing, status = :fixed)
    if auto_l || auto_k
        est = auto_lmax_estimate(a, p, e, x; tol)
        if est.status == :probe_failed
            @warn "Automatic lmax/kmax: the equatorial circular probe at r = p / (1 + e) failed; using $(AUTO_CAP_MIN)." a p e x
        elseif est.status == :no_falloff
            @warn "Automatic lmax/kmax: the probe shows no l falloff by l = 10; using $(AUTO_CAP_MAX)." a p e x
        end
    end
    L = auto_l ? est.cap : lmax
    nest = (cap = nothing, status = :fixed)
    if auto_n
        nest = auto_nmax_estimate(a, p, e, x, tol, L)
        nest.status == :scale_failed && @warn "Automatic nmax: the radial harmonic scale could not be evaluated; using $(AUTO_NMAX_MAX)." a p e x
    end
    return (lmax = L, kmax = auto_k ? est.cap : kmax, nmax = auto_n ? nest.cap : nmax, auto_l, auto_k, auto_n,
        estimate = est.cap, n_estimate = nest.cap, status = est.status)
end

_result_lmax_hit(result, lmax) = get(result, :lmax_hit_inf, false) || get(result, :lmax_hit_hor, false) ||
    (haskey(result, :l_reached_inf) && !haskey(result, :lmax_hit_inf) && (something(result.l_reached_inf, 0) >= lmax || something(get(result, :l_reached_hor, 0), 0) >= lmax))
_result_kmax_hit(result) = get(result, :kmax_hit_inf, false) || get(result, :kmax_hit_hor, false)
_result_nmax_hit(result) = get(result, :nmax_hit_inf, false) || get(result, :nmax_hit_hor, false)

# Run run(lmax, kmax, nmax, warn) with the start caps; when a cap came from :auto and the sum hit it, raise that cap
# by half (lmax and kmax at most AUTO_CAP_MAX, nmax at most AUTO_NMAX_MAX) and sum again, at most
# AUTO_CAP_MAX_RAISES times. Adds the caps to the result.
function _with_auto_caps(run, caps)
    L, Kc, Nc = caps.lmax, caps.kmax, caps.nmax
    auto = caps.auto_l || caps.auto_k || caps.auto_n
    raises = 0
    result = run(L, Kc, Nc, !auto)
    while auto && raises < AUTO_CAP_MAX_RAISES
        raise_l = caps.auto_l && L < AUTO_CAP_MAX && _result_lmax_hit(result, L)
        raise_k = caps.auto_k && Kc isa Integer && Kc < AUTO_CAP_MAX && _result_kmax_hit(result)
        raise_n = caps.auto_n && Nc isa Integer && Nc < AUTO_NMAX_MAX && _result_nmax_hit(result)
        (raise_l || raise_k || raise_n) || break
        raise_l && (L = min(AUTO_CAP_MAX, ceil(Int, 1.5L)))
        raise_k && (Kc = min(AUTO_CAP_MAX, ceil(Int, 1.5Kc)))
        raise_n && (Nc = min(AUTO_NMAX_MAX, ceil(Int, 1.5Nc)))
        raises += 1
        result = run(L, Kc, Nc, false)
    end
    if auto
        _result_lmax_hit(result, L) && @warn "Mode summation hit lmax = $(L) with modes still above the cutoff after $(raises) automatic raise(s): the l sum is truncated." lmax = L
        _result_kmax_hit(result) && @warn "Mode summation ran to kmax = $(Kc) without converged k shells after $(raises) automatic raise(s): the k sum may be truncated." kmax = Kc
        _result_nmax_hit(result) && @warn "Mode summation ran to nmax = $(Nc) without converged n shells after $(raises) automatic raise(s): the n sum may be truncated." nmax = Nc
    end
    return merge(result, (lmax_used = L, kmax_used = Kc, nmax_used = Nc, lmax_auto = caps.auto_l, kmax_auto = caps.auto_k,
        nmax_auto = caps.auto_n, cap_estimate = caps.estimate, nmax_estimate = caps.n_estimate, cap_status = caps.status,
        cap_raises = raises))
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

    shells_inf = Float64[]
    shells_hor = Float64[]
    last_inf = (0.0, 0.0, 0.0)
    last_hor = (0.0, 0.0, 0.0)
    stopped_inf = false
    stopped_hor = false
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

        push!(shells_inf, shell_infinity_energy)
        last_inf = (shell_infinity_energy, shell_infinity_angular, shell_infinity_carter)
        threshold = tol * abs(total_infinity_energy)
        below_count = abs(shell_infinity_energy) <= threshold ? below_count + 1 : 0
        if below_count >= min_consecutive
            stopped_inf = true
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

        push!(shells_hor, shell_horizon_energy)
        last_hor = (shell_horizon_energy, shell_horizon_angular, shell_horizon_carter)
        threshold = tol * abs(total_horizon_energy)
        below_count = abs(shell_horizon_energy) <= threshold ? below_count + 1 : 0
        if below_count >= min_consecutive
            stopped_hor = true
            break
        end
    end

    # Geometric l tail beyond the stop (as the n tail of eccentric orbits): the l shells of a circular orbit fall off
    # geometrically, and min_consecutive shells below tol can still leave a few tol in the shells not evaluated.
    l_tail_factor_inf, l_tail_ratio_inf = stopped_inf ? _geometric_shell_tail(shells_inf, min_consecutive) : (0.0, 0.0)
    l_tail_factor_hor, l_tail_ratio_hor = stopped_hor ? _geometric_shell_tail(shells_hor, min_consecutive) : (0.0, 0.0)
    l_tail_inf = l_tail_factor_inf .* last_inf
    l_tail_hor = l_tail_factor_hor .* last_hor
    total_infinity_energy += l_tail_inf[1]
    total_infinity_angular += l_tail_inf[2]
    total_infinity_carter += l_tail_inf[3]
    total_horizon_energy += l_tail_hor[1]
    total_horizon_angular += l_tail_hor[2]
    total_horizon_carter += l_tail_hor[3]

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
        l_tail_infinity = (energy = l_tail_inf[1], angular_momentum = l_tail_inf[2], carter_constant = l_tail_inf[3], ratio = l_tail_ratio_inf),
        l_tail_horizon = (energy = l_tail_hor[1], angular_momentum = l_tail_hor[2], carter_constant = l_tail_hor[3], ratio = l_tail_ratio_hor),
    )
end

# Total circular-orbit fluxes: sum m = 1..l (doubled for -m) over l shells until min_consecutive shells are
# below tol times the running total, separately at infinity and at the horizon.
function circular_mode_summation(a, p; submission = nothing, x=1.0, tol = 1e-8, lmax = :auto, min_consecutive = 5, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing, Nmax = submission === nothing ? 2^14 : submission.Nmax, Kmax = submission === nothing ? 2^12 : submission.Kmax)
    caps = _auto_caps_start(a, p, 0.0, x, tol, lmax, nothing)
    scope, owned = ensure_public_submission(submission, a, p, 0.0, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        _with_auto_caps(caps) do L, _, _, _
            with_submission_operation(scope) do
                _circular_mode_summation_scoped(scope, a, p; x = x, tol = tol, lmax = L, min_consecutive = min_consecutive)
            end
        end
    end
end

_record_tag(x) = replace(replace(@sprintf("%.16g", x), "." => "p"), "-" => "m")

function _eccentric_mode_summation_scoped(scope, a, p, e; x=1.0, N = 64, N0 = N, Nmax = 2^14, tol = 1e-8, lmax = 30, nmax = 500, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, tail_levin = :auto, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_local_n::Int = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, levin_max_depth::Int = 8)
    rules = orbit_rules(a, p, e, x)
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
    # The n shells repeat about every Ups_phi / Ups_r shells (lowering m by one and raising n by Ups_phi / Ups_r gives the
    # same omega). Near the separatrix (zoom-whirl) that is ~ 30 shells, with deep dips between the peaks: at
    # (a, p, e, x) = (0.254, 8.61, 0.821, -1) the walk stopped in the dip at n = 250 and missed the next peaks (1e-6 of
    # E_inf at tol = 1e-8). The n stop needs a run of below-threshold shells spanning one such period.
    n_stop_run = _n_stop_run(minimum_consecutive, radial_over_azimuthal)
    F_orbit = KG["Frequencies"]
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
    last_shell_inf = (0.0, 0.0, 0.0)
    last_shell_hor = (0.0, 0.0, 0.0)
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
    # Absolute energy tolerance of an adaptive Levin mode: LEVIN_ABS_FRACTION tol times the branch total so far.
    branch_levin_abs_inf = Ref(0.0)
    branch_levin_abs_hor = Ref(0.0)
    levin_floor_tracks_mode_floor = levin_mode_abs_floor == mode_abs_floor
    branch_mode_floor = s -> s == -2 ? branch_mode_floor_inf[] : branch_mode_floor_hor[]
    branch_levin_floor = s -> levin_floor_tracks_mode_floor ? branch_mode_floor(s) : Float64(levin_mode_abs_floor)
    ci_kwargs = fast ? (Nmax = Nmax, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling) : (Nmax = Nmax, zero_low_flux = zero_low_flux, threaded_sampling = false)
    tail_levin_infinity_active = Ref(false)
    tail_levin_horizon_active = Ref(false)
    tail_levin_infinity_start_n = Ref{Union{Nothing, Int}}(nothing)
    tail_levin_horizon_start_n = Ref{Union{Nothing, Int}}(nothing)
    # tail_levin = true: every mode from the first |n| >= nmin on. :auto: each mode whose trapezoid grid,
    # estimated from its radial phase rate, passes the cost rule of _ecc_levin_auto.
    use_eccentric_tail_levin = (s, m, n) -> begin
        branch_policy = s == -2 ? tail_policy_inf : tail_policy_hor
        _tail_enabled(branch_policy) && fast || return false
        active = s == -2 ? tail_levin_infinity_active : tail_levin_horizon_active
        start_n = s == -2 ? tail_levin_infinity_start_n : tail_levin_horizon_start_n
        branch_policy === true && active[] && return true
        _use_levin_tail(n; nmin = levin_nmin) || return false
        if branch_policy === true || _ecc_levin_auto(
                _mode_radial_trapezoid_n(KG_sample, m, n, 0, N, Nmax), n,
                adaptive_local_n, levin_max_depth, levin_nmin)
            active[] = true
            start_n[] === nothing && (start_n[] = n)
            return true
        end
        return false
    end
    cap_trackers = (CapTracker(lmax, 0), CapTracker(lmax, 0))
    failed_modes = [0, 0]
    near_static_zero = [0, 0]
    ecc_eval = (s, l, m, n, max_flux, tail_values) -> begin
        use_levin = use_eccentric_tail_levin(s, m, n)
        floor = MODE_REFINE_FLOOR_FRACTION * branch_mode_floor(s)
        ω_mode = (m * F_orbit["ϒϕ"] + n * F_orbit["ϒr"]) / F_orbit["ϒt"]
        mode = if _negligible_low_frequency_mode(s, l, ω_mode, F_orbit["ϒϕ"] / F_orbit["ϒt"], a, p, e)
            res = ConvolutionIntegrals._zero_radiative_mode(ω_mode, nothing; reason = "negligible_low_frequency")
            res["NearStaticOverflow"] = true
            res
        else
        let; ConvolutionIntegrals._with_coverage_retry(s, l, m, n, 0; omega = ω_mode, Ωφ = F_orbit["ϒϕ"] / F_orbit["ϒt"], rmax = p / (1 - e), negligible_if_failed = _negligible_low_frequency_mode(s, l, ω_mode, F_orbit["ϒϕ"] / F_orbit["ϒt"], a, p, e; bound = ConvolutionIntegrals.NEGLIGIBLE_FAILED_MODE_BOUND)) do; if use_levin
            with_submission_mode(scope; mode = (s, l, m, n, 0)) do mode_cache
                ConvolutionIntegrals.convolution_integral_eccentric_adaptive_levin_isem(KG_sample, s, l, m, n; local_n = adaptive_local_n, max_depth = levin_max_depth, tol = tol, mode_abs_floor = _tail_levin_mode_floor(floor, branch_levin_floor(s)), energy_abs_tol = s == -2 ? branch_levin_abs_inf[] : branch_levin_abs_hor[], zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, cache = mode_cache)
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
        end; end; end
        end
        get(mode, "FailedMode", false) && (failed_modes[s == -2 ? 1 : 2] += 1)
        get(mode, "NearStaticOverflow", false) && (near_static_zero[s == -2 ? 1 : 2] += 1)
        ConvolutionIntegrals._check_finite_mode_flux(mode, s, l, m, n, 0)
        _note_cap_mode!(cap_trackers[s == -2 ? 1 : 2], l, m, n, 0, mode["EnergyFlux"])
        mode
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

    n0_reach_inf = 0
    for m in 1:lmax
        n0_reach_inf = m
        lmin = max(2, abs(m))
        shell_m_infinity_energy = 0.0
        shell_m_infinity_angular = 0.0
        shell_m_infinity_carter = 0.0
        below_count_l = 0

        for l in lmin:lmax
            mode = let; with_submission_mode(scope; mode = (-2, l, m, 0, 0)) do mode_cache
                ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, -2, l, m, 0, N; tol = tol, sample_tol = sample_tol, mode_abs_floor = _adaptive_branch_mode_floor(shell_n_infinity_energy + shell_m_infinity_energy, tol, mode_abs_floor), ci_kwargs..., cache = mode_cache)
            end; end
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
    n0_reach_hor = 0
    for m in 1:lmax
        n0_reach_hor = m
        lmin = max(2, abs(m))
        shell_m_horizon_energy = 0.0
        shell_m_horizon_angular = 0.0
        shell_m_horizon_carter = 0.0
        below_count_l = 0

        for l in lmin:lmax
            mode = let; with_submission_mode(scope; mode = (2, l, m, 0, 0)) do mode_cache
                ConvolutionIntegrals.convolution_integral_eccentric_trapezoidal_isem(KG_sample, 2, l, m, 0, N; tol = tol, sample_tol = sample_tol, mode_abs_floor = _adaptive_branch_mode_floor(shell_n_horizon_energy + shell_m_horizon_energy, tol, mode_abs_floor), ci_kwargs..., cache = mode_cache)
            end; end
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
    branch_levin_abs_inf[] = LEVIN_ABS_FRACTION * tol * abs(total_infinity_energy)
    branch_levin_abs_hor[] = LEVIN_ABS_FRACTION * tol * abs(total_horizon_energy)

    prev_peak_m_pos = 1
    prev_peak_m_neg = 1
    # |m| reached by the previous n shell's walks; an m shell counts toward the stop only past MWALK_REACH_FRACTION of it
    prev_reach_pos = n0_reach_inf
    prev_reach_neg = n0_reach_inf
    for n in 1:nmax
        m_omega_zero = -n * radial_over_azimuthal
        omega_zero_gate = !rules.high_e || abs(m_omega_zero) < lmax
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
        reach_pos = 0
        reach_neg = 0
        last_big_pos = 0
        last_big_neg = 0
        for m in 1:lmax
            reach_pos = m
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
            abs(shell_m_infinity_energy) > layer_threshold && (last_big_pos = abs(m))
            # Near the ISSO Ups_r << Ups_phi and the n shells of even a nearly circular orbit peak at large m (at
            # (a, p, e) = (0.932, 2.095, 0.009) the n = 6 shell sits at m = l = 15..27 while m <= 2 are tiny, and the walk
            # dropped 4e-6 of E_inf there); a shell counts toward the stop only past the peak m of the previous n shell.
            below_count_m = !(omega_zero_gate && _before_omega_zero(m, m_omega_zero)) && (abs(shell_m_infinity_energy) <= layer_threshold) && (abs(shell_n_infinity_energy) > min(0.1, exp(1 - 1 / e)) * abs(shell_n_infinity_energy_last)) &&
                (m >= max(prev_peak_m_pos, floor(Int, MWALK_REACH_FRACTION * prev_reach_pos))) ? below_count_m + 1 : 0
            if below_count_m >= minimum_consecutive
                break
            end
        end

        below_count_m = 0
        Max_flux_neg = 0.0
        Max_m_neg = 1
        for m in -1:-1:-lmax
            reach_neg = -m
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
                if abs(mode["EnergyFlux"]) > Max_flux_neg
                    Max_flux_neg = abs(mode["EnergyFlux"])
                    Max_m_neg = abs(m)
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
            abs(shell_m_infinity_energy) > layer_threshold && (last_big_neg = abs(m))
            # As in the positive branch, the streak only counts once the n shell holds a fair part of the previous one:
            # for retrograde orbits (x = -1) the negative m carry the flux, and at high n it sits at |m| ~ 5..10, past
            # m = -1, -2, which are tiny (at (a, p, e) = (0.525, 9.52, 0.531) the walk stopped there for n = 22, 25,
            # 27.., dropping 7.6e-4 of E_inf).
            below_count_m = !(omega_zero_gate && _before_omega_zero(m, m_omega_zero)) && (abs(shell_m_infinity_energy) <= layer_threshold) && (abs(m) >= max(Max_m, prev_peak_m_neg, floor(Int, MWALK_REACH_FRACTION * prev_reach_neg))) &&
                (abs(shell_n_infinity_energy) > min(0.1, exp(1 - 1 / e)) * abs(shell_n_infinity_energy_last)) ? below_count_m + 1 : 0
            if below_count_m >= minimum_consecutive
                break
            end
        end
        shell_n_infinity_energy_last = shell_n_infinity_energy
        prev_peak_m_pos = Max_m
        prev_peak_m_neg = Max_m_neg
        prev_reach_pos = rules.high_e ? min(lmax, last_big_pos + minimum_consecutive) : reach_pos
        prev_reach_neg = rules.high_e ? min(lmax, last_big_neg + minimum_consecutive) : reach_neg

        push!(Energy_flux_inf, shell_n_infinity_energy)
        push!(n_list_inf, n)
        last_shell_inf = (shell_n_infinity_energy, shell_n_infinity_angular, shell_n_infinity_carter)
        total_infinity_energy += shell_n_infinity_energy
        total_infinity_angular += shell_n_infinity_angular
        total_infinity_carter += shell_n_infinity_carter
        branch_mode_floor_inf[] = _adaptive_branch_mode_floor(total_infinity_energy, tol, mode_abs_floor)
        branch_levin_abs_inf[] = LEVIN_ABS_FRACTION * tol * abs(total_infinity_energy)
        n_threshold = _mode_cutoff_threshold(tol * abs(total_infinity_energy), branch_mode_floor_inf[])
        below_count_n = abs(shell_n_infinity_energy) <= n_threshold ? below_count_n + 1 : 0
        if below_count_n >= n_stop_run
            n_reached_inf = n
            break
        end
    end

    below_count_n = 0
    prev_peak_m_pos = 1
    prev_peak_m_neg = 1
    prev_reach_pos = n0_reach_hor
    prev_reach_neg = n0_reach_hor
    for n in 1:nmax
        m_omega_zero = -n * radial_over_azimuthal
        omega_zero_gate = !rules.high_e || abs(m_omega_zero) < lmax
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
        reach_pos = 0
        reach_neg = 0
        last_big_pos = 0
        last_big_neg = 0
        for m in 1:lmax
            reach_pos = m
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

            # Horizon shells are negative in the superradiant regime, so the shell ratio compares magnitudes.
            layer_threshold = _mode_cutoff_threshold(tol * abs(shell_n_horizon_energy_last), branch_mode_floor_hor[])
            abs(shell_m_horizon_energy) > layer_threshold && (last_big_pos = abs(m))
            # Near the ISSO Ups_r << Ups_phi and the n shells of even a nearly circular orbit peak at large m (at
            # (a, p, e) = (0.932, 2.095, 0.009) the n = 6 shell sits at m = l = 15..27 while m <= 2 are tiny, and the walk
            # dropped 4e-6 of E_inf there); a shell counts toward the stop only past the peak m of the previous n shell.
            below_count_m = !(omega_zero_gate && _before_omega_zero(m, m_omega_zero)) && (abs(shell_m_horizon_energy) <= layer_threshold) && (abs(shell_n_horizon_energy) > min(0.1, exp(1 - 1 / e)) * abs(shell_n_horizon_energy_last)) &&
                (m >= max(prev_peak_m_pos, floor(Int, MWALK_REACH_FRACTION * prev_reach_pos))) ? below_count_m + 1 : 0
            if below_count_m >= minimum_consecutive

                break
            end
        end

        below_count_m = 0
        Max_flux_neg = 0.0
        Max_m_neg = 1
        for m in -1:-1:-lmax
            reach_neg = -m
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
                if abs(mode["EnergyFlux"]) > Max_flux_neg
                    Max_flux_neg = abs(mode["EnergyFlux"])
                    Max_m_neg = abs(m)
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
            abs(shell_m_horizon_energy) > layer_threshold && (last_big_neg = abs(m))
            # As in the positive branch, the streak only counts once the n shell holds a fair part of the previous one:
            # for retrograde orbits (x = -1) the negative m carry the flux, and at high n it sits at |m| ~ 5..10, past
            # m = -1, -2, which are tiny (at (a, p, e) = (0.525, 9.52, 0.531) the walk stopped there for n = 22, 25,
            # 27.., dropping 7.6e-4 of E_inf).
            below_count_m = !(omega_zero_gate && _before_omega_zero(m, m_omega_zero)) && (abs(shell_m_horizon_energy) <= layer_threshold) && (abs(m) >= max(Max_m, prev_peak_m_neg, floor(Int, MWALK_REACH_FRACTION * prev_reach_neg))) &&
                (abs(shell_n_horizon_energy) > min(0.1, exp(1 - 1 / e)) * abs(shell_n_horizon_energy_last)) ? below_count_m + 1 : 0
            if below_count_m >= minimum_consecutive
                break
            end
        end
        shell_n_horizon_energy_last = shell_n_horizon_energy
        prev_peak_m_pos = Max_m
        prev_peak_m_neg = Max_m_neg
        prev_reach_pos = rules.high_e ? min(lmax, last_big_pos + minimum_consecutive) : reach_pos
        prev_reach_neg = rules.high_e ? min(lmax, last_big_neg + minimum_consecutive) : reach_neg

        push!(Energy_flux_hor, shell_n_horizon_energy)
        push!(n_list_hor, n)
        last_shell_hor = (shell_n_horizon_energy, shell_n_horizon_angular, shell_n_horizon_carter)
        total_horizon_energy += shell_n_horizon_energy
        total_horizon_angular += shell_n_horizon_angular
        total_horizon_carter += shell_n_horizon_carter
        branch_mode_floor_hor[] = _adaptive_branch_mode_floor(total_horizon_energy, tol, mode_abs_floor)
        branch_levin_abs_hor[] = LEVIN_ABS_FRACTION * tol * abs(total_horizon_energy)
        n_threshold = _mode_cutoff_threshold(tol * abs(total_horizon_energy), branch_mode_floor_hor[])
        below_count_n = abs(shell_n_horizon_energy) <= n_threshold ? below_count_n + 1 : 0
        if below_count_n >= n_stop_run
            n_reached_hor = n
            break
        end
    end

    tail_factor_inf, tail_ratio_inf = _geometric_n_tail_factor(Energy_flux_inf, e, n_reached_inf, n_stop_run)
    tail_factor_hor, tail_ratio_hor = _geometric_n_tail_factor(Energy_flux_hor, e, n_reached_hor, n_stop_run)
    n_tail_inf = tail_factor_inf .* last_shell_inf
    n_tail_hor = tail_factor_hor .* last_shell_hor
    total_infinity_energy += n_tail_inf[1]
    total_infinity_angular += n_tail_inf[2]
    total_infinity_carter += n_tail_inf[3]
    total_horizon_energy += n_tail_hor[1]
    total_horizon_angular += n_tail_hor[2]
    total_horizon_carter += n_tail_hor[3]

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
        record_h5["meta/n_tail_infinity"] = collect(n_tail_inf)
        record_h5["meta/n_tail_horizon"] = collect(n_tail_hor)
        record_h5["meta/n_tail_ratio_infinity"] = tail_ratio_inf
        record_h5["meta/n_tail_ratio_horizon"] = tail_ratio_hor
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
        lmax_hit_inf = _cap_hits(cap_trackers[1], tol, total_infinity_energy)[1],
        lmax_hit_hor = _cap_hits(cap_trackers[2], tol, total_horizon_energy)[1],
        nmax_hit_inf = _nmax_hit(n_reached_inf, Energy_flux_inf, nmax, tol, total_infinity_energy),
        nmax_hit_hor = _nmax_hit(n_reached_hor, Energy_flux_hor, nmax, tol, total_horizon_energy),
        failed_modes_inf = failed_modes[1],
        failed_modes_hor = failed_modes[2],
        near_static_zero_inf = near_static_zero[1],
        near_static_zero_hor = near_static_zero[2],
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
        n_tail_infinity = (energy = n_tail_inf[1], angular_momentum = n_tail_inf[2], carter_constant = n_tail_inf[3], ratio = tail_ratio_inf),
        n_tail_horizon = (energy = n_tail_hor[1], angular_momentum = n_tail_hor[2], carter_constant = n_tail_hor[3], ratio = tail_ratio_hor),
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
function eccentric_mode_summation(a, p, e; submission = nothing, x=1.0, N = 64, N0 = N, Nmax = submission === nothing ? 2^14 : submission.Nmax, tol = 1e-8, lmax = :auto, nmax = :auto, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, tail_levin = :auto, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_local_n::Int = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, levin_max_depth::Int = 8, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing, Kmax = submission === nothing ? 2^12 : submission.Kmax)
    caps = _auto_caps_start(a, p, e, x, tol, lmax, nothing, nmax)
    scope, owned = ensure_public_submission(submission, a, p, e, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        _with_auto_caps(caps) do L, Kc, Nc, warn
            with_submission_operation(scope) do
                _eccentric_mode_summation_scoped(scope, a, p, e; x = x, N = N, N0 = N0, Nmax = Nmax, tol = tol, lmax = L, nmax = Nc, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, tail_levin = tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_local_n = levin_local_n, levin_max_depth = levin_max_depth)
            end
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
        mode = circular_mode_flux(a, p, l, m; x = x, submission = scope, lifecycle_observer = scope.lifecycle_observer, geometry_observer = scope.geometry_observer)
        return (l = l, m = m, k = k, infinity = mode.infinity, horizon = mode.horizon)
    end

    KG = submission_orbit!(scope)
    KG_sample = submission_master!(scope)
    flux_inf = let; with_submission_mode(scope; mode = (-2, l, m, 0, k)) do mode_cache
        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; Kmax = Kmax, sample_tol = sample_tol, cache = mode_cache)
    end; end
    flux_hor = let; with_submission_mode(scope; mode = (2, l, m, 0, k)) do mode_cache
        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; Kmax = Kmax, sample_tol = sample_tol, cache = mode_cache)
    end; end
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

function _inclined_mode_summation_scoped(scope, a, p, x; K = 16, K0 = K, Kmax = 2^12, tol = 1e-8, lmax = 30, kmax = 20, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, warn_truncation::Bool = true)
    lmax < 2 && throw(ArgumentError("lmax must be at least 2"))
    minimum_consecutive < 1 && throw(ArgumentError("minimum_consecutive must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))

    if x == 1.0 || x == -1.0
        return circular_mode_summation(a, p; x = x, tol = tol, lmax = lmax, min_consecutive = minimum_consecutive, submission = scope, lifecycle_observer = scope.lifecycle_observer, geometry_observer = scope.geometry_observer)
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
    # Highest l evaluated per boundary, and whether an l (or m) walk hit lmax with modes still above the cutoff.
    l_reached_inf = 0
    l_reached_hor = 0
    lmax_hit_inf = false
    lmax_hit_hor = false
    # Near-static modes recorded as zero after a matching overflow (see ConvolutionIntegrals.NEAR_STATIC_RATIO).
    near_static_zero_inf = 0
    near_static_zero_hor = 0
    failed_modes_inf = 0
    failed_modes_hor = 0

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

        # The k = 0 m walks stop per boundary once minimum_consecutive m shells past the largest one are below tol
        # times the running k = 0 total, so modes far below the cutoff (which can underflow to NaN, l = m ~ 20 at
        # p ~ 2e5) are never evaluated.
        done_inf = false
        done_hor = false
        below_m_inf = 0
        below_m_hor = 0
        peak_m_inf = 0.0
        peak_m_hor = 0.0
        m_trail_inf = Float64[]
        m_trail_hor = Float64[]
        for m in 1:lmax
            lmin = max(2, abs(m))
            shell_m_infinity_energy = 0.0
            shell_m_infinity_angular = 0.0
            shell_m_infinity_carter = 0.0
            shell_m_horizon_energy = 0.0
            shell_m_horizon_angular = 0.0
            shell_m_horizon_carter = 0.0
            below_count_l = 0
            if !done_inf
                l_trail = Float64[]

                for l in lmin:lmax
                    mode_inf = let; with_submission_mode(scope; mode = (-2, l, m, 0, 0)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, 0, K; tol = tol, sample_tol = sample_tol, mode_abs_floor = _adaptive_branch_mode_floor(shell_k_infinity_energy + shell_m_infinity_energy, tol, mode_abs_floor), ci_kwargs..., cache = mode_cache)
                    end; end
                    _record_current_mode_inclined!(record_h5, record, "I", 0, m, l, mode_index, mode_inf, K)
                    get(mode_inf, "NearStaticOverflow", false) && (near_static_zero_inf += 1)
                    get(mode_inf, "FailedMode", false) && (failed_modes_inf += 1)
                    shell_m_infinity_energy += 2 * mode_inf["EnergyFlux"]
                    shell_m_infinity_angular += 2 * mode_inf["AngularMomentumFlux"]
                    shell_m_infinity_carter += 2 * mode_inf["CarterConstantFlux"]
                    total_modes += 2
                    l_reached_inf = max(l_reached_inf, l)
                    _push_trail!(l_trail, mode_inf["EnergyFlux"])
                    layer_threshold = _mode_cutoff_threshold(tol * abs(shell_m_infinity_energy), _adaptive_branch_mode_floor(shell_k_infinity_energy + shell_m_infinity_energy, tol, mode_abs_floor))
                    below_count_l = abs(mode_inf["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end
                # A walk that ran to lmax is truncated when its extrapolated tail matters for the boundary total; the cutoff
                # inside the walk is relative to this m shell, so short walks at large m always end above it.
                below_count_l < minimum_consecutive && _tail_hit(2 * _geometric_tail(l_trail, lmax), tol, shell_k_infinity_energy + shell_m_infinity_energy) && (lmax_hit_inf = true)

                shell_k_infinity_energy += shell_m_infinity_energy
                shell_k_infinity_angular += shell_m_infinity_angular
                shell_k_infinity_carter += shell_m_infinity_carter
                _push_trail!(m_trail_inf, shell_m_infinity_energy)
                peak_m_inf = max(peak_m_inf, abs(shell_m_infinity_energy))
                below_m_inf = abs(shell_m_infinity_energy) < peak_m_inf && abs(shell_m_infinity_energy) <= tol * abs(shell_k_infinity_energy) ? below_m_inf + 1 : 0
                done_inf = below_m_inf >= minimum_consecutive
            end

            if !done_hor
                below_count_l = 0
                l_trail = Float64[]
                for l in lmin:lmax
                    mode_hor = let; with_submission_mode(scope; mode = (2, l, m, 0, 0)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, 0, K; tol = tol, sample_tol = sample_tol, mode_abs_floor = _adaptive_branch_mode_floor(shell_k_horizon_energy + shell_m_horizon_energy, tol, mode_abs_floor), ci_kwargs..., cache = mode_cache)
                    end; end
                    _record_current_mode_inclined!(record_h5, record, "H", 0, m, l, mode_index, mode_hor, K)
                    get(mode_hor, "NearStaticOverflow", false) && (near_static_zero_hor += 1)
                    get(mode_hor, "FailedMode", false) && (failed_modes_hor += 1)
                    shell_m_horizon_energy += 2 * mode_hor["EnergyFlux"]
                    shell_m_horizon_angular += 2 * mode_hor["AngularMomentumFlux"]
                    shell_m_horizon_carter += 2 * mode_hor["CarterConstantFlux"]
                    total_modes += 2
                    l_reached_hor = max(l_reached_hor, l)
                    _push_trail!(l_trail, mode_hor["EnergyFlux"])
                    layer_threshold = _mode_cutoff_threshold(tol * abs(shell_m_horizon_energy), _adaptive_branch_mode_floor(shell_k_horizon_energy + shell_m_horizon_energy, tol, mode_abs_floor))
                    below_count_l = abs(mode_hor["EnergyFlux"]) <= layer_threshold ? below_count_l + 1 : 0
                    if below_count_l >= minimum_consecutive
                        break
                    end
                end
                # A walk that ran to lmax is truncated when its extrapolated tail matters for the boundary total; the cutoff
                # inside the walk is relative to this m shell, so short walks at large m always end above it.
                below_count_l < minimum_consecutive && _tail_hit(2 * _geometric_tail(l_trail, lmax), tol, shell_k_horizon_energy + shell_m_horizon_energy) && (lmax_hit_hor = true)

                shell_k_horizon_energy += shell_m_horizon_energy
                shell_k_horizon_angular += shell_m_horizon_angular
                shell_k_horizon_carter += shell_m_horizon_carter
                _push_trail!(m_trail_hor, shell_m_horizon_energy)
                peak_m_hor = max(peak_m_hor, abs(shell_m_horizon_energy))
                below_m_hor = abs(shell_m_horizon_energy) < peak_m_hor && abs(shell_m_horizon_energy) <= tol * abs(shell_k_horizon_energy) ? below_m_hor + 1 : 0
                done_hor = below_m_hor >= minimum_consecutive
            end
            done_inf && done_hor && break
        end

        done_inf || _tail_hit(_geometric_tail(m_trail_inf, lmax), tol, shell_k_infinity_energy) && (lmax_hit_inf = true)
        done_hor || _tail_hit(_geometric_tail(m_trail_hor, lmax), tol, shell_k_horizon_energy) && (lmax_hit_hor = true)

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
            l_tail = _branch_l_walk(2, _branch_peak_l(2, 0, k, x), lmax, layer_threshold, minimum_consecutive, false) do l
                mode = let; with_submission_mode(scope; mode = (-2, l, 0, 0, k)) do mode_cache
                    ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, 0, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, mode_abs_floor = MODE_REFINE_FLOOR_FRACTION * branch_mode_floor_inf[], ci_kwargs..., cache = mode_cache)
                end; end
                _record_current_mode_inclined!(record_h5, record, "I", k, 0, l, mode_index, mode, K)
                get(mode, "NearStaticOverflow", false) && (near_static_zero_inf += 1)
                get(mode, "FailedMode", false) && (failed_modes_inf += 1)
                l_reached_inf = max(l_reached_inf, l)
                shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                mode["EnergyFlux"]
            end
            lmax_hit_inf = lmax_hit_inf || _tail_hit(l_tail, tol, total_infinity_energy + shell_k_infinity_energy + shell_m_infinity_energy)

            shell_k_infinity_energy += shell_m_infinity_energy
            shell_k_infinity_angular += shell_m_infinity_angular
            shell_k_infinity_carter += shell_m_infinity_carter

            Max_flux = 0.0
            Max_m = 1
            m_trail = Float64[]
            # an m shell counts towards the stop only past the largest m shell of this k walk
            peak_m_infinity = 0.0
            for m in 1:lmax
                lmin = max(2, abs(m))
                shell_m_infinity_energy = 0.0
                shell_m_infinity_angular = 0.0
                shell_m_infinity_carter = 0.0
                below_count_l = 0

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
                l_tail = _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, layer_threshold, minimum_consecutive, false) do l
                    mode = let; with_submission_mode(scope; mode = (-2, l, m, 0, k)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, mode_abs_floor = MODE_REFINE_FLOOR_FRACTION * branch_mode_floor_inf[], ci_kwargs..., cache = mode_cache)
                    end; end
                    _record_current_mode_inclined!(record_h5, record, "I", k, m, l, mode_index, mode, K)
                    get(mode, "NearStaticOverflow", false) && (near_static_zero_inf += 1)
                    get(mode, "FailedMode", false) && (failed_modes_inf += 1)
                    l_reached_inf = max(l_reached_inf, l)
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
                lmax_hit_inf = lmax_hit_inf || _tail_hit(l_tail, tol, total_infinity_energy + shell_k_infinity_energy + shell_m_infinity_energy)

                shell_k_infinity_energy += shell_m_infinity_energy
                shell_k_infinity_angular += shell_m_infinity_angular
                shell_k_infinity_carter += shell_m_infinity_carter
                _push_trail!(m_trail, shell_m_infinity_energy)

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
                below_count_m = !_before_m_peak(m, k, x) && (abs(shell_m_infinity_energy) <= layer_threshold) && (abs(shell_m_infinity_energy) < peak_m_infinity) ? below_count_m + 1 : 0
                peak_m_infinity = max(peak_m_infinity, abs(shell_m_infinity_energy))
                if below_count_m >= minimum_consecutive
                    break
                end
            end
            below_count_m < minimum_consecutive && _tail_hit(_geometric_tail(m_trail, lmax), tol, total_infinity_energy + shell_k_infinity_energy) && (lmax_hit_inf = true)

            below_count_m = 0
            m_trail = Float64[]
            for m in -1:-1:-lmax
                lmin = max(2, abs(m))
                shell_m_infinity_energy = 0.0
                shell_m_infinity_angular = 0.0
                shell_m_infinity_carter = 0.0
                below_count_l = 0

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
                l_tail = _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, layer_threshold, minimum_consecutive, false) do l
                    mode = let; with_submission_mode(scope; mode = (-2, l, m, 0, k)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, -2, l, m, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_inf, mode_abs_floor = MODE_REFINE_FLOOR_FRACTION * branch_mode_floor_inf[], ci_kwargs..., cache = mode_cache)
                    end; end
                    _record_current_mode_inclined!(record_h5, record, "I", k, m, l, mode_index, mode, K)
                    get(mode, "NearStaticOverflow", false) && (near_static_zero_inf += 1)
                    get(mode, "FailedMode", false) && (failed_modes_inf += 1)
                    l_reached_inf = max(l_reached_inf, l)
                    shell_m_infinity_energy += 2 * mode["EnergyFlux"]
                    shell_m_infinity_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_infinity_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    mode["EnergyFlux"]
                end
                lmax_hit_inf = lmax_hit_inf || _tail_hit(l_tail, tol, total_infinity_energy + shell_k_infinity_energy + shell_m_infinity_energy)

                shell_k_infinity_energy += shell_m_infinity_energy
                shell_k_infinity_angular += shell_m_infinity_angular
                shell_k_infinity_carter += shell_m_infinity_carter
                _push_trail!(m_trail, shell_m_infinity_energy)

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_infinity_energy_last), branch_mode_floor_inf[])
                below_count_m = !_before_m_peak(m, k, x) && (abs(shell_m_infinity_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
                if below_count_m >= minimum_consecutive
                    shell_k_infinity_energy_last = shell_k_infinity_energy
                    break
                end
            end
            below_count_m < minimum_consecutive && _tail_hit(_geometric_tail(m_trail, lmax), tol, total_infinity_energy + shell_k_infinity_energy) && (lmax_hit_inf = true)

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
            l_tail = _branch_l_walk(2, _branch_peak_l(2, 0, k, x), lmax, layer_threshold, minimum_consecutive, true) do l
                mode = let; with_submission_mode(scope; mode = (2, l, 0, 0, k)) do mode_cache
                    ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, 0, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, mode_abs_floor = MODE_REFINE_FLOOR_FRACTION * branch_mode_floor_hor[], ci_kwargs..., cache = mode_cache)
                end; end
                _record_current_mode_inclined!(record_h5, record, "H", k, 0, l, mode_index, mode, K)
                get(mode, "NearStaticOverflow", false) && (near_static_zero_hor += 1)
                get(mode, "FailedMode", false) && (failed_modes_hor += 1)
                l_reached_hor = max(l_reached_hor, l)
                shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                total_modes += 2
                mode["EnergyFlux"]
            end
            lmax_hit_hor = lmax_hit_hor || _tail_hit(l_tail, tol, total_horizon_energy + shell_k_horizon_energy + shell_m_horizon_energy)

            shell_k_horizon_energy += shell_m_horizon_energy
            shell_k_horizon_angular += shell_m_horizon_angular
            shell_k_horizon_carter += shell_m_horizon_carter

            Max_flux = 0.0
            Max_m = 1
            m_trail = Float64[]
            # an m shell counts towards the stop only past the largest m shell of this k walk
            peak_m_horizon = 0.0
            for m in 1:lmax
                lmin = max(2, abs(m))
                shell_m_horizon_energy = 0.0
                shell_m_horizon_angular = 0.0
                shell_m_horizon_carter = 0.0
                below_count_l = 0

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
                l_tail = _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, layer_threshold, minimum_consecutive, true) do l
                    mode = let; with_submission_mode(scope; mode = (2, l, m, 0, k)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, mode_abs_floor = MODE_REFINE_FLOOR_FRACTION * branch_mode_floor_hor[], ci_kwargs..., cache = mode_cache)
                    end; end
                    _record_current_mode_inclined!(record_h5, record, "H", k, m, l, mode_index, mode, K)
                    get(mode, "NearStaticOverflow", false) && (near_static_zero_hor += 1)
                    get(mode, "FailedMode", false) && (failed_modes_hor += 1)
                    l_reached_hor = max(l_reached_hor, l)
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
                lmax_hit_hor = lmax_hit_hor || _tail_hit(l_tail, tol, total_horizon_energy + shell_k_horizon_energy + shell_m_horizon_energy)

                shell_k_horizon_energy += shell_m_horizon_energy
                shell_k_horizon_angular += shell_m_horizon_angular
                shell_k_horizon_carter += shell_m_horizon_carter
                _push_trail!(m_trail, shell_m_horizon_energy)

                # Horizon shells are negative in the superradiant regime, so the shell ratio compares magnitudes.
                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
                below_count_m = !_before_m_peak(m, k, x) && (abs(shell_m_horizon_energy) <= layer_threshold) && (abs(shell_m_horizon_energy) < peak_m_horizon) ? below_count_m + 1 : 0
                peak_m_horizon = max(peak_m_horizon, abs(shell_m_horizon_energy))
                if below_count_m >= minimum_consecutive
                    break
                end
            end
            below_count_m < minimum_consecutive && _tail_hit(_geometric_tail(m_trail, lmax), tol, total_horizon_energy + shell_k_horizon_energy) && (lmax_hit_hor = true)

            below_count_m = 0
            m_trail = Float64[]
            for m in -1:-1:-lmax
                lmin = max(2, abs(m))
                shell_m_horizon_energy = 0.0
                shell_m_horizon_angular = 0.0
                shell_m_horizon_carter = 0.0
                below_count_l = 0

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
                l_tail = _branch_l_walk(lmin, _branch_peak_l(lmin, m, k, x), lmax, layer_threshold, minimum_consecutive, true) do l
                    mode = let; with_submission_mode(scope; mode = (2, l, m, 0, k)) do mode_cache
                        ConvolutionIntegrals.convolution_integral_inclined_trapezoidal_isem(KG_sample, 2, l, m, k, K; tol = tol, sample_tol = sample_tol, max_flux = 2 * E_estimate_hor, mode_abs_floor = MODE_REFINE_FLOOR_FRACTION * branch_mode_floor_hor[], ci_kwargs..., cache = mode_cache)
                    end; end
                    _record_current_mode_inclined!(record_h5, record, "H", k, m, l, mode_index, mode, K)
                    get(mode, "NearStaticOverflow", false) && (near_static_zero_hor += 1)
                    get(mode, "FailedMode", false) && (failed_modes_hor += 1)
                    l_reached_hor = max(l_reached_hor, l)
                    shell_m_horizon_energy += 2 * mode["EnergyFlux"]
                    shell_m_horizon_angular += 2 * mode["AngularMomentumFlux"]
                    shell_m_horizon_carter += 2 * mode["CarterConstantFlux"]
                    total_modes += 2
                    mode["EnergyFlux"]
                end
                lmax_hit_hor = lmax_hit_hor || _tail_hit(l_tail, tol, total_horizon_energy + shell_k_horizon_energy + shell_m_horizon_energy)

                shell_k_horizon_energy += shell_m_horizon_energy
                shell_k_horizon_angular += shell_m_horizon_angular
                shell_k_horizon_carter += shell_m_horizon_carter
                _push_trail!(m_trail, shell_m_horizon_energy)

                layer_threshold = _mode_cutoff_threshold(tol * abs(shell_k_horizon_energy_last), branch_mode_floor_hor[])
                below_count_m = !_before_m_peak(m, k, x) && (abs(shell_m_horizon_energy) <= layer_threshold) && (abs(m) >= Max_m) ? below_count_m + 1 : 0
                if below_count_m >= minimum_consecutive
                    shell_k_horizon_energy_last = shell_k_horizon_energy
                    break
                end
            end
            below_count_m < minimum_consecutive && _tail_hit(_geometric_tail(m_trail, lmax), tol, total_horizon_energy + shell_k_horizon_energy) && (lmax_hit_hor = true)

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

        # A k sum that ran to kmax is truncated when the geometric tail of its last three k shells matters.
        kmax_hit_inf = k_reached_inf === nothing && _tail_hit(_geometric_tail(abs.(Energy_flux_inf[max(1, end - 2):end]), kmax), tol, total_infinity_energy)
        kmax_hit_hor = k_reached_hor === nothing && _tail_hit(_geometric_tail(abs.(Energy_flux_hor[max(1, end - 2):end]), kmax), tol, total_horizon_energy)
        if warn_truncation && (lmax_hit_inf || lmax_hit_hor)
            @warn "Inclined mode summation hit lmax = $(lmax) with modes still above the cutoff: the l sum is truncated. Raise lmax." a p x lmax l_reached_inf l_reached_hor lmax_hit_inf lmax_hit_hor
        end
        if warn_truncation && (kmax_hit_inf || kmax_hit_hor)
            @warn "Inclined mode summation ran to kmax = $(kmax) without $(minimum_consecutive + 1) consecutive negligible k shells: the k sum may be truncated. Raise kmax." a p x kmax kmax_hit_inf kmax_hit_hor
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
            record_h5["meta/l_reached_inf"] = l_reached_inf
            record_h5["meta/l_reached_hor"] = l_reached_hor
            record_h5["meta/lmax_hit_inf"] = lmax_hit_inf
            record_h5["meta/lmax_hit_hor"] = lmax_hit_hor
            record_h5["meta/near_static_zero_inf"] = near_static_zero_inf
            record_h5["meta/near_static_zero_hor"] = near_static_zero_hor
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
            l_reached_inf = l_reached_inf,
            l_reached_hor = l_reached_hor,
            lmax_hit_inf = lmax_hit_inf,
            lmax_hit_hor = lmax_hit_hor,
            kmax_hit_inf = kmax_hit_inf,
            kmax_hit_hor = kmax_hit_hor,
            near_static_zero_inf = near_static_zero_inf,
            near_static_zero_hor = near_static_zero_hor,
            failed_modes_inf = failed_modes_inf,
            failed_modes_hor = failed_modes_hor,
            infinity_energy_flux_list = Energy_flux_inf,
            horizon_energy_flux_list = Energy_flux_hor,
            record_file = record ? record_path : nothing,
        )
    finally
        record && _public_close_resource!(scope, record_h5)
    end
end

# Total inclined circular fluxes, summed over k shells as the eccentric summation sums n shells.
function inclined_mode_summation(a, p, x; submission = nothing, K = 16, K0 = K, Kmax = submission === nothing ? 2^12 : submission.Kmax, tol = 1e-8, lmax = :auto, kmax = :auto, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing, Nmax = submission === nothing ? 2^14 : submission.Nmax)
    caps = _auto_caps_start(a, p, 0.0, x, tol, lmax, kmax)
    scope, owned = ensure_public_submission(submission, a, p, 0.0, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        _with_auto_caps(caps) do L, Kc, Nc, warn
            with_submission_operation(scope) do
                _inclined_mode_summation_scoped(scope, a, p, x; K = K, K0 = K0, Kmax = Kmax, tol = tol, lmax = L, kmax = Kc, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, warn_truncation = warn)
            end
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
        return eccentric_mode_summation(a, p, e; x = x, N = N0, N0 = N0, Nmax = Nmax, tol = tol, lmax = lmax, nmax = nmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, tail_levin = tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_local_n = levin_local_n, levin_max_depth = levin_max_depth, submission = scope, lifecycle_observer = scope.lifecycle_observer, geometry_observer = scope.geometry_observer)
    end
    if e == 0.0
        return inclined_mode_summation(a, p, x; K = K0, K0 = K0, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, submission = isequal(e, 0.0) ? scope : nothing, stop_requested = isequal(e, 0.0) ? nothing : (() -> _mode_flux_parent_stop_requested(scope)), lifecycle_observer = scope.lifecycle_observer, geometry_observer = scope.geometry_observer)
    end

    record_path = record_path === nothing ? "generic_mode_data_a_$(_record_tag(a))_p_$(_record_tag(p))_e_$(_record_tag(e))_x_$(_record_tag(x)).h5" : record_path

    KG = submission_orbit!(scope)
    radial_over_azimuthal = KG["Frequencies"]["ϒr"] / KG["Frequencies"]["ϒϕ"]
    F_orbit = KG["Frequencies"]
    polar_over_azimuthal = KG["Frequencies"]["ϒθ"] / KG["Frequencies"]["ϒϕ"]

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
        # n-shell trigger of the generic tail (forced Levin for the shell and the m-walk tail shortcut). It is decided
        # per shell from the radial trapezoid grid the shell's (m = 2, k = 0) mode needs, and is not latched: a shell
        # whose phase is slow again goes back to the trapezoid. start_n records the first shell that triggered.
        use_generic_tail_levin = (s, n) -> begin
            branch_policy = s == -2 ? tail_policy_inf : tail_policy_hor
            _tail_enabled(branch_policy) && fast || return false
            active = s == -2 ? tail_levin_infinity_active : tail_levin_horizon_active
            start_n = s == -2 ? tail_levin_infinity_start_n : tail_levin_horizon_start_n
            triggered = _use_levin_tail(n; nmin = levin_nmin) &&
                (branch_policy === true || _gen_levin_auto_mode(
                    _mode_radial_trapezoid_n(KG_master, 2, n, 0, N0, Nmax), n, K0, Kmax,
                    levin_local_n, levin_max_depth, levin_nmin))
            if triggered
                active[] = true
                start_n[] === nothing && (start_n[] = n)
            end
            return triggered
        end
        # one tracker per boundary: the two branches may run on different tasks
        cap_trackers = (CapTracker(lmax, kmax), CapTracker(lmax, kmax))
        failed_modes = [0, 0]
        near_static_zero = [0, 0]
        eval_mode = (s, l, m, n, k, max_flux, use_levin=false) -> begin
            auto_levin = use_generic_tail_levin(s, n)
            use_levin = use_levin || auto_levin ||
                (fast && (s == -2 ? tail_policy_inf : tail_policy_hor) === :auto &&
                    _gen_levin_auto_mode(_mode_radial_trapezoid_n(KG_master, m, n, k, N0, Nmax), n, K0, Kmax,
                        levin_local_n, levin_max_depth, levin_nmin))
            floor = MODE_REFINE_FLOOR_FRACTION * branch_mode_floor(s)
            ω_mode = (m * F_orbit["ϒϕ"] + n * F_orbit["ϒr"] + k * F_orbit["ϒθ"]) / F_orbit["ϒt"]
            mode = if _negligible_low_frequency_mode(s, l, ω_mode, F_orbit["ϒϕ"] / F_orbit["ϒt"], a, p, e)
                res = ConvolutionIntegrals._zero_radiative_mode(ω_mode, nothing; reason = "negligible_low_frequency")
                res["NearStaticOverflow"] = true
                res
            else
            let; ConvolutionIntegrals._with_coverage_retry(s, l, m, n, k; omega = ω_mode, Ωφ = F_orbit["ϒϕ"] / F_orbit["ϒt"], rmax = p / (1 - e), negligible_if_failed = _negligible_low_frequency_mode(s, l, ω_mode, F_orbit["ϒϕ"] / F_orbit["ϒt"], a, p, e; bound = ConvolutionIntegrals.NEGLIGIBLE_FAILED_MODE_BOUND)) do; if fast
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
            end; end; end
            end
            get(mode, "FailedMode", false) && (failed_modes[s == -2 ? 1 : 2] += 1)
            get(mode, "NearStaticOverflow", false) && (near_static_zero[s == -2 ? 1 : 2] += 1)
            ConvolutionIntegrals._check_finite_mode_flux(mode, s, l, m, n, k)
            _note_cap_mode!(cap_trackers[s == -2 ? 1 : 2], l, m, n, k, mode["EnergyFlux"])
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
            last_shell = (0.0, 0.0, 0.0)
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
            peak_k_shell = 0.0
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

                # a k shell counts towards the stop only past the largest k shell (shells can rise before they fall)
                below_count_k = abs(shell_k) <= initial_k_shell_threshold && abs(shell_k) < peak_k_shell ? below_count_k + 1 : 0
                peak_k_shell = max(peak_k_shell, abs(shell_k))
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
                walk_max_l = 0.0
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
                    if below_count_l >= minimum_consecutive || _deep_stop(abs(mode["EnergyFlux"]), deep, threshold, walk_max_l)
                        break
                    end
                    walk_max_l = max(walk_max_l, abs(mode["EnergyFlux"]))
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
                    walk_max_l = 0.0
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
                        if below_count_l >= minimum_consecutive || _deep_stop(abs(mode["EnergyFlux"]), deep, threshold, walk_max_l)
                            break
                        end
                        walk_max_l = max(walk_max_l, abs(mode["EnergyFlux"]))
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
                    walk_max_l = 0.0
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
                        if below_count_l >= minimum_consecutive || _deep_stop(abs(mode["EnergyFlux"]), deep, threshold, walk_max_l)
                            break
                        end
                        walk_max_l = max(walk_max_l, abs(mode["EnergyFlux"]))
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
                        below_count_m = !growing_m && !_before_m_peak(m, k, x) && !_before_omega_zero(m, -(n * radial_over_azimuthal + k * polar_over_azimuthal)) && _positive_m_tail_below(shell_m, shell_k, threshold, tol, m, Max_m, mode_abs_floor, tail_cutoff) ? below_count_m + 1 : 0
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
                        below_count_m = !growing_m && !_before_m_peak(m, k, x) && !_before_omega_zero(m, -(n * radial_over_azimuthal + k * polar_over_azimuthal)) && _negative_m_tail_below(shell_m, threshold, neg_branch_scale, m, Max_m, tail_cutoff) ? below_count_m + 1 : 0
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
                        below_count_m = !growing_m && !_before_m_peak(m, k, x) && !_before_omega_zero(m, -(n * radial_over_azimuthal + k * polar_over_azimuthal)) && _positive_m_tail_below(shell_m, shell_k, threshold, tol, m, Max_m, mode_abs_floor, tail_cutoff) ? below_count_m + 1 : 0
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
                        below_count_m = !growing_m && !_before_m_peak(m, k, x) && !_before_omega_zero(m, -(n * radial_over_azimuthal + k * polar_over_azimuthal)) && _negative_m_tail_below(shell_m, threshold, neg_branch_scale, m, Max_m, tail_cutoff) ? below_count_m + 1 : 0
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
                last_shell = (shell_n, shell_n_angular, shell_n_carter)
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

            # Geometric n tail at high eccentricity, as for eccentric equatorial orbits.
            tail_factor, tail_ratio = _geometric_n_tail_factor(shell_list, e, n_reached, minimum_consecutive + 1)
            n_tail = tail_factor .* last_shell
            total_energy += n_tail[1]
            total_angular += n_tail[2]
            total_carter += n_tail[3]

            return (
                energy_flux = total_energy,
                angular_momentum_flux = total_angular,
                carter_constant_flux = total_carter,
                total_modes = total_modes,
                n_reached = n_reached,
                shell_list = shell_list,
                truncation_floor = branch_mode_floor(s),
                n_tail = (energy = n_tail[1], angular_momentum = n_tail[2], carter_constant = n_tail[3], ratio = tail_ratio),
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

        cap_hits_inf = _cap_hits(cap_trackers[1], tol, inf.energy_flux)
        cap_hits_hor = _cap_hits(cap_trackers[2], tol, hor.energy_flux)
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
            n_tail_infinity = inf.n_tail,
            n_tail_horizon = hor.n_tail,
            lmax_hit_inf = cap_hits_inf[1],
            lmax_hit_hor = cap_hits_hor[1],
            kmax_hit_inf = cap_hits_inf[2],
            kmax_hit_hor = cap_hits_hor[2],
            nmax_hit_inf = _nmax_hit(inf.n_reached, inf.shell_list, nmax, tol, inf.energy_flux),
            nmax_hit_hor = _nmax_hit(hor.n_reached, hor.shell_list, nmax, tol, hor.energy_flux),
            failed_modes_inf = failed_modes[1],
            failed_modes_hor = failed_modes[2],
            near_static_zero_inf = near_static_zero[1],
            near_static_zero_hor = near_static_zero[2],
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
function generic_mode_summation(a, p, e, x; submission = nothing, N0 = 64, K0 = 16, Nmax = submission === nothing ? 2^14 : submission.Nmax, Kmax = submission === nothing ? 2^12 : submission.Kmax, tol = 1e-8, lmax = :auto, kmax = :auto, nmax = :auto, minimum_consecutive = 2, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, neg_branch_scale = 0.1, beta = DEFAULT_GENERIC_BETA, deep = DEFAULT_GENERIC_DEEP, tail_levin = :auto, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_local_n::Int = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, levin_max_depth::Int = 8, progress_interval::Int = 0, progress_path = nothing, stop_requested = nothing, lifecycle_observer = nothing, geometry_observer = nothing)
    caps = _auto_caps_start(a, p, e, x, tol, lmax, kmax, nmax)
    scope, owned = ensure_public_submission(submission, a, p, e, x; Nmax, Kmax, stop_requested, lifecycle_observer, geometry_observer)
    return _with_public_scope(scope, owned) do
        _with_auto_caps(caps) do L, Kc, Nc, warn
            with_submission_operation(scope) do
                _generic_mode_summation_scoped(scope, a, p, e, x; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, tol = tol, lmax = L, kmax = Kc, nmax = Nc, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, neg_branch_scale = neg_branch_scale, beta=beta, deep=deep, tail_levin = tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_local_n = levin_local_n, levin_max_depth = levin_max_depth, progress_interval = progress_interval, progress_path = progress_path)
            end
        end
    end
end

end
