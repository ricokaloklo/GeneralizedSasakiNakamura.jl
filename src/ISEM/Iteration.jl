module Iteration

using ..CoefficientTables: CoefficientSet
using ..Parameters: horizon_tail
using ..LargeFrequencyExpansion: LargeFrequencyExpansionPlan, large_frequency_expansion_plan, large_frequency_expansion_ab!
using ..NearExtreme: near_extreme_selected
using ..Eikonal:
    EikonalPatch,
    EIKONAL_TERMS,
    EIKONAL_ORDER,
    LARGE_FREQUENCY_EXPANSION_EIKONAL_ORDER,
    EikonalScratch,
    eikonal_candidate,
    eikonal_pair,
    eikonal_path,
    eikonal_split,
    eikonal_split_pair,
    eikonal_state,
    eikonal_value
using ..HorizonLargeFrequencyExpansion:
    HorizonLargeFrequencyExpansionPatch,
    large_frequency_expansion_horizon_pair,
    large_frequency_expansion_horizon_single,
    large_frequency_expansion_horizon_state,
    large_frequency_expansion_horizon_value
using ..InfinityLargeFrequencyExpansion:
    InfinityLargeFrequencyExpansionPatch,
    large_frequency_expansion_infinity_anchor,
    large_frequency_expansion_infinity_anchor_pair,
    large_frequency_expansion_infinity_pair,
    large_frequency_expansion_infinity_reach,
    large_frequency_expansion_infinity_reach_pair,
    large_frequency_expansion_infinity_single,
    large_frequency_expansion_infinity_state,
    large_frequency_expansion_infinity_value
using ..LocalSolutionAtInfinity:
    infinity_guarded,
    infinity_exponents,
    infinity_series_triple
using ..LocalSolutionAtZero:
    is_horizon_threshold,
    horizon_series_pair,
    horizon_series_triple,
    horizon_series_value,
    zero_log_solution,
    zero_local_solution
using ..OrdinaryPointExpansion:
    endpoint_ab_series,
    ordinary_ab_series!,
    ordinary_scaled_y_ab_series!,
    poly_pair,
    poly_triple,
    poly_value,
    series_issue

export Patch, ScaledPatch, LogScaledPatch
export LogScaledState, LogScaledBasis
export Basis, InfinitySelection, IterationScratch, Truncation
export EndpointCertificateError
export endpoint_seed_scale, select_horizon_seed, select_infinity_endpoint
export infinity_seed_wronskian
export iterate_from_zero, iterate_pair_from_zero
export iterate_from_infinity, iterate_pair_from_infinity
export iterate_from_state, iterate_pair_from_state
export iterate_logscaled_from_state
export iterate_scaled_y
export basis_value, basis_state, basis_patch_count
export logscaled_state, logscaled_basis_state
export materialize_logscaled_state

const _RADIUS_SAFETY = 0.8
const _ORDINARY_STEP_SAFETY = 0.9
const _MIN_STEP = 1e-5
const _MIN_INFINITY_ENDPOINT = 1024eps(Float64)
const _MAX_REJECTIONS = 160
const _ENDPOINT_SEED_TARGET = 1e-13
const _HORIZON_REACH_ORDER = 64
const _HORIZON_REACH_COST_RATIO = 0.75
const _HORIZON_REACH_MIN_SAVING = 2.0
# Load-time environment switches (read once when the module loads):
# ISEM_HORIZON_REACH_AUTO=1 (default) lets the horizon seed be re-selected with an
# order-_HORIZON_REACH_ORDER series when that shortens the march (see
# `_horizon_seed`); ISEM_ORDINARY_EIKONAL=1 (default) allows eikonal patches in
# the ordinary regime; ISEM_HORIZON_LARGE_FREQUENCY_EXPANSION=1 (default 0)
# compiles in the horizon large-frequency carrier route (`_large_frequency_expansion_horizon_basis`);
# ISEM_INFINITY_LARGE_FREQUENCY_EXPANSION=1 (default) compiles in the infinity
# large-frequency carrier of the hybrid marches and `iterate_from_infinity`.
const _HORIZON_REACH_AUTO = get(
    ENV, "ISEM_HORIZON_REACH_AUTO", "1") == "1"
const _INFINITY_ENDPOINT_TARGET = 1e-10
const _INFINITY_EIKONAL_LIMIT = 1e-12 + 256eps(Float64)
const _LARGE_FREQUENCY_EXPANSION_EIKONAL_TERMS = 3
const _EIKONAL_SLOW_VARIATION = 0.015
const _EIKONAL_CERTIFIED_VARIATION = 0.02
const _ORDINARY_EIKONAL = get(
    ENV, "ISEM_ORDINARY_EIKONAL", "1") == "1"
const _HORIZON_LARGE_FREQUENCY_EXPANSION = get(ENV, "ISEM_HORIZON_LARGE_FREQUENCY_EXPANSION", "0") == "1"
const _INFINITY_LARGE_FREQUENCY_EXPANSION = get(ENV, "ISEM_INFINITY_LARGE_FREQUENCY_EXPANSION", "1") == "1"
const _SCALED_Y_CANDIDATES = (0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.48, 0.64, 0.8, 1.0)
const _INFINITY_ENDPOINT_PHASES = (
    1.0,
    0.8408964152537145,
    0.7071067811865476,
    0.5946035575013605,
)
const _SMALL_FREQUENCY_EXPANSION_BRIDGE_ORDER = 32
const _SMALL_FREQUENCY_EXPANSION_BRIDGE_TARGET_Y = 1.0e-3
const _SMALL_FREQUENCY_EXPANSION_BRIDGE_STEP_SAFETY = 0.8
const _SMALL_FREQUENCY_EXPANSION_BRIDGE_TOLERANCE = 1.0e-13

@inline _eikonal_allowed(regime::Symbol) =
    regime != :ordinary || _ORDINARY_EIKONAL

# Taylor patch about center_x, valid through next_x: the solution is scale times
# sum_k coeffs[k+1] (x - center_x)^k.
struct Patch
    center_x::Float64
    next_x::Float64
    coeffs::Vector{ComplexF64}
    scale::Float64
end

Patch(center_x, next_x, coeffs) =
    Patch(center_x, next_x, coeffs, 1.0)

# Taylor patch in Y = (1 - x) / omega_scale about center_y, valid through next_x;
# the solution is scale times the series in Y - center_y.
struct ScaledPatch
    center_x::Float64
    next_x::Float64
    center_y::Float64
    omega_scale::Float64
    coeffs::Vector{ComplexF64}
    scale::Float64
end

# State (X, dXdx) normalized so that max(|X|, |dXdx|) = 1; the true state is
# exp(log_scale) times the stored one.
struct LogScaledState
    X::ComplexF64
    dXdx::ComplexF64
    log_scale::Float64
end

# Taylor patch like `Patch` whose amplitude is stored as its natural log, for
# propagation through regions of exponential growth or decay.
struct LogScaledPatch
    center_x::Float64
    next_x::Float64
    coeffs::Vector{ComplexF64}
    log_scale::Float64
end

"""
    Truncation

Record of a series used below its requested order: where (`location`), for
which branch, the requested and effective orders, and the reason.
"""
struct Truncation
    location::Symbol
    branch::Symbol
    requested_order::Int
    effective_order::Int
    reason::Symbol
end

# Result of `iterate_logscaled_from_state`: the seed state, the log-scaled
# patches to the match point, the state there, step statistics and truncations.
struct LogScaledBasis
    kind::Symbol
    seed_x::Float64
    seed_state::LogScaledState
    patches::Vector{LogScaledPatch}
    match_state::LogScaledState
    step_count::Int
    max_step::Float64
    min_step::Float64
    truncations::Vector{Truncation}
end

"""
    Basis

One propagated solution: its kind (`:horizon_in`, `:horizon_out`,
`:infinity_out` or `:infinity_in`), the endpoint series and seed, the patches
from the seed to the match point, and the state there (`match_X`,
`match_dXdx`). `endpoint_valid` is false when the basis was started from a
state rather than an endpoint series. At the superradiant threshold the outgoing
horizon solution carries a logarithm, X = log(x) X_log + x^mu sum c_n x^n; the
endpoint series of X_log is then stored in `endpoint_log_coeffs` (empty otherwise)
with exponent `endpoint_log_exponent`.
"""
struct Basis{P}
    kind::Symbol
    scale::ComplexF64
    exponent1::ComplexF64
    exponent2::ComplexF64
    seed_x::Float64
    seed_coord::Float64
    seed_X::ComplexF64
    seed_dXdx::ComplexF64
    endpoint_coeffs::Vector{ComplexF64}
    patches::Vector{P}
    match_X::ComplexF64
    match_dXdx::ComplexF64
    step_count::Int
    max_step::Float64
    min_step::Float64
    endpoint_variable_scale::ComplexF64
    truncations::Vector{Truncation}
    endpoint_valid::Bool
    endpoint_log_coeffs::Vector{ComplexF64}
    endpoint_log_exponent::ComplexF64
end

Basis(
    kind,
    scale,
    exponent1,
    exponent2,
    seed_x,
    seed_coord,
    seed_X,
    seed_dXdx,
    endpoint_coeffs,
    patches,
    match_X,
    match_dXdx,
    step_count,
    max_step,
    min_step,
    endpoint_variable_scale,
    truncations,
    endpoint_valid,
) = Basis(
    kind,
    scale,
    exponent1,
    exponent2,
    seed_x,
    seed_coord,
    seed_X,
    seed_dXdx,
    endpoint_coeffs,
    patches,
    match_X,
    match_dXdx,
    step_count,
    max_step,
    min_step,
    endpoint_variable_scale,
    truncations,
    endpoint_valid,
    ComplexF64[],
    0.0 + 0.0im,
)

Basis(
    kind,
    scale,
    exponent1,
    exponent2,
    seed_x,
    seed_coord,
    seed_X,
    seed_dXdx,
    endpoint_coeffs,
    patches,
    match_X,
    match_dXdx,
    step_count,
    max_step,
    min_step,
    endpoint_variable_scale,
) = Basis(
    kind,
    scale,
    exponent1,
    exponent2,
    seed_x,
    seed_coord,
    seed_X,
    seed_dXdx,
    endpoint_coeffs,
    patches,
    match_X,
    match_dXdx,
    step_count,
    max_step,
    min_step,
    endpoint_variable_scale,
    Truncation[],
    true,
)

Basis(
    kind,
    scale,
    exponent1,
    exponent2,
    seed_x,
    seed_coord,
    seed_X,
    seed_dXdx,
    endpoint_coeffs,
    patches,
    match_X,
    match_dXdx,
    step_count,
    max_step,
    min_step,
    endpoint_variable_scale,
    truncations,
) = Basis(
    kind,
    scale,
    exponent1,
    exponent2,
    seed_x,
    seed_coord,
    seed_X,
    seed_dXdx,
    endpoint_coeffs,
    patches,
    match_X,
    match_dXdx,
    step_count,
    max_step,
    min_step,
    endpoint_variable_scale,
    truncations,
    true,
)

# Chosen infinity endpoint: endpoint_y (y = 1 - x), the truncated :out and :in
# series with their orders, the selection score, the variable scale of the
# series variable y / variable_scale, the requested order and the truncation
# reasons of the two series.
struct InfinitySelection
    endpoint_y::Float64
    out_coeffs::Vector{ComplexF64}
    in_coeffs::Vector{ComplexF64}
    out_order::Int
    in_order::Int
    score::Float64
    variable_scale::ComplexF64
    requested_order::Int
    out_reason::Symbol
    in_reason::Symbol
end

# Thrown when no endpoint coordinate meets the seed target: `endpoint` is
# :horizon or :infinity, with the target, the best score and coordinate tried,
# and the series order in use.
struct EndpointCertificateError <: Exception
    endpoint::Symbol
    target::Float64
    best_score::Float64
    best_coord::Float64
    requested_order::Int
end

function Base.showerror(io::IO, err::EndpointCertificateError)
    print(io,
        "ISEM GSN ", err.endpoint,
        " endpoint has no certified seed: best score=", err.best_score,
        ", target=", err.target,
        ", coordinate=", err.best_coord,
        ", requested_order=", err.requested_order)
end


InfinitySelection(
    endpoint_y,
    out_coeffs,
    in_coeffs,
    out_order,
    in_order,
    score,
    variable_scale,
) = InfinitySelection(
    endpoint_y,
    out_coeffs,
    in_coeffs,
    out_order,
    in_order,
    score,
    variable_scale,
    max(out_order, in_order),
    :none,
    :none,
)

# Work buffers for the ordinary-point marches: `pq` for the coefficient-table
# values, `avec`/`bvec` for the local A/B series, `local_coeffs`/`local_coeffs2`
# for the Taylor coefficients of up to two solutions, `radius_buffer`/
# `radius_buffer2` for the radius estimate, `denom_inv[n+1] = 1/((n+2)(n+1))`,
# and the optional eikonal scratch, large-frequency plan and near-extreme
# double-double plan (`ordinary_dd`).
struct IterationScratch
    pq::Vector{ComplexF64}
    avec::Vector{ComplexF64}
    bvec::Vector{ComplexF64}
    local_coeffs::Vector{ComplexF64}
    local_coeffs2::Vector{ComplexF64}
    radius_buffer::Vector{Float64}
    radius_buffer2::Vector{Float64}
    denom_inv::Vector{Float64}
    eikonal::Union{Nothing,EikonalScratch}
    large_frequency_expansion::Union{Nothing,LargeFrequencyExpansionPlan}
    ordinary_dd::Union{Nothing,LargeFrequencyExpansionPlan}
end

function IterationScratch(
    coefficients::CoefficientSet,
    order::Integer;
    force_eikonal::Bool=false,
    large_frequency_expansion::Bool=false,
    small_frequency_expansion::Bool=false,
    route_branch::Union{Nothing,Symbol}=nothing,
    ordinary_dd::Bool=near_extreme_selected(coefficients),
)
    order_int = Int(order)
    order_int >= 0 || throw(ArgumentError("ordinary order must be non-negative."))
    regime = large_frequency_expansion ? :large_frequency_expansion : small_frequency_expansion ? :small_frequency_expansion : :ordinary
    use_eikonal = force_eikonal ||
        (_eikonal_allowed(regime) &&
         eikonal_candidate(coefficients, regime, route_branch))
    eikonal_order = large_frequency_expansion ? LARGE_FREQUENCY_EXPANSION_EIKONAL_ORDER : EIKONAL_ORDER
    buffer_order = use_eikonal ? max(order_int, eikonal_order) : order_int
    return IterationScratch(
        Vector{ComplexF64}(undef, coefficients.ordinary.value_count),
        Vector{ComplexF64}(undef, buffer_order + 1),
        Vector{ComplexF64}(undef, buffer_order + 1),
        Vector{ComplexF64}(undef, buffer_order + 1),
        Vector{ComplexF64}(undef, buffer_order + 1),
        Vector{Float64}(undef, buffer_order + 1),
        Vector{Float64}(undef, buffer_order + 1),
        buffer_order >= 2 ? [1.0 / ((n + 2) * (n + 1)) for n in 0:(buffer_order - 2)] : Float64[],
        use_eikonal ? EikonalScratch(eikonal_order) : nothing,
        large_frequency_expansion ? large_frequency_expansion_plan(coefficients.params,
            coefficients.ordinary.basis, buffer_order) : nothing,
        ordinary_dd ? large_frequency_expansion_plan(coefficients.params,
            coefficients.ordinary.basis, buffer_order) : nothing,
    )
end

function _eikonal_enabled(controls, coefficients, scratch=nothing)
    regime = getproperty(controls, :large_frequency_expansion) ? :large_frequency_expansion :
        getproperty(controls, :small_frequency_expansion) ? :small_frequency_expansion : :ordinary
    _eikonal_allowed(regime) || return false
    scratch === nothing || return scratch.eikonal !== nothing
    return eikonal_candidate(coefficients, regime)
end

_scratch(coefficients, order, controls) = IterationScratch(
    coefficients, order;
    force_eikonal=getproperty(controls, :large_frequency_expansion),
    large_frequency_expansion=getproperty(controls, :large_frequency_expansion),
    small_frequency_expansion=getproperty(controls, :small_frequency_expansion))

function _local_ab!(work, coefficients, x, order)
    work.large_frequency_expansion === nothing || return large_frequency_expansion_ab!(
        work.avec, work.bvec, work.large_frequency_expansion, x, order)
    work.ordinary_dd === nothing || return large_frequency_expansion_ab!(
        work.avec, work.bvec, work.ordinary_dd, x, order)
    return ordinary_ab_series!(
        work.avec, work.bvec, work.pq, coefficients, x, order)
end

_param_a(params) = hasproperty(params, :a) ? getproperty(params, :a) : sin(getproperty(params, :nu))
_param_kappa(params) = hasproperty(params, :kappa) ? getproperty(params, :kappa) : cos(getproperty(params, :nu))
_param_omega(params) = getproperty(params, :omega)
_param_m(params) = getproperty(params, :m)

# Slow-variation indicator at a point from f = a'/2 + a^2/4 - b (the invariant
# of the normal form u'' = f u): max(|f'|/|f|^(3/2), sqrt(|f''/2|)/|f|).
# Small values mean the eikonal form is accurate.
function _eikonal_indicator(a_coeffs, b_coeffs)
    length(a_coeffs) >= 4 && length(b_coeffs) >= 3 || return Inf
    @inbounds begin
        a0, a1, a2, a3 = a_coeffs[1], a_coeffs[2], a_coeffs[3], a_coeffs[4]
        b0, b1, b2 = b_coeffs[1], b_coeffs[2], b_coeffs[3]
    end
    f0 = a1 / 2 + a0^2 / 4 - b0
    isfinite(real(f0)) && isfinite(imag(f0)) && !iszero(f0) || return Inf
    f1 = a2 + a0 * a1 / 2 - b1
    f2 = 3a3 / 2 + a0 * a2 / 2 + a1^2 / 4 - b2
    e1 = abs(f1) / max(abs(f0)^(3 / 2), eps(Float64))
    e2 = sqrt(abs(f2) / max(abs(f0)^2, eps(Float64)))
    return max(e1, e2)
end

function _eikonal_score(coefficients, x, work)
    a_terms, b_terms = _local_ab!(work, coefficients, x, 4)
    return _eikonal_indicator(
        @view(work.avec[1:a_terms]), @view(work.bvec[1:b_terms]))
end

function _eikonal_ready(coefficients, x, work)
    work.eikonal === nothing && return false
    return _eikonal_score(coefficients, x, work) <=
        _EIKONAL_SLOW_VARIATION
end

function _normalized_horizon_branch(branch::Symbol)
    branch in (:in, :ingoing, :IN) && return :in
    branch in (:out, :outgoing, :OUT) && return :out
    throw(ArgumentError("horizon branch must be :in or :out."))
end

function _normalized_infinity_branch(branch::Symbol)
    branch in (:out, :outgoing, :up, :UP) && return :out
    branch in (:in, :ingoing, :down, :DOWN) && return :in
    throw(ArgumentError("infinity branch must be :out/:up or :in/:down."))
end

function _kind(endpoint::Symbol, branch::Symbol)
    endpoint == :horizon && return branch == :in ? :horizon_in : :horizon_out
    endpoint == :infinity && return branch == :out ? :infinity_out : :infinity_in
    throw(ArgumentError("endpoint must be :horizon or :infinity."))
end

_is_horizon(kind::Symbol) = kind in (:horizon_in, :horizon_out)

# Constant phase that normalizes each endpoint branch: exp(-/+ i p (r_+ + 2 log kappa))
# for horizon in/out with p = omega - m a / (2 r_+), and
# exp(+/- i omega (r_- + 2 log kappa)) for infinity out/in, with r_+/- = 1 +/- kappa.
function endpoint_seed_scale(kind::Symbol, params)
    a = _param_a(params)
    kappa = _param_kappa(params)
    omega = _param_omega(params)
    m = _param_m(params)
    r_plus = 1 + kappa
    r_minus = 1 - kappa
    p = omega - m * a / (2 * r_plus)
    horizon_constant = r_plus + 2 * log(kappa)
    infinity_constant = r_minus + 2 * log(kappa)
    kind == :horizon_in && return ComplexF64(exp(-1im * p * horizon_constant))
    kind == :horizon_out && return ComplexF64(exp(1im * p * horizon_constant))
    kind == :infinity_out && return ComplexF64(exp(1im * omega * infinity_constant))
    kind == :infinity_in && return ComplexF64(exp(-1im * omega * infinity_constant))
    throw(ArgumentError("unknown ISEM GSN endpoint kind $kind."))
end

# Infinity series exp(Sigma/Y) Y^rho sum_n c_n Y^n with Y = y / variable_scale,
# times variable_scale^rho; the pair and triple forms add x-derivatives
# (dx = -dy).
function _scaled_infinity_series_value(coeffs, rho, Sigma, y, variable_scale)
    Y = y / variable_scale
    scale = variable_scale^rho
    return scale * exp(Sigma / Y) * Y^rho * poly_value(coeffs, Y)
end

function _scaled_infinity_series_pair(coeffs, rho, Sigma, y, variable_scale)
    Y = y / variable_scale
    series, dseries = poly_pair(coeffs, Y)
    scale = variable_scale^rho
    prefactor = exp(Sigma / Y) * Y^rho
    value = scale * prefactor * series
    derivative_Y = prefactor * (dseries + (-Sigma / Y^2 + rho / Y) * series)
    return value, -scale * derivative_Y / variable_scale
end

function _scaled_infinity_series_triple(coeffs, rho, Sigma, y, variable_scale)
    Y = y / variable_scale
    series, dseries, ddseries = poly_triple(coeffs, Y)
    scale = variable_scale^rho
    prefactor = exp(Sigma / Y) * Y^rho
    log_derivative = -Sigma / Y^2 + rho / Y
    log_second = 2 * Sigma / Y^3 - rho / Y^2
    value = scale * prefactor * series
    derivative_Y = prefactor * (dseries + log_derivative * series)
    second_Y = prefactor * (
        ddseries + 2 * log_derivative * dseries +
        (log_derivative^2 + log_second) * series
    )
    inv_scale = inv(variable_scale)
    return value, -scale * derivative_Y * inv_scale, scale * second_Y * inv_scale^2
end

function _endpoint_value(kind::Symbol, coeffs, exponent1, exponent2, x, variable_scale=1.0 + 0.0im)
    if _is_horizon(kind)
        return horizon_series_value(coeffs, exponent1, x)
    end
    y = 1 - x
    return _scaled_infinity_series_value(coeffs, exponent1, exponent2, y, variable_scale)
end

function _endpoint_pair(kind::Symbol, coeffs, exponent1, exponent2, coord, variable_scale=1.0 + 0.0im)
    if _is_horizon(kind)
        return horizon_series_pair(coeffs, exponent1, coord)
    end
    return _scaled_infinity_series_pair(coeffs, exponent1, exponent2, coord, variable_scale)
end

function _endpoint_state(kind::Symbol, coeffs, exponent1, exponent2, x, variable_scale=1.0 + 0.0im)
    if _is_horizon(kind)
        value, derivative, _ = horizon_series_triple(coeffs, exponent1, x)
        return value, derivative
    end
    value, derivative, _ = _scaled_infinity_series_triple(
        coeffs,
        exponent1,
        exponent2,
        1 - x,
        variable_scale,
    )
    return value, derivative
end

function _endpoint_tail_radius(coeffs::Vector{ComplexF64}, target::Float64,
        cap::Float64)
    order = length(coeffs) - 1
    order <= 0 && return cap
    radius = cap
    found = false
    first_degree = max(1, div(order, 2))
    @inbounds for degree in first_degree:order
        magnitude = abs(coeffs[degree + 1])
        iszero(magnitude) && continue
        candidate = (target / magnitude)^(1 / degree)
        if isfinite(candidate) && candidate > 0
            radius = min(radius, candidate)
            found = true
        end
    end
    return found ? radius : cap
end

function _endpoint_last_term_radius(coeffs::Vector{ComplexF64},
        target::Float64, cap::Float64)
    order = length(coeffs) - 1
    order <= 0 && return cap
    magnitude = abs(coeffs[end])
    iszero(magnitude) && return cap
    radius = (target / magnitude)^(1 / order)
    return isfinite(radius) && radius > 0 ? min(cap, radius) : cap
end

# Median of the ratios |c_(n-1)/c_n| over the upper half of the series, an
# estimate of its radius of convergence, capped at `cap`.
function _endpoint_effective_radius!(buffer::Vector{Float64}, coeffs, cap::Float64)
    count = 0
    first_idx = max(2, div(length(coeffs), 2))
    @inbounds for idx in first_idx:(length(coeffs) - 1)
        denom = abs(coeffs[idx + 1])
        numer = abs(coeffs[idx])
        if denom > eps(Float64) && numer > 0
            r = numer / denom
            if isfinite(r) && r > 0
                count += 1
                buffer[count] = r
            end
        end
    end
    count == 0 && return cap
    @inbounds for i in 2:count
        value = buffer[i]
        j = i - 1
        while j >= 1 && buffer[j] > value
            buffer[j + 1] = buffer[j]
            j -= 1
        end
        buffer[j + 1] = value
    end
    radius = isodd(count) ? buffer[cld(count, 2)] :
        0.5 * (buffer[div(count, 2)] + buffer[div(count, 2) + 1])
    return min(cap, radius)
end

# Taylor coefficients about an ordinary point of X'' + a X' + b X = 0 with
# X(0) = c0, X'(0) = c1, given the series of a and b.
function _ordinary_point_coeffs!(
    coeffs::Vector{ComplexF64},
    a_coeffs,
    b_coeffs,
    c0,
    c1,
    order::Integer,
    denom_inv::Vector{Float64},
)
    coeffs[1] = ComplexF64(c0)
    order >= 1 && (coeffs[2] = ComplexF64(c1))
    order < 2 && return coeffs
    @inbounds for n in 0:(order - 2)
        numerator = 0.0 + 0.0im
        for j in 0:n
            numerator += a_coeffs[j + 1] * (n - j + 1) * coeffs[n - j + 2]
            numerator += b_coeffs[j + 1] * coeffs[n - j + 1]
        end
        coeffs[n + 3] = -numerator * denom_inv[n + 1]
    end
    return coeffs
end

function _ordinary_point_pair!(
    coeffs1::Vector{ComplexF64},
    coeffs2::Vector{ComplexF64},
    a_coeffs,
    b_coeffs,
    c01,
    c11,
    c02,
    c12,
    order::Integer,
    denom_inv::Vector{Float64},
)
    coeffs1[1] = ComplexF64(c01)
    coeffs2[1] = ComplexF64(c02)
    if order >= 1
        coeffs1[2] = ComplexF64(c11)
        coeffs2[2] = ComplexF64(c12)
    end
    order < 2 && return coeffs1, coeffs2
    @inbounds for n in 0:(order - 2)
        numerator1 = 0.0 + 0.0im
        numerator2 = 0.0 + 0.0im
        for j in 0:n
            aw = a_coeffs[j + 1] * (n - j + 1)
            b = b_coeffs[j + 1]
            numerator1 += aw * coeffs1[n - j + 2]
            numerator1 += b * coeffs1[n - j + 1]
            numerator2 += aw * coeffs2[n - j + 2]
            numerator2 += b * coeffs2[n - j + 1]
        end
        factor = -denom_inv[n + 1]
        coeffs1[n + 3] = numerator1 * factor
        coeffs2[n + 3] = numerator2 * factor
    end
    return coeffs1, coeffs2
end

@inline _state_scale(value, derivative) =
    max(abs(value), abs(derivative), floatmin(Float64))

# Checks of a horizon seed x: relative ODE residual of x^mu sum c_n x^n with
# A = a(x)/x, B = b(x)/x^2 (`_horizon_residual_metric`); round-off condition
# number of the series value and derivative (`_horizon_eval_condition`); and the
# relative size of the upper half of the series (`_horizon_tail_metric`).
function _horizon_residual_metric(coeffs, mu, a_coeffs, b_coeffs, x)
    value, derivative, second = horizon_series_triple(coeffs, mu, x)
    A = poly_value(a_coeffs, x) / x
    B = poly_value(b_coeffs, x) / x^2
    residual = second + A * derivative + B * value
    scale = abs(second) + abs(A * derivative) + abs(B * value)
    return abs(residual) / max(scale, eps(Float64))
end

function _horizon_eval_condition(coeffs, mu, x)
    value = 0.0 + 0.0im
    derivative = 0.0 + 0.0im
    value_scale = 0.0
    derivative_scale = 0.0
    power = 1.0 + 0.0im
    derivative_power = inv(ComplexF64(x))
    @inbounds for n in 0:(length(coeffs) - 1)
        value_term = coeffs[n + 1] * power
        derivative_term = (mu + n) * coeffs[n + 1] * derivative_power
        value += value_term
        derivative += derivative_term
        value_scale += abs(value_term)
        derivative_scale += abs(derivative_term)
        power *= x
        derivative_power *= x
    end
    value_condition = value_scale / max(abs(value), floatmin(Float64))
    derivative_condition = derivative_scale /
        max(abs(derivative), floatmin(Float64))
    return max(value_condition, derivative_condition)
end

function _horizon_tail_metric(coeffs, x)
    order = length(coeffs) - 1
    order <= 0 && return 0.0
    first_degree = max(1, div(order, 2))
    absolute_x = abs(x)
    power = absolute_x^first_degree
    tail = 0.0
    @inbounds for degree in first_degree:order
        tail += abs(coeffs[degree + 1]) * power
        power *= absolute_x
    end
    value = poly_value(coeffs, x)
    return tail / max(abs(value), eps(Float64))
end

# Seed point for a horizon Frobenius series with exponent mu: start from the
# radius allowed by the tail and shrink until max(tail, ODE residual, roundoff)
# <= tol. Returns (seed_x, score), the best seed if tol is never met.
function select_horizon_seed(coeffs, mu, a_coeffs, b_coeffs,
        match_x; tol=1e-14, tail_aware=false)
    radius = tail_aware ? _endpoint_tail_radius(coeffs, tol, match_x) :
        _endpoint_last_term_radius(coeffs, tol, match_x)
    seed_x = max(_MIN_STEP, min(match_x, _RADIUS_SAFETY * radius))
    best_seed = seed_x
    best_score = Inf
    for _ in 0:32
        tail = tail_aware ? _horizon_tail_metric(coeffs, seed_x) :
            abs(coeffs[end] * seed_x^(length(coeffs) - 1)) /
            max(abs(poly_value(coeffs, seed_x)), eps(Float64))
        residual = _horizon_residual_metric(coeffs, mu, a_coeffs, b_coeffs, seed_x)
        roundoff = 8eps(Float64) * _horizon_eval_condition(coeffs, mu, seed_x)
        score = max(tail, residual, roundoff)
        if score < best_score
            best_score = score
            best_seed = seed_x
        end
        score <= tol && return seed_x, score
        seed_x *= tail_aware ? 0.8 : 0.5
        seed_x < _MIN_STEP && break
    end
    return best_seed, best_score
end

@inline _tail_seed(coefficients::CoefficientSet) =
    horizon_tail(coefficients.params)

function _tail_metric(coeffs, value, h, order::Integer)
    start_k = max(0, order - 2)
    tail = 0.0
    power = abs(h)^start_k
    denom = max(abs(value), eps(Float64))
    @inbounds for k in start_k:order
        tail = max(tail, abs(coeffs[k + 1]) * power / denom)
        power *= abs(h)
    end
    return tail
end

# Step shrink factor from the tail score: the order-based estimate of the step
# that would meet tol, clamped to [0.05, 0.9]; 1 when tol is already met.
@inline function _step_factor(score, tol, order)
    score <= tol && return 1.0
    (!isfinite(score) || score <= 0 || tol <= 0) && return 0.5
    exponent = max(Int(order) - 2, 1)
    estimate = exp((log(tol) - log(score)) / exponent)
    return clamp(0.98 * estimate, 0.05, 0.9)
end

# Shrink a Taylor step from h_start until the last-terms tail metric meets tol
# (or min_step is reached); returns the best attempt and its state.
function _single_step(coeffs, direction, h_start, min_step, tol, order)
    h_abs = h_start
    best_score = Inf
    best_h = direction * h_abs
    best_value = coeffs[1]
    best_derivative = length(coeffs) > 1 ? coeffs[2] : zero(eltype(coeffs))
    for _ in 0:_MAX_REJECTIONS
        h = direction * h_abs
        value, derivative = poly_pair(coeffs, h)
        score = _tail_metric(coeffs, value, h, order)
        if score < best_score
            best_score = score
            best_h = h
            best_value = value
            best_derivative = derivative
        end
        score <= tol && break
        next_h = max(min_step, _step_factor(score, tol, order) * h_abs)
        next_h < h_abs * (1 - 8eps(Float64)) || break
        h_abs = next_h
    end
    return (score=best_score, h=best_h, value=best_value,
        derivative=best_derivative)
end

function _pair_step(coeffs1, coeffs2, direction, h_start, min_step, tol,
        order1, order2)
    h_abs = h_start
    best_score = Inf
    best_h = direction * h_abs
    best_value1 = coeffs1[1]
    best_deriv1 = length(coeffs1) > 1 ? coeffs1[2] : zero(eltype(coeffs1))
    best_value2 = coeffs2[1]
    best_deriv2 = length(coeffs2) > 1 ? coeffs2[2] : zero(eltype(coeffs2))
    for _ in 0:_MAX_REJECTIONS
        h = direction * h_abs
        value1, deriv1 = poly_pair(coeffs1, h)
        value2, deriv2 = poly_pair(coeffs2, h)
        score1 = _tail_metric(coeffs1, value1, h, order1)
        score2 = _tail_metric(coeffs2, value2, h, order2)
        score = max(score1, score2)
        if score < best_score
            best_score = score
            best_h = h
            best_value1 = value1
            best_deriv1 = deriv1
            best_value2 = value2
            best_deriv2 = deriv2
        end
        score <= tol && break
        factor = min(
            _step_factor(score1, tol, order1),
            _step_factor(score2, tol, order2),
        )
        next_h = max(min_step, factor * h_abs)
        next_h < h_abs * (1 - 8eps(Float64)) || break
        h_abs = next_h
    end
    return (
        score=best_score,
        h=best_h,
        value1=best_value1,
        deriv1=best_deriv1,
        value2=best_value2,
        deriv2=best_deriv2,
    )
end

# Highest order whose coefficients pass `series_issue`, and the reason
# for stopping (:none if all pass).
function _credible_solution_order(coeffs, order::Integer)
    max_order = min(Int(order), length(coeffs) - 1)
    scale = max(abs(coeffs[1]), max_order >= 1 ? abs(coeffs[2]) : 0.0, 1.0)
    @inbounds for n in 0:max_order
        reason = series_issue(coeffs, n + 1, scale)
        reason == :none || return n - 1, reason
    end
    return max_order, :none
end

function _record_truncation!(records, location, branch, requested, effective, reason)
    reason == :none && return records
    push!(records, Truncation(
        Symbol(location),
        Symbol(branch),
        Int(requested),
        Int(effective),
        Symbol(reason),
    ))
    return records
end

function _infinity_tail_metric(coeffs, y, order::Int=length(coeffs) - 1)
    partial = poly_value(@view(coeffs[1:(order + 1)]), y)
    term = coeffs[order + 1] * y^order
    return abs(term) / max(abs(partial), eps(Float64))
end

function _infinity_adjacent_metric(coeffs, y, order::Int=length(coeffs) - 1)
    order <= 1 && return Inf
    current = poly_value(@view(coeffs[1:(order + 1)]), y)
    previous = poly_value(@view(coeffs[1:order]), y)
    return abs(current - previous) / max(abs(current), eps(Float64))
end

# Order and trial y with the smallest max(tail, adjacent-order) metric of an
# infinity series; used as a first guess for the endpoint search.
function _asymptotic_order_and_seed_radius(coeffs::Vector{ComplexF64}, target::Float64)
    max_order = length(coeffs) - 1
    best_order = min(12, max_order)
    best_seed = _MIN_STEP
    best_score = Inf
    stale = 0
    for order in max(2, best_order):max_order
        for y in (0.002, 0.005, 0.01, 0.02, 0.04, 0.08)
            score = max(_infinity_tail_metric(coeffs, y, order),
                _infinity_adjacent_metric(coeffs, y, order))
            if score < best_score
                best_order = order
                best_seed = y
                best_score = score
                stale = 0
            else
                stale += 1
            end
        end
        stale >= 20 && order > best_order && break
    end
    best_score <= target || return best_order, best_seed
    return best_order, best_seed
end

function _infinity_residual_metric(coeffs, rho, sigma, y, A, B)
    value, derivative, second = infinity_series_triple(coeffs, rho, sigma, y)
    residual = second + A * derivative + B * value
    scale = abs(second) + abs(A * derivative) + abs(B * value)
    return abs(residual) / max(scale, eps(Float64))
end

function _infinity_order_metrics(a_coeffs, b_coeffs, coeffs, rho, sigma, y, order::Int)
    A = poly_value(a_coeffs, y) / y
    B = poly_value(b_coeffs, y) / y^4
    truncated = @view(coeffs[1:(order + 1)])
    tail = _infinity_tail_metric(coeffs, y, order)
    adjacent = _infinity_adjacent_metric(coeffs, y, order)
    residual = _infinity_residual_metric(truncated, rho, sigma, y, A, B)
    return (order=order, score=max(tail, adjacent, residual))
end

function _best_infinity_order_metrics(a_coeffs, b_coeffs, coeffs, rho, sigma, y)
    max_order = length(coeffs) - 1
    min_order = min(12, max_order)
    best = nothing
    for order in min_order:max_order
        metrics = _infinity_order_metrics(a_coeffs, b_coeffs, coeffs, rho, sigma, y, order)
        if best === nothing || metrics.score < best.score
            best = metrics
        end
    end
    return best
end

# Laurent coefficients of A (from y^-1) and B (from y^-4) rewritten for the
# scaled variable Y = y / omega.
function _scale_infinity_ab(a_coeffs, b_coeffs, omega)
    a_scaled = similar(a_coeffs)
    b_scaled = similar(b_coeffs)
    @inbounds for i in eachindex(a_coeffs)
        p = i - 2
        a_scaled[i] = a_coeffs[i] * omega^(p + 1)
    end
    @inbounds for i in eachindex(b_coeffs)
        p = i - 5
        b_scaled[i] = b_coeffs[i] * omega^(p + 2)
    end
    return a_scaled, b_scaled
end

function _scaled_y_series_triple(coeffs, rho, Sigma, Y)
    series, dseries, ddseries = poly_triple(coeffs, Y)
    prefactor = exp(Sigma / Y) * Y^rho
    log_derivative = -Sigma / Y^2 + rho / Y
    log_second = 2 * Sigma / Y^3 - rho / Y^2
    value = prefactor * series
    derivative_Y = prefactor * (dseries + log_derivative * series)
    second_Y = prefactor * (
        ddseries + 2 * log_derivative * dseries +
        (log_derivative^2 + log_second) * series
    )
    return value, derivative_Y, second_Y
end

function _scaled_y_tail_metric(coeffs, Y, order::Int)
    partial = poly_value(@view(coeffs[1:(order + 1)]), Y)
    term = coeffs[order + 1] * Y^order
    return abs(term) / max(abs(partial), eps(Float64))
end

function _scaled_y_adjacent_metric(coeffs, Y, order::Int)
    order <= 1 && return Inf
    current = poly_value(@view(coeffs[1:(order + 1)]), Y)
    previous = poly_value(@view(coeffs[1:order]), Y)
    return abs(current - previous) / max(abs(current), eps(Float64))
end

function _scaled_y_residual_metric(a_scaled, b_scaled, coeffs, rho, Sigma, Y, order::Int)
    truncated = @view coeffs[1:(order + 1)]
    value, derivative_Y, second_Y = _scaled_y_series_triple(truncated, rho, Sigma, Y)
    A = poly_value(a_scaled, Y) / Y
    B = poly_value(b_scaled, Y) / Y^4
    residual = second_Y - A * derivative_Y + B * value
    scale = abs(second_Y) + abs(A * derivative_Y) + abs(B * value)
    return abs(residual) / max(scale, eps(Float64))
end

# For one branch: the widest Y in _SCALED_Y_CANDIDATES (over all orders) whose
# max(tail, adjacent-order, residual) score meets _INFINITY_ENDPOINT_TARGET, or
# the best-scoring (order, Y) if none does.
function _best_scaled_y_branch_metric(a_scaled, b_scaled, coeffs, rho, Sigma)
    max_order = length(coeffs) - 1
    min_order = min(12, max_order)
    best = (score=Inf, order=0, Y=NaN)
    widest_accepted = nothing
    for order in min_order:max_order
        for Y in _SCALED_Y_CANDIDATES
            tail = _scaled_y_tail_metric(coeffs, Y, order)
            adjacent = _scaled_y_adjacent_metric(coeffs, Y, order)
            residual = _scaled_y_residual_metric(a_scaled, b_scaled, coeffs, rho, Sigma, Y, order)
            score = max(tail, adjacent, residual)
            if score < best.score
                best = (score=score, order=order, Y=Y)
            end
            candidate = (score=score, order=order, Y=Y)
            if score <= _INFINITY_ENDPOINT_TARGET &&
                    (widest_accepted === nothing ||
                     Y > widest_accepted.Y ||
                     (Y == widest_accepted.Y && score < widest_accepted.score))
                widest_accepted = candidate
            end
        end
    end
    return widest_accepted === nothing ? best : widest_accepted
end

# Low-frequency infinity selection in the scaled coordinate Y = y/omega (real
# omega > 0); picks for each branch the widest Y meeting
# _INFINITY_ENDPOINT_TARGET, or the best score if none does.
function _scaled_y_infinity_selection(coefficients::CoefficientSet, order::Int)
    params = coefficients.params
    omega = _param_omega(params)
    iszero(omega) && return nothing
    omega isa Real && omega > 0 || return nothing
    a_coeffs, b_coeffs = endpoint_ab_series(coefficients, :I, order)
    a_scaled, b_scaled = _scale_infinity_ab(a_coeffs, b_coeffs, omega)
    candidates = Dict{Symbol,Any}()
    for branch in (:out, :in)
        rho, sigma = infinity_exponents(params, branch)
        Sigma = sigma / omega
        guarded = infinity_guarded(
            a_scaled,
            b_scaled,
            rho,
            Sigma,
            order;
            a_first_power=-1,
            b_first_power=-4,
        )
        scaled_coeffs = guarded.coefficients
        best = _best_scaled_y_branch_metric(a_scaled, b_scaled, scaled_coeffs, rho, Sigma)
        candidates[branch] = (
            coeffs=ComplexF64.(scaled_coeffs[1:(best.order + 1)]),
            y=abs(omega) * best.Y,
            order=best.order,
            score=best.score,
            reason=guarded.truncation_reason,
        )
    end
    out = candidates[:out]
    inn = candidates[:in]
    endpoint_y = max(100eps(Float64), min(out.y, inn.y))
    return InfinitySelection(
        endpoint_y,
        out.coeffs,
        inn.coeffs,
        out.order,
        inn.order,
        max(out.score, inn.score),
        ComplexF64(omega),
        order,
        out.reason,
        inn.reason,
    )
end

function _infinity_work_order(controls)
    value = getproperty(controls, :infinity_order)
    return Int(value)
end

# Wronskian in * out' - out * in' (x-derivatives) of the two normalized infinity
# series at endpoint_y.
function infinity_seed_wronskian(
    coefficients::CoefficientSet,
    endpoint_y::Float64,
    out_coeffs,
    out_order::Int,
    in_coeffs,
    in_order::Int,
    variable_scale=1.0 + 0.0im,
)
    function endpoint_state(branch, coeffs, order)
        kind = _kind(:infinity, branch)
        rho, sigma = infinity_exponents(coefficients.params, branch)
        raw, derivative = _endpoint_pair(
            kind,
            @view(coeffs[1:(order + 1)]),
            rho,
            sigma / variable_scale,
            endpoint_y,
            variable_scale,
        )
        scale = endpoint_seed_scale(kind, coefficients.params)
        return scale * raw, scale * derivative
    end

    out_value, out_derivative = endpoint_state(:out, out_coeffs, out_order)
    in_value, in_derivative = endpoint_state(:in, in_coeffs, in_order)
    return in_value * out_derivative - out_value * in_derivative
end

function infinity_seed_wronskian(
    coefficients::CoefficientSet,
    selection::InfinitySelection,
)
    return infinity_seed_wronskian(
        coefficients,
        selection.endpoint_y,
        selection.out_coeffs,
        selection.out_order,
        selection.in_coeffs,
        selection.in_order,
        selection.variable_scale,
    )
end

# Search y = 1 - x by octaves (then intermediate phases) for an endpoint where
# both infinity series pass tail, adjacent-order and residual checks and the
# optional `certificate` is met; throws EndpointCertificateError if none.
function _select_infinity_endpoint_order(
    coefficients::CoefficientSet,
    match_x::Float64,
    order::Int;
    controls,
    certificate=nothing,
)
    a_coeffs, b_coeffs = endpoint_ab_series(coefficients, :I, order)
    rho_out, sigma_out = infinity_exponents(coefficients.params, :out)
    rho_in, sigma_in = infinity_exponents(coefficients.params, :in)
    out_guarded = infinity_guarded(a_coeffs, b_coeffs, rho_out, sigma_out, order)
    in_guarded = infinity_guarded(a_coeffs, b_coeffs, rho_in, sigma_in, order)
    out_full = out_guarded.coefficients
    in_full = in_guarded.coefficients
    seed_target = getproperty(controls, :large_frequency_expansion) ?
        max(getproperty(controls, :tolerance), 256eps(Float64)) :
        max(getproperty(controls, :tolerance), _ENDPOINT_SEED_TARGET)
    endpoint_target = seed_target
    certificate_target = max(endpoint_target, _INFINITY_ENDPOINT_TARGET)
    _, out_seed = _asymptotic_order_and_seed_radius(out_full, seed_target)
    _, in_seed = _asymptotic_order_and_seed_radius(in_full, seed_target)
    coefficient_seed = _RADIUS_SAFETY * min(out_seed, in_seed, 1.0 - match_x)
    interval_seed = _RADIUS_SAFETY * (1.0 - match_x)
    interval_seed > eps(Float64) || throw(EndpointCertificateError(
        :infinity, endpoint_target, Inf, interval_seed, order))
    initial_y = min(interval_seed, max(coefficient_seed, 0.25 * interval_seed))
    endpoint_floor = max(eps(Float64), min(_MIN_INFINITY_ENDPOINT, 0.25 * interval_seed))
    best = nothing
    for phase_pass in 1:2
        octave_scale = 1.0
        for _ in 0:64
            exhausted = true
            phase_indices = phase_pass == 1 ? (1:1) : (2:length(_INFINITY_ENDPOINT_PHASES))
            for phase_index in phase_indices
                phase = _INFINITY_ENDPOINT_PHASES[phase_index]
                y = initial_y * octave_scale * phase
                y < endpoint_floor && continue
                exhausted = false
                out_best = _best_infinity_order_metrics(
                    a_coeffs, b_coeffs, out_full, rho_out, sigma_out, y)
                in_best = _best_infinity_order_metrics(
                    a_coeffs, b_coeffs, in_full, rho_in, sigma_in, y)
                certificate_score = certificate === nothing ? 0.0 : Float64(certificate(
                    coefficients,
                    y,
                    out_full,
                    out_best.order,
                    in_full,
                    in_best.order,
                    1.0 + 0.0im,
                ))
                branch_score = max(out_best.score, in_best.score)
                score = max(branch_score, certificate_score)
                if best === nothing || score < best.score
                    best = (
                        score=score,
                        y=y,
                        out_order=out_best.order,
                        in_order=in_best.order,
                    )
                end
                if branch_score <= endpoint_target &&
                        certificate_score <= certificate_target
                    return InfinitySelection(
                        y,
                        ComplexF64.(out_full[1:(out_best.order + 1)]),
                        ComplexF64.(in_full[1:(in_best.order + 1)]),
                        out_best.order,
                        in_best.order,
                        score,
                        1.0 + 0.0im,
                        order,
                        out_guarded.truncation_reason,
                        in_guarded.truncation_reason,
                    )
                end
            end
            exhausted && break
            octave_scale *= 0.5
        end
    end
    throw(EndpointCertificateError(
        :infinity,
        endpoint_target,
        best === nothing ? Inf : best.score,
        best === nothing ? initial_y : best.y,
        order,
    ))
end

# Infinity endpoint selection: the scaled-y selection for the low-frequency
# basis, else `_select_infinity_endpoint_order`, retried at order 32 (non-large_frequency_expansion,
# requested order < 32) when it fails or gives an endpoint closer than _MIN_STEP.
function select_infinity_endpoint(
    coefficients::CoefficientSet,
    match_x::Float64;
    controls,
    certificate=nothing,
)
    order = _infinity_work_order(controls)
    if getproperty(controls, :basis) == :low_frequency_scaled_y
        scaled = _scaled_y_infinity_selection(coefficients, order)
        scaled === nothing || return scaled
    end

    can_escalate = !getproperty(controls, :large_frequency_expansion) && order < 32
    selected = try
        _select_infinity_endpoint_order(
            coefficients, match_x, order; controls, certificate)
    catch err
        can_escalate && err isa EndpointCertificateError || rethrow()
        _select_infinity_endpoint_order(
            coefficients, match_x, 32; controls, certificate)
    end
    if can_escalate && selected.requested_order == order &&
            selected.endpoint_y < _MIN_STEP
        return _select_infinity_endpoint_order(
            coefficients, match_x, 32; controls, certificate)
    end
    return selected
end

# Basis seeded from an endpoint series at seed_coord, scaled by
# `endpoint_seed_scale` and propagated to match_x. With certified_steps it
# returns nothing when a step misses the tolerance; step_budget caps the steps.
function _build_basis_series(
    coefficients::CoefficientSet,
    kind::Symbol,
    endpoint_coeffs::Vector{ComplexF64},
    exponent1,
    exponent2,
    match_x::Float64,
    seed_coord::Float64,
    order::Integer;
    controls,
    scratch=nothing,
    endpoint_variable_scale=1.0 + 0.0im,
    endpoint_requested_order=length(endpoint_coeffs) - 1,
    endpoint_effective_order=length(endpoint_coeffs) - 1,
    endpoint_reason=:none,
    certified_steps=false,
    step_budget=nothing,
)
    work = scratch === nothing ? _scratch(coefficients, order, controls) : scratch
    seed_x = _is_horizon(kind) ? seed_coord : 1.0 - seed_coord
    raw_seed, raw_dseed = _endpoint_pair(
        kind,
        endpoint_coeffs,
        exponent1,
        exponent2,
        seed_coord,
        endpoint_variable_scale,
    )
    scale = endpoint_seed_scale(kind, coefficients.params)
    current_x = seed_x
    current_value = scale * raw_seed
    current_derivative = scale * raw_dseed
    if _eikonal_enabled(controls, coefficients, work) && !_is_horizon(kind)
        hybrid = _hybrid_state(
            coefficients,
            kind,
            seed_x,
            current_value,
            current_derivative,
        match_x;
        controls=controls,
        scratch=work,
    )
        truncations = Truncation[]
        _record_truncation!(
            truncations,
            :infinity,
            kind,
            endpoint_requested_order,
            endpoint_effective_order,
            endpoint_reason,
        )
        append!(truncations, hybrid.truncations)
        return Basis(
            kind,
            scale,
            ComplexF64(exponent1),
            ComplexF64(exponent2),
            seed_x,
            seed_coord,
            current_value,
            current_derivative,
            endpoint_coeffs,
            hybrid.patches,
            hybrid.match_X,
            hybrid.match_dXdx,
            hybrid.step_count,
            hybrid.max_step,
            hybrid.min_step,
            ComplexF64(endpoint_variable_scale),
            truncations,
        )
    end
    direction = match_x >= current_x ? 1.0 : -1.0
    patches = Patch[]
    sizehint!(patches, 192)
    step_count = 0
    max_step = 0.0
    min_step = Inf
    tol = getproperty(controls, :tolerance)
    truncations = Truncation[]
    _record_truncation!(
        truncations,
        _is_horizon(kind) ? :horizon : :infinity,
        kind,
        endpoint_requested_order,
        endpoint_effective_order,
        endpoint_reason,
    )

    while direction * (match_x - current_x) > 100eps(Float64)
        step_count += 1
        step_budget !== nothing && step_count > step_budget && return nothing
        a_terms, b_terms = _local_ab!(
            work, coefficients, current_x, order)
        a_coeffs = @view work.avec[1:a_terms]
        b_coeffs = @view work.bvec[1:b_terms]
        state_scale = _state_scale(current_value, current_derivative)
        local_coeffs = _ordinary_point_coeffs!(
            work.local_coeffs,
            a_coeffs,
            b_coeffs,
            current_value / state_scale,
            current_derivative / state_scale,
            order,
            work.denom_inv,
        )
        effective_order, reason = _credible_solution_order(local_coeffs, order)
        certified_steps && effective_order < 2 && return nothing
        effective_order >= 1 || error(
            "ISEM GSN ordinary recurrence failed at x=$(current_x): " *
            "reason=$(reason), state_scale=$(max(abs(current_value), abs(current_derivative))), " *
            "c0=$(local_coeffs[1]), c1=$(local_coeffs[2])",
        )
        _record_truncation!(truncations, :ordinary, kind, order, effective_order, reason)
        finite_coeffs = @view local_coeffs[1:(effective_order + 1)]
        remaining = abs(match_x - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_cap = min(0.20, singular_distance)
        local_radius = min(radius_cap, _endpoint_effective_radius!(work.radius_buffer, finite_coeffs, radius_cap))
        min_step = 100eps(Float64) * max(1.0, abs(current_x))
        h_abs = min(remaining, max(min_step, min(remaining, _ORDINARY_STEP_SAFETY * local_radius)))
        best = _single_step(
            finite_coeffs, direction, h_abs, min_step, tol, effective_order)
        certified_steps && !(best.score <= tol) && return nothing
        next_x = current_x + best.h
        push!(patches,
            Patch(current_x, next_x, copy(finite_coeffs), state_scale))
        current_x = next_x
        current_value = state_scale * best.value
        current_derivative = state_scale * best.derivative
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))
    end
    isempty(patches) && (min_step = 0.0)
    return Basis(
        kind,
        scale,
        ComplexF64(exponent1),
        ComplexF64(exponent2),
        seed_x,
        seed_coord,
        scale * raw_seed,
        scale * raw_dseed,
        endpoint_coeffs,
        patches,
        current_value,
        current_derivative,
        step_count,
        max_step,
        isfinite(min_step) ? min_step : 0.0,
        ComplexF64(endpoint_variable_scale),
        truncations,
    )
end

# Propagation that mixes methods: a large-frequency-expansion carrier patch from infinity when
# available, eikonal segments where the coefficients vary slowly, and Taylor
# patches elsewhere. Returns a Basis with endpoint_valid = false.
function _hybrid_state(
    coefficients::CoefficientSet,
    kind::Symbol,
    seed_x::Real,
    seed_X,
    seed_dXdx,
    match_x::Real;
    controls,
    scratch=nothing,
    partial_eikonal=false,
)
    order = getproperty(controls, :ordinary_order)
    work = scratch === nothing ? _scratch(coefficients, order, controls) : scratch
    current_x = Float64(seed_x)
    target_x = Float64(match_x)
    current_value = ComplexF64(seed_X)
    current_derivative = ComplexF64(seed_dXdx)
    direction = target_x >= current_x ? 1.0 : -1.0
    patches = Union{Patch,EikonalPatch,InfinityLargeFrequencyExpansionPatch}[]
    sizehint!(patches, 48)
    step_count = 0
    max_step = 0.0
    min_step = Inf
    tol = getproperty(controls, :tolerance)
    eikonal_limit = _INFINITY_EIKONAL_LIMIT
    truncations = Truncation[]
    eikonal_attempted = false
    carrier_anchored = false

    @static if _INFINITY_LARGE_FREQUENCY_EXPANSION
        if work.large_frequency_expansion !== nothing && kind in (:infinity_in, :infinity_out)
            branch = kind == :infinity_in ? :in : :out
            scale = endpoint_seed_scale(kind, coefficients.params)
            carrier = large_frequency_expansion_infinity_anchor(coefficients, branch,
                1.0 - current_x, target_x, scale)
            if carrier !== nothing
                push!(patches, carrier.patch)
                current_x = carrier.patch.next_x
                current_value = carrier.state.X
                current_derivative = carrier.state.dXdx
                step_count = 1
                max_step = abs(carrier.patch.next_x -
                    carrier.patch.center_x)
                min_step = max_step
                carrier_anchored = true
            end
        end
    end

    if current_x != target_x && _eikonal_ready(coefficients, current_x, work)
        suffix = eikonal_split(
            coefficients,
            current_x,
            current_value,
            current_derivative,
            target_x;
            avec=work.avec,
            bvec=work.bvec,
            pq=work.pq,
            scratch=work.eikonal,
            residual_limit=eikonal_limit,
            check_limit=eikonal_limit,
            terms=partial_eikonal ? EIKONAL_TERMS : _LARGE_FREQUENCY_EXPANSION_EIKONAL_TERMS,
            large_frequency_expansion=work.large_frequency_expansion,
            allow_partial=partial_eikonal,
        )
        if suffix !== nothing
            append!(patches, suffix.patches)
            current_value = suffix.value
            current_derivative = suffix.derivative
            current_x = suffix.x
            step_count += length(suffix.patches)
            max_step = suffix.max_step
            min_step = suffix.min_step
            eikonal_attempted = true
        end
    end

    if !carrier_anchored && current_x != target_x && work.large_frequency_expansion !== nothing &&
            kind in (:infinity_in, :infinity_out)
        branch = kind == :infinity_in ? :in : :out
        carrier = large_frequency_expansion_infinity_reach(coefficients, branch, 1 - current_x,
            target_x, current_value, current_derivative)
        if carrier !== nothing
            push!(patches, carrier.patch)
            current_x = carrier.patch.next_x
            current_value = carrier.state.X
            current_derivative = carrier.state.dXdx
            step_count = 1
            max_step = abs(carrier.patch.next_x - carrier.patch.center_x)
            min_step = max_step
            if _eikonal_ready(coefficients, current_x, work)
                suffix = eikonal_split(
                    coefficients,
                    current_x,
                    current_value,
                    current_derivative,
                    target_x;
                    avec=work.avec,
                    bvec=work.bvec,
                    pq=work.pq,
                    scratch=work.eikonal,
                    residual_limit=eikonal_limit,
                    check_limit=eikonal_limit,
                    terms=partial_eikonal ? EIKONAL_TERMS : _LARGE_FREQUENCY_EXPANSION_EIKONAL_TERMS,
                    large_frequency_expansion=work.large_frequency_expansion,
                    allow_partial=partial_eikonal,
                )
                if suffix !== nothing
                    append!(patches, suffix.patches)
                    current_value = suffix.value
                    current_derivative = suffix.derivative
                    current_x = suffix.x
                    step_count += length(suffix.patches)
                    max_step = max(max_step, suffix.max_step)
                    min_step = min(min_step, suffix.min_step)
                    eikonal_attempted = true
                end
            end
        end
    end

    while direction * (target_x - current_x) > 100eps(Float64) * max(1.0, abs(current_x))
        step_count += 1
        a_terms, b_terms = _local_ab!(
            work, coefficients, current_x, order)
        a_coeffs = @view work.avec[1:a_terms]
        b_coeffs = @view work.bvec[1:b_terms]
        eikonal_ready = work.eikonal !== nothing &&
            _eikonal_indicator(a_coeffs, b_coeffs) <=
                _EIKONAL_SLOW_VARIATION
        state_scale = _state_scale(current_value, current_derivative)
        local_coeffs = _ordinary_point_coeffs!(
            work.local_coeffs,
            a_coeffs,
            b_coeffs,
            current_value / state_scale,
            current_derivative / state_scale,
            order,
            work.denom_inv,
        )
        effective_order, reason = _credible_solution_order(local_coeffs, order)
        effective_order >= 1 || error(
            "ISEM GSN hybrid recurrence failed at x=$(current_x): reason=$(reason)")
        _record_truncation!(truncations, :ordinary, kind, order, effective_order, reason)
        finite_coeffs = @view local_coeffs[1:(effective_order + 1)]
        remaining = abs(target_x - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_cap = min(0.20, singular_distance)
        local_radius = min(radius_cap,
            _endpoint_effective_radius!(work.radius_buffer, finite_coeffs, radius_cap))
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        h_abs = min(remaining,
            max(min_abs_step, min(remaining, _ORDINARY_STEP_SAFETY * local_radius)))
        best = _single_step(
            finite_coeffs, direction, h_abs, min_abs_step, tol, effective_order)
        next_x = current_x + best.h
        push!(patches,
            Patch(current_x, next_x, copy(finite_coeffs), state_scale))
        current_x = next_x
        current_value = state_scale * best.value
        current_derivative = state_scale * best.derivative
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))

        remaining_after = direction * (target_x - current_x)
        remaining_after > 100eps(Float64) || continue
        if eikonal_ready && !eikonal_attempted
            eikonal_attempted = true
            suffix = eikonal_split(
                coefficients,
                current_x,
                current_value,
                current_derivative,
                target_x;
                avec=work.avec,
                bvec=work.bvec,
                pq=work.pq,
                scratch=work.eikonal,
                residual_limit=eikonal_limit,
                check_limit=eikonal_limit,
                large_frequency_expansion=work.large_frequency_expansion,
                allow_partial=partial_eikonal,
            )
            if suffix !== nothing
                append!(patches, suffix.patches)
                current_value = suffix.value
                current_derivative = suffix.derivative
                current_x = suffix.x
                step_count += length(suffix.patches)
                max_step = max(max_step, suffix.max_step)
                min_step = min(min_step, suffix.min_step)
                direction * (target_x - current_x) <=
                    100eps(Float64) * max(1.0, abs(current_x)) && break
            end
        end
    end
    isempty(patches) && (min_step = 0.0)
    seed_coord = _is_horizon(kind) ? Float64(seed_x) : 1.0 - Float64(seed_x)
    return Basis(
        kind,
        1.0 + 0.0im,
        0.0 + 0.0im,
        0.0 + 0.0im,
        Float64(seed_x),
        seed_coord,
        ComplexF64(seed_X),
        ComplexF64(seed_dXdx),
        ComplexF64[ComplexF64(seed_X)],
        patches,
        current_value,
        current_derivative,
        step_count,
        max_step,
        isfinite(min_step) ? min_step : 0.0,
        1.0 + 0.0im,
        truncations,
        false,
    )
end

# Paired form of `_hybrid_state` with common steps for both solutions.
function _hybrid_pair(
    coefficients::CoefficientSet,
    kind1::Symbol,
    kind2::Symbol,
    seed_x::Real,
    seed_X1,
    seed_dXdx1,
    seed_X2,
    seed_dXdx2,
    match_x::Real;
    controls,
    scratch=nothing,
)
    order = getproperty(controls, :ordinary_order)
    work = scratch === nothing ? _scratch(coefficients, order, controls) : scratch
    current_x = Float64(seed_x)
    target_x = Float64(match_x)
    value1 = ComplexF64(seed_X1)
    deriv1 = ComplexF64(seed_dXdx1)
    value2 = ComplexF64(seed_X2)
    deriv2 = ComplexF64(seed_dXdx2)
    direction = target_x >= current_x ? 1.0 : -1.0
    patches1 = Union{Patch,EikonalPatch,InfinityLargeFrequencyExpansionPatch}[]
    patches2 = Union{Patch,EikonalPatch,InfinityLargeFrequencyExpansionPatch}[]
    sizehint!(patches1, 48)
    sizehint!(patches2, 48)
    local1 = work.local_coeffs
    local2 = work.local_coeffs2
    step_count = 0
    max_step = 0.0
    min_step = Inf
    tol = getproperty(controls, :tolerance)
    eikonal_limit = _INFINITY_EIKONAL_LIMIT
    truncations1 = Truncation[]
    truncations2 = Truncation[]
    eikonal_attempted = false
    carrier_anchored = false

    @static if _INFINITY_LARGE_FREQUENCY_EXPANSION
        if work.large_frequency_expansion !== nothing &&
                kind1 in (:infinity_in, :infinity_out) &&
                kind2 in (:infinity_in, :infinity_out)
            branch1 = kind1 == :infinity_in ? :in : :out
            branch2 = kind2 == :infinity_in ? :in : :out
            scale1 = endpoint_seed_scale(kind1, coefficients.params)
            scale2 = endpoint_seed_scale(kind2, coefficients.params)
            carrier = large_frequency_expansion_infinity_anchor_pair(coefficients, branch1, branch2,
                1.0 - current_x, target_x, scale1, scale2)
            if carrier !== nothing
                push!(patches1, carrier.patch1)
                push!(patches2, carrier.patch2)
                current_x = carrier.patch1.next_x
                value1 = carrier.state1.X
                deriv1 = carrier.state1.dXdx
                value2 = carrier.state2.X
                deriv2 = carrier.state2.dXdx
                step_count = 1
                max_step = abs(carrier.patch1.next_x -
                    carrier.patch1.center_x)
                min_step = max_step
                carrier_anchored = true
            end
        end
    end

    if current_x != target_x && work.large_frequency_expansion !== nothing &&
            _eikonal_ready(coefficients, current_x, work)
        suffix = eikonal_split_pair(
            coefficients,
            current_x,
            value1,
            deriv1,
            value2,
            deriv2,
            target_x;
            avec=work.avec,
            bvec=work.bvec,
            pq=work.pq,
            scratch=work.eikonal,
            residual_limit=eikonal_limit,
            check_limit=eikonal_limit,
            terms=_LARGE_FREQUENCY_EXPANSION_EIKONAL_TERMS,
            large_frequency_expansion=work.large_frequency_expansion,
        )
        if suffix !== nothing
            append!(patches1, suffix.patches1)
            append!(patches2, suffix.patches2)
            value1 = suffix.value1
            deriv1 = suffix.deriv1
            value2 = suffix.value2
            deriv2 = suffix.deriv2
            current_x = target_x
            step_count += length(suffix.patches1)
            max_step = suffix.max_step
            min_step = suffix.min_step
            eikonal_attempted = true
        end
    end

    if !carrier_anchored && current_x != target_x && work.large_frequency_expansion !== nothing &&
            kind1 in (:infinity_in, :infinity_out) &&
            kind2 in (:infinity_in, :infinity_out)
        branch1 = kind1 == :infinity_in ? :in : :out
        branch2 = kind2 == :infinity_in ? :in : :out
        carrier = large_frequency_expansion_infinity_reach_pair(coefficients, branch1, branch2,
            1 - current_x, target_x,
            value1, deriv1, value2, deriv2)
        if carrier !== nothing
            push!(patches1, carrier.patch1)
            push!(patches2, carrier.patch2)
            current_x = carrier.patch1.next_x
            value1 = carrier.state1.X
            deriv1 = carrier.state1.dXdx
            value2 = carrier.state2.X
            deriv2 = carrier.state2.dXdx
            step_count = 1
            max_step = abs(carrier.patch1.next_x -
                carrier.patch1.center_x)
            min_step = max_step
            if _eikonal_ready(coefficients, current_x, work)
                suffix = eikonal_split_pair(
                    coefficients,
                    current_x,
                    value1,
                    deriv1,
                    value2,
                    deriv2,
                    target_x;
                    avec=work.avec,
                    bvec=work.bvec,
                    pq=work.pq,
                    scratch=work.eikonal,
                    residual_limit=eikonal_limit,
                    check_limit=eikonal_limit,
                    terms=_LARGE_FREQUENCY_EXPANSION_EIKONAL_TERMS,
                    large_frequency_expansion=work.large_frequency_expansion,
                )
                if suffix !== nothing
                    append!(patches1, suffix.patches1)
                    append!(patches2, suffix.patches2)
                    value1 = suffix.value1
                    deriv1 = suffix.deriv1
                    value2 = suffix.value2
                    deriv2 = suffix.deriv2
                    current_x = target_x
                    step_count += length(suffix.patches1)
                    max_step = max(max_step, suffix.max_step)
                    min_step = min(min_step, suffix.min_step)
                    eikonal_attempted = true
                end
            end
        end
    end

    while direction * (target_x - current_x) > 100eps(Float64) * max(1.0, abs(current_x))
        step_count += 1
        a_terms, b_terms = _local_ab!(
            work, coefficients, current_x, order)
        a_coeffs = @view work.avec[1:a_terms]
        b_coeffs = @view work.bvec[1:b_terms]
        eikonal_ready = work.eikonal !== nothing &&
            _eikonal_indicator(a_coeffs, b_coeffs) <=
                _EIKONAL_SLOW_VARIATION
        scale1 = _state_scale(value1, deriv1)
        scale2 = _state_scale(value2, deriv2)
        _ordinary_point_pair!(local1, local2, a_coeffs, b_coeffs,
            value1 / scale1, deriv1 / scale1,
            value2 / scale2, deriv2 / scale2, order, work.denom_inv)
        eff1, reason1 = _credible_solution_order(local1, order)
        eff2, reason2 = _credible_solution_order(local2, order)
        eff1 >= 1 || error(
            "ISEM GSN hybrid pair failed for first solution at x=$(current_x)")
        eff2 >= 1 || error(
            "ISEM GSN hybrid pair failed for second solution at x=$(current_x)")
        _record_truncation!(truncations1, :ordinary, kind1, order, eff1, reason1)
        _record_truncation!(truncations2, :ordinary, kind2, order, eff2, reason2)
        coeffs1 = @view local1[1:(eff1 + 1)]
        coeffs2 = @view local2[1:(eff2 + 1)]
        remaining = abs(target_x - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_cap = min(0.20, singular_distance)
        r1 = _endpoint_effective_radius!(work.radius_buffer, coeffs1, radius_cap)
        r2 = _endpoint_effective_radius!(work.radius_buffer2, coeffs2, radius_cap)
        local_radius = min(radius_cap, r1, r2)
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        h_abs = min(remaining,
            max(min_abs_step, min(remaining, _ORDINARY_STEP_SAFETY * local_radius)))
        best = _pair_step(coeffs1, coeffs2, direction, h_abs,
            min_abs_step, tol, eff1, eff2)
        next_x = current_x + best.h
        push!(patches1, Patch(current_x, next_x, copy(coeffs1), scale1))
        push!(patches2, Patch(current_x, next_x, copy(coeffs2), scale2))
        current_x = next_x
        value1 = scale1 * best.value1
        deriv1 = scale1 * best.deriv1
        value2 = scale2 * best.value2
        deriv2 = scale2 * best.deriv2
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))

        remaining_after = direction * (target_x - current_x)
        remaining_after > 100eps(Float64) || continue
        if eikonal_ready && !eikonal_attempted
            eikonal_attempted = true
            suffix = eikonal_pair(
                coefficients,
                current_x,
                value1,
                deriv1,
                value2,
                deriv2,
                target_x;
                avec=work.avec,
                bvec=work.bvec,
                pq=work.pq,
                scratch=work.eikonal,
                residual_limit=eikonal_limit,
                check_limit=eikonal_limit,
                large_frequency_expansion=work.large_frequency_expansion,
            )
            if suffix !== nothing
                append!(patches1, suffix.patches1)
                append!(patches2, suffix.patches2)
                value1 = suffix.value1
                deriv1 = suffix.deriv1
                value2 = suffix.value2
                deriv2 = suffix.deriv2
                current_x = target_x
                step_count += length(suffix.patches1)
                max_step = max(max_step, suffix.max_step)
                min_step = min(min_step, suffix.min_step)
                break
            end
        end
    end
    isempty(patches1) && (min_step = 0.0)
    seed_coord = _is_horizon(kind1) ? Float64(seed_x) : 1.0 - Float64(seed_x)
    common_min = isfinite(min_step) ? min_step : 0.0
    basis1 = Basis(
        kind1, 1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im,
        Float64(seed_x), seed_coord, ComplexF64(seed_X1), ComplexF64(seed_dXdx1),
        ComplexF64[ComplexF64(seed_X1)], patches1, value1, deriv1,
        step_count, max_step, common_min, 1.0 + 0.0im, truncations1, false)
    basis2 = Basis(
        kind2, 1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im,
        Float64(seed_x), seed_coord, ComplexF64(seed_X2), ComplexF64(seed_dXdx2),
        ComplexF64[ComplexF64(seed_X2)], patches2, value2, deriv2,
        step_count, max_step, common_min, 1.0 + 0.0im, truncations2, false)
    return basis1, basis2
end

"""
    iterate_from_state(coefficients, kind, seed_x, seed_X, seed_dXdx, match_x; controls, scratch=nothing, allow_eikonal=false, allow_partial=false)

Propagate the state `(seed_X, seed_dXdx)` from `seed_x` to `match_x` and return
a `Basis` of kind `kind` without an endpoint series. Eikonal segments are
used when enabled for these parameters and either `kind` is not a horizon kind
or `allow_eikonal` is set; Taylor patches are used otherwise.
"""
function iterate_from_state(
    coefficients::CoefficientSet,
    kind::Symbol,
    seed_x::Real,
    seed_X,
    seed_dXdx,
    match_x::Real;
    controls,
    scratch=nothing,
    allow_eikonal=false,
    allow_partial=false,
)
    if _eikonal_enabled(controls, coefficients, scratch) &&
            (!_is_horizon(kind) || allow_eikonal)
        return _hybrid_state(
            coefficients,
            kind,
            seed_x,
            seed_X,
            seed_dXdx,
            match_x;
            controls=controls,
            scratch=scratch,
            partial_eikonal=allow_partial,
        )
    end
    order = getproperty(controls, :ordinary_order)
    work = scratch === nothing ? _scratch(coefficients, order, controls) : scratch
    current_x = Float64(seed_x)
    target_x = Float64(match_x)
    current_value = ComplexF64(seed_X)
    current_derivative = ComplexF64(seed_dXdx)
    direction = target_x >= current_x ? 1.0 : -1.0
    patches = Patch[]
    sizehint!(patches, 192)
    step_count = 0
    max_step = 0.0
    min_step = Inf
    tol = getproperty(controls, :tolerance)
    truncations = Truncation[]

    while direction * (target_x - current_x) > 100eps(Float64) * max(1.0, abs(current_x))
        step_count += 1
        a_terms, b_terms = _local_ab!(
            work, coefficients, current_x, order)
        a_coeffs = @view work.avec[1:a_terms]
        b_coeffs = @view work.bvec[1:b_terms]
        state_scale = _state_scale(current_value, current_derivative)
        local_coeffs = _ordinary_point_coeffs!(
            work.local_coeffs,
            a_coeffs,
            b_coeffs,
            current_value / state_scale,
            current_derivative / state_scale,
            order,
            work.denom_inv,
        )
        effective_order, reason = _credible_solution_order(local_coeffs, order)
        effective_order >= 1 || error(
            "ISEM GSN state propagation failed at x=$(current_x): " *
            "reason=$(reason), state_scale=$(max(abs(current_value), abs(current_derivative))), " *
            "c0=$(local_coeffs[1]), c1=$(local_coeffs[2])",
        )
        _record_truncation!(truncations, :ordinary, kind, order, effective_order, reason)
        finite_coeffs = @view local_coeffs[1:(effective_order + 1)]
        remaining = abs(target_x - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_cap = min(0.20, singular_distance)
        local_radius = min(radius_cap, _endpoint_effective_radius!(work.radius_buffer, finite_coeffs, radius_cap))
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        h_abs = min(remaining, max(min_abs_step, min(remaining, _ORDINARY_STEP_SAFETY * local_radius)))
        best = _single_step(finite_coeffs, direction, h_abs,
            min_abs_step, tol, effective_order)
        next_x = current_x + best.h
        push!(patches,
            Patch(current_x, next_x, copy(finite_coeffs), state_scale))
        current_x = next_x
        current_value = state_scale * best.value
        current_derivative = state_scale * best.derivative
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))
    end
    isempty(patches) && (min_step = 0.0)
    seed_coord = _is_horizon(kind) ? Float64(seed_x) : 1.0 - Float64(seed_x)
    return Basis(
        kind,
        1.0 + 0.0im,
        0.0 + 0.0im,
        0.0 + 0.0im,
        Float64(seed_x),
        seed_coord,
        ComplexF64(seed_X),
        ComplexF64(seed_dXdx),
        ComplexF64[ComplexF64(seed_X)],
        patches,
        current_value,
        current_derivative,
        step_count,
        max_step,
        isfinite(min_step) ? min_step : 0.0,
        1.0 + 0.0im,
        truncations,
        false,
    )
end

@inline _scaled_finite(value) =
    isfinite(real(value)) && isfinite(imag(value))

# State (X, dXdx) times exp(log_scale), stored as a unit-size pair and the
# accumulated log of its magnitude.
function logscaled_state(X, dXdx, log_scale::Real=0.0)
    value = ComplexF64(X)
    derivative = ComplexF64(dXdx)
    _scaled_finite(value) && _scaled_finite(derivative) ||
        error("ISEM GSN log-scaled state is nonfinite")
    scale = max(abs(value), abs(derivative))
    isfinite(scale) && scale > 0 ||
        error("ISEM GSN log-scaled state is zero or nonfinite")
    total_log_scale = Float64(log_scale) + log(scale)
    isfinite(total_log_scale) ||
        error("ISEM GSN log-scaled state scale is nonfinite")
    return LogScaledState(
        ComplexF64(value / scale),
        ComplexF64(derivative / scale),
        total_log_scale,
    )
end

function _materialize_logscaled_component(component, log_scale, phase)
    iszero(component) && return ComplexF64(0)
    component_log = log_scale + log(abs(component))
    component_log <= log(floatmax(Float64)) ||
        error("ISEM GSN log-scaled state exceeds Float64 range")
    value = exp(component_log) * cis(angle(component) + phase)
    _scaled_finite(value) ||
        error("ISEM GSN log-scaled state materialization is nonfinite")
    return ComplexF64(value)
end

# Named tuple (X, dXdx) equal to `scale * exp(state.log_scale)` times the stored
# components; errors when the result is not finite.
function materialize_logscaled_state(
    state::LogScaledState;
    scale=1.0 + 0.0im,
)
    applied = ComplexF64(scale)
    _scaled_finite(applied) ||
        error("ISEM GSN log-scaled materialization scale is nonfinite")
    iszero(applied) && return (X=ComplexF64(0), dXdx=ComplexF64(0))
    applied_log = state.log_scale + log(abs(applied))
    applied_phase = angle(applied)
    return (
        X=_materialize_logscaled_component(
            state.X, applied_log, applied_phase),
        dXdx=_materialize_logscaled_component(
            state.dXdx, applied_log, applied_phase),
    )
end

# Taylor propagation of a log-scaled state, renormalizing after every patch so
# that states beyond the Float64 range stay representable.
function iterate_logscaled_from_state(
    coefficients::CoefficientSet,
    kind::Symbol,
    seed_x::Real,
    seed_state::LogScaledState,
    match_x::Real;
    controls,
    scratch=nothing,
)
    order = getproperty(controls, :ordinary_order)
    work = scratch === nothing ? _scratch(coefficients, order, controls) : scratch
    current_x = Float64(seed_x)
    target_x = Float64(match_x)
    state = seed_state
    direction = target_x >= current_x ? 1.0 : -1.0
    patches = LogScaledPatch[]
    sizehint!(patches, 192)
    step_count = 0
    max_step = 0.0
    min_step = Inf
    tol = getproperty(controls, :tolerance)
    truncations = Truncation[]

    while direction * (target_x - current_x) >
            100eps(Float64) * max(1.0, abs(current_x))
        step_count += 1
        a_terms, b_terms = _local_ab!(work, coefficients, current_x, order)
        a_coeffs = @view work.avec[1:a_terms]
        b_coeffs = @view work.bvec[1:b_terms]
        local_coeffs = _ordinary_point_coeffs!(
            work.local_coeffs,
            a_coeffs,
            b_coeffs,
            state.X,
            state.dXdx,
            order,
            work.denom_inv,
        )
        effective_order, reason = _credible_solution_order(local_coeffs, order)
        effective_order >= 1 || error(
            "ISEM GSN log-scaled propagation failed at x=$(current_x): " *
            "reason=$(reason)",
        )
        _record_truncation!(
            truncations, :ordinary, kind, order, effective_order, reason)
        finite_coeffs = @view local_coeffs[1:(effective_order + 1)]
        remaining = abs(target_x - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_cap = min(0.20, singular_distance)
        local_radius = min(radius_cap, _endpoint_effective_radius!(
            work.radius_buffer, finite_coeffs, radius_cap))
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        h_abs = min(remaining, max(min_abs_step,
            min(remaining, _ORDINARY_STEP_SAFETY * local_radius)))
        best = _single_step(
            finite_coeffs, direction, h_abs, min_abs_step, tol, effective_order)
        _scaled_finite(best.value) && _scaled_finite(best.derivative) ||
            error("ISEM GSN log-scaled propagation produced a nonfinite state")
        next_x = current_x + best.h
        push!(patches, LogScaledPatch(
            current_x,
            next_x,
            copy(finite_coeffs),
            state.log_scale,
        ))
        state = logscaled_state(
            best.value, best.derivative, state.log_scale)
        current_x = next_x
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))
    end
    isempty(patches) && (min_step = 0.0)
    return LogScaledBasis(
        kind,
        Float64(seed_x),
        seed_state,
        patches,
        state,
        step_count,
        max_step,
        isfinite(min_step) ? min_step : 0.0,
        truncations,
    )
end

# Small-frequency propagation in the scaled coordinate y = (1 - x)/|omega| from
# seed_x toward target_y, stopping early when stop_check(x, X, dXdx) is true or
# a step fails its checks. Returns a basis with ScaledPatch patches, or
# nothing if no step was taken.
function iterate_scaled_y(
    coefficients::CoefficientSet,
    kind::Symbol,
    seed_x::Real,
    seed_X,
    seed_dXdx;
    target_y=nothing,
    order::Integer=_SMALL_FREQUENCY_EXPANSION_BRIDGE_ORDER,
    tolerance::Real=_SMALL_FREQUENCY_EXPANSION_BRIDGE_TOLERANCE,
    stop_check=nothing,
)
    omega_scale = abs(_param_omega(coefficients.params))
    omega_scale > 0 || return nothing
    start_y = (1.0 - Float64(seed_x)) / omega_scale
    final_y = target_y === nothing ?
        max(_SMALL_FREQUENCY_EXPANSION_BRIDGE_TARGET_Y, _MIN_STEP / omega_scale) :
        Float64(target_y)
    final_y == start_y && return nothing
    direction = final_y > start_y ? 1.0 : -1.0
    work = IterationScratch(coefficients, order; ordinary_dd=false)
    current_y = start_y
    current_value = ComplexF64(seed_X)
    current_derivative_y = ComplexF64(-omega_scale * seed_dXdx)
    patches = ScaledPatch[]
    sizehint!(patches, 128)
    max_step = 0.0
    min_step = Inf

    while direction * (final_y - current_y) > 100eps(Float64) * max(1.0, current_y)
        current_x = 1.0 - omega_scale * current_y
        a_terms, b_terms = ordinary_scaled_y_ab_series!(
            work.avec, work.bvec, work.pq, coefficients, current_y, omega_scale, order)
        a_coeffs = @view work.avec[1:a_terms]
        b_coeffs = @view work.bvec[1:b_terms]
        all(_scaled_finite, a_coeffs) && all(_scaled_finite, b_coeffs) || break

        state_scale = _state_scale(current_value, current_derivative_y)
        isfinite(state_scale) || break
        local_coeffs = _ordinary_point_coeffs!(
            work.local_coeffs,
            a_coeffs,
            b_coeffs,
            current_value / state_scale,
            current_derivative_y / state_scale,
            order,
            work.denom_inv,
        )
        all(_scaled_finite, local_coeffs) || break

        radius_cap = min(current_y, inv(omega_scale) - current_y)
        radius_cap > 0 && isfinite(radius_cap) || break
        radius = _endpoint_effective_radius!(
            work.radius_buffer, local_coeffs, radius_cap)
        isfinite(radius) && radius > 0 || break
        remaining = direction * (final_y - current_y)
        h_start = min(remaining, _SMALL_FREQUENCY_EXPANSION_BRIDGE_STEP_SAFETY * radius)
        min_step_y = 100eps(Float64) * max(1.0, current_y)
        h_start > min_step_y || break
        step = _single_step(local_coeffs, direction, h_start,
            min_step_y, Float64(tolerance), order)
        step.score <= tolerance || break
        next_y = current_y + step.h
        direction * (next_y - current_y) > 0 || break
        next_value = ComplexF64(state_scale * step.value)
        next_derivative_y = ComplexF64(state_scale * step.derivative)
        _scaled_finite(next_value) && _scaled_finite(next_derivative_y) || break

        next_x = 1.0 - omega_scale * next_y
        direction * (current_x - next_x) > 0 || break
        push!(patches, ScaledPatch(
            current_x,
            next_x,
            current_y,
            omega_scale,
            copy(local_coeffs),
            state_scale,
        ))
        x_step = abs(next_x - current_x)
        max_step = max(max_step, x_step)
        min_step = min(min_step, x_step)
        current_y = next_y
        current_value = next_value
        current_derivative_y = next_derivative_y
        if stop_check !== nothing && stop_check(
                next_x,
                current_value,
                ComplexF64(-current_derivative_y / omega_scale),
            )
            break
        end
    end
    if direction > 0 && final_y - current_y > 100eps(Float64) * max(1.0, current_y)
        error("ISEM GSN scaled-y propagation did not reach the matching point")
    end
    isempty(patches) && return nothing

    return Basis(
        kind,
        1.0 + 0.0im,
        0.0 + 0.0im,
        0.0 + 0.0im,
        Float64(seed_x),
        start_y,
        ComplexF64(seed_X),
        ComplexF64(seed_dXdx),
        ComplexF64[ComplexF64(seed_X)],
        patches,
        current_value,
        ComplexF64(-current_derivative_y / omega_scale),
        length(patches),
        max_step,
        isfinite(min_step) ? min_step : 0.0,
        ComplexF64(-omega_scale),
        Truncation[],
        false,
    )
end

# Two states propagated from seed_x to match_x with common steps (the hybrid
# march for infinity kinds when eikonal patches are enabled); returns
# (basis1, basis2) with endpoint_valid = false.
function iterate_pair_from_state(
    coefficients::CoefficientSet,
    kind1::Symbol,
    kind2::Symbol,
    seed_x::Real,
    seed_X1,
    seed_dXdx1,
    seed_X2,
    seed_dXdx2,
    match_x::Real;
    controls,
    scratch=nothing,
)
    if _eikonal_enabled(controls, coefficients, scratch) && !_is_horizon(kind1)
        return _hybrid_pair(
            coefficients,
            kind1,
            kind2,
            seed_x,
            seed_X1,
            seed_dXdx1,
            seed_X2,
            seed_dXdx2,
            match_x;
            controls=controls,
            scratch=scratch,
        )
    end
    order = getproperty(controls, :ordinary_order)
    work = scratch === nothing ? _scratch(coefficients, order, controls) : scratch
    current_x = Float64(seed_x)
    target_x = Float64(match_x)
    value1 = ComplexF64(seed_X1)
    deriv1 = ComplexF64(seed_dXdx1)
    value2 = ComplexF64(seed_X2)
    deriv2 = ComplexF64(seed_dXdx2)
    direction = target_x >= current_x ? 1.0 : -1.0
    patches1 = Patch[]
    patches2 = Patch[]
    sizehint!(patches1, 192)
    sizehint!(patches2, 192)
    local1 = work.local_coeffs
    local2 = work.local_coeffs2
    radius1 = work.radius_buffer
    radius2 = work.radius_buffer2
    step_count = 0
    max_step = 0.0
    min_step = Inf
    tol = getproperty(controls, :tolerance)
    truncations1 = Truncation[]
    truncations2 = Truncation[]

    while direction * (target_x - current_x) > 100eps(Float64) * max(1.0, abs(current_x))
        step_count += 1
        a_terms, b_terms = _local_ab!(
            work, coefficients, current_x, order)
        a_coeffs = @view work.avec[1:a_terms]
        b_coeffs = @view work.bvec[1:b_terms]
        state_scale1 = _state_scale(value1, deriv1)
        state_scale2 = _state_scale(value2, deriv2)
        _ordinary_point_pair!(local1, local2, a_coeffs, b_coeffs,
            value1 / state_scale1, deriv1 / state_scale1,
            value2 / state_scale2, deriv2 / state_scale2, order, work.denom_inv)
        eff1, reason1 = _credible_solution_order(local1, order)
        eff2, reason2 = _credible_solution_order(local2, order)
        eff1 >= 1 || error(
            "ISEM GSN paired state recurrence failed for first solution at x=$(current_x): " *
            "reason=$(reason1), state_scale=$(max(abs(value1), abs(deriv1))), " *
            "c0=$(local1[1]), c1=$(local1[2])",
        )
        eff2 >= 1 || error(
            "ISEM GSN paired state recurrence failed for second solution at x=$(current_x): " *
            "reason=$(reason2), state_scale=$(max(abs(value2), abs(deriv2))), " *
            "c0=$(local2[1]), c1=$(local2[2])",
        )
        _record_truncation!(truncations1, :ordinary, kind1, order, eff1, reason1)
        _record_truncation!(truncations2, :ordinary, kind2, order, eff2, reason2)
        coeffs1 = @view local1[1:(eff1 + 1)]
        coeffs2 = @view local2[1:(eff2 + 1)]
        remaining = abs(target_x - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_cap = min(0.20, singular_distance)
        r1 = _endpoint_effective_radius!(radius1, coeffs1, radius_cap)
        r2 = _endpoint_effective_radius!(radius2, coeffs2, radius_cap)
        local_radius = min(radius_cap, r1, r2)
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        h_abs = min(remaining, max(min_abs_step, min(remaining, _ORDINARY_STEP_SAFETY * local_radius)))
        best = _pair_step(coeffs1, coeffs2, direction, h_abs,
            min_abs_step, tol, eff1, eff2)
        next_x = current_x + best.h
        push!(patches1,
            Patch(current_x, next_x, copy(coeffs1), state_scale1))
        push!(patches2,
            Patch(current_x, next_x, copy(coeffs2), state_scale2))
        current_x = next_x
        value1 = state_scale1 * best.value1
        deriv1 = state_scale1 * best.deriv1
        value2 = state_scale2 * best.value2
        deriv2 = state_scale2 * best.deriv2
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))
    end
    isempty(patches1) && (min_step = 0.0)
    seed_coord = _is_horizon(kind1) ? Float64(seed_x) : 1.0 - Float64(seed_x)
    common_min_step = isfinite(min_step) ? min_step : 0.0
    basis1 = Basis(
        kind1,
        1.0 + 0.0im,
        0.0 + 0.0im,
        0.0 + 0.0im,
        Float64(seed_x),
        seed_coord,
        ComplexF64(seed_X1),
        ComplexF64(seed_dXdx1),
        ComplexF64[ComplexF64(seed_X1)],
        patches1,
        value1,
        deriv1,
        step_count,
        max_step,
        common_min_step,
        1.0 + 0.0im,
        truncations1,
        false,
    )
    basis2 = Basis(
        kind2,
        1.0 + 0.0im,
        0.0 + 0.0im,
        0.0 + 0.0im,
        Float64(seed_x),
        seed_coord,
        ComplexF64(seed_X2),
        ComplexF64(seed_dXdx2),
        ComplexF64[ComplexF64(seed_X2)],
        patches2,
        value2,
        deriv2,
        step_count,
        max_step,
        common_min_step,
        1.0 + 0.0im,
        truncations2,
        false,
    )
    return basis1, basis2
end

# Value and x-derivative of a horizon local solution at x, including the
# log(x) term of a threshold (logarithmic) solution when present.
function _horizon_solution_pair(solution, x)
    value, derivative = horizon_series_pair(
        solution.coefficients, solution.exponent, x)
    if hasproperty(solution, :log_coefficients) &&
            !isempty(solution.log_coefficients) &&
            !iszero(solution.log_coefficient)
        log_value, log_derivative = horizon_series_pair(
            solution.log_coefficients, solution.log_exponent, x)
        logarithm = log(x)
        value += solution.log_coefficient * logarithm * log_value
        derivative += solution.log_coefficient * (
            logarithm * log_derivative + log_value / x)
    end
    return ComplexF64(value), ComplexF64(derivative)
end

# Taylor march of both horizon solutions from the common seed_x to match_x with
# shared steps; the plain route of `iterate_pair_from_zero`.
function _build_pair_from_horizon_seed(
    coefficients::CoefficientSet,
    solution1,
    solution2,
    branch1::Symbol,
    branch2::Symbol,
    seed_x::Float64,
    match_x::Float64;
    controls,
    scratch=nothing,
)
    order = getproperty(controls, :ordinary_order)
    work = scratch === nothing ? _scratch(coefficients, order, controls) : scratch
    kind1 = _kind(:horizon, branch1)
    kind2 = _kind(:horizon, branch2)
    raw1, draw1 = _horizon_solution_pair(solution1, seed_x)
    raw2, draw2 = _horizon_solution_pair(solution2, seed_x)
    scale1 = endpoint_seed_scale(kind1, coefficients.params)
    scale2 = endpoint_seed_scale(kind2, coefficients.params)
    current_x = seed_x
    value1 = scale1 * raw1
    deriv1 = scale1 * draw1
    value2 = scale2 * raw2
    deriv2 = scale2 * draw2
    direction = match_x >= current_x ? 1.0 : -1.0
    patches1 = Patch[]
    patches2 = Patch[]
    sizehint!(patches1, 192)
    sizehint!(patches2, 192)
    local1 = work.local_coeffs
    local2 = work.local_coeffs2
    radius1 = work.radius_buffer
    radius2 = work.radius_buffer2
    step_count = 0
    max_step = 0.0
    min_step = Inf
    tol = getproperty(controls, :tolerance)
    truncations1 = Truncation[]
    truncations2 = Truncation[]
    _record_truncation!(truncations1, :horizon, kind1,
        solution1.requested_order, solution1.effective_order, solution1.truncation_reason)
    _record_truncation!(truncations2, :horizon, kind2,
        solution2.requested_order, solution2.effective_order, solution2.truncation_reason)

    while direction * (match_x - current_x) > 100eps(Float64)
        step_count += 1
        a_terms, b_terms = _local_ab!(
            work, coefficients, current_x, order)
        a_coeffs = @view work.avec[1:a_terms]
        b_coeffs = @view work.bvec[1:b_terms]
        state_scale1 = _state_scale(value1, deriv1)
        state_scale2 = _state_scale(value2, deriv2)
        _ordinary_point_pair!(local1, local2, a_coeffs, b_coeffs,
            value1 / state_scale1, deriv1 / state_scale1,
            value2 / state_scale2, deriv2 / state_scale2, order, work.denom_inv)
        eff1, reason1 = _credible_solution_order(local1, order)
        eff2, reason2 = _credible_solution_order(local2, order)
        eff1 >= 1 || error(
            "ISEM GSN paired horizon recurrence failed for first solution at x=$(current_x): " *
            "reason=$(reason1), state_scale=$(max(abs(value1), abs(deriv1))), " *
            "c0=$(local1[1]), c1=$(local1[2])",
        )
        eff2 >= 1 || error(
            "ISEM GSN paired horizon recurrence failed for second solution at x=$(current_x): " *
            "reason=$(reason2), state_scale=$(max(abs(value2), abs(deriv2))), " *
            "c0=$(local2[1]), c1=$(local2[2])",
        )
        _record_truncation!(truncations1, :ordinary, kind1, order, eff1, reason1)
        _record_truncation!(truncations2, :ordinary, kind2, order, eff2, reason2)
        coeffs1 = @view local1[1:(eff1 + 1)]
        coeffs2 = @view local2[1:(eff2 + 1)]
        remaining = abs(match_x - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_cap = min(0.20, singular_distance)
        r1 = _endpoint_effective_radius!(radius1, coeffs1, radius_cap)
        r2 = _endpoint_effective_radius!(radius2, coeffs2, radius_cap)
        local_radius = min(radius_cap, r1, r2)
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        h_abs = min(remaining, max(min_abs_step, min(remaining, _ORDINARY_STEP_SAFETY * local_radius)))
        best = _pair_step(coeffs1, coeffs2, direction, h_abs,
            min_abs_step, tol, eff1, eff2)
        next_x = current_x + best.h
        push!(patches1,
            Patch(current_x, next_x, copy(coeffs1), state_scale1))
        push!(patches2,
            Patch(current_x, next_x, copy(coeffs2), state_scale2))
        current_x = next_x
        value1 = state_scale1 * best.value1
        deriv1 = state_scale1 * best.deriv1
        value2 = state_scale2 * best.value2
        deriv2 = state_scale2 * best.deriv2
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))
    end
    isempty(patches1) && (min_step = 0.0)
    common_min_step = isfinite(min_step) ? min_step : 0.0
    basis1 = Basis(
        kind1,
        scale1,
        ComplexF64(solution1.exponent),
        0.0 + 0.0im,
        seed_x,
        seed_x,
        scale1 * raw1,
        scale1 * draw1,
        solution1.coefficients,
        patches1,
        value1,
        deriv1,
        step_count,
        max_step,
        common_min_step,
        1.0 + 0.0im,
        truncations1,
    )
    basis2 = Basis(
        kind2,
        scale2,
        ComplexF64(solution2.exponent),
        0.0 + 0.0im,
        seed_x,
        seed_x,
        scale2 * raw2,
        scale2 * draw2,
        solution2.coefficients,
        patches2,
        value2,
        deriv2,
        step_count,
        max_step,
        common_min_step,
        1.0 + 0.0im,
        truncations2,
    )
    return _with_horizon_log(basis1, solution1), _with_horizon_log(basis2, solution2)
end

# Attach the logarithmic part of a threshold horizon solution to its basis, so that
# evaluation between the horizon and the seed includes log(x) X_log.
function _with_horizon_log(basis::Basis, solution)
    hasproperty(solution, :log_coefficients) &&
        !isempty(solution.log_coefficients) &&
        !iszero(solution.log_coefficient) || return basis
    return Basis(
        basis.kind, basis.scale, basis.exponent1, basis.exponent2,
        basis.seed_x, basis.seed_coord, basis.seed_X, basis.seed_dXdx,
        basis.endpoint_coeffs, basis.patches, basis.match_X, basis.match_dXdx,
        basis.step_count, basis.max_step, basis.min_step,
        basis.endpoint_variable_scale, basis.truncations, basis.endpoint_valid,
        ComplexF64.(solution.log_coefficient .* solution.log_coefficients),
        ComplexF64(solution.log_exponent),
    )
end

# Logarithmic part log(x) X_log of a threshold horizon basis and its x-derivative.
function _horizon_log_pair(basis::Basis, x)
    log_value, log_derivative = horizon_series_pair(
        basis.endpoint_log_coeffs, basis.endpoint_log_exponent, x)
    logarithm = log(x)
    return logarithm * log_value, logarithm * log_derivative + log_value / x
end

function _patch_range(patches, seed_x)
    isempty(patches) && return 0.0, 0.0
    maximum_step = 0.0
    minimum_step = Inf
    current_x = seed_x
    for patch in patches
        step = abs(patch.next_x - current_x)
        maximum_step = max(maximum_step, step)
        minimum_step = min(minimum_step, step)
        current_x = patch.next_x
    end
    return maximum_step, isfinite(minimum_step) ? minimum_step : 0.0
end

# Horizon large-frequency route: one horizon-series patch from seed_x, then an
# eikonal path (or a Taylor march if that fails) to match_x; nothing when no
# patch is available. Called only when ISEM_HORIZON_LARGE_FREQUENCY_EXPANSION=1.
# `_large_frequency_expansion_horizon_basis_pair` is the two-branch version.
function _large_frequency_expansion_horizon_basis(
    coefficients::CoefficientSet,
    solution,
    branch::Symbol,
    seed_x::Float64,
    match_x::Float64,
    work::IterationScratch,
    controls,
)
    kind = _kind(:horizon, branch)
    raw, draw = _endpoint_pair(
        kind, solution.coefficients, solution.exponent,
        0.0 + 0.0im, seed_x)
    scale = endpoint_seed_scale(kind, coefficients.params)
    seed_value = scale * raw
    seed_derivative = scale * draw
    candidate = large_frequency_expansion_horizon_single(
        coefficients, branch, seed_x, scale, match_x)
    candidate === nothing && return nothing
    PatchType = Union{HorizonLargeFrequencyExpansionPatch,EikonalPatch,Patch}
    patches = PatchType[candidate.patch]
    if candidate.target < match_x * (1 - 100eps(Float64))
        suffix = eikonal_path(
            coefficients, candidate.target,
            candidate.state.X, candidate.state.dXdx, match_x;
            avec=work.avec, bvec=work.bvec, pq=work.pq,
            scratch=work.eikonal, large_frequency_expansion=work.large_frequency_expansion,
            residual_limit=_INFINITY_EIKONAL_LIMIT,
            check_limit=_INFINITY_EIKONAL_LIMIT)
        if suffix === nothing
            basis = iterate_from_state(
                coefficients, kind, candidate.target,
                candidate.state.X, candidate.state.dXdx, match_x;
                controls, scratch=work)
            append!(patches, basis.patches)
            match_value = basis.match_X
            match_derivative = basis.match_dXdx
        else
            append!(patches, suffix.patches)
            match_value = suffix.value
            match_derivative = suffix.derivative
        end
    else
        match_value = candidate.state.X
        match_derivative = candidate.state.dXdx
    end
    max_step, min_step = _patch_range(patches, seed_x)
    truncations = Truncation[]
    _record_truncation!(truncations, :horizon, kind,
        solution.requested_order, solution.effective_order,
        solution.truncation_reason)
    return Basis(
        kind, scale, ComplexF64(solution.exponent), 0.0 + 0.0im,
        seed_x, seed_x, seed_value, seed_derivative,
        solution.coefficients, patches,
        match_value, match_derivative, length(patches),
        max_step, min_step, 1.0 + 0.0im, truncations)
end

function _large_frequency_expansion_horizon_basis_pair(
    coefficients::CoefficientSet,
    solution1,
    solution2,
    branch1::Symbol,
    branch2::Symbol,
    seed_x::Float64,
    match_x::Float64,
    work::IterationScratch,
    controls,
)
    kind1 = _kind(:horizon, branch1)
    kind2 = _kind(:horizon, branch2)
    raw1, draw1 = _endpoint_pair(
        kind1, solution1.coefficients, solution1.exponent,
        0.0 + 0.0im, seed_x)
    raw2, draw2 = _endpoint_pair(
        kind2, solution2.coefficients, solution2.exponent,
        0.0 + 0.0im, seed_x)
    scale1 = endpoint_seed_scale(kind1, coefficients.params)
    scale2 = endpoint_seed_scale(kind2, coefficients.params)
    value1 = scale1 * raw1
    deriv1 = scale1 * draw1
    value2 = scale2 * raw2
    deriv2 = scale2 * draw2
    candidate = large_frequency_expansion_horizon_pair(coefficients,
        branch1, scale1, branch2, scale2, seed_x, match_x)
    candidate === nothing && return nothing
    PatchType = Union{HorizonLargeFrequencyExpansionPatch,EikonalPatch,Patch}
    patches1 = PatchType[candidate.patch1]
    patches2 = PatchType[candidate.patch2]
    if candidate.target < match_x * (1 - 100eps(Float64))
        suffix = eikonal_pair(
            coefficients, candidate.target,
            candidate.state1.X, candidate.state1.dXdx,
            candidate.state2.X, candidate.state2.dXdx, match_x;
            avec=work.avec, bvec=work.bvec, pq=work.pq,
            scratch=work.eikonal, large_frequency_expansion=work.large_frequency_expansion,
            residual_limit=_INFINITY_EIKONAL_LIMIT,
            check_limit=_INFINITY_EIKONAL_LIMIT)
        if suffix === nothing
            basis1, basis2 = iterate_pair_from_state(
                coefficients, kind1, kind2, candidate.target,
                candidate.state1.X, candidate.state1.dXdx,
                candidate.state2.X, candidate.state2.dXdx, match_x;
                controls, scratch=work)
            append!(patches1, basis1.patches)
            append!(patches2, basis2.patches)
            match_value1 = basis1.match_X
            match_deriv1 = basis1.match_dXdx
            match_value2 = basis2.match_X
            match_deriv2 = basis2.match_dXdx
        else
            append!(patches1, suffix.patches1)
            append!(patches2, suffix.patches2)
            match_value1 = suffix.value1
            match_deriv1 = suffix.deriv1
            match_value2 = suffix.value2
            match_deriv2 = suffix.deriv2
        end
    else
        match_value1 = candidate.state1.X
        match_deriv1 = candidate.state1.dXdx
        match_value2 = candidate.state2.X
        match_deriv2 = candidate.state2.dXdx
    end
    max_step1, min_step1 = _patch_range(patches1, seed_x)
    max_step2, min_step2 = _patch_range(patches2, seed_x)
    truncations1 = Truncation[]
    truncations2 = Truncation[]
    _record_truncation!(truncations1, :horizon, kind1,
        solution1.requested_order, solution1.effective_order,
        solution1.truncation_reason)
    _record_truncation!(truncations2, :horizon, kind2,
        solution2.requested_order, solution2.effective_order,
        solution2.truncation_reason)
    basis1 = Basis(
        kind1, scale1, ComplexF64(solution1.exponent), 0.0 + 0.0im,
        seed_x, seed_x, value1, deriv1, solution1.coefficients,
        patches1, match_value1, match_deriv1, length(patches1),
        max_step1, min_step1, 1.0 + 0.0im, truncations1)
    basis2 = Basis(
        kind2, scale2, ComplexF64(solution2.exponent), 0.0 + 0.0im,
        seed_x, seed_x, value2, deriv2, solution2.coefficients,
        patches2, match_value2, match_deriv2, length(patches2),
        max_step2, min_step2, 1.0 + 0.0im, truncations2)
    return basis1, basis2
end

function _hybrid_horizon(
    coefficients::CoefficientSet,
    solution,
    branch::Symbol,
    seed_x::Float64,
    match_x::Float64,
    work::IterationScratch,
    controls,
)
    kind = _kind(:horizon, branch)
    raw, draw = _endpoint_pair(
        kind, solution.coefficients, solution.exponent, 0.0 + 0.0im, seed_x)
    scale = endpoint_seed_scale(kind, coefficients.params)
    seed_value = scale * raw
    seed_derivative = scale * draw
    result = _hybrid_state(
        coefficients, kind, seed_x, seed_value, seed_derivative, match_x;
        controls, scratch=work)
    truncations = Truncation[]
    _record_truncation!(truncations, :horizon, kind,
        solution.requested_order, solution.effective_order,
        solution.truncation_reason)
    append!(truncations, result.truncations)
    return Basis(
        kind,
        scale,
        ComplexF64(solution.exponent),
        0.0 + 0.0im,
        seed_x,
        seed_x,
        seed_value,
        seed_derivative,
        solution.coefficients,
        result.patches,
        result.match_X,
        result.match_dXdx,
        result.step_count,
        result.max_step,
        result.min_step,
        1.0 + 0.0im,
        truncations,
    )
end

# Eikonal route from the horizon seed: `eikonal_path` when the slow-variation
# score is at most _EIKONAL_SLOW_VARIATION, `eikonal_split` when only
# _EIKONAL_CERTIFIED_VARIATION is met and a residual `limit` is given; nothing
# otherwise. `_eikonal_horizon_pair` is the two-branch version.
function _eikonal_horizon(
    coefficients::CoefficientSet,
    solution,
    branch::Symbol,
    seed_x::Float64,
    match_x::Float64,
    work::IterationScratch,
    limit,
)
    score = _eikonal_score(coefficients, seed_x, work)
    fast = score <= _EIKONAL_SLOW_VARIATION
    certified = limit !== nothing &&
        score <= _EIKONAL_CERTIFIED_VARIATION
    (fast || certified) || return nothing
    kind = _kind(:horizon, branch)
    raw, draw = _endpoint_pair(
        kind, solution.coefficients, solution.exponent, 0.0 + 0.0im, seed_x)
    scale = endpoint_seed_scale(kind, coefficients.params)
    seed_value = scale * raw
    seed_derivative = scale * draw
    result = limit === nothing ?
        eikonal_path(
            coefficients, seed_x, seed_value, seed_derivative, match_x;
            avec=work.avec, bvec=work.bvec, pq=work.pq,
            scratch=work.eikonal, large_frequency_expansion=work.large_frequency_expansion) :
        fast ? eikonal_path(
            coefficients, seed_x, seed_value, seed_derivative, match_x;
            avec=work.avec, bvec=work.bvec, pq=work.pq,
            scratch=work.eikonal, residual_limit=limit,
            check_limit=limit, large_frequency_expansion=work.large_frequency_expansion) :
        eikonal_split(
            coefficients, seed_x, seed_value, seed_derivative, match_x;
            avec=work.avec, bvec=work.bvec, pq=work.pq,
            scratch=work.eikonal, residual_limit=limit,
            check_limit=limit, large_frequency_expansion=work.large_frequency_expansion)
    result === nothing && return nothing
    truncations = Truncation[]
    _record_truncation!(truncations, :horizon, kind,
        solution.requested_order, solution.effective_order, solution.truncation_reason)
    return Basis(
        kind,
        scale,
        ComplexF64(solution.exponent),
        0.0 + 0.0im,
        seed_x,
        seed_x,
        seed_value,
        seed_derivative,
        solution.coefficients,
        result.patches,
        result.value,
        result.derivative,
        length(result.patches),
        result.max_step,
        result.min_step,
        1.0 + 0.0im,
        truncations,
    )
end

function _hybrid_horizon_pair(
    coefficients::CoefficientSet,
    solution1,
    solution2,
    branch1::Symbol,
    branch2::Symbol,
    seed_x::Float64,
    match_x::Float64,
    work::IterationScratch,
    controls,
)
    kind1 = _kind(:horizon, branch1)
    kind2 = _kind(:horizon, branch2)
    raw1, draw1 = _endpoint_pair(
        kind1, solution1.coefficients, solution1.exponent,
        0.0 + 0.0im, seed_x)
    raw2, draw2 = _endpoint_pair(
        kind2, solution2.coefficients, solution2.exponent,
        0.0 + 0.0im, seed_x)
    scale1 = endpoint_seed_scale(kind1, coefficients.params)
    scale2 = endpoint_seed_scale(kind2, coefficients.params)
    seed_value1 = scale1 * raw1
    seed_derivative1 = scale1 * draw1
    seed_value2 = scale2 * raw2
    seed_derivative2 = scale2 * draw2
    result1, result2 = _hybrid_pair(
        coefficients, kind1, kind2, seed_x,
        seed_value1, seed_derivative1, seed_value2, seed_derivative2,
        match_x; controls, scratch=work)
    truncations1 = Truncation[]
    truncations2 = Truncation[]
    _record_truncation!(truncations1, :horizon, kind1,
        solution1.requested_order, solution1.effective_order,
        solution1.truncation_reason)
    _record_truncation!(truncations2, :horizon, kind2,
        solution2.requested_order, solution2.effective_order,
        solution2.truncation_reason)
    append!(truncations1, result1.truncations)
    append!(truncations2, result2.truncations)
    basis1 = Basis(
        kind1,
        scale1,
        ComplexF64(solution1.exponent),
        0.0 + 0.0im,
        seed_x,
        seed_x,
        seed_value1,
        seed_derivative1,
        solution1.coefficients,
        result1.patches,
        result1.match_X,
        result1.match_dXdx,
        result1.step_count,
        result1.max_step,
        result1.min_step,
        1.0 + 0.0im,
        truncations1,
    )
    basis2 = Basis(
        kind2,
        scale2,
        ComplexF64(solution2.exponent),
        0.0 + 0.0im,
        seed_x,
        seed_x,
        seed_value2,
        seed_derivative2,
        solution2.coefficients,
        result2.patches,
        result2.match_X,
        result2.match_dXdx,
        result2.step_count,
        result2.max_step,
        result2.min_step,
        1.0 + 0.0im,
        truncations2,
    )
    return basis1, basis2
end

function _eikonal_horizon_pair(
    coefficients::CoefficientSet,
    solution1,
    solution2,
    branch1::Symbol,
    branch2::Symbol,
    seed_x::Float64,
    match_x::Float64,
    work::IterationScratch,
    limit,
)
    score = _eikonal_score(coefficients, seed_x, work)
    fast = score <= _EIKONAL_SLOW_VARIATION
    certified = limit !== nothing &&
        score <= _EIKONAL_CERTIFIED_VARIATION
    (fast || certified) || return nothing
    kind1 = _kind(:horizon, branch1)
    kind2 = _kind(:horizon, branch2)
    raw1, draw1 = _endpoint_pair(
        kind1, solution1.coefficients, solution1.exponent, 0.0 + 0.0im, seed_x)
    raw2, draw2 = _endpoint_pair(
        kind2, solution2.coefficients, solution2.exponent, 0.0 + 0.0im, seed_x)
    scale1 = endpoint_seed_scale(kind1, coefficients.params)
    scale2 = endpoint_seed_scale(kind2, coefficients.params)
    seed_value1 = scale1 * raw1
    seed_derivative1 = scale1 * draw1
    seed_value2 = scale2 * raw2
    seed_derivative2 = scale2 * draw2
    result = limit === nothing ?
        eikonal_pair(
            coefficients, seed_x, seed_value1, seed_derivative1,
            seed_value2, seed_derivative2, match_x;
            avec=work.avec, bvec=work.bvec, pq=work.pq,
            scratch=work.eikonal, large_frequency_expansion=work.large_frequency_expansion) :
        fast ? eikonal_pair(
            coefficients, seed_x, seed_value1, seed_derivative1,
            seed_value2, seed_derivative2, match_x;
            avec=work.avec, bvec=work.bvec, pq=work.pq,
            scratch=work.eikonal, residual_limit=limit,
            check_limit=limit, large_frequency_expansion=work.large_frequency_expansion) :
        eikonal_split_pair(
            coefficients, seed_x, seed_value1, seed_derivative1,
            seed_value2, seed_derivative2, match_x;
            avec=work.avec, bvec=work.bvec, pq=work.pq,
            scratch=work.eikonal, residual_limit=limit,
            check_limit=limit, large_frequency_expansion=work.large_frequency_expansion)
    result === nothing && return nothing
    truncations1 = Truncation[]
    truncations2 = Truncation[]
    _record_truncation!(truncations1, :horizon, kind1,
        solution1.requested_order, solution1.effective_order, solution1.truncation_reason)
    _record_truncation!(truncations2, :horizon, kind2,
        solution2.requested_order, solution2.effective_order, solution2.truncation_reason)
    basis1 = Basis(
        kind1,
        scale1,
        ComplexF64(solution1.exponent),
        0.0 + 0.0im,
        seed_x,
        seed_x,
        seed_value1,
        seed_derivative1,
        solution1.coefficients,
        result.patches1,
        result.value1,
        result.deriv1,
        length(result.patches1),
        result.max_step,
        result.min_step,
        1.0 + 0.0im,
        truncations1,
    )
    basis2 = Basis(
        kind2,
        scale2,
        ComplexF64(solution2.exponent),
        0.0 + 0.0im,
        seed_x,
        seed_x,
        seed_value2,
        seed_derivative2,
        solution2.coefficients,
        result.patches2,
        result.value2,
        result.deriv2,
        length(result.patches2),
        result.max_step,
        result.min_step,
        1.0 + 0.0im,
        truncations2,
    )
    return basis1, basis2
end

# Seed point of a horizon solution from `select_horizon_seed` (using the
# logarithmic part's series for a threshold solution).
@inline function _horizon_seed_value(solution, coefficients, controls, match_x)
    seed_coefficients = hasproperty(solution, :log_coefficients) &&
        !isempty(solution.log_coefficients) ?
        solution.log_coefficients : solution.coefficients
    seed, _ = select_horizon_seed(
        seed_coefficients,
        solution.exponent,
        solution.A,
        solution.B,
        match_x;
        tol=max(getproperty(controls, :tolerance), _ENDPOINT_SEED_TARGET),
        tail_aware=_tail_seed(coefficients),
    )
    return seed
end

@inline function _horizon_reach_cost(seed, match_x)
    ratio = match_x / max(seed, floatmin(Float64))
    return log(max(ratio, 1.0))
end

function _horizon_reach_promising(
        solution, coefficients, seed, match_x)
    current_cost = _horizon_reach_cost(seed, match_x)
    current_cost >= _HORIZON_REACH_MIN_SAVING || return false, current_cost

    tail = _tail_seed(coefficients) ?
        _horizon_tail_metric(solution.coefficients, seed) :
        abs(solution.coefficients[end] *
            seed^(length(solution.coefficients) - 1)) /
        max(abs(poly_value(solution.coefficients, seed)), eps(Float64))
    residual = _horizon_residual_metric(
        solution.coefficients, solution.exponent,
        solution.A, solution.B, seed)
    roundoff = 8eps(Float64) * _horizon_eval_condition(
        solution.coefficients, solution.exponent, seed)
    return roundoff < max(tail, residual), current_cost
end

# Whether the order-_HORIZON_REACH_ORDER seed extension is tried: switched on by
# ISEM_HORIZON_REACH_AUTO, off in the small- and large-frequency regimes and
# when the horizon order is already at least _HORIZON_REACH_ORDER.
@inline function _horizon_reach_auto(coefficients, controls)
    _HORIZON_REACH_AUTO || return false
    getproperty(controls, :small_frequency_expansion) && return false
    getproperty(controls, :large_frequency_expansion) && return false
    getproperty(controls, :horizon_order) < _HORIZON_REACH_ORDER || return false
    return true
end

# Horizon series and seed point. When the seed is limited by truncation rather
# than roundoff, an order-_HORIZON_REACH_ORDER series is tried and kept if it
# moves the seed out enough to save propagation (cost = log(match_x / seed)).
function _horizon_seed(coefficients, branch, controls, match_x)
    order = getproperty(controls, :horizon_order)
    solution = zero_local_solution(coefficients, branch, order)
    seed = _horizon_seed_value(
        solution, coefficients, controls, match_x)
    _horizon_reach_auto(coefficients, controls) || return solution, seed
    promising, current_cost = _horizon_reach_promising(
        solution, coefficients, seed, match_x)
    promising || return solution, seed

    candidate = zero_local_solution(
        coefficients, branch, _HORIZON_REACH_ORDER)
    candidate_seed = _horizon_seed_value(
        candidate, coefficients, controls, match_x)
    candidate_cost = _horizon_reach_cost(candidate_seed, match_x)
    (candidate_seed > seed &&
        candidate_cost <= _HORIZON_REACH_COST_RATIO * current_cost &&
        current_cost - candidate_cost >= _HORIZON_REACH_MIN_SAVING) ||
        return solution, seed
    return candidate, candidate_seed
end

function _horizon_pair(coefficients, branch1, branch2, controls, match_x)
    order = getproperty(controls, :horizon_order)
    threshold = is_horizon_threshold(coefficients.params)
    solution1 = threshold && branch1 == :out ?
        zero_log_solution(coefficients, :out, order) :
        zero_local_solution(coefficients, branch1, order)
    solution2 = threshold && branch2 == :out ?
        zero_log_solution(coefficients, :out, order) :
        zero_local_solution(coefficients, branch2, order)
    seed1 = _horizon_seed_value(solution1, coefficients, controls, match_x)
    seed2 = _horizon_seed_value(solution2, coefficients, controls, match_x)
    seed = min(seed1, seed2)
    threshold && return solution1, solution2, seed
    _horizon_reach_auto(coefficients, controls) ||
        return solution1, solution2, seed
    limiting_solution = seed1 <= seed2 ? solution1 : solution2
    promising, current_cost = _horizon_reach_promising(
        limiting_solution, coefficients, seed, match_x)
    promising || return solution1, solution2, seed

    candidate1 = zero_local_solution(
        coefficients, branch1, _HORIZON_REACH_ORDER)
    candidate2 = zero_local_solution(
        coefficients, branch2, _HORIZON_REACH_ORDER)
    candidate_seed1 = _horizon_seed_value(
        candidate1, coefficients, controls, match_x)
    candidate_seed2 = _horizon_seed_value(
        candidate2, coefficients, controls, match_x)
    candidate_seed = min(candidate_seed1, candidate_seed2)
    candidate_cost = _horizon_reach_cost(candidate_seed, match_x)
    (candidate_seed > seed &&
        candidate_cost <= _HORIZON_REACH_COST_RATIO * current_cost &&
        current_cost - candidate_cost >= _HORIZON_REACH_MIN_SAVING) ||
        return solution1, solution2, seed
    return candidate1, candidate2, candidate_seed
end

# Horizon basis of one branch (:in or :out) from the Frobenius series at the
# horizon, propagated outward to match_x: eikonal or hybrid paths when enabled,
# otherwise Taylor patches.
function iterate_from_zero(
    coefficients::CoefficientSet,
    branch::Symbol;
    controls,
    match_x::Float64=getproperty(controls, :match_x),
    scratch=nothing,
)
    normalized = _normalized_horizon_branch(branch)
    endpoint_solution, seed_x = _horizon_seed(
        coefficients, normalized, controls, match_x)
    if _eikonal_enabled(controls, coefficients, scratch)
        work = scratch === nothing ? _scratch(coefficients,
            getproperty(controls, :ordinary_order), controls) : scratch
        @static if _HORIZON_LARGE_FREQUENCY_EXPANSION
            if getproperty(controls, :large_frequency_expansion)
                basis = _large_frequency_expansion_horizon_basis(
                    coefficients, endpoint_solution, normalized,
                    seed_x, match_x, work, controls)
                basis === nothing || return basis
            end
        end
        limit = (getproperty(controls, :large_frequency_expansion) || _ORDINARY_EIKONAL) ?
            _INFINITY_EIKONAL_LIMIT : nothing
        basis = _eikonal_horizon(
            coefficients, endpoint_solution, normalized,
            seed_x, match_x, work, limit)
        if basis === nothing && getproperty(controls, :large_frequency_expansion)
            basis = _hybrid_horizon(
                coefficients, endpoint_solution, normalized,
                seed_x, match_x, work, controls)
        end
        basis === nothing || return basis
    end
    return _build_basis_series(
        coefficients,
        _kind(:horizon, normalized),
        endpoint_solution.coefficients,
        endpoint_solution.exponent,
        0.0 + 0.0im,
        match_x,
        seed_x,
        getproperty(controls, :ordinary_order);
        controls=controls,
        scratch=scratch,
        endpoint_requested_order=endpoint_solution.requested_order,
        endpoint_effective_order=endpoint_solution.effective_order,
        endpoint_reason=endpoint_solution.truncation_reason,
    )
end

# Both horizon branches from one common seed with shared steps, returning
# (basis1, basis2). At the superradiant threshold the outgoing branch uses the
# logarithmic solution and the eikonal routes are skipped.
function iterate_pair_from_zero(
    coefficients::CoefficientSet,
    branch1::Symbol,
    branch2::Symbol;
    controls,
    match_x::Float64=getproperty(controls, :match_x),
    scratch=nothing,
)
    normalized1 = _normalized_horizon_branch(branch1)
    normalized2 = _normalized_horizon_branch(branch2)
    solution1, solution2, seed_x = _horizon_pair(
        coefficients, normalized1, normalized2, controls, match_x)
    threshold_log_pair =
        hasproperty(solution1, :representation) &&
        solution1.representation == :threshold_log_frobenius ||
        hasproperty(solution2, :representation) &&
        solution2.representation == :threshold_log_frobenius
    if !threshold_log_pair &&
            _eikonal_enabled(controls, coefficients, scratch)
        work = scratch === nothing ? _scratch(coefficients,
            getproperty(controls, :ordinary_order), controls) : scratch
        @static if _HORIZON_LARGE_FREQUENCY_EXPANSION
            if getproperty(controls, :large_frequency_expansion)
                bases = _large_frequency_expansion_horizon_basis_pair(
                    coefficients, solution1, solution2,
                    normalized1, normalized2,
                    seed_x, match_x, work, controls)
                bases === nothing || return bases
            end
        end
        limit = (getproperty(controls, :large_frequency_expansion) || _ORDINARY_EIKONAL) ?
            _INFINITY_EIKONAL_LIMIT : nothing
        bases = _eikonal_horizon_pair(
            coefficients,
            solution1,
            solution2,
            normalized1,
            normalized2,
            seed_x,
            match_x,
            work,
            limit,
        )
        if bases === nothing && getproperty(controls, :large_frequency_expansion)
            bases = _hybrid_horizon_pair(
                coefficients, solution1, solution2, normalized1, normalized2,
                seed_x, match_x, work, controls)
        end
        bases === nothing || return bases
    end
    return _build_pair_from_horizon_seed(
        coefficients,
        solution1,
        solution2,
        normalized1,
        normalized2,
        seed_x,
        match_x;
        controls=controls,
        scratch=scratch,
    )
end

function _large_frequency_expansion_infinity_data(coefficients, selected, branch)
    kind = _kind(:infinity, branch)
    rho, sigma = infinity_exponents(coefficients.params, branch)
    endpoint_coeffs = branch == :out ?
        selected.out_coeffs : selected.in_coeffs
    endpoint_order = branch == :out ?
        selected.out_order : selected.in_order
    endpoint_reason = branch == :out ?
        selected.out_reason : selected.in_reason
    scale = endpoint_seed_scale(kind, coefficients.params)
    return (; kind, rho, sigma, endpoint_coeffs,
        endpoint_order, endpoint_reason, scale)
end

# Large-frequency infinity route: one infinity-series patch from endpoint_y
# straight to match_x; nothing when it is not certified or the selection uses a
# scaled variable. `_large_frequency_expansion_infinity_pair` is the two-branch version.
function _large_frequency_expansion_infinity_basis(coefficients, selected, branch, match_x)
    selected.variable_scale == 1 || return nothing
    data = _large_frequency_expansion_infinity_data(coefficients, selected, branch)
    candidate = large_frequency_expansion_infinity_single(coefficients, branch,
        selected.endpoint_y, match_x, data.scale)
    candidate === nothing && return nothing
    seed_x = 1 - selected.endpoint_y
    seed = large_frequency_expansion_infinity_state(candidate.patch, seed_x)
    truncations = Truncation[]
    _record_truncation!(truncations, :infinity, data.kind,
        selected.requested_order, data.endpoint_order,
        data.endpoint_reason)
    step = abs(seed_x - match_x)
    return Basis(
        data.kind, data.scale, ComplexF64(data.rho),
        ComplexF64(data.sigma), seed_x, selected.endpoint_y,
        seed.X, seed.dXdx, data.endpoint_coeffs,
        InfinityLargeFrequencyExpansionPatch[candidate.patch],
        candidate.state.X, candidate.state.dXdx,
        1, step, step, 1.0 + 0.0im, truncations)
end

function _large_frequency_expansion_infinity_pair(coefficients, selected, branch1, branch2, match_x)
    selected.variable_scale == 1 || return nothing
    data1 = _large_frequency_expansion_infinity_data(coefficients, selected, branch1)
    data2 = _large_frequency_expansion_infinity_data(coefficients, selected, branch2)
    candidate = large_frequency_expansion_infinity_pair(coefficients, branch1, branch2,
        selected.endpoint_y, match_x, data1.scale, data2.scale)
    candidate === nothing && return nothing
    seed_x = 1 - selected.endpoint_y
    seed1 = large_frequency_expansion_infinity_state(candidate.patch1, seed_x)
    seed2 = large_frequency_expansion_infinity_state(candidate.patch2, seed_x)
    function wrap(data, patch, seed, state)
        truncations = Truncation[]
        _record_truncation!(truncations, :infinity, data.kind,
            selected.requested_order, data.endpoint_order,
            data.endpoint_reason)
        step = abs(seed_x - match_x)
        return Basis(
            data.kind, data.scale, ComplexF64(data.rho),
            ComplexF64(data.sigma), seed_x, selected.endpoint_y,
            seed.X, seed.dXdx, data.endpoint_coeffs,
            InfinityLargeFrequencyExpansionPatch[patch], state.X, state.dXdx,
            1, step, step, 1.0 + 0.0im, truncations)
    end
    return wrap(data1, candidate.patch1, seed1, candidate.state1),
        wrap(data2, candidate.patch2, seed2, candidate.state2)
end

# Infinity basis of one branch (:out/:up or :in/:down) from the selected
# endpoint series, propagated inward to match_x (single large-frequency-expansion carrier patch when
# large_frequency_expansion is on and it reaches). certified_steps and step_budget are passed to
# `_build_basis_series`.
function iterate_from_infinity(
    coefficients::CoefficientSet,
    branch::Symbol;
    controls,
    match_x::Float64=getproperty(controls, :match_x),
    selection=nothing,
    scratch=nothing,
    certified_steps=false,
    step_budget=nothing,
)
    normalized = _normalized_infinity_branch(branch)
    selected = selection === nothing ?
        select_infinity_endpoint(coefficients, match_x; controls=controls) :
        selection
    @static if _INFINITY_LARGE_FREQUENCY_EXPANSION
        if getproperty(controls, :large_frequency_expansion)
            basis = _large_frequency_expansion_infinity_basis(
                coefficients, selected, normalized, match_x)
            basis === nothing || return basis
        end
    end
    rho, sigma = infinity_exponents(coefficients.params, normalized)
    endpoint_coeffs = normalized == :out ? selected.out_coeffs : selected.in_coeffs
    endpoint_order = normalized == :out ? selected.out_order : selected.in_order
    endpoint_reason = normalized == :out ? selected.out_reason : selected.in_reason
    variable_scale = selected.variable_scale
    endpoint_sigma = sigma / variable_scale
    return _build_basis_series(
        coefficients,
        _kind(:infinity, normalized),
        endpoint_coeffs,
        rho,
        endpoint_sigma,
        match_x,
        selected.endpoint_y,
        getproperty(controls, :ordinary_order);
        controls=controls,
        scratch=scratch,
        endpoint_variable_scale=variable_scale,
        endpoint_requested_order=selected.requested_order,
        endpoint_effective_order=endpoint_order,
        endpoint_reason=endpoint_reason,
        certified_steps,
        step_budget,
    )
end

# Both infinity branches from the common selected endpoint, marched together by
# `_hybrid_pair` (or the large-frequency patch pair when enabled); returns
# (basis1, basis2).
function iterate_pair_from_infinity(
    coefficients::CoefficientSet,
    branch1::Symbol,
    branch2::Symbol;
    controls,
    match_x::Float64=getproperty(controls, :match_x),
    selection=nothing,
    scratch=nothing,
)
    normalized1 = _normalized_infinity_branch(branch1)
    normalized2 = _normalized_infinity_branch(branch2)
    selected = selection === nothing ?
        select_infinity_endpoint(coefficients, match_x;
            controls=controls) : selection
    @static if _INFINITY_LARGE_FREQUENCY_EXPANSION
        if getproperty(controls, :large_frequency_expansion)
            bases = _large_frequency_expansion_infinity_pair(
                coefficients, selected, normalized1, normalized2, match_x)
            bases === nothing || return bases
        end
    end
    variable_scale = selected.variable_scale
    seed_coord = selected.endpoint_y
    seed_x = 1.0 - seed_coord

    function endpoint_data(branch)
        kind = _kind(:infinity, branch)
        rho, sigma = infinity_exponents(coefficients.params, branch)
        endpoint_coeffs = branch == :out ?
            selected.out_coeffs : selected.in_coeffs
        endpoint_order = branch == :out ?
            selected.out_order : selected.in_order
        endpoint_reason = branch == :out ?
            selected.out_reason : selected.in_reason
        endpoint_sigma = sigma / variable_scale
        raw, derivative = _endpoint_pair(
            kind, endpoint_coeffs, rho, endpoint_sigma,
            seed_coord, variable_scale)
        scale = endpoint_seed_scale(kind, coefficients.params)
        return (; kind, rho, endpoint_sigma, endpoint_coeffs,
            endpoint_order, endpoint_reason, scale,
            value=scale * raw, derivative=scale * derivative)
    end

    data1 = endpoint_data(normalized1)
    data2 = endpoint_data(normalized2)
    basis1, basis2 = _hybrid_pair(
        coefficients,
        data1.kind,
        data2.kind,
        seed_x,
        data1.value,
        data1.derivative,
        data2.value,
        data2.derivative,
        match_x;
        controls,
        scratch,
    )

    function wrap(data, basis)
        truncations = Truncation[]
        _record_truncation!(truncations, :infinity, data.kind,
            selected.requested_order, data.endpoint_order,
            data.endpoint_reason)
        append!(truncations, basis.truncations)
        return Basis(
            data.kind,
            data.scale,
            ComplexF64(data.rho),
            ComplexF64(data.endpoint_sigma),
            seed_x,
            seed_coord,
            data.value,
            data.derivative,
            data.endpoint_coeffs,
            basis.patches,
            basis.match_X,
            basis.match_dXdx,
            basis.step_count,
            basis.max_step,
            basis.min_step,
            ComplexF64(variable_scale),
            truncations,
        )
    end
    return wrap(data1, basis1), wrap(data2, basis2)
end

function _patch_contains(patch::Patch, x)
    lo = min(patch.center_x, patch.next_x)
    hi = max(patch.center_x, patch.next_x)
    return lo - 20eps(Float64) <= x <= hi + 20eps(Float64)
end

function _patch_contains(patch::EikonalPatch, x)
    lo = min(patch.center_x, patch.next_x)
    hi = max(patch.center_x, patch.next_x)
    return lo - 20eps(Float64) <= x <= hi + 20eps(Float64)
end

function _patch_contains(patch::HorizonLargeFrequencyExpansionPatch, x)
    lo = min(patch.center_x, patch.next_x)
    hi = max(patch.center_x, patch.next_x)
    return lo - 20eps(Float64) <= x <= hi + 20eps(Float64)
end

function _patch_contains(patch::InfinityLargeFrequencyExpansionPatch, x)
    lo = min(patch.center_x, patch.next_x)
    hi = max(patch.center_x, patch.next_x)
    return lo - 20eps(Float64) <= x <= hi + 20eps(Float64)
end

function _patch_contains(patch::ScaledPatch, x)
    lo = min(patch.center_x, patch.next_x)
    hi = max(patch.center_x, patch.next_x)
    return lo - 20eps(Float64) <= x <= hi + 20eps(Float64)
end

function _patch_contains(patch::LogScaledPatch, x)
    lo = min(patch.center_x, patch.next_x)
    hi = max(patch.center_x, patch.next_x)
    return lo - 20eps(Float64) <= x <= hi + 20eps(Float64)
end

# Binary search for the first patch whose interval (with a 20 eps margin)
# contains x, for patches ordered in either direction; 0 if none does.
function _patch_index(patches, x)
    isempty(patches) && return 0
    forward = patches[1].center_x <= patches[end].center_x
    lo = 1
    hi = length(patches)
    while lo <= hi
        mid = (lo + hi) >>> 1
        patch = patches[mid]
        xmin = min(patch.center_x, patch.next_x) - 20eps(Float64)
        xmax = max(patch.center_x, patch.next_x) + 20eps(Float64)
        if x < xmin
            forward ? (hi = mid - 1) : (lo = mid + 1)
        elseif x > xmax
            forward ? (lo = mid + 1) : (hi = mid - 1)
        else
            while mid > 1 && _patch_contains(patches[mid - 1], x)
                mid -= 1
            end
            return mid
        end
    end
    return 0
end

_patch_value(patch::Patch, x) = patch.scale *
    poly_value(patch.coeffs, x - patch.center_x)
_patch_value(patch::EikonalPatch, x) = eikonal_value(patch, x)
_patch_value(patch::HorizonLargeFrequencyExpansionPatch, x) = large_frequency_expansion_horizon_value(patch, x)
_patch_value(patch::InfinityLargeFrequencyExpansionPatch, x) = large_frequency_expansion_infinity_value(patch, x)
function _patch_value(patch::ScaledPatch, x)
    y = (1.0 - Float64(x)) / patch.omega_scale
    return patch.scale * poly_value(patch.coeffs, y - patch.center_y)
end

function _patch_state(patch::Patch, x)
    value, derivative = poly_pair(patch.coeffs, x - patch.center_x)
    return patch.scale * value, patch.scale * derivative
end

_patch_state(patch::EikonalPatch, x) = eikonal_state(patch, x)
function _patch_state(patch::HorizonLargeFrequencyExpansionPatch, x)
    state = large_frequency_expansion_horizon_state(patch, x)
    return state.X, state.dXdx
end
function _patch_state(patch::InfinityLargeFrequencyExpansionPatch, x)
    state = large_frequency_expansion_infinity_state(patch, x)
    return state.X, state.dXdx
end
function _patch_state(patch::ScaledPatch, x)
    y = (1.0 - Float64(x)) / patch.omega_scale
    value, derivative_y = poly_pair(
        patch.coeffs, y - patch.center_y)
    return patch.scale * value,
        -patch.scale * derivative_y / patch.omega_scale
end

# Value of a basis at x: the endpoint series between the endpoint and the seed
# (when valid), otherwise the patch covering x; errors outside the basis.
function basis_value(basis::Basis, x)
    x == basis.seed_x && return basis.seed_X
    if basis.endpoint_valid && _is_horizon(basis.kind)
        if x <= basis.seed_x
            value = _endpoint_value(basis.kind, basis.endpoint_coeffs, basis.exponent1, basis.exponent2, x)
            isempty(basis.endpoint_log_coeffs) || (value += _horizon_log_pair(basis, x)[1])
            return basis.scale * value
        end
    elseif basis.endpoint_valid
        x >= basis.seed_x && return basis.scale *
            _endpoint_value(
                basis.kind,
                basis.endpoint_coeffs,
                basis.exponent1,
                basis.exponent2,
                x,
                basis.endpoint_variable_scale,
            )
    end
    patch_index = _patch_index(basis.patches, x)
    patch_index > 0 && return _patch_value(basis.patches[patch_index], x)
    isempty(basis.patches) && error("ISEM GSN basis has no patch for x=$x")
    error("ISEM GSN basis does not cover x=$x")
end

# (X, dX/dx) of a basis at x: the endpoint series between the endpoint and the
# seed (when valid), otherwise the patch covering x; errors outside the basis.
function basis_state(basis::Basis, x)
    x == basis.seed_x && return basis.seed_X, basis.seed_dXdx
    if basis.endpoint_valid && _is_horizon(basis.kind)
        if x <= basis.seed_x
            value, derivative = _endpoint_state(basis.kind, basis.endpoint_coeffs,
                basis.exponent1, basis.exponent2, x)
            if !isempty(basis.endpoint_log_coeffs)
                log_value, log_derivative = _horizon_log_pair(basis, x)
                value += log_value
                derivative += log_derivative
            end
            return basis.scale * value, basis.scale * derivative
        end
    elseif basis.endpoint_valid
        if x >= basis.seed_x
            value, derivative = _endpoint_state(basis.kind, basis.endpoint_coeffs,
                basis.exponent1, basis.exponent2, x, basis.endpoint_variable_scale)
            return basis.scale * value, basis.scale * derivative
        end
    end
    patch_index = _patch_index(basis.patches, x)
    patch_index > 0 && return _patch_state(basis.patches[patch_index], x)
    isempty(basis.patches) && error("ISEM GSN basis has no patch for x=$x")
    error("ISEM GSN basis does not cover x=$x")
end

# State at x when the caller also supplies infinity_y = 1 - x computed from the
# physical radius; near infinity this keeps the digits that forming 1 - x loses.
function basis_state(basis::Basis{Patch}, x, infinity_y)
    if basis.endpoint_valid && !_is_horizon(basis.kind) && x >= basis.seed_x
        value, derivative, _ = _scaled_infinity_series_triple(
            basis.endpoint_coeffs, basis.exponent1, basis.exponent2,
            infinity_y, basis.endpoint_variable_scale)
        return basis.scale * value, basis.scale * derivative
    end
    index = _patch_index(basis.patches, x)
    index > 0 || return basis_state(basis, x)
    patch = basis.patches[index]
    # Offset from the patch centre, x - center_x = (1 - center_x) - infinity_y.
    value, derivative = poly_pair(
        patch.coeffs, (1 - patch.center_x) - infinity_y)
    return patch.scale * value, patch.scale * derivative
end

# LogScaledState of a log-scaled basis at x: the seed state at seed_x, otherwise
# the covering patch renormalized; errors outside the basis.
function logscaled_basis_state(basis::LogScaledBasis, x)
    x == basis.seed_x && return basis.seed_state
    patch_index = _patch_index(basis.patches, x)
    patch_index > 0 || error("ISEM GSN log-scaled basis does not cover x=$x")
    patch = basis.patches[patch_index]
    value, derivative = poly_pair(
        patch.coeffs, Float64(x) - patch.center_x)
    return logscaled_state(value, derivative, patch.log_scale)
end

# Number of patches in a basis.
basis_patch_count(basis::Basis) = length(basis.patches)
basis_patch_count(basis::LogScaledBasis) = length(basis.patches)

end
