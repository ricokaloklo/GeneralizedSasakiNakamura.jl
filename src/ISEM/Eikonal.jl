module Eikonal

using ..CoefficientTables: CoefficientSet
using ..LargeFrequencyExpansion: LargeFrequencyExpansionPlan, large_frequency_expansion_ab!
using ..OrdinaryPointExpansion:
    ordinary_ab_series!,
    poly_pair,
    poly_triple,
    poly_value

export EikonalPatch, EikonalSelection, EikonalScratch
export eikonal_candidate, eikonal_preflight, eikonal_path, eikonal_pair
export eikonal_split, eikonal_split_pair
export eikonal_value, eikonal_state
export EIKONAL_TERMS, EIKONAL_ORDER, LARGE_FREQUENCY_EXPANSION_EIKONAL_ORDER
export reset_eikonal_profile!, eikonal_profile

const EIKONAL_TERMS = 8
const EIKONAL_ORDER = 24
const LARGE_FREQUENCY_EXPANSION_EIKONAL_ORDER = 24
const EIKONAL_CHECK_GAP = 4
const EIKONAL_RESIDUAL_LIMIT = 2e-11 + 256eps(Float64)
const EIKONAL_CHECK_LIMIT = 2e-11 + 256eps(Float64)
const EIKONAL_SPLIT_LIMIT = 1e-10 + 256eps(Float64)
const EIKONAL_MAX_STEP = 0.1
const EIKONAL_RADIUS_SAFETY = 0.8
const EIKONAL_STEP_GROWTH = 1.55
const EIKONAL_PARTIAL_STALL =
    get(ENV, "ISEM_EIKONAL_PARTIAL_STALL", "0") != "0"
const EIKONAL_PARTIAL_STALL_RATIO = sqrt(eps(Float64))
const EIKONAL_COMPRESSION_MIN = 1536.0
const EIKONAL_SMALL_FREQUENCY_EXPANSION_COMPRESSION_MIN = 8192.0
const EIKONAL_ROUTE_COMPRESSION_MIN = 192.0
const EIKONAL_IN_PHASE_MIN = 5.0
# Optional counters of accepted and rejected eikonal steps, compiled in only when
# ISEM_EIKONAL_PROFILE=1 is set at load time (read with `eikonal_profile`).
const EIKONAL_PROFILE = get(ENV, "ISEM_EIKONAL_PROFILE", "0") == "1"
const _PROFILE_ACCEPTED = Threads.Atomic{Int}(0)
const _PROFILE_REJECTED = Threads.Atomic{Int}(0)
const _PROFILE_RESIDUAL_REJECTED = Threads.Atomic{Int}(0)
const _PROFILE_CHECK_REJECTED = Threads.Atomic{Int}(0)
const _PROFILE_SPLIT_REJECTED = Threads.Atomic{Int}(0)
const _PROFILE_NONFINITE_REJECTED = Threads.Atomic{Int}(0)
const _PROFILE_OTHER_REJECTED = Threads.Atomic{Int}(0)
const _PROFILE_PLANS = Threads.Atomic{Int}(0)
const _PROFILE_PLAN_NS = Threads.Atomic{UInt64}(0)

function reset_eikonal_profile!()
    _PROFILE_ACCEPTED[] = 0
    _PROFILE_REJECTED[] = 0
    _PROFILE_RESIDUAL_REJECTED[] = 0
    _PROFILE_CHECK_REJECTED[] = 0
    _PROFILE_SPLIT_REJECTED[] = 0
    _PROFILE_NONFINITE_REJECTED[] = 0
    _PROFILE_OTHER_REJECTED[] = 0
    _PROFILE_PLANS[] = 0
    _PROFILE_PLAN_NS[] = 0
    return nothing
end

function eikonal_profile()
    accepted = _PROFILE_ACCEPTED[]
    plans = _PROFILE_PLANS[]
    return (
        enabled=EIKONAL_PROFILE,
        accepted_patch_count=accepted,
        rejected_attempt_count=_PROFILE_REJECTED[],
        rejected_residual=_PROFILE_RESIDUAL_REJECTED[],
        rejected_check=_PROFILE_CHECK_REJECTED[],
        rejected_split=_PROFILE_SPLIT_REJECTED[],
        rejected_nonfinite=_PROFILE_NONFINITE_REJECTED[],
        rejected_other=_PROFILE_OTHER_REJECTED[],
        plan_build_count=plans,
        plans_per_accepted_patch=plans / max(accepted, 1),
        plan_build_us=Float64(_PROFILE_PLAN_NS[]) / 1.0e3,
    )
end

@inline function _profile_accept!(count)
    @static if EIKONAL_PROFILE
        Threads.atomic_add!(_PROFILE_ACCEPTED, count)
    end
    return nothing
end

@inline function _profile_reject!(reason::Symbol=:other)
    @static if EIKONAL_PROFILE
        Threads.atomic_add!(_PROFILE_REJECTED, 1)
        reason === :residual ? (Threads.atomic_add!(_PROFILE_RESIDUAL_REJECTED, 1)) :
        reason === :check ? (Threads.atomic_add!(_PROFILE_CHECK_REJECTED, 1)) :
        reason === :split ? (Threads.atomic_add!(_PROFILE_SPLIT_REJECTED, 1)) :
        reason === :nonfinite ? (Threads.atomic_add!(_PROFILE_NONFINITE_REJECTED, 1)) :
        (Threads.atomic_add!(_PROFILE_OTHER_REJECTED, 1))
    end
    return nothing
end

@inline function _reject_reason(residual, check, residual_limit, check_limit)
    (!isfinite(residual) || !isfinite(check)) && return :nonfinite
    residual > residual_limit && return :residual
    check > check_limit && return :check
    return :other
end

@inline function _reject_reason(
        residual, check, split, residual_limit, check_limit, split_limit)
    (!isfinite(residual) || !isfinite(check) || !isfinite(split)) &&
        return :nonfinite
    residual > residual_limit && return :residual
    check > check_limit && return :check
    split > split_limit && return :split
    return :other
end

# One eikonal patch on [center_x, next_x]. With h = x - center_x the solution is
# scale * exp(-aint/2) * Z * (cplus * exp(xi*eta) * plus + cminus * exp(-xi*eta) * minus),
# where A, B, eta, Z, aint, plus and minus are Taylor coefficients in h.
struct EikonalPatch
    center_x::Float64
    next_x::Float64
    xi::ComplexF64
    A::Vector{ComplexF64}
    B::Vector{ComplexF64}
    eta::Vector{ComplexF64}
    Z::Vector{ComplexF64}
    aint::Vector{ComplexF64}
    plus::Vector{ComplexF64}
    minus::Vector{ComplexF64}
    cplus::ComplexF64
    cminus::ComplexF64
    scale::Float64
end

struct EikonalPlan
    A::Vector{ComplexF64}
    B::Vector{ComplexF64}
    xi::ComplexF64
    eta::Vector{ComplexF64}
    Z::Vector{ComplexF64}
    aint::Vector{ComplexF64}
    plus::Vector{ComplexF64}
    minus::Vector{ComplexF64}
end

struct EikonalScratch
    f::Vector{ComplexF64}
    ft::Vector{ComplexF64}
    root::Vector{ComplexF64}
    invroot::Vector{ComplexF64}
    fp::Vector{ComplexF64}
    fpp::Vector{ComplexF64}
    work1::Vector{ComplexF64}
    work2::Vector{ComplexF64}
    curvroot::Vector{ComplexF64}
    previous::Vector{ComplexF64}
    next::Vector{ComplexF64}
    derivative::Vector{ComplexF64}
    product::Vector{ComplexF64}
    plan_A::Vector{ComplexF64}
    plan_B::Vector{ComplexF64}
    plan_eta::Vector{ComplexF64}
    plan_Z::Vector{ComplexF64}
    plan_aint::Vector{ComplexF64}
    plan_plus::Vector{ComplexF64}
    plan_minus::Vector{ComplexF64}
end

function EikonalScratch(order::Integer=EIKONAL_ORDER)
    order >= 2 || throw(ArgumentError("eikonal order must be at least two"))
    buffer() = zeros(ComplexF64, Int(order) + 1)
    return EikonalScratch(ntuple(_ -> buffer(), 20)...)
end

struct EikonalKernel
    A::Vector{ComplexF64}
    B::Vector{ComplexF64}
    xi::ComplexF64
    eta::Vector{ComplexF64}
    Z::Vector{ComplexF64}
    aint::Vector{ComplexF64}
    plus::Vector{ComplexF64}
    minus::Vector{ComplexF64}
    cplus::ComplexF64
    cminus::ComplexF64
end

struct EikonalSelection
    candidate::Bool
    reason::Symbol
    angular_compression::Float64
    phase_compression::Float64
    m_ratio::Float64
end

_finite(z) = isfinite(real(z)) && isfinite(imag(z))
_finite_vector(values) = all(_finite, values)

# Decide whether eikonal propagation is worth trying for these parameters, from
# the angular compression (l + 1/2)/kappa, the phase compression
# max(|omega|, |k_H|, |m Omega_H|)/kappa, the frequency regime and the route
# branch. Returns a EikonalSelection with the reason for the decision.
function eikonal_preflight(
    params,
    regime::Symbol=:ordinary,
    route_branch::Union{Nothing,Symbol}=nothing,
)
    l = getproperty(params, :l)
    m = getproperty(params, :m)
    a = getproperty(params, :a)
    omega = getproperty(params, :omega)
    kappa = max(getproperty(params, :kappa), eps(Float64))
    angular_compression = (l + 0.5) / kappa
    horizon_frequency = a / (2 * (1 + kappa))
    horizon_wavenumber = omega - m * horizon_frequency
    phase_scale = max(abs(omega), abs(horizon_wavenumber),
        abs(m * horizon_frequency))
    phase_compression = phase_scale / kappa
    m_ratio = abs(m) / max(l, 1)

    candidate, reason = if regime == :large_frequency_expansion
        true, :large_frequency
    elseif regime == :small_frequency_expansion
        if l >= 30
            true, :high_l_small_frequency
        elseif angular_compression >= EIKONAL_SMALL_FREQUENCY_EXPANSION_COMPRESSION_MIN
            true, :small_frequency_angular_compression
        elseif phase_compression >= EIKONAL_SMALL_FREQUENCY_EXPANSION_COMPRESSION_MIN
            true, :small_frequency_phase_compression
        else
            false, :small_frequency_taylor_cost
        end
    elseif angular_compression >= EIKONAL_COMPRESSION_MIN
        true, :angular_compression
    elseif phase_compression >= EIKONAL_COMPRESSION_MIN
        true, :phase_compression
    elseif l >= 30 && route_branch !== nothing &&
            max(angular_compression, phase_compression) >=
                EIKONAL_ROUTE_COMPRESSION_MIN
        true, :route_compression
    elseif l >= 30 && route_branch in (:UP, :up)
        false, :up_taylor_cost
    elseif l >= 30 && route_branch in (:IN, :in) &&
            phase_compression < EIKONAL_IN_PHASE_MIN
        false, :in_low_phase_taylor_cost
    elseif l >= 30
        true, :high_l
    else
        false, :taylor_cost_low
    end
    return EikonalSelection(candidate, reason,
        angular_compression, phase_compression, m_ratio)
end

function eikonal_candidate(
        coefficients::CoefficientSet,
        regime::Symbol=:ordinary,
        route_branch::Union{Nothing,Symbol}=nothing,
)
    return eikonal_preflight(
        coefficients.params, regime, route_branch).candidate
end

function _poly_copy(source, count)
    result = zeros(ComplexF64, count)
    count = min(length(source), length(result))
    @inbounds for index in 1:count
        result[index] = source[index]
    end
    return result
end

function _poly_copy!(result, source)
    fill!(result, 0)
    count = min(length(source), length(result))
    copyto!(result, 1, source, 1, count)
    return result
end

function _poly_mul!(out, left, right)
    fill!(out, 0)
    last = length(out)
    @inbounds for i in eachindex(left)
        li = left[i]
        iszero(li) && continue
        stop = min(length(right), last - i + 1)
        for j in 1:stop
            out[i + j - 1] += li * right[j]
        end
    end
    return out
end

function _poly_derivative!(out, source)
    fill!(out, 0)
    @inbounds for degree in 1:(min(length(source), length(out) + 1) - 1)
        out[degree] = degree * source[degree + 1]
    end
    return out
end

# Taylor series of g^(1/2), g^(-1/2), g^(-3) and g^(-1/4) from that of g, by
# the recurrence n g0 u_n = sum_k (p(n-k) - k) g_(n-k) u_k for u = g^p; false
# if g(0) is zero or a coefficient is not finite.
function _poly_wkb_powers!(root, invroot, inverse_cube, quarter,
        source)
    for output in (root, invroot, inverse_cube, quarter)
        fill!(output, 0)
    end
    g0 = source[1]
    _finite(g0) && !iszero(g0) || return false
    root[1] = g0^0.5
    invroot[1] = g0^-0.5
    inverse_cube[1] = g0^-3
    quarter[1] = g0^-0.25
    all(_finite, (root[1], invroot[1], inverse_cube[1], quarter[1])) ||
        return false
    @inbounds for n in 1:(length(root) - 1)
        rhs_root = ComplexF64(0)
        rhs_invroot = ComplexF64(0)
        rhs_cube = ComplexF64(0)
        rhs_quarter = ComplexF64(0)
        for k in 0:(n - 1)
            weighted_source = (n - k) * source[n - k + 1]
            rhs_root += root[k + 1] * weighted_source
            rhs_invroot += invroot[k + 1] * weighted_source
            rhs_cube += inverse_cube[k + 1] * weighted_source
            rhs_quarter += quarter[k + 1] * weighted_source
        end
        lhs_root = ComplexF64(0)
        lhs_invroot = ComplexF64(0)
        lhs_cube = ComplexF64(0)
        lhs_quarter = ComplexF64(0)
        for k in 1:(n - 1)
            weighted_source = source[k + 1] * (n - k)
            lhs_root += weighted_source * root[n - k + 1]
            lhs_invroot += weighted_source * invroot[n - k + 1]
            lhs_cube += weighted_source * inverse_cube[n - k + 1]
            lhs_quarter += weighted_source * quarter[n - k + 1]
        end
        denominator = n * g0
        root[n + 1] = (0.5 * rhs_root - lhs_root) / denominator
        invroot[n + 1] = (-0.5 * rhs_invroot - lhs_invroot) /
            denominator
        inverse_cube[n + 1] = (-3 * rhs_cube - lhs_cube) /
            denominator
        quarter[n + 1] = (-0.25 * rhs_quarter - lhs_quarter) /
            denominator
        all(_finite, (root[n + 1], invroot[n + 1],
            inverse_cube[n + 1], quarter[n + 1])) || return false
    end
    return true
end

function _poly_integral(source)
    result = zeros(ComplexF64, length(source))
    @inbounds for degree in 1:(length(source) - 1)
        result[degree + 1] = source[degree] / degree
    end
    return result
end

function _poly_integral!(result, source)
    fill!(result, 0)
    @inbounds for degree in 1:(min(length(source), length(result)) - 1)
        result[degree + 1] = source[degree] / degree
    end
    return result
end

# Local eikonal (WKB) plan of X'' + A X' + B X = 0 from the Taylor series of A
# and B at the current point. With f = A'/2 + A^2/4 - B and ft = f / f0,
# xi = sqrt(f0), eta = integral of sqrt(ft), Z = ft^(-1/4), and plus/minus are
# the amplitude corrections to `terms` orders in 1/(+xi) and 1/(-xi). Writes
# into `scratch` (patches must copy it); nothing if a series is not finite.
function _eikonal_plan!(
    a_coeffs,
    b_coeffs,
    scratch::EikonalScratch,
    terms::Integer=EIKONAL_TERMS,
)
    terms >= 0 || throw(ArgumentError("eikonal term count must be non-negative."))
    A = _poly_copy!(scratch.plan_A, a_coeffs)
    B = _poly_copy!(scratch.plan_B, b_coeffs)
    order = length(A) - 1
    f = scratch.f
    fill!(f, 0)
    @inbounds for n in 0:order
        derivative = n < order ? (n + 1) * A[n + 2] / 2 : 0
        square = ComplexF64(0)
        for k in 0:n
            square += A[k + 1] * A[n - k + 1]
        end
        f[n + 1] = derivative + square / 4 - B[n + 1]
    end
    f0 = f[1]
    _finite(f0) && !iszero(f0) || return nothing
    xi = ComplexF64(sqrt(f0))
    _finite(xi) && !iszero(xi) || return nothing
    ft = scratch.ft
    @inbounds for index in eachindex(ft)
        ft[index] = f[index] / f0
    end
    root = scratch.root
    invroot = scratch.invroot
    inverse_cube = scratch.product
    Z = scratch.plan_Z
    _poly_wkb_powers!(root, invroot, inverse_cube, Z, ft) ||
        return nothing
    fp = _poly_derivative!(scratch.fp, ft)
    fpp = _poly_derivative!(scratch.fpp, fp)
    work1 = _poly_mul!(scratch.work1, ft, fpp)
    work2 = _poly_mul!(scratch.work2, fp, fp)
    @inbounds for index in eachindex(work1)
        work1[index] = (4 * work1[index] - 5 * work2[index]) / 16
    end
    curvroot = _poly_mul!(scratch.curvroot, work1, inverse_cube)
    _poly_mul!(work1, curvroot, root)
    @inbounds for index in eachindex(curvroot)
        curvroot[index] = work1[index] / 2
    end

    plus = scratch.plan_plus
    minus = scratch.plan_minus
    fill!(plus, 0)
    fill!(minus, 0)
    previous = scratch.previous
    next = scratch.next
    fill!(previous, 0)
    fill!(next, 0)
    previous[1] = 1
    plus[1] = 1
    minus[1] = 1
    plus_weight = ComplexF64(1)
    minus_weight = ComplexF64(1)
    for _ in 1:terms
        derivative = _poly_derivative!(scratch.derivative, previous)
        product = _poly_mul!(scratch.product, derivative, invroot)
        @inbounds for index in eachindex(next)
            next[index] = -product[index] / 2
        end
        _poly_mul!(product, curvroot, previous)
        @inbounds for degree in 1:order
            next[degree + 1] += product[degree] / degree
        end
        plus_weight /= xi
        minus_weight /= -xi
        @inbounds for index in eachindex(next)
            plus[index] += plus_weight * next[index]
            minus[index] += minus_weight * next[index]
        end
        previous, next = next, previous
    end

    eta = _poly_integral!(scratch.plan_eta, root)
    aint = _poly_integral!(scratch.plan_aint, A)
    vectors = (A, B, eta, Z, aint, plus, minus)
    all(_finite_vector, vectors) || return nothing
    return EikonalPlan(A, B, xi, eta, Z, aint, plus, minus)
end

@inline function _profile_plan!(a_coeffs, b_coeffs, scratch, terms)
    @static if EIKONAL_PROFILE
        start = time_ns()
        plan = _eikonal_plan!(a_coeffs, b_coeffs, scratch, terms)
        Threads.atomic_add!(_PROFILE_PLAN_NS, time_ns() - start)
        Threads.atomic_add!(_PROFILE_PLANS, 1)
        return plan
    else
        return _eikonal_plan!(a_coeffs, b_coeffs, scratch, terms)
    end
end

# Allocating form of `_eikonal_plan!`: the plan vectors are fresh, so patches
# built with `_patch` (no copy) stay valid after the next plan.
function _eikonal_plan(
    a_coeffs,
    b_coeffs,
    scratch::EikonalScratch,
    terms::Integer=EIKONAL_TERMS,
)
    terms >= 0 || throw(ArgumentError("eikonal term count must be non-negative."))
    count = length(scratch.f)
    order = count - 1
    A = _poly_copy(a_coeffs, count)
    B = _poly_copy(b_coeffs, count)
    f = scratch.f
    fill!(f, 0)
    @inbounds for n in 0:order
        derivative = n < order ? (n + 1) * A[n + 2] / 2 : 0
        square = ComplexF64(0)
        for k in 0:n
            square += A[k + 1] * A[n - k + 1]
        end
        f[n + 1] = derivative + square / 4 - B[n + 1]
    end
    f0 = f[1]
    _finite(f0) && !iszero(f0) || return nothing
    xi = ComplexF64(sqrt(f0))
    _finite(xi) && !iszero(xi) || return nothing
    ft = scratch.ft
    @inbounds for index in eachindex(ft)
        ft[index] = f[index] / f0
    end
    root = scratch.root
    invroot = scratch.invroot
    inverse_cube = scratch.product
    Z = zeros(ComplexF64, count)
    _poly_wkb_powers!(root, invroot, inverse_cube, Z, ft) ||
        return nothing
    fp = _poly_derivative!(scratch.fp, ft)
    fpp = _poly_derivative!(scratch.fpp, fp)
    work1 = _poly_mul!(scratch.work1, ft, fpp)
    work2 = _poly_mul!(scratch.work2, fp, fp)
    @inbounds for index in eachindex(work1)
        work1[index] = (4 * work1[index] - 5 * work2[index]) / 16
    end
    curvroot = _poly_mul!(scratch.curvroot, work1, inverse_cube)
    _poly_mul!(work1, curvroot, root)
    @inbounds for index in eachindex(curvroot)
        curvroot[index] = work1[index] / 2
    end

    plus = zeros(ComplexF64, count)
    minus = zeros(ComplexF64, count)
    previous = scratch.previous
    next = scratch.next
    fill!(previous, 0)
    fill!(next, 0)
    previous[1] = 1
    plus[1] = 1
    minus[1] = 1
    plus_weight = ComplexF64(1)
    minus_weight = ComplexF64(1)
    for _ in 1:terms
        derivative = _poly_derivative!(scratch.derivative, previous)
        product = _poly_mul!(scratch.product, derivative, invroot)
        @inbounds for index in eachindex(next)
            next[index] = -product[index] / 2
        end
        _poly_mul!(product, curvroot, previous)
        @inbounds for degree in 1:order
            next[degree + 1] += product[degree] / degree
        end
        plus_weight /= xi
        minus_weight /= -xi
        @inbounds for index in eachindex(next)
            plus[index] += plus_weight * next[index]
            minus[index] += minus_weight * next[index]
        end
        previous, next = next, previous
    end

    eta = _poly_integral(root)
    aint = _poly_integral(A)
    vectors = (A, B, eta, Z, aint, plus, minus)
    all(_finite_vector, vectors) || return nothing
    return EikonalPlan(A, B, xi, eta, Z, aint, plus, minus)
end

@inline function _profile_plan(a_coeffs, b_coeffs, scratch, terms)
    @static if EIKONAL_PROFILE
        start = time_ns()
        plan = _eikonal_plan(a_coeffs, b_coeffs, scratch, terms)
        Threads.atomic_add!(_PROFILE_PLAN_NS, time_ns() - start)
        Threads.atomic_add!(_PROFILE_PLANS, 1)
        return plan
    else
        return _eikonal_plan(a_coeffs, b_coeffs, scratch, terms)
    end
end

# Fix the coefficients cplus/cminus of the two eikonal modes so that the
# kernel reproduces (X0, dX0) at h = 0; nothing if the modes are degenerate.
function _eikonal_kernel(plan::EikonalPlan, X0, dX0)
    eta0, deta0 = poly_pair(plan.eta, 0.0)
    z0, dz0 = poly_pair(plan.Z, 0.0)
    p0, dp0 = poly_pair(plan.plus, 0.0)
    m0, dm0 = poly_pair(plan.minus, 0.0)
    a0 = plan.A[1]
    vp = z0 * p0
    vm = z0 * m0
    dp = dz0 * p0 + z0 * (dp0 + plan.xi * deta0 * p0) - a0 * vp / 2
    dm = dz0 * m0 + z0 * (dm0 - plan.xi * deta0 * m0) - a0 * vm / 2
    denominator = vp * dm - vm * dp
    _finite(denominator) && !iszero(denominator) || return nothing
    cplus = (X0 * dm - dX0 * vm) / denominator
    cminus = (dX0 * vp - X0 * dp) / denominator
    _finite(cplus) && _finite(cminus) || return nothing
    return EikonalKernel(
        plan.A,
        plan.B,
        plan.xi,
        plan.eta,
        plan.Z,
        plan.aint,
        plan.plus,
        plan.minus,
        ComplexF64(cplus),
        ComplexF64(cminus),
    )
end

function _poly_pair(coeffs, h, degree)
    last = min(length(coeffs), degree + 1)
    return poly_pair(@view(coeffs[1:last]), h)
end

function _kernel_pair(kernel::EikonalKernel, h, degree=length(kernel.A) - 1)
    eta, deta = _poly_pair(kernel.eta, h, degree)
    Z, dZ = _poly_pair(kernel.Z, h, degree)
    plus, dplus = _poly_pair(kernel.plus, h, degree)
    minus, dminus = _poly_pair(kernel.minus, h, degree)
    plus *= kernel.cplus
    dplus *= kernel.cplus
    minus *= kernel.cminus
    dminus *= kernel.cminus
    A = poly_value(@view(kernel.A[1:min(end, degree + 1)]), h)
    aint = poly_value(@view(kernel.aint[1:min(end, degree + 1)]), h)
    factor = exp(-aint / 2)
    phase_plus = exp(kernel.xi * eta)
    phase_minus = exp(-kernel.xi * eta)
    chi = phase_plus * plus + phase_minus * minus
    dchi = phase_plus * (dplus + kernel.xi * deta * plus) +
        phase_minus * (dminus - kernel.xi * deta * minus)
    value = factor * Z * chi
    derivative = factor * (dZ * chi + Z * dchi - A * Z * chi / 2)
    return value, derivative
end

function _kernel_triple(kernel::EikonalKernel, h)
    eta, deta, d2eta = poly_triple(kernel.eta, h)
    Z, dZ, d2Z = poly_triple(kernel.Z, h)
    plus, dplus, d2plus = poly_triple(kernel.plus, h)
    minus, dminus, d2minus = poly_triple(kernel.minus, h)
    plus *= kernel.cplus
    dplus *= kernel.cplus
    d2plus *= kernel.cplus
    minus *= kernel.cminus
    dminus *= kernel.cminus
    d2minus *= kernel.cminus
    A, dA = poly_pair(kernel.A, h)
    factor = exp(-poly_value(kernel.aint, h) / 2)
    phase_plus = exp(kernel.xi * eta)
    phase_minus = exp(-kernel.xi * eta)
    chi = phase_plus * plus + phase_minus * minus
    dchi = phase_plus * (dplus + kernel.xi * deta * plus) +
        phase_minus * (dminus - kernel.xi * deta * minus)
    d2chi = phase_plus * (d2plus + 2kernel.xi * deta * dplus +
        (kernel.xi * d2eta + kernel.xi^2 * deta^2) * plus) +
        phase_minus * (d2minus - 2kernel.xi * deta * dminus +
        (-kernel.xi * d2eta + kernel.xi^2 * deta^2) * minus)
    psi = Z * chi
    dpsi = dZ * chi + Z * dchi
    d2psi = d2Z * chi + 2dZ * dchi + Z * d2chi
    value = factor * psi
    derivative = factor * (dpsi - A * psi / 2)
    second = factor * (d2psi - A * dpsi + (-dA / 2 + A^2 / 4) * psi)
    return value, derivative, second
end

@inline function _mode_pair(cplus, cminus, plus, dplus, minus, dminus,
        phase_plus, phase_minus, xi, deta)
    plus_value = cplus * plus
    plus_derivative = cplus * dplus
    minus_value = cminus * minus
    minus_derivative = cminus * dminus
    chi = phase_plus * plus_value + phase_minus * minus_value
    dchi = phase_plus * (plus_derivative + xi * deta * plus_value) +
        phase_minus * (minus_derivative - xi * deta * minus_value)
    return chi, dchi
end

@inline function _mode_triple(cplus, cminus,
        plus, dplus, d2plus, minus, dminus, d2minus,
        phase_plus, phase_minus, xi, deta, d2eta)
    plus_value = cplus * plus
    plus_derivative = cplus * dplus
    plus_second = cplus * d2plus
    minus_value = cminus * minus
    minus_derivative = cminus * dminus
    minus_second = cminus * d2minus
    chi = phase_plus * plus_value + phase_minus * minus_value
    dchi = phase_plus * (plus_derivative + xi * deta * plus_value) +
        phase_minus * (minus_derivative - xi * deta * minus_value)
    d2chi = phase_plus * (plus_second + 2xi * deta * plus_derivative +
        (xi * d2eta + xi^2 * deta^2) * plus_value) +
        phase_minus * (minus_second - 2xi * deta * minus_derivative +
        (-xi * d2eta + xi^2 * deta^2) * minus_value)
    return chi, dchi, d2chi
end

# Two kernels built from the same plan: the series are evaluated once from
# kernel1 and only the mode coefficients cplus/cminus differ.
function _kernel_pair2(kernel1::EikonalKernel, kernel2::EikonalKernel,
        h, degree=length(kernel1.A) - 1)
    eta, deta = _poly_pair(kernel1.eta, h, degree)
    Z, dZ = _poly_pair(kernel1.Z, h, degree)
    plus, dplus = _poly_pair(kernel1.plus, h, degree)
    minus, dminus = _poly_pair(kernel1.minus, h, degree)
    A = poly_value(
        @view(kernel1.A[1:min(end, degree + 1)]), h)
    aint = poly_value(
        @view(kernel1.aint[1:min(end, degree + 1)]), h)
    factor = exp(-aint / 2)
    phase_plus = exp(kernel1.xi * eta)
    phase_minus = exp(-kernel1.xi * eta)
    chi1, dchi1 = _mode_pair(
        kernel1.cplus, kernel1.cminus, plus, dplus, minus, dminus,
        phase_plus, phase_minus, kernel1.xi, deta)
    chi2, dchi2 = _mode_pair(
        kernel2.cplus, kernel2.cminus, plus, dplus, minus, dminus,
        phase_plus, phase_minus, kernel1.xi, deta)
    value1 = factor * Z * chi1
    derivative1 = factor * (dZ * chi1 + Z * dchi1 - A * Z * chi1 / 2)
    value2 = factor * Z * chi2
    derivative2 = factor * (dZ * chi2 + Z * dchi2 - A * Z * chi2 / 2)
    return (value1, derivative1), (value2, derivative2)
end

function _kernel_triple2(kernel1::EikonalKernel, kernel2::EikonalKernel, h)
    eta, deta, d2eta = poly_triple(kernel1.eta, h)
    Z, dZ, d2Z = poly_triple(kernel1.Z, h)
    plus, dplus, d2plus = poly_triple(kernel1.plus, h)
    minus, dminus, d2minus = poly_triple(kernel1.minus, h)
    A, dA = poly_pair(kernel1.A, h)
    factor = exp(-poly_value(kernel1.aint, h) / 2)
    phase_plus = exp(kernel1.xi * eta)
    phase_minus = exp(-kernel1.xi * eta)
    chi1, dchi1, d2chi1 = _mode_triple(
        kernel1.cplus, kernel1.cminus,
        plus, dplus, d2plus, minus, dminus, d2minus,
        phase_plus, phase_minus, kernel1.xi, deta, d2eta)
    chi2, dchi2, d2chi2 = _mode_triple(
        kernel2.cplus, kernel2.cminus,
        plus, dplus, d2plus, minus, dminus, d2minus,
        phase_plus, phase_minus, kernel1.xi, deta, d2eta)
    psi1 = Z * chi1
    dpsi1 = dZ * chi1 + Z * dchi1
    d2psi1 = d2Z * chi1 + 2dZ * dchi1 + Z * d2chi1
    value1 = factor * psi1
    derivative1 = factor * (dpsi1 - A * psi1 / 2)
    second1 = factor * (
        d2psi1 - A * dpsi1 + (-dA / 2 + A^2 / 4) * psi1)
    psi2 = Z * chi2
    dpsi2 = dZ * chi2 + Z * dchi2
    d2psi2 = d2Z * chi2 + 2dZ * dchi2 + Z * d2chi2
    value2 = factor * psi2
    derivative2 = factor * (dpsi2 - A * psi2 / 2)
    second2 = factor * (
        d2psi2 - A * dpsi2 + (-dA / 2 + A^2 / 4) * psi2)
    return (value1, derivative1, second1),
        (value2, derivative2, second2)
end

function _pair_error(candidate, reference)
    scale = max(abs(candidate[1]), abs(candidate[2]),
        abs(reference[1]), abs(reference[2]), floatmin(Float64))
    return max(abs(candidate[1] - reference[1]),
        abs(candidate[2] - reference[2])) / scale
end

# State at offset h with two error measures: the relative residual of
# X'' + A X' + B X = 0, and the change when the series are cut
# EIKONAL_CHECK_GAP orders lower.
function _kernel_step(kernel::EikonalKernel, h)
    value, derivative, second = _kernel_triple(kernel, h)
    A = poly_value(kernel.A, h)
    B = poly_value(kernel.B, h)
    residual = second + A * derivative + B * value
    residual_scale = abs(second) + abs(A * derivative) + abs(B * value)
    score = abs(residual) / max(residual_scale, eps(Float64))
    check_order = max(4, length(kernel.A) - 1 - EIKONAL_CHECK_GAP)
    checked = _kernel_pair(kernel, h, check_order)
    check_error = _pair_error((value, derivative), checked)
    return (; value, derivative, score, check_error)
end

function _kernel_steps(kernel1::EikonalKernel, kernel2::EikonalKernel, h)
    triple1, triple2 = _kernel_triple2(kernel1, kernel2, h)
    value1, derivative1, second1 = triple1
    value2, derivative2, second2 = triple2
    A = poly_value(kernel1.A, h)
    B = poly_value(kernel1.B, h)
    residual1 = second1 + A * derivative1 + B * value1
    residual2 = second2 + A * derivative2 + B * value2
    scale1 = abs(second1) + abs(A * derivative1) + abs(B * value1)
    scale2 = abs(second2) + abs(A * derivative2) + abs(B * value2)
    score1 = abs(residual1) / max(scale1, eps(Float64))
    score2 = abs(residual2) / max(scale2, eps(Float64))
    check_order = max(4, length(kernel1.A) - 1 - EIKONAL_CHECK_GAP)
    checked1, checked2 = _kernel_pair2(
        kernel1, kernel2, h, check_order)
    check1 = _pair_error((value1, derivative1), checked1)
    check2 = _pair_error((value2, derivative2), checked2)
    return (
        (; value=value1, derivative=derivative1,
            score=score1, check_error=check1),
        (; value=value2, derivative=derivative2,
            score=score2, check_error=check2),
    )
end

# Shrink the step by 0.8 until the residual and truncation checks pass; returns
# nothing below min_step. `rejected` records whether any shrink happened.
function _accepted_single(kernel, direction, h_abs, min_step,
    residual_limit, check_limit)
    rejected = false
    while h_abs >= min_step
        h = direction * h_abs
        step = _kernel_step(kernel, h)
        if isfinite(step.score) && isfinite(step.check_error) &&
                step.score <= residual_limit && step.check_error <= check_limit
            return (; h, step, residual=step.score,
                check_error=step.check_error, rejected)
        end
        @static if EIKONAL_PROFILE
            _profile_reject!(_reject_reason(
                step.score, step.check_error, residual_limit, check_limit))
        end
        rejected = true
        h_abs *= 0.8
    end
    return nothing
end

function _accepted_pair(kernel1, kernel2, direction, h_abs, min_step,
    residual_limit, check_limit)
    rejected = false
    while h_abs >= min_step
        h = direction * h_abs
        step1, step2 = _kernel_steps(kernel1, kernel2, h)
        residual = max(step1.score, step2.score)
        check_error = max(step1.check_error, step2.check_error)
        if isfinite(residual) && isfinite(check_error) &&
                residual <= residual_limit && check_error <= check_limit
            return (; h, step1, step2, residual, check_error, rejected)
        end
        @static if EIKONAL_PROFILE
            _profile_reject!(_reject_reason(
                residual, check_error, residual_limit, check_limit))
        end
        rejected = true
        h_abs *= 0.8
    end
    return nothing
end

function _patch(kernel::EikonalKernel, center_x, next_x, scale)
    return EikonalPatch(
        center_x,
        next_x,
        kernel.xi,
        kernel.A,
        kernel.B,
        kernel.eta,
        kernel.Z,
        kernel.aint,
        kernel.plus,
        kernel.minus,
        kernel.cplus,
        kernel.cminus,
        scale,
    )
end

function _copy_patch(kernel::EikonalKernel, center_x, next_x, scale)
    return EikonalPatch(
        center_x,
        next_x,
        kernel.xi,
        copy(kernel.A),
        copy(kernel.B),
        copy(kernel.eta),
        copy(kernel.Z),
        copy(kernel.aint),
        copy(kernel.plus),
        copy(kernel.minus),
        kernel.cplus,
        kernel.cminus,
        scale,
    )
end

function _copy_patch_pair(kernel1::EikonalKernel, kernel2::EikonalKernel,
        center_x, next_x, scale1, scale2)
    A = copy(kernel1.A)
    B = copy(kernel1.B)
    eta = copy(kernel1.eta)
    Z = copy(kernel1.Z)
    aint = copy(kernel1.aint)
    plus = copy(kernel1.plus)
    minus = copy(kernel1.minus)
    patch1 = EikonalPatch(
        center_x, next_x, kernel1.xi, A, B, eta, Z, aint, plus, minus,
        kernel1.cplus, kernel1.cminus, scale1)
    patch2 = EikonalPatch(
        center_x, next_x, kernel2.xi, A, B, eta, Z, aint, plus, minus,
        kernel2.cplus, kernel2.cminus, scale2)
    return patch1, patch2
end

# Propagate (seed_X, seed_dXdx) from seed_x to target_x with eikonal patches,
# re-planning at each patch start. Returns patches, the final state and step
# statistics, or nothing as soon as a plan or an accepted step is unavailable.
function eikonal_path(
    coefficients::CoefficientSet,
    seed_x,
    seed_X,
    seed_dXdx,
    target_x;
    avec,
    bvec,
    pq,
    scratch=EikonalScratch(),
    residual_limit=EIKONAL_RESIDUAL_LIMIT,
    check_limit=EIKONAL_CHECK_LIMIT,
    terms::Integer=EIKONAL_TERMS,
    large_frequency_expansion::Union{Nothing,LargeFrequencyExpansionPlan}=nothing,
)
    current_x = Float64(seed_x)
    target = Float64(target_x)
    value = ComplexF64(seed_X)
    derivative = ComplexF64(seed_dXdx)
    direction = target >= current_x ? 1.0 : -1.0
    patches = EikonalPatch[]
    sizehint!(patches, 12)
    max_step = 0.0
    min_step = Inf
    max_residual = 0.0
    max_check = 0.0
    order = length(scratch.f) - 1
    next_h_abs = Inf
    while direction * (target - current_x) > 100eps(Float64)
        a_terms, b_terms = large_frequency_expansion === nothing ?
            ordinary_ab_series!(
                avec, bvec, pq, coefficients, current_x, order) :
            large_frequency_expansion_ab!(avec, bvec, large_frequency_expansion, current_x, order)
        @static if EIKONAL_PROFILE
            plan = _profile_plan(@view(avec[1:a_terms]),
                @view(bvec[1:b_terms]), scratch, terms)
        else
            plan = _eikonal_plan(@view(avec[1:a_terms]),
                @view(bvec[1:b_terms]), scratch, terms)
        end
        plan === nothing && return nothing
        state_scale = max(abs(value), abs(derivative), floatmin(Float64))
        kernel = _eikonal_kernel(plan, value / state_scale, derivative / state_scale)
        kernel === nothing && return nothing
        remaining = abs(target - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_step = EIKONAL_RADIUS_SAFETY * singular_distance
        step_cap = large_frequency_expansion === nothing ?
            min(remaining, radius_step, EIKONAL_MAX_STEP) :
            min(remaining, radius_step)
        h_abs = min(step_cap, next_h_abs)
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        best = _accepted_single(kernel, direction, h_abs, min_abs_step,
            residual_limit, check_limit)
        best === nothing && return nothing
        next_x = current_x + best.h
        push!(patches, _patch(kernel, current_x, next_x, state_scale))
        @static if EIKONAL_PROFILE
            _profile_accept!(1)
        end
        value = state_scale * best.step.value
        derivative = state_scale * best.step.derivative
        current_x = next_x
        next_h_abs = (best.rejected ? 1.0 : EIKONAL_STEP_GROWTH) * abs(best.h)
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))
        max_residual = max(max_residual, best.residual)
        max_check = max(max_check, best.check_error)
    end
    return (; patches, value, derivative, max_step,
        min_step=isfinite(min_step) ? min_step : 0.0,
        max_residual, max_check)
end

# Like `eikonal_path`, but each step of length 2h is checked against two
# half-steps, the second re-planned at the midpoint; the step is accepted only
# when the full and split states agree within `split_limit`. With
# `allow_partial`, a failure returns the patches built so far (with `x`).
function eikonal_split(
    coefficients::CoefficientSet,
    seed_x,
    seed_X,
    seed_dXdx,
    target_x;
    avec,
    bvec,
    pq,
    scratch=EikonalScratch(),
    residual_limit=EIKONAL_RESIDUAL_LIMIT,
    check_limit=EIKONAL_CHECK_LIMIT,
    split_limit=EIKONAL_SPLIT_LIMIT,
    terms::Integer=EIKONAL_TERMS,
    large_frequency_expansion::Union{Nothing,LargeFrequencyExpansionPlan}=nothing,
    allow_partial::Bool=false,
)
    current_x = Float64(seed_x)
    target = Float64(target_x)
    value = ComplexF64(seed_X)
    derivative = ComplexF64(seed_dXdx)
    direction = target >= current_x ? 1.0 : -1.0
    patches = EikonalPatch[]
    sizehint!(patches, 24)
    max_step = 0.0
    min_step = Inf
    max_residual = 0.0
    max_check = 0.0
    max_split = 0.0
    order = length(scratch.f) - 1
    middle_scratch = EikonalScratch(order)
    next_h_abs = Inf
    partial_result() = isempty(patches) ? nothing : (
        patches=patches,
        value=value,
        derivative=derivative,
        x=current_x,
        max_step=max_step,
        min_step=isfinite(min_step) ? min_step : 0.0,
        max_residual=max_residual,
        max_check=max_check,
        max_split=max_split,
    )
    while direction * (target - current_x) > 100eps(Float64)
        remaining = abs(target - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_step = EIKONAL_RADIUS_SAFETY * singular_distance
        step_cap = large_frequency_expansion === nothing ?
            min(remaining, radius_step, EIKONAL_MAX_STEP) :
            min(remaining, radius_step)
        h_abs = min(step_cap, next_h_abs)
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        a_terms, b_terms = large_frequency_expansion === nothing ?
            ordinary_ab_series!(
                avec, bvec, pq, coefficients, current_x, order) :
            large_frequency_expansion_ab!(avec, bvec, large_frequency_expansion, current_x, order)
        @static if EIKONAL_PROFILE
            plan1 = _profile_plan!(@view(avec[1:a_terms]),
                @view(bvec[1:b_terms]), scratch, terms)
        else
            plan1 = _eikonal_plan!(@view(avec[1:a_terms]),
                @view(bvec[1:b_terms]), scratch, terms)
        end
        plan1 === nothing && return allow_partial ? partial_result() : nothing
        scale1 = max(abs(value), abs(derivative), floatmin(Float64))
        kernel1 = _eikonal_kernel(
            plan1, value / scale1, derivative / scale1)
        kernel1 === nothing && return nothing
        accepted = nothing
        rejected = false
        while h_abs >= min_abs_step
                half_h = direction * h_abs / 2
                full = _kernel_step(kernel1, 2 * half_h)
                if !isfinite(full.score) || !isfinite(full.check_error) ||
                        full.score > residual_limit ||
                        full.check_error > check_limit
                    @static if EIKONAL_PROFILE
                        _profile_reject!(_reject_reason(full.score,
                            full.check_error, residual_limit, check_limit))
                    end
                    rejected = true
                    h_abs *= 0.8
                    continue
                end
                step1 = _kernel_step(kernel1, half_h)
                if !isfinite(step1.score) || !isfinite(step1.check_error) ||
                        step1.score > residual_limit ||
                        step1.check_error > check_limit
                    @static if EIKONAL_PROFILE
                        _profile_reject!(_reject_reason(step1.score,
                            step1.check_error, residual_limit, check_limit))
                    end
                    rejected = true
                    h_abs *= 0.8
                    continue
                end
                midpoint_x = current_x + half_h
                midpoint_value = scale1 * step1.value
                midpoint_derivative = scale1 * step1.derivative

                middle_a, middle_b = large_frequency_expansion === nothing ?
                    ordinary_ab_series!(avec, bvec, pq,
                        coefficients, midpoint_x, order) :
                    large_frequency_expansion_ab!(avec, bvec, large_frequency_expansion, midpoint_x, order)
                @static if EIKONAL_PROFILE
                    plan2 = _profile_plan!(@view(avec[1:middle_a]),
                        @view(bvec[1:middle_b]), middle_scratch, terms)
                else
                    plan2 = _eikonal_plan!(@view(avec[1:middle_a]),
                        @view(bvec[1:middle_b]), middle_scratch, terms)
                end
                if plan2 === nothing
                    @static if EIKONAL_PROFILE
                        _profile_reject!(:other)
                    end
                    rejected = true
                    h_abs *= 0.8
                    continue
                end
                scale2 = max(abs(midpoint_value), abs(midpoint_derivative),
                    floatmin(Float64))
                kernel2 = _eikonal_kernel(plan2,
                    midpoint_value / scale2, midpoint_derivative / scale2)
                if kernel2 === nothing
                    @static if EIKONAL_PROFILE
                        _profile_reject!(:other)
                    end
                    rejected = true
                    h_abs *= 0.8
                    continue
                end
                step2 = _kernel_step(kernel2, half_h)
                full_state = (scale1 * full.value, scale1 * full.derivative)
                split_state = (scale2 * step2.value, scale2 * step2.derivative)
                split_error = _pair_error(full_state, split_state)
                residual = max(full.score, step1.score, step2.score)
                check_error = max(
                    full.check_error, step1.check_error, step2.check_error)
                if isfinite(residual) && isfinite(check_error) &&
                        isfinite(split_error) &&
                        residual <= residual_limit && check_error <= check_limit &&
                        split_error <= split_limit
                    accepted = (; half_h, midpoint_x, scale1, scale2,
                        kernel1, kernel2, step2, residual, check_error, split_error)
                    break
                end
                @static if EIKONAL_PROFILE
                    _profile_reject!(_reject_reason(residual, check_error,
                        split_error, residual_limit, check_limit, split_limit))
                end
                rejected = true
                h_abs *= 0.8
        end
        accepted === nothing && return allow_partial ? partial_result() : nothing
        if EIKONAL_PARTIAL_STALL && allow_partial && 2abs(accepted.half_h) <=
                EIKONAL_PARTIAL_STALL_RATIO * step_cap
            return partial_result()
        end
        next_x = current_x + 2 * accepted.half_h
        push!(patches, _copy_patch(accepted.kernel1, current_x,
            accepted.midpoint_x, accepted.scale1))
        push!(patches, _copy_patch(accepted.kernel2, accepted.midpoint_x,
            next_x, accepted.scale2))
        @static if EIKONAL_PROFILE
            _profile_accept!(2)
        end
        value = accepted.scale2 * accepted.step2.value
        derivative = accepted.scale2 * accepted.step2.derivative
        current_x = next_x
        next_h_abs = 2 * (rejected ? 1.0 : EIKONAL_STEP_GROWTH) *
            abs(accepted.half_h)
        max_step = max(max_step, abs(accepted.half_h))
        min_step = min(min_step, abs(accepted.half_h))
        max_residual = max(max_residual, accepted.residual)
        max_check = max(max_check, accepted.check_error)
        max_split = max(max_split, accepted.split_error)
    end
    return (; patches, value, derivative, x=current_x, max_step,
        min_step=isfinite(min_step) ? min_step : 0.0,
        max_residual, max_check, max_split)
end

# Paired form of `eikonal_split` with a common step size for both solutions.
# Returns nothing (no partial result) when any step fails.
function eikonal_split_pair(
    coefficients::CoefficientSet,
    seed_x,
    seed_X1,
    seed_dXdx1,
    seed_X2,
    seed_dXdx2,
    target_x;
    avec,
    bvec,
    pq,
    scratch=EikonalScratch(),
    residual_limit=EIKONAL_RESIDUAL_LIMIT,
    check_limit=EIKONAL_CHECK_LIMIT,
    split_limit=EIKONAL_SPLIT_LIMIT,
    terms::Integer=EIKONAL_TERMS,
    large_frequency_expansion::Union{Nothing,LargeFrequencyExpansionPlan}=nothing,
)
    current_x = Float64(seed_x)
    target = Float64(target_x)
    value1 = ComplexF64(seed_X1)
    deriv1 = ComplexF64(seed_dXdx1)
    value2 = ComplexF64(seed_X2)
    deriv2 = ComplexF64(seed_dXdx2)
    direction = target >= current_x ? 1.0 : -1.0
    patches1 = EikonalPatch[]
    patches2 = EikonalPatch[]
    sizehint!(patches1, 24)
    sizehint!(patches2, 24)
    max_step = 0.0
    min_step = Inf
    max_residual = 0.0
    max_check = 0.0
    max_split = 0.0
    order = length(scratch.f) - 1
    middle_scratch = EikonalScratch(order)
    next_h_abs = Inf
    while direction * (target - current_x) > 100eps(Float64)
        remaining = abs(target - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_step = EIKONAL_RADIUS_SAFETY * singular_distance
        step_cap = large_frequency_expansion === nothing ?
            min(remaining, radius_step, EIKONAL_MAX_STEP) :
            min(remaining, radius_step)
        h_abs = min(step_cap, next_h_abs)
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        a_terms, b_terms = large_frequency_expansion === nothing ?
            ordinary_ab_series!(
                avec, bvec, pq, coefficients, current_x, order) :
            large_frequency_expansion_ab!(avec, bvec, large_frequency_expansion, current_x, order)
        @static if EIKONAL_PROFILE
            plan1 = _profile_plan!(@view(avec[1:a_terms]),
                @view(bvec[1:b_terms]), scratch, terms)
        else
            plan1 = _eikonal_plan!(@view(avec[1:a_terms]),
                @view(bvec[1:b_terms]), scratch, terms)
        end
        plan1 === nothing && return nothing
        scale11 = max(abs(value1), abs(deriv1), floatmin(Float64))
        scale12 = max(abs(value2), abs(deriv2), floatmin(Float64))
        kernel11 = _eikonal_kernel(
            plan1, value1 / scale11, deriv1 / scale11)
        kernel12 = _eikonal_kernel(
            plan1, value2 / scale12, deriv2 / scale12)
        (kernel11 === nothing || kernel12 === nothing) && return nothing
        accepted = nothing
        rejected = false
        while h_abs >= min_abs_step
                half_h = direction * h_abs / 2
                full1, full2 = _kernel_steps(
                    kernel11, kernel12, 2 * half_h)
                full_residual = max(full1.score, full2.score)
                full_check = max(full1.check_error, full2.check_error)
                if !isfinite(full_residual) || !isfinite(full_check) ||
                        full_residual > residual_limit || full_check > check_limit
                    @static if EIKONAL_PROFILE
                        _profile_reject!(_reject_reason(full_residual,
                            full_check, residual_limit, check_limit))
                    end
                    rejected = true
                    h_abs *= 0.8
                    continue
                end
                first1, first2 = _kernel_steps(kernel11, kernel12, half_h)
                first_residual = max(first1.score, first2.score)
                first_check = max(first1.check_error, first2.check_error)
                if !isfinite(first_residual) || !isfinite(first_check) ||
                        first_residual > residual_limit || first_check > check_limit
                    @static if EIKONAL_PROFILE
                        _profile_reject!(_reject_reason(first_residual,
                            first_check, residual_limit, check_limit))
                    end
                    rejected = true
                    h_abs *= 0.8
                    continue
                end
                midpoint_x = current_x + half_h
                midpoint_value1 = scale11 * first1.value
                midpoint_deriv1 = scale11 * first1.derivative
                midpoint_value2 = scale12 * first2.value
                midpoint_deriv2 = scale12 * first2.derivative

                middle_a, middle_b = large_frequency_expansion === nothing ?
                    ordinary_ab_series!(avec, bvec, pq,
                        coefficients, midpoint_x, order) :
                    large_frequency_expansion_ab!(avec, bvec, large_frequency_expansion, midpoint_x, order)
                @static if EIKONAL_PROFILE
                    plan2 = _profile_plan!(@view(avec[1:middle_a]),
                        @view(bvec[1:middle_b]), middle_scratch, terms)
                else
                    plan2 = _eikonal_plan!(@view(avec[1:middle_a]),
                        @view(bvec[1:middle_b]), middle_scratch, terms)
                end
                if plan2 === nothing
                    @static if EIKONAL_PROFILE
                        _profile_reject!(:other)
                    end
                    rejected = true
                    h_abs *= 0.8
                    continue
                end
                scale21 = max(abs(midpoint_value1), abs(midpoint_deriv1),
                    floatmin(Float64))
                scale22 = max(abs(midpoint_value2), abs(midpoint_deriv2),
                    floatmin(Float64))
                kernel21 = _eikonal_kernel(plan2,
                    midpoint_value1 / scale21, midpoint_deriv1 / scale21)
                kernel22 = _eikonal_kernel(plan2,
                    midpoint_value2 / scale22, midpoint_deriv2 / scale22)
                if kernel21 === nothing || kernel22 === nothing
                    @static if EIKONAL_PROFILE
                        _profile_reject!(:other)
                    end
                    rejected = true
                    h_abs *= 0.8
                    continue
                end
                second1, second2 = _kernel_steps(kernel21, kernel22, half_h)
                split1 = _pair_error(
                    (scale11 * full1.value, scale11 * full1.derivative),
                    (scale21 * second1.value, scale21 * second1.derivative))
                split2 = _pair_error(
                    (scale12 * full2.value, scale12 * full2.derivative),
                    (scale22 * second2.value, scale22 * second2.derivative))
                split_error = max(split1, split2)
                residual = max(full1.score, full2.score,
                    first1.score, first2.score, second1.score, second2.score)
                check_error = max(full1.check_error, full2.check_error,
                    first1.check_error, first2.check_error,
                    second1.check_error, second2.check_error)
                if isfinite(residual) && isfinite(check_error) &&
                        isfinite(split_error) &&
                        residual <= residual_limit && check_error <= check_limit &&
                        split_error <= split_limit
                    accepted = (; half_h, midpoint_x,
                        scale11, scale12, scale21, scale22,
                        kernel11, kernel12, kernel21, kernel22,
                        second1, second2, residual, check_error, split_error)
                    break
                end
                @static if EIKONAL_PROFILE
                    _profile_reject!(_reject_reason(residual, check_error,
                        split_error, residual_limit, check_limit, split_limit))
                end
                rejected = true
                h_abs *= 0.8
        end
        accepted === nothing && return nothing
        next_x = current_x + 2 * accepted.half_h
        first1, first2 = _copy_patch_pair(
            accepted.kernel11, accepted.kernel12, current_x,
            accepted.midpoint_x, accepted.scale11, accepted.scale12)
        second1, second2 = _copy_patch_pair(
            accepted.kernel21, accepted.kernel22, accepted.midpoint_x,
            next_x, accepted.scale21, accepted.scale22)
        push!(patches1, first1)
        push!(patches2, first2)
        push!(patches1, second1)
        push!(patches2, second2)
        @static if EIKONAL_PROFILE
            _profile_accept!(2)
        end
        value1 = accepted.scale21 * accepted.second1.value
        deriv1 = accepted.scale21 * accepted.second1.derivative
        value2 = accepted.scale22 * accepted.second2.value
        deriv2 = accepted.scale22 * accepted.second2.derivative
        current_x = next_x
        next_h_abs = 2 * (rejected ? 1.0 : EIKONAL_STEP_GROWTH) *
            abs(accepted.half_h)
        max_step = max(max_step, abs(accepted.half_h))
        min_step = min(min_step, abs(accepted.half_h))
        max_residual = max(max_residual, accepted.residual)
        max_check = max(max_check, accepted.check_error)
        max_split = max(max_split, accepted.split_error)
    end
    return (; patches1, patches2, value1, deriv1, value2, deriv2,
        max_step, min_step=isfinite(min_step) ? min_step : 0.0,
        max_residual, max_check, max_split)
end

# Paired form of `eikonal_path`: two solutions share one plan per step (they
# differ only in cplus/cminus) and advance with a common step size.
function eikonal_pair(
    coefficients::CoefficientSet,
    seed_x,
    seed_X1,
    seed_dXdx1,
    seed_X2,
    seed_dXdx2,
    target_x;
    avec,
    bvec,
    pq,
    scratch=EikonalScratch(),
    residual_limit=EIKONAL_RESIDUAL_LIMIT,
    check_limit=EIKONAL_CHECK_LIMIT,
    terms::Integer=EIKONAL_TERMS,
    large_frequency_expansion::Union{Nothing,LargeFrequencyExpansionPlan}=nothing,
)
    current_x = Float64(seed_x)
    target = Float64(target_x)
    value1 = ComplexF64(seed_X1)
    deriv1 = ComplexF64(seed_dXdx1)
    value2 = ComplexF64(seed_X2)
    deriv2 = ComplexF64(seed_dXdx2)
    direction = target >= current_x ? 1.0 : -1.0
    patches1 = EikonalPatch[]
    patches2 = EikonalPatch[]
    sizehint!(patches1, 12)
    sizehint!(patches2, 12)
    max_step = 0.0
    min_step = Inf
    max_residual = 0.0
    max_check = 0.0
    order = length(scratch.f) - 1
    next_h_abs = Inf
    while direction * (target - current_x) > 100eps(Float64)
        a_terms, b_terms = large_frequency_expansion === nothing ?
            ordinary_ab_series!(
                avec, bvec, pq, coefficients, current_x, order) :
            large_frequency_expansion_ab!(avec, bvec, large_frequency_expansion, current_x, order)
        @static if EIKONAL_PROFILE
            plan = _profile_plan(@view(avec[1:a_terms]),
                @view(bvec[1:b_terms]), scratch, terms)
        else
            plan = _eikonal_plan(@view(avec[1:a_terms]),
                @view(bvec[1:b_terms]), scratch, terms)
        end
        plan === nothing && return nothing
        scale1 = max(abs(value1), abs(deriv1), floatmin(Float64))
        scale2 = max(abs(value2), abs(deriv2), floatmin(Float64))
        kernel1 = _eikonal_kernel(plan, value1 / scale1, deriv1 / scale1)
        kernel2 = _eikonal_kernel(plan, value2 / scale2, deriv2 / scale2)
        (kernel1 === nothing || kernel2 === nothing) && return nothing
        remaining = abs(target - current_x)
        singular_distance = min(abs(current_x), abs(1.0 - current_x))
        radius_step = EIKONAL_RADIUS_SAFETY * singular_distance
        step_cap = large_frequency_expansion === nothing ?
            min(remaining, radius_step, EIKONAL_MAX_STEP) :
            min(remaining, radius_step)
        h_abs = min(step_cap, next_h_abs)
        min_abs_step = 100eps(Float64) * max(1.0, abs(current_x))
        best = _accepted_pair(kernel1, kernel2, direction, h_abs, min_abs_step,
            residual_limit, check_limit)
        best === nothing && return nothing
        next_x = current_x + best.h
        push!(patches1, _patch(kernel1, current_x, next_x, scale1))
        push!(patches2, _patch(kernel2, current_x, next_x, scale2))
        @static if EIKONAL_PROFILE
            _profile_accept!(1)
        end
        value1 = scale1 * best.step1.value
        deriv1 = scale1 * best.step1.derivative
        value2 = scale2 * best.step2.value
        deriv2 = scale2 * best.step2.derivative
        current_x = next_x
        next_h_abs = (best.rejected ? 1.0 : EIKONAL_STEP_GROWTH) * abs(best.h)
        max_step = max(max_step, abs(best.h))
        min_step = min(min_step, abs(best.h))
        max_residual = max(max_residual, best.residual)
        max_check = max(max_check, best.check_error)
    end
    return (; patches1, patches2, value1, deriv1, value2, deriv2, max_step,
        min_step=isfinite(min_step) ? min_step : 0.0,
        max_residual, max_check)
end

function _patch_kernel(patch::EikonalPatch)
    return EikonalKernel(
        patch.A,
        patch.B,
        patch.xi,
        patch.eta,
        patch.Z,
        patch.aint,
        patch.plus,
        patch.minus,
        patch.cplus,
        patch.cminus,
    )
end

function eikonal_value(patch::EikonalPatch, x)
    h = x - patch.center_x
    kernel = _patch_kernel(patch)
    eta = poly_value(kernel.eta, h)
    Z = poly_value(kernel.Z, h)
    plus = kernel.cplus * poly_value(kernel.plus, h)
    minus = kernel.cminus * poly_value(kernel.minus, h)
    factor = exp(-poly_value(kernel.aint, h) / 2)
    value = factor * Z * (exp(kernel.xi * eta) * plus +
        exp(-kernel.xi * eta) * minus)
    return patch.scale * value
end

function eikonal_state(patch::EikonalPatch, x)
    value, derivative = _kernel_pair(_patch_kernel(patch), x - patch.center_x)
    return patch.scale * value, patch.scale * derivative
end

end
