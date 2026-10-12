module HorizonLargeFrequencyExpansion

using ..CoefficientTables: CoefficientSet
using ..LocalSolutionAtZero:
    horizon_series_triple,
    zero_local_solution
using ..OrdinaryPointExpansion:
    ordinary_ab_series,
    poly_value

export HorizonLargeFrequencyExpansionPatch
export lfeh_single, lfeh_pair, lfeh_value, lfeh_state

# Horizon patches for the large-frequency route: the horizon Frobenius series is
# used directly from the seed point out to the largest x below the matching point
# where it passes the tail, ODE-residual and nested-order checks.

const LARGE_FREQUENCY_EXPANSION_HORIZON_MIN_ORDER = 24
const LARGE_FREQUENCY_EXPANSION_HORIZON_MAX_ORDER = 96
const LARGE_FREQUENCY_EXPANSION_HORIZON_ORDER_STEP = 8
const LARGE_FREQUENCY_EXPANSION_HORIZON_TAIL = 1e-10
const LARGE_FREQUENCY_EXPANSION_HORIZON_RESIDUAL = 1e-14
const LARGE_FREQUENCY_EXPANSION_HORIZON_CHECK = 1e-10
const LARGE_FREQUENCY_EXPANSION_HORIZON_SAFETY = 0.8

struct HorizonLargeFrequencyExpansionPatch
    center_x::Float64
    next_x::Float64
    exponent::ComplexF64
    coeffs::Vector{ComplexF64}
    scale::ComplexF64
end

function _state(solution, x, order)
    coeffs = @view solution.coefficients[1:(order + 1)]
    X, dXdx, d2Xdx2 = horizon_series_triple(
        coeffs, solution.exponent, x)
    return (; X, dXdx, d2Xdx2)
end

function _error(left, right)
    scale = max(abs(left.X), abs(left.dXdx),
        abs(right.X), abs(right.dXdx), floatmin(Float64))
    return max(abs(left.X - right.X),
        abs(left.dXdx - right.dXdx)) / scale
end

# Size of the last `width` series terms relative to the series value.
function _tail(solution, x, order; width=4)
    value = 0.0
    @inbounds for degree in max(0, order - width + 1):order
        value += abs(solution.coefficients[degree + 1] * x^degree)
    end
    series = poly_value(
        @view(solution.coefficients[1:(order + 1)]), x)
    return value / max(abs(series), eps(Float64))
end

# Relative residual of X'' + A X' + B X = 0 at one point.
function _residual(candidate, A, B)
    second = candidate.d2Xdx2
    first = A * candidate.dXdx
    zeroth = B * candidate.X
    return abs(second + first + zeroth) /
        max(abs(second) + abs(first) + abs(zeroth), eps(Float64))
end

function _ab_value(coefficients, x)
    avec, bvec = ordinary_ab_series(coefficients, x, 2)
    return avec[1], bvec[1]
end

# Lowest order (from LARGE_FREQUENCY_EXPANSION_HORIZON_MIN_ORDER in steps of LARGE_FREQUENCY_EXPANSION_HORIZON_ORDER_STEP) at which the
# series passes all three checks at x, or `nothing`.
function _candidate(solution, x, A, B)
    highest = solution.effective_order - mod(
        solution.effective_order - LARGE_FREQUENCY_EXPANSION_HORIZON_MIN_ORDER, LARGE_FREQUENCY_EXPANSION_HORIZON_ORDER_STEP)
    highest >= LARGE_FREQUENCY_EXPANSION_HORIZON_MIN_ORDER || return nothing
    for order in LARGE_FREQUENCY_EXPANSION_HORIZON_MIN_ORDER:LFEH_ORDER_STEP:highest
        full = _state(solution, x, order)
        nested = _state(solution, x, order - LARGE_FREQUENCY_EXPANSION_HORIZON_ORDER_STEP)
        tail = _tail(solution, x, order)
        residual = _residual(full, A, B)
        check = _error(full, nested)
        if tail <= LARGE_FREQUENCY_EXPANSION_HORIZON_TAIL && residual <= LARGE_FREQUENCY_EXPANSION_HORIZON_RESIDUAL &&
                check <= LARGE_FREQUENCY_EXPANSION_HORIZON_CHECK
            return (; state=full, order, tail, residual, check)
        end
    end
    return nothing
end

# Step the target down from match_x by LARGE_FREQUENCY_EXPANSION_HORIZON_SAFETY until the series is accepted,
# then accept a patch only if it also passes one further safety step inward.
function _single(coefficients, solution, seed_x, scale, match_x)
    target = Float64(match_x)
    while target > seed_x * (1 + 100eps(Float64))
        A, B = _ab_value(coefficients, target)
        candidate = _candidate(solution, target, A, B)
        if candidate !== nothing
            safe_target = LARGE_FREQUENCY_EXPANSION_HORIZON_SAFETY * target
            safe_target > seed_x * (1 + 100eps(Float64)) || return nothing
            A, B = _ab_value(coefficients, safe_target)
            safe = _candidate(solution, safe_target, A, B)
            if safe !== nothing
                patch = HorizonLargeFrequencyExpansionPatch(
                    Float64(seed_x), safe_target,
                    ComplexF64(solution.exponent),
                    copy(@view(solution.coefficients[1:(safe.order + 1)])),
                    ComplexF64(scale))
                return (; patch, state=lfeh_state(patch, safe_target),
                    target=safe_target, order=safe.order, tail=safe.tail,
                    residual=safe.residual, check=safe.check)
            end
        end
        target *= LARGE_FREQUENCY_EXPANSION_HORIZON_SAFETY
    end
    return nothing
end

function lfeh_single(coefficients::CoefficientSet, branch,
        seed_x, scale, match_x)
    solution = zero_local_solution(
        coefficients, branch, LARGE_FREQUENCY_EXPANSION_HORIZON_MAX_ORDER)
    return _single(
        coefficients, solution, Float64(seed_x), scale, match_x)
end

# Two-branch version of lfeh_single; both branches must pass at the same end point.
function lfeh_pair(coefficients::CoefficientSet,
        branch1, scale1, branch2, scale2, seed_x, match_x)
    solution1 = zero_local_solution(
        coefficients, branch1, LARGE_FREQUENCY_EXPANSION_HORIZON_MAX_ORDER)
    solution2 = zero_local_solution(
        coefficients, branch2, LARGE_FREQUENCY_EXPANSION_HORIZON_MAX_ORDER)
    target = Float64(match_x)
    while target > seed_x * (1 + 100eps(Float64))
        A, B = _ab_value(coefficients, target)
        first = _candidate(solution1, target, A, B)
        second = _candidate(solution2, target, A, B)
        if first !== nothing && second !== nothing
            safe_target = LARGE_FREQUENCY_EXPANSION_HORIZON_SAFETY * target
            safe_target > seed_x * (1 + 100eps(Float64)) || return nothing
            A, B = _ab_value(coefficients, safe_target)
            first = _candidate(solution1, safe_target, A, B)
            second = _candidate(solution2, safe_target, A, B)
            if first !== nothing && second !== nothing
                patch1 = HorizonLargeFrequencyExpansionPatch(
                    Float64(seed_x), safe_target,
                    ComplexF64(solution1.exponent),
                    copy(@view(solution1.coefficients[1:(first.order + 1)])),
                    ComplexF64(scale1))
                patch2 = HorizonLargeFrequencyExpansionPatch(
                    Float64(seed_x), safe_target,
                    ComplexF64(solution2.exponent),
                    copy(@view(solution2.coefficients[1:(second.order + 1)])),
                    ComplexF64(scale2))
                return (
                    patch1,
                    patch2,
                    state1=lfeh_state(patch1, safe_target),
                    state2=lfeh_state(patch2, safe_target),
                    target=safe_target,
                    tail=max(first.tail, second.tail),
                    residual=max(first.residual, second.residual),
                    check=max(first.check, second.check),
                )
            end
        end
        target *= LARGE_FREQUENCY_EXPANSION_HORIZON_SAFETY
    end
    return nothing
end

function lfeh_state(patch::HorizonLargeFrequencyExpansionPatch, x)
    X, dXdx, d2Xdx2 = horizon_series_triple(
        patch.coeffs, patch.exponent, x)
    return (
        X=patch.scale * X,
        dXdx=patch.scale * dXdx,
        d2Xdx2=patch.scale * d2Xdx2,
    )
end

lfeh_value(patch::HorizonLargeFrequencyExpansionPatch, x) = lfeh_state(patch, x).X

end
