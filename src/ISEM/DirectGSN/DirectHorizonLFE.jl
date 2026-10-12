module DirectHorizonLFE

using ..DirectCoefficientTables: DirectCoefficientSet
using ..DirectLocalSolutionAtZero:
    direct_horizon_series_triple,
    direct_zero_local_solution
using ..DirectOrdinaryPointExpansion:
    direct_ordinary_ab_series,
    direct_poly_value

export HorizonLFEPatch
export lfeh_single, lfeh_pair, lfeh_value, lfeh_state

# Horizon patches for the large-frequency route: the horizon Frobenius series is
# used directly from the seed point out to the largest x below the matching point
# where it passes the tail, ODE-residual and nested-order checks.

const LFEH_MIN_ORDER = 24
const LFEH_MAX_ORDER = 96
const LFEH_ORDER_STEP = 8
const LFEH_TAIL = 1e-10
const LFEH_RESIDUAL = 1e-14
const LFEH_CHECK = 1e-10
const LFEH_SAFETY = 0.8

struct HorizonLFEPatch
    center_x::Float64
    next_x::Float64
    exponent::ComplexF64
    coeffs::Vector{ComplexF64}
    scale::ComplexF64
end

function _state(solution, x, order)
    coeffs = @view solution.coefficients[1:(order + 1)]
    X, dXdx, d2Xdx2 = direct_horizon_series_triple(
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
    series = direct_poly_value(
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
    avec, bvec = direct_ordinary_ab_series(coefficients, x, 2)
    return avec[1], bvec[1]
end

# Lowest order (from LFEH_MIN_ORDER in steps of LFEH_ORDER_STEP) at which the
# series passes all three checks at x, or `nothing`.
function _candidate(solution, x, A, B)
    highest = solution.effective_order - mod(
        solution.effective_order - LFEH_MIN_ORDER, LFEH_ORDER_STEP)
    highest >= LFEH_MIN_ORDER || return nothing
    for order in LFEH_MIN_ORDER:LFEH_ORDER_STEP:highest
        full = _state(solution, x, order)
        nested = _state(solution, x, order - LFEH_ORDER_STEP)
        tail = _tail(solution, x, order)
        residual = _residual(full, A, B)
        check = _error(full, nested)
        if tail <= LFEH_TAIL && residual <= LFEH_RESIDUAL &&
                check <= LFEH_CHECK
            return (; state=full, order, tail, residual, check)
        end
    end
    return nothing
end

# Step the target down from match_x by LFEH_SAFETY until the series is accepted,
# then accept a patch only if it also passes one further safety step inward.
function _single(coefficients, solution, seed_x, scale, match_x)
    target = Float64(match_x)
    while target > seed_x * (1 + 100eps(Float64))
        A, B = _ab_value(coefficients, target)
        candidate = _candidate(solution, target, A, B)
        if candidate !== nothing
            safe_target = LFEH_SAFETY * target
            safe_target > seed_x * (1 + 100eps(Float64)) || return nothing
            A, B = _ab_value(coefficients, safe_target)
            safe = _candidate(solution, safe_target, A, B)
            if safe !== nothing
                patch = HorizonLFEPatch(
                    Float64(seed_x), safe_target,
                    ComplexF64(solution.exponent),
                    copy(@view(solution.coefficients[1:(safe.order + 1)])),
                    ComplexF64(scale))
                return (; patch, state=lfeh_state(patch, safe_target),
                    target=safe_target, order=safe.order, tail=safe.tail,
                    residual=safe.residual, check=safe.check)
            end
        end
        target *= LFEH_SAFETY
    end
    return nothing
end

function lfeh_single(coefficients::DirectCoefficientSet, branch,
        seed_x, scale, match_x)
    solution = direct_zero_local_solution(
        coefficients, branch, LFEH_MAX_ORDER)
    return _single(
        coefficients, solution, Float64(seed_x), scale, match_x)
end

# Two-branch version of lfeh_single; both branches must pass at the same end point.
function lfeh_pair(coefficients::DirectCoefficientSet,
        branch1, scale1, branch2, scale2, seed_x, match_x)
    solution1 = direct_zero_local_solution(
        coefficients, branch1, LFEH_MAX_ORDER)
    solution2 = direct_zero_local_solution(
        coefficients, branch2, LFEH_MAX_ORDER)
    target = Float64(match_x)
    while target > seed_x * (1 + 100eps(Float64))
        A, B = _ab_value(coefficients, target)
        first = _candidate(solution1, target, A, B)
        second = _candidate(solution2, target, A, B)
        if first !== nothing && second !== nothing
            safe_target = LFEH_SAFETY * target
            safe_target > seed_x * (1 + 100eps(Float64)) || return nothing
            A, B = _ab_value(coefficients, safe_target)
            first = _candidate(solution1, safe_target, A, B)
            second = _candidate(solution2, safe_target, A, B)
            if first !== nothing && second !== nothing
                patch1 = HorizonLFEPatch(
                    Float64(seed_x), safe_target,
                    ComplexF64(solution1.exponent),
                    copy(@view(solution1.coefficients[1:(first.order + 1)])),
                    ComplexF64(scale1))
                patch2 = HorizonLFEPatch(
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
        target *= LFEH_SAFETY
    end
    return nothing
end

function lfeh_state(patch::HorizonLFEPatch, x)
    X, dXdx, d2Xdx2 = direct_horizon_series_triple(
        patch.coeffs, patch.exponent, x)
    return (
        X=patch.scale * X,
        dXdx=patch.scale * dXdx,
        d2Xdx2=patch.scale * d2Xdx2,
    )
end

lfeh_value(patch::HorizonLFEPatch, x) = lfeh_state(patch, x).X

end
