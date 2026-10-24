module Matching

include("Parameters.jl")
using .Parameters
using ...TeukolskyTransformation
include("Iteration.jl")
using .Iteration
include("AsymptoticAmplitudes.jl")
using .AsymptoticAmplitudes

export _Pin, _Pup

_TruncatioN = 40
_TOLERANCE = 1e-13
_xm = - 1.0
_rhom = 1.0
_MATCHING_ACCEPT_TOL = 1e-12
_MATCHING_WARNING_TOL = 1e-8
_MIN_MATCHING_N = 10

function _warn_matching_mismatch(context, mismatch, split, N, requested_tol)
    mismatch <= _MATCHING_WARNING_TOL && return nothing
    @warn "ISEM adaptive matching mismatch above tolerance" context=context mismatch=mismatch split=split N=N warning_tol=_MATCHING_WARNING_TOL requested_tol=requested_tol
    return nothing
end

function _components_with_N_metadata(comp, N_used, N_candidates)
    return merge(comp, (N = N_used, N_candidates = N_candidates))
end

function _normalize_solution(eval_P, x, norm)
    P, dP, d2P, error = eval_P(x)
    return (P / norm, dP / norm, d2P / norm, error)
end

# P(x) from lower_eval for x < split and from upper_eval for split <= x < 0, divided by `norm`.
function _make_piecewise_solution(split, lower_eval, upper_eval, norm)
    function _P(x)
        x_selector = real(x)
        if x_selector < split
            return _normalize_solution(lower_eval, x, norm)
        elseif split <= x_selector < 0.0
            return _normalize_solution(upper_eval, x, norm)
        end
    end
    return _P
end

# Moves the junction within [split_min, split_max] by interval halving to reduce the relative
# mismatch |P_lower - P_upper|/min(|P_lower|, |P_upper|); returns the solution, junction and mismatch.
function _adaptive_piecewise_solution(split, lower_eval, upper_eval, norm, split_min, split_max; mismatch_tol = 1e-12, maxiter = 16, tol = _TOLERANCE)
    lo = min(split_min, split_max)
    hi = max(split_min, split_max)
    x0 = clamp(split, lo, hi)

    mismatch_at = x -> begin
        P1, _, _, _ = lower_eval(x)
        P2, _, _, _ = upper_eval(x)
        d0 = min(abs(P1), abs(P2))
        return abs(P1 - P2) / d0
    end

    best_x = x0
    best_m = mismatch_at(x0)
    if best_m <= mismatch_tol
        return _make_piecewise_solution(best_x, lower_eval, upper_eval, norm), best_x, best_m
    end

    for _ in 1:maxiter
        left_x = (lo + x0) / 2
        right_x = (x0 + hi) / 2
        left_m = mismatch_at(left_x)
        right_m = mismatch_at(right_x)
        if left_m <= best_m
            best_x = left_x
            best_m = left_m
        end
        if right_m <= best_m
            best_x = right_x
            best_m = right_m
        end
        if left_m <= right_m
            hi = x0
            x0 = left_x
        else
            lo = x0
            x0 = right_x
        end
        if best_m <= mismatch_tol || abs(hi - lo) <= tol
            break
        end
    end

    return _make_piecewise_solution(best_x, lower_eval, upper_eval, norm), best_x, best_m
end

# Junction in [xm_min, xm_max] for the UP and complex-frequency matches: `nothing`
# takes _xm when it lies in the range and the midpoint otherwise; a requested xm
# outside the range is rejected.
function _junction(xm, xm_min, xm_max)
    xm === nothing && return xm_min <= _xm <= xm_max ? _xm : (xm_min + xm_max) / 2
    xm_min <= xm <= xm_max || throw(ArgumentError(
        "matching needs xm in [$xm_min, $xm_max] (got xm = $xm); omit xm for the automatic junction."))
    return xm
end

# Complex-frequency IN solution: the amplitudes come from matching on the rotated contour
# x = -rho/z (z = omega/|omega|) at rhom; the returned P on the real x axis is matched separately.
function _Pin_contour(s, l, m, a, omega, rhom = _rhom, xm = nothing, N = _TruncatioN, tol = _TOLERANCE; components = false)
    params = rho_parameters(TeukolskyParameters(s, l, m, a, omega))
    z = params.z
    kappa = params.kappa
    epsilon = params.epsilon
    tau = params.tau
    lambda = params.lambda
    Btrans = B_trans(s, epsilon, tau, kappa)

    function build_zero_in_rho(target_rho)
        rho_list, _, evaluator = iterate_zero_in_rho(target_rho, s, epsilon, tau, kappa, lambda, z, N, tol)
        return rho_list, evaluator
    end

    rho_list_zero_in_rho, P_zero_in_rho = build_zero_in_rho(rhom)
    # All bases are matched at rhom_effective = max(rhom, horizon seed radius); a
    # horizon series that stopped short of it is rebuilt up to it.
    rhom_effective = max(Float64(rhom), rho_list_zero_in_rho[1])
    if rho_list_zero_in_rho[end] < rhom_effective
        rho_list_zero_in_rho, P_zero_in_rho = build_zero_in_rho(rhom_effective)
    end

    _, _, P_inf_in_rho = iterate_inf_in_rho(rhom_effective, s, epsilon, tau, kappa, lambda, z, N, tol)

    _, _, P_inf_out_rho = iterate_inf_out_rho(rhom_effective, s, epsilon, tau, kappa, lambda, z, N, tol)

    CInfIn = - C_inf(rhom_effective, P_zero_in_rho, P_inf_out_rho, s, epsilon, tau, kappa, z)
    CInfOut = C_inf(rhom_effective, P_zero_in_rho, P_inf_in_rho, s, epsilon, tau, kappa, z)
    Binc = B_inc(CInfIn, epsilon, tau, kappa)
    Bref = B_ref(CInfOut, s, epsilon, tau, kappa)

    coeffs_zero = Iteration.an_zero(s, epsilon, tau, kappa, lambda, N)
    xm_max = -Iteration.convergence_radius_an_zero(coeffs_zero, tol)

    x_list_inf_in, coe_list_inf_in, P_inf_in = iterate_inf_in(xm_max, s, epsilon, tau, kappa, lambda, N, tol)

    x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, tol)

    xm_min = min(x_list_inf_in[1], x_list_inf_out[1])

    x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol)

    lower_eval(x) = begin
        Pin, dPin, d2Pin, error_in = P_inf_in(x)
        Pout, dPout, d2Pout, error_out = P_inf_out(x)
        P = Pin * CInfIn + Pout * CInfOut
        dP = dPin * CInfIn + dPout * CInfOut
        d2P = d2Pin * CInfIn + d2Pout * CInfOut
        return (P, dP, d2P, max(error_in, error_out))
    end
    upper_eval(x) = P_zero_in(x)
    xm_match = _junction(xm, xm_min, xm_max)
    _P, xsplit, split_mismatch = _adaptive_piecewise_solution(xm_match, lower_eval, upper_eval, Btrans, xm_min, xm_max)
    if components
        return _P, (Binc / Btrans, Bref / Btrans), (
            P_zero_in = P_zero_in,
            P_inf_in = P_inf_in,
            P_inf_out = P_inf_out,
            P_zero_in_rho = P_zero_in_rho,
            P_inf_in_rho = P_inf_in_rho,
            P_inf_out_rho = P_inf_out_rho,
            CInfIn = CInfIn,
            CInfOut = CInfOut,
            Btrans = Btrans,
            xm = xm,
            xmin = xm_min,
            xmax = xm_max,
            xm_match = xm_match,
            xsplit = xsplit,
            split_mismatch = split_mismatch,
            rhom = rhom_effective,
            horizon_in_patch_count = length(x_list_zero_in),
            infinity_in_patch_count = length(x_list_inf_in),
            infinity_out_patch_count = length(x_list_inf_out),
        )
    end

    return _P, (Binc / Btrans, Bref / Btrans)
end

# Complex-frequency UP solution, structured as _Pin_contour.
function _Pup_contour(s, l, m, a, omega, rhom = _rhom, xm = nothing, N = _TruncatioN, tol = _TOLERANCE; components = false)
    params = rho_parameters(TeukolskyParameters(s, l, m, a, omega))
    z = params.z
    kappa = params.kappa
    epsilon = params.epsilon
    tau = params.tau
    lambda = params.lambda
    Ctrans = _f4(s, epsilon, kappa)

    function build_zero_in_rho(target_rho)
        rho_list, _, evaluator = iterate_zero_in_rho(target_rho, s, epsilon, tau, kappa, lambda, z, N, tol)
        return rho_list, evaluator
    end

    function build_zero_out_rho(target_rho)
        rho_list, _, evaluator = iterate_zero_out_rho(target_rho, s, epsilon, tau, kappa, lambda, z, N, tol)
        return rho_list, evaluator
    end

    rho_list_zero_in_rho, P_zero_in_rho = build_zero_in_rho(rhom)
    rho_list_zero_out_rho, P_zero_out_rho = build_zero_out_rho(rhom)
    # All bases are matched at rhom_effective = max(rhom, horizon seed radii); a
    # horizon series that stopped short of it is rebuilt up to it.
    rhom_effective = max(Float64(rhom), rho_list_zero_in_rho[1], rho_list_zero_out_rho[1])
    if rho_list_zero_in_rho[end] < rhom_effective
        rho_list_zero_in_rho, P_zero_in_rho = build_zero_in_rho(rhom_effective)
    end
    if rho_list_zero_out_rho[end] < rhom_effective
        rho_list_zero_out_rho, P_zero_out_rho = build_zero_out_rho(rhom_effective)
    end

    _, _, P_inf_out_rho = iterate_inf_out_rho(rhom_effective, s, epsilon, tau, kappa, lambda, z, N, tol)

    CZeroIn = - C_zero(rhom_effective, P_inf_out_rho, P_zero_out_rho, s, epsilon, tau, kappa, z)
    CZeroOut = C_zero(rhom_effective, P_inf_out_rho, P_zero_in_rho, s, epsilon, tau, kappa, z)
    Cinc = C_inc(CZeroOut, epsilon, tau, kappa)
    Cref = C_ref(CZeroIn, s, epsilon, tau, kappa)

    coeffs_zero_in = Iteration.an_zero(s, epsilon, tau, kappa, lambda, N)
    xm_a = -Iteration.convergence_radius_an_zero(coeffs_zero_in, tol)
    coeffs_zero_out = Iteration.bn_zero(s, epsilon, tau, kappa, lambda, N)
    xm_b = -Iteration.convergence_radius_bn_zero(coeffs_zero_out, tol)
    xm_max = max(xm_a, xm_b)

    x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, tol)

    xm_min = x_list_inf_out[1]

    x_list_zero_in, _, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol)

    x_list_zero_out, _, P_zero_out = iterate_zero_out(xm_min, s, epsilon, tau, kappa, lambda, N, tol)

    lower_eval(x) = P_inf_out(x)
    upper_eval(x) = begin
        Pin, dPin, d2Pin, error_in = P_zero_in(x)
        Pout, dPout, d2Pout, error_out = P_zero_out(x)
        P = Pin * CZeroIn + Pout * CZeroOut
        dP = dPin * CZeroIn + dPout * CZeroOut
        d2P = d2Pin * CZeroIn + d2Pout * CZeroOut
        return (P, dP, d2P, max(error_in, error_out))
    end
    xm_match = _junction(xm, xm_min, xm_max)
    _P, xsplit, split_mismatch = _adaptive_piecewise_solution(xm_match, lower_eval, upper_eval, Ctrans, xm_min, xm_max)
    if components
        return _P, (Cinc / Ctrans, Cref / Ctrans), (
            P_zero_in = P_zero_in,
            P_zero_out = P_zero_out,
            P_inf_out = P_inf_out,
            P_zero_in_rho = P_zero_in_rho,
            P_zero_out_rho = P_zero_out_rho,
            P_inf_out_rho = P_inf_out_rho,
            CZeroIn = CZeroIn,
            CZeroOut = CZeroOut,
            Ctrans = Ctrans,
            xm = xm,
            xm_min = xm_min,
            xm_max = xm_max,
            xm_match = xm_match,
            xsplit = xsplit,
            split_mismatch = split_mismatch,
            rhom = rhom_effective,
            horizon_in_patch_count = length(x_list_zero_in),
            horizon_out_patch_count = length(x_list_zero_out),
            infinity_out_patch_count = length(x_list_inf_out),
        )
    end

    return _P, (Cinc / Ctrans, Cref / Ctrans)
end

# IN solution P and (B_inc/B_trans, B_ref/B_trans) for complex omega; `components = true` also
# returns the matching metadata, and then a fixed N warns if the matching mismatch is above tolerance.
function _Pin(s, l, m, a, omega; xm = nothing, rhom = _rhom, N = nothing, tol = _TOLERANCE, components = false)
    solve_at_N = Ncand -> _Pin_contour(s, l, m, a, omega, rhom, xm, Ncand === nothing ? _TruncatioN : Ncand, tol; components = true)
    if N === nothing
        result = solve_at_N(nothing)
        if components
            return result[1], result[2], _components_with_N_metadata(result[3], nothing, [nothing])
        end
        return result[1], result[2]
    end
    Nfixed = max(Int(N), _MIN_MATCHING_N)
    result = solve_at_N(Nfixed)
    if components
        comp = _components_with_N_metadata(result[3], Nfixed, [Nfixed])
        _warn_matching_mismatch((branch = "Pin", s = s, l = l, m = m, a = a, omega = omega), comp.split_mismatch, comp.xsplit, Nfixed, _MATCHING_ACCEPT_TOL)
        return result[1], result[2], comp
    end
    return result[1], result[2]
end

# UP analogue of _Pin, returning (C_inc/C_trans, C_ref/C_trans).
function _Pup(s, l, m, a, omega; xm = nothing, rhom = _rhom, N = nothing, tol = _TOLERANCE, components = false)
    solve_at_N = Ncand -> _Pup_contour(s, l, m, a, omega, rhom, xm, Ncand === nothing ? _TruncatioN : Ncand, tol; components = true)
    if N === nothing
        result = solve_at_N(nothing)
        if components
            return result[1], result[2], _components_with_N_metadata(result[3], nothing, [nothing])
        end
        return result[1], result[2]
    end
    Nfixed = max(Int(N), _MIN_MATCHING_N)
    result = solve_at_N(Nfixed)
    if components
        comp = _components_with_N_metadata(result[3], Nfixed, [Nfixed])
        _warn_matching_mismatch((branch = "Pup", s = s, l = l, m = m, a = a, omega = omega), comp.split_mismatch, comp.xsplit, Nfixed, _MATCHING_ACCEPT_TOL)
        return result[1], result[2], comp
    end
    return result[1], result[2]
end
end
