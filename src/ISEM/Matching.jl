module Matching

include("Parameters.jl")
using .Parameters
include("TeukolskyTransformation.jl")
using .TeukolskyTransformation
include("Iteration.jl")
using .Iteration
include("AsymptoticAmplitudes.jl")
using .AsymptoticAmplitudes

export _Pin, _Pup, _P

_TruncatioN = 40
_TOLERANCE = 1e-13
_xm = - 1.0
_rhom = 1.0
_MATCHING_ACCEPT_TOL = 1e-12
_MATCHING_WARNING_TOL = 1e-8
_MAX_N_MATCHING_ATTEMPTS = 20
_MIN_MATCHING_N = 10

@inline _is_omega_complex(omega) = !iszero(imag(complex(omega)))

@inline function _warn_matching_mismatch(context, mismatch, split, N, requested_tol)
    mismatch <= _MATCHING_WARNING_TOL && return nothing
    @warn "ISEM adaptive matching mismatch above tolerance" context=context mismatch=mismatch split=split N=N warning_tol=_MATCHING_WARNING_TOL requested_tol=requested_tol
    return nothing
end

function _matching_N_candidates(N; max_attempts = _MAX_N_MATCHING_ATTEMPTS)
    N0 = max(Int(N), _MIN_MATCHING_N)
    raw_candidates = (
        N0,
        N0 + 1,
        max(N0 - 1, 10),
        N0 + 2,
        max(N0 - 2, 10),
        N0 + 3,
        max(N0 - 3, 10),
        N0 + 5,
        max(N0 - 5, 10),
        N0 + 8,
        max(N0 - 8, 10),
        N0 + 10,
        max(N0 - 10, 10),
        N0 + 15,
        max(N0 - 15, 10),
        20,
        14,
        13,
        11,
        10,
    )
    candidates = Int[]
    for Ncand0 in raw_candidates
        Ncand = max(Ncand0, _MIN_MATCHING_N)
        if Ncand >= _MIN_MATCHING_N && !(Ncand in candidates)
            push!(candidates, Ncand)
        end
        length(candidates) >= max_attempts && break
    end
    return candidates
end

@inline function _components_with_N_metadata(comp, N_used, N_candidates)
    return merge(comp, (N = N_used, N_candidates = N_candidates))
end

function _adaptive_N_solution(context, solve_at_N, N, components; mismatch_tol = _MATCHING_ACCEPT_TOL, max_attempts = _MAX_N_MATCHING_ATTEMPTS)
    candidates = _matching_N_candidates(N; max_attempts = max_attempts)
    best_result = nothing
    best_components = nothing
    best_mismatch = Inf
    best_N = first(candidates)

    for Ncand in candidates
        result = solve_at_N(Ncand)
        comp = result[3]
        mismatch = hasproperty(comp, :split_mismatch) ? comp.split_mismatch : 0.0
        if mismatch < best_mismatch
            best_result = result
            best_components = comp
            best_mismatch = mismatch
            best_N = Ncand
        end
        mismatch <= mismatch_tol && break
    end

    best_components = _components_with_N_metadata(best_components, best_N, candidates)
    _warn_matching_mismatch(context, best_mismatch, best_components.xsplit, best_N, mismatch_tol)

    if components
        return best_result[1], best_result[2], best_components
    end
    return best_result[1], best_result[2]
end

@inline function _normalize_solution(eval_P, x, norm)
    P, dP, d2P, error = eval_P(x)
    return (P / norm, dP / norm, d2P / norm, error)
end

@inline function _make_piecewise_solution(split, lower_eval, upper_eval, norm)
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

@inline _is_horizon_superradiance_threshold(s, epsilon, tau) = abs(epsilon + tau) < 1e-12

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

function _iterative_adaptive_refit_pin(xm_match, xm_min, xm_max, P_inf_in, P_inf_out, P_zero_in, (Btrans, s, epsilon, tau, kappa); mismatch_tol = 1e-12, maxiter = 5, N = nothing)
    CInfIn = 0im
    CInfOut = 0im
    xsplit = xm_match
    split_mismatch = Inf
    _P = nothing
    for _ in 1:maxiter
        CInfIn = - C_inf(xm_match, P_zero_in, P_inf_out, s, epsilon, tau, kappa)
        CInfOut = C_inf(xm_match, P_zero_in, P_inf_in, s, epsilon, tau, kappa)
        lower_eval(x) = begin
            Pin, dPin, d2Pin, error_in = P_inf_in(x)
            Pout, dPout, d2Pout, error_out = P_inf_out(x)
            P = Pin * CInfIn + Pout * CInfOut
            dP = dPin * CInfIn + dPout * CInfOut
            d2P = d2Pin * CInfIn + d2Pout * CInfOut
            return (P, dP, d2P, max(error_in, error_out))
        end
        upper_eval(x) = P_zero_in(x)
        _P, xsplit, split_mismatch = _adaptive_piecewise_solution(xm_match, lower_eval, upper_eval, Btrans, xm_min, xm_max)
        if split_mismatch <= mismatch_tol
            break
        end
        if xsplit < xm_match
            xm_max = xm_match
            xm_match = xsplit
        else
            xm_min = xm_match
            xm_match = xsplit
        end
    end
    return _P, (CInfIn, CInfOut), xsplit, split_mismatch
end

function _iterative_adaptive_refit_pup(xm_match, xm_min, xm_max, P_inf_out, P_zero_in, P_zero_out, (Ctrans, s, epsilon, tau, kappa); mismatch_tol = 1e-12, maxiter = 5, N = nothing)
    CZeroIn = 0im
    CZeroOut = 0im
    xsplit = xm_match
    split_mismatch = Inf
    _P = nothing
    for _ in 1:maxiter
        CZeroIn = - C_zero(xm_match, P_inf_out, P_zero_out, s, epsilon, tau, kappa)
        CZeroOut = C_zero(xm_match, P_inf_out, P_zero_in, s, epsilon, tau, kappa)
        lower_eval(x) = P_inf_out(x)
        upper_eval(x) = begin
            Pin, dPin, d2Pin, error_in = P_zero_in(x)
            Pout, dPout, d2Pout, error_out = P_zero_out(x)
            P = Pin * CZeroIn + Pout * CZeroOut
            dP = dPin * CZeroIn + dPout * CZeroOut
            d2P = d2Pin * CZeroIn + d2Pout * CZeroOut
            return (P, dP, d2P, max(error_in, error_out))
        end
        _P, xsplit, split_mismatch = _adaptive_piecewise_solution(xm_match, lower_eval, upper_eval, Ctrans, xm_min, xm_max)
        if split_mismatch <= mismatch_tol
            break
        end
        if xsplit < xm_match
            xm_max = xm_match
            xm_match = xsplit
        else
            xm_min = xm_match
            xm_match = xsplit
        end
    end
    return _P, (CZeroIn, CZeroOut), xsplit, split_mismatch
end

function _Pin_real(s, l, m, a, omega, xm = _xm, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0, TSinInf = 0, TSoutInf = 0, TSinHor = 0; components = false)
    params = isem_parameters(TeukolskyParameters(s, l, m, a, omega))
    kappa = params.kappa
    epsilon = params.epsilon
    tau = params.tau
    lambda = params.lambda
    Btrans = B_trans(s, epsilon, tau, kappa)

    if TSinHor == 1
        if s == -2
            x_list_zero_in, coe_list_zero_in, P_zero_in_plus = iterate_zero_in(xm, -s, epsilon, tau, kappa, lambda - 4, N, tol, lfe)
            Rplus, P_zero_in = PlusTwo_to_MinusTwo_Hor(P_zero_in_plus, m, a, omega, lambda - 4)
        else
            x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    else
        x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    if TSinInf == 1
        if s == 2
            x_list_inf_in, coe_list_inf_in, P_inf_in_minus = iterate_inf_in(xm, -s, epsilon, tau, kappa, lambda + 4, N, l, tol, sfe)
            Rminus, P_inf_in = MinusTwo_to_PlusTwo_Inf(P_inf_in_minus, m, a, omega, lambda + 4)
        else
            x_list_inf_in, coe_list_inf_in, P_inf_in = iterate_inf_in(xm, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_in, coe_list_inf_in, P_inf_in = iterate_inf_in(xm, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    if TSoutInf == 1
        if s == -2
            x_list_inf_out, coe_list_inf_out, P_inf_out_plus = iterate_inf_out(xm, -s, epsilon, tau, kappa, lambda - 4, N, l, tol, sfe)
            Rplus, P_inf_out = PlusTwo_to_MinusTwo_Inf(P_inf_out_plus, m, a, omega, lambda - 4)
        else
            x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    CInfIn = - C_inf(xm, P_zero_in, P_inf_out, s, epsilon, tau, kappa)
    CInfOut = C_inf(xm, P_zero_in, P_inf_in, s, epsilon, tau, kappa)
    Binc = B_inc(CInfIn, epsilon, tau, kappa)
    Bref = B_ref(CInfOut, s, epsilon, tau, kappa)

    lower_eval(x) = begin
        Pin, dPin, d2Pin, error_in = P_inf_in(x)
        Pout, dPout, d2Pout, error_out = P_inf_out(x)
        P = Pin * CInfIn + Pout * CInfOut
        dP = dPin * CInfIn + dPout * CInfOut
        d2P = d2Pin * CInfIn + d2Pout * CInfOut
        return (P, dP, d2P, max(error_in, error_out))
    end
    upper_eval(x) = P_zero_in(x)
    _P = _make_piecewise_solution(xm, lower_eval, upper_eval, Btrans)

    if components
        return _P, (Binc / Btrans, Bref / Btrans), (
            P_zero_in = P_zero_in,
            P_inf_in = P_inf_in,
            P_inf_out = P_inf_out,
            CInfIn = CInfIn,
            CInfOut = CInfOut,
            Btrans = Btrans,
            xm = xm,
            xsplit = xm,
        )
    end

    return _P, (Binc / Btrans, Bref / Btrans)
end

function _Pin_real_adaptive(s, l, m, a, omega, xm = _xm, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0, TSinInf = 0, TSoutInf = 0, TSinHor = 0; components = false)
    params = isem_parameters(TeukolskyParameters(s, l, m, a, omega))
    kappa = params.kappa
    epsilon = params.epsilon
    tau = params.tau
    lambda = params.lambda
    Btrans = B_trans(s, epsilon, tau, kappa)

    coeffs_zero = N === nothing ?
        Iteration.LocalSolutionAtZero.an_zero_auto(s, epsilon, tau, kappa, lambda, Iteration._auto_nmax(lambda, epsilon, kappa); tol = tol)[1] :
        Iteration.an_zero(s, epsilon, tau, kappa, lambda, N)
    xm_max = -Iteration.convergence_radius_an_zero(coeffs_zero, tol)

    if TSinInf == 1
        if s == 2
            x_list_inf_in, coe_list_inf_in, P_inf_in_minus = iterate_inf_in(xm_max, -s, epsilon, tau, kappa, lambda + 4, N, l, tol, sfe)
            Rminus, P_inf_in = MinusTwo_to_PlusTwo_Inf(P_inf_in_minus, m, a, omega, lambda + 4)
        else
            x_list_inf_in, coe_list_inf_in, P_inf_in = iterate_inf_in(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_in, coe_list_inf_in, P_inf_in = iterate_inf_in(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    if TSoutInf == 1
        if s == -2
            x_list_inf_out, coe_list_inf_out, P_inf_out_plus = iterate_inf_out(xm_max, -s, epsilon, tau, kappa, lambda - 4, N, l, tol, sfe)
            Rplus, P_inf_out = PlusTwo_to_MinusTwo_Inf(P_inf_out_plus, m, a, omega, lambda - 4)
        else
            x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    xm_min = min(x_list_inf_in[1], x_list_inf_out[1])

    if TSinHor == 1
        if s == -2
            x_list_zero_in, coe_list_zero_in, P_zero_in_plus = iterate_zero_in(xm_min, -s, epsilon, tau, kappa, lambda - 4, N, tol, lfe)
            Rplus, P_zero_in = PlusTwo_to_MinusTwo_Hor(P_zero_in_plus, m, a, omega, lambda - 4)
        else
            x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    else
        x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    xm_match = (xm < xm_min || xm > xm_max) ? (xm_min + xm_max) / 2 : xm
    _P, (CInfIn, CInfOut), xsplit, split_mismatch = _iterative_adaptive_refit_pin(xm_match, xm_min, xm_max, P_inf_in, P_inf_out, P_zero_in, (Btrans, s, epsilon, tau, kappa); mismatch_tol = 1e-12, maxiter = 5, N = N)
    Binc = B_inc(CInfIn, epsilon, tau, kappa)
    Bref = B_ref(CInfOut, s, epsilon, tau, kappa)

    if components
        return _P, (Binc / Btrans, Bref / Btrans), (
            P_zero_in = P_zero_in,
            P_inf_in = P_inf_in,
            P_inf_out = P_inf_out,
            CInfIn = CInfIn,
            CInfOut = CInfOut,
            Btrans = Btrans,
            xm = xm,
            xm_min = xm_min,
            xm_max = xm_max,
            xm_match = xm_match,
            xsplit = xsplit,
            split_mismatch = split_mismatch,
        )
    end

    return _P, (Binc / Btrans, Bref / Btrans)
end

function _Pup_real(s, l, m, a, omega, xm = _xm, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0, TSoutInf = 0, TSinHor = 0, TSoutHor = 0; components = false)
    params = isem_parameters(TeukolskyParameters(s, l, m, a, omega))
    kappa = params.kappa
    epsilon = params.epsilon
    tau = params.tau
    lambda = params.lambda
    Ctrans = C_trans(s, epsilon, kappa)

    if TSinHor == 1
        if s == -2
            x_list_zero_in, coe_list_zero_in, P_zero_in_plus = iterate_zero_in(xm, -s, epsilon, tau, kappa, lambda - 4, N, tol, lfe)
            Rplus, P_zero_in = PlusTwo_to_MinusTwo_Hor(P_zero_in_plus, m, a, omega, lambda - 4)
        else
            x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    else
        x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    if TSoutHor == 1
        if s == 2
            x_list_zero_out, coe_list_zero_out, P_zero_out_minus = iterate_zero_out(xm, -s, epsilon, tau, kappa, lambda + 4, N, tol, lfe)
            Rminus, P_zero_out = MinusTwo_to_PlusTwo_Hor(P_zero_out_minus, m, a, omega, lambda + 4)
        else
            x_list_zero_out, coe_list_zero_out, P_zero_out = iterate_zero_out(xm, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    else
        x_list_zero_out, coe_list_zero_out, P_zero_out = iterate_zero_out(xm, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    if TSoutInf == 1
        if s == -2
            x_list_inf_out, coe_list_inf_out, P_inf_out_plus = iterate_inf_out(xm, -s, epsilon, tau, kappa, lambda - 4, N, l, tol, sfe)
            Rplus, P_inf_out = PlusTwo_to_MinusTwo_Inf(P_inf_out_plus, m, a, omega, lambda - 4)
        else
            x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    CZeroIn = - C_zero(xm, P_inf_out, P_zero_out, s, epsilon, tau, kappa)
    CZeroOut = C_zero(xm, P_inf_out, P_zero_in, s, epsilon, tau, kappa)
    Cinc = C_inc(CZeroOut, epsilon, tau, kappa)
    Cref = C_ref(CZeroIn, s, epsilon, tau, kappa)

    lower_eval(x) = P_inf_out(x)
    upper_eval(x) = begin
        Pin, dPin, d2Pin, error_in = P_zero_in(x)
        Pout, dPout, d2Pout, error_out = P_zero_out(x)
        P = Pin * CZeroIn + Pout * CZeroOut
        dP = dPin * CZeroIn + dPout * CZeroOut
        d2P = d2Pin * CZeroIn + d2Pout * CZeroOut
        return (P, dP, d2P, max(error_in, error_out))
    end
    _P = _make_piecewise_solution(xm, lower_eval, upper_eval, Ctrans)

    if components
        return _P, (Cinc / Ctrans, Cref / Ctrans), (
            P_zero_in = P_zero_in,
            P_zero_out = P_zero_out,
            P_inf_out = P_inf_out,
            CZeroIn = CZeroIn,
            CZeroOut = CZeroOut,
            Ctrans = Ctrans,
            xm = xm,
            xsplit = xm,
        )
    end

    return _P, (Cinc / Ctrans, Cref / Ctrans)
end

function _Pup_real_adaptive(s, l, m, a, omega, xm = _xm, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0, TSoutInf = 0, TSinHor = 0, TSoutHor = 0; components = false)
    params = isem_parameters(TeukolskyParameters(s, l, m, a, omega))
    kappa = params.kappa
    epsilon = params.epsilon
    tau = params.tau
    lambda = params.lambda
    Ctrans = C_trans(s, epsilon, kappa)
    use_resonant_zero_in = s > 0 && _is_horizon_superradiance_threshold(s, epsilon, tau) && TSinHor == 0
    use_resonant_zero_out = s < 0 && _is_horizon_superradiance_threshold(s, epsilon, tau) && TSoutHor == 0

    xm_a = if use_resonant_zero_in
        -0.25
    else
        coeffs_zero_in = N === nothing ?
            Iteration.LocalSolutionAtZero.an_zero_auto(s, epsilon, tau, kappa, lambda, Iteration._auto_nmax(lambda, epsilon, kappa); tol = tol)[1] :
            Iteration.an_zero(s, epsilon, tau, kappa, lambda, N)
        -Iteration.convergence_radius_an_zero(coeffs_zero_in, tol)
    end
    xm_b = if use_resonant_zero_out
        -0.25
    else
        coeffs_zero_out = N === nothing ?
            Iteration.LocalSolutionAtZero.bn_zero_auto(s, epsilon, tau, kappa, lambda, Iteration._auto_nmax(lambda, epsilon, kappa); tol = tol)[1] :
            Iteration.bn_zero(s, epsilon, tau, kappa, lambda, N)
        -Iteration.convergence_radius_bn_zero(coeffs_zero_out, tol)
    end
    xm_max = max(xm_a, xm_b)

    if TSoutInf == 1
        if s == -2
            x_list_inf_out, coe_list_inf_out, P_inf_out_plus = iterate_inf_out(xm_max, -s, epsilon, tau, kappa, lambda - 4, N, l, tol, sfe)
            Rplus, P_inf_out = PlusTwo_to_MinusTwo_Inf(P_inf_out_plus, m, a, omega, lambda - 4)
        else
            x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    xm_min = x_list_inf_out[1]

    if TSinHor == 1
        if s == -2
            x_list_zero_in, coe_list_zero_in, P_zero_in_plus = iterate_zero_in(xm_min, -s, epsilon, tau, kappa, lambda - 4, N, tol, lfe)
            Rplus, P_zero_in = PlusTwo_to_MinusTwo_Hor(P_zero_in_plus, m, a, omega, lambda - 4)
        else
            x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    elseif use_resonant_zero_in
        x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in_resonant_log(xm_min, s, epsilon, kappa, lambda, N, tol)
    else
        x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    if TSoutHor == 1
        if s == 2
            x_list_zero_out, coe_list_zero_out, P_zero_out_minus = iterate_zero_out(xm_min, -s, epsilon, tau, kappa, lambda + 4, N, tol, lfe)
            Rminus, P_zero_out = MinusTwo_to_PlusTwo_Hor(P_zero_out_minus, m, a, omega, lambda + 4)
        else
            x_list_zero_out, coe_list_zero_out, P_zero_out = iterate_zero_out(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    elseif use_resonant_zero_out
        x_list_zero_out, coe_list_zero_out, P_zero_out = iterate_zero_out_resonant_log(xm_min, s, epsilon, kappa, lambda, N, tol)
    else
        x_list_zero_out, coe_list_zero_out, P_zero_out = iterate_zero_out(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    xm_match = (xm < xm_min || xm > xm_max) ? (xm_min + xm_max) / 2 : xm
    _P, (CZeroIn, CZeroOut), xsplit, split_mismatch = _iterative_adaptive_refit_pup(xm_match, xm_min, xm_max, P_inf_out, P_zero_in, P_zero_out, (Ctrans, s, epsilon, tau, kappa); mismatch_tol = 1e-12, maxiter = 5, N = N)
    Cinc = C_inc(CZeroOut, epsilon, tau, kappa)
    Cref = C_ref(CZeroIn, s, epsilon, tau, kappa)

    if components
        return _P, (Cinc / Ctrans, Cref / Ctrans), (
            P_zero_in = P_zero_in,
            P_zero_out = P_zero_out,
            P_inf_out = P_inf_out,
            CZeroIn = CZeroIn,
            CZeroOut = CZeroOut,
            Ctrans = Ctrans,
            xm = xm,
            xm_min = xm_min,
            xm_max = xm_max,
            xm_match = xm_match,
            xsplit = xsplit,
            split_mismatch = split_mismatch,
        )
    end

    return _P, (Cinc / Ctrans, Cref / Ctrans)
end

function _Pin_contour(s, l, m, a, omega, rhom = _rhom, xm = _xm, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0, TSinInf = 0, TSoutInf = 0, TSinHor = 0; components = false)
    params = rho_parameters(TeukolskyParameters(s, l, m, a, omega))
    z = params.z
    kappa = params.kappa
    epsilon = params.epsilon
    tau = params.tau
    lambda = params.lambda
    Btrans = B_trans(s, epsilon, tau, kappa)

    if TSinHor == 1
        if s == -2
            _, _, P_zero_in_plus_rho = iterate_zero_in_rho(rhom, -s, epsilon, tau, kappa, lambda - 4, z, N, tol)
            _, P_zero_in_rho = contour_PlusTwo_to_MinusTwo_Hor(P_zero_in_plus_rho, m, a, omega, lambda - 4, z)
        else
            _, _, P_zero_in_rho = iterate_zero_in_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
        end
    else
        _, _, P_zero_in_rho = iterate_zero_in_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
    end

    if TSinInf == 1
        if s == 2
            _, _, P_inf_in_minus_rho = iterate_inf_in_rho(rhom, -s, epsilon, tau, kappa, lambda + 4, z, N, tol)
            _, P_inf_in_rho = contour_MinusTwo_to_PlusTwo_Inf(P_inf_in_minus_rho, m, a, omega, lambda + 4, z)
        else
            _, _, P_inf_in_rho = iterate_inf_in_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
        end
    else
        _, _, P_inf_in_rho = iterate_inf_in_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
    end

    if TSoutInf == 1
        if s == -2
            _, _, P_inf_out_plus_rho = iterate_inf_out_rho(rhom, -s, epsilon, tau, kappa, lambda - 4, z, N, tol)
            _, P_inf_out_rho = contour_PlusTwo_to_MinusTwo_Inf(P_inf_out_plus_rho, m, a, omega, lambda - 4, z)
        else
            _, _, P_inf_out_rho = iterate_inf_out_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
        end
    else
        _, _, P_inf_out_rho = iterate_inf_out_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
    end

    CInfIn = - C_inf(rhom, P_zero_in_rho, P_inf_out_rho, s, epsilon, tau, kappa, z)
    CInfOut = C_inf(rhom, P_zero_in_rho, P_inf_in_rho, s, epsilon, tau, kappa, z)
    Binc = B_inc(CInfIn, epsilon, tau, kappa)
    Bref = B_ref(CInfOut, s, epsilon, tau, kappa)

    coeffs_zero = Iteration.an_zero(s, epsilon, tau, kappa, lambda, N)
    xm_max = -Iteration.convergence_radius_an_zero(coeffs_zero, tol)

    if TSinInf == 1
        if s == 2
            x_list_inf_in, coe_list_inf_in, P_inf_in_minus = iterate_inf_in(xm_max, -s, epsilon, tau, kappa, lambda + 4, N, l, tol, sfe)
            Rminus, P_inf_in = MinusTwo_to_PlusTwo_Inf(P_inf_in_minus, m, a, omega, lambda + 4)
        else
            x_list_inf_in, coe_list_inf_in, P_inf_in = iterate_inf_in(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_in, coe_list_inf_in, P_inf_in = iterate_inf_in(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    if TSoutInf == 1
        if s == -2
            x_list_inf_out, coe_list_inf_out, P_inf_out_plus = iterate_inf_out(xm_max, -s, epsilon, tau, kappa, lambda - 4, N, l, tol, sfe)
            Rplus, P_inf_out = PlusTwo_to_MinusTwo_Inf(P_inf_out_plus, m, a, omega, lambda - 4)
        else
            x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    xm_min = min(x_list_inf_in[1], x_list_inf_out[1])

    if TSinHor == 1
        if s == -2
            x_list_zero_in, coe_list_zero_in, P_zero_in_plus = iterate_zero_in(xm_min, -s, epsilon, tau, kappa, lambda - 4, N, tol, lfe)
            Rplus, P_zero_in = PlusTwo_to_MinusTwo_Hor(P_zero_in_plus, m, a, omega, lambda - 4)
        else
            x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    else
        x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    lower_eval(x) = begin
        Pin, dPin, d2Pin, error_in = P_inf_in(x)
        Pout, dPout, d2Pout, error_out = P_inf_out(x)
        P = Pin * CInfIn + Pout * CInfOut
        dP = dPin * CInfIn + dPout * CInfOut
        d2P = d2Pin * CInfIn + d2Pout * CInfOut
        return (P, dP, d2P, max(error_in, error_out))
    end
    upper_eval(x) = P_zero_in(x)
    xm_match = (xm < xm_min || xm > xm_max) ? (xm_min + xm_max) / 2 : xm
    _P, xsplit, split_mismatch = _adaptive_piecewise_solution(xm_match, lower_eval, upper_eval, Btrans, xm_min, xm_max)
    if components
        return _P, (Binc / Btrans, Bref / Btrans), (
            P_zero_in = P_zero_in,
            P_inf_in = P_inf_in,
            P_inf_out = P_inf_out,
            CInfIn = CInfIn,
            CInfOut = CInfOut,
            Btrans = Btrans,
            xm = xm,
            xmin = xm_min,
            xmax = xm_max,
            xm_match = xm_match,
            xsplit = xsplit,
            split_mismatch = split_mismatch,
        )
    end

    return _P, (Binc / Btrans, Bref / Btrans)
end

function _Pup_contour(s, l, m, a, omega, rhom = _rhom, xm = _xm, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0, TSoutInf = 0, TSinHor = 0, TSoutHor = 0; components = false)
    params = rho_parameters(TeukolskyParameters(s, l, m, a, omega))
    z = params.z
    kappa = params.kappa
    epsilon = params.epsilon
    tau = params.tau
    lambda = params.lambda
    Ctrans = _f4(s, epsilon, kappa)

    if TSinHor == 1
        if s == -2
            _, _, P_zero_in_plus_rho = iterate_zero_in_rho(rhom, -s, epsilon, tau, kappa, lambda - 4, z, N, tol)
            _, P_zero_in_rho = contour_PlusTwo_to_MinusTwo_Hor(P_zero_in_plus_rho, m, a, omega, lambda - 4, z)
        else
            _, _, P_zero_in_rho = iterate_zero_in_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
        end
    else
        _, _, P_zero_in_rho = iterate_zero_in_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
    end

    if TSoutHor == 1
        if s == 2
            _, _, P_zero_out_minus_rho = iterate_zero_out_rho(rhom, -s, epsilon, tau, kappa, lambda + 4, z, N, tol)
            _, P_zero_out_rho = contour_MinusTwo_to_PlusTwo_Hor(P_zero_out_minus_rho, m, a, omega, lambda + 4, z)
        else
            _, _, P_zero_out_rho = iterate_zero_out_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
        end
    else
        _, _, P_zero_out_rho = iterate_zero_out_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
    end

    if TSoutInf == 1
        if s == -2
            _, _, P_inf_out_plus_rho = iterate_inf_out_rho(rhom, -s, epsilon, tau, kappa, lambda - 4, z, N, tol)
            _, P_inf_out_rho = contour_PlusTwo_to_MinusTwo_Inf(P_inf_out_plus_rho, m, a, omega, lambda - 4, z)
        else
            _, _, P_inf_out_rho = iterate_inf_out_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
        end
    else
        _, _, P_inf_out_rho = iterate_inf_out_rho(rhom, s, epsilon, tau, kappa, lambda, z, N, tol)
    end

    CZeroIn = - C_zero(rhom, P_inf_out_rho, P_zero_out_rho, s, epsilon, tau, kappa, z)
    CZeroOut = C_zero(rhom, P_inf_out_rho, P_zero_in_rho, s, epsilon, tau, kappa, z)
    Cinc = C_inc(CZeroOut, epsilon, tau, kappa)
    Cref = C_ref(CZeroIn, s, epsilon, tau, kappa)

    coeffs_zero_in = Iteration.an_zero(s, epsilon, tau, kappa, lambda, N)
    xm_a = -Iteration.convergence_radius_an_zero(coeffs_zero_in, tol)
    coeffs_zero_out = Iteration.bn_zero(s, epsilon, tau, kappa, lambda, N)
    xm_b = -Iteration.convergence_radius_bn_zero(coeffs_zero_out, tol)
    xm_max = max(xm_a, xm_b)

    if TSoutInf == 1
        if s == -2
            x_list_inf_out, coe_list_inf_out, P_inf_out_plus = iterate_inf_out(xm_max, -s, epsilon, tau, kappa, lambda - 4, N, l, tol, sfe)
            Rplus, P_inf_out = PlusTwo_to_MinusTwo_Inf(P_inf_out_plus, m, a, omega, lambda - 4)
        else
            x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
        end
    else
        x_list_inf_out, coe_list_inf_out, P_inf_out = iterate_inf_out(xm_max, s, epsilon, tau, kappa, lambda, N, l, tol, sfe)
    end

    xm_min = x_list_inf_out[1]

    if TSinHor == 1
        if s == -2
            x_list_zero_in, coe_list_zero_in, P_zero_in_plus = iterate_zero_in(xm_min, -s, epsilon, tau, kappa, lambda - 4, N, tol, lfe)
            Rplus, P_zero_in = PlusTwo_to_MinusTwo_Hor(P_zero_in_plus, m, a, omega, lambda - 4)
        else
            x_list_zero_in, coe_list_zero_in, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    else
        x_list_zero_in, _, P_zero_in = iterate_zero_in(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    if TSoutHor == 1
        if s == 2
            x_list_zero_out, coe_list_zero_out, P_zero_out_minus = iterate_zero_out(xm_min, -s, epsilon, tau, kappa, lambda + 4, N, tol, lfe)
            Rminus, P_zero_out = MinusTwo_to_PlusTwo_Hor(P_zero_out_minus, m, a, omega, lambda + 4)
        else
            x_list_zero_out, coe_list_zero_out, P_zero_out = iterate_zero_out(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
        end
    else
        x_list_zero_out, _, P_zero_out = iterate_zero_out(xm_min, s, epsilon, tau, kappa, lambda, N, tol, lfe)
    end

    lower_eval(x) = P_inf_out(x)
    upper_eval(x) = begin
        Pin, dPin, d2Pin, error_in = P_zero_in(x)
        Pout, dPout, d2Pout, error_out = P_zero_out(x)
        P = Pin * CZeroIn + Pout * CZeroOut
        dP = dPin * CZeroIn + dPout * CZeroOut
        d2P = d2Pin * CZeroIn + d2Pout * CZeroOut
        return (P, dP, d2P, max(error_in, error_out))
    end
    xm_match = (xm < xm_min || xm > xm_max) ? (xm_min + xm_max) / 2 : xm
    _P, xsplit, split_mismatch = _adaptive_piecewise_solution(xm_match, lower_eval, upper_eval, Ctrans, xm_min, xm_max)
    if components
        return _P, (Cinc / Ctrans, Cref / Ctrans), (
            P_zero_in = P_zero_in,
            P_zero_out = P_zero_out,
            P_inf_out = P_inf_out,
            CZeroIn = CZeroIn,
            CZeroOut = CZeroOut,
            Ctrans = Ctrans,
            xm = xm,
            xm_min = xm_min,
            xm_max = xm_max,
            xm_match = xm_match,
            xsplit = xsplit,
            split_mismatch = split_mismatch,
        )
    end

    return _P, (Cinc / Ctrans, Cref / Ctrans)
end

function _Pin(s, l, m, a, omega; xm = _xm, rhom = _rhom, N = nothing, tol = _TOLERANCE, sfe = 0, lfe = 0, TSinInf = 0, TSoutInf = 0, TSinHor = 0, components = false)
    solve_at_N = if _is_omega_complex(omega)
        Ncand -> _Pin_contour(s, l, m, a, omega, rhom, xm, Ncand === nothing ? _TruncatioN : Ncand, tol, sfe, lfe, TSinInf, TSoutInf, TSinHor; components = true)
    else
        Ncand -> _Pin_real_adaptive(s, l, m, a, omega, xm, Ncand, tol, sfe, lfe, TSinInf, TSoutInf, TSinHor; components = true)
    end
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

function _Pup(s, l, m, a, omega; xm = _xm, rhom = _rhom, N = nothing, tol = _TOLERANCE, sfe = 0, lfe = 0, TSoutInf = 0, TSinHor = 0, TSoutHor = 0, components = false)
    solve_at_N = if _is_omega_complex(omega)
        Ncand -> _Pup_contour(s, l, m, a, omega, rhom, xm, Ncand === nothing ? _TruncatioN : Ncand, tol, sfe, lfe, TSoutInf, TSinHor, TSoutHor; components = true)
    else
        Ncand -> _Pup_real_adaptive(s, l, m, a, omega, xm, Ncand, tol, sfe, lfe, TSoutInf, TSinHor, TSoutHor; components = true)
    end
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

function _P(s, l, m, a, omega; xm = _xm, rhom = _rhom, N = nothing, tol = _TOLERANCE, sfe = 0, lfe = 0, TSinInf = 0, TSoutInf = 0, TSinHor = 0, TSoutHor = 0, components = false)
    Pin, (Binc, Bref) = _Pin(s, l, m, a, omega; xm = xm, rhom = rhom, N = N, tol = tol, sfe = sfe, lfe = lfe, TSinInf = TSinInf, TSoutInf = TSoutInf, TSinHor = TSinHor, components = components)
    Pup, (Cinc, Cref) = _Pup(s, l, m, a, omega; xm = xm, rhom = rhom, N = N, tol = tol, sfe = sfe, lfe = lfe, TSoutInf = TSoutInf, TSinHor = TSinHor, TSoutHor = TSoutHor, components = components)
    return Pin, Pup, (Binc, Bref, Cinc, Cref)
end

end
