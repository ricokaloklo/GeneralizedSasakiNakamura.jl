module Iteration

include("LocalSolutionAtZero.jl")
using .LocalSolutionAtZero
include("LocalSolutionAtInfinity.jl")
using .LocalSolutionAtInfinity
include("OrdinaryPointExpansion.jl")
using .OrdinaryPointExpansion
include("SmallFrequencyExpansion.jl")
using .SmallFrequencyExpansion
include("LargeFrequencyExpansion.jl")
using .LargeFrequencyExpansion

export iterate_inf_in, iterate_inf_out, iterate_zero_in, iterate_zero_out,
       iterate_inf_in_rho, iterate_inf_out_rho, iterate_zero_in_rho, iterate_zero_out_rho,
       iterate_zero_out_resonant_log, iterate_zero_in_resonant_log

_TruncatioN = 40
_TOLERANCE = 1e-13
_LFE_MIN_RADIUS = 1e-6
_LFE_MAX_STEPS = 100
_LFE_MIN_EK = 3.0
_LFE_INTERNAL_NMAX = 16
_AUTO_NMAX = 100
_AUTO_PATIENCE = 12
_AUTO_SFE_N = 30
_AUTO_STEP_SAFETY = 0.25

@inline function _fixed_or_auto_N(N)
    return N === nothing ? _AUTO_NMAX : N
end

@inline function _auto_nmax(lambda, epsilon, kappa)
    if abs(epsilon) < 0.01 && abs(kappa) < 0.15
        return 60
    end
    return _AUTO_NMAX
end

@inline function _auto_order(lambda)
    return min(_AUTO_NMAX, max(10, ceil(Int, sqrt(abs(lambda))) + 5))
end

@inline function _sfe_order(N)
    return N === nothing ? _AUTO_SFE_N : N
end

@inline function _ordinary_coeffs(s, epsilon, tau, kappa, lambda, x0, P, P_prime, N, tol)
    if N === nothing
        return ordinary_point_coeffs_auto(s, epsilon, tau, kappa, lambda, x0, P, P_prime, _auto_nmax(lambda, epsilon, kappa); tol = tol, patience = _AUTO_PATIENCE)[1]
    end
    return ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, x0, P, P_prime, N)
end

@inline function _ordinary_step_size(left_radius, right_radius, N)
    Δx = min(left_radius, right_radius)
    return N === nothing ? _AUTO_STEP_SAFETY * Δx : Δx
end

function iterate_inf_in(xm, s, epsilon, tau, kappa, lambda, N = _TruncatioN, l = nothing, tol = _TOLERANCE, sfe = 0)
    use_sfe = sfe == 1
    if use_sfe
        nu_in, coeffs_inf, Psfe_in, Rin = sfe_in(s, epsilon, tau, kappa, lambda, _sfe_order(N), l)
        x_min = - Rin
        P, P_prime, _, _ = Psfe_in(x_min)
    else
        coeffs_inf, Nopt, Rin = an_inf(s, epsilon, tau, kappa, lambda, tol)
        x_min = - convergence_radius_an_inf(Rin, s, epsilon, tau, kappa, lambda, coeffs_inf, tol)
        P, P_prime, _ = evaluate_inf_in(x_min, s, epsilon, tau, coeffs_inf)
    end
    coeffs_list = Vector{typeof(coeffs_inf)}()
    x_list = Vector{Float64}()
    push!(x_list, x_min)

    while x_min < xm
            coeffs_ord = _ordinary_coeffs(s, epsilon, tau, kappa, lambda, x_min, P, P_prime, N, tol)
        Δx = _ordinary_step_size(convergence_radius_x0_right(x_min, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), convergence_radius_x0_left(x_min, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), N)
        if !(isfinite(Δx)) || Δx <= 0
            Base.error("Step size Δx invalid: $Δx")
        end
        P, P_prime, _ = evaluate_ordinary_point_series(x_min + Δx, x_min, coeffs_ord)
        x_min += Δx
        push!(coeffs_list, coeffs_ord)
        push!(x_list, x_min)
    end

    function evaluate_P(x)
        if x <= x_list[1]
            if use_sfe
                return Psfe_in(x)
            else
                Px, Px_prime, Px_primeprime = evaluate_inf_in(x, s, epsilon, tau, coeffs_inf)
                err = error_inf_in(x, s, epsilon, tau, kappa, lambda, coeffs_inf)
                return (Px, Px_prime, Px_primeprime, err)
            end
        elseif x_list[end] < x < 0.0
            Base.error("x = $x is outside computed domain x <= $(x_list[end]).")
        elseif x >= 0.0
            Base.error("x should be a negative real number.")
        end

        i = searchsortedlast(x_list, x)
        i > length(coeffs_list) && Base.error("x = $x is outside computed ordinary-point domain.")
        x0 = x_list[i]
        coeffs = coeffs_list[i]
        Px, Px_prime, Px_primeprime = evaluate_ordinary_point_series(x, x0, coeffs)
        err = error_ordinary((x, Px, Px_prime, Px_primeprime), s, epsilon, tau, kappa, lambda)
        return (Px, Px_prime, Px_primeprime, err)
    end

    return x_list, coeffs_list, evaluate_P
end

function iterate_inf_out(xm, s, epsilon, tau, kappa, lambda, N = _TruncatioN, l = nothing, tol = _TOLERANCE, sfe = 0)
    use_sfe = sfe == 1
    if use_sfe
        nu_out, coeffs_inf, Psfe_out, Rout = sfe_out(s, epsilon, tau, kappa, lambda, _sfe_order(N), l)
        x_min = - Rout
        P, P_prime, _, _ = Psfe_out(x_min)
    else
        coeffs_inf, Nopt, Rout = bn_inf(s, epsilon, tau, kappa, lambda, tol)
        x_min = - convergence_radius_bn_inf(Rout, s, epsilon, tau, kappa, lambda, coeffs_inf, tol)
        P, P_prime, _ = evaluate_inf_out(x_min, s, epsilon, tau, kappa, coeffs_inf)
    end
    coeffs_list = Vector{typeof(coeffs_inf)}()
    x_list = Vector{Float64}()
    push!(x_list, x_min)

    while x_min < xm
            coeffs_ord = _ordinary_coeffs(s, epsilon, tau, kappa, lambda, x_min, P, P_prime, N, tol)
        Δx = _ordinary_step_size(convergence_radius_x0_right(x_min, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), convergence_radius_x0_left(x_min, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), N)
        if !(isfinite(Δx)) || Δx <= 0
            Base.error("Step size Δx invalid: $Δx")
        end
        P, P_prime, _ = evaluate_ordinary_point_series(x_min + Δx, x_min, coeffs_ord)
        x_min += Δx
        push!(coeffs_list, coeffs_ord)
        push!(x_list, x_min)
    end

    function evaluate_P(x)
        if x <= x_list[1]
            if use_sfe
                return Psfe_out(x)
            else
                Px, Px_prime, Px_primeprime = evaluate_inf_out(x, s, epsilon, tau, kappa, coeffs_inf)
                err = error_inf_out(x, s, epsilon, tau, kappa, lambda, coeffs_inf)
                return (Px, Px_prime, Px_primeprime, err)
            end
        elseif x_list[end] < x < 0.0
            Base.error("x = $x is outside computed domain x <= $(x_list[end]).")
        elseif x >= 0.0
            Base.error("x should be a negative real number.")
        end

        i = searchsortedlast(x_list, x)
        i > length(coeffs_list) && Base.error("x = $x is outside computed ordinary-point domain.")
        x0 = x_list[i]
        coeffs = coeffs_list[i]
        Px, Px_prime, Px_primeprime = evaluate_ordinary_point_series(x, x0, coeffs)
        err = error_ordinary((x, Px, Px_prime, Px_primeprime), s, epsilon, tau, kappa, lambda)
        return (Px, Px_prime, Px_primeprime, err)
    end

    return x_list, coeffs_list, evaluate_P
end

function iterate_zero_in(xm, s, epsilon, tau, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE, lfe = 0)
    if lfe == 1
        Nlfe = min(_fixed_or_auto_N(N), _LFE_INTERNAL_NMAX)
        if abs(epsilon * kappa) < _LFE_MIN_EK
            @warn("LFE zero-in skipped because |epsilon*kappa| = $(abs(epsilon * kappa)) is below the threshold $_LFE_MIN_EK.")
            return iterate_zero_in(xm, s, epsilon, tau, kappa, lambda, N, tol, 0)
        end
        Xin, Rin = lfe_in_zero(s, epsilon, tau, kappa, lambda, Nlfe)
        if !isfinite(Rin) || Rin < _LFE_MIN_RADIUS || ceil(Int, max(0.0, abs(xm) - Rin) / Rin) > _LFE_MAX_STEPS
            @warn("LFE zero-in radius $Rin is too small; falling back to the standard zero expansion.")
            return iterate_zero_in(xm, s, epsilon, tau, kappa, lambda, N, tol, 0)
        end
        x_max = - Rin
        X, X_prime, _, _ = Xin(x_max)
        func_list = Vector{Any}()
        P_list = Vector{Any}()
        x_list = Vector{Float64}()
        push!(x_list, x_max)
        push!(func_list, Xin)
        push!(P_list, X_to_P(Xin, s, epsilon, tau, kappa))

        while x_max > xm
            X0, Δx = lfe_ordinary_point(s, epsilon, tau, kappa, lambda, x_max, X, X_prime, Nlfe, tol)
            if !(isfinite(Δx)) || Δx <= 0
                Base.error("Step size Δx invalid: $Δx")
            end
            X, X_prime, _, _ = X0(x_max - Δx)
            x_max -= Δx
            push!(func_list, X0)
            push!(P_list, X_to_P(X0, s, epsilon, tau, kappa))
            push!(x_list, x_max)
        end

        function lfe_P(x)
            if x_list[1] <= x < 0.0
                return P_list[1](x)
            elseif x < x_list[end]
                Base.error("x = $x is outside computed domain x >= $(x_list[end]).")
            elseif x >= 0.0
                Base.error("x = $x should be a negative real number.")
            end
            i = searchsortedlast(x_list, x; rev=true)
            i > length(P_list) && Base.error("x = $x is outside computed ordinary-point domain.")
            return P_list[i](x)
        end

        return x_list, func_list, lfe_P
    else
        coeffs_zero = if N === nothing
            an_zero_auto(s, epsilon, tau, kappa, lambda, _auto_nmax(lambda, epsilon, kappa); tol = tol, patience = _AUTO_PATIENCE)[1]
        else
            an_zero(s, epsilon, tau, kappa, lambda, N)
        end
        coeffs_list = Vector{typeof(coeffs_zero)}()
        x_list = Vector{Float64}()
        x_max = - convergence_radius_an_zero(coeffs_zero, tol) / sqrt(abs(lambda))
        push!(x_list, x_max)
        P, P_prime, _ = evaluate_zero_in(x_max, coeffs_zero)

        while x_max > xm
            coeffs_ord = _ordinary_coeffs(s, epsilon, tau, kappa, lambda, x_max, P, P_prime, N, tol)
            Δx = _ordinary_step_size(convergence_radius_x0_left(x_max, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), convergence_radius_x0_right(x_max, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), N)
            if !(isfinite(Δx)) || Δx <= 0
                Base.error("Step size Δx invalid: $Δx")
            end
            P, P_prime, _ = evaluate_ordinary_point_series(x_max - Δx, x_max, coeffs_ord)
            x_max -= Δx
            push!(coeffs_list, coeffs_ord)
            push!(x_list, x_max)
        end

        function evaluate_P(x)
            if x_list[1] <= x < 0.0
                Px, Px_prime, Px_primeprime = evaluate_zero_in(x, coeffs_zero)
                error = error_zero_in((x, Px, Px_prime, Px_primeprime), s, epsilon, tau, kappa, lambda)
                return (Px, Px_prime, Px_primeprime, error)
            elseif x < x_list[end]
                Base.error("x = $x is outside computed domain x >= $(x_list[end]).")
            elseif x >= 0.0
                Base.error("x = $x should be a negative real number.")
            end
            i = searchsortedlast(x_list, x; rev=true)
            i > length(coeffs_list) && Base.error("x = $x is outside computed ordinary-point domain.")
            x0 = x_list[i]
            coeffs = coeffs_list[i]
            Px, Px_prime, Px_primeprime = evaluate_ordinary_point_series(x, x0, coeffs)
            error = error_ordinary((x, Px, Px_prime, Px_primeprime), s, epsilon, tau, kappa, lambda)
            return (Px, Px_prime, Px_primeprime, error)
        end

        return x_list, coeffs_list, evaluate_P
    end
end

function iterate_zero_out(xm, s, epsilon, tau, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE, lfe = 0)
    if lfe == 1
        Nlfe = min(_fixed_or_auto_N(N), _LFE_INTERNAL_NMAX)
        if abs(epsilon * kappa) < _LFE_MIN_EK
            @warn("LFE zero-out skipped because |epsilon*kappa| = $(abs(epsilon * kappa)) is below the threshold $_LFE_MIN_EK.")
            return iterate_zero_out(xm, s, epsilon, tau, kappa, lambda, N, tol, 0)
        end
        Xout, Rout = lfe_out_zero(s, epsilon, tau, kappa, lambda, Nlfe)
        if !isfinite(Rout) || Rout < _LFE_MIN_RADIUS || ceil(Int, max(0.0, abs(xm) - Rout) / Rout) > _LFE_MAX_STEPS
            @warn("LFE zero-out radius $Rout is too small; falling back to the standard zero expansion.")
            return iterate_zero_out(xm, s, epsilon, tau, kappa, lambda, N, tol, 0)
        end
        x_max = - Rout
        X, X_prime, _, _ = Xout(x_max)
        func_list = Vector{Any}()
        P_list = Vector{Any}()
        x_list = Vector{Float64}()
        push!(x_list, x_max)
        push!(func_list, Xout)
        push!(P_list, X_to_P(Xout, s, epsilon, tau, kappa))

        while x_max > xm
            X0, Δx = lfe_ordinary_point(s, epsilon, tau, kappa, lambda, x_max, X, X_prime, Nlfe, tol)
            if !(isfinite(Δx)) || Δx <= 0
                Base.error("Step size Δx invalid: $Δx")
            end
            X, X_prime, _, _ = X0(x_max - Δx)
            x_max -= Δx
            push!(func_list, X0)
            push!(P_list, X_to_P(X0, s, epsilon, tau, kappa))
            push!(x_list, x_max)
        end

        function lfe_P(x)
            if x_list[1] <= x < 0.0
                return P_list[1](x)
            elseif x < x_list[end]
                Base.error("x = $x is outside computed domain x >= $(x_list[end]).")
            elseif x >= 0.0
                Base.error("x = $x should be a negative real number.")
            end
            i = searchsortedlast(x_list, x; rev=true)
            i > length(P_list) && Base.error("x = $x is outside computed ordinary-point domain.")
            return P_list[i](x)
        end

        return x_list, func_list, lfe_P
    else
        coeffs_zero = if N === nothing
            bn_zero_auto(s, epsilon, tau, kappa, lambda, _auto_nmax(lambda, epsilon, kappa); tol = tol, patience = _AUTO_PATIENCE)[1]
        else
            bn_zero(s, epsilon, tau, kappa, lambda, N)
        end
        coeffs_list = Vector{typeof(coeffs_zero)}()
        x_list = Vector{Float64}()
        x_max = - convergence_radius_bn_zero(coeffs_zero, tol) / sqrt(abs(lambda))
        push!(x_list, x_max)
        P, P_prime, _ = evaluate_zero_out(x_max, s, epsilon, tau, coeffs_zero)

        while x_max > xm
            coeffs_ord = _ordinary_coeffs(s, epsilon, tau, kappa, lambda, x_max, P, P_prime, N, tol)
            Δx = _ordinary_step_size(convergence_radius_x0_left(x_max, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), convergence_radius_x0_right(x_max, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), N)
            if !(isfinite(Δx)) || Δx <= 0
                Base.error("Step size Δx invalid: $Δx")
            end
            P, P_prime, _ = evaluate_ordinary_point_series(x_max - Δx, x_max, coeffs_ord)
            x_max -= Δx
            push!(coeffs_list, coeffs_ord)
            push!(x_list, x_max)
        end

        function evaluate_P(x)
            if x_list[1] <= x < 0.0
                Px, Px_prime, Px_primeprime = evaluate_zero_out(x, s, epsilon, tau, coeffs_zero)
                error = error_zero_out((x, Px, Px_prime, Px_primeprime), s, epsilon, tau, kappa, lambda)
                return (Px, Px_prime, Px_primeprime, error)
            elseif x < x_list[end]
                Base.error("x = $x is outside computed domain x >= $(x_list[end]).")
            elseif x >= 0.0
                Base.error("x = $x should be a negative real number.")
            end
            i = searchsortedlast(x_list, x; rev=true)
            i > length(coeffs_list) && Base.error("x = $x is outside computed ordinary-point domain.")
            x0 = x_list[i]
            coeffs = coeffs_list[i]
            Px, Px_prime, Px_primeprime = evaluate_ordinary_point_series(x, x0, coeffs)
            error = error_ordinary((x, Px, Px_prime, Px_primeprime), s, epsilon, tau, kappa, lambda)
            return (Px, Px_prime, Px_primeprime, error)
        end

        return x_list, coeffs_list, evaluate_P
    end
end

function iterate_inf_in_rho(rho_m, s, epsilon, tau, kappa, lambda, z, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0)
    coeffs_inf, Nopt, R0 = inf_in_coeffs(s, epsilon, tau, kappa, lambda, z, N; tol = tol)
    rho = R0
    P, Pp, _ = evaluate_inf_in(rho, s, epsilon, tau, coeffs_inf, z)
    coeffs_list = Vector{typeof(coeffs_inf)}()
    rho_list = Float64[]
    push!(rho_list, rho)
    push!(coeffs_list, coeffs_inf)
    while rho > rho_m
        coeffs_ord = ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, z, rho, P, Pp, N)
        Δρ = min(
            convergence_radius_rho_left(rho, s, epsilon, tau, kappa, lambda, z, coeffs_ord, tol),
            convergence_radius_rho_right(rho, s, epsilon, tau, kappa, lambda, z, coeffs_ord, tol),
            rho - rho_m,
        )
        Δρ > 0 || error("Step size Δρ invalid: $Δρ")
        rho -= Δρ
        P, Pp, _ = evaluate_ordinary_point_series(rho, rho + Δρ, coeffs_ord)
        push!(rho_list, rho)
        push!(coeffs_list, coeffs_ord)
    end
    function evaluate_P(rho_eval)
        if rho_eval > rho_list[1]
            if rho_eval >= R0
                P0, P0p, P0pp = evaluate_inf_in_rho(rho_eval, s, epsilon, tau, coeffs_inf, z)
                err = error_inf_in_rho(rho_eval, s, epsilon, tau, kappa, lambda, coeffs_inf, z)
                return (P0, P0p, P0pp, err)
            end
            error("rho = $rho_eval is outside the computed domain rho <= $(rho_list[1]).")
        elseif rho_eval < rho_list[end]
            error("rho = $rho_eval is outside the computed domain rho >= $(rho_list[end]).")
        end
        i = searchsortedlast(rho_list, rho_eval; rev = true)
        i == 0 && error("rho = $rho_eval is outside the computed inf_in patch list.")
        if i == length(rho_list)
            rho0 = rho_list[end - 1]
            coeffs = coeffs_list[end]
        else
            rho0 = rho_list[i]
            coeffs = coeffs_list[i + 1]
        end
        Pr, Prp, Prpp = evaluate_ordinary_point_series(rho_eval, rho0, coeffs)
        err = error_ordinary(rho_eval, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
        return (Pr, Prp, Prpp, err)
    end
    return rho_list, coeffs_list, evaluate_P
end

function iterate_inf_out_rho(rho_m, s, epsilon, tau, kappa, lambda, z, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0)
    coeffs_inf, Nopt, R0 = inf_out_coeffs(s, epsilon, tau, kappa, lambda, z, N; tol = tol)
    rho = R0
    P, Pp, _ = evaluate_inf_out(rho, s, epsilon, tau, kappa, coeffs_inf, z)
    coeffs_list = Vector{typeof(coeffs_inf)}()
    rho_list = Float64[]
    push!(rho_list, rho)
    push!(coeffs_list, coeffs_inf)
    while rho > rho_m
        coeffs_ord = ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, z, rho, P, Pp, N)
        Δρ = min(
            convergence_radius_rho_left(rho, s, epsilon, tau, kappa, lambda, z, coeffs_ord, tol),
            convergence_radius_rho_right(rho, s, epsilon, tau, kappa, lambda, z, coeffs_ord, tol),
            rho - rho_m,
        )
        Δρ > 0 || error("Step size Δρ invalid: $Δρ")
        rho -= Δρ
        P, Pp, _ = evaluate_ordinary_point_series(rho, rho + Δρ, coeffs_ord)
        push!(rho_list, rho)
        push!(coeffs_list, coeffs_ord)
    end
    function evaluate_P(rho_eval)
        if rho_eval > rho_list[1]
            if rho_eval >= R0
                P0, P0p, P0pp = evaluate_inf_out_rho(rho_eval, s, epsilon, tau, kappa, coeffs_inf, z)
                err = error_inf_out_rho(rho_eval, s, epsilon, tau, kappa, lambda, coeffs_inf, z)
                return (P0, P0p, P0pp, err)
            end
            error("rho = $rho_eval is outside the computed domain rho <= $(rho_list[1]).")
        elseif rho_eval < rho_list[end]
            error("rho = $rho_eval is outside the computed domain rho >= $(rho_list[end]).")
        end
        i = searchsortedlast(rho_list, rho_eval; rev = true)
        i == 0 && error("rho = $rho_eval is outside the computed inf_out patch list.")
        if i == length(rho_list)
            rho0 = rho_list[end - 1]
            coeffs = coeffs_list[end]
        else
            rho0 = rho_list[i]
            coeffs = coeffs_list[i + 1]
        end
        Pr, Prp, Prpp = evaluate_ordinary_point_series(rho_eval, rho0, coeffs)
        err = error_ordinary(rho_eval, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
        return (Pr, Prp, Prpp, err)
    end
    return rho_list, coeffs_list, evaluate_P
end

function iterate_zero_in_rho(rho_m, s, epsilon, tau, kappa, lambda, z, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0)
    coeffs_zero, R0 = zero_in_coeffs(s, epsilon, tau, kappa, lambda, z, N; tol = tol)
    lambda_scale = iszero(abs(lambda)) ? one(abs(lambda)) : sqrt(abs(lambda))
    rho = R0 / lambda_scale
    P, Pp, _ = evaluate_zero_in_rho(rho, coeffs_zero)
    coeffs_list = Vector{typeof(coeffs_zero)}()
    rho_list = Float64[]
    push!(rho_list, rho)
    push!(coeffs_list, coeffs_zero)
    while rho < rho_m
        coeffs_ord = ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, z, rho, P, Pp, N)
        Δρ = min(
            convergence_radius_rho_right(rho, s, epsilon, tau, kappa, lambda, z, coeffs_ord, tol),
            convergence_radius_rho_left(rho, s, epsilon, tau, kappa, lambda, z, coeffs_ord, tol),
            rho_m - rho,
        )
        Δρ > 0 || error("Step size Δρ invalid: $Δρ")
        rho += Δρ
        P, Pp, _ = evaluate_ordinary_point_series(rho, rho - Δρ, coeffs_ord)
        push!(rho_list, rho)
        push!(coeffs_list, coeffs_ord)
    end
    function evaluate_P(rho_eval)
        if rho_eval < rho_list[1]
            error("rho = $rho_eval is outside the computed domain rho >= $(rho_list[1]).")
        elseif rho_eval <= rho_list[end]
            if rho_eval <= R0
                P0, P0p, P0pp = evaluate_zero_in_rho(rho_eval, coeffs_zero)
                err = error_zero_in_rho(rho_eval, s, epsilon, tau, kappa, lambda, z, coeffs_zero)
                return (P0, P0p, P0pp, err)
            end
            i = searchsortedlast(rho_list, rho_eval)
            i == 0 && error("rho = $rho_eval is outside the computed zero_in patch list.")
            if i == length(rho_list)
                rho0 = rho_list[end - 1]
                coeffs = coeffs_list[end]
            else
                rho0 = rho_list[i]
                coeffs = coeffs_list[i + 1]
            end
            Pr, Prp, Prpp = evaluate_ordinary_point_series(rho_eval, rho0, coeffs)
            err = error_ordinary(rho_eval, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
            return (Pr, Prp, Prpp, err)
        else
            error("rho = $rho_eval is outside the computed domain rho <= $(rho_list[end]).")
        end
    end
    return rho_list, coeffs_list, evaluate_P
end

function iterate_zero_out_resonant_log(xm, s, epsilon, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE)
    Nuse = N === nothing ? _AUTO_NMAX : N
    h, b, log_coefficient = zero_out_resonant_log_coeffs(s, epsilon, kappa, lambda, Nuse)
    x_max = -0.25
    P, P_prime, _ = evaluate_zero_out_resonant_log(x_max, s, h, b, log_coefficient)
    coeffs_list = Vector{Vector{ComplexF64}}()
    x_list = Float64[x_max]

    while x_max > xm
        coeffs_ord = _ordinary_coeffs(s, epsilon, -epsilon, kappa, lambda, x_max, P, P_prime, N, tol)
        Δx = _ordinary_step_size(
            convergence_radius_x0_left(x_max, s, epsilon, -epsilon, kappa, lambda, coeffs_ord, tol),
            convergence_radius_x0_right(x_max, s, epsilon, -epsilon, kappa, lambda, coeffs_ord, tol),
            N,
        )
        if !(isfinite(Δx)) || Δx <= 0
            Base.error("Step size Δx invalid: $Δx")
        end
        P, P_prime, _ = evaluate_ordinary_point_series(x_max - Δx, x_max, coeffs_ord)
        x_max -= Δx
        push!(coeffs_list, coeffs_ord)
        push!(x_list, x_max)
    end

    function evaluate_P(x)
        if x_list[1] <= real(x) < 0.0
            Px, Px_prime, Px_primeprime = evaluate_zero_out_resonant_log(x, s, h, b, log_coefficient)
            error = error_zero_out_resonant_log(x, s, epsilon, kappa, lambda, h, b, log_coefficient)
            return (Px, Px_prime, Px_primeprime, error)
        elseif real(x) < x_list[end]
            Base.error("x = $x is outside computed domain x >= $(x_list[end]).")
        elseif real(x) >= 0.0
            Base.error("x = $x should be a negative real number.")
        end
        i = searchsortedlast(x_list, real(x); rev=true)
        i > length(coeffs_list) && Base.error("x = $x is outside computed ordinary-point domain.")
        x0 = x_list[i]
        coeffs = coeffs_list[i]
        Px, Px_prime, Px_primeprime = evaluate_ordinary_point_series(x, x0, coeffs)
        error = error_ordinary((x, Px, Px_prime, Px_primeprime), s, epsilon, -epsilon, kappa, lambda)
        return (Px, Px_prime, Px_primeprime, error)
    end

    return x_list, coeffs_list, evaluate_P
end

function iterate_zero_in_resonant_log(xm, s, epsilon, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE)
    Nuse = N === nothing ? _AUTO_NMAX : N
    a, b, log_coefficient = zero_in_resonant_log_coeffs(s, epsilon, kappa, lambda, Nuse)
    x_max = -0.25
    P, P_prime, _ = evaluate_zero_in_resonant_log(x_max, s, a, b, log_coefficient)
    coeffs_list = Vector{Vector{ComplexF64}}()
    x_list = Float64[x_max]

    while x_max > xm
        coeffs_ord = _ordinary_coeffs(s, epsilon, -epsilon, kappa, lambda, x_max, P, P_prime, N, tol)
        Δx = _ordinary_step_size(
            convergence_radius_x0_left(x_max, s, epsilon, -epsilon, kappa, lambda, coeffs_ord, tol),
            convergence_radius_x0_right(x_max, s, epsilon, -epsilon, kappa, lambda, coeffs_ord, tol),
            N,
        )
        if !(isfinite(Δx)) || Δx <= 0
            Base.error("Step size Δx invalid: $Δx")
        end
        P, P_prime, _ = evaluate_ordinary_point_series(x_max - Δx, x_max, coeffs_ord)
        x_max -= Δx
        push!(coeffs_list, coeffs_ord)
        push!(x_list, x_max)
    end

    function evaluate_P(x)
        if x_list[1] <= real(x) < 0.0
            Px, Px_prime, Px_primeprime = evaluate_zero_in_resonant_log(x, s, a, b, log_coefficient)
            error = error_zero_in_resonant_log(x, s, epsilon, kappa, lambda, a, b, log_coefficient)
            return (Px, Px_prime, Px_primeprime, error)
        elseif real(x) < x_list[end]
            Base.error("x = $x is outside computed domain x >= $(x_list[end]).")
        elseif real(x) >= 0.0
            Base.error("x = $x should be a negative real number.")
        end
        i = searchsortedlast(x_list, real(x); rev=true)
        i > length(coeffs_list) && Base.error("x = $x is outside computed ordinary-point domain.")
        x0 = x_list[i]
        coeffs = coeffs_list[i]
        Px, Px_prime, Px_primeprime = evaluate_ordinary_point_series(x, x0, coeffs)
        error = error_ordinary((x, Px, Px_prime, Px_primeprime), s, epsilon, -epsilon, kappa, lambda)
        return (Px, Px_prime, Px_primeprime, error)
    end

    return x_list, coeffs_list, evaluate_P
end

function iterate_zero_out_rho(rho_m, s, epsilon, tau, kappa, lambda, z, N = _TruncatioN, tol = _TOLERANCE, sfe = 0, lfe = 0)
    coeffs_zero, R0 = zero_out_coeffs(s, epsilon, tau, kappa, lambda, z, N; tol = tol)
    lambda_scale = iszero(abs(lambda)) ? one(abs(lambda)) : sqrt(abs(lambda))
    rho = R0 / lambda_scale
    P, Pp, _ = evaluate_zero_out_rho(rho, s, epsilon, tau, coeffs_zero, z)
    coeffs_list = Vector{typeof(coeffs_zero)}()
    rho_list = Float64[]
    push!(rho_list, rho)
    push!(coeffs_list, coeffs_zero)
    while rho < rho_m
        coeffs_ord = ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, z, rho, P, Pp, N)
        Δρ = min(
            convergence_radius_rho_right(rho, s, epsilon, tau, kappa, lambda, z, coeffs_ord, tol),
            convergence_radius_rho_left(rho, s, epsilon, tau, kappa, lambda, z, coeffs_ord, tol),
            rho_m - rho,
        )
        Δρ > 0 || error("Step size Δρ invalid: $Δρ")
        rho += Δρ
        P, Pp, _ = evaluate_ordinary_point_series(rho, rho - Δρ, coeffs_ord)
        push!(rho_list, rho)
        push!(coeffs_list, coeffs_ord)
    end
    function evaluate_P(rho_eval)
        if rho_eval < rho_list[1]
            error("rho = $rho_eval is outside the computed domain rho >= $(rho_list[1]).")
        elseif rho_eval <= rho_list[end]
            if rho_eval <= R0
                P0, P0p, P0pp = evaluate_zero_out_rho(rho_eval, s, epsilon, tau, coeffs_zero, z)
                err = error_zero_out_rho(rho_eval, s, epsilon, tau, kappa, lambda, z, coeffs_zero)
                return (P0, P0p, P0pp, err)
            end
            i = searchsortedlast(rho_list, rho_eval)
            i == 0 && error("rho = $rho_eval is outside the computed zero_out patch list.")
            if i == length(rho_list)
                rho0 = rho_list[end - 1]
                coeffs = coeffs_list[end]
            else
                rho0 = rho_list[i]
                coeffs = coeffs_list[i + 1]
            end
            Pr, Prp, Prpp = evaluate_ordinary_point_series(rho_eval, rho0, coeffs)
            err = error_ordinary(rho_eval, rho0, s, epsilon, tau, kappa, lambda, z, coeffs)
            return (Pr, Prp, Prpp, err)
        else
            error("rho = $rho_eval is outside the computed domain rho <= $(rho_list[end]).")
        end
    end
    return rho_list, coeffs_list, evaluate_P
end

end
