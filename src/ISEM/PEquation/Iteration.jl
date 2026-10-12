module Iteration

include("LocalSolutionAtZero.jl")
using .LocalSolutionAtZero
include("LocalSolutionAtInfinity.jl")
using .LocalSolutionAtInfinity
include("OrdinaryPointExpansion.jl")
using .OrdinaryPointExpansion

export iterate_inf_in, iterate_inf_out, iterate_zero_in, iterate_zero_out,
       iterate_inf_in_rho, iterate_inf_out_rho, iterate_zero_in_rho, iterate_zero_out_rho

_TruncatioN = 40
_TOLERANCE = 1e-13

# The contour horizon series is started inside its estimated convergence radius,
# closer in for large |lambda|; for |lambda| < 1 dividing by sqrt|lambda| would
# start it outside the radius (s = +2, l = 2 at small a).
@inline _seed_scale(lambda) = max(one(abs(lambda)), sqrt(abs(lambda)))

# Infinity IN solution of the P equation on the real x axis: start from the asymptotic series and
# step towards the horizon with order-N ordinary-point Taylor patches until x >= xm. Returns the patch
# start points, the patch coefficients and x -> (P, dP/dx, d2P/dx2, residual). iterate_inf_out is the
# OUT analogue.
function iterate_inf_in(xm, s, epsilon, tau, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE)
    coeffs_inf, Nopt, Rin = an_inf(s, epsilon, tau, kappa, lambda, tol)
    x_min = - convergence_radius_an_inf(Rin, s, epsilon, tau, kappa, lambda, coeffs_inf, tol)
    P, P_prime, _ = evaluate_inf_in(x_min, s, epsilon, tau, coeffs_inf)
    coeffs_list = Vector{typeof(coeffs_inf)}()
    x_list = Vector{Float64}()
    push!(x_list, x_min)

    while x_min < xm
        coeffs_ord = ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, x_min, P, P_prime, N)
        Δx = min(convergence_radius_x0_right(x_min, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), convergence_radius_x0_left(x_min, s, epsilon, tau, kappa, lambda, coeffs_ord, tol))
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
            Px, Px_prime, Px_primeprime = evaluate_inf_in(x, s, epsilon, tau, coeffs_inf)
            err = error_inf_in(x, s, epsilon, tau, kappa, lambda, coeffs_inf)
            return (Px, Px_prime, Px_primeprime, err)
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

function iterate_inf_out(xm, s, epsilon, tau, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE)
    coeffs_inf, Nopt, Rout = bn_inf(s, epsilon, tau, kappa, lambda, tol)
    x_min = - convergence_radius_bn_inf(Rout, s, epsilon, tau, kappa, lambda, coeffs_inf, tol)
    P, P_prime, _ = evaluate_inf_out(x_min, s, epsilon, tau, kappa, coeffs_inf)
    coeffs_list = Vector{typeof(coeffs_inf)}()
    x_list = Vector{Float64}()
    push!(x_list, x_min)

    while x_min < xm
        coeffs_ord = ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, x_min, P, P_prime, N)
        Δx = min(convergence_radius_x0_right(x_min, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), convergence_radius_x0_left(x_min, s, epsilon, tau, kappa, lambda, coeffs_ord, tol))
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
            Px, Px_prime, Px_primeprime = evaluate_inf_out(x, s, epsilon, tau, kappa, coeffs_inf)
            err = error_inf_out(x, s, epsilon, tau, kappa, lambda, coeffs_inf)
            return (Px, Px_prime, Px_primeprime, err)
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

# Horizon IN solution: start from the series at x = 0 and step away from the horizon with
# ordinary-point patches until x <= xm. iterate_zero_out is the OUT analogue.
function iterate_zero_in(xm, s, epsilon, tau, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE)
    coeffs_zero = an_zero(s, epsilon, tau, kappa, lambda, N)
    coeffs_list = Vector{typeof(coeffs_zero)}()
    x_list = Vector{Float64}()
    x_max = - convergence_radius_an_zero(coeffs_zero, tol) / sqrt(abs(lambda))
    push!(x_list, x_max)
    P, P_prime, _ = evaluate_zero_in(x_max, coeffs_zero)

    while x_max > xm
        coeffs_ord = ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, x_max, P, P_prime, N)
        Δx = min(convergence_radius_x0_left(x_max, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), convergence_radius_x0_right(x_max, s, epsilon, tau, kappa, lambda, coeffs_ord, tol))
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

function iterate_zero_out(xm, s, epsilon, tau, kappa, lambda, N = _TruncatioN, tol = _TOLERANCE)
    coeffs_zero = bn_zero(s, epsilon, tau, kappa, lambda, N)
    coeffs_list = Vector{typeof(coeffs_zero)}()
    x_list = Vector{Float64}()
    x_max = - convergence_radius_bn_zero(coeffs_zero, tol) / sqrt(abs(lambda))
    push!(x_list, x_max)
    P, P_prime, _ = evaluate_zero_out(x_max, s, epsilon, tau, coeffs_zero)

    while x_max > xm
        coeffs_ord = ordinary_point_coeffs(s, epsilon, tau, kappa, lambda, x_max, P, P_prime, N)
        Δx = min(convergence_radius_x0_left(x_max, s, epsilon, tau, kappa, lambda, coeffs_ord, tol), convergence_radius_x0_right(x_max, s, epsilon, tau, kappa, lambda, coeffs_ord, tol))
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

# Rotated-contour versions (x = -rho/z): step from the infinity series inward to rho_m, or from the
# horizon series outward; the seed series is stored as the first patch.
function iterate_inf_in_rho(rho_m, s, epsilon, tau, kappa, lambda, z, N = _TruncatioN, tol = _TOLERANCE)
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

function iterate_inf_out_rho(rho_m, s, epsilon, tau, kappa, lambda, z, N = _TruncatioN, tol = _TOLERANCE)
    coeffs_inf, Nopt, canonical_R0 = inf_out_coeffs(
        s, epsilon, tau, kappa, lambda, z, N; tol = tol)
    R0 = canonical_R0
    seed_target = max(0.05 * Float64(tol), 2.5 * eps(Float64))
    seed_error = error_inf_out(
        R0, s, epsilon, tau, kappa, lambda, coeffs_inf, z)
    if !isfinite(seed_error) || seed_error > seed_target
        for step in 1:10
            factor = Float64(100 + 5 * step) / 100.0
            candidate_R0 = factor * canonical_R0
            candidate_error = error_inf_out(
                candidate_R0, s, epsilon, tau, kappa, lambda,
                coeffs_inf, z)
            if isfinite(candidate_error) && candidate_error <= seed_target
                R0 = candidate_R0
                break
            end
        end
    end
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

function iterate_zero_in_rho(rho_m, s, epsilon, tau, kappa, lambda, z, N = _TruncatioN, tol = _TOLERANCE)
    coeffs_zero, R0 = zero_in_coeffs(s, epsilon, tau, kappa, lambda, z, N; tol = tol)
    rho = R0 / _seed_scale(lambda)
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
        if rho_eval == rho_list[1]
            P0, P0p, P0pp = evaluate_zero_in_rho(rho_eval, coeffs_zero)
            err = error_zero_in_rho(rho_eval, s, epsilon, tau, kappa, lambda, z, coeffs_zero)
            return (P0, P0p, P0pp, err)
        elseif rho_eval < rho_list[1]
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

function iterate_zero_out_rho(rho_m, s, epsilon, tau, kappa, lambda, z, N = _TruncatioN, tol = _TOLERANCE)
    coeffs_zero, R0 = zero_out_coeffs(s, epsilon, tau, kappa, lambda, z, N; tol = tol)
    rho = R0 / _seed_scale(lambda)
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
        if rho_eval == rho_list[1]
            P0, P0p, P0pp = evaluate_zero_out_rho(rho_eval, s, epsilon, tau, coeffs_zero, z)
            err = error_zero_out_rho(rho_eval, s, epsilon, tau, kappa, lambda, z, coeffs_zero)
            return (P0, P0p, P0pp, err)
        elseif rho_eval < rho_list[1]
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
