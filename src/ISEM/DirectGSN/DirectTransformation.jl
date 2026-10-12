module DirectTransformation

using ....Coordinates: r_from_rstar, rstar_from_r
using ....ConversionFactors
using ....Solutions: Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix
using ..DirectAsymptoticAmplitudes: _try_unit_pair
using ..DirectMatching:
    DirectRoute,
    DirectConjugatedRoute,
    DirectScaledBasis,
    direct_evaluate,
    direct_state
using ..DirectMatching: direct_route_plan, direct_route_truncations
using ..DirectIteration: _endpoint_last_term_radius, _ENDPOINT_SEED_TARGET
using ..DirectMatching: _mst_evaluation_region

export direct_kappa, direct_r_plus, direct_r_minus
export direct_y_from_x, direct_x_from_y
export direct_x_to_r, direct_r_to_x, direct_y_to_r, direct_r_to_y
export direct_dx_dr, direct_dr_dx, direct_dx_drstar
export direct_evaluate_x, direct_evaluate_y, direct_evaluate_r
export direct_state_x, direct_state_y, direct_state_r, direct_state_rstar
export direct_gsn_solution_rstar, direct_gsn_radial_function
export direct_teukolsky_solution, direct_teukolsky_radial_function
export direct_y_branch_supported, direct_y_solution, direct_y_radial_function

# Teukolsky horizon solutions R_in -> Delta^{-s} e^{-i k r*} and R_out -> e^{+i k r*}
# (leading coefficient one, r* as in Coordinates) at real r near r_+, from the
# Frobenius series of z^2 p2 R'' + z p1 R' + p0 R = 0 in z = r - r_+ (the
# Teukolsky equation times Delta), exponents -s - i u and i u, u = r_+ k / kappa.
# Returns (R_in, dR_in/dr, R_out, dR_out/dr, relative truncation error); the
# series converges for |z| < 2 kappa (the next singular point is r_-).
function direct_teukolsky_horizon_pair(s, m, a, omega, lambda, r::Real)
    kappa = sqrt(1 - a^2)
    rplus = 1 + kappa
    rminus = 1 - kappa
    z = Float64(r) - rplus
    k = omega - m * a / (2 * rplus)
    u = rplus * k / kappa
    K0 = 2 * rplus * omega - a * m
    L0 = lambda - 4im * s * omega * rplus
    # p2 = (z + 2 kappa)^2, p1 = 2 (s + 1)(z + kappa)(z + 2 kappa),
    # p0 = K^2 - 2 i s (z + kappa) K - z (z + 2 kappa)(lambda - 4 i s omega r),
    # K = K0 + 2 r_+ omega z + omega z^2.
    Kc = (K0, 2 * rplus * omega, omega)
    K2 = (Kc[1]^2, 2 * Kc[1] * Kc[2], Kc[2]^2 + 2 * Kc[1] * Kc[3],
        2 * Kc[2] * Kc[3], Kc[3]^2)
    sK = (kappa * Kc[1], Kc[1] + kappa * Kc[2], Kc[2] + kappa * Kc[3], Kc[3])
    zz = (0.0, 2 * kappa, 1.0)                      # z (z + 2 kappa)
    L = (L0, -4im * s * omega)                       # lambda - 4 i s omega r
    zzL = (0.0im, zz[2] * L[1], zz[2] * L[2] + zz[3] * L[1], zz[3] * L[2])
    p2 = (complex(4 * kappa^2), complex(4 * kappa), 1.0 + 0.0im, 0.0im, 0.0im)
    p1 = (complex(4 * (s + 1) * kappa^2), complex(6 * (s + 1) * kappa),
        complex(2 * (s + 1)), 0.0im, 0.0im)
    p0 = ntuple(j -> K2[j] - 2im * s * (j <= 4 ? sK[j] : 0.0im) -
        (j <= 4 ? zzL[j] : 0.0im), 5)
    function series(rho)
        T(j, n) = p2[j + 1] * (n + rho) * (n + rho - 1) + p1[j + 1] * (n + rho) +
            p0[j + 1]
        coefficients = ComplexF64[1.0]
        value = 1.0 + 0.0im
        derivative = rho / z
        zn = 1.0
        n = 0
        tail = Inf
        while n < 400
            n += 1
            acc = 0.0im
            for j in 1:min(n, 4)
                acc += T(j, n - j) * coefficients[n - j + 1]
            end
            c = -acc / T(0, n)
            push!(coefficients, c)
            zn *= z
            term = c * zn
            value += term
            derivative += (n + rho) * term / z
            tail = abs(term) / max(abs(value), floatmin(Float64))
            n >= 4 && tail <= eps(Float64) && break
        end
        power = exp(rho * log(z))
        return power * value, power * derivative, tail
    end
    rho_in = -s - im * u
    rho_out = im * u
    norm_in = (2 * kappa)^(-s) * exp(-im * k * rplus) * exp(im * u * log(2.0)) *
        exp(im * k * rminus / kappa * log(kappa))
    norm_out = exp(im * k * rplus) * exp(-im * u * log(2.0)) *
        exp(-im * k * rminus / kappa * log(kappa))
    Rin, dRin, tin = series(rho_in)
    Rout, dRout, tout = series(rho_out)
    # Delta^{-s} e^{-i k r*} = norm_in z^{rho_in} (1 + O(z)), e^{i k r*} likewise.
    return (Rin * norm_in, dRin * norm_in, Rout * norm_out, dRout * norm_out,
        max(tin, tout))
end

# Teukolsky solution callables: r -> (R, R', R'', error), r -> (R, R') and r -> R.
struct DirectTeukolskySolution{F,G,H}
    tuple_value::F
    pair_value::G
    radial_value::H
end

(solution::DirectTeukolskySolution)(r) = solution.tuple_value(r)

# Asymptotic Teukolsky series at infinity in z = 1 / (omega r), valid for |z| <= radius.
struct DirectTeukolskyInfinitySeries{P}
    params::P
    coefficients::Vector{ComplexF64}
    radius::Float64
    tolerance::Float64
end

function _teukolsky_infinity_coefficients(s, m, a, w, L, order)
    a2 = a*a
    a4 = a2*a2
    # Equation for R=exp(i*w*rstar)*r^(-2s-1)*f(1/r), multiplied by r^-4.
    L2 = ComplexF64.((0im,0im,0im,0im,1+0im,-4+0im,2(a2+2),-4a2,a4))
    L1 = ComplexF64.((0im,0im,-2im*w,2(s+2im*w+1),-2(2im*a2*w+5s+5),
        2(3a2*s+2im*a2*w+3a2+6s+6),-2a2*(im*a2*w+7s+7),4a4*(s+1),0im))
    L0 = ComplexF64.((0im,0im,-L-2a*m*w,
        -2(-L+2im*a2*s*w+im*a2*w-im*a*m*s+2s*s+3s+1),
        -L*a2-2a*a2*m*w+a2*m*m+4a2*s*s+8im*a2*s*w+6a2*s+
            4im*a2*w+2a2-2im*a*m*s+8s*s+12s+4,
        -2a2*(2s+1)*(im*a2*w+3s+3),2a4*(s+1)*(2s+1),0im,0im))
    coefficient = zeros(ComplexF64, order+1)
    coefficient[1] = 1
    at(v,n) = 0<=n<length(v) ? v[n+1] : 0im
    # Store coefficients in 1/(w*r); cancel w in the indicial denominator first.
    for n in 1:order
        acc = 0im
        for j in max(0,n-8):n-1
            acc += (at(L2,n-j+3)*j*(j-1)+at(L1,n-j+2)*j+at(L0,n-j+1))*
                w^(n-j-1)*coefficient[j+1]
        end
        coefficient[n+1] = acc/(2im*n)
    end
    return coefficient
end

function direct_teukolsky_infinity_series(p, order, tolerance)
    coefficient = _teukolsky_infinity_coefficients(p.s, p.m, p.a, p.omega,
        p.lambda, order)
    radius = _endpoint_last_term_radius(coefficient,Float64(tolerance),Inf)
    return DirectTeukolskyInfinitySeries(p,coefficient,radius,Float64(tolerance))
end

# (R, dR/dr, error estimate) from the series, or `nothing` if r is outside the
# radius or the truncation error exceeds the tolerance.
function (series::DirectTeukolskyInfinitySeries)(r)
    p = series.params
    z = inv(p.omega*r)
    abs(z)<=series.radius || return nothing
    c = series.coefficients
    value = c[end]
    derivative = zero(value)
    for j in length(c)-1:-1:1
        derivative = derivative*z+value
        value = value*z+c[j]
    end
    power = -2p.s-1
    radial_factor = power+im*p.omega*(r*r+p.a*p.a)*r/((r-direct_r_plus(p))*(r-direct_r_minus(p)))
    slope = radial_factor*value-z*derivative
    term, norm_value, norm_slope = 1+0im, 0.0, 0.0
    for j in eachindex(c)
        contribution = c[j]*term
        norm_value += abs(contribution)
        norm_slope += abs((radial_factor-(j-1))*contribution)
        term *= z
    end
    ratio = abs(c[end]*z/c[end-1])
    ratio<1 || return nothing
    tail = abs(c[end]*z^(length(c)-1))/(1-ratio)
    error = max((tail+eps(Float64)*norm_value)/abs(value),
        (tail*abs(radial_factor-length(c))+eps(Float64)*norm_slope)/abs(slope))
    error<=series.tolerance || return nothing
    phase = p.omega*rstar_from_r(p.a,r)
    scale = exp(im*phase)*r^power
    return (scale*value,scale*slope/r,error+eps(Float64)*abs(phase))
end

struct DirectTeukolskyPairEvaluator{R,C,T,E}
    route::R
    converter::C
    inv_scale::T
    infinity_series::E
    infinity_scale::ComplexF64
end

# Use the infinity series where it applies (UP, real omega); otherwise convert
# the GSN state at r to (R, R').
function _teukolsky_pair_with_error(e::DirectTeukolskyPairEvaluator,r)
    if e.infinity_series !== nothing
        endpoint = e.infinity_series(r)
        endpoint !== nothing && return e.infinity_scale*endpoint[1],e.infinity_scale*endpoint[2],endpoint[3]
    end
    x = direct_r_to_x(e.route.params, r)
    state = direct_state(e.route, x; radius=r)
    R, Rp = e.converter(r, state.X, state.dXdx)
    return e.inv_scale * R, e.inv_scale * Rp, 0.0
end

function (e::DirectTeukolskyPairEvaluator)(r)
    R,Rp,_ = _teukolsky_pair_with_error(e,r)
    return R,Rp
end

struct DirectTeukolskyRadialEvaluator{P}
    pair::P
end

function (e::DirectTeukolskyRadialEvaluator)(r)
    R, _ = e.pair(r)
    return R
end

struct DirectTeukolskyTupleEvaluator{P,E,F}
    params::P
    pair::E
    d2R::F
end

function (e::DirectTeukolskyTupleEvaluator)(r)
    R, Rp, error = _teukolsky_pair_with_error(e.pair,r)
    p = e.params
    Rpp = e.d2R(r, R, Rp, p.s, p.a, p.omega, p.m, p.lambda)
    return (R, Rp, Rpp, error)
end

# Compactified coordinate x = (r - r_+) / (r - r_-) in [0, 1), and y = 1 - x.
direct_y_from_x(x) = one(x) - x
direct_x_from_y(y) = one(y) - y

function _a_value(a_or_params)
    return hasproperty(a_or_params, :a) ? getproperty(a_or_params, :a) : a_or_params
end

function direct_kappa(a_or_params)
    if hasproperty(a_or_params, :kappa)
        return getproperty(a_or_params, :kappa)
    end
    a = _a_value(a_or_params)
    return sqrt((one(a) - a) * (one(a) + a))
end

direct_r_plus(a_or_params) = one(direct_kappa(a_or_params)) + direct_kappa(a_or_params)
direct_r_minus(a_or_params) = one(direct_kappa(a_or_params)) - direct_kappa(a_or_params)

function direct_x_to_r(a_or_params, x)
    rp = direct_r_plus(a_or_params)
    rm = direct_r_minus(a_or_params)
    return (rp - rm * x) / (one(x) - x)
end

function direct_r_to_x(a_or_params, r)
    rp = direct_r_plus(a_or_params)
    rm = direct_r_minus(a_or_params)
    return (r - rp) / (r - rm)
end

direct_y_to_r(a_or_params, y) = direct_x_to_r(a_or_params, direct_x_from_y(y))
direct_r_to_y(a_or_params, r) = direct_y_from_x(direct_r_to_x(a_or_params, r))

function direct_dr_dx(a_or_params, x)
    kappa = direct_kappa(a_or_params)
    return 2 * kappa / (one(x) - x)^2
end

function direct_dx_dr(a_or_params, x)
    return inv(direct_dr_dx(a_or_params, x))
end

function direct_dx_drstar(a_or_params, x)
    kappa = direct_kappa(a_or_params)
    return kappa * (x - one(x))^2 * x /
        (kappa + (x - one(x))^2 + 2 * kappa^2 * x - kappa * x^2)
end

# x = 1 / (1 + exp(-t)), y = 1 - x, log x and log y, evaluated without overflow.
@inline function _logit_parts(t)
    if t >= 0.0
        inverse = exp(-t)
        log_norm = log1p(inverse)
        x = 1.0 / (1.0 + inverse)
        y = inverse / (1.0 + inverse)
        return x, y, -log_norm, -t - log_norm
    end
    forward = exp(t)
    log_norm = log1p(forward)
    x = forward / (1.0 + forward)
    y = 1.0 / (1.0 + forward)
    return x, y, t - log_norm, -log_norm
end

@inline function _rstar_logit(kappa, rplus, t)
    x, y, log_x, log_y = _logit_parts(t)
    return rplus + 2kappa * exp(t) + 2log(kappa) +
        (rplus / kappa) * log_x - 2log_y
end

@inline function _rstar_slope(kappa, rplus, t)
    x, y, _, _ = _logit_parts(t)
    return 2kappa * exp(t) + (rplus / kappa) * y + 2x
end

# r*(t) and dr*/dt for x = logistic(t).
@inline function _rstar_parts(kappa, rplus, logkappa2,
        rplus_over_kappa, t)
    x, y, log_x, log_y = _logit_parts(t)
    radial = 2kappa * exp(t)
    value = rplus + radial + logkappa2 +
        rplus_over_kappa * log_x - 2log_y
    slope = radial + rplus_over_kappa * y + 2x
    return value, slope
end

# Invert r*(x) with a safeguarded Newton iteration in t = logit(x), starting from
# the horizon or far-zone asymptotic estimate; falls back to r_from_rstar.
function _x_from_rstar(params, target)
    kappa = params.kappa
    a = params.a
    fallback() = direct_r_to_x(params, r_from_rstar(a, target))
    (isfinite(target) && isfinite(kappa) && kappa > 0.0) || return fallback()

    rplus = 1.0 + kappa
    logkappa2 = 2log(kappa)
    rplus_over_kappa = rplus / kappa
    lower = log(nextfloat(0.0))
    upper_x = prevfloat(1.0)
    upper = log(upper_x) - log1p(-upper_x)
    horizon_t = clamp(
        (target - rplus - logkappa2) / rplus_over_kappa,
        lower,
        upper,
    )
    far_scale = target - rplus - logkappa2
    infinity_t = clamp(
        far_scale > 0.0 ? log(far_scale / (2kappa)) : 0.0,
        lower,
        upper,
    )
    horizon_value, _ = _rstar_parts(kappa, rplus, logkappa2,
        rplus_over_kappa, horizon_t)
    infinity_value, _ = _rstar_parts(kappa, rplus, logkappa2,
        rplus_over_kappa, infinity_t)
    horizon_error = abs(horizon_value - target)
    infinity_error = abs(infinity_value - target)
    t = horizon_error <= infinity_error ? horizon_t : infinity_t
    converged = false

    for _ in 1:24
        rstar_value, slope = _rstar_parts(kappa, rplus, logkappa2,
            rplus_over_kappa, t)
        value = rstar_value - target
        if abs(value) <= 8eps(max(abs(target), 1.0))
            converged = true
            break
        end
        value < 0.0 ? (lower = t) : (upper = t)
        step = value / slope
        candidate = t - step
        used_newton = isfinite(candidate) && lower < candidate < upper
        if !used_newton
            candidate = 0.5 * (lower + upper)
        end
        candidate == t && break
        t = candidate
        if used_newton && abs(step) <= 4eps(max(abs(t), 1.0))
            converged = true
            break
        end
    end

    x, _, _, _ = _logit_parts(t)
    if !converged || !isfinite(x) || x <= 1.0e-10 || x >= 1.0 - 1.0e-10
        return fallback()
    end
    return x
end

direct_evaluate_x(route::DirectRoute, x) = direct_evaluate(route, x)
direct_evaluate_y(route::DirectRoute, y) = direct_evaluate(route, direct_x_from_y(y))
direct_evaluate_r(route::DirectRoute, r) = direct_evaluate(route, direct_r_to_x(route.params, r))
direct_evaluate_x(route::DirectConjugatedRoute, x) =
    direct_evaluate(route, x)
direct_evaluate_y(route::DirectConjugatedRoute, y) =
    direct_evaluate(route, direct_x_from_y(y))
direct_evaluate_r(route::DirectConjugatedRoute, r) =
    direct_evaluate(route, direct_r_to_x(route.params, r))

direct_state_x(route::DirectRoute, x) = direct_state(route, x)
direct_state_x(route::DirectConjugatedRoute, x) = direct_state(route, x)

function direct_state_y(route::DirectRoute, y)
    state = direct_state(route, direct_x_from_y(y))
    return (X=state.X, dXdy=-state.dXdx)
end

function direct_state_y(route::DirectConjugatedRoute, y)
    state = direct_state(route, direct_x_from_y(y))
    return (X=state.X, dXdy=-state.dXdx)
end

function direct_state_r(route::DirectRoute, r)
    x = direct_r_to_x(route.params, r)
    state = direct_state(route, x)
    return (X=state.X, dXdr=state.dXdx * direct_dx_dr(route.params, x))
end

function direct_state_r(route::DirectConjugatedRoute, r)
    x = direct_r_to_x(route.params, r)
    state = direct_state(route, x)
    return (X=state.X, dXdr=state.dXdx * direct_dx_dr(route.params, x))
end

# Note: takes x (not r*) and returns X with dX/dr*.
function direct_state_rstar(route::DirectRoute, x)
    state = direct_state(route, x)
    return (X=state.X, dXdrstar=state.dXdx * direct_dx_drstar(route.params, x))
end

function direct_state_rstar(route::DirectConjugatedRoute, x)
    state = direct_state(route, x)
    return (
        X=state.X,
        dXdrstar=state.dXdx * direct_dx_drstar(route.params, x),
    )
end

function _route_boundary_condition(route::DirectRoute)
    root = _root_module()
    route.branch == :IN && return getfield(root, :IN)
    route.branch == :UP && return getfield(root, :UP)
    throw(ArgumentError("direct route branch must be :IN or :UP."))
end

function _root_module()
    return parentmodule(parentmodule(parentmodule(@__MODULE__)))
end

function _route_mode(route::DirectRoute)
    mode_type = getfield(_root_module(), :Mode)
    return mode_type(
        route.params.s,
        route.params.l,
        route.params.m,
        route.params.a,
        route.params.omega,
        route.params.lambda,
    )
end

function _isem_module()
    return parentmodule(parentmodule(@__MODULE__))
end

function _teukolsky_transformation_module()
    isem = _isem_module()
    isdefined(isem, :Matching) || error("ISEM.Matching is not loaded; cannot build direct GSN Teukolsky/Y transformations.")
    matching = getfield(isem, :Matching)
    isdefined(matching, :TeukolskyTransformation) || error("ISEM.Matching.TeukolskyTransformation is not loaded.")
    return getfield(matching, :TeukolskyTransformation)
end

# X and dX/dr* at r; dx/dr* = 2 kappa x / (r^2 + a^2).
function _direct_r_state(route::DirectRoute, r)
    x = direct_r_to_x(route.params, r)
    state = direct_state(route, x; radius=r)
    return state.X, state.dXdx * (2route.params.kappa*x/(r^2+route.params.a^2))
end

# r* -> (X, dX/dr*, 0.0). Points in the MST region, or beyond the matching point
# when the route has an endpoint basis or plan, are evaluated through r.
function direct_gsn_solution_rstar(route::DirectRoute)
    return rs -> begin
        x = _x_from_rstar(route.params, rs)
        if _mst_evaluation_region(route, x) ||
                (x > route.match_x && (route.endpoint_bridge isa DirectScaledBasis ||
                    route.endpoint_plan !== nothing))
            r = r_from_rstar(route.params.a, rs)
            X, dXdrstar = _direct_r_state(route, r)
            return (X, dXdrstar, 0.0)
        end
        state = direct_state(route, x)
        dXdrstar = state.dXdx * direct_dx_drstar(route.params, x)
        return (state.X, dXdrstar, 0.0)
    end
end

function direct_gsn_solution_rstar(route::DirectConjugatedRoute)
    positive_solution = direct_gsn_solution_rstar(route.route)
    return rs -> begin
        X, dXdrstar, residual = positive_solution(rs)
        return (conj(X), conj(dXdrstar), residual)
    end
end

function _unit_gsn_solution(route::DirectRoute)
    scale = route.solution_scale
    scaled_solution = direct_gsn_solution_rstar(route)
    scale == 1.0 && return scaled_solution
    return rs -> begin
        X, dXdrstar, residual = scaled_solution(rs)
        return (X / scale, dXdrstar / scale, residual)
    end
end

# The unit-transmission normalisation divides by solution_scale; it is used only
# if the rescaled state stays at least exp(32) below floatmax at sample points.
const _UNIT_STATE_LOG_MARGIN = 32.0

@inline function _unit_point_safe(route::DirectRoute, x, limit)
    state = direct_state(route, x)
    dXdrstar = state.dXdx * direct_dx_drstar(route.params, x)
    return isfinite(real(state.X)) && isfinite(imag(state.X)) &&
        isfinite(real(dXdrstar)) && isfinite(imag(dXdrstar)) &&
        abs(state.X) <= limit && abs(dXdrstar) <= limit
end

function _unit_basis_safe(route::DirectRoute, basis, limit)
    basis === nothing && return true
    local_basis = hasproperty(basis, :basis) ? basis.basis : basis
    for patch in local_basis.patches
        center = patch.center_x
        endpoint = patch.next_x
        _unit_point_safe(route, center, limit) || return false
        _unit_point_safe(route, (center + endpoint) / 2, limit) || return false
        _unit_point_safe(route, endpoint, limit) || return false
    end
    return true
end

function _unit_state_safe(route::DirectRoute)
    scale = route.solution_scale
    scale == 1.0 && return true
    limit = exp(log(floatmax(Float64)) + log(scale) -
        _UNIT_STATE_LOG_MARGIN)
    isfinite(limit) && limit > 0 || return false
    _unit_point_safe(route, 1.0e-5, limit) || return false
    _unit_point_safe(route, route.match_x, limit) || return false
    _unit_point_safe(route, 1.0 - 1.0e-5, limit) || return false
    _unit_basis_safe(route, route.bridge, limit) || return false
    return _unit_basis_safe(route, route.endpoint_bridge, limit)
end

# Wrap a route as the package's GSNRadialFunction, with unit GSN transmission when
# representable and the route diagnostics in the metadata.
function direct_gsn_radial_function(route::DirectRoute)
    controls = route.controls
    root = _root_module()
    plan = direct_route_plan(route)
    threshold_log = any(
        truncation -> truncation.reason == :threshold_log_companion,
        direct_route_truncations(route),
    )
    unit_pair = _try_unit_pair(route)
    unit_pair_available = unit_pair !== nothing
    unit_state_available = unit_pair_available && _unit_state_safe(route)
    unit_available = unit_pair_available && unit_state_available
    transmission = unit_available ? 1.0 + 0.0im :
        ComplexF64(route.solution_scale)
    incidence, reflection = unit_available ? unit_pair :
        (route.incidence, route.reflection)
    solution = unit_available ? _unit_gsn_solution(route) :
        direct_gsn_solution_rstar(route)
    metadata = (
        xm = route.match_x,
        tol = controls.tolerance,
        horizon_order = controls.horizon_order,
        ordinary_order = controls.ordinary_order,
        infinity_out_order = route.infinity_out_order,
        infinity_in_order = route.infinity_in_order,
        infinity_endpoint_y = route.infinity_endpoint_y,
        patch_count = plan.ordinary_patches + plan.eikonal_patches +
            plan.endpoint_patches,
        sfe = controls.sfe,
        lfe = controls.lfe,
        frequency_regime = controls.frequency_selection.regime,
        frequency_reason = plan.frequency_reason,
        spin_regime = plan.spin,
        spin_reason = plan.spin_reason,
        endpoint_representation = plan.endpoint,
        horizon_representation = threshold_log ?
            :threshold_log_frobenius : :power_frobenius,
        propagation_representation = plan.propagation,
        matching_representation = plan.matching,
        amplitude_observable = threshold_log ?
            :threshold_log_connection_coefficients : :plane_wave_amplitudes,
        eikonal_candidate = plan.eikonal_candidate,
        ordinary_patches = plan.ordinary_patches,
        eikonal_patches = plan.eikonal_patches,
        endpoint_patches = plan.endpoint_patches,
        horizon_patches = plan.horizon_patches,
        infinity_patches = plan.infinity_patches,
        angular_compression = plan.angular_compression,
        phase_compression = plan.phase_compression,
        solution_scale = route.solution_scale,
        unit_pair_available = unit_pair_available,
        unit_state_available = unit_state_available,
        unit_normalization_available = unit_available,
        endpoint_handoff_x = route.endpoint_plan !== nothing ?
            route.endpoint_plan.seed_x :
            plan.endpoint == :direct_endpoint_basis ? route.match_x : missing,
        endpoint_handoff_scale = route.endpoint_plan !== nothing ?
            route.endpoint_plan.scale :
            plan.endpoint == :direct_endpoint_basis ?
                route.endpoint_bridge.scale : missing,
        normalization = unit_available ?
            :unit_gsn_transmission : :scaled_gsn_transmission,
    )
    return getfield(root, :GSNRadialFunction)(
        _route_mode(route),
        _route_boundary_condition(route),
        missing,
        missing,
        rstar_from_r(route.params.a, direct_x_to_r(route.params, route.match_x)),
        controls.horizon_order,
        max(route.infinity_out_order, route.infinity_in_order),
        transmission,
        incidence,
        reflection,
        metadata,
        missing,
        solution,
        getfield(root, :UNIT_GSN_TRANS),
        "GSN-ISEM",
    )
end

function direct_gsn_radial_function(route::DirectConjugatedRoute)
    positive = direct_gsn_radial_function(route.route)
    isem = _isem_module()
    conjugate_function = getfield(isem, :_conjugate_gsn_radial_function)
    return conjugate_function(
        positive,
        route.params.m,
        route.params.omega,
    )
end

function _teukolsky_amplitude_factors(route::DirectRoute)
    p = route.params
    if route.branch == :IN
        return (
            trans = ConversionFactors.Btrans(p.s, p.m, p.a, p.omega, p.lambda),
            incidence = ConversionFactors.Binc(p.s, p.m, p.a, p.omega, p.lambda),
            reflection = ConversionFactors.Bref(p.s, p.m, p.a, p.omega, p.lambda),
        )
    elseif route.branch == :UP
        return (
            trans = ConversionFactors.Ctrans(p.s, p.m, p.a, p.omega, p.lambda),
            incidence = ConversionFactors.Cinc(p.s, p.m, p.a, p.omega, p.lambda),
            reflection = ConversionFactors.Cref(p.s, p.m, p.a, p.omega, p.lambda),
        )
    end
    throw(ArgumentError("direct route branch must be :IN or :UP."))
end

function _normalization_scale_from_factors(route::DirectRoute, normalize, factors)
    root = _root_module()
    if normalize == :teukolsky
        return route.solution_scale * factors.trans,
            getfield(root, :UNIT_TEUKOLSKY_TRANS)
    elseif normalize == :gsn
        return route.solution_scale, getfield(root, :UNIT_GSN_TRANS)
    end
    throw(ArgumentError("normalize must be :teukolsky or :gsn."))
end

function _normalization_scale(route::DirectRoute, normalize)
    return _normalization_scale_from_factors(route, normalize, _teukolsky_amplitude_factors(route))
end

function direct_teukolsky_solution(route::DirectRoute; normalize=:teukolsky, scale=nothing, converter=nothing)
    p = route.params
    scale_value = scale === nothing ? first(_normalization_scale(route, normalize)) : scale
    pair_converter = converter === nothing ? _direct_teukolsky_converter(p) : converter
    return _direct_teukolsky_solution(route, scale_value, pair_converter)
end

function _direct_teukolsky_solution(route::DirectRoute, scale_value, pair_converter)
    p = route.params
    tt = _teukolsky_transformation_module()
    inv_scale = inv(scale_value)
    endpoint = route.branch == :UP && isreal(p.omega) ?
        direct_teukolsky_infinity_series(p,max(route.infinity_out_order,route.controls.ordinary_order),
            max(route.controls.tolerance,_ENDPOINT_SEED_TARGET)) : nothing
    endpoint_scale = ComplexF64(first(_normalization_scale(route,:teukolsky))/scale_value)
    pair_value = DirectTeukolskyPairEvaluator(route,pair_converter,inv_scale,endpoint,endpoint_scale)
    tuple_value = DirectTeukolskyTupleEvaluator(p, pair_value, getfield(tt, :d2R))
    radial_value = DirectTeukolskyRadialEvaluator(pair_value)
    return DirectTeukolskySolution(tuple_value, pair_value, radial_value)
end

# Converter (r, X, dX/dx) -> (R, R'): closed forms for s = -2, 0, +2, otherwise
# the generic GSN-to-Teukolsky conversion matrix.
function _direct_teukolsky_converter(p)
    if p.s == -2
        return DirectTeukolskySMinus2Converter(p.m, p.a, p.omega, p.lambda)
    elseif p.s == 2
        return DirectTeukolskySPlus2Converter(p.m, p.a, p.omega, p.lambda)
    elseif p.s == 0
        return DirectTeukolskySZeroConverter(p.a)
    end
    return DirectTeukolskyGenericConverter(p)
end

struct DirectTeukolskyGenericConverter{P}
    params::P
end

function (c::DirectTeukolskyGenericConverter)(r, X, dXdx)
    p = c.params
    vals = Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(
        p.s,
        p.m,
        p.a,
        p.omega,
        p.lambda,
        r,
    )
    x = direct_r_to_x(p, r)
    dXdrstar = dXdx * direct_dx_drstar(p, x)
    return vals[1, 1] * X + vals[1, 2] * dXdrstar,
        vals[2, 1] * X + vals[2, 2] * dXdrstar
end

struct DirectTeukolskySZeroConverter{T}
    kappa::T
    rminus::T
    a2::T
end

function DirectTeukolskySZeroConverter(a)
    aT = float(a)
    a2 = aT * aT
    kappa = sqrt((one(aT) - aT) * (one(aT) + aT))
    return DirectTeukolskySZeroConverter(kappa, one(kappa) - kappa, a2)
end

function (c::DirectTeukolskySZeroConverter)(r, X, dXdx)
    rho2 = c.a2 + r * r
    invrho = inv(sqrt(rho2))
    offset = r - c.rminus
    R = invrho * X
    Rp = invrho * (
        -(r / rho2) * X + (2 * c.kappa / (offset * offset)) * dXdx)
    return R, Rp
end

struct DirectTeukolskySMinus2Converter{T}
    m::Int
    a::T
    omega::T
    lambda::T
    a2::T
    a3::T
    a4::T
    a5::T
    a6::T
    a7::T
    a8::T
    m2::Int
    lambda2::T
    omega2::T
    omega3::T
    kappa::T
    rminus::T
end

function DirectTeukolskySMinus2Converter(m, a, omega, lambda)
    aT, omegaT, lambdaT = promote(float(a), omega, lambda)
    a2 = aT * aT
    a3 = a2 * aT
    a4 = a2 * a2
    a5 = a4 * aT
    a6 = a3 * a3
    a7 = a6 * aT
    a8 = a4 * a4
    lambda2 = lambdaT * lambdaT
    omega2 = omegaT * omegaT
    omega3 = omega2 * omegaT
    kappa = sqrt(one(aT) - a2)
    rminus = one(kappa) - kappa
    return DirectTeukolskySMinus2Converter{typeof(aT)}(
        Int(m),
        aT,
        omegaT,
        lambdaT,
        a2,
        a3,
        a4,
        a5,
        a6,
        a7,
        a8,
        Int(m) * Int(m),
        lambda2,
        omega2,
        omega3,
        kappa,
        rminus,
    )
end

function (c::DirectTeukolskySMinus2Converter)(r, X, dXdx)
    return _direct_teukolsky_pair_sminus2(
        c.m,
        c.a,
        c.omega,
        c.lambda,
        r,
        X,
        dXdx,
        c.a2,
        c.a3,
        c.a4,
        c.a5,
        c.a6,
        c.a7,
        c.a8,
        c.m2,
        c.lambda2,
        c.omega2,
        c.omega3,
        c.kappa,
        c.rminus,
    )
end

struct DirectTeukolskySPlus2Converter{T}
    m::Int
    a::T
    omega::T
    lambda::T
    a2::T
    a3::T
    a4::T
    a5::T
    a6::T
    a7::T
    a8::T
    m2::Int
    lambda2::T
    omega2::T
    omega3::T
    kappa::T
    rminus::T
end

function DirectTeukolskySPlus2Converter(m, a, omega, lambda)
    aT, omegaT, lambdaT = promote(float(a), omega, lambda)
    a2 = aT * aT
    a3 = a2 * aT
    a4 = a2 * a2
    a5 = a4 * aT
    a6 = a3 * a3
    a7 = a6 * aT
    a8 = a4 * a4
    lambda2 = lambdaT * lambdaT
    omega2 = omegaT * omegaT
    omega3 = omega2 * omegaT
    kappa = sqrt(one(aT) - a2)
    rminus = one(kappa) - kappa
    return DirectTeukolskySPlus2Converter{typeof(aT)}(
        Int(m),
        aT,
        omegaT,
        lambdaT,
        a2,
        a3,
        a4,
        a5,
        a6,
        a7,
        a8,
        Int(m) * Int(m),
        lambda2,
        omega2,
        omega3,
        kappa,
        rminus,
    )
end

function (c::DirectTeukolskySPlus2Converter)(r, X, dXdx)
    return _direct_teukolsky_pair_splus2(
        c.m,
        c.a,
        c.omega,
        c.lambda,
        r,
        X,
        dXdx,
        c.a2,
        c.a3,
        c.a4,
        c.a5,
        c.a6,
        c.a7,
        c.a8,
        c.m2,
        c.lambda2,
        c.omega2,
        c.omega3,
        c.kappa,
        c.rminus,
    )
end

function _direct_teukolsky_matrix_sminus2(m, a, omega, lambda, r)
    a2 = a^2
    a3 = a2 * a
    a4 = a2^2
    a5 = a4 * a
    a6 = a3^2
    a7 = a6 * a
    a8 = a4^2
    m2 = m * m
    lambda2 = lambda * lambda
    omega2 = omega * omega
    omega3 = omega2 * omega
    return _direct_teukolsky_matrix_sminus2(
        m,
        a,
        omega,
        lambda,
        r,
        a2,
        a3,
        a4,
        a5,
        a6,
        a7,
        a8,
        m2,
        lambda2,
        omega2,
        omega3,
    )
end

function _direct_teukolsky_matrix_sminus2(
    m,
    a,
    omega,
    lambda,
    r,
    a2,
    a3,
    a4,
    a5,
    a6,
    a7,
    a8,
    m2,
    lambda2,
    omega2,
    omega3,
)
    delta = a2 + (-2 + r) * r
    r2 = r^2
    r3 = r2 * r
    r4 = r2^2
    r5 = r4 * r
    r6 = r3^2
    r7 = r6 * r
    r8 = r4^2
    denom = (-r4) * (2 * lambda + lambda2 - 12im * omega) +
        24 * a3 * m * r * (im - 2 * r * omega) -
        4 * a * m * r2 * (6im + 2im * r * lambda + 3 * r2 * omega) +
        12 * a4 * (-1 - 2im * r * omega + 2 * r2 * omega2) +
        4 * a2 * r * (6 + r * (-3 + 6 * m2 + 6im * omega) +
        2im * r2 * (-3 + lambda) * omega + 3 * r3 * omega2)
    invdenom = inv(denom)
    m11_num = r2 * (
        -4 * a5 * m * r * (-im + r * omega) -
        2 * a * m * r4 * (im + 2 * r2 * omega) -
        2 * a3 * m * r2 * (3im - 2im * r + 4 * r2 * omega) +
        2 * a6 * (-5 - 2im * r * omega + r2 * omega2) +
        a4 * r * (32 + r * (-20 + 2 * m2 - lambda + 6im * omega) -
        10im * r2 * omega + 6 * r3 * omega2) +
        r4 * (-12 + 2 * r * (5 + lambda) -
        r2 * (2 + lambda - 6im * omega) - 2im * r3 * omega +
        2 * r4 * omega2) +
        2 * a2 * r2 * (-12 + r * (19 + lambda) +
        r2 * (-6 + m2 - lambda + 6im * omega) - 4im * r3 * omega +
        3 * r4 * omega2)
    )
    sqrtfac = sqrt((a2 + r2) / delta^2)
    m11 = m11_num * invdenom / (delta^3 * sqrtfac^3)
    m12 = (
        2 * r3 * delta * sqrtfac *
        (im * a * m * r + a2 * (-2 - im * r * omega) +
        r * (3 - r - im * r2 * omega))
    ) * invdenom
    m21_num = 2 * a7 * m * r * (2 + 4im * r * omega - 3 * r2 * omega2) -
        2 * a5 * m * r2 * (4 + r * (-1 + m2 - lambda + 6im * omega) -
        12im * r2 * omega + 9 * r3 * omega2) -
        2 * a3 * m * r4 * (-5 + 2 * lambda + r * (3 + m2 - 2 * lambda +
        16im * omega) - 12im * r2 * omega + 9 * r3 * omega2) -
        2 * a * m * r5 * (6 + r * (-7 + 2 * lambda) +
        r2 * (2 - lambda + 10im * omega) - 4im * r3 * omega +
        3 * r4 * omega2) +
        2 * a8 * (6im - 2 * r * omega - 2im * r2 * omega2 + r3 * omega3) +
        2 * a6 * r * (-36im + r2 * (-12 + 3 * m2 - lambda + 3im * omega) *
        omega - 10im * r3 * omega2 + 4 * r4 * omega3 +
        r * (21im - 2im * m2 + 2im * lambda + 4 * omega)) +
        r5 * (-48im + 12im * r * (6 + lambda) -
        12im * r2 * (3 + lambda - 3im * omega) -
        2 * r4 * (4 + lambda - 9im * omega) * omega -
        8im * r5 * omega2 + 2 * r6 * omega3 +
        r3 * (6im + 3im * lambda + 34 * omega + 4 * lambda * omega)) +
        2 * a2 * r3 * (-48im + 4im * r * (27 + 2 * lambda) +
        im * r2 * (-72 + m2 - 14 * lambda + 30im * omega) +
        r4 * (-16 + 3 * m2 - 3 * lambda + 21im * omega) * omega -
        14im * r5 * omega2 + 4 * r6 * omega3 +
        r3 * (15im + 5im * lambda + 48 * omega + 4 * lambda * omega)) +
        a4 * r2 * (144im + 2im * r * (-90 + 3 * m2 - 8 * lambda) +
        2 * r3 * (-22 + 6 * m2 - 3 * lambda + 15im * omega) * omega -
        36im * r4 * omega2 + 12 * r5 * omega3 +
        r2 * (54im - 4im * m2 + 11im * lambda + 70 * omega + 4 * lambda * omega))
    m21 = im * r * sqrtfac * m21_num * invdenom / ((a2 + r2)^2)
    m22 = (
        r2 * sqrtfac *
        (-4 * a3 * m * r * (-im + r * omega) -
        2 * a * m * r2 * (3im - im * r + 2 * r2 * omega) +
        2 * a4 * (-3 - 2im * r * omega + r2 * omega2) +
        a2 * r * (24 + r * (-12 + 2 * m2 - lambda + 6im * omega) -
        12im * r2 * omega + 4 * r3 * omega2) +
        r2 * (-24 + 2 * r * (12 + lambda) -
        r2 * (6 + lambda - 18im * omega) - 8im * r3 * omega +
        2 * r4 * omega2))
    ) * invdenom
    return m11, m12, m21, m22
end

# Explicit s = -2 GSN-to-Teukolsky conversion applied to (X, dX/dx).
function _direct_teukolsky_pair_sminus2(
    m,
    a,
    omega,
    lambda,
    r,
    X,
    dXdx,
    a2,
    a3,
    a4,
    a5,
    a6,
    a7,
    a8,
    m2,
    lambda2,
    omega2,
    omega3,
    kappa,
    rminus,
)
    delta = a2 + (-2 + r) * r
    r2 = r * r
    r3 = r2 * r
    r4 = r2 * r2
    r5 = r4 * r
    r6 = r3 * r3
    r7 = r6 * r
    r8 = r4 * r4
    rshift = r - rminus
    dxdrstar = 2 * kappa * delta / (rshift * rshift * (a2 + r2))
    denom = (-r4) * (2 * lambda + lambda2 - 12im * omega) +
        24 * a3 * m * r * (im - 2 * r * omega) -
        4 * a * m * r2 * (6im + 2im * r * lambda + 3 * r2 * omega) +
        12 * a4 * (-1 - 2im * r * omega + 2 * r2 * omega2) +
        4 * a2 * r * (6 + r * (-3 + 6 * m2 + 6im * omega) +
        2im * r2 * (-3 + lambda) * omega + 3 * r3 * omega2)
    invdenom = inv(denom)
    m11_num = r2 * (
        -4 * a5 * m * r * (-im + r * omega) -
        2 * a * m * r4 * (im + 2 * r2 * omega) -
        2 * a3 * m * r2 * (3im - 2im * r + 4 * r2 * omega) +
        2 * a6 * (-5 - 2im * r * omega + r2 * omega2) +
        a4 * r * (32 + r * (-20 + 2 * m2 - lambda + 6im * omega) -
        10im * r2 * omega + 6 * r3 * omega2) +
        r4 * (-12 + 2 * r * (5 + lambda) -
        r2 * (2 + lambda - 6im * omega) - 2im * r3 * omega +
        2 * r4 * omega2) +
        2 * a2 * r2 * (-12 + r * (19 + lambda) +
        r2 * (-6 + m2 - lambda + 6im * omega) - 4im * r3 * omega +
        3 * r4 * omega2)
    )
    sqrtfac = sqrt((a2 + r2) / delta^2)
    m11 = m11_num * invdenom / (delta^3 * sqrtfac^3)
    m12_dx = (
        2 * r3 * delta * sqrtfac *
        (im * a * m * r + a2 * (-2 - im * r * omega) +
        r * (3 - r - im * r2 * omega))
    ) * invdenom * dxdrstar
    m21_num = 2 * a7 * m * r * (2 + 4im * r * omega - 3 * r2 * omega2) -
        2 * a5 * m * r2 * (4 + r * (-1 + m2 - lambda + 6im * omega) -
        12im * r2 * omega + 9 * r3 * omega2) -
        2 * a3 * m * r4 * (-5 + 2 * lambda + r * (3 + m2 - 2 * lambda +
        16im * omega) - 12im * r2 * omega + 9 * r3 * omega2) -
        2 * a * m * r5 * (6 + r * (-7 + 2 * lambda) +
        r2 * (2 - lambda + 10im * omega) - 4im * r3 * omega +
        3 * r4 * omega2) +
        2 * a8 * (6im - 2 * r * omega - 2im * r2 * omega2 + r3 * omega3) +
        2 * a6 * r * (-36im + r2 * (-12 + 3 * m2 - lambda + 3im * omega) *
        omega - 10im * r3 * omega2 + 4 * r4 * omega3 +
        r * (21im - 2im * m2 + 2im * lambda + 4 * omega)) +
        r5 * (-48im + 12im * r * (6 + lambda) -
        12im * r2 * (3 + lambda - 3im * omega) -
        2 * r4 * (4 + lambda - 9im * omega) * omega -
        8im * r5 * omega2 + 2 * r6 * omega3 +
        r3 * (6im + 3im * lambda + 34 * omega + 4 * lambda * omega)) +
        2 * a2 * r3 * (-48im + 4im * r * (27 + 2 * lambda) +
        im * r2 * (-72 + m2 - 14 * lambda + 30im * omega) +
        r4 * (-16 + 3 * m2 - 3 * lambda + 21im * omega) * omega -
        14im * r5 * omega2 + 4 * r6 * omega3 +
        r3 * (15im + 5im * lambda + 48 * omega + 4 * lambda * omega)) +
        a4 * r2 * (144im + 2im * r * (-90 + 3 * m2 - 8 * lambda) +
        2 * r3 * (-22 + 6 * m2 - 3 * lambda + 15im * omega) * omega -
        36im * r4 * omega2 + 12 * r5 * omega3 +
        r2 * (54im - 4im * m2 + 11im * lambda + 70 * omega + 4 * lambda * omega))
    m21 = im * r * sqrtfac * m21_num * invdenom / ((a2 + r2)^2)
    m22_dx = (
        r2 * sqrtfac *
        (-4 * a3 * m * r * (-im + r * omega) -
        2 * a * m * r2 * (3im - im * r + 2 * r2 * omega) +
        2 * a4 * (-3 - 2im * r * omega + r2 * omega2) +
        a2 * r * (24 + r * (-12 + 2 * m2 - lambda + 6im * omega) -
        12im * r2 * omega + 4 * r3 * omega2) +
        r2 * (-24 + 2 * r * (12 + lambda) -
        r2 * (6 + lambda - 18im * omega) - 8im * r3 * omega +
        2 * r4 * omega2))
    ) * invdenom * dxdrstar
    return m11 * X + m12_dx * dXdx, m21 * X + m22_dx * dXdx
end

# Explicit s = +2 GSN-to-Teukolsky conversion applied to (X, dX/dx).
function _direct_teukolsky_pair_splus2(
    m,
    a,
    omega,
    lambda,
    r,
    X,
    dXdx,
    a2,
    a3,
    a4,
    a5,
    a6,
    a7,
    a8,
    m2,
    lambda2,
    omega2,
    omega3,
    kappa,
    rminus,
)
    delta = a2 + (-2 + r) * r
    r2 = r * r
    r3 = r2 * r
    r4 = r2 * r2
    r5 = r4 * r
    r6 = r3 * r3
    r7 = r6 * r
    r8 = r4 * r4
    sumr2 = a2 + r2
    rshift = r - rminus
    dxdrstar = 2 * kappa * delta / (rshift * rshift * sumr2)
    denom = (-r4) * (24 + 10 * lambda + lambda2 + 12im * omega) -
        24 * a3 * m * r * (im + 2 * r * omega) +
        4 * a * m * r2 * (6im + 2im * r * (4 + lambda) - 3 * r2 * omega) +
        12 * a4 * (-1 + 2im * r * omega + 2 * r2 * omega2) +
        4 * a2 * r * (6 + r * (-3 + 6 * m2 - 6im * omega) -
        2im * r2 * (1 + lambda) * omega + 3 * r3 * omega2)
    common_inv = inv(delta^3 * sumr2 * sqrt(sumr2) * denom)
    m11_num = r2 * delta * (
        -4 * a5 * m * r * (im + r * omega) -
        2 * a3 * m * r2 * (-3im + 2im * r + 4 * r2 * omega) +
        a * m * (2im * r4 - 4 * r6 * omega) +
        2 * a6 * (-5 + 2im * r * omega + r2 * omega2) +
        a4 * r * (32 + r * (-24 + 2 * m2 - lambda - 6im * omega) +
        10im * r2 * omega + 6 * r3 * omega2) +
        r4 * (-12 + 2 * r * (9 + lambda) -
        r2 * (6 + lambda + 6im * omega) + 2im * r3 * omega +
        2 * r4 * omega2) +
        2 * a2 * r2 * (-12 + r * (23 + lambda) +
        r2 * (-10 + m2 - lambda - 6im * omega) + 4im * r3 * omega +
        3 * r4 * omega2)
    )
    m11 = m11_num * common_inv
    m12_dx = (
        2 * r3 * delta * sumr2^2 *
        (-im * a * m * r + a2 * (-2 + im * r * omega) +
        r * (3 - r + im * r2 * omega))
    ) * common_inv * dxdrstar
    m21_num =
        -2 * a7 * m * r * (-2 + 4im * r * omega + 3 * r2 * omega2) -
        2 * a5 * m * r2 * (-4 + r * (3 + m2 - lambda + 2im * omega) +
        4im * r2 * omega + 9 * r3 * omega2) +
        2 * a * m * r5 * (-10 + r * (3 - 2 * lambda) +
        r2 * (2 + lambda + 2im * omega) + 4im * r3 * omega -
        3 * r4 * omega2) -
        2 * a3 * m * r3 *
        (12 + r2 * (3 + m2 - 2 * lambda) + r * (-17 + 2 * lambda) -
        4im * r3 * omega + 9 * r4 * omega2) +
        2 * a8 * (-6im - 2 * r * omega + 2im * r2 * omega2 + r3 * omega3) +
        2 * a6 * r * (16im +
        im * r * (-9 + 2 * m2 - 2 * lambda + 4im * omega) +
        r2 * (-8 + 3 * m2 - lambda + im * omega) * omega +
        6im * r3 * omega2 + 4 * r4 * omega3) +
        r6 * (-4im * lambda - 2 * r3 * (4 + lambda + 5im * omega) * omega +
        2 * r5 * omega3 - 12 * r * (im + omega) +
        r2 * (6im + im * lambda + 18 * omega + 4 * lambda * omega)) +
        a4 * r2 * (-16im + 2im * r * (10 + m2 + 6 * lambda -
        12im * omega) + 2 * r3 * (-14 + 6 * m2 - 3 * lambda -
        3im * omega) * omega + 12im * r4 * omega2 + 12 * r5 * omega3 +
        r2 * (-2im - 4im * m2 - 7im * lambda + 22 * omega +
        4 * lambda * omega)) +
        2 * a2 * r4 * (-4im * lambda +
        im * r * (-4 + 3 * m2 + 6 * lambda + 6im * omega) +
        3 * r3 * (-4 + m2 - lambda - 3im * omega) * omega +
        2im * r4 * omega2 + 4 * r5 * omega3 +
        r2 * (5im - 4im * m2 - im * lambda + 24 * omega +
        4 * lambda * omega))
    m21 = -im * r * m21_num * common_inv
    m22_dx = (
        r2 * sumr2^2 *
        (-4 * a3 * m * r * (im + r * omega) -
        2 * a * m * r2 * (im - 3im * r + 2 * r2 * omega) +
        2 * a4 * (-3 + 2im * r * omega + r2 * omega2) +
        r3 * (2 * lambda - r * (2 + lambda + 10im * omega) + 2 * r3 * omega2) +
        a2 * r * (8 + 2 * m2 * r - r * lambda + 2im * r * omega +
        4im * r2 * omega + 4 * r3 * omega2))
    ) * common_inv * dxdrstar
    return m11 * X + m12_dx * dXdx, m21 * X + m22_dx * dXdx
end

# Wrap a route as the package's TeukolskyRadialFunction, normalised to unit
# Teukolsky (:teukolsky) or unit GSN (:gsn) transmission.
function direct_teukolsky_radial_function(route::DirectRoute; normalize=:teukolsky)
    root = _root_module()
    factors = _teukolsky_amplitude_factors(route)
    scale, convention = _normalization_scale_from_factors(route, normalize, factors)
    solution = direct_teukolsky_solution(route; normalize=normalize, scale=scale)
    return getfield(root, :TeukolskyRadialFunction)(
        _route_mode(route),
        _route_boundary_condition(route),
        route.solution_scale * factors.trans / scale,
        route.incidence * factors.incidence / scale,
        route.reflection * factors.reflection / scale,
        missing,
        direct_gsn_radial_function(route),
        solution,
        convention,
    )
end

function direct_teukolsky_radial_function(
    route::DirectConjugatedRoute;
    normalize=:teukolsky,
)
    positive = direct_teukolsky_radial_function(
        route.route;
        normalize=normalize,
    )
    isem = _isem_module()
    conjugate_function =
        getfield(isem, :_conjugate_teukolsky_radial_function)
    return conjugate_function(
        positive,
        route.params.m,
        route.params.omega,
    )
end

direct_y_branch_supported(route::DirectRoute) =
    (route.params.s == -2 && route.branch == :IN) || (route.params.s == 2 && route.branch == :UP)
direct_y_branch_supported(route::DirectConjugatedRoute) =
    direct_y_branch_supported(route.route)

function _check_y_branch(route::DirectRoute)
    direct_y_branch_supported(route) && return nothing
    error("direct Y transformation supports only s = -2 with IN and s = +2 with UP.")
end

function _direct_y_amplitudes(route::DirectRoute)
    p = route.params
    if p.s == -2 && route.branch == :IN
        return (
            trans = route.solution_scale *
                ConversionFactors.Btrans(p.s, p.m, p.a, p.omega, p.lambda),
            incidence = ConversionFactors.Binc(p.s, p.m, p.a, p.omega, p.lambda) * route.incidence,
            reflection = ConversionFactors.Bref(p.s, p.m, p.a, p.omega, p.lambda) * route.reflection,
        )
    elseif p.s == 2 && route.branch == :UP
        return (
            trans = route.solution_scale *
                ConversionFactors.Ctrans(p.s, p.m, p.a, p.omega, p.lambda),
            incidence = ConversionFactors.Cinc(p.s, p.m, p.a, p.omega, p.lambda) * route.incidence,
            reflection = ConversionFactors.Cref(p.s, p.m, p.a, p.omega, p.lambda) * route.reflection,
        )
    end
    error("unsupported direct Y branch.")
end

# Y-function solution built from the GSN solution; supported for s = -2 IN and s = +2 UP.
function direct_y_solution(route::DirectRoute; amplitudes=nothing)
    _check_y_branch(route)
    iszero(route.params.omega) && error("direct Y transformation requires nonzero omega.")
    p = route.params
    tt = _teukolsky_transformation_module()
    amps = amplitudes === nothing ? _direct_y_amplitudes(route) : amplitudes
    teukolsky_from_gsn_matrix = r -> Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(
        p.s,
        p.m,
        p.a,
        p.omega,
        p.lambda,
        r,
    )
    return getfield(tt, :GSN_to_Y_solution_from_matrix)(
        r -> rstar_from_r(p.a, r),
        teukolsky_from_gsn_matrix,
        direct_gsn_solution_rstar(route),
        p.s,
        p.m,
        p.a,
        p.omega,
        p.lambda,
        amps.trans,
    )
end

function direct_y_radial_function(route::DirectRoute)
    _check_y_branch(route)
    isem = _isem_module()
    isdefined(isem, :YRadialFunction) || error("ISEM.YRadialFunction is not loaded.")
    y_type = getfield(isem, :YRadialFunction)
    amps = _direct_y_amplitudes(route)
    ysoln = direct_y_solution(route; amplitudes=amps)
    xsoln = r -> direct_gsn_solution_rstar(route)(rstar_from_r(route.params.a, r))[1]
    yscalar = r -> ysoln(r)[1]
    return y_type(
        _route_mode(route),
        _route_boundary_condition(route),
        amps.trans,
        amps.incidence,
        amps.reflection,
        missing,
        missing,
        xsoln,
        yscalar,
        ysoln,
        getfield(_root_module(), :UNIT_TEUKOLSKY_TRANS),
    )
end

function direct_y_radial_function(route::DirectConjugatedRoute)
    positive = direct_y_radial_function(route.route)
    isem = _isem_module()
    conjugate_function = getfield(isem, :_conjugate_y_radial_function)
    return conjugate_function(
        positive,
        route.params.m,
        route.params.omega,
    )
end

end
