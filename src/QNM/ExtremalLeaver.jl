# Exact-extremal Kerr (a = M = 1) radial QNM condition.
#
# Richartz, Phys. Rev. D 93, 064062 (2016) expands the radial Teukolsky function as
#   R = e^{i ω r} e^{J0/(r-1)} (r-1)^{J1} r^{J2} Σ a_n x^n,   x = (r-2)/r,
# with J0 = i(2ω - m), J1 = -2s - 2iω and J2 = -1 + 4iω. The a_n obey the five-term recurrence
# (3.2)-(3.4). The P3 coefficient below contains -2im(1 + 6iω); the printed -2im(1 - 6iω) does not
# satisfy the Teukolsky equation for m != 0 (checked symbolically and numerically).
#
# A QNM is a solution whose series converges absolutely at both ends x = -1 (horizon) and x = 1
# (infinity). Of the four asymptotic solutions, the two that decay forward [the minus sign in (3.5)
# and the plus sign in (3.6)] span the admissible subspace. Backward recursion amplifies exactly
# that subspace; it is tracked with a re-orthonormalised pair of solutions, and the QNM condition
# is that some vector of it also satisfies the boundary rows (3.2) and (3.3). Evaluating the
# decoupled continued fractions of the paper from one truncated elimination chain satisfies (3.11)
# for every ω, so it is not used.
#
# Near ω = m/2 (J0 -> 0) the truncated condition also has zeros that approximate the branch cut of
# the zero-damping accumulation; they move with the truncation depth. Roots are therefore accepted
# only when they are stable under doubling of the depth.

function _extremal_five_term(s, m, omega, lambda, n)
    P1 = 4 * (im * m - im * omega - s)
    P2 = 2 * (1 - 4im * m + 16im * omega)
    P3 = -2 * (1 + 2s + 2lambda + 2im * omega * (4 + 23im * omega) -
        2im * m * (1 + 6im * omega))
    P4 = 4 * (im * m - im * omega + s)
    P5 = -4 * (im * m - im * omega + s) * (1 + 4im * omega)
    P6 = -3 - 8im * omega
    P7 = 2 * (1 + 6im * omega - 8omega^2)
    return (n^2 + n, P1 * n, -2n^2 + P2 * n + P3, P4 * n + P5,
        n^2 + P6 * n + P7)
end

function _orthonormal_pair(u, v)
    u = u / sqrt(sum(abs2, u))
    v = v - sum(conj.(u) .* v) * u
    return u, v / sqrt(sum(abs2, v))
end

"""
    extremal_radial_residual(s, m, omega, angular_A; depth=4000)

Normalised boundary determinant of the exact-extremal (`a = 1`, `M = 1`) radial recurrence at
truncation depth `depth`. It vanishes at a quasinormal frequency.
"""
function extremal_radial_residual(s::Integer, m::Integer, omega::Complex{T},
        angular_A; depth::Int=4000) where {T<:AbstractFloat}
    depth >= 8 || throw(ArgumentError("depth must be at least 8."))
    lambda = Complex{T}(angular_A)
    one_c = one(Complex{T})
    u = Complex{T}[one_c, 0, T(0.3), T(-0.1)]
    v = Complex{T}[0, one_c, T(-0.2), T(0.4)]
    for k in depth:-1:3
        al, be, ga, de, ep = _extremal_five_term(s, m, omega, lambda, T(k))
        un = -(al * u[1] + be * u[2] + ga * u[3] + de * u[4]) / ep
        vn = -(al * v[1] + be * v[2] + ga * v[3] + de * v[4]) / ep
        u = Complex{T}[u[2], u[3], u[4], un]
        v = Complex{T}[v[2], v[3], v[4], vn]
        u, v = _orthonormal_pair(u, v)
    end
    # After the loop the windows hold (a3, a2, a1, a0).
    al1, be1, ga1, _, _ = _extremal_five_term(s, m, omega, lambda, one(T))
    al2, be2, ga2, de2, _ = _extremal_five_term(s, m, omega, lambda, T(2))
    row1(w) = al1 * w[2] + be1 * w[3] + ga1 * w[4]
    row2(w) = al2 * w[1] + be2 * w[2] + ga2 * w[3] + de2 * w[4]
    e11, e12, e21, e22 = row1(u), row1(v), row2(u), row2(v)
    scale = sqrt(abs2(e11) + abs2(e12)) * sqrt(abs2(e21) + abs2(e22))
    return (e11 * e22 - e12 * e21) / scale
end

# Secant iteration on `extremal_radial_residual` at one truncation depth, with
# steps capped at `trust_radius`.
function _extremal_secant(mode::QNMMode, seed::Complex{T}; depth::Int,
        angular_order::Int, sheet_id::Symbol, root_tolerance::T,
        maximum_iterations::Int, trust_radius::T) where {T<:AbstractFloat}
    angular_ok(angular) = angular.status in (
        :continued, :predictor_corrected, :spherical_anchor,
        :high_precision_refined)
    function evaluate(omega)
        angular = angular_branch(mode, one(T), omega;
            truncation_order=angular_order, sheet_id)
        angular_ok(angular) || return (value=Complex{T}(NaN, NaN), angular)
        value = extremal_radial_residual(mode.s, mode.m, omega,
            angular.angular_A; depth)
        return (value=value, angular=angular)
    end
    x0 = seed
    x1 = seed + Complex{T}(T(1e-4) * max(one(T), abs(seed)))
    f0 = evaluate(x0)
    f1 = evaluate(x1)
    rows = NamedTuple[]
    converged = false
    for iteration in 1:maximum_iterations
        (isfinite(f0.value) && isfinite(f1.value)) || break
        f1.value == f0.value && break
        step = -f1.value * (x1 - x0) / (f1.value - f0.value)
        abs(step) > trust_radius && (step *= trust_radius / abs(step))
        x0, f0 = x1, f1
        x1 = x1 + step
        f1 = evaluate(x1)
        push!(rows, (iteration, omega=x1, residual=abs(f1.value),
            step=abs(step)))
        if abs(step) <= root_tolerance * max(one(T), abs(x1))
            converged = isfinite(f1.value)
            break
        end
    end
    return (omega=x1, residual=abs(f1.value), angular=f1.angular,
        converged, rows)
end

"""
    extremal_qnm_root(mode, seed; kwargs...)

Polish a quasinormal frequency of exact-extremal Kerr (`a = 1`) from `seed`. The root is solved
at increasing truncation depths; it is accepted only when two successive depths agree to
`root_tolerance` (relative) and the boundary determinant is below `residual_tolerance`.
"""
function extremal_qnm_root(mode::QNMMode, seed::Complex{T};
        angular_order::Int=40,
        sheet_id::Symbol=:qnm_extremal_exact,
        initial_depth::Int=2000,
        maximum_depth::Int=64000,
        root_tolerance=nothing,
        residual_tolerance=nothing,
        maximum_iterations::Int=60,
        trust_radius=nothing) where {T<:AbstractFloat}
    rtol = root_tolerance === nothing ?
        (T === Float64 ? T(1e-11) : T(1e-20)) : T(root_tolerance)
    ftol = residual_tolerance === nothing ?
        (T === Float64 ? T(1e-11) : T(1e-25)) : T(residual_tolerance)
    radius = trust_radius === nothing ? T(0.01) : T(trust_radius)
    depth = initial_depth
    previous = nothing
    current = seed
    depth_rows = NamedTuple[]
    status = :failed
    stop_reason = :extremal_depth_not_stable
    last = nothing
    while depth <= maximum_depth
        solve = _extremal_secant(mode, current; depth, angular_order,
            sheet_id=Symbol(sheet_id, :_depth_, depth),
            root_tolerance=rtol / 10, maximum_iterations,
            trust_radius=radius)
        last = solve
        push!(depth_rows, (depth, omega=solve.omega,
            residual=solve.residual, converged=solve.converged,
            iterations=length(solve.rows)))
        if !solve.converged
            stop_reason = :extremal_secant_not_converged
            break
        end
        if previous !== nothing
            drift = abs(solve.omega - previous) /
                max(one(T), abs(solve.omega))
            if drift <= rtol && solve.residual <= ftol
                status = :accepted
                stop_reason = :extremal_depth_stable
                break
            end
        end
        previous = solve.omega
        current = solve.omega
        depth *= 2
    end
    drift = length(depth_rows) >= 2 ?
        abs(depth_rows[end].omega - depth_rows[end-1].omega) /
        max(one(T), abs(depth_rows[end].omega)) : T(Inf)
    return (omega=last.omega, angular=last.angular, status, stop_reason,
        residual=T(last.residual), depth_drift=T(drift),
        depth=depth_rows[end].depth, depth_rows,
        seed_shift=T(abs(last.omega - seed)))
end
