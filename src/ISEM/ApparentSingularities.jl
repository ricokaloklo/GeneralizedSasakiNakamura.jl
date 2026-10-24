# Apparent singularities of the GSN equation next to the real x path.
#
# The GSN potentials carry eta'/eta and eta''/eta, so every zero r1 of eta(r) = c0 + c1/r + ... + c4/r^4 is
# a regular singular point of the equation in x, with indices 0 and 2 and both solutions analytic there: an
# apparent singularity. When r1 sits next to the real axis the Taylor patches about real points see a radius
# of only |Im x1|, the steps collapse and the march mixes in the other solution: Y_radial(-2, 2, 2, 0.99,
# 0.75885, IN) (x1 = 0.70896 + 5.4e-7 i, inside the infinity march to match_x = 0.7) was off by 1e-5 for
# r >= 2.5, and the Double64 route overflowed there ("noninvertible q0"). A march that reaches such a point
# stops at 0.35 of the distance rho to the next singular point, takes the local series
# X = sum a_k (x - x1)^k (a0, a2 free, a1 and a_k>=3 from the Laurent coefficients of A and B on a circle
# about x1; the k = 2 resonance is checked) and resumes past x1. Real frequencies only; the same code serves
# the binary64 marches (Iteration, A and B from the closed-form potentials, as the tables lose 1e-12 near the
# horizon) and the arbitrary-precision route (its own tables, in its own number type).
module ApparentSingularities

using ...Transformation: eta_coefficient
using LinearAlgebra: eigvals

export ApparentPoint, eta_roots_x, fixed_singular_points_x, apparent_points, apparent_point_in, apparent_basis, apparent_jump, next_apparent,
    APPARENT_STEP_FRACTION

# The plain march loses accuracy only when the point is very close to the axis: Y of (-2, 2, 2, 0.99) was off by
# 6e-10 at |Im x1| = 2.6e-5 and 1e-5 at 5.4e-7, but good to 1.2e-11 at 6e-4 (omega = 0.76).
const APPARENT_IM_MAX = 1e-3
const APPARENT_STEP_FRACTION = 0.35
const APPARENT_CHECK_TOL = 1e-8
# Binary64 series order and contour points; wider types take more (see _apparent_sizes).
const APPARENT_SERIES_ORDER = 60
const APPARENT_CONTOUR_POINTS = 64

struct ApparentPoint{C}
    x1::C
    radius::Float64
end

_eta_coefficients(p, ::Type{T}) where {T} =
    T[eta_coefficient(p.s, p.m, p.a, p.omega, p.lambda, -k) for k in 0:4]

# The roots of eta r^4 (companion matrix, binary64) as points x = (r - r+) / (r - r-); empty for complex
# frequencies, a outside [0, 1) or a constant eta.
function eta_roots_x(p)
    xs = ComplexF64[]
    (hasproperty(p, :a) && isreal(p.omega) && isreal(p.lambda)) || return xs
    a, w = Float64(real(p.a)), Float64(real(p.omega))
    (0 <= a < 1 && w != 0) || return xs
    c = ComplexF64[eta_coefficient(p.s, p.m, a, w, Float64(real(p.lambda)), -k) for k in 0:4]
    (all(isfinite, c) && abs(c[1]) > 0) || return xs
    degree = findlast(!iszero, c) - 1
    degree >= 1 || return xs
    companion = zeros(ComplexF64, degree, degree)
    companion[1, :] = -c[2:degree+1] ./ c[1]
    for i in 2:degree
        companion[i, i-1] = 1
    end
    kappa = sqrt(1 - a^2)
    for r in eigvals(companion)
        push!(xs, (r - (1 + kappa)) / (r - (1 - kappa)))
    end
    return xs
end

# The other singular points of the GSN equation in x: the horizon, infinity and x(+-i a).
function fixed_singular_points_x(a)
    kappa = sqrt(1 - a^2)
    to_x(r) = (r - (1 + kappa)) / (r - (1 - kappa))
    return a > 0 ? ComplexF64[0.0, 1.0, to_x(im * a), to_x(-im * a)] : ComplexF64[0.0, 1.0]
end

# Apparent points in 0 < Re x < 1 within APPARENT_IM_MAX of the real axis; radius is the distance to the
# nearest other singular point.
function apparent_points(p)
    points = ApparentPoint{ComplexF64}[]
    xs = eta_roots_x(p)
    isempty(xs) && return points
    others = fixed_singular_points_x(Float64(real(p.a)))
    for (i, x1) in enumerate(xs)
        (0 < real(x1) < 1 && abs(imag(x1)) < APPARENT_IM_MAX) || continue
        candidates = vcat(others, ComplexF64[xs[j] for j in eachindex(xs) if j != i])
        push!(points, ApparentPoint(x1, minimum(abs(x1 - z) for z in candidates)))
    end
    return points
end

# The point in Complex{R}: the binary64 root refined by Newton steps on eta r^4 in the type R.
function apparent_point_in(::Type{R}, point::ApparentPoint, p) where {R}
    R === Float64 && return point
    C = Complex{R}
    c = _eta_coefficients(p, C)
    kappa = sqrt(one(R) - R(p.a)^2)
    rp, rm = one(R) + kappa, one(R) - kappa
    r = C((rp - rm * point.x1) / (1 - point.x1))
    for _ in 1:8
        f = zero(C); df = zero(C)
        for k in 0:4
            f = f * r + c[k+1]
            k < 4 && (df = df * r + (4 - k) * c[k+1])
        end
        iszero(df) && break
        dr = f / df
        r -= dr
        abs(dr) <= 4 * eps(R) * abs(r) && break
    end
    return ApparentPoint{C}((r - rp) / (r - rm), point.radius)
end

function _apparent_sizes(::Type{R}) where {R}
    bits = -log2(Float64(eps(R)))
    R === Float64 && return APPARENT_CONTOUR_POINTS, APPARENT_SERIES_ORDER
    # Aliasing of the contour sums goes as 2^-N; the series terms at 0.35 rho as 0.35^K.
    return max(APPARENT_CONTOUR_POINTS, ceil(Int, 1.2 * bits)),
        max(APPARENT_SERIES_ORDER, ceil(Int, bits * log10(2) / 0.45) + 10)
end

# Two local solutions about x1 (a0 = 1, a2 = 0 and a0 = 0, a2 = 1) from the Laurent coefficients of A and B on
# a circle of radius rho / 2; `ab(x)` returns (A(x), B(x)). Nothing if x1 does not behave as an apparent point
# with indices 0 and 2.
function apparent_basis(ab, point::ApparentPoint{C}) where {C}
    R = real(C)
    N, K = _apparent_sizes(R)
    rc = R(point.radius) / 2
    twopi = 2 * convert(R, pi)
    roots = C[cis(twopi * k / N) for k in 0:N-1]
    Avals = Vector{C}(undef, N)
    Bvals = Vector{C}(undef, N)
    for k in 1:N
        Avals[k], Bvals[k] = ab(point.x1 + rc * roots[k])
    end
    # Laurent coefficient j of a function sampled on the circle; offsets: alpha[j + 3], beta[j + 4].
    laurent(vals, j) = sum(vals[k+1] * conj(roots[mod(j * k, N)+1]) for k in 0:N-1) / (N * rc^j)
    alpha = C[laurent(Avals, j) for j in -2:K]
    beta = C[laurent(Bvals, j) for j in -3:K]
    al(j) = alpha[j+3]
    be(j) = beta[j+4]
    scaleA = Float64(maximum(abs, Avals) * rc)
    scaleB = Float64(maximum(abs, Bvals) * rc^2)
    (abs(al(-2)) <= APPARENT_CHECK_TOL * scaleA && abs(be(-3)) <= APPARENT_CHECK_TOL * scaleB &&
        abs(al(-1) + 1) <= APPARENT_CHECK_TOL * max(1.0, scaleA) &&
        abs(be(-2)) <= APPARENT_CHECK_TOL * max(1.0, scaleB)) || return nothing
    L(n) = n * (n - 1) + al(-1) * n + be(-2)
    basis = Vector{Vector{C}}()
    resonance = 0.0
    for start in (0, 2)
        c = zeros(C, K + 1)
        c[start+1] = 1
        for n in 1:K
            n == start && continue
            rhs = zero(C)
            size = zero(R)
            for j in 0:n-1
                term = al(j) * (n - 1 - j) * c[n-j]
                rhs -= term
                size += abs(term)
            end
            for j in -1:n-2
                term = be(j) * c[n-1-j]
                rhs -= term
                size += abs(term)
            end
            if n == 2
                # Resonance: the k = 2 equation has no a2 term at an apparent point and must hold by itself.
                resonance = max(resonance, Float64(abs(rhs) / max(size, floatmin(R))))
                continue
            end
            c[n+1] = rhs / L(n)
        end
        push!(basis, c)
    end
    resonance <= APPARENT_CHECK_TOL || return nothing
    return basis
end

function _series_pair(c, t)
    value = zero(eltype(c))
    derivative = zero(eltype(c))
    power = one(eltype(c))
    for n in 0:length(c)-1
        value += c[n+1] * power
        n < length(c) - 1 && (derivative += (n + 1) * c[n+2] * power)
        power *= t
    end
    return value, derivative
end

# Taylor coefficients about xc of the series with coefficients c about x1.
function _recenter(c, x1, xc)
    K = length(c) - 1
    delta = xc - x1
    out = zeros(eltype(c), K + 1)
    for n in 0:K
        cn = c[n+1]
        iszero(cn) && continue
        term = cn
        for k in n:-1:0
            out[k+1] += term
            term = term * delta * k / (n - k + 1)
        end
    end
    return out
end

# The solution with state (X_a, dX_a) at the real point x_a, as Taylor coefficients about Re x1, and its state
# at x_b; nothing if the local basis cannot represent it.
function apparent_jump(basis, point::ApparentPoint, x_a, X_a, dX_a, x_b)
    xc = real(point.x1)
    p0, d0 = _series_pair(basis[1], x_a - point.x1)
    p2, d2 = _series_pair(basis[2], x_a - point.x1)
    det = p0 * d2 - p2 * d0
    iszero(det) && return nothing
    c0 = (X_a * d2 - p2 * dX_a) / det
    c2 = (p0 * dX_a - X_a * d0) / det
    coeffs = _recenter(c0 .* basis[1] .+ c2 .* basis[2], point.x1, xc)
    X_b, dX_b = _series_pair(coeffs, x_b - xc)
    (isfinite(X_b) && isfinite(dX_b)) || return nothing
    return (center = xc, coeffs = coeffs, X = X_b, dXdx = dX_b)
end

# The next apparent point a march from current_x towards target_x must cross, or nothing.
function next_apparent(points, current_x, target_x, direction)
    best = nothing
    for point in points
        xr = real(point.x1)
        reach = APPARENT_STEP_FRACTION * point.radius
        (direction * (xr - current_x) > -reach && direction * (target_x - xr) > -reach) || continue
        (best === nothing || direction * (xr - real(best.x1)) < 0) && (best = point)
    end
    return best
end

end
