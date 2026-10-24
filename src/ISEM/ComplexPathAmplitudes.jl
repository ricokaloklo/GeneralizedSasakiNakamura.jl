module ComplexPathAmplitudes

# Asymptotic amplitudes of the Teukolsky IN and UP solutions at real frequency from the Teukolsky equation
# continued along a path in the complex r plane, in Float64.
#
# At high frequency the reflected part of a basis solution is exponentially small next to the incident part
# (|B_ref / B_inc| = 2e-17 at (s, l, m, a, omega) = (-2, 2, 2, 0.9, 3), 1e-56 at (-2, 4, 4, 0.9, 10)). On the real
# axis it is then carried by the last digits of numbers of the size of B_inc, and any Float64 route, ISEM or a
# real-axis march, leaves it wrong by about eps |B_inc|. Below the real axis (Im r* < 0 for omega > 0) the
# outgoing wave e^{i omega r*} grows and the ingoing one decays by e^{-2 omega |Im r*|}, so there the reflected
# part dominates and can be projected out at its own relative precision. Nothing here is a large-frequency
# expansion: the horizon Frobenius series, an ordinary-point Taylor march of the Teukolsky equation and the
# series at infinity are the exact ones, only evaluated at complex r.
#
# IN (unit transmission R -> Delta^{-s} e^{-i k r*} at the horizon): B_inc from a march along the real axis,
# B_ref from a march that leaves the horizon series at r_+ + kappa e^{-i theta} (theta near pi), drops to
# Im r = -depth and runs out to X - i depth, where both are projected on the infinity series.
# UP (unit transmission R -> r^{-2s-1} e^{i omega r*} at infinity): C_up (horizon e^{i k r*}) along the real axis,
# C_ref (horizon Delta^{-s} e^{-i k r*}) along the mirror path above the axis, run inward from X + i depth to
# r_+ + kappa e^{i theta} (the path below the axis leaves C_ref wrong by 1e11 - 1e42 on the 896-bit grid).
#
# Measured against BigFloat(256) references (B_ref / B_inc relative error): 2e-14 at omega = 2.35 and 3, 1e-13 at
# omega = 6, 2e-11 at omega = 10, where the real-axis Float64 routes give 0.3, 8e3, 2e21 and 3e44.
#
# All state is local to one call; nothing is cached, so calls are thread safe.

export complex_path_amplitudes

const TAYLOR_ORDER = 40
const SEED_THETA = 3.0
const PATH_DEPTH = 10.0
const SERIES_TOL = 1e-15
const SEED_CONDITION = 1e2
const STEP_WAVELENGTHS = 0.5
const PATH_LEFT = 1.5
const NEAR_EXTREMAL_KAPPA = 0.1
const PASS_KAPPAS = 5.0
const PROJECTION_RADII = (4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 90.0, 120.0, 180.0, 250.0, 400.0, 600.0)

const C = ComplexF64

_horizons(a) = (k = sqrt(1 - a^2); (k, 1 + k, 1 - k))

# Tortoise coordinate in the package convention (Coordinates.rstar_from_r), continued to complex r.
function _rstar(a, r)
    a == 1 && return r + 2log(r - 1) - 2 / (r - 1) - log(4.0)
    kappa, rp, rm = _horizons(a)
    return r + 2rp / (rp - rm) * log((r - rp) / 2) - 2rm / (rp - rm) * log((r - rm) / 2)
end

@inline _at(v, k) = 0 <= k < length(v) ? v[k+1] : zero(C)

function _pmul(p, q)
    out = zeros(C, length(p) + length(q) - 1)
    for i in eachindex(p), j in eachindex(q)
        out[i+j-1] += p[i] * q[j]
    end
    return out
end

function _padd(ps...)
    out = zeros(C, maximum(length, ps))
    for p in ps
        out[1:length(p)] .+= p
    end
    return out
end

# Horizon Frobenius solution R = norm z^rho sum c_n z^n, z = r - r_+, at (complex) r; :in is Delta^{-s} e^{-i k r*}
# (rho = -s - i u), :out is e^{i k r*} (rho = i u), u = r_+ k / kappa. Returns (R e^{-L}, dR/dr e^{-L}, last term / sum, L) with L the real log magnitude.
function _horizon_solution(s, m, a, omega, lambda, r, branch::Symbol)
    kappa, rp, rm = _horizons(a)
    z = C(r) - rp
    k = omega - m * a / (2rp)
    u = rp * k / kappa
    K0 = 2rp * omega - a * m
    Kc = (C(K0), C(2rp * omega), C(omega))
    K2 = (Kc[1]^2, 2Kc[1] * Kc[2], Kc[2]^2 + 2Kc[1] * Kc[3], 2Kc[2] * Kc[3], Kc[3]^2)
    sK = (kappa * Kc[1], Kc[1] + kappa * Kc[2], Kc[2] + kappa * Kc[3], Kc[3])
    L = (C(lambda) - 4im * s * omega * rp, C(-4im * s * omega))
    zzL = (zero(C), 2kappa * L[1], 2kappa * L[2] + L[1], L[2])
    p2 = (C(4kappa^2), C(4kappa), one(C), zero(C), zero(C))
    p1 = (C(4(s + 1) * kappa^2), C(6(s + 1) * kappa), C(2(s + 1)), zero(C), zero(C))
    p0 = ntuple(j -> K2[j] - 2im * s * (j <= 4 ? sK[j] : zero(C)) - (j <= 4 ? zzL[j] : zero(C)), 5)
    rho = branch === :in ? C(-s - im * u) : C(im * u)
    # Run on the terms t_n = c_n z^n: near extremality (u ~ 1 / kappa) the c_n alone overflow before the terms decay.
    terms = C[one(C)]
    value = one(C)
    derivative = rho
    tail = Inf
    magnitude = 1.0
    for n in 1:4000
        acc = zero(C)
        zj = one(C)
        for j in 1:min(n, 4)
            zj *= z
            nn = (n - j) + rho
            acc += (p2[j+1] * nn * (nn - 1) + p1[j+1] * nn + p0[j+1]) * terms[n-j+1] * zj
        end
        nn = n + rho
        term = -acc / (p2[1] * nn * (nn - 1) + p1[1] * nn + p0[1])
        push!(terms, term)
        value += term
        derivative += nn * term
        magnitude += abs(term)
        tail = abs(term) / abs(value)
        n >= 4 && tail <= 1e-17 && break
    end
    derivative /= z
    lognorm = if branch === :in
        -s * log(2kappa) - im * k * rp + im * u * log(2.0) + im * k * rm / kappa * log(kappa)
    else
        im * k * rp - im * u * log(2.0) - im * k * rm / kappa * log(kappa)
    end
    # R = e^{L} (value, derivative) with the magnitude of e^{L} returned separately (it under- or overflows off the
    # real axis, |z^{-i u}| = e^{u arg z}).
    L = lognorm + rho * log(z)
    phase = cis(imag(L))
    return phase * value, phase * derivative, max(tail, magnitude / abs(value) * eps(Float64)), real(L)
end

# Exact extremal spin: r = 1 is an irregular singular point. R = e^{i beta / z} z^rho sum c_n z^n (z = r - 1) is an
# asymptotic series, summed to its smallest term; :in has beta = 2 omega - m (e^{-i k r*}) and rho = -2s - 2 i omega,
# :out has beta = -(2 omega - m) and rho = 2 i omega. Its normalization is not the package's, so only ratios of
# amplitudes built from one seed are used. Returns (R, dR/dr, smallest term / sum, log magnitude) as for _horizon_solution.
function _extremal_horizon_solution(s, m, omega, lambda, r, branch::Symbol)
    z = C(r) - 1
    K0 = 2omega - m
    beta = branch === :in ? C(K0) : C(-K0)
    p01 = 4omega * K0 - 2im * s * K0
    q1 = -2im * beta * s + p01
    q2 = C(4omega^2 + 2omega * K0 - lambda)
    q3 = C(4omega^2 + 2im * s * omega)
    q4 = C(omega^2)
    rho = q1 / (2im * beta)
    c = C[one(C)]
    value = one(C)
    dsum = rho / z
    zn = one(C)
    t1 = Inf
    t2 = Inf
    for n in 1:400
        nn = n - 1 + rho
        acc = (nn * (nn - 1) + 2(s + 1) * nn + q2) * c[n]
        n >= 2 && (acc += q3 * c[n-1])
        n >= 3 && (acc += q4 * c[n-2])
        cn = acc / (2im * beta * n)
        zn *= z
        term = cn * zn
        # Stop where the terms start to grow, as in _asymptotic_sum.
        (isfinite(term) && !(n >= 2 && abs(term) > t1 && abs(term) > t2)) || break
        push!(c, cn)
        value += term
        dsum += (n + rho) * term / z
        t2, t1 = t1, abs(term)
        max(t1, t2) <= 1e-18 * abs(value) && break
    end
    L = im * beta / z + rho * log(z)
    phase = cis(imag(L))
    return phase * value, phase * (dsum - im * beta / z^2 * value), min(t1, t2) / abs(value), real(L)
end

# Infinity series: R_in = e^{-i omega r*} sum c_n r^{-1-n}, R_out = e^{i omega r*} r^{-2s-1} sum d_n (omega r)^{-n}.
function _infinity_in_coefficients(s, m, a, w, lambda, order)
    D = C[a^2, -2, 1]
    rho = C[a^2, 0, 1]
    P2 = _pmul(D, D)
    P1 = _pmul(D, _padd(-2im * w .* rho, C[-2(s + 1), 2(s + 1)]))
    P0 = _padd(_pmul(D, C[4im * s * w - lambda, -2im * w]), C[0, 8im * s * w, -8im * s * w], -2a * m * w .* rho,
        C[a^2 * m^2 - 2im * s * a * m, 2im * s * a * m])
    c = zeros(C, order + 1)
    c[1] = 1
    for n in 1:order
        acc = zero(C)
        for d in 1:5
            j = n - d
            j < 0 && break
            acc += (_at(P2, 5 - d) * (j + 1) * (j + 2) - _at(P1, 4 - d) * (j + 1) + _at(P0, 3 - d)) * c[j+1]
        end
        c[n+1] = -acc / (2im * w * n)
    end
    return c
end

function _infinity_out_coefficients(s, m, a, w, lambda, order)
    a2 = a * a
    a4 = a2 * a2
    L2 = C.((0, 0, 0, 0, 1, -4, 2(a2 + 2), -4a2, a4))
    L1 = C.((0, 0, -2im * w, 2(s + 2im * w + 1), -2(2im * a2 * w + 5s + 5), 2(3a2 * s + 2im * a2 * w + 3a2 + 6s + 6),
        -2a2 * (im * a2 * w + 7s + 7), 4a4 * (s + 1), 0))
    L0 = C.((0, 0, -lambda - 2a * m * w, -2(-lambda + 2im * a2 * s * w + im * a2 * w - im * a * m * s + 2s * s + 3s + 1),
        -lambda * a2 - 2a * a2 * m * w + a2 * m * m + 4a2 * s * s + 8im * a2 * s * w + 6a2 * s + 4im * a2 * w + 2a2 -
            2im * a * m * s + 8s * s + 12s + 4,
        -2a2 * (2s + 1) * (im * a2 * w + 3s + 3), 2a4 * (s + 1) * (2s + 1), 0, 0))
    d = zeros(C, order + 1)
    d[1] = 1
    for n in 1:order
        acc = zero(C)
        for j in max(0, n - 8):n-1
            acc += (_at(L2, n - j + 3) * j * (j - 1) + _at(L1, n - j + 2) * j + _at(L0, n - j + 1)) * w^(n - j - 1) * d[j+1]
        end
        d[n+1] = acc / (2im * n)
    end
    return d
end

# sum c_n u^n and its u-derivative, stopped where the terms start to grow (a term above both of its two
# predecessors: a single small term, as d_1 ~ lambda ~ 1e-4 at s = l = 2, a = 1e-4, is not the optimal truncation);
# (h, dh/du, last kept term / |h|).
function _asymptotic_sum(c, u)
    h = zero(C)
    hu = zero(C)
    t1 = Inf
    t2 = Inf
    p = one(C)
    for n in 0:length(c)-1
        t = c[n+1] * p
        isfinite(t) || break
        n >= 2 && abs(t) > t1 && abs(t) > t2 && break
        h += t
        n > 0 && (hu += n * c[n+1] * p / u)
        t2, t1 = t1, abs(t)
        n >= 2 && max(t1, t2) < 1e-18 * abs(h) && break
        p *= u
    end
    return h, hu, max(t1, t2) / abs(h)
end

function _infinity_solutions(s, a, w, r, cin, cout)
    delta = r^2 - 2r + a^2
    phase = _rstar(a, r)
    h, hu, e_in = _asymptotic_sum(cin, 1 / r)
    g = h / r
    gp = -hu / r^3 - h / r^2
    Rin = cis(-w * phase) * g
    dRin = cis(-w * phase) * (gp - im * w * (r^2 + a^2) / delta * g)
    h, hu, e_out = _asymptotic_sum(cout, 1 / (w * r))
    p = -2s - 1
    g = r^p * h
    gp = p * r^(p - 1) * h - r^(p - 2) * hu / w
    Rout = cis(w * phase) * g
    dRout = cis(w * phase) * (gp + im * w * (r^2 + a^2) / delta * g)
    return Rin, dRin, Rout, dRout, max(e_in, e_out)
end

# One Taylor step of Delta^2 R'' + 2(s + 1)(r - 1) Delta R' + P0 R = 0 from c to c + h, with
# P0 = K^2 - 2 i s (r - 1) K - Delta (lambda - 4 i s omega r). Returns (R, dR/dr, last term / R).
function _taylor_step(s, m, a, w, lambda, R, dR, c, h, coefficients)
    # Local polynomials in t = r - c, as tuples (no allocation per step).
    d0, d1 = c^2 - 2c + a^2, 2c - 2
    k0, k1, k2 = (c^2 + a^2) * w - a * m, 2c * w, C(w)
    r0 = c - 1
    p2 = (d0 * d0, 2d0 * d1, d1 * d1 + 2d0, 2d1, one(C))
    q = 2(s + 1)
    p1 = (q * r0 * d0, q * (r0 * d1 + d0), q * (r0 + d1), C(q))
    KK = (k0 * k0, 2k0 * k1, k1 * k1 + 2k0 * k2, 2k1 * k2, k2 * k2)
    RK = (r0 * k0, r0 * k1 + k0, r0 * k2 + k1, k2)
    cD = (c * d0, c * d1 + d0, c + d1, one(C))
    is2 = 2im * s
    iw4 = 4im * s * w
    p0 = (KK[1] - is2 * RK[1] - lambda * d0 + iw4 * cD[1],
        KK[2] - is2 * RK[2] - lambda * d1 + iw4 * cD[2],
        KK[3] - is2 * RK[3] - lambda + iw4 * cD[3],
        KK[4] - is2 * RK[4] + iw4 * cD[4],
        KK[5])
    # Recur on b_n = c_n h^n: near the horizon the c_n alone overflow (|u / z|^n at order 40) while b_n decay.
    N = length(coefficients) - 1
    b = coefficients
    b[1] = R
    b[2] = dR * h
    hk = (one(C), h, h^2, h^3, h^4, h^5, h^6)
    for n in 0:N-2
        acc = zero(C)
        for k in 1:min(n + 2, 4)
            acc += p2[k+1] * hk[k+1] * (n - k + 2) * (n - k + 1) * b[n-k+3]
        end
        for k in 0:min(n + 1, 3)
            acc += _at(p1, k) * hk[k+2] * (n - k + 1) * b[n-k+2]
        end
        for k in 0:min(n, 4)
            acc += _at(p0, k) * hk[k+3] * b[n-k+1]
        end
        b[n+3] = -acc / (p2[1] * (n + 2) * (n + 1))
    end
    v = zero(C)
    d = zero(C)
    for n in 0:N
        v += b[n+1]
        n > 0 && (d += n * b[n+1])
    end
    return v, d / h, abs(b[N+1]) / abs(v)
end

# March (R, R') through the vertices of a path; steps are kept to 0.4 of the distance to r_+- and to one
# local wavelength of the tortoise phase, 1 / (omega |dr*/dr|). Returns (R, R', last Taylor term, log scale): the
# solution is (R, R') e^{log scale}.
function _march(s, m, a, w, lambda, R, dR, path; step_wavelengths = STEP_WAVELENGTHS)
    kappa, rp, rm = _horizons(a)
    coefficients = zeros(C, TAYLOR_ORDER + 1)
    worst = 0.0
    logscale = 0.0
    for q in 1:length(path)-1
        c = C(path[q])
        target = C(path[q+1])
        while c != target
            dist = min(abs(c - rp), abs(c - rm))
            hmax = min(0.4dist, step_wavelengths / (abs(w) * abs((c^2 + a^2) / (c^2 - 2c + a^2))))
            remaining = target - c
            last = abs(remaining) <= hmax
            h = last ? remaining : hmax * remaining / abs(remaining)
            R, dR, t = _taylor_step(s, m, a, w, lambda, R, dR, c, h, coefficients)
            worst = max(worst, t)
            c = last ? target : c + h
            # Rescale: near extremality |R| changes by e^{u |arg z|} (e^{2000} at a = 1 - 1e-6) along the path.
            size = abs(R)
            if !(1e-100 < size < 1e100)
                R /= size
                dR /= size
                logscale += log(size)
            end
        end
    end
    return R, dR, worst, logscale
end

# Projection of (R, R') on the two infinity solutions at the first radius in PROJECTION_RADII (shifted by
# -i depth) where both series reach SERIES_TOL.
function _projection_point(s, a, w, depth, cin, cout)
    for X in PROJECTION_RADII
        r = C(X, -depth)
        sols = _infinity_solutions(s, a, w, r, cin, cout)
        sols[5] <= SERIES_TOL && return r, sols
    end
    return nothing, nothing
end

_wronskian_split(R, dR, A, dA, B, dB) = (det = A * dB - B * dA; ((R * dB - B * dR) / det, (A * dR - R * dA) / det))

_seed(s, m, a, w, lambda, r, branch) = a == 1 ? _extremal_horizon_solution(s, m, w, lambda, r, branch) :
    _horizon_solution(s, m, a, w, lambda, r, branch)

"""
    complex_path_amplitudes(s, m, a, omega, lambda, branch; theta = SEED_THETA, depth = PATH_DEPTH, order = 400)

Teukolsky amplitudes of the IN (`branch = :IN`: B_inc, B_ref) or UP (`:UP`: C_up, C_ref) solution for real
`omega > 0` and `0 <= a <= 1`. The incident amplitude comes from the real axis, the reflected one from the complex
path (see the module notes). Returns `ratio` = reflection / incidence, and, for `a < 1`, `incidence` and
`reflection` as ratios to the unit transmission amplitude (at `a = 1` the seed is not normalized and only `ratio` is
meaningful), with `diagnostics` (last Taylor term, series truncation, projection radius, seed radius).
"""
function complex_path_amplitudes(s::Int, m::Int, a::Real, omega::Real, lambda::Real, branch::Symbol;
        theta::Real = SEED_THETA, depth::Real = PATH_DEPTH, order::Int = 400, seed_distance = nothing, pass_height = nothing, path = :left,
        reflection_only::Bool = false, step_wavelengths::Real = STEP_WAVELENGTHS)
    a = Float64(a); w = Float64(omega); lambda = Float64(lambda)
    (0 <= a <= 1) || throw(ArgumentError("complex_path_amplitudes needs 0 <= a <= 1."))
    w > 0 || throw(ArgumentError("complex_path_amplitudes needs real omega > 0."))
    branch in (:IN, :UP) || throw(ArgumentError("branch must be :IN or :UP."))
    extremal = a == 1
    extremal && 2w == m && throw(DomainError(w, "omega = m / 2 at a = 1: the horizon solutions are degenerate."))
    extremal && reflection_only && throw(ArgumentError("reflection_only needs a < 1 (the a = 1 seed is not normalized)."))
    kappa, rp, rm = _horizons(a)
    # Seed distance from r_+: kappa (inside the Frobenius radius 2 kappa), or at a = 1 a distance where the asymptotic
    # horizon series still reaches 1e-17 (|2 omega - m| / |z| >= 40).
    delta = seed_distance !== nothing ? Float64(seed_distance) : extremal ? min(0.3, abs(2w - m) / 40) : kappa
    cin = _infinity_in_coefficients(s, m, a, w, lambda, order)
    cout = _infinity_out_coefficients(s, m, a, w, lambda, order)
    r_real, real_sols = reflection_only ? (C(PROJECTION_RADII[1]), nothing) : _projection_point(s, a, w, 0.0, cin, cout)
    # IN runs below the axis; UP, read off at the horizon, runs along the mirror path above it.
    side = branch === :IN ? 1.0 : -1.0
    r_cplx, cplx_sols = _projection_point(s, a, w, side * Float64(depth), cin, cout)
    (r_real === nothing || r_cplx === nothing) && error("infinity series did not reach $(SERIES_TOL) by r = $(last(PROJECTION_RADII)).")
    # Near extremality (u = r_+ k / kappa large) the Frobenius terms grow far beyond the sum before they decay; the
    # seed moves in until the summed term sizes stay within SEED_CONDITION of the value (the error reported as
    # condition * eps), at both seeds and for both horizon solutions.
    if !extremal && seed_distance === nothing
        resolved = false
        while delta > 1e-3 * kappa
            worst = maximum(_horizon_solution(s, m, a, w, lambda, r, br)[3]
                for r in (rp + delta, rp + delta * cis(-side * theta)), br in (:in, :out))
            (resolved = worst <= SEED_CONDITION * eps(Float64)) && break
            delta /= 2
        end
        # 1 - a below ~1e-7 (u = r_+ k / kappa ~ 1e4): no seed distance resolves the series; callers fall back.
        resolved || throw(DomainError(a, "complex_path_amplitudes: horizon series unresolved this close to extremality."))
    end
    r_seed_real = C(rp + delta)
    r_seed = rp + delta * cis(-side * theta)
    path_real = (r_seed_real, r_real)
    # The path passes under r_- before it descends: straight down from the seed, Im r* rises to about -5 near
    # r ~ 1 - 2 (r_- enters r* with the opposite sign), and the roundoff there sets the floor of B_ref; left of r_-
    # both logarithms add up to -2 pi. At (s, l, m, a, omega) = (-2, 2, 2, 0.5, 10), B_ref / B_inc = 4e-63: 2e7 -> 1e-7.
    # Near extremality (kappa < NEAR_EXTREMAL_KAPPA) the path passes r_- at PASS_KAPPAS kappa: closer, the two r_-
    # solutions differ by e^{pi u_-} along the half turn and the march loses the recessive one (a = 0.9999, omega = 10:
    # 1e83 at the seed height, 3e-8 at 5 kappa). Away from extremality the seed height is best (a = 0, omega = 10:
    # 1e-8, against 5e3 at height 0.5).
    lift = pass_height !== nothing ? Float64(pass_height) :
        kappa < NEAR_EXTREMAL_KAPPA ? max(abs(imag(r_seed)), PASS_KAPPAS * kappa, extremal ? delta : 0.0) : abs(imag(r_seed))
    path_cplx = path === :down ? (r_seed, r_seed - im * side * depth, r_cplx) :
        (r_seed, C(real(r_seed), -side * max(lift, abs(imag(r_seed)))), C(rm - PATH_LEFT, -side * lift),
            C(rm - PATH_LEFT, -side * depth), r_cplx)
    if reflection_only
        # The seed is normalized to unit transmission: the reflection alone, without the real-axis march.
        if branch === :IN
            R, dR, _, L_cplx = _seed(s, m, a, w, lambda, r_seed, :in)
            R, dR, t_cplx, M = _march(s, m, a, w, lambda, R, dR, path_cplx; step_wavelengths)
            _, reflection_scaled = _wronskian_split(R, dR, cplx_sols[1:4]...)
            reflection = reflection_scaled * exp(L_cplx + M)
        else
            _, _, R, dR, _ = cplx_sols
            R, dR, t_cplx, M_cplx = _march(s, m, a, w, lambda, R, dR, reverse(path_cplx); step_wavelengths)
            A, dA, _, _ = _seed(s, m, a, w, lambda, r_seed, :out)
            B, dB, _, LB = _seed(s, m, a, w, lambda, r_seed, :in)
            _, reflection_scaled = _wronskian_split(R, dR, A, dA, B, dB)
            reflection = reflection_scaled * exp(M_cplx - LB)
        end
        return (ratio = NaN + NaN * im, incidence = NaN + NaN * im, reflection = reflection,
            diagnostics = (taylor_tail = t_cplx, series = cplx_sols[5], projection_real = NaN,
                projection_complex = r_cplx, seed_distance = delta))
    end
    if branch === :IN
        R, dR, _, L_real = _seed(s, m, a, w, lambda, r_seed_real, :in)
        R, dR, t_real, M = _march(s, m, a, w, lambda, R, dR, path_real; step_wavelengths)
        L_real += M
        incidence_scaled, _ = _wronskian_split(R, dR, real_sols[1:4]...)
        R, dR, _, L_cplx = _seed(s, m, a, w, lambda, r_seed, :in)
        R, dR, t_cplx, M = _march(s, m, a, w, lambda, R, dR, path_cplx; step_wavelengths)
        L_cplx += M
        _, reflection_scaled = _wronskian_split(R, dR, cplx_sols[1:4]...)
        ratio = reflection_scaled / incidence_scaled * exp(L_cplx - L_real)
        incidence = incidence_scaled * exp(L_real)
    else
        # UP: start from the outgoing series at infinity and march inward; split on the horizon solutions
        # (:out = e^{i k r*}, coefficient C_up; :in = Delta^{-s} e^{-i k r*}, coefficient C_ref).
        _, _, R, dR, _ = real_sols
        R, dR, t_real, M_real = _march(s, m, a, w, lambda, R, dR, reverse(path_real); step_wavelengths)
        A, dA, _, LA = _seed(s, m, a, w, lambda, r_seed_real, :out)
        B, dB, _, _ = _seed(s, m, a, w, lambda, r_seed_real, :in)
        incidence_scaled, _ = _wronskian_split(R, dR, A, dA, B, dB)
        _, _, R, dR, _ = cplx_sols
        R, dR, t_cplx, M_cplx = _march(s, m, a, w, lambda, R, dR, reverse(path_cplx); step_wavelengths)
        A, dA, _, _ = _seed(s, m, a, w, lambda, r_seed, :out)
        B, dB, _, LB = _seed(s, m, a, w, lambda, r_seed, :in)
        _, reflection_scaled = _wronskian_split(R, dR, A, dA, B, dB)
        ratio = reflection_scaled / incidence_scaled * exp(LA - LB + M_cplx - M_real)
        incidence = incidence_scaled * exp(M_real - LA)
    end
    reflection = ratio * incidence
    return (ratio = ratio, incidence = extremal ? NaN + NaN * im : incidence,
        reflection = extremal ? NaN + NaN * im : reflection,
        diagnostics = (taylor_tail = max(t_real, t_cplx), series = max(real_sols[5], cplx_sols[5]),
            projection_real = real(r_real), projection_complex = r_cplx, seed_distance = delta))
end

end
