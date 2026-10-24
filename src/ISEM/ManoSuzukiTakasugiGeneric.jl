# Type-generic Mano-Suzuki-Takasugi (MST) representation of the Teukolsky UP and DOWN solutions.
#
# The series of Sasaki & Tagoshi, Living Rev. Relativ. 6, 6 (2003), Sec. 4.4 (R_-^nu, the outgoing
# "up" solution, and R_+^nu, the ingoing "down" solution at infinity) are evaluated in the floating
# type R of the mode parameters, with the conventions of the Black Hole Perturbation Toolkit
# Teukolsky paclet (Kernel/MST/MST.m) so that the unit-transmission UP solution is directly comparable
# with TeukolskyRadial[..., "BoundaryConditions" -> "Up"]:
#
#   R_up(r)   = prefacUp(zhat) sum_n fUp_n (-2 i zhat)^n U(n + nu + s + 1 - i eps, 2n + 2nu + 2, -2 i zhat) / UpTrans,
#   fUp_n     = (-1)^n (nu + 1 + s - i eps)_n / (nu + 1 - s + i eps)_n f_n,
#   UpTrans   = (eps/2)^(-1-2s) exp(i eps (log eps - (1 - kappa)/2)) 2^(-s-1+i eps) exp(-pi eps/2 - i pi (nu+1+s)/2) sum_n fUp_n,
#
# with zhat = omega (r - r_-), eps = 2 omega, tau = (eps - m a)/kappa, f_n the MST series coefficients
# (f_0 = 1) and nu the renormalized angular momentum. The DOWN solution (ingoing at infinity) is the
# companion Coulomb-type series
#
#   R_down(r) = e^{-i zhat} zhat^{nu + i eps_+} (zhat - eps kappa)^{-s - i eps_+} sum_n f_n (2 i zhat)^n U(n + nu + 1 - s + i eps, 2n + 2nu + 2, 2 i zhat) / DownTrans,
#   DownTrans = (sum_n f_n) (2 i)^{-nu-1+s-i eps} omega^{-1-i eps} e^{i omega r_-} e^{-i eps log 2},
#
# eps_+ = (eps + tau)/2, normalized to a unit ingoing wave r^-1 e^{-i omega r*}; its form was fixed by
# comparison with the asymptotic-series DOWN solution at real frequency (constant ratio at every r),
# and DownTrans is the leading large-r coefficient of that series. Both solutions converge at every
# finite r > r_+, so they are free of the Stokes ambiguity of the asymptotic series at complex
# frequencies (see ArbitraryPrecisionRoute.jl).
#
# The confluent hypergeometric U is built from Kummer's M by the connection formula at a working
# precision raised by the cancellation of the Kummer series (about |z| / ln 2 bits) and by the
# closeness of the second parameter to an integer; its Gamma factors come from `_loggamma_big`.
module ManoSuzukiTakasugiGeneric

using SpecialFunctions: loggamma
using ..Parameters: GSNParameters

export MSTGenericData, mst_generic_data
export mst_nu_generic, mst_up_state, mst_down_state, mst_up_transmission, mst_down_transmission

const NU_CF_DEPTH_INITIAL = 40
const NU_CF_DEPTH_MAX = 20000
const NU_MAX_ITERATIONS = 100
const COEFFICIENT_CAP = 2000
const SERIES_MIN_TERMS = 4
const COEFFICIENT_MIN_TERMS = 32
const COEFFICIENT_EXTENSION = 64
const U_GUARD_BITS = 64
const KUMMER_MAX_TERMS = 200000

_real_type(::GSNParameters{R}) where {R} = R

# pi in the real type of a complex working value (a bare pi next to im is rounded to Float64).
_pi(z) = convert(real(typeof(z)), pi)
_complex_type(::GSNParameters{R}) where {R} = Complex{R}

"""
    MSTGenericData

MST data of one mode in the floating type of its parameters: the renormalized angular momentum `nu`,
the series coefficients `f_n` for `nmin <= n <= nmax` (`coefficient(data, n)`), the tolerance of the
sums, and the derived quantities eps = 2 omega and tau = (eps - m a)/kappa.
"""
mutable struct MSTGenericData{P<:GSNParameters,R,C}
    const params::P
    const epsilon::C
    const tau::C
    const nu::C
    coefficients::Vector{C}
    up_factors::Vector{C}       # (-1)^n (nu+1+s-i eps)_n / (nu+1-s+i eps)_n, aligned with coefficients
    nmin::Int
    nmax::Int
    const tolerance::R
    const nu_residual::R
    const nu_depth::Int
    up_normalization::Union{Nothing,C}      # cached mst_up_transmission
    down_normalization::Union{Nothing,C}    # cached mst_down_transmission
end

coefficient(data::MSTGenericData, n::Integer) = data.coefficients[n - data.nmin + 1]

# The Pochhammer ratios of the UP series, (-1)^n (x)_n / (y)_n with x = nu+1+s-i eps, y = nu+1-s+i eps,
# for nmin <= n <= nmax by the recurrence g(n)/g(n-1) = (x+n-1)/(y+n-1) in both directions from n = 0
# (two complex loggammas per term otherwise, the whole cost of an UP evaluation).
function _up_factors(p::GSNParameters, nu, eps_, nmin::Int, nmax::Int)
    C = typeof(nu)
    x = nu + 1 + p.s - im * eps_
    y = nu + 1 - p.s + im * eps_
    factors = Vector{C}(undef, nmax - nmin + 1)
    factors[1 - nmin] = one(C)
    g = one(C)
    for n in 1:nmax
        g *= -(x + n - 1) / (y + n - 1)
        factors[n - nmin + 1] = g
    end
    g = one(C)
    for n in -1:-1:nmin
        g *= -(y + n) / (x + n)
        factors[n - nmin + 1] = g
    end
    return factors
end

# Extend the stored coefficients to |n| <= N (the Tricomi sums at small zhat need more terms than the
# decay of f_n alone suggests, because the U factors grow with n); the ratios are recomputed from a
# deeper continued fraction, which also refreshes the stored ones.
function _extend_coefficients!(data::MSTGenericData, N::Int)
    N <= max(-data.nmin, data.nmax) && return data
    N <= COEFFICIENT_CAP || error("MST generic: more than $COEFFICIENT_CAP series coefficients would be needed.")
    p = data.params
    C = typeof(data.nu)
    ratios_R, ratios_L = _ratios(p, data.epsilon, p.kappa, data.tau, data.nu, N, 60)
    positive = Vector{C}(undef, N)
    negative = Vector{C}(undef, N)
    f = one(C)
    for n in 1:N
        f *= ratios_R[n]
        positive[n] = f
    end
    f = one(C)
    for n in 1:N
        f *= ratios_L[n]
        negative[n] = f
    end
    data.coefficients = vcat(reverse(negative), [one(C)], positive)
    data.nmin = -N
    data.nmax = N
    data.up_factors = _up_factors(p, data.nu, data.epsilon, -N, N)
    data.up_normalization = nothing
    data.down_normalization = nothing
    return data
end

# ---------------------------------------------------------------------------------------------
# Three-term recurrence alpha_n f_(n+1) + beta_n f_n + gamma_n f_(n-1) = 0 (ST Eq. 124)
# ---------------------------------------------------------------------------------------------
@inline function _alpha(p::GSNParameters, eps, kappa, tau, nu, n)
    s = p.s
    npnu1 = n + nu + 1
    return im * eps * kappa * (npnu1 + s + im * eps) * (npnu1 + s - im * eps) * (npnu1 + im * tau) /
        (npnu1 * (2 * npnu1 + 1))
end
@inline function _beta(p::GSNParameters, eps, kappa, tau, nu, n)
    s = p.s
    product = (n + nu) * (n + nu + 1)
    linear = eps * (eps - p.m * p.a)
    return -p.lambda - s * (s + 1) + product + eps^2 + linear + linear * (s^2 + eps^2) / product
end
@inline function _gamma(p::GSNParameters, eps, kappa, tau, nu, n)
    s = p.s
    npnu = n + nu
    return -im * eps * kappa * (npnu - s + im * eps) * (npnu - s - im * eps) * (npnu - im * tau) /
        (npnu * (2 * npnu - 1))
end

# R_n = f_n / f_(n-1) for n = 1..N by backward recursion from depth N + extra (ST Eq. 125), and
# L_n = f_n / f_(n+1) for n = -1..-N by the mirrored recursion. Returned as vectors R[n], L[-n].
function _ratios(p::GSNParameters, eps, kappa, tau, nu, N::Int, depth::Int)
    C = typeof(nu)
    top = N + depth
    R = Vector{C}(undef, N)
    r = -_gamma(p, eps, kappa, tau, nu, top) / _beta(p, eps, kappa, tau, nu, top)
    for k in (top - 1):-1:1
        r = -_gamma(p, eps, kappa, tau, nu, k) / (_beta(p, eps, kappa, tau, nu, k) + _alpha(p, eps, kappa, tau, nu, k) * r)
        k <= N && (R[k] = r)
    end
    L = Vector{C}(undef, N)
    l = -_alpha(p, eps, kappa, tau, nu, -top) / _beta(p, eps, kappa, tau, nu, -top)
    for k in (-top + 1):-1
        l = -_alpha(p, eps, kappa, tau, nu, k) / (_beta(p, eps, kappa, tau, nu, k) + _gamma(p, eps, kappa, tau, nu, k) * l)
        k >= -N && (L[-k] = l)
    end
    return R, L
end

# Characteristic function whose zero is nu: beta_0 + alpha_0 R_1 + gamma_0 L_(-1) (continued fractions
# of depth `depth`).
function _characteristic(p::GSNParameters, eps, kappa, tau, nu, depth::Int)
    r = -_gamma(p, eps, kappa, tau, nu, depth) / _beta(p, eps, kappa, tau, nu, depth)
    for k in (depth - 1):-1:1
        r = -_gamma(p, eps, kappa, tau, nu, k) / (_beta(p, eps, kappa, tau, nu, k) + _alpha(p, eps, kappa, tau, nu, k) * r)
    end
    l = -_alpha(p, eps, kappa, tau, nu, -depth) / _beta(p, eps, kappa, tau, nu, -depth)
    for k in (-depth + 1):-1
        l = -_alpha(p, eps, kappa, tau, nu, k) / (_beta(p, eps, kappa, tau, nu, k) + _gamma(p, eps, kappa, tau, nu, k) * l)
    end
    return _beta(p, eps, kappa, tau, nu, 0) + _alpha(p, eps, kappa, tau, nu, 0) * r + _gamma(p, eps, kappa, tau, nu, 0) * l
end

"""
    mst_nu_generic(params, seed; tol)

Renormalized angular momentum in the floating type of `params`, by secant iteration on the
continued-fraction characteristic equation from `seed` (a binary64 value, e.g. from the Float64 MST
module). The continued-fraction depth is doubled until the characteristic function at the seed is
stable to `tol`. Returns `(nu, residual, depth)`.
"""
function mst_nu_generic(p::GSNParameters, seed; tol=nothing)
    R = _real_type(p)
    C = _complex_type(p)
    tolerance = tol === nothing ? 64 * eps(R) : R(tol)
    eps_ = C(2 * p.omega)
    kappa = p.kappa
    tau = (eps_ - p.m * p.a) / kappa
    nu0 = C(seed)
    depth = NU_CF_DEPTH_INITIAL
    g0 = _characteristic(p, eps_, kappa, tau, nu0, depth)
    while depth < NU_CF_DEPTH_MAX
        g1 = _characteristic(p, eps_, kappa, tau, nu0, 2depth)
        abs(g1 - g0) <= tolerance * max(one(R), abs(g1)) && break
        depth *= 2
        g0 = g1
    end
    depth = 2depth
    scale = max(one(R), abs(p.lambda), R(p.l * (p.l + 1)))
    x0 = nu0
    x1 = nu0 + C(1e-6)
    f0 = _characteristic(p, eps_, kappa, tau, x0, depth)
    f1 = _characteristic(p, eps_, kappa, tau, x1, depth)
    for _ in 1:NU_MAX_ITERATIONS
        f1 == f0 && break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        f2 = _characteristic(p, eps_, kappa, tau, x2, depth)
        x0, f0, x1, f1 = x1, f1, x2, f2
        if abs(x1 - x0) <= tolerance * max(one(R), abs(x1)) && abs(f1) <= sqrt(tolerance) * scale
            # alpha_n and gamma_n have poles where n + nu = -1/2 or 1/2, and the characteristic function
            # tends to zero there: at omega = 2 - 6i (a = 0.7, l = m = 2) every seed converged to a
            # half-integer, whose MST sums then gave amplitudes that changed by orders of magnitude with
            # the precision. Such a limit is not a renormalized angular momentum.
            abs(2 * x1 - round(real(2 * x1))) <= sqrt(tolerance) &&
                error("MST generic: the iteration from the seed $seed converged to the half-integer pole nu = $(ComplexF64(x1)) of the continued fraction, not to a root.")
            return x1, abs(f1), depth
        end
    end
    error("MST generic: the renormalized angular momentum iteration did not converge from the seed $seed (residual $(abs(f1))).")
end

"""
    mst_generic_data(params, nu; tol, N)

MST coefficients f_n for |n| <= N (default: until |f_n| falls below tol relative to f_0 by a wide
margin, at most $(COEFFICIENT_CAP)) from the continued-fraction ratios, for the given renormalized
angular momentum `nu` (from `mst_nu_generic`).
"""
function mst_generic_data(p::GSNParameters, nu; tol=nothing, N=nothing, nu_residual=zero(_real_type(p)), nu_depth=0)
    R = _real_type(p)
    C = _complex_type(p)
    tolerance = tol === nothing ? 64 * eps(R) : R(tol)
    eps_ = C(2 * p.omega)
    kappa = p.kappa
    tau = (eps_ - p.m * p.a) / kappa
    nu_c = C(nu)
    # f_n by products of the continued-fraction ratios; each side stops once the coefficients are far
    # below the tolerance (the U factors and powers of 2 zhat can enhance a term by many orders of
    # magnitude, hence the margin). The ratios come from a backward recursion that must start above the
    # last term wanted, so the count is doubled from a small start until both sides have converged
    # (computing them for the cap at once cost more than the rest of the MST data).
    margin = tolerance * ldexp(one(R), -40)
    local positive, negative
    count = N === nothing ? 2 * COEFFICIENT_MIN_TERMS : Int(N)
    while true
        ratios_R, ratios_L = _ratios(p, eps_, kappa, tau, nu_c, count, 60)
        positive = C[]
        f = one(C)
        for n in 1:count
            f *= ratios_R[n]
            push!(positive, f)
            N === nothing && abs(f) < margin && n >= COEFFICIENT_MIN_TERMS && break
        end
        negative = C[]
        f = one(C)
        for n in 1:count
            f *= ratios_L[n]
            push!(negative, f)
            N === nothing && abs(f) < margin && n >= COEFFICIENT_MIN_TERMS && break
        end
        (N !== nothing || (length(positive) < count && length(negative) < count) || count >= COEFFICIENT_CAP) && break
        count = min(2count, COEFFICIENT_CAP)
    end
    nmin = -length(negative)
    nmax = length(positive)
    coefficients = vcat(reverse(negative), [one(C)], positive)
    return MSTGenericData{typeof(p),R,C}(p, eps_, tau, nu_c, coefficients, _up_factors(p, nu_c, eps_, nmin, nmax),
        nmin, nmax, tolerance, R(nu_residual), nu_depth, nothing, nothing)
end

# ---------------------------------------------------------------------------------------------
# Confluent hypergeometric functions in the working type
# ---------------------------------------------------------------------------------------------
# Kummer M(a, b, z) by its power series (entire in z), summed at the current precision; also returns
# the sum of the term magnitudes, whose ratio to |M| is the cancellation of the series.
function _kummer_M(a::C, b::C, z::C) where {C<:Complex}
    T = real(C)
    term = one(C)
    total = one(C)
    absolute = one(T)
    threshold = eps(T) / 4
    for k in 0:KUMMER_MAX_TERMS
        term *= (a + k) / (b + k) * z / (k + 1)
        total += term
        absolute += abs(term)
        k > 2 * abs(z) && abs(term) <= threshold * abs(total) && return total, absolute
    end
    error("MST generic: the Kummer series did not converge for |z| = $(Float64(abs(z))).")
end

# Complex log-gamma at the current BigFloat precision, up to a multiple of 2 pi i (it enters the
# connection formula only through exp). SpecialFunctions 2.8 recomputes the Bernoulli numbers in rational
# arithmetic on every complex BigFloat call (quadratic in the precision, with growing integers): at the
# strongly damped and imaginary-axis frequencies, whose MST sums re-anchor their recurrences often, that
# was three quarters of the time of a 288-bit evaluation, and its Stirling truncation leaves about 7% of
# the bits (2^-1873 at 2000 bits). The same Stirling series here takes the coefficients B_2k / (2k (2k - 1))
# from a cache filled once per process, shifts the argument to |z| >= the precision in bits and sums
# until the terms are below 2^-(bits + 16): exp of the result agrees with a 300-bit-higher evaluation
# to the last bit at 128..2000 bits, at 1/7 of the time at 2000 bits.
const _BERNOULLI_NUMBERS = Rational{BigInt}[1 // 1, -1 // 2]      # B_0, B_1, B_2, ...
const _STIRLING_COEFFICIENTS = Rational{BigInt}[]
const _STIRLING_LOCK = ReentrantLock()

# B_m from sum_{k=0}^{m} binomial(m + 1, k) B_k = 0 (B_m = 0 for odd m >= 3); caller holds the lock.
function _bernoulli_number(m::Int)
    B = _BERNOULLI_NUMBERS
    while length(B) <= m
        n = length(B)
        if isodd(n)
            push!(B, 0 // 1)
            continue
        end
        total = zero(Rational{BigInt})
        binomial_ = one(BigInt)                 # binomial(n + 1, k)
        for k in 0:n-1
            iszero(B[k+1]) || (total += binomial_ * B[k+1])
            binomial_ = binomial_ * (n + 1 - k) ÷ (k + 1)
        end
        push!(B, -total // (n + 1))
    end
    return B[m+1]
end

function _stirling_coefficients(n::Int)
    lock(_STIRLING_LOCK) do
        while length(_STIRLING_COEFFICIENTS) < n
            k = length(_STIRLING_COEFFICIENTS) + 1
            push!(_STIRLING_COEFFICIENTS, _bernoulli_number(2k) // (2k * (2k - 1)))
        end
        return _STIRLING_COEFFICIENTS[1:n]
    end
end

function _loggamma_big(z::Complex{BigFloat})
    bits = precision(BigFloat)
    value = setprecision(BigFloat, bits + 16) do
        zz = Complex{BigFloat}(z)
        real(zz) < 1 // 2 && return log(big(pi) / sinpi(zz)) - _loggamma_big(1 - zz)     # reflection
        shift = max(0, ceil(Int, bits - Float64(abs(zz))))
        zr = zz + shift
        w = inv(zr)
        w2 = w * w
        threshold = ldexp(one(BigFloat), -(bits + 16))
        total = (zr - one(BigFloat) / 2) * log(zr) - zr + log(2 * big(pi)) / 2
        coefficients = _stirling_coefficients(max(16, bits ÷ 10))
        power = w
        k = 0
        while true
            k += 1
            # (|z + shift| >= bits reaches the threshold within about bits / 12 terms; a non-finite argument
            # never does)
            isfinite(total) && k <= bits || return loggamma(zz)
            k > length(coefficients) && (coefficients = _stirling_coefficients(min(2 * length(coefficients), bits)))
            term = BigFloat(coefficients[k]) * power
            total += term
            abs(term) <= threshold * (1 + abs(total)) && break
            power *= w2
        end
        product = one(Complex{BigFloat})
        for i in 0:shift-1
            product *= zz + i
        end
        total - log(product)
    end
    return Complex{BigFloat}(value)
end

# Distance of a complex number to the nearest integer, in the real type.
_integer_distance(b) = abs(b - round(real(b)))

# Tricomi U(a, b, z) for non-integer b by the connection formula
#   U = Gamma(1-b)/Gamma(a-b+1) M(a,b,z) + Gamma(b-1)/Gamma(a) z^(1-b) M(a-b+1,2-b,z),
# evaluated at a working precision raised by the cancellation of the Kummer series (|z| / ln 2 bits
# plus a half) and by the closeness of b to an integer, then rounded to the caller's precision.
# The cancellation of the connection formula (the two terms against their sum) and of the two Kummer
# series is measured, and the evaluation is repeated at a higher precision until it is covered with a
# margin of U_CANCELLATION_MARGIN_BITS.
const U_CANCELLATION_MARGIN_BITS = 32
const U_MAX_ATTEMPTS = 6

# Large argument: the asymptotic expansion U(a, b, z) ~ z^(-a) sum_k (a)_k (a-b+1)_k / (k! (-z)^k)
# (DLMF 13.7.3, |ph z| < 3 pi / 2), whose smallest term is about exp(-|z|) times the parameter growth.
# The Kummer series would need about |z| terms at |z| / ln 2 extra bits (four million terms at r = 1e6),
# so the expansion is tried first when |z| exceeds the bits of the precision; it is accepted when a
# term falls below the tolerance, and the Kummer route is used when the terms grow again past the
# smallest one (k > |z| / 2) first. While |(a + k)(a - b + 1 + k)| > (k + 1)|z| the terms grow before they
# decrease: the sums of the MST series reach |n| ~ 400 at r = 1000, where that holds, and refusing the
# expansion there sent single values to a Kummer series at 9000 extra bits (two minutes for one point
# of omega = 1e-6 - 3i). The growth costs log2(largest term / |sum|) bits, which the evaluation measures
# and covers by repeating it at a precision raised by as much.
const U_ASYMPTOTIC_MIN_Z_PER_BIT = 0.8
function _tricomi_U_asymptotic(a::Complex{BigFloat}, b::Complex{BigFloat}, z::Complex{BigFloat}, bits::Int)
    extra = U_GUARD_BITS
    for _ in 1:2
        value, lost = setprecision(BigFloat, bits + extra) do
            aa = Complex{BigFloat}(a); bb = Complex{BigFloat}(b); zz = Complex{BigFloat}(z)
            threshold = ldexp(one(BigFloat), -(bits + 8))
            c = aa - bb + 1
            divergence = Float64(abs(zz)) / 2
            term = one(Complex{BigFloat})
            total = term
            previous = abs(term)
            largest = previous
            for k in 0:200000
                term *= (aa + k) * (c + k) / ((k + 1) * (-zz))
                magnitude = abs(term)
                magnitude > previous && k > divergence && return nothing, 0   # diverging before the tolerance
                total += term
                largest = max(largest, magnitude)
                if magnitude <= threshold * abs(total)
                    return zz^(-aa) * total, Float64(log2(largest / abs(total)))
                end
                previous = magnitude
            end
            nothing, 0
        end
        (value === nothing || !isfinite(lost)) && return nothing
        lost <= extra - 16 && return Complex{BigFloat}(value)
        extra = ceil(Int, lost) + U_GUARD_BITS
    end
    return nothing
end

function _tricomi_U(a::Complex{BigFloat}, b::Complex{BigFloat}, z::Complex{BigFloat})
    bits = precision(BigFloat)
    # (valid for |ph z| < 3 pi / 2, so on the whole principal sheet; near the negative real axis the
    # recessive exp(z) part it leaves out is exp(-|z|) small. Excluding |ph z| > 0.95 pi sent the far points
    # of near-axis frequencies, where z = -2 i zhat is nearly negative real, to a Kummer series of 1e5 terms.)
    if Float64(abs(z)) >= U_ASYMPTOTIC_MIN_Z_PER_BIT * bits
        asymptotic = _tricomi_U_asymptotic(a, b, z, bits)
        asymptotic === nothing || return asymptotic
    end
    distance = Float64(_integer_distance(b))
    distance > 0 || error("MST generic: Tricomi U with an integer second parameter is not implemented.")
    extra = ceil(Int, 1.5 * Float64(abs(z))) + U_GUARD_BITS + ceil(Int, max(0.0, -log2(distance)))
    local value
    for _ in 1:U_MAX_ATTEMPTS
        value, cancellation = setprecision(BigFloat, bits + extra) do
            aa = Complex{BigFloat}(a); bb = Complex{BigFloat}(b); zz = Complex{BigFloat}(z)
            m1, abs1 = _kummer_M(aa, bb, zz)
            m2, abs2 = _kummer_M(aa - bb + 1, 2 - bb, zz)
            f1 = exp(_loggamma_big(1 - bb) - _loggamma_big(aa - bb + 1))
            f2 = exp(_loggamma_big(bb - 1) - _loggamma_big(aa) + (1 - bb) * log(zz))
            t1 = f1 * m1
            t2 = f2 * m2
            total = t1 + t2
            magnitude = abs(f1) * abs1 + abs(f2) * abs2       # every term of both series, in magnitude
            total, Float64(magnitude / max(abs(total), floatmin(BigFloat)))
        end
        lost = isfinite(cancellation) && cancellation > 0 ? ceil(Int, log2(cancellation)) : extra
        lost + U_CANCELLATION_MARGIN_BITS <= extra && break
        extra = lost + U_CANCELLATION_MARGIN_BITS + U_GUARD_BITS
    end
    return Complex{BigFloat}(value)
end
function _tricomi_U(a::C, b::C, z::C) where {C<:Complex}
    # Other floating types: evaluate through BigFloat at a precision matching the type plus the guard.
    T = real(C)
    bits = precision(T)
    value = setprecision(BigFloat, bits + U_GUARD_BITS) do
        _tricomi_U(Complex{BigFloat}(a), Complex{BigFloat}(b), Complex{BigFloat}(z))
    end
    return C(value)
end

# log of the Pochhammer symbol (x)_n = Gamma(x + n) / Gamma(x) for any integer n.
_log_pochhammer(x, n::Integer) = loggamma(x + n) - loggamma(x)

# ---------------------------------------------------------------------------------------------
# UP solution (ST Eq. 155 / BHPT prefacUp, fUp, UpTrans)
# ---------------------------------------------------------------------------------------------
_zhat(data::MSTGenericData, r) = data.params.omega * (r - (one(data.params.kappa) - data.params.kappa))

# fUp_n = (-1)^n (nu + 1 + s - i eps)_n / (nu + 1 - s + i eps)_n f_n (the ratio from the stored recurrence)
_fup(data::MSTGenericData, n::Integer) = data.up_factors[n - data.nmin + 1] * coefficient(data, n)

"""
    mst_up_transmission(data)

The amplitude UpTrans of the raw UP series (BHPT prefacUpTrans * A_-), i.e. the Teukolsky
transmission amplitude of the raw series; dividing by it gives the unit-transmission UP solution
R_up -> r^(-2s-1) exp(i omega r*).
"""
function mst_up_transmission(data::MSTGenericData)
    data.up_normalization === nothing || return data.up_normalization
    data.up_normalization = _mst_up_transmission(data)
    return data.up_normalization
end
function _mst_up_transmission(data::MSTGenericData)
    p = data.params
    e = data.epsilon
    k = p.kappa
    nu = data.nu
    C = typeof(nu)
    total = zero(C)
    for n in data.nmin:data.nmax
        total += _fup(data, n)
    end
    pi_ = _pi(nu)
    aminus = 2^(-p.s - 1 + im * e) * exp(-pi_ * e / 2 - im * pi_ * (nu + 1 + p.s) / 2) * total
    prefactor = (e / 2)^(-1 - 2 * p.s) * exp(im * e * (log(e) - (1 - k) / 2))
    return prefactor * aminus
end

# prefacUp(zhat) and d log prefacUp / d zhat
function _up_prefactor(data::MSTGenericData, zhat)
    p = data.params
    e = data.epsilon
    k = p.kappa
    nu = data.nu
    eplus = (e + data.tau) / 2
    pi_ = _pi(nu)
    value = 2^nu * exp(-pi_ * e) * exp(-im * pi_ * (nu + 1 + p.s)) * exp(im * zhat) *
        zhat^(nu + im * eplus) * (zhat - e * k)^(-p.s - im * eplus)
    dlog = im + (nu + im * eplus) / zhat + (-p.s - im * eplus) / (zhat - e * k)
    return value, dlog
end

# h_n = c^n U(a + n, b + 2n, c) and dh_n/dc by direct evaluation (two Tricomi functions).
function _h_direct(a, b, c, n::Integer)
    U0 = _tricomi_U(a + n, b + 2n, c)
    U1 = _tricomi_U(a + n + 1, b + 2n + 1, c)
    cn = c^n
    return cn * U0, n * cn / c * U0 - cn * (a + n) * U1
end

# Three-term recurrences of h_n in n (DLMF 13.3.7 arranged as in the BHPT paclet), with the
# c-derivatives of the c-dependent coefficients:
#   increasing n:  h_n = A_n h_(n-2) + B_n(c) h_(n-1)
#   decreasing n:  h_n = C_n(c) h_(n+1) + D_n h_(n+2)
function _forward_coefficients(a, b, c, n)
    den = (a + n - 1) * (b + 2n - 4)
    A = (b + n - a - 2) * (b + 2n - 2) / den
    bn = b + 2n
    B = (bn - 3) * ((8 + bn^2 - 6bn) / c + (2(a + n) - bn)) / den
    dB = -(bn - 3) * (8 + bn^2 - 6bn) / (c^2 * den)
    return A, B, dB
end
function _backward_coefficients(a, b, c, n)
    den = (b + n - a) * (b + 2n + 2)
    C = -(b + 2n + 1) * ((b^2 + 4n * (n + 1) + b * (2 + 4n)) / c + (2a - b)) / den
    dC = (b + 2n + 1) * (b^2 + 4n * (n + 1) + b * (2 + 4n)) / (c^2 * den)
    D = (a + n + 1) * (b + 2n) / den
    return C, dC, D
end

# The recurrences lose about log2(condition) bits per step, where the condition is the ratio of the
# sum of the term magnitudes to the result; they run at a precision raised by RECURRENCE_EXTRA_BITS
# and are re-anchored by two direct evaluations once the accumulated loss exceeds the re-anchor budget.
const RECURRENCE_EXTRA_BITS = 160
const RECURRENCE_REANCHOR_BITS = 120

@inline function _step_condition(h, t1, t2, dh, s1, s2, s3)
    R = real(typeof(h))
    cond_value = (abs(t1) + abs(t2)) / max(abs(h), floatmin(R))
    cond_deriv = (abs(s1) + abs(s2) + abs(s3)) / max(abs(dh), floatmin(R))
    return Float64(max(cond_value, cond_deriv, one(R)))
end

# Sum over n of weight_n h_n and weight_n dh_n/dzhat, adding terms on each side of n = 0 until they fall
# below the tolerance relative to the sum (the stored coefficients are extended as needed). The h_n
# come from the recurrences above with direct anchors. The sum runs in BigFloat at the caller's
# precision raised by RECURRENCE_EXTRA_BITS whatever the working type (a fixed-precision type such as
# Double64 cannot carry the extra bits of the recurrence), and the totals are returned in that type.
function _tricomi_sum(data::MSTGenericData, weight, a, b, c, dc_dzhat)
    C = typeof(c)
    R = real(C)
    CB = Complex{BigFloat}
    outer_bits = R === BigFloat ? precision(BigFloat) : precision(R)
    work_bits = outer_bits + RECURRENCE_EXTRA_BITS
    total, dtotal, absolute, nused = setprecision(BigFloat, work_bits) do
        aa = CB(a); bb = CB(b); cc = CB(c)
        total = zero(CB)
        dtotal = zero(CB)
        absolute = zero(BigFloat)
        nused = (0, 0)
        # anchors of both directions: h_(n - direction) and h_n are h_-1 and h_0 going up, h_0 and h_-1 going down
        h_minus, dh_minus = _h_direct(aa, bb, cc, -1)
        h_zero, dh_zero = _h_direct(aa, bb, cc, 0)
        for direction in (1, -1)
            n = direction == 1 ? 0 : -1
            h_prev1, dh_prev1 = direction == 1 ? (h_minus, dh_minus) : (h_zero, dh_zero)   # h_(n-dir)
            h_curr, dh_curr = direction == 1 ? (h_zero, dh_zero) : (h_minus, dh_minus)     # h_n
            loss = 0.0
            streak = 0
            first_step = true
            while true
                if direction == 1 ? n > data.nmax : n < data.nmin
                    _extend_coefficients!(data, max(-data.nmin, data.nmax) + COEFFICIENT_EXTENSION)
                end
                if !first_step
                    # advance the recurrence from (h_prev1, h_curr) = (h_(n-2dir), h_(n-dir)) to h_n
                    if direction == 1
                        A, B, dB = _forward_coefficients(aa, bb, cc, n)
                        t1 = A * h_prev1; t2 = B * h_curr
                        h_next = t1 + t2
                        s1 = A * dh_prev1; s2 = dB * h_curr; s3 = B * dh_curr
                        dh_next = s1 + s2 + s3
                    else
                        Cn, dCn, Dn = _backward_coefficients(aa, bb, cc, n)
                        t1 = Cn * h_curr; t2 = Dn * h_prev1
                        h_next = t1 + t2
                        s1 = dCn * h_curr; s2 = Cn * dh_curr; s3 = Dn * dh_prev1
                        dh_next = s1 + s2 + s3
                    end
                    loss += log2(_step_condition(h_next, t1, t2, dh_next, s1, s2, s3))
                    if loss > RECURRENCE_REANCHOR_BITS || !isfinite(loss)
                        # re-anchor with two direct values
                        h_prev1, dh_prev1 = _h_direct(aa, bb, cc, n - direction)
                        h_next, dh_next = _h_direct(aa, bb, cc, n)
                        loss = 0.0
                    else
                        h_prev1, dh_prev1 = h_curr, dh_curr
                    end
                    h_curr, dh_curr = h_next, dh_next
                end
                first_step = false
                w = CB(weight(n))
                term = w * h_curr
                total += term
                dtotal += w * dh_curr * CB(dc_dzhat)
                absolute += abs(term)
                if abs(n) >= SERIES_MIN_TERMS && abs(term) <= data.tolerance * abs(total)
                    streak += 1
                    streak >= 2 && break
                else
                    streak = 0
                end
                n += direction
            end
            nused = direction == 1 ? (nused[1], n) : (n, nused[2])
        end
        total, dtotal, absolute, nused
    end
    if R === BigFloat
        round_(x) = Complex{BigFloat}(BigFloat(real(x), precision=outer_bits), BigFloat(imag(x), precision=outer_bits))
        return round_(total), round_(dtotal), BigFloat(absolute, precision=outer_bits), nused
    end
    return C(total), C(dtotal), R(absolute), nused
end

"""
    mst_up_state(data, r)

(R_up, dR_up/dr) at the Boyer-Lindquist radius r, unit transmission at infinity. Returns a named
tuple with the state, the relative cancellation of the sum and the terms used.
"""
function mst_up_state(data::MSTGenericData, r)
    p = data.params
    R = _real_type(p)
    zhat = _zhat(data, R(r))
    c = -2im * zhat
    a = data.nu + p.s + 1 - im * data.epsilon
    b = 2 * data.nu + 2
    prefactor, dlog = _up_prefactor(data, zhat)
    total, dtotal, absolute, nused = _tricomi_sum(data, n -> _fup(data, n), a, b, c, -2im)
    norm = mst_up_transmission(data)
    value = prefactor * total / norm
    derivative = prefactor * (dlog * total + dtotal) * p.omega / norm   # dzhat/dr = omega
    return (R=value, dR=derivative, cancellation=Float64(absolute / max(abs(total), floatmin(R))), terms=nused)
end

# ---------------------------------------------------------------------------------------------
# DOWN solution: the ingoing branch at infinity (see the module header)
#   R_down^raw = e^{-i zhat} zhat^{nu + i eps_+} (zhat - eps kappa)^{-s - i eps_+}
#                sum_n f_n (2 i zhat)^n U(n+nu+1-s+i eps, 2n+2nu+2, 2 i zhat).
# ---------------------------------------------------------------------------------------------
function _down_prefactor(data::MSTGenericData, zhat)
    p = data.params
    e = data.epsilon
    k = p.kappa
    nu = data.nu
    eplus = (e + data.tau) / 2
    value = exp(-im * zhat) * zhat^(nu + im * eplus) * (zhat - e * k)^(-p.s - im * eplus)
    dlog = -im + (nu + im * eplus) / zhat + (-p.s - im * eplus) / (zhat - e * k)
    return value, dlog
end

"""
    mst_down_transmission(data)

Normalization of the raw DOWN series to a unit ingoing wave r^-1 exp(-i omega r*): its leading
large-r coefficient, (sum_n f_n) (2 i)^(-nu-1+s-i eps) omega^(-1-i eps) e^{i omega r_-} e^{-i eps log 2}
(with U(a, b, z) -> z^-a, zhat -> omega r and r* -> r + 2 log r - 2 log 2).
"""
function mst_down_transmission(data::MSTGenericData)
    data.down_normalization === nothing || return data.down_normalization
    data.down_normalization = _mst_down_transmission(data)
    return data.down_normalization
end
function _mst_down_transmission(data::MSTGenericData)
    p = data.params
    e = data.epsilon
    k = p.kappa
    nu = data.nu
    C = typeof(nu)
    R = real(C)
    total = zero(C)
    for n in data.nmin:data.nmax
        total += coefficient(data, n)
    end
    pi_ = _pi(nu)
    omega = e / 2
    power = -nu - 1 + p.s - im * e
    two_i_power = exp(power * (log(R(2)) + im * pi_ / 2))      # (2 i)^power, principal branch
    return total * two_i_power * omega^(-1 - im * e) * exp(im * omega * (1 - k)) * exp(-im * e * log(R(2)))
end

"""
    mst_down_state(data, r)

(R_down, dR_down/dr) at the Boyer-Lindquist radius r, unit ingoing amplitude at infinity.
"""
function mst_down_state(data::MSTGenericData, r)
    p = data.params
    R = _real_type(p)
    zhat = _zhat(data, R(r))
    c = 2im * zhat
    a = data.nu + 1 - p.s + im * data.epsilon
    b = 2 * data.nu + 2
    prefactor, dlog = _down_prefactor(data, zhat)
    total, dtotal, absolute, nused = _tricomi_sum(data, n -> coefficient(data, n), a, b, c, 2im)
    norm = mst_down_transmission(data)
    value = prefactor * total / norm
    derivative = prefactor * (dlog * total + dtotal) * p.omega / norm
    return (R=value, dR=derivative, cancellation=Float64(absolute / max(abs(total), floatmin(R))), terms=nused)
end

end
