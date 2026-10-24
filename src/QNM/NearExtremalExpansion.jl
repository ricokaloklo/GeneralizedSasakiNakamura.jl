# Near-extremal expansion of damped Kerr QNM frequencies, omega(x) = sum_k c_k x^k with
# x = 1 - a^2, solved order by order from the Teukolsky equation at a = 1 in ComplexF64.
#
# The radial function is written as R = P(r) sum_j x^j (r/(r-1))^(3j) sum_n a_n^(j) z^n,
# z = (r-2)/r, with Richartz's exact-extremal prefactor P; the angular function is expanded in
# spin-weighted spherical harmonics, where the spectral matrix is M0 + c M1 + c^2 M2 in
# c = a omega, a = sqrt(1 - x). Order 0 is the exact-extremal problem (Newton on the truncated,
# bordered radial-angular system, seeded by the a = 1 root); every order k >= 1 is one bordered
# linear solve for (a^(k), v^(k), c_k, A_k). Nothing is fitted. The operator blocks are generated
# in NearExtremalOperators.jl.

struct _NEBlock
    j::Int
    kk::Int
    cden::Int
    pole::Int
    g0::Function   # z-coefficients multiplying g
    g1::Function   # ... g'
    g2::Function   # ... g''
end

# Truncated power series in x (coefficients 0..K).
struct _NEJet{T<:Number} <: Number
    c::Vector{T}
end
_ne_const(x, K) = (v = zeros(typeof(x), K + 1); v[1] = x; _NEJet(v))
Base.:+(a::_NEJet, b::_NEJet) = _NEJet(a.c .+ b.c)
Base.:-(a::_NEJet, b::_NEJet) = _NEJet(a.c .- b.c)
Base.:-(a::_NEJet) = _NEJet(-a.c)
function Base.:*(a::_NEJet, b::_NEJet)
    K = length(a.c) - 1
    T = promote_type(eltype(a.c), eltype(b.c))
    out = zeros(T, K + 1)
    @inbounds for i in 0:K, j in 0:(K - i)
        out[i + j + 1] += a.c[i + 1] * b.c[j + 1]
    end
    return _NEJet(out)
end
Base.:+(a::_NEJet, b::Number) = (c = promote_type(eltype(a.c), typeof(b)).(a.c); c[1] += b; _NEJet(c))
Base.:+(b::Number, a::_NEJet) = a + b
Base.:-(a::_NEJet, b::Number) = a + (-b)
Base.:-(b::Number, a::_NEJet) = (-a) + b
Base.:*(a::_NEJet, b::Number) = _NEJet(a.c .* b)
Base.:*(b::Number, a::_NEJet) = _NEJet(a.c .* b)
Base.:/(a::_NEJet, b::Number) = _NEJet(a.c ./ b)
Base.zero(a::_NEJet) = _NEJet(zero.(a.c))
Base.iszero(a::_NEJet) = all(iszero, a.c)
Base.:^(a::_NEJet, p::Integer) =
    p == 0 ? _ne_const(one(eltype(a.c)), length(a.c) - 1) : p == 1 ? a : a * a^(p - 1)

include("NearExtremalOperators.jl")

# Taylor coefficients of sqrt(1 - x).
const _NE_SQRT1MX = (1.0, -1 / 2, -1 / 8, -1 / 16, -5 / 128, -7 / 256, -21 / 1024,
    -33 / 2048, -429 / 32768)
const _NE_ORDER = 8
const _NE_POLE = 3

_ne_falling(n, d) = d == 0 ? one(n) : prod(n - t for t in 0:(d - 1))

# polynomial (z-coefficients) times (1 + z)^e
function _ne_times_onepz(p::Vector, e::Int)
    isempty(p) && return p
    out = p
    for _ in 1:e
        z = zero(out[1])
        out = [(i <= length(out) ? out[i] : z) + (i > 1 ? out[i - 1] : z)
            for i in 1:(length(out) + 1)]
    end
    return out
end

# z-polynomials of a block at (w, A) (numbers, or jets of length K + 1), times factor (1+z)^e
function _ne_block_polys(b::_NEBlock, w, A, m, s, factor, e, K)
    usejet = w isa _NEJet || A isa _NEJet
    conv(x) = usejet ? (x isa _NEJet ? _NEJet(ComplexF64.(x.c)) :
        _ne_const(ComplexF64(x), K)) : ComplexF64(x)
    T = usejet ? _NEJet{ComplexF64} : ComplexF64
    return [_ne_times_onepz(T[conv(factor * x) for x in f(w, A, m, s)], e)
        for f in (b.g0, b.g1, b.g2)]
end

_ne_jetcoef(Q, t) = [[x isa _NEJet ? x.c[t + 1] : (t == 0 ? x : zero(x)) for x in q] for q in Q]

# rows n = 0..N of sum_d Q_d(z) g^(d)(z) applied to the coefficient vector a
function _ne_apply_rows(Q, a::AbstractVector, N::Int)
    out = zeros(ComplexF64, N + 1)
    for d in 0:2, (ii, q) in enumerate(Q[d + 1])
        iszero(q) && continue
        i = ii - 1
        for n in 0:N
            idx = n - i + d
            (0 <= idx <= N) || continue
            out[n + 1] += q * _ne_falling(idx, d) * a[idx + 1]
        end
    end
    return out
end

# Banded matrix of the same rows acting on a_0..a_N: kl subdiagonals and ku superdiagonals,
# stored LAPACK-style with kl extra rows for the fill-in of the LU factorisation.
struct _NEBand
    AB::Matrix{ComplexF64}
    kl::Int
    ku::Int
    n::Int
end
_ne_bget(B::_NEBand, i, j) = (i - j > B.kl || j - i > B.ku) ? zero(ComplexF64) :
    B.AB[B.kl + B.ku + 1 + i - j, j]

function _ne_rows_band(Q, N::Int)
    kl = max(0, maximum(length(Q[d + 1]) - 1 - d for d in 0:2)); ku = 2
    AB = zeros(ComplexF64, 2kl + ku + 1, N + 1)
    for d in 0:2, (ii, q) in enumerate(Q[d + 1])
        iszero(q) && continue
        i = ii - 1
        for n in 0:N
            idx = n - i + d
            (0 <= idx <= N) || continue
            AB[kl + ku + 1 + n - idx, idx + 1] += q * _ne_falling(idx, d)
        end
    end
    return _NEBand(AB, kl, ku, N + 1)
end

function _ne_band_mul(B::_NEBand, x::AbstractVector)
    y = zeros(ComplexF64, B.n)
    for j in 1:B.n, i in max(1, j - B.ku):min(B.n, j + B.kl)
        y[i] += B.AB[B.kl + B.ku + 1 + i - j, j] * x[j]
    end
    return y
end

# banded LU with partial pivoting (fill-in up to kl + ku superdiagonals), several right-hand sides
function _ne_banded_solve(B::_NEBand, rhs::Matrix{ComplexF64})
    n, kl, ku = B.n, B.kl, B.ku
    AB = copy(B.AB)
    get(i, j) = AB[kl + ku + 1 + i - j, j]
    set!(i, j, v) = (AB[kl + ku + 1 + i - j, j] = v)
    X = copy(rhs)
    for k in 1:n
        lastr = min(n, k + kl)
        p = k; best = abs(get(k, k))
        for i in (k + 1):lastr
            a = abs(get(i, k)); a > best && (best = a; p = i)
        end
        lastc = min(n, k + kl + ku)
        if p != k
            for j in k:lastc
                t1 = get(k, j); set!(k, j, get(p, j)); set!(p, j, t1)
            end
            X[k, :], X[p, :] = X[p, :], X[k, :]
        end
        piv = get(k, k)
        for i in (k + 1):lastr
            f = get(i, k) / piv
            iszero(f) && continue
            set!(i, k, zero(ComplexF64))
            for j in (k + 1):lastc
                set!(i, j, get(i, j) - f * get(k, j))
            end
            @views X[i, :] .-= f .* X[k, :]
        end
    end
    for k in n:-1:1
        lastc = min(n, k + kl + ku)
        for j in (k + 1):lastc
            @views X[k, :] .-= get(k, j) .* X[j, :]
        end
        @views X[k, :] ./= get(k, k)
    end
    return X
end

# Solve Md g + c colw + d colA = rr, g[1] = r0, Ma v + c angw + d vcol = ar, v[lidx] = rv. The
# truncated radial block is nearly singular at a root; on the constraint g[1] = r0 it acts like
# Md + alpha e0 e0^T, which is well conditioned and banded.
function _ne_bordered_solve(Md, colw, colA, rr, r0, Ma, angw, vcol, ar, rv, lidx)
    T = eltype(rr)
    alpha = one(T) + maximum(abs(_ne_bget(Md, 1, j)) for j in 1:min(Md.n, 1 + Md.ku))
    AB = copy(Md.AB); AB[Md.kl + Md.ku + 1, 1] += alpha
    Mdef = _NEBand(AB, Md.kl, Md.ku, Md.n)
    rhs1 = copy(rr); rhs1[1] += alpha * r0
    Y = _ne_banded_solve(Mdef, hcat(rhs1, colw, colA))
    Na = length(ar)
    Mbig = zeros(T, Na + 2, Na + 2)
    Mbig[1:Na, 1:Na] .= Ma; Mbig[1:Na, Na + 1] .= angw; Mbig[1:Na, Na + 2] .= vcol
    Mbig[Na + 1, lidx] = one(T)
    Mbig[Na + 2, Na + 1] = Y[1, 2]; Mbig[Na + 2, Na + 2] = Y[1, 3]
    sol = Mbig \ vcat(ar, rv, Y[1, 1] - r0)
    v = sol[1:Na]; c = sol[Na + 1]; d = sol[Na + 2]
    g = Y[:, 1] .- c .* Y[:, 2] .- d .* Y[:, 3]
    return g, v, c, d
end

struct _NESetup
    s::Int; l::Int; m::Int
    N::Int; Na::Int
    M0::Matrix{ComplexF64}; M1::Matrix{ComplexF64}; M2::Matrix{ComplexF64}
    lidx::Int
end

function _NESetup(s, l, m; N, Na)
    sm(c) = ComplexF64.(_big_angular_spectral_matrix(ComplexF64(c), s, m, Na))
    M0 = sm(0); Mp = sm(1); Mm = sm(-1)
    return _NESetup(s, l, m, N, Na, M0, (Mp - Mm) / 2, (Mp + Mm) / 2 - M0,
        l - max(abs(m), abs(s)) + 1)
end

_ne_order0_polys(S, w, A, K) = _ne_block_polys(_NE_OPERATOR_BLOCKS[(0, 0)], w, A, S.m, S.s, 1.0, 0, K)

# Order 0 at the given exact-extremal root: c_0 = omega0 is not iterated. The radial null vector
# a of the truncated block solves (Md + alpha e0 e0^T) a = alpha e0 (a_0 = 1 when Md a = 0), and
# the angular eigenpair (A_0, v) is the eigenvalue of M(omega0) nearest to `A0` with v[lidx] = 1.
# Returns the relative radial residual |Md a| / (|Md| |a|) with it.
function _ne_order0_at_root(S::_NESetup, w0, A0)
    N, Na = S.N, S.Na
    w = ComplexF64(w0)
    Ma = S.M0 + w * S.M1 + w^2 * S.M2
    F = eigen(Ma)
    j = argmin(abs.(F.values .- A0))
    A = F.values[j]
    v = F.vectors[:, j] ./ F.vectors[S.lidx, j]
    Mr = _ne_rows_band(_ne_order0_polys(S, w, A, 1), N)
    alpha = one(ComplexF64) + maximum(abs(_ne_bget(Mr, 1, k)) for k in 1:min(Mr.n, 1 + Mr.ku))
    AB = copy(Mr.AB); AB[Mr.kl + Mr.ku + 1, 1] += alpha
    rhs = zeros(ComplexF64, N + 1, 1); rhs[1, 1] = alpha
    a = vec(_ne_banded_solve(_NEBand(AB, Mr.kl, Mr.ku, Mr.n), rhs))
    scale = maximum(abs, Mr.AB)
    residual = maximum(abs, _ne_band_mul(Mr, a)) / (scale * maximum(abs, a))
    return (; w, A, a, v, residual)
end

function _ne_orders(S::_NESetup, o0, K)
    N, Na, s, m = S.N, S.Na, S.s, S.m
    cs = ComplexF64[o0.w]; ds = ComplexF64[o0.A]; gs = [o0.a]; vs = [o0.v]
    ajet = _NEJet(ComplexF64.(collect(_NE_SQRT1MX[1:(K + 1)])))
    for k in 1:K
        wj = _NEJet(vcat(cs, zeros(ComplexF64, K + 1 - k)))
        Aj = _NEJet(vcat(ds, zeros(ComplexF64, K + 1 - k)))
        rk = zeros(ComplexF64, N + 1)
        for j in 0:(k - 1), kk in 0:(k - j)
            b = _NE_OPERATOR_BLOCKS[(j, kk)]
            factor = 4.0 / b.cden / 2.0^(_NE_POLE * (k - j))
            e = 1 + _NE_POLE * (k - j) - b.pole
            Q = _ne_jetcoef(_ne_block_polys(b, wj, Aj, m, s, factor, e, K), k - j - kk)
            rk .+= _ne_apply_rows(Q, gs[j + 1], N)
        end
        b00 = _NE_OPERATOR_BLOCKS[(0, 0)]
        f00 = 4.0 / b00.cden / 2.0^(_NE_POLE * k); e00 = 1 + _NE_POLE * k - b00.pole
        e1 = vcat([o0.w, 1.0 + 0im], zeros(ComplexF64, K - 1))
        e2 = vcat([o0.A, 1.0 + 0im], zeros(ComplexF64, K - 1))
        Qw = _ne_jetcoef(_ne_block_polys(b00, _NEJet(e1), _ne_const(o0.A, K), m, s, f00, e00, K), 1)
        QA = _ne_jetcoef(_ne_block_polys(b00, _ne_const(o0.w, K), _NEJet(e2), m, s, f00, e00, K), 1)
        bkk = _NE_OPERATOR_BLOCKS[(k, 0)]
        Md = _ne_rows_band(_ne_block_polys(bkk, o0.w, o0.A, m, s, 4.0 / bkk.cden,
            1 - bkk.pole, K), N)
        cjet = ajet * wj
        c2 = cjet * cjet
        Mjet(t) = (t == 0 ? S.M0 : zero(S.M0)) + cjet.c[t + 1] * S.M1 + c2.c[t + 1] * S.M2
        ak = zeros(ComplexF64, Na)
        for i in 1:k; ak .+= Mjet(i) * vs[k - i + 1]; end
        for i in 1:(k - 1); ak .-= ds[i + 1] * vs[k - i + 1]; end
        Ma0 = S.M0 + o0.w * S.M1 + o0.w^2 * S.M2 - o0.A * I
        g, v, c, d = _ne_bordered_solve(Md, _ne_apply_rows(Qw, gs[1], N),
            _ne_apply_rows(QA, gs[1], N), -rk, zero(ComplexF64), Ma0,
            (S.M1 + 2 * o0.w * S.M2) * vs[1], -vs[1], -ak, zero(ComplexF64), S.lidx)
        push!(gs, g); push!(vs, v); push!(cs, c); push!(ds, d)
    end
    return cs
end

"""
    near_extremal_expansion_coefficients(mode, omega0; maximum_N=20000)

Coefficients `c_0, ..., c_8` of `omega(x) = sum_k c_k x^k`, `x = 1 - a^2`, for the damped
exact-extremal root `omega0` of `mode`, in ComplexF64, with `c_0 = omega0` taken as given (it is
not re-solved, so the expansion stays on the root it is handed). The radial truncation starts at
300 and grows by 1.5 until the order-0 coefficients at the cut are below 1e-17 of their maximum.
Returns `nothing` if the truncated radial system does not annihilate its null vector to 1e-10
(`omega0` is not a root at that truncation) or the truncation exceeds `maximum_N`.
"""
function near_extremal_expansion_coefficients(mode::QNMMode, omega0; maximum_N::Int=20000)
    s, l, m = mode.s, mode.l, mode.m
    w0 = ComplexF64(omega0)
    Na = l - max(abs(m), abs(s)) + 25
    A0 = angular_branch(mode, 1.0, w0; truncation_order=Na,
        sheet_id=:near_extremal_expansion).angular_A
    N = 300
    while true
        S = _NESetup(s, l, m; N, Na)
        o = _ne_order0_at_root(S, w0, A0)
        tail = maximum(abs, @view o.a[(end - 20):end]) / maximum(abs, o.a)
        if tail < 1.0e-17
            o.residual <= 1.0e-10 || return nothing
            return _ne_orders(S, o, _NE_ORDER)
        end
        N >= maximum_N && return nothing
        N = min(maximum_N, round(Int, 1.5N))
    end
end

# Per-mode coefficient cache (coefficients or `nothing`), shared across threads.
const _NE_CACHE = Dict{NTuple{4,Int},Any}()
const _NE_CACHE_LOCK = ReentrantLock()

function _cached_near_extremal_coefficients(mode::QNMMode, omega0)
    key = (mode.s, mode.l, mode.m, mode.n)
    lock(_NE_CACHE_LOCK) do
        haskey(_NE_CACHE, key) && return _NE_CACHE[key]
    end
    c = near_extremal_expansion_coefficients(mode, omega0)
    lock(_NE_CACHE_LOCK) do
        _NE_CACHE[key] = c === nothing ? nothing : Tuple(c)
    end
    return c === nothing ? nothing : Tuple(c)
end
