module DirectParameters

using LinearAlgebra
using LinearAlgebra.BLAS: @blasfunc, BlasInt
using SpinWeightedSpheroidalHarmonics

export DirectGSNParameters, DirectGSNControls, DirectFrequencySelection
export DirectSpinSelection
export direct_gsn_parameters, direct_gsn_controls
export direct_sfe_regime, direct_lfe_regime, direct_frequency_regime
export direct_spin_regime
export direct_horizon_tail
export direct_swsh_eigenvalue
export supported_spin_weight

const DEFAULT_HORIZON_ORDER = 24
const DEFAULT_ORDINARY_ORDER = 24
const DEFAULT_REGULAR_ORDINARY_ORDER = 40
const DEFAULT_INFINITY_ORDER = 24
const DEFAULT_MATCH_X = 0.5
const DEFAULT_TOLERANCE = 1e-14
const SFE_FREQUENCY_LIMIT = 0.1
const SFE_TRANSITION_LIMIT = 0.15
const SFE_KAPPA_OMEGA_LIMIT = 1e-3
const SFE_HORIZON_WAVENUMBER_MIN = 4.0
const HORIZON_TAIL_DETUNING_LIMIT = 0.1
const SFE_ZERO_LAMBDA_ORDINARY_MIN = 1e-4
const SFE_ZERO_LAMBDA_ULPS = 128.0
const LFE_PHASE_THRESHOLD = 20.0
const LFE_NEAR_EXTREME_FREQUENCY = 10.0
const LFE_KAPPA_OMEGA_LIMIT = 0.1
const NEAR_EXTREME_KAPPA_LIMIT = 0.02
const NEAR_EXTREME_COMPRESSION_MIN = 512.0
const SMALL_C_SWSH_LIMIT = 1e-6
const SWSH_EIGENVALUE_ATOL = 5e-14
const SWSH_EIGENVALUE_RTOL = 5e-15
const SWSH_MIN_BUFFER = 10
const SWSH_MAX_REFINEMENTS = 20
const SWSH_LAPACK = LinearAlgebra.LAPACK.liblapack

supported_spin_weight(s::Integer) = s in (-2, -1, 0, 1, 2)

# Mode parameters in units M = 1, with a = sin(nu) and kappa = cos(nu) = sqrt(1 - a^2).
# omega and lambda are Float64 for real frequencies and ComplexF64 otherwise.
struct DirectGSNParameters{T<:Number}
    s::Int
    l::Int
    m::Int
    a::Float64
    omega::T
    lambda::T
    nu::Float64
    kappa::Float64
end

# Outcome of the small/large-frequency regime selection, with the features it used.
struct DirectFrequencySelection
    sfe_request::Symbol
    lfe_request::Symbol
    regime::Symbol
    reason::Symbol
    omega_abs::Float64
    kappa_omega_abs::Float64
    tau_abs::Float64
    horizon_wavenumber_abs::Float64
    m_ratio::Float64
    a_abs::Float64
    mst_triggered::Bool
end

# Outcome of the regular / near-extremal / exact-extremal spin classification.
struct DirectSpinSelection
    regime::Symbol
    reason::Symbol
    kappa::Float64
    angular_compression::Float64
    phase_compression::Float64
    tau_abs::Float64
    horizon_detuning_abs::Float64
    m_ratio::Float64
end

struct DirectGSNControls
    horizon_order::Int
    ordinary_order::Int
    infinity_order::Int
    match_x::Float64
    tolerance::Float64
    basis::Symbol
    endpoint_basis::NamedTuple
    ordinary_basis::Symbol
    sfe::Bool
    lfe::Bool
    frequency_selection::DirectFrequencySelection
    source::Symbol
end

function Base.show(io::IO, ::MIME"text/plain", p::DirectGSNParameters)
    println(io, "DirectGSNParameters:")
    println(io, "  s      = ", p.s)
    println(io, "  l      = ", p.l)
    println(io, "  m      = ", p.m)
    println(io, "  a      = ", p.a)
    println(io, "  omega  = ", p.omega)
    println(io, "  lambda = ", p.lambda)
    println(io, "  nu     = ", p.nu)
    print(io, "  kappa  = ", p.kappa)
end

function Base.show(io::IO, ::MIME"text/plain", c::DirectGSNControls)
    println(io, "DirectGSNControls:")
    println(io, "  horizon_order  = ", c.horizon_order)
    println(io, "  ordinary_order = ", c.ordinary_order)
    println(io, "  infinity_order = ", c.infinity_order)
    println(io, "  match_x        = ", c.match_x)
    println(io, "  tolerance      = ", c.tolerance)
    println(io, "  basis          = ", c.basis)
    println(io, "  endpoint_basis = ", c.endpoint_basis)
    println(io, "  ordinary_basis = ", c.ordinary_basis)
    println(io, "  sfe            = ", c.sfe)
    println(io, "  lfe            = ", c.lfe)
    println(io, "  regime         = ", c.frequency_selection.regime)
    println(io, "  regime_reason  = ", c.frequency_selection.reason)
    println(io, "  mst_triggered  = ", c.frequency_selection.mst_triggered)
    print(io, "  source         = ", c.source)
end

@inline _binary64_compatible(value::Real) =
    value isa Integer || value isa Rational ||
    value isa Float16 || value isa Float32 || value isa Float64

@inline _binary64_compatible(value::Complex) =
    _binary64_compatible(real(value)) && _binary64_compatible(imag(value))

@inline _binary64_compatible(value) = false

function _require_binary64_input(name::AbstractString, value)
    _binary64_compatible(value) && return value
    throw(ArgumentError(
        "$name has type $(typeof(value)); Direct GSN currently supports " *
        "binary64-compatible inputs only and will not silently downcast " *
        "arbitrary-precision floating values."))
end

function _real_float(name::AbstractString, value)
    _require_binary64_input(name, value)
    z = complex(value)
    iszero(imag(z)) || throw(ArgumentError("$name must be real for the first direct GSN implementation."))
    return Float64(real(z))
end

function _checked_spin(s, l, m)
    supported_spin_weight(s) || throw(ArgumentError("direct GSN coefficient tables support s = -2, -1, 0, 1, 2."))
    l >= abs(s) || throw(ArgumentError("l must satisfy l >= abs(s)."))
    abs(m) <= l || throw(ArgumentError("m must satisfy abs(m) <= l."))
    return Int(s), Int(l), Int(m)
end

# Coupling coefficients of the spherical-harmonic expansion of the spin-weighted
# spheroidal eigenproblem; the matrix below is pentadiagonal in l.
@inline function _swsh_f(s::Int, l::Int, m::Int)
    (l == -1 && iszero(m) && iszero(s)) && return 0.0
    lp1 = l + 1
    return sqrt(((lp1^2 - m^2) / ((2l + 3) * (2l + 1))) *
        ((lp1^2 - s^2) / lp1^2))
end

@inline function _swsh_g(s::Int, l::Int, m::Int)
    iszero(l) && return 0.0
    return sqrt(((l^2 - m^2) / (4l^2 - 1)) * ((l^2 - s^2) / l^2))
end

@inline _swsh_h(s::Int, l::Int, m::Int) =
    (iszero(l) || iszero(s)) ? 0.0 : -m * s / (l * (l + 1))

# Second-order small-c expansion of lambda for |c| <= SMALL_C_SWSH_LIMIT.
@inline function _small_c_swsh_eigenvalue(s::Int, l::Int, m::Int, c::Float64)
    g_lower = _swsh_g(s, l, m)
    g_upper = _swsh_g(s, l + 1, m)
    h = _swsh_h(s, l, m)
    b = g_lower^2 + g_upper^2 + h^2
    lower_shift = iszero(l) ? 0.0 : g_lower^2 / l
    lambda0 = Float64(l * (l + 1) - s * (s + 1))
    lambda1 = 2s * h - 2m
    lambda2 = 1 - b + 2s^2 * (lower_shift - g_upper^2 / (l + 1))
    return muladd(c, muladd(c, lambda2, lambda1), lambda0)
end

@inline function _swsh_matrix_coefficient(c, s::Int, m::Int, l::Int, lprime::Int)
    if lprime == l
        f = _swsh_f(s, l, m)
        g = _swsh_g(s, l, m)
        h = _swsh_h(s, l, m)
        b = f * _swsh_g(s, l + 1, m) + g * _swsh_f(s, l - 1, m) + h^2
        return l * (l + 1) - s * (s + 1) - c^2 * b + 2c * s * h
    elseif lprime == l + 1
        e = _swsh_g(s, lprime, m) *
            (_swsh_h(s, lprime - 1, m) + _swsh_h(s, lprime, m))
        return -c^2 * e + 2c * s * _swsh_g(s, lprime, m)
    elseif lprime == l + 2
        return -c^2 * _swsh_g(s, lprime, m) * _swsh_g(s, lprime - 1, m)
    end
    return 0.0
end

# Matrix element whose eigenvalue is lambda directly (diagonal shifted by c^2 - 2 m c).
@inline function _swsh_lambda_matrix_coefficient(
    c::Float64,
    s::Int,
    m::Int,
    l::Int,
    lprime::Int,
)
    lprime == l || return _swsh_matrix_coefficient(c, s, m, l, lprime)
    f = _swsh_f(s, l, m)
    g = _swsh_g(s, l, m)
    h = _swsh_h(s, l, m)
    b = f * _swsh_g(s, l + 1, m) +
        g * _swsh_f(s, l - 1, m) + h^2
    spherical = Float64(l * (l + 1) - s * (s + 1))
    return muladd(c^2, 1 - b, muladd(2c, s * h - m, spherical))
end

# The `index`-th smallest eigenvalue of a symmetric band matrix in LAPACK upper
# band storage, via dsbevx (eigenvalues only).
function _selected_banded_eigenvalue!(band::Matrix{Float64}, index::Int)
    n = BlasInt(size(band, 2))
    kd = BlasInt(size(band, 1) - 1)
    leading_band = BlasInt(size(band, 1))
    q = zeros(Float64, 1)
    leading_q = BlasInt(1)
    lower_value = 0.0
    upper_value = 0.0
    lower_index = BlasInt(index)
    upper_index = BlasInt(index)
    # LAPACK recommends twice the safe minimum for high relative accuracy.
    absolute_tolerance = 2floatmin(Float64)
    found = Ref{BlasInt}()
    eigenvalues = zeros(Float64, Int(n))
    eigenvectors = zeros(Float64, 1)
    leading_eigenvectors = BlasInt(1)
    work = zeros(Float64, 7 * Int(n))
    integer_work = zeros(BlasInt, 5 * Int(n))
    failed = zeros(BlasInt, Int(n))
    info = Ref{BlasInt}()

    ccall((@blasfunc(dsbevx_), SWSH_LAPACK), Cvoid,
        (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
            Ref{BlasInt}, Clong, Clong, Clong),
        'N', 'I', 'U', n, kd, band, leading_band, q, leading_q,
        lower_value, upper_value, lower_index, upper_index,
        absolute_tolerance, found, eigenvalues, eigenvectors,
        leading_eigenvectors, work, integer_work, failed, info, 1, 1, 1)

    info[] == 0 || error("LAPACK dsbevx failed with info=$(info[]).")
    found[] == 1 || error("LAPACK dsbevx returned $(found[]) eigenvalues.")
    return eigenvalues[1]
end

# Eigenvalue for (s, l, m, c) from a spectral matrix truncated to `size` l-modes.
function _swsh_eigenvalue_at_size(
    s::Int,
    l::Int,
    m::Int,
    c::Float64,
    size::Int,
)
    lmin = max(abs(m), abs(s))
    index = l - lmin + 1
    size >= index || throw(ArgumentError("SWSH matrix does not contain l=$l."))
    bandwidth = 2
    band = zeros(Float64, bandwidth + 1, size)
    @inbounds for row_l in lmin:(lmin + size - 1)
        i = row_l - lmin + 1
        for col_l in row_l:min(row_l + bandwidth, lmin + size - 1)
            j = col_l - lmin + 1
            band[bandwidth + 1 + i - j, j] =
                _swsh_lambda_matrix_coefficient(c, s, m, row_l, col_l)
        end
    end
    return _selected_banded_eigenvalue!(band, index)
end

# Spin-weighted spheroidal eigenvalue lambda for real c = a omega. Exact at c = 0,
# series for tiny |c|, otherwise the banded spectral matrix is enlarged until the
# eigenvalue converges. `diagnostics=true` returns the size and convergence data too.
function direct_swsh_eigenvalue(
    s::Int,
    l::Int,
    m::Int,
    c::Real;
    diagnostics::Bool=false,
)
    c_float = Float64(c)
    if iszero(c)
        value = Float64(l * (l + 1) - s * (s + 1))
        result = (value=value, size=0, refinement=0, delta=0.0,
            method=:spherical)
        return diagnostics ? result : result.value
    end
    if abs(c_float) <= SMALL_C_SWSH_LIMIT
        value = _small_c_swsh_eigenvalue(s, l, m, c_float)
        result = (value=value, size=0, refinement=0, delta=0.0,
            method=:small_c)
        return diagnostics ? result : result.value
    end

    lmin = max(abs(m), abs(s))
    index = l - lmin + 1
    size = max(index + SWSH_MIN_BUFFER,
        index + ceil(Int, abs(c_float)) + 8)
    step = max(8, ceil(Int, abs(c_float) / 8))
    previous = _swsh_eigenvalue_at_size(s, l, m, c_float, size)

    for refinement in 1:SWSH_MAX_REFINEMENTS
        next_size = size + step
        current = _swsh_eigenvalue_at_size(s, l, m, c_float, next_size)
        delta = abs(current - previous)
        threshold = SWSH_EIGENVALUE_ATOL +
            SWSH_EIGENVALUE_RTOL * max(abs(previous), abs(current))
        if delta <= threshold
            result = (value=current, size=next_size, refinement=refinement,
                delta=delta, method=:adaptive_spectral)
            return diagnostics ? result : result.value
        end
        size = next_size
        previous = current
    end

    throw(ErrorException("SWSH eigenvalue failed to converge after $(SWSH_MAX_REFINEMENTS) matrix refinements."))
end

_direct_swsh_eigenvalue(s::Int, l::Int, m::Int, c::Float64) =
    direct_swsh_eigenvalue(s, l, m, c)

# Validate inputs and build DirectGSNParameters. lambda defaults to the
# spheroidal eigenvalue (real omega) or the SpinWeightedSpheroidalHarmonics value
# (complex omega); nu, if given, must satisfy sin(nu) = a.
function direct_gsn_parameters(s::Integer, l::Integer, m::Integer, a, omega; lambda=nothing, nu=nothing)
    ss, ll, mm = _checked_spin(s, l, m)
    a_float = _real_float("a", a)
    _require_binary64_input("omega", omega)
    omega_complex = ComplexF64(omega)
    isfinite(real(omega_complex)) && isfinite(imag(omega_complex)) ||
        throw(ArgumentError("omega must be finite."))
    abs(a_float) < 1 || throw(ArgumentError("direct GSN requires abs(a) < 1."))
    kappa = sqrt((1 - a_float) * (1 + a_float))

    nu_float = if nu === nothing
        atan(a_float, kappa)
    else
        value = _real_float("nu", nu)
        isapprox(sin(value), a_float; atol=100eps(Float64), rtol=100eps(Float64)) ||
            throw(ArgumentError("provided nu is inconsistent with a = sin(nu)."))
        value
    end
    kappa > 0 || throw(ArgumentError("direct GSN requires a non-extremal positive kappa."))

    if iszero(imag(omega_complex))
        omega_float = Float64(real(omega_complex))
        lambda_float = if lambda === nothing
            value = _direct_swsh_eigenvalue(ss, ll, mm, a_float * omega_float)
            _real_float("lambda", value)
        else
            _real_float("lambda", lambda)
        end
        return DirectGSNParameters(
            ss, ll, mm, a_float, omega_float, lambda_float, nu_float, kappa)
    end

    lambda_complex = if lambda === nothing
        ComplexF64(
            SpinWeightedSpheroidalHarmonics.spin_weighted_spheroidal_eigenvalue(
                ss, ll, mm, a_float * omega_complex))
    else
        _require_binary64_input("lambda", lambda)
        ComplexF64(lambda)
    end
    isfinite(real(lambda_complex)) && isfinite(imag(lambda_complex)) ||
        throw(ArgumentError("lambda must be finite."))

    return DirectGSNParameters(
        ss, ll, mm, a_float, omega_complex, lambda_complex, nu_float, kappa)
end

# Real-frequency large-phase test used by lfe = :auto.
function direct_lfe_regime(params::DirectGSNParameters)
    isreal(params.omega) || return false
    horizon_frequency = params.a / (2 * (1 + params.kappa))
    horizon_wavenumber = params.omega - params.m * horizon_frequency
    phase_scale = max(abs(params.omega), abs(horizon_wavenumber))
    near_extreme = abs(params.omega) >= LFE_NEAR_EXTREME_FREQUENCY &&
        abs(params.kappa * params.omega) < LFE_KAPPA_OMEGA_LIMIT
    return phase_scale >= LFE_PHASE_THRESHOLD || near_extreme
end

# True when the horizon wavenumber omega - m Omega_H is small relative to max(|omega|, 1).
@inline function direct_horizon_tail(params::DirectGSNParameters)
    horizon_frequency = params.a / (2 * (1 + params.kappa))
    wavenumber = params.omega - params.m * horizon_frequency
    return abs(wavenumber) /
        max(abs(params.omega), one(abs(params.omega))) <
        HORIZON_TAIL_DETUNING_LIMIT
end

# Real-frequency small-frequency test used by sfe = :auto.
@inline function direct_sfe_regime(params::DirectGSNParameters)
    isreal(params.omega) || return false
    _sfe_zero_lambda_ordinary(params) && return false
    features = _frequency_features(params)
    0 < features.omega_abs <= SFE_FREQUENCY_LIMIT && return true
    return features.omega_abs <= SFE_TRANSITION_LIMIT &&
        features.kappa_omega_abs <= SFE_KAPPA_OMEGA_LIMIT &&
        features.horizon_wavenumber_abs >= SFE_HORIZON_WAVENUMBER_MIN
end

# Small-frequency regular-spin modes whose lambda vanishes to rounding level;
# these are sent to the ordinary route instead of the small-frequency one.
@inline function _sfe_zero_lambda_ordinary(params::DirectGSNParameters)
    omega = abs(params.omega)
    SFE_ZERO_LAMBDA_ORDINARY_MIN <= omega <= SFE_FREQUENCY_LIMIT ||
        return false
    direct_spin_regime(params).regime == :regular || return false
    c = params.a * params.omega
    angular_scale = max(
        1.0,
        abs(Float64(params.l * (params.l + 1))),
        abs(Float64(params.s * (params.s + 1))),
        2 * abs(params.m * c),
        abs2(c),
    )
    return abs(params.lambda) <=
        SFE_ZERO_LAMBDA_ULPS * eps(Float64) * angular_scale
end

# Normalise an sfe/lfe keyword to :on, :off or :auto.
@inline function _frequency_request(name::AbstractString, value)
    value === nothing && return :off
    value === true && return :on
    value === false && return :off
    value === :auto && return :auto
    value isa Integer && value == 1 && return :on
    value isa Integer && value == 0 && return :off
    if value isa AbstractString || value isa Symbol
        text = lowercase(String(value))
        text in ("true", "on") && return :on
        text in ("false", "off") && return :off
        text == "auto" && return :auto
    end
    throw(ArgumentError("$name must be true, false, :auto, or \"auto\"."))
end

@inline function _frequency_features(params::DirectGSNParameters)
    epsilon = 2 * params.omega
    horizon_frequency = params.a / (2 * (1 + params.kappa))
    horizon_wavenumber = params.omega - params.m * horizon_frequency
    return (
        omega_abs=abs(params.omega),
        kappa_omega_abs=abs(params.kappa * params.omega),
        tau_abs=abs((epsilon - params.m * params.a) / params.kappa),
        horizon_wavenumber_abs=abs(horizon_wavenumber),
        m_ratio=abs(params.m) / params.l,
        a_abs=abs(params.a),
    )
end

# Classify the spin as :exact_extremal (kappa = 0), :near_extremal (small kappa
# with large angular, phase or horizon-detuning compression), or :regular.
@inline function direct_spin_regime(params::DirectGSNParameters)
    features = _frequency_features(params)
    angular_scale = params.l + 0.5
    kappa_scale = max(params.kappa, eps(Float64))
    angular_compression = angular_scale / kappa_scale
    phase_scale = max(
        features.omega_abs,
        features.horizon_wavenumber_abs,
        abs(params.m * params.a) / (2 * (1 + params.kappa)),
    )
    phase_compression = phase_scale / kappa_scale

    if iszero(params.kappa)
        regime, reason = :exact_extremal, :zero_kappa
    elseif params.kappa <= NEAR_EXTREME_KAPPA_LIMIT &&
            max(angular_compression, phase_compression,
                features.tau_abs) >= NEAR_EXTREME_COMPRESSION_MIN
        reason = angular_compression >= max(phase_compression,
            features.tau_abs) ? :angular_compression :
            phase_compression >= features.tau_abs ? :phase_compression :
            :horizon_detuning
        regime = :near_extremal
    else
        regime, reason = :regular, :regular_spin
    end

    return DirectSpinSelection(
        regime,
        reason,
        params.kappa,
        angular_compression,
        phase_compression,
        features.tau_abs,
        features.horizon_wavenumber_abs,
        features.m_ratio,
    )
end

# Choose the frequency regime (:sfe, :lfe or :ordinary) from the sfe/lfe requests
# and the mode parameters, recording the reason.
@inline function direct_frequency_regime(
    params::DirectGSNParameters;
    sfe=false,
    lfe=false,
)
    sfe_request = _frequency_request("sfe", sfe)
    lfe_request = _frequency_request("lfe", lfe)
    sfe_request == :on && lfe_request == :on &&
        throw(ArgumentError("sfe=true and lfe=true are mutually exclusive."))
    complex_frequency = !isreal(params.omega)
    principal_mst_frequency =
        iszero(real(params.omega)) && imag(params.omega) < 0 &&
        abs(params.omega) <= SFE_FREQUENCY_LIMIT
    complex_frequency && sfe_request == :on && !principal_mst_frequency &&
        throw(ArgumentError(
            "sfe=true for complex omega is restricted to the negative " *
            "pure-imaginary principal-MST axis."))
    complex_frequency && lfe_request == :on && throw(ArgumentError(
        "lfe=true is available only for real omega; use the ordinary complex-frequency route."))

    regime, reason = if principal_mst_frequency && sfe_request == :on
        (:sfe, :forced_principal_mst)
    elseif complex_frequency
        (:ordinary, :complex_frequency)
    elseif sfe_request == :on
        (:sfe, :forced_sfe)
    elseif lfe_request == :on
        (:lfe, :forced_lfe)
    elseif sfe_request == :auto && _sfe_zero_lambda_ordinary(params)
        (:ordinary, :zero_lambda_ordinary)
    elseif sfe_request == :auto && direct_sfe_regime(params)
        (:sfe, :validated_small_frequency)
    elseif lfe_request == :auto && direct_lfe_regime(params)
        (:lfe, :validated_large_frequency)
    elseif iszero(params.omega) && sfe_request == :auto
        (:ordinary, :static_frequency)
    else
        (:ordinary, :ordinary_frequency)
    end

    features = _frequency_features(params)
    return DirectFrequencySelection(
        sfe_request,
        lfe_request,
        regime,
        reason,
        features.omega_abs,
        features.kappa_omega_abs,
        features.tau_abs,
        features.horizon_wavenumber_abs,
        features.m_ratio,
        features.a_abs,
        regime == :sfe,
    )
end

# Copy of `ctrl` switched to the ordinary regime after the small-frequency (MST)
# route was rejected; mst_triggered is set to true to record that attempt.
function _ordinary_fallback(ctrl::DirectGSNControls, reason::Symbol)
    selection = ctrl.frequency_selection
    fallback = DirectFrequencySelection(
        selection.sfe_request,
        selection.lfe_request,
        :ordinary,
        reason,
        selection.omega_abs,
        selection.kappa_omega_abs,
        selection.tau_abs,
        selection.horizon_wavenumber_abs,
        selection.m_ratio,
        selection.a_abs,
        true,
    )
    return DirectGSNControls(
        ctrl.horizon_order,
        ctrl.ordinary_order,
        ctrl.infinity_order,
        ctrl.match_x,
        ctrl.tolerance,
        ctrl.basis,
        ctrl.endpoint_basis,
        ctrl.ordinary_basis,
        false,
        false,
        fallback,
        ctrl.source,
    )
end

# Solver controls: series orders, matching point xm in (0, 1), tolerance, bases
# and the selected frequency regime. The ordinary order is raised to
# DEFAULT_REGULAR_ORDINARY_ORDER in the ordinary regime unless N or ordinary_order is given.
function direct_gsn_controls(
    params::DirectGSNParameters;
    N=nothing,
    horizon_order=nothing,
    ordinary_order=nothing,
    infinity_order=nothing,
    xm=nothing,
    tol=nothing,
    basis=nothing,
    endpoint_basis=nothing,
    ordinary_basis=nothing,
    sfe=false,
    lfe=false,
)
    h_order = Int(horizon_order === nothing ? (N === nothing ? DEFAULT_HORIZON_ORDER : N) : horizon_order)
    o_order = Int(ordinary_order === nothing ? (N === nothing ? DEFAULT_ORDINARY_ORDER : N) : ordinary_order)
    h_order > 0 || throw(ArgumentError("horizon_order must be positive."))
    o_order > 0 || throw(ArgumentError("ordinary_order must be positive."))

    match_x = Float64(xm === nothing ? DEFAULT_MATCH_X : xm)
    0 < match_x < 1 || throw(ArgumentError("xm must lie in (0, 1)."))

    tolerance = Float64(tol === nothing ? DEFAULT_TOLERANCE : tol)
    tolerance > 0 || throw(ArgumentError("tol must be positive."))

    endpoint = endpoint_basis === nothing ? (horizon=:hc, infinity=:ic) : endpoint_basis
    ordinary = ordinary_basis === nothing ? :auto : Symbol(ordinary_basis)
    selected_basis = basis === nothing ? :stable_default : Symbol(basis)
    frequency_selection = direct_frequency_regime(params; sfe, lfe)
    if N === nothing && ordinary_order === nothing &&
            frequency_selection.regime == :ordinary
        o_order = max(o_order, DEFAULT_REGULAR_ORDINARY_ORDER)
    end
    use_sfe = frequency_selection.regime == :sfe
    use_lfe = frequency_selection.regime == :lfe

    control_source = xm === nothing && basis === nothing &&
        endpoint_basis === nothing && ordinary_basis === nothing ? :default : :user

    return DirectGSNControls(
        h_order,
        o_order,
        Int(infinity_order === nothing ? (N === nothing ? DEFAULT_INFINITY_ORDER : N) : infinity_order),
        match_x,
        tolerance,
        selected_basis,
        endpoint,
        ordinary,
        use_sfe,
        use_lfe,
        frequency_selection,
        control_source,
    )
end

direct_gsn_controls(s::Integer, l::Integer, m::Integer, a, omega; kwargs...) =
    direct_gsn_controls(direct_gsn_parameters(s, l, m, a, omega); kwargs...)

end
