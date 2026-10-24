# Arbitrary-precision ISEM route of the GSN equation.
#
# For mode parameters beyond binary64 (a BigFloat spin, frequency or separation constant) the IN and
# UP solutions are built from the same ingredients as the Float64 routes -- the exact generated
# coefficient tables, the horizon Frobenius series, the asymptotic series at infinity and Taylor
# patches between them, matched at an interior point of the compact coordinate x -- with every number
# in the floating type R of the inputs. The series orders, the seed radii and every truncation test
# follow the requested tolerance (default: 64 ulps of the inputs' precision), and the BigFloat
# working precision is raised internally by a guard and by the documented continuation loss.
#
# The regimes of the Float64 routes that rest on binary64 data (large-frequency tables, near-extremal
# endpoint factors, the eikonal patches and the rotated complex contour) are not used here. Complex
# frequencies are seeded by the asymptotic series at infinity, whose error law is recorded in the
# diagnostics (`complex_loss_digits`); the MST representation of the infinity side replaces it in
# ManoSuzukiTakasugiGeneric.
module ArbitraryPrecisionRoute

using ...Coordinates: r_from_rstar, rstar_from_r, horizon_distance_from_rstar
using ...ConversionFactors
using ...Solutions: Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix
using ..Parameters: GSNParameters, GSNControls, gsn_controls, is_binary64_parameters
using ..Parameters: DEFAULT_TOLERANCE, DEFAULT_MATCH_X
using ..CoefficientTables: CoefficientSet, gsn_coefficients
using ..OrdinaryPointExpansion: ab_term_counts, ordinary_ab_series!, ordinary_point_coeffs!, poly_pair
using ..LocalSolutionAtZero: zero_local_solution, is_horizon_threshold
using ..ApparentSingularities: apparent_points, apparent_point_in, apparent_basis, apparent_jump, next_apparent,
    APPARENT_STEP_FRACTION
using ..LocalSolutionAtInfinity: infinity_local_solution
using ..Matching: ConjugatedRoute
using ..RadialFunctions: TeukolskySolution, TeukolskyRadialEvaluator, TeukolskyTupleEvaluator,
    _teukolsky_converter
import ..RadialFunctions: gsn_radial_function, teukolsky_radial_function, y_radial_function,
    _teukolsky_pair_with_error
using ..TeukolskyTransformation: d2R, GSN_to_Y_solution_from_matrix
using ..ManoSuzukiTakasugiInfinity: mano_suzuki_takasugi_nu, mano_suzuki_takasugi_nu_complex
using ..ManoSuzukiTakasugiGeneric: MSTGenericData, mst_generic_data, mst_nu_generic,
    mst_up_state, mst_down_state
using StaticArrays: SA

export ArbitraryPrecisionRadialRoute, arbitrary_precision_route
export route_state, route_state_r

# Default series tolerance in ulps of the inputs' precision, Taylor order per decimal digit, guard
# bits added to the working precision, and the largest seed coordinates at the two ends.
const TOLERANCE_ULPS = 64
const ORDER_PER_DIGIT = 2.2
const MINIMUM_ORDER = 24
const GUARD_BITS = 24
const MATCHING_GUARD_BITS = 8
# The MST side (nu iteration, series coefficients, recurrences) costs about two decimal digits
# relative to the asymptotic series at the same working precision; this covers them.
const MST_GUARD_BITS = 16
const MAX_REBUILDS = 4
const HORIZON_SEED_XMAX = 1 // 2
const INFINITY_SEED_YMAX = 1 // 2
const MAX_PATCHES = 2000

_root_module() = parentmodule(parentmodule(@__MODULE__))
_isem_module() = parentmodule(@__MODULE__)

_real_type(::GSNParameters{R}) where {R} = R
_complex_type(::GSNParameters{R}) where {R} = Complex{R}

# Compact coordinate x = (r - r_+) / (r - r_-), r_+- = 1 +- kappa; dx/dr* = 2 kappa x / (r^2 + a^2).
r_plus(p::GSNParameters) = one(p.kappa) + p.kappa
r_minus(p::GSNParameters) = one(p.kappa) - p.kappa
x_from_r(p::GSNParameters, r) = (r - r_plus(p)) / (r - r_minus(p))
r_from_x(p::GSNParameters, x) = (r_plus(p) - r_minus(p) * x) / (one(x) - x)
dx_drstar(p::GSNParameters, x, r) = 2 * p.kappa * x / (r^2 + p.a^2)
horizon_wavenumber(p::GSNParameters) = p.omega - p.m * p.a / (2 * r_plus(p))

# Constant phases of the unit-transmission normalization (Iteration.endpoint_seed_scale in R):
# horizon in/out = exp(-/+ i k (r_+ + 2 log kappa)) with k = omega - m a / (2 r_+),
# infinity up/down = exp(+/- i omega (r_- + 2 log kappa)).
function _horizon_scale(p::GSNParameters, branch::Symbol)
    phase = horizon_wavenumber(p) * (r_plus(p) + 2 * log(p.kappa))
    return branch == :in ? exp(-im * phase) : exp(im * phase)
end
function _infinity_scale(p::GSNParameters, branch::Symbol)
    phase = p.omega * (r_minus(p) + 2 * log(p.kappa))
    return branch == :up ? exp(im * phase) : exp(-im * phase)
end

# Taylor patch about x0 valid between x0 and xend.
struct SeriesPatch{R,C}
    x0::R
    xend::R
    coeffs::Vector{C}
end

# Horizon Frobenius series X = scale x^mu sum c_n x^n, used for x <= x_seed.
struct HorizonExpansion{R,C}
    exponent::C
    scale::C
    coeffs::Vector{C}
    x_seed::R
end

# Representations of the two infinity branches (:up outgoing, :down ingoing).
abstract type InfinityBranch end

# Asymptotic series at infinity X = scale exp(sigma / y) y^rho sum d_n y^n, y = 1 - x, used for
# y <= y_seed (optimally truncated at the seed); Taylor patches continue it inward.
struct InfinityExpansion{R,C} <: InfinityBranch
    rho::C
    sigma::C
    scale::C
    coeffs::Vector{C}
    y_seed::R
end

# MST representation (ManoSuzukiTakasugiGeneric): the Teukolsky UP or DOWN solution, convergent at
# every r > r_+, converted to the GSN state at the requested point. `gsn_scale` is the Teukolsky
# amplitude of the unit-GSN-transmission branch (ConversionFactors.Ctrans for UP, Binc for DOWN).
struct MSTInfinityBranch{D<:MSTGenericData,C} <: InfinityBranch
    data::D
    branch::Symbol
    gsn_scale::C
end

"""
    ArbitraryPrecisionRadialRoute

IN or UP solution of the GSN equation in the floating type of the mode parameters, in the GSN
unit-transmission normalization. Stores the endpoint expansions, the Taylor patches to the matching
point and the connection coefficients: `incidence`, `reflection` are B_inc, B_ref for IN
(X_in = B_inc X_down + B_ref X_up beyond the matching point) and C_inc, C_ref for UP
(X_up = C_inc X_out + C_ref X_in inside it). `route_state(route, x)` evaluates (X, dX/dx).
"""
struct ArbitraryPrecisionRadialRoute{P<:GSNParameters,R,C}
    branch::Symbol
    params::P
    controls::GSNControls
    coefficients::Any           # CoefficientSet of the mode, for extending the MST chains on demand
    extension_lock::ReentrantLock
    order::Int
    tolerance::R
    working_bits::Int
    match_x::R
    horizon_in::HorizonExpansion{R,C}
    horizon_in_patches::Vector{SeriesPatch{R,C}}
    horizon_out::Union{Nothing,HorizonExpansion{R,C}}
    horizon_out_patches::Vector{SeriesPatch{R,C}}
    infinity_up::InfinityBranch
    infinity_up_patches::Vector{SeriesPatch{R,C}}
    infinity_down::Union{Nothing,InfinityBranch}
    infinity_down_patches::Vector{SeriesPatch{R,C}}
    incidence::C
    reflection::C
    diagnostics::NamedTuple
end

function _normalized_branch(branch)
    branch in (:IN, :in, :ingoing) && return :IN
    branch in (:UP, :up, :outgoing) && return :UP
    throw(ArgumentError("the arbitrary-precision ISEM route builds :IN or :UP; DOWN and OUT are combined from them."))
end

# Precision of the inputs in bits (BigFloat values carry their own; other types their fixed one).
_input_bits(p::GSNParameters{BigFloat}) = max(precision(p.a), precision(real(p.omega)), precision(imag(p.omega)))
_input_bits(p::GSNParameters{R}) where {R} = precision(R)

# Smallest relative spacing at `bits` bits in the type R.
_ulp(::Type{BigFloat}, bits::Int) = ldexp(one(BigFloat), 1 - bits)
_ulp(::Type{R}, bits::Int) where {R} = eps(R)

# Bits of working precision consumed by the cancellation of a Taylor step: the tail criterion lets a
# step reach k h ~ N tol^(1/N) / e for a solution of local wavenumber k, whose series terms grow to
# about exp(k h) before they cancel to the value. The endpoint series are seeded where their own
# cancellation stays below the same bound.
function _taylor_cancellation_bits(order::Int, tolerance)
    return ceil(Int, order * Float64(tolerance)^(1 / order) / (ℯ * log(2)))
end

# The parameters with their derived quantities (kappa, nu) recomputed at the current precision. The
# coefficient tables take nu, the local solutions kappa and a; rounded independently at the input
# precision they describe slightly different equations, and every solution carries a generic error of
# that size, which the matching amplifies by its cancellation. a, omega and lambda keep their values.
_working_parameters(p::GSNParameters) = p
function _working_parameters(p::GSNParameters{BigFloat})
    a = p.a
    kappa = sqrt((one(a) - a) * (one(a) + a))
    return GSNParameters(p.s, p.l, p.m, a, p.omega, p.lambda, atan(a, kappa), kappa)
end

# sum |terms| / |sum| of a series evaluated at h, as a Float64 (1 when there is none).
function _series_cancellation(coeffs, h, value)
    R = typeof(h)
    total = zero(R)
    hk = one(R)
    @inbounds for c in coeffs
        total += abs(c) * hk
        hk *= h
    end
    return Float64(max(total / max(abs(value), floatmin(R)), one(R)))
end

# Bits of working precision lost by continuing the inward-decaying branch from the seed radius r0 to
# the matching radius r_m at complex frequency: exp(2 |Im omega| (r0 - r_m)); the seed radius of the
# asymptotic series is about ln(1 / tol) / (2 |omega|) (Sec. 5.3 of the design note).
function _complex_seed_radius(p::GSNParameters, tol)
    return log(1 / Float64(tol)) / (2 * Float64(abs(p.omega)))
end

"""
    arbitrary_precision_route(params, branch; N=nothing, xm=nothing, tol=nothing, infinity=:auto)

IN or UP route in the floating type of `params` (see `ArbitraryPrecisionRadialRoute`). `tol` is the
relative truncation tolerance of every series (default: $(TOLERANCE_ULPS) ulps of the inputs'
precision), `N` the Taylor order of the patches (default: $(ORDER_PER_DIGIT) per requested decimal
digit, at least $(MINIMUM_ORDER); the endpoint series use 2N and 3N terms), `xm` the matching point in
x (default 1/2). `infinity` selects the representation of the infinity side: `:series` (asymptotic
series seeds continued by Taylor patches; exact only for real frequencies), `:mst` (the convergent
MST series of ManoSuzukiTakasugiGeneric, evaluated at every requested point), or `:auto` (`:mst`,
falling back to `:series` at a real frequency if the MST data cannot be built). Negative real frequencies are built from the mirror mode
(m -> -m, omega -> -omega) and returned as a `ConjugatedRoute`. `max_rebuilds` caps the rebuilds that
cover the measured cancellation of the matching (default $(MAX_REBUILDS)). `root_search=true` leaves
the cancellation of the incidence numerator out of that measure: a root search evaluates the
incidence amplitude where it vanishes by nature, and needs it accurate in absolute terms only.
"""
function arbitrary_precision_route(params::GSNParameters, branch=:IN; N=nothing, xm=nothing, tol=nothing, controls=nothing, infinity::Symbol=:auto,
        max_rebuilds::Integer=MAX_REBUILDS, root_search::Bool=false)
    route_branch = _normalized_branch(branch)
    iszero(params.omega) && throw(ArgumentError("the arbitrary-precision ISEM route requires a nonzero frequency."))
    infinity in (:auto, :series, :mst) || throw(ArgumentError("infinity must be :auto, :series or :mst."))
    # :auto takes the MST side at every frequency for BigFloat (since the UP-series factors are stored,
    # it is faster than the asymptotic series with its inward Taylor chain, 1.5-4 times at 256 bits, with
    # the same digits); at a real frequency the asymptotic series remains the fallback should the MST
    # data fail. For a fixed-precision type (Double64) the MST data live in BigFloat and the asymptotic
    # series, in the type itself, is 20 times faster, so it is taken at real frequencies.
    R = _real_type(params)
    infinity_mode = infinity === :auto ? ((R === BigFloat || !isreal(params.omega)) ? :mst : :series) : infinity
    if infinity === :auto && isreal(params.omega) && R === BigFloat
        try
            return arbitrary_precision_route(params, branch; N, xm, tol, controls, infinity=:mst, max_rebuilds, root_search)
        catch err
            err isa ArgumentError && rethrow(err)
            @debug "arbitrary-precision route: MST side failed at a real frequency, using the asymptotic series" exception=err
            return arbitrary_precision_route(params, branch; N, xm, tol, controls, infinity=:series, max_rebuilds, root_search)
        end
    end
    # Re omega < 0 (real or complex): the principal-branch MST formulas are not on the physical sheet
    # there, so the route is built for the mirror mode (m, omega) -> (-m, -conj(omega)), whose
    # separation constant is conj(lambda), and complex-conjugated on evaluation.
    if real(params.omega) < 0
        positive = GSNParameters(params.s, params.l, -params.m, params.a, -conj(params.omega),
            conj(params.lambda), params.nu, params.kappa)
        return ConjugatedRoute(
            arbitrary_precision_route(positive, route_branch; N, xm, tol, controls, infinity, max_rebuilds, root_search), params)
    end
    if controls !== nothing
        # Controls built by gsn_controls carry binary64 defaults; only explicit values are taken.
        tol === nothing && controls.tolerance != DEFAULT_TOLERANCE && (tol = controls.tolerance)
        xm === nothing && controls.match_x != DEFAULT_MATCH_X && (xm = controls.match_x)
    end
    R = _real_type(params)
    bits = _input_bits(params)
    tolerance = tol === nothing ? TOLERANCE_ULPS * _ulp(R, bits) : R(tol)
    tolerance > 0 || throw(ArgumentError("tol must be positive."))
    digits = -log10(Float64(tolerance))
    order = N === nothing ? max(MINIMUM_ORDER, ceil(Int, ORDER_PER_DIGIT * digits)) : Int(N)
    order > 0 || throw(ArgumentError("N must be positive."))
    loss_bits = 0
    if !isreal(params.omega) && infinity_mode == :series
        r0 = _complex_seed_radius(params, tolerance)
        loss_bits = ceil(Int, Float64(abs(imag(params.omega))) * r0 / log(2))
    end
    tolerance_bits = ceil(Int, -log2(Float64(tolerance)))
    extra_bits = GUARD_BITS + loss_bits + (infinity_mode == :mst ? MST_GUARD_BITS : 0)
    working_bits = R === BigFloat ? max(bits, tolerance_bits + extra_bits + _taylor_cancellation_bits(order, tolerance)) : precision(R)
    if R !== BigFloat
        # A fixed-precision type is built once at its own precision; when the measured cancellation of
        # the matching or of the series sums exceeds the guard, the route is rebuilt in BigFloat at the
        # type's precision with the full loop below, and that route is returned (its values are BigFloat).
        route = _build_route(params, route_branch, order, tolerance, xm, working_bits, infinity_mode)
        d = route.diagnostics
        matching = root_search ? max(d.determinant_condition, d.reflection_condition) : d.matching_condition
        lost_bits = ceil(Int, log2(matching)) +
            max(0, ceil(Int, log2(d.mst_condition)) - (infinity_mode == :mst ? MST_GUARD_BITS : 0))
        series_bits = ceil(Int, log2(d.series_condition))
        # Without guard bits every bit of cancellation is a bit of the result: the type's slack is the
        # 64-ulp tolerance plus half the guard, beyond which BigFloat takes over.
        lost_bits + series_bits <= GUARD_BITS ÷ 2 && return route
        return arbitrary_precision_route(_bigfloat_parameters(params, precision(R)), route_branch;
            N, xm, tol, controls, infinity, max_rebuilds, root_search)
    end
    # The connection coefficients come from 2 x 2 determinants whose terms cancel (by
    # exp(2 |Im omega| r*_m) at complex frequency, by |B_inc / B_ref| when the reflection is
    # exponentially small), and the MST states cancel in their sums and in the Teukolsky -> GSN
    # conversion. The cancellation amplifies the truncation errors of the series as well as the
    # rounding errors of their sums, so the route is rebuilt with the series tolerance tightened by the
    # measured bits (plus MATCHING_GUARD_BITS) and the working precision raised to cover the measured
    # cancellation of the sums, until both are covered. A build cannot measure a cancellation beyond
    # its own noise (about 1 / tol), so the measurement is repeated after each rebuild.
    tolerance_k, order_k, working_k = tolerance, order, working_bits
    route = setprecision(BigFloat, working_k) do
        _build_route(params, route_branch, order_k, tolerance_k, xm, working_k, infinity_mode)
    end
    for rebuild in 1:max_rebuilds
        d = route.diagnostics
        # The MST sums are already truncated MST_GUARD_BITS below the tolerance (and summed at a
        # precision raised by as much), so only their cancellation beyond that costs anything.
        # A root search on the incidence amplitude evaluates it where it vanishes by nature; there the
        # cancellation of its numerator is not a loss, and only the determinant's and the other
        # coefficient's are covered.
        matching = root_search ? max(d.determinant_condition, d.reflection_condition) : d.matching_condition
        lost_bits = ceil(Int, log2(matching)) +
            max(0, ceil(Int, log2(d.mst_condition)) - (infinity_mode == :mst ? MST_GUARD_BITS : 0))
        series_bits = ceil(Int, log2(d.series_condition))
        tolerance_bits_k = ceil(Int, -log2(Float64(tolerance_k)))
        tolerance_ok = tolerance_bits_k >= tolerance_bits + lost_bits - MATCHING_GUARD_BITS
        rounding_ok = working_k >= tolerance_bits_k + series_bits - MATCHING_GUARD_BITS
        tolerance_ok && rounding_ok && break
        if !tolerance_ok
            tolerance_k = tolerance * ldexp(one(R), -(lost_bits + MATCHING_GUARD_BITS))
            order_k = N === nothing ? max(MINIMUM_ORDER, ceil(Int, ORDER_PER_DIGIT * -log10(Float64(tolerance_k)))) : order
            tolerance_bits_k = ceil(Int, -log2(Float64(tolerance_k)))
        end
        working_k = max(bits, working_k + 1,
            tolerance_bits_k + extra_bits + max(_taylor_cancellation_bits(order_k, tolerance_k), series_bits))
        route = setprecision(BigFloat, working_k) do
            _build_route(params, route_branch, order_k, tolerance_k, xm, working_k, infinity_mode; rebuilds=rebuild)
        end
    end
    return route
end

# The matching point stays at x = 1/2 for the MST side too: at complex frequency the two infinity
# branches differ in size by exp(2 |Im omega| r*) at the matching radius, and extracting the smaller
# coefficient loses that factor, which grows with r_m (the MST sums need more terms at small zhat,
# but the recurrence makes them cheap). The remaining loss is measured and absorbed by the working
# precision (see `arbitrary_precision_route`).
function _default_match_x(p::GSNParameters, tolerance, infinity_mode::Symbol)
    R = _real_type(p)
    (isreal(p.omega) || infinity_mode == :mst) && return R(1 // 2)
    # Asymptotic seed at complex frequency: match at half the seed radius so that both sides share
    # the continuation loss.
    r0 = _complex_seed_radius(p, tolerance)
    return max(R(1 // 2), x_from_r(p, R(max(r0 / 2, Float64(r_plus(p)) + 2))))
end

# Binary64 seeds of the renormalized angular momentum from the Float64 MST module, in the order to try:
# at the frequency itself, then (the Float64 routine rejects some frequencies within 1e-15 of the
# positive imaginary axis) at the frequency with its real part moved to 1e-6 and 1e-3 of |omega|,
# which changes nu by that much and still seeds the generic iteration; `l` last.
function _nu_seeds(p::GSNParameters)
    a = Float64(p.a)
    omega = ComplexF64(p.omega)
    lambda = ComplexF64(p.lambda)
    seeds = ComplexF64[]
    shifts = isreal(p.omega) ? (omega,) :
        (omega, ComplexF64(1e-6 * abs(omega), imag(omega)), ComplexF64(1e-3 * abs(omega), imag(omega)))
    for w in shifts
        try
            seed = isreal(p.omega) ? ComplexF64(mano_suzuki_takasugi_nu(p.s, p.l, p.m, a, real(w), lambda)) :
                ComplexF64(mano_suzuki_takasugi_nu_complex(p.s, p.l, p.m, a, w, lambda))
            isfinite(seed) && push!(seeds, seed)
        catch
        end
    end
    push!(seeds, ComplexF64(p.l))
    return seeds
end
_nu_seed(p::GSNParameters) = first(_nu_seeds(p))

# MST data and the two infinity branches in GSN unit-transmission normalization. The MST sums are
# truncated MST_GUARD_BITS below the route tolerance: their Teukolsky states pass through the inverse
# conversion matrix, whose conditioning would otherwise cost about two digits (the working precision
# is raised by the same bits).
# The parameters in BigFloat at `bits` (exact values of a, omega, lambda; kappa and nu recomputed).
function _bigfloat_parameters(p::GSNParameters, bits::Int)
    return setprecision(BigFloat, bits) do
        omega = isreal(p.omega) ? BigFloat(real(p.omega)) : Complex{BigFloat}(p.omega)
        lambda = isreal(p.lambda) ? BigFloat(real(p.lambda)) : Complex{BigFloat}(p.lambda)
        _working_parameters(GSNParameters(p.s, p.l, p.m, BigFloat(p.a), omega, lambda, BigFloat(p.nu), BigFloat(p.kappa)))
    end
end

# Precision of the MST data of a branch whose GSN states are of type C: that of the type for BigFloat
# routes (set by the caller), the type's precision plus the MST guard for fixed-precision types, whose
# MST data live in BigFloat.
_mst_bits(::Type{Complex{BigFloat}}) = precision(BigFloat)
_mst_bits(::Type{C}) where {C<:Complex} = precision(real(C)) + MST_GUARD_BITS

function _mst_branches(p::GSNParameters, tolerance)
    R = _real_type(p)
    if R !== BigFloat
        # A fixed-precision type cannot carry the guard bits of the MST sums: the data are built in
        # BigFloat at the type's precision plus the guard, and the states are converted on evaluation.
        C = _complex_type(p)
        data, up, down = setprecision(BigFloat, _mst_bits(C)) do
            _mst_branches(_bigfloat_parameters(p, _mst_bits(C)), BigFloat(tolerance))
        end
        return data, MSTInfinityBranch(data, :up, C(up.gsn_scale)), MSTInfinityBranch(data, :down, C(down.gsn_scale))
    end
    mst_tolerance = tolerance * ldexp(one(tolerance), -MST_GUARD_BITS)
    seeds = _nu_seeds(p)
    local nu, residual, depth
    for (k, seed) in enumerate(seeds)
        try
            nu, residual, depth = mst_nu_generic(p, seed; tol=mst_tolerance)
            break
        catch err
            k == length(seeds) && rethrow(err)
        end
    end
    data = mst_generic_data(p, nu; tol=mst_tolerance, nu_residual=residual, nu_depth=depth)
    # (the conversion factors are the integer 1 for s = 0: the scale is held in the complex working type)
    C = _complex_type(p)
    up = MSTInfinityBranch(data, :up, C(ConversionFactors.Ctrans(p.s, p.m, p.a, p.omega, p.lambda)))
    down = MSTInfinityBranch(data, :down, C(ConversionFactors.Binc(p.s, p.m, p.a, p.omega, p.lambda)))
    return data, up, down
end

function _build_route(p::GSNParameters, branch::Symbol, order::Int, tolerance, xm, working_bits::Int, infinity_mode::Symbol;
        rebuilds::Int=0)
    p = _working_parameters(p)
    R = _real_type(p)
    C = _complex_type(p)
    cancellation_limit = 2.0^_taylor_cancellation_bits(order, tolerance)
    controls = gsn_controls(p; N=order, tol=Float64(tolerance),
        xm=xm === nothing ? nothing : Float64(xm),
        small_frequency_expansion=false, large_frequency_expansion=false)
    coefficients = gsn_coefficients(p; controls=controls)
    match_x = xm === nothing ? _default_match_x(p, tolerance, infinity_mode) : R(xm)
    0 < match_x < 1 || throw(ArgumentError("xm must lie in (0, 1)."))
    work = _Workspace(C, coefficients, order)

    t_horizon = time()
    horizon_in = _horizon_expansion(coefficients, :in, 2order, tolerance, cancellation_limit)
    Xh, dXh = _horizon_state(horizon_in, horizon_in.x_seed)
    horizon_in_patches, Xin, dXin, chain_condition = _continue(work, coefficients, horizon_in.x_seed, Xh, dXh, match_x, order, tolerance)
    series_condition = max(_seed_cancellation(horizon_in), chain_condition)
    t_horizon = time() - t_horizon

    t_infinity = time()
    mst = nothing
    infinity_up_patches = SeriesPatch{R,C}[]
    infinity_down_patches = SeriesPatch{R,C}[]
    local infinity_up, Xup, dXup, mst_down
    mst_condition = 1.0
    if infinity_mode == :mst
        mst, infinity_up, mst_down = _mst_branches(p, tolerance)
        Xup, dXup, up_condition = _mst_gsn_state(infinity_up, match_x)
        mst_condition = max(mst_condition, up_condition)
    else
        infinity_up = _infinity_expansion(coefficients, :up, 3order, tolerance, cancellation_limit)
        Xu, dXu = _infinity_state(infinity_up, one(R) - infinity_up.y_seed)
        infinity_up_patches, Xup, dXup, chain_condition = _continue(work, coefficients, one(R) - infinity_up.y_seed, Xu, dXu, match_x, order, tolerance)
        series_condition = max(series_condition, _seed_cancellation(infinity_up), chain_condition)
    end

    horizon_out = nothing
    horizon_out_patches = SeriesPatch{R,C}[]
    infinity_down = nothing
    if branch == :IN
        if infinity_mode == :mst
            infinity_down = mst_down
            Xdn, dXdn, down_condition = _mst_gsn_state(infinity_down, match_x)
            mst_condition = max(mst_condition, down_condition)
        else
            infinity_down = _infinity_expansion(coefficients, :down, 3order, tolerance, cancellation_limit)
            Xd, dXd = _infinity_state(infinity_down, one(R) - infinity_down.y_seed)
            infinity_down_patches, Xdn, dXdn, chain_condition = _continue(work, coefficients, one(R) - infinity_down.y_seed, Xd, dXd, match_x, order, tolerance)
            series_condition = max(series_condition, _seed_cancellation(infinity_down), chain_condition)
        end
        t_infinity = time() - t_infinity
        # X_in = B_inc X_down + B_ref X_up at the matching point.
        det = Xdn * dXup - Xup * dXdn
        incidence = (Xin * dXup - Xup * dXin) / det
        reflection = (Xdn * dXin - Xin * dXdn) / det
        conditions = (_cancellation(Xdn * dXup, Xup * dXdn), _cancellation(Xin * dXup, Xup * dXin),
            _cancellation(Xdn * dXin, Xin * dXdn))
    else
        t_infinity = time() - t_infinity
        is_horizon_threshold(p) && throw(ArgumentError(
            "the arbitrary-precision UP route does not build the logarithmic horizon companion at the superradiant threshold."))
        t_out = time()
        horizon_out = _horizon_expansion(coefficients, :out, 2order, tolerance, cancellation_limit)
        Xo, dXo = _horizon_state(horizon_out, horizon_out.x_seed)
        horizon_out_patches, Xout, dXout, chain_condition = _continue(work, coefficients, horizon_out.x_seed, Xo, dXo, match_x, order, tolerance)
        series_condition = max(series_condition, _seed_cancellation(horizon_out), chain_condition)
        t_horizon += time() - t_out
        # X_up = C_inc X_out + C_ref X_in at the matching point.
        det = Xout * dXin - Xin * dXout
        incidence = (Xup * dXin - Xin * dXup) / det
        reflection = (Xout * dXup - Xup * dXout) / det
        conditions = (_cancellation(Xout * dXin, Xin * dXout), _cancellation(Xup * dXin, Xin * dXup),
            _cancellation(Xout * dXup, Xup * dXout))
    end

    # Evaluation outside the matching point: the MST branches are continued outward by Taylor patches
    # on demand (_extend_mst_chain!), at a precision raised by the growth exp(2 |Im omega| (r*_end - r*_m))
    # of the branch that decays outward (the admixture of the growing branch that rounding leaves), so
    # that both chains carry the tolerance. Nothing is built here: the amplitudes and the horizon side
    # are what most callers need, and a chain costs as much as the rest of the route.

    r_m = r_from_x(p, match_x)
    series = infinity_up isa InfinityExpansion
    r0 = series ? Float64(2 * p.kappa / infinity_up.y_seed) : Float64(r_m)
    loss = (isreal(p.omega) || !series) ? 0.0 :
        Float64(2 * abs(imag(p.omega)) * max(r0 - Float64(r_m), 0) / log(10))
    diagnostics = (
        working_bits = working_bits,
        input_bits = _input_bits(p),
        tolerance = Float64(tolerance),
        order = order,
        horizon_order = 2order,
        infinity_order = series ? 3order : 0,
        match_x = Float64(match_x),
        r_m = Float64(r_m),
        horizon_seed_x = Float64(horizon_in.x_seed),
        infinity_seed_y = series ? Float64(infinity_up.y_seed) : 0.0,
        infinity_seed_r = r0,
        infinity_terms = series ? length(infinity_up.coeffs) : (mst.nmax - mst.nmin + 1),
        horizon_patches = length(horizon_in_patches) + length(horizon_out_patches),
        infinity_patches = length(infinity_up_patches) + length(infinity_down_patches),
        infinity_chain_end_x = isempty(infinity_up_patches) ? Float64(match_x) : Float64(maximum(max(patch.x0, patch.xend) for patch in infinity_up_patches)),
        matching_determinant = det,
        matching_condition = maximum(conditions),
        determinant_condition = conditions[1],
        incidence_condition = conditions[2],
        reflection_condition = conditions[3],
        mst_condition = mst_condition,
        series_condition = series_condition,
        rebuilds = rebuilds,
        complex_loss_digits = loss,
        frequency_regime = :arbitrary_precision_series,
        infinity_representation = series ? :asymptotic_series : :mst,
        horizon_seconds = t_horizon,
        infinity_seconds = t_infinity,
        mst_nu = mst === nothing ? missing : ComplexF64(mst.nu),
        mst_nu_residual = mst === nothing ? missing : Float64(mst.nu_residual),
        mst_depth = mst === nothing ? missing : mst.nu_depth,
    )
    return ArbitraryPrecisionRadialRoute{typeof(p),R,C}(
        branch, p, controls, coefficients, ReentrantLock(), order, tolerance, working_bits, match_x,
        horizon_in, horizon_in_patches, horizon_out, horizon_out_patches,
        infinity_up, infinity_up_patches, infinity_down, infinity_down_patches,
        incidence, reflection, diagnostics)
end

# Guard bits of the outward MST chains above the loss of the decaying branch along them.
const MST_CHAIN_GUARD_BITS = 16
# Reach of the outward chains: a chain costs a patch per fraction of a wavelength and, at complex
# frequency, 2 |Im omega| / ln 2 extra bits per unit of r*; beyond MST_CHAIN_MAX_PHASE radians of phase
# from the matching radius or MST_CHAIN_MAX_LOSS_BITS extra bits the MST series are summed directly at the
# requested point instead (about 0.1 s per point at 256 bits, against minutes for a chain to r = 1000 at
# omega = 0.5 - 0.3i, whose precision would grow by 900 bits).
const MST_CHAIN_MAX_PHASE = 60.0
const MST_CHAIN_MAX_LOSS_BITS = 64.0
# Segments of an UP chain before the points beyond are summed directly instead (each starts from an MST sum).
const MST_CHAIN_MAX_UP_SEGMENTS = 16

# The raised precision of a chain covers its own rounding, but not the error of its starting state, the
# MST state (truncated MST_GUARD_BITS below the route tolerance), which grows along the chain with the
# ratio of the other solution to the branch. In the IN solution, B_inc X_down + B_ref X_up, that admixture
# keeps its size relative to the solution (the solution is the regular one in the centrifugal region and
# holds the growing branch at complex frequency), so an IN chain is bounded by its rounding only. The UP
# solution is the branch itself: its chain is cut into segments that each start from the MST state at
# their first point and lose at most MST_GUARD_BITS (`_up_chain_loss_bits`). A single chain from the
# matching radius kept 18 of the 34 digits at 113 bits for s = 1, l = 4, omega = 1/20 at r = 1000.
function _chain_worthwhile(route::ArbitraryPrecisionRadialRoute, x)
    p = route.params
    r_m = r_from_x(p, route.match_x)
    r = r_from_x(p, x)
    Float64(abs(p.omega)) * Float64(r - r_m) <= MST_CHAIN_MAX_PHASE || return false
    route.branch == :UP && return _up_chain_segments(p, r_m, r) <= MST_CHAIN_MAX_UP_SEGMENTS
    return _chain_loss_bits(p, r_m, r) <= MST_CHAIN_MAX_LOSS_BITS
end

# Bits lost by continuing a branch outward from r0 to r1: the branch that decays outward picks up the
# growing one through rounding, so its relative error grows by the ratio of the two, at complex frequency
# exp(2 |Im omega| (r*_1 - r*_0)).
function _chain_loss_bits(p::GSNParameters, r0, r1)
    damping = Float64(abs(imag(p.omega)))
    iszero(damping) && return 0.0
    a = Float64(p.a)
    return 2 * damping * (rstar_from_r(a, Float64(r1)) - rstar_from_r(a, Float64(r0))) / log(2)
end

# Growth of the other solution relative to the UP branch from r0 to r1: exp(2 Im omega (r*_1 - r*_0)) at
# Im omega > 0 only (the branch grows outward at Im omega < 0); in the centrifugal region
# r < (l + 1/2) / |omega| the static ratio v^(2l+1) in v = (r - r_+) / (r_+ - r_-) (the UP solution is the
# irregular one there; near the horizon of a = 0.99 the measured loss of a chain to three matching radii
# at l = 5, 35 bits, followed v, where r gives 17); beyond the turning point r^(2s) for s > 0 (the UP
# solution goes as r^(-1-2s) against r^(-1) in the Teukolsky function that R and Y are made of).
function _up_chain_loss_bits(p::GSNParameters, r0, r1)
    growth = Float64(imag(p.omega))
    loss = growth > 0 ? 2 * growth * (rstar_from_r(Float64(p.a), Float64(r1)) - rstar_from_r(Float64(p.a), Float64(r0))) / log(2) : 0.0
    rp = Float64(r_plus(p))
    turning = max((p.l + 0.5) / max(Float64(abs(p.omega)), floatmin(Float64)), rp)
    lo = min(Float64(r0), turning)
    hi = min(Float64(r1), turning)
    hi > lo && (loss += (2 * p.l + 1) * log2((hi - rp) / (lo - rp)))
    lo = max(Float64(r0), turning)
    hi = max(Float64(r1), turning)
    hi > lo && p.s > 0 && (loss += 2 * p.s * log2(hi / lo))
    return loss
end

# Largest radius in [r0, r1] whose loss from r0 stays within `budget` (the losses grow with r).
function _chain_reach(loss, r0, r1, budget)
    lo, hi = Float64(r0), Float64(r1)
    loss(hi) <= budget && return hi
    for _ in 1:48
        mid = (lo + hi) / 2
        loss(mid) <= budget ? (lo = mid) : (hi = mid)
    end
    return lo
end

# Number of UP segments from r0 to r1, each within MST_GUARD_BITS of loss.
function _up_chain_segments(p::GSNParameters, r0, r1)
    n = 0
    r = Float64(r0)
    while r < Float64(r1) && n <= MST_CHAIN_MAX_UP_SEGMENTS
        r = max(_chain_reach(x -> _up_chain_loss_bits(p, r, x), r, r1, MST_GUARD_BITS), nextfloat(r))
        n += 1
    end
    return n
end

# Extend the outward Taylor chain of an MST branch so that it covers x_target: from the end state of the
# chain (or the MST state at the matching point), to at least twice the present distance from the matching
# radius but no further than a chain may reach, at a precision raised by the growth of the decaying branch
# over the new stretch. An UP chain grows by segments that start from the MST state at their first point.
# The patch vectors of the route are appended under the route's lock.
function _extend_mst_chain!(route::ArbitraryPrecisionRadialRoute, which::Symbol, x_target)
    lock(route.extension_lock) do
        p = route.params
        R = _real_type(p)
        C = _complex_type(p)
        patches = which == :up ? route.infinity_up_patches : route.infinity_down_patches
        expansion = which == :up ? route.infinity_up : route.infinity_down
        segmented = route.branch == :UP
        r_m = r_from_x(p, route.match_x)
        r_target = r_from_x(p, R(x_target))
        for _ in 1:MST_CHAIN_MAX_UP_SEGMENTS + 1
            if !isempty(patches) && x_target <= _chain_range(patches)[2]
                return nothing                  # covered (also by another evaluation meanwhile)
            end
            if isempty(patches)
                x0 = route.match_x
                X0, dX0 = _infinity_state(expansion, x0)
            else
                last = patches[end]
                x0 = last.xend
                X0, dX0 = segmented ? _infinity_state(expansion, x0) : poly_pair(last.coeffs, x0 - last.x0)
            end
            r0 = r_from_x(p, x0)
            r_new = max(r_target, r_m + 2 * (r0 - r_m), r0 + one(R))
            # (no further than a chain or, for UP, the present segment may reach)
            phase_reach = Float64(r_m) + MST_CHAIN_MAX_PHASE / max(Float64(abs(p.omega)), floatmin(Float64))
            reach = segmented ?
                _chain_reach(r -> _up_chain_loss_bits(p, r0, r), r0, min(Float64(r_new), phase_reach), MST_GUARD_BITS) :
                _chain_reach(r -> _chain_loss_bits(p, r_m, r), r_target, min(Float64(r_new), max(phase_reach, Float64(r_target))), MST_CHAIN_MAX_LOSS_BITS)
            r_new = min(r_new, R(reach))
            # (ending exactly at the target when it is reached: the round trip through r can fall short of it)
            x_new = r_new >= r_target ? max(x_from_r(p, r_new), R(x_target)) : x_from_r(p, r_new)
            x_new > x0 || return nothing
            loss = ceil(Int, segmented ? _up_chain_loss_bits(p, r0, r_new) : _chain_loss_bits(p, r0, r_new)) + MST_CHAIN_GUARD_BITS
            build(coefficients, work) = _continue(work, coefficients, x0, X0, dX0, x_new, route.order, route.tolerance)[1]
            new_patches = if R === BigFloat
                setprecision(BigFloat, route.working_bits + loss) do
                    coefficients = gsn_coefficients(_working_parameters(p); controls=route.controls)
                    build(coefficients, _Workspace(C, coefficients, route.order))
                end
            else
                build(route.coefficients, _Workspace(C, route.coefficients, route.order))
            end
            append!(patches, new_patches)
        end
        return nothing
    end
end

# Cancellation of a difference u - v: (|u| + |v|) / |u - v|, as a Float64 (1 when there is none).
function _cancellation(u, v)
    R = real(typeof(u))
    return Float64(max((abs(u) + abs(v)) / max(abs(u - v), floatmin(R)), one(R)))
end

# Buffers of the ordinary-point expansions: packed P/Q values and the A, B Taylor coefficients.
struct _Workspace{C}
    pq::Vector{C}
    a::Vector{C}
    b::Vector{C}
end
function _Workspace(::Type{C}, coefficients::CoefficientSet, order::Int) where {C}
    a_terms, b_terms = ab_term_counts(:O, order)
    return _Workspace{C}(
        Vector{C}(undef, coefficients.ordinary.value_count),
        Vector{C}(undef, a_terms),
        Vector{C}(undef, b_terms))
end

# Relative size of the last four terms of a polynomial at distance h, against the value.
function _tail_small(coeffs, h, value, tolerance)
    R = typeof(h)
    tail = zero(R)
    hk = h^(length(coeffs) - 4)
    @inbounds for k in max(1, length(coeffs) - 3):length(coeffs)
        tail += abs(coeffs[k]) * hk
        hk *= h
    end
    return tail <= tolerance * max(abs(value), floatmin(R))
end

# Horizon Frobenius series and its seed point: the largest x <= HORIZON_SEED_XMAX at which the last
# terms are below the tolerance and the cancellation of the sum is below `cancellation_limit`.
function _horizon_expansion(coefficients::CoefficientSet, branch::Symbol, order::Int, tolerance, cancellation_limit)
    p = coefficients.params
    R = _real_type(p)
    local_solution = zero_local_solution(coefficients, branch, order)
    coeffs = local_solution.coefficients
    x = R(HORIZON_SEED_XMAX)
    for _ in 1:400
        value, _ = poly_pair(coeffs, x)
        _tail_small(coeffs, x, value, tolerance) && _series_cancellation(coeffs, x, value) <= cancellation_limit && break
        x *= R(4 // 5)
    end
    return HorizonExpansion{R,eltype(coeffs)}(local_solution.exponent, _horizon_scale(p, branch), coeffs, x)
end

_seed_cancellation(h::HorizonExpansion) = _series_cancellation(h.coeffs, h.x_seed, poly_pair(h.coeffs, h.x_seed)[1])
_seed_cancellation(e::InfinityExpansion) = _series_cancellation(e.coeffs, e.y_seed, poly_pair(e.coeffs, e.y_seed)[1])

# (X, dX/dx) of a horizon expansion at x.
function _horizon_state(h::HorizonExpansion, x)
    value, derivative = poly_pair(h.coeffs, x)
    prefactor = h.scale * x^h.exponent
    return prefactor * value, prefactor * (derivative + h.exponent * value / x)
end

# Asymptotic series at infinity, optimally truncated: the seed y is the largest y <= INFINITY_SEED_YMAX
# at which a term falls below the tolerance relative to the partial sum (with its predecessor within
# a factor 10 of it) and the cancellation of the sum is below `cancellation_limit`; the series is cut
# there.
function _infinity_expansion(coefficients::CoefficientSet, branch::Symbol, order::Int, tolerance, cancellation_limit)
    p = coefficients.params
    R = _real_type(p)
    local_solution = infinity_local_solution(coefficients, branch, order)
    d = local_solution.coefficients
    C = eltype(d)
    y = R(INFINITY_SEED_YMAX)
    nopt = length(d)
    for _ in 1:400
        partial = zero(C)
        total = zero(R)
        previous = zero(R)
        yk = one(R)
        found = false
        @inbounds for k in 1:length(d)
            term = abs(d[k]) * yk
            partial += d[k] * yk
            total += term
            scale = tolerance * max(abs(partial), floatmin(R))
            if k >= 3 && term <= scale && previous <= 10 * scale
                nopt = k
                found = total <= cancellation_limit * max(abs(partial), floatmin(R))
                break
            end
            previous = term
            yk *= y
        end
        found && break
        y *= R(4 // 5)
    end
    return InfinityExpansion{R,C}(local_solution.rho, local_solution.sigma,
        _infinity_scale(p, branch), d[1:nopt], y)
end

# (X, dX/dx) of an infinity expansion at x (y = 1 - x; dx = -dy).
function _infinity_state(e::InfinityExpansion, x)
    y = one(x) - x
    value, derivative = poly_pair(e.coeffs, y)
    prefactor = e.scale * exp(e.sigma / y) * y^e.rho
    dXdy = prefactor * (derivative + (-e.sigma / y^2 + e.rho / y) * value)
    return prefactor * value, -dXdy
end

# (X, dX/dx) of an MST branch at x: the Teukolsky state at r, scaled to the GSN unit-transmission
# normalization, and converted with the inverse of the GSN -> Teukolsky conversion matrix. The third
# value is the cancellation met on the way (the Tricomi sums and the conversion), by which the
# truncation error of the series is amplified in the GSN state.
function _mst_gsn_state(e::MSTInfinityBranch{<:MSTGenericData,C}, x) where {C}
    real(C) === BigFloat && return _mst_gsn_state_at_precision(e, x)
    X, dXdx, condition = setprecision(BigFloat, _mst_bits(C)) do
        _mst_gsn_state_at_precision(e, x)
    end
    return C(X), C(dXdx), condition
end
function _mst_gsn_state_at_precision(e::MSTInfinityBranch, x)
    p = e.data.params
    r = r_from_x(p, x)
    state = e.branch == :up ? mst_up_state(e.data, r) : mst_down_state(e.data, r)
    matrix = Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(
        p.s, p.m, p.a, p.omega, p.lambda, r)
    inverse = inv(matrix)
    R_ = e.gsn_scale * state.R
    dR_ = e.gsn_scale * state.dR
    X = inverse[1, 1] * R_ + inverse[1, 2] * dR_
    dXdrstar = inverse[2, 1] * R_ + inverse[2, 2] * dR_
    conversion = max(_cancellation(inverse[1, 1] * R_, -inverse[1, 2] * dR_),
        _cancellation(inverse[2, 1] * R_, -inverse[2, 2] * dR_))
    return X, dXdrstar / dx_drstar(p, x, r), state.cancellation * conversion
end
function _infinity_state(e::MSTInfinityBranch, x)
    X, dXdx, _ = _mst_gsn_state(e, x)
    return X, dXdx
end

# One Taylor patch from x0 (state X, dX) towards `target`: the step is at most half the empirical
# radius of the solution series (ratio of its last coefficients), then shrunk until the last four
# terms are below the tolerance relative to the value. Returns the patch, the new state and the
# cancellation sum |terms| / |value| of the step.
function _taylor_step(work::_Workspace, coefficients::CoefficientSet, x0, X, dX, target, order::Int, tolerance)
    R = typeof(x0)
    C = eltype(work.pq)
    ordinary_ab_series!(work.a, work.b, work.pq, coefficients, x0, order)
    coeffs = Vector{C}(undef, order + 1)
    ordinary_point_coeffs!(coeffs, work.a, work.b, X, dX, order)
    radius = R(Inf)
    @inbounds for k in max(1, length(coeffs) - 6):(length(coeffs) - 1)
        (iszero(coeffs[k + 1]) || iszero(coeffs[k])) && continue
        radius = min(radius, abs(coeffs[k]) / abs(coeffs[k + 1]))
    end
    distance = abs(target - x0)
    direction = target >= x0 ? one(R) : -one(R)
    h = min(distance, radius / 2)
    xh = h == distance ? target : x0 + direction * h
    for _ in 1:200
        value, _ = poly_pair(coeffs, xh - x0)
        _tail_small(coeffs, abs(xh - x0), value, tolerance) && break
        h *= R(2 // 3)
        xh = x0 + direction * h
    end
    value, derivative = poly_pair(coeffs, xh - x0)
    # The cancellation is taken relative to the larger of the two end values, so that a step ending
    # near a zero of the solution is not reported as a loss of precision.
    cancellation = _series_cancellation(coeffs, abs(xh - x0), max(abs(value), abs(coeffs[1])))
    return SeriesPatch{R,C}(x0, xh, coeffs), xh, value, derivative, cancellation
end

function _continue(work::_Workspace, coefficients::CoefficientSet, x0, X, dX, target, order::Int, tolerance)
    R = typeof(x0)
    C = eltype(work.pq)
    patches = SeriesPatch{R,C}[]
    cancellation = 1.0
    x = x0
    direction = target >= x0 ? 1 : -1
    apparent = apparent_points(coefficients.params)
    for _ in 1:MAX_PATCHES
        # An apparent point next to the path is crossed with its local series (ApparentSingularities).
        stop = target
        point = isempty(apparent) ? nothing : next_apparent(apparent, x, target, direction)
        if point !== nothing
            reach = R(APPARENT_STEP_FRACTION * point.radius)
            local_point = apparent_point_in(R, point, coefficients.params)
            xr = real(local_point.x1)
            if direction * (xr - x) <= reach * (1 + 1e-6)
                filter!(q -> q !== point, apparent)
                x_b = direction * (target - xr) <= reach ? target : xr + direction * reach
                basis = apparent_basis(local_point) do z
                    ordinary_ab_series!(work.a, work.b, work.pq, coefficients, z, 2)
                    return work.a[1], work.b[1]
                end
                jump = basis === nothing ? nothing : apparent_jump(basis, local_point, x, X, dX, x_b)
                if jump !== nothing
                    push!(patches, SeriesPatch{R,C}(jump.center, x, jump.coeffs), SeriesPatch{R,C}(jump.center, x_b, jump.coeffs))
                    x, X, dX = x_b, jump.X, jump.dXdx
                    x == target && return patches, X, dX, cancellation
                    continue
                end
            else
                stop = xr - direction * reach
            end
        end
        patch, x, X, dX, step_cancellation = _taylor_step(work, coefficients, x, X, dX, stop, order, tolerance)
        cancellation = max(cancellation, step_cancellation)
        push!(patches, patch)
        x == target && return patches, X, dX, cancellation
        abs(x - target) <= 4 * eps(R) * max(one(R), abs(target)) && return patches, X, dX, cancellation
    end
    error("arbitrary-precision ISEM: the Taylor continuation did not reach x = $target from x = $x0 in $MAX_PATCHES patches.")
end

# (X, dX/dx) from a chain of patches at x; `nothing` outside the chain.
function _patch_state(patches::Vector{<:SeriesPatch}, x)
    isempty(patches) && return nothing
    R = typeof(x)
    slack = 8 * eps(R) * max(one(R), abs(x))
    @inbounds for patch in patches
        lo, hi = minmax(patch.x0, patch.xend)
        if lo - slack <= x <= hi + slack
            value, derivative = poly_pair(patch.coeffs, x - patch.x0)
            return value, derivative
        end
    end
    return nothing
end

function _chain_range(patches::Vector{<:SeriesPatch})
    lo = minimum(min(patch.x0, patch.xend) for patch in patches)
    hi = maximum(max(patch.x0, patch.xend) for patch in patches)
    return lo, hi
end

# State of the infinity branch `which` (:up or :down) at x: patches inside their chain, the
# asymptotic series beyond the seed point.
function _infinity_branch_state(route::ArbitraryPrecisionRadialRoute, which::Symbol, x)
    expansion = which == :up ? route.infinity_up : route.infinity_down
    patches = which == :up ? route.infinity_up_patches : route.infinity_down_patches
    expansion === nothing && error("the $(route.branch) route has no $which expansion.")
    if expansion isa MSTInfinityBranch
        if isempty(patches) || x > _chain_range(patches)[2] + 8 * eps(typeof(x))
            x < one(x) - _real_type(route.params)(1 // 2^20) || return _infinity_state(expansion, x)
            _chain_worthwhile(route, x) || return _infinity_state(expansion, x)
            _extend_mst_chain!(route, which, x)
        end
        state = _patch_state(patches, x)
        return state === nothing ? _infinity_state(expansion, x) : state
    end
    x >= one(x) - expansion.y_seed && return _infinity_state(expansion, x)
    state = _patch_state(patches, x)
    state === nothing || return state
    error("arbitrary-precision ISEM: x = $x is outside the $which chain of the $(route.branch) route.")
end

function _horizon_branch_state(route::ArbitraryPrecisionRadialRoute, which::Symbol, x)
    expansion = which == :in ? route.horizon_in : route.horizon_out
    patches = which == :in ? route.horizon_in_patches : route.horizon_out_patches
    expansion === nothing && error("the $(route.branch) route has no horizon $which expansion.")
    x <= expansion.x_seed && return _horizon_state(expansion, x)
    state = _patch_state(patches, x)
    state === nothing || return state
    error("arbitrary-precision ISEM: x = $x is outside the horizon $which chain of the $(route.branch) route.")
end

# Evaluations of a BigFloat route run at its working precision: the stored series carry that
# precision and their sums cancel (by up to the Taylor-step bound), as does the GSN -> Teukolsky
# conversion at large frequency, so an evaluation at the caller's precision would lose those bits.
# The results are returned at the working precision, like the amplitudes.
_at_working_precision(f, route::ArbitraryPrecisionRadialRoute{<:GSNParameters,BigFloat}) =
    setprecision(f, BigFloat, route.working_bits)
_at_working_precision(f, route::ArbitraryPrecisionRadialRoute) = f()

"""
    route_state(route, x)

(X, dX/dx) of the route's solution at the compact coordinate x in (0, 1), at the working precision
of the route.
"""
function route_state(route::ArbitraryPrecisionRadialRoute, x)
    return _at_working_precision(route) do
        _route_state(route, x)
    end
end
function _route_state(route::ArbitraryPrecisionRadialRoute, x)
    R = _real_type(route.params)
    xx = R(x)
    0 < xx < 1 || throw(DomainError(x, "the compact coordinate must lie in (0, 1)."))
    if route.branch == :IN
        xx <= route.match_x && return _horizon_branch_state(route, :in, xx)
        down = _infinity_branch_state(route, :down, xx)
        up = _infinity_branch_state(route, :up, xx)
        return route.incidence * down[1] + route.reflection * up[1],
            route.incidence * down[2] + route.reflection * up[2]
    end
    xx >= route.match_x && return _infinity_branch_state(route, :up, xx)
    out = _horizon_branch_state(route, :out, xx)
    inn = _horizon_branch_state(route, :in, xx)
    return route.incidence * out[1] + route.reflection * inn[1],
        route.incidence * out[2] + route.reflection * inn[2]
end

"""
    route_state_r(route, r)

(X, dX/dr*) at the Boyer-Lindquist radius r.
"""
function route_state_r(route::ArbitraryPrecisionRadialRoute, r)
    return _at_working_precision(route) do
        p = route.params
        rr = _real_type(p)(r)
        x = x_from_r(p, rr)
        X, dXdx = _route_state(route, x)
        X, dXdx * dx_drstar(p, x, rr)
    end
end

# ---------------------------------------------------------------------------------------------
# Public wrappers
# ---------------------------------------------------------------------------------------------
function _route_mode(route::ArbitraryPrecisionRadialRoute)
    p = route.params
    return getfield(_root_module(), :Mode)(p.s, p.l, p.m, p.a, p.omega, p.lambda)
end

function _route_boundary_condition(route::ArbitraryPrecisionRadialRoute)
    root = _root_module()
    return route.branch == :IN ? getfield(root, :IN) : getfield(root, :UP)
end

# r* -> (X, dX/dr*, 0) in the floating type of the route.
struct GSNStateEvaluator{Rt}
    route::Rt
end
# The compact coordinate is formed from the distance to the horizon, x = h / (h + r_+ - r_-), so that
# deep in the horizon region (r* -> -infinity) x keeps its relative precision; through r = r_+ + h it
# would keep only the digits of r beyond those of r_+.
function (e::GSNStateEvaluator)(rs)
    return _at_working_precision(e.route) do
        p = e.route.params
        R = _real_type(p)
        h = horizon_distance_from_rstar(p.a, R(rs))
        x = h / (h + r_plus(p) - r_minus(p))
        X, dXdx = _route_state(e.route, x)
        (X, dXdx * dx_drstar(p, x, r_plus(p) + h), zero(R))
    end
end

function gsn_radial_function(route::ArbitraryPrecisionRadialRoute)
    p = route.params
    root = _root_module()
    C = _complex_type(p)
    rsmp = rstar_from_r(p.a, r_from_x(p, route.match_x))
    metadata = merge(route.diagnostics, (
        xm = Float64(route.match_x),
        tol = Float64(route.tolerance),
        ordinary_order = route.order,
        unit_normalization_available = true,
        normalization = :unit_gsn_transmission,
    ))
    return getfield(root, :GSNRadialFunction)(
        _route_mode(route),
        _route_boundary_condition(route),
        missing,
        missing,
        rsmp,
        2 * route.order,
        3 * route.order,
        one(C),
        route.incidence,
        route.reflection,
        metadata,
        missing,
        GSNStateEvaluator(route),
        getfield(root, :UNIT_GSN_TRANS),
        "ISEM",
    )
end

function _teukolsky_amplitude_factors(route::ArbitraryPrecisionRadialRoute)
    p = route.params
    R = _real_type(p)
    if R !== BigFloat
        # The closed-form factors cancel by several digits at some frequencies (Cref at a = 0.9 near the
        # positive imaginary axis); a fixed-precision type has no guard, so they are evaluated in
        # BigFloat at the type's precision plus the guard and converted.
        C = _complex_type(p)
        factors = setprecision(BigFloat, precision(R) + GUARD_BITS) do
            _teukolsky_amplitude_factors_of(_bigfloat_parameters(p, precision(R) + GUARD_BITS), route.branch)
        end
        return (trans=C(factors.trans), incidence=C(factors.incidence), reflection=C(factors.reflection))
    end
    return _at_working_precision(route) do
        _teukolsky_amplitude_factors_of(p, route.branch)
    end
end
function _teukolsky_amplitude_factors_of(p::GSNParameters, branch::Symbol)
    begin
        if branch == :IN
            (
                trans = ConversionFactors.Btrans(p.s, p.m, p.a, p.omega, p.lambda),
                incidence = ConversionFactors.Binc(p.s, p.m, p.a, p.omega, p.lambda),
                reflection = ConversionFactors.Bref(p.s, p.m, p.a, p.omega, p.lambda),
            )
        else
            (
                trans = ConversionFactors.Ctrans(p.s, p.m, p.a, p.omega, p.lambda),
                incidence = ConversionFactors.Cinc(p.s, p.m, p.a, p.omega, p.lambda),
                reflection = ConversionFactors.Cref(p.s, p.m, p.a, p.omega, p.lambda),
            )
        end
    end
end

# r -> (R, dR/dr) divided by `scale`, through the closed-form GSN-to-Teukolsky converter.
struct TeukolskyPairEvaluator{Rt,Cv,T}
    route::Rt
    converter::Cv
    inv_scale::T
end
function _teukolsky_pair_with_error(e::TeukolskyPairEvaluator, r)
    return _at_working_precision(e.route) do
        p = e.route.params
        rr = _real_type(p)(r)
        X, dXdx = _route_state(e.route, x_from_r(p, rr))
        R, Rp = e.converter(rr, X, dXdx)
        e.inv_scale * R, e.inv_scale * Rp, zero(_real_type(p))
    end
end
function (e::TeukolskyPairEvaluator)(r)
    R, Rp, _ = _teukolsky_pair_with_error(e, r)
    return R, Rp
end

function teukolsky_radial_function(route::ArbitraryPrecisionRadialRoute; normalize=:teukolsky)
    p = route.params
    root = _root_module()
    factors = _teukolsky_amplitude_factors(route)
    normalize in (:teukolsky, :gsn) || throw(ArgumentError("normalize must be :teukolsky or :gsn."))
    convention = getfield(root, normalize == :teukolsky ? :UNIT_TEUKOLSKY_TRANS : :UNIT_GSN_TRANS)
    # The converter precomputes its constants, whose combinations cancel at large frequency, so it is
    # built at the working precision like everything the evaluators use.
    pair, trans, incidence, reflection = _at_working_precision(route) do
        scale = normalize == :teukolsky ? factors.trans : one(factors.trans)
        TeukolskyPairEvaluator(route, _teukolsky_converter(p), inv(scale)),
            factors.trans / scale, route.incidence * factors.incidence / scale,
            route.reflection * factors.reflection / scale
    end
    solution = TeukolskySolution(
        TeukolskyTupleEvaluator(p, pair, d2R), pair, TeukolskyRadialEvaluator(pair))
    return getfield(root, :TeukolskyRadialFunction)(
        _route_mode(route),
        _route_boundary_condition(route),
        trans,
        incidence,
        reflection,
        missing,
        gsn_radial_function(route),
        solution,
        convention,
    )
end

function y_radial_function(route::ArbitraryPrecisionRadialRoute)
    p = route.params
    (p.s == -2 && route.branch == :IN) || (p.s == 2 && route.branch == :UP) ||
        error("ISEM Y transformation supports only s = -2 with IN and s = +2 with UP.")
    root = _root_module()
    isem = _isem_module()
    factors = _teukolsky_amplitude_factors(route)
    amplitudes = (
        trans = factors.trans,
        incidence = factors.incidence * route.incidence,
        reflection = factors.reflection * route.reflection,
    )
    gsn_solution = GSNStateEvaluator(route)
    conversion_matrix = r -> Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(
        p.s, p.m, p.a, p.omega, p.lambda, r)
    ysoln_inner = GSN_to_Y_solution_from_matrix(
        r -> rstar_from_r(p.a, r), conversion_matrix, gsn_solution,
        p.s, p.m, p.a, p.omega, p.lambda, amplitudes.trans)
    ysoln = r -> _at_working_precision(() -> ysoln_inner(r), route)
    xsoln = r -> _at_working_precision(() -> gsn_solution(rstar_from_r(p.a, r))[1], route)
    yscalar = r -> ysoln(r)[1]
    return getfield(isem, :YRadialFunction)(
        _route_mode(route),
        _route_boundary_condition(route),
        amplitudes.trans,
        amplitudes.incidence,
        amplitudes.reflection,
        missing,
        missing,
        xsoln,
        yscalar,
        ysoln,
        getfield(root, :UNIT_TEUKOLSKY_TRANS),
    )
end

end
