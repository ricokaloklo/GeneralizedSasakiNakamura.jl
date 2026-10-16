module GeneralizedSasakiNakamura

include("Homogeneous/Kerr.jl")
include("Homogeneous/Coordinates.jl")
include("Homogeneous/AsymptoticExpansionCoefficients.jl")
include("Homogeneous/InitialConditions.jl")
include("Homogeneous/ConversionFactors.jl")
include("Homogeneous/Transformation.jl")
include("Homogeneous/Potentials.jl")
include("Homogeneous/Solutions.jl")
include("Homogeneous/ExtremalRadial.jl")
include("Homogeneous/ComplexFrequencies.jl")

using .Coordinates
export r_from_rstar, rstar_from_r # Useful to be exposed
using .Solutions

using SpinWeightedSpheroidalHarmonics
using DifferentialEquations # Should have been compiled by now
using Logging, LoggingExtras
using KerrGeodesics

export GSN_radial, Teukolsky_radial # Homogeneous solutions
export GSN_pointparticle_mode, Teukolsky_pointparticle_mode, Teukolsky_pointparticle_flux # Inhomogeneous solutions
export BoundaryCondition, NormalizationConvention, Mode, GSNRadialFunction, TeukolskyRadialFunction, YRadialFunction
export Y_radial
export ISEM

# Default values
_DEFAULT_rsin = -50
_DEFAULT_rsout = 1000
_DEFAULT_horizon_expansion_order = 3
_DEFAULT_infinity_expansion_order = 6

_DEFAULT_rhoin = -5000
_DEFAULT_rhoout = 5000
_DEFAULT_horizon_expansion_order_for_cplx_freq = 25
_DEFAULT_infinity_expansion_order_for_cplx_freq = 25

# Horizon series order grows like 10 log10(1/(1 - |a|)) for near-extremal spins (50 at |a| >= 1).
# Returns (scale, order); the scale also pushes rsin inwards in the legacy ODE route.
function _adaptive_horizon_expansion_order(a, order::Int)
    spin = abs(float(real(a)))
    spin >= 1 && return (1.0, max(order, 50))
    scale = log10(inv(max(1 - spin, eps(Float64))))
    return scale, max(order, _DEFAULT_horizon_expansion_order, ceil(Int, 10 * max(scale, 0) - 1e-12))
end

# Outer boundary for the numerical integration: placed where omega * rstar reaches
# _INFINITY_BOUNDARY_PHASE, i.e. rsout = max(rsout_default, phase / |omega|), so the
# large-r asymptotic series is always evaluated well inside its region of validity
# (with omega * rsout below the series order the asymptotic series diverges).
# The returned order is a CAP: fansatz/dfansatz_dr truncate the series at its
# smallest term (optimal truncation of an asymptotic series), so a generous cap is
# harmless and the actual number of terms is chosen by the series itself.
const _INFINITY_BOUNDARY_PHASE = 30.0
const _MAX_infinity_expansion_order = 40
function _adaptive_infinity_expansion_order(omega, order::Int)
    freq = abs(omega)
    freq == 0 && return order
    # GSN_radial uses rsout = max(10 * scale * rsout, rsout) with the default rsout
    scale = _INFINITY_BOUNDARY_PHASE / (freq * 10 * _DEFAULT_rsout)
    return max(scale, 0.0), max(order, _DEFAULT_infinity_expansion_order, _MAX_infinity_expansion_order)
end

# IN for purely-ingoing at the horizon and UP for purely-outgoing at infinity
# OUT for purely-outgoing at the horizon and DOWN for purely-ingoing at infinity
@enum BoundaryCondition begin
    IN = 1
    UP = 2
    OUT = 3
    DOWN = 4
end
export IN, UP, OUT, DOWN # Use these to specify the BC

# Normalization convention, UNIT_GSN_TRANS means that the transmission amplitude for GSN functions is normalized to 1, vice versa
@enum NormalizationConvention begin
    UNIT_GSN_TRANS = 1
    UNIT_TEUKOLSKY_TRANS = 2
end
export UNIT_GSN_TRANS, UNIT_TEUKOLSKY_TRANS

struct Mode
    s::Int # spin weight
    l::Int # harmonic index
    m::Int # azimuthal index
    a # Kerr spin parameter
    omega::Union{Real, Complex} # frequency
    lambda # SWSH eigenvalue
end

# Implement pretty-printing for Mode type
# REPL
function Base.show(io::IO, ::MIME"text/plain", mode::Mode)
    print(io, "Mode(s = $(mode.s), l = $(mode.l), m = $(mode.m), a = $(mode.a), omega = $(mode.omega), lambda = $(mode.lambda))")
end

struct GSNRadialFunction{F}
    mode::Mode # Information about the mode
    boundary_condition::BoundaryCondition # The boundary condition that this radial function statisfies
    rsin # The numerical inner boundary where the GSN equation is numerically evolved
    rsout # The numerical outer boundary where the GSN equation is numerically evolved
    rsmp # The matching point in tortoise coordinate if used
    horizon_expansion_order::Union{Int, Missing} # The order of the asymptotic expansion at the horizon
    infinity_expansion_order::Union{Int, Missing} # The order of the asymptotic expansion at infinity
    transmission_amplitude # In GSN formalism
    incidence_amplitude # In GSN formalism
    reflection_amplitude # In GSN formalism
    numerical_GSN_solution # Store the numerical solution to the GSN equation in [rsin, rsout]
    numerical_Riccati_solution # Store the numerical solution to the GSN equation in the Riccati form if applicable
    GSN_solution::F # Store the *full* GSN solution where asymptotic solutions are smoothly attached
    normalization_convention::NormalizationConvention # The normalization convention used for the *stored* GSN solution
    method # The method used to solve the GSN equation
end

# Implement pretty-printing for GSNRadialFunction
# Mostly to suppress the printing of the numerical solution
function Base.show(io::IO, ::MIME"text/plain", gsn_func::GSNRadialFunction)
    println(io, "GSNRadialFunction(")
    print(io, "    mode = "); show(io, "text/plain", gsn_func.mode); println(io, ",")
    println(io, "    boundary_condition = $(gsn_func.boundary_condition),")
    println(io, "    rsin = $(gsn_func.rsin),")
    println(io, "    rsout = $(gsn_func.rsout),")
    println(io, "    rsmp = $(gsn_func.rsmp),")
    println(io, "    horizon_expansion_order = $(gsn_func.horizon_expansion_order),")
    println(io, "    infinity_expansion_order = $(gsn_func.infinity_expansion_order),")
    println(io, "    transmission_amplitude = $(gsn_func.transmission_amplitude),")
    println(io, "    incidence_amplitude = $(gsn_func.incidence_amplitude),")
    println(io, "    reflection_amplitude = $(gsn_func.reflection_amplitude),")
    println(io, "    normalization_convention = $(gsn_func.normalization_convention),")
    println(io, "    method = $(gsn_func.method)")
    print(io, ")")
end

function Base.show(io::IO, gsn_func::GSNRadialFunction)
    print(io, "GSNRadialFunction(mode = Mode(s = $(gsn_func.mode.s), l = $(gsn_func.mode.l), m = $(gsn_func.mode.m), a = $(gsn_func.mode.a), omega = $(gsn_func.mode.omega), lambda = $(gsn_func.mode.lambda)), boundary_condition = $(gsn_func.boundary_condition))")
end

function Base.show(io::IO, ::MIME"text/plain", gsn_funcs::Tuple{GSNRadialFunction, GSNRadialFunction})
    println(io, "(")
    show(io, "text/plain", gsn_funcs[1]); println(io, ",")
    show(io, "text/plain", gsn_funcs[2])
    print(io, ")")
end

struct TeukolskyRadialFunction{P,G,T}
    mode::Mode # Information about the mode
    boundary_condition::BoundaryCondition # The boundary condition that this radial function statisfies
    transmission_amplitude # In Teukolsky formalism
    incidence_amplitude # In Teukolsky formalism
    reflection_amplitude # In Teukolsky formalism
    P_solution::P # Store the full P solution if available
    GSN_solution::G # Store the full GSN solution if available
    Teukolsky_solution::T # Store the full Teukolsky solution
    normalization_convention::NormalizationConvention # The normalization convention used for the *stored* Teukolsky solution
end

# Implement pretty-printing for TeukolskyRadialFunction
# Mostly to suppress the printing of the numerical solution
function Base.show(io::IO, ::MIME"text/plain", teuk_func::TeukolskyRadialFunction)
    println(io, "TeukolskyRadialFunction(")
    print(io, "    mode = "); show(io, "text/plain", teuk_func.mode); println(io, ",")
    println(io, "    boundary_condition = $(teuk_func.boundary_condition),")
    println(io, "    transmission_amplitude = $(teuk_func.transmission_amplitude),")
    println(io, "    incidence_amplitude = $(teuk_func.incidence_amplitude),")
    println(io, "    reflection_amplitude = $(teuk_func.reflection_amplitude),")
    println(io, "    normalization_convention = $(teuk_func.normalization_convention)")
    print(io, ")")
end

function Base.show(io::IO, teuk_func::TeukolskyRadialFunction)
    print(io, "TeukolskyRadialFunction(mode = Mode(s = $(teuk_func.mode.s), l = $(teuk_func.mode.l), m = $(teuk_func.mode.m), a = $(teuk_func.mode.a), omega = $(teuk_func.mode.omega), lambda = $(teuk_func.mode.lambda)), boundary_condition = $(teuk_func.boundary_condition))")
end

function Base.show(io::IO, ::MIME"text/plain", teuk_funcs::Tuple{TeukolskyRadialFunction, TeukolskyRadialFunction})
    println(io, "(")
    show(io, "text/plain", teuk_funcs[1]); println(io, ",")
    show(io, "text/plain", teuk_funcs[2])
    print(io, ")")
end

include("ISEM/ISEM.jl")

const YRadialFunction = ISEM.YRadialFunction
const Y_radial = ISEM.Y_radial

_is_auto_method(method) = method == "auto"
_is_isem_method(method) = method == "ISEM"
_use_isem_method(method) = method == "auto" || _is_isem_method(method)
_is_static_frequency(omega) = iszero(omega)
_is_horizon_superradiance_frequency(a, m, omega) = isreal(omega) &&
    !iszero(omega) && omega == m * Kerr.omega_horizon(a)

# Numerical boundaries for the exact-extremal solver. For real omega the inner start
# follows the scale of the irregular horizon expansion, h = r - 1 = 0.05 |2 omega - m a|^0.7
# (at most 0.25, at r* >= -400): starting deeper amplifies the error of the asymptotic
# horizon state for s > 0, starting shallower leaves its region of use.
function _extremal_default_bounds(m, a, omega)
    if isreal(omega)
        h = min(0.05 * abs(2 * real(omega) - m * a)^0.7, 0.25)
        rsin = h > 0 ? max(Coordinates.rstar_from_r(a, 1 + h), -400.0) : -400.0
        return (rsin, 160.0)
    end
    damping = abs(imag(omega))
    damping == 0 && return (-80.0, 160.0)
    extent = clamp(12 / float(damping), 24.0, 80.0)
    return (-extent, extent)
end

function _warn_solve_failure(method, retcode)
    alternatives = method == "Riccati" ?
        "method=\"linear\" or method=\"ISEM\"" :
        "method=\"Riccati\" or method=\"ISEM\""
    @warn "The GSN solve using method=\"$(method)\" failed with retcode=$(retcode). Consider using $(alternatives) instead."
end

function _warn_complex_solve_failures(method, solutions)
    for solution in solutions
        if !SciMLBase.successful_retcode(solution)
            _warn_solve_failure(method, solution.retcode)
        end
    end
end

function _swsh_eigenvalue(s::Int, l::Int, m::Int, c)
    if c isa Real && s in (-2, -1, 0, 1, 2) && l >= abs(s) && abs(m) <= l
        return ISEM.swsh_eigenvalue(s, l, m, c)
    end
    return spin_weighted_spheroidal_eigenvalue(s, l, m, c)
end

# ISEM takes the matching point xm in the compact coordinate, strictly inside (0, 1).
function _check_radial_options(method, xm)
    if (_is_auto_method(method) || _is_isem_method(method)) && xm !== nothing
        xm isa Real && 0 < xm < 1 || throw(ArgumentError(
            "ISEM uses xm in (0,1)."))
    end
    return nothing
end

const _LEGACY_FALLBACK_OMEGA_WARN = 1e-8

# Called when method="auto" falls back from ISEM to the ODE routes: below |omega| = 1e-8
# the inward integration loses the UP solution, so the result should not be trusted.
function _warn_low_frequency_fallback(omega)
    abs(omega) >= _LEGACY_FALLBACK_OMEGA_WARN && return nothing
    @warn "ISEM failed at |omega| = $(abs(omega)) < $(_LEGACY_FALLBACK_OMEGA_WARN); falling back to the legacy ODE routes, which are unreliable at such low frequency (the inward integration loses the UP solution). Treat this mode's result with suspicion." maxlog=10
    return nothing
end

function _combine_isem_down(Rin, Rup)
    Binc = Rin.incidence_amplitude
    Bref = Rin.reflection_amplitude
    Cinc = Rup.incidence_amplitude
    Cref = Rup.reflection_amplitude

    function teuk_func(coord)
        Rin_vals = Rin.Teukolsky_solution(coord)
        Rup_vals = Rup.Teukolsky_solution(coord)
        R = Binc^-1 * (Rin_vals[1] - Bref * Rup_vals[1])
        dR = Binc^-1 * (Rin_vals[2] - Bref * Rup_vals[2])
        d2R = Binc^-1 * (Rin_vals[3] - Bref * Rup_vals[3])
        err = max(Rin_vals[4], Rup_vals[4])
        return (R, dR, d2R, err)
    end

    return TeukolskyRadialFunction(
        Rin.mode,
        DOWN,
        one(Binc),
        Binc^-1 - Bref * Cref / Binc,
        -Bref * Cinc / Binc,
        missing,
        missing,
        teuk_func,
        UNIT_TEUKOLSKY_TRANS,
    )
end

function _combine_isem_out(Rin, Rup)
    Binc = Rin.incidence_amplitude
    Bref = Rin.reflection_amplitude
    Cinc = Rup.incidence_amplitude
    Cref = Rup.reflection_amplitude

    function teuk_func(coord)
        Rin_vals = Rin.Teukolsky_solution(coord)
        Rup_vals = Rup.Teukolsky_solution(coord)
        R = Cinc^-1 * (Rup_vals[1] - Cref * Rin_vals[1])
        dR = Cinc^-1 * (Rup_vals[2] - Cref * Rin_vals[2])
        d2R = Cinc^-1 * (Rup_vals[3] - Cref * Rin_vals[3])
        err = max(Rin_vals[4], Rup_vals[4])
        return (R, dR, d2R, err)
    end

    return TeukolskyRadialFunction(
        Rin.mode,
        OUT,
        one(Cinc),
        Cinc^-1 - Bref * Cref / Cinc,
        -Binc * Cref / Cinc,
        missing,
        missing,
        teuk_func,
        UNIT_TEUKOLSKY_TRANS,
    )
end

function _isem_branch(boundary_condition)
    boundary_condition == IN && return :IN
    boundary_condition == UP && return :UP
    error("method = \"ISEM\" constructs DOWN and OUT from IN and UP; this boundary condition is not supported here.")
end

function _isem_route(s, l, m, a, omega, boundary_condition; xm=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing)
    return ISEM.radial_route(
        s,
        l,
        m,
        a,
        omega,
        _isem_branch(boundary_condition);
        N=N,
        xm=xm,
        tol=tol,
        small_frequency_expansion=sfe === nothing ? :auto : sfe,
        large_frequency_expansion=lfe === nothing ? :auto : lfe,
    )
end

function _combine_gsn_down(Xin, Xup)
    Btrans = Xin.transmission_amplitude
    Binc = Xin.incidence_amplitude
    Bref = Xin.reflection_amplitude
    Ctrans = Xup.transmission_amplitude
    Cinc = Xup.incidence_amplitude
    Cref = Xup.reflection_amplitude

    function solution(rs)
        xin = Xin.GSN_solution(rs)
        xup = Xup.GSN_solution(rs)
        scale = inv(Binc)
        return (
            scale * (xin[1] - Bref / Ctrans * xup[1]),
            scale * (xin[2] - Bref / Ctrans * xup[2]),
            max(xin[3], xup[3]),
        )
    end

    return GSNRadialFunction(
        Xin.mode,
        DOWN,
        missing,
        missing,
        missing,
        max(Xin.horizon_expansion_order, Xup.horizon_expansion_order),
        max(Xin.infinity_expansion_order, Xup.infinity_expansion_order),
        one(Binc),
        Btrans / Binc - Bref * Cref / (Binc * Ctrans),
        -Bref * Cinc / (Binc * Ctrans),
        (fundamental_solutions = (Xin, Xup),),
        missing,
        solution,
        UNIT_GSN_TRANS,
        "ISEM",
    )
end

function _combine_gsn_out(Xin, Xup)
    Btrans = Xin.transmission_amplitude
    Binc = Xin.incidence_amplitude
    Bref = Xin.reflection_amplitude
    Ctrans = Xup.transmission_amplitude
    Cinc = Xup.incidence_amplitude
    Cref = Xup.reflection_amplitude

    function solution(rs)
        xin = Xin.GSN_solution(rs)
        xup = Xup.GSN_solution(rs)
        scale = inv(Cinc)
        return (
            scale * (xup[1] - Cref / Btrans * xin[1]),
            scale * (xup[2] - Cref / Btrans * xin[2]),
            max(xin[3], xup[3]),
        )
    end

    return GSNRadialFunction(
        Xin.mode,
        OUT,
        missing,
        missing,
        missing,
        max(Xin.horizon_expansion_order, Xup.horizon_expansion_order),
        max(Xin.infinity_expansion_order, Xup.infinity_expansion_order),
        one(Cinc),
        Ctrans / Cinc - Bref * Cref / (Btrans * Cinc),
        -Binc * Cref / (Btrans * Cinc),
        (fundamental_solutions = (Xin, Xup),),
        missing,
        solution,
        UNIT_GSN_TRANS,
        "ISEM",
    )
end

# IN solutions at |a| >= 0.99 with Re(omega) <= -10 and |Im(omega)| <= 1e-3 |Re(omega)|
# (and no explicit matching controls) are built from the conjugated P-equation route of
# the positive-frequency partner (s, l, -m, -conj(omega)); see `_gsn_from_isem`.
function _use_lhp_highspin_p_fallback(a, omega, boundary_condition; xm, N, tol, sfe, lfe)
    boundary_condition == IN || return false
    xm === nothing && N === nothing && tol === nothing &&
        sfe === nothing && lfe === nothing || return false
    omega isa Complex || return false
    real_part = real(omega)
    real_part <= -10 || return false
    abs(imag(omega)) <= 1.0e-3 * abs(real_part) || return false
    return abs(real(a)) >= 0.99
end

function _gsn_from_isem(s, l, m, a, omega, boundary_condition; xm=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing)
    if boundary_condition == DOWN || boundary_condition == OUT
        Xin = _gsn_from_isem(s, l, m, a, omega, IN;
            xm=xm, N=N, tol=tol, sfe=sfe, lfe=lfe)
        Xup = _gsn_from_isem(s, l, m, a, omega, UP;
            xm=xm, N=N, tol=tol, sfe=sfe, lfe=lfe)
        return boundary_condition == DOWN ?
            _combine_gsn_down(Xin, Xup) :
            _combine_gsn_out(Xin, Xup)
    end
    if _use_lhp_highspin_p_fallback(
            a, omega, boundary_condition; xm, N, tol, sfe, lfe)
        positive_omega = -conj(omega)
        route = ISEM.radial_route(
            s, l, -m, a, positive_omega,
            _isem_branch(boundary_condition);
            backend=:p_equation,
            N=40,
            tol=1.0e-14,
        )
        positive = ISEM.gsn_radial_function(route)
        return ISEM._conjugate_gsn_radial_function(
            positive, m, omega)
    end
    route = _isem_route(s, l, m, a, omega, boundary_condition; xm=xm, N=N, tol=tol, sfe=sfe, lfe=lfe)
    return ISEM.gsn_radial_function(route)
end

function _extremal_method_symbol(method)
    (_is_auto_method(method) || _is_isem_method(method)) &&
        return :gsn_isem
    method == "Riccati" && return :riccati
    method == "linear" && return :linear
    throw(ArgumentError(
        "exact-extremal Kerr supports method=\"auto\", \"ISEM\", " *
        "\"Riccati\", or \"linear\""))
end

function _extremal_method_label(method)
    symbol = _extremal_method_symbol(method)
    symbol === :gsn_isem && return "ISEM"
    symbol === :riccati && return "Riccati"
    return "linear"
end

function _combine_extremal_gsn(Xin, Xup, boundary_condition, method_label)
    Btrans = Xin.transmission_amplitude
    Binc = Xin.incidence_amplitude
    Bref = Xin.reflection_amplitude
    Ctrans = Xup.transmission_amplitude
    Cinc = Xup.incidence_amplitude
    Cref = Xup.reflection_amplitude

    if boundary_condition == DOWN
        solution = rs -> begin
            xin = Xin.GSN_solution(rs)
            xup = Xup.GSN_solution(rs)
            inv(Binc) .* (xin .- Bref / Ctrans .* xup)
        end
        return GSNRadialFunction(
            Xin.mode, DOWN, Xin.rsin, Xin.rsout, missing,
            Xin.horizon_expansion_order, Xin.infinity_expansion_order,
            one(Binc), Btrans / Binc - Bref * Cref / (Binc * Ctrans),
            -Bref * Cinc / (Binc * Ctrans),
            (fundamental_solutions=(Xin, Xup),), missing, solution,
            UNIT_GSN_TRANS, method_label)
    end

    solution = rs -> begin
        xin = Xin.GSN_solution(rs)
        xup = Xup.GSN_solution(rs)
        inv(Cinc) .* (xup .- Cref / Btrans .* xin)
    end
    return GSNRadialFunction(
        Xin.mode, OUT, Xin.rsin, Xin.rsout, missing,
        Xin.horizon_expansion_order, Xin.infinity_expansion_order,
        one(Cinc), Ctrans / Cinc - Bref * Cref / (Btrans * Cinc),
        -Binc * Cref / (Btrans * Cinc),
        (fundamental_solutions=(Xin, Xup),), missing, solution,
        UNIT_GSN_TRANS, method_label)
end

function _gsn_from_extremal(s, l, m, a, omega, boundary_condition,
        rsin, rsout; horizon_expansion_order, infinity_expansion_order,
        method, data_type, ODE_algorithm, tolerance, lambda=nothing)
    method_symbol = _extremal_method_symbol(method)
    method_label = _extremal_method_label(method)
    lambda_value = lambda === nothing ?
        _swsh_eigenvalue(s, l, m, a * omega) : lambda

    if boundary_condition == DOWN || boundary_condition == OUT
        Xin = _gsn_from_extremal(
            s, l, m, a, omega, IN, rsin, rsout;
            horizon_expansion_order, infinity_expansion_order,
            method, data_type, ODE_algorithm, tolerance,
            lambda=lambda_value)
        Xup = _gsn_from_extremal(
            s, l, m, a, omega, UP, rsin, rsout;
            horizon_expansion_order, infinity_expansion_order,
            method, data_type, ODE_algorithm, tolerance,
            lambda=lambda_value)
        return _combine_extremal_gsn(
            Xin, Xup, boundary_condition, method_label)
    end

    branch = boundary_condition == IN ? :IN :
        boundary_condition == UP ? :UP :
        throw(ArgumentError("unknown boundary condition $boundary_condition"))
    solved = ExtremalRadial.solve_extremal_gsn(
        s, m, a, omega, lambda_value, branch, rsin, rsout;
        method=method_symbol,
        horizon_order=max(horizon_expansion_order, 48),
        infinity_order=max(infinity_expansion_order,
            isreal(omega) ? 12 : 18),
        series_order=36,
        data_type, ode_algorithm=ODE_algorithm, tolerance)
    mode = Mode(s, l, m, a, omega, lambda_value)
    return GSNRadialFunction(
        mode, boundary_condition, rsin, rsout, missing,
        max(horizon_expansion_order, 48),
        max(infinity_expansion_order, isreal(omega) ? 12 : 18),
        solved.transmission, solved.incidence, solved.reflection,
        solved.numerical, solved.riccati, solved.solution,
        UNIT_GSN_TRANS, method_label)
end

function _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing)
    if boundary_condition == DOWN || boundary_condition == OUT
        Rin = _teukolsky_from_isem(s, l, m, a, omega, IN;
            xm=xm, N=N, tol=tol, sfe=sfe, lfe=lfe)
        Rup = _teukolsky_from_isem(s, l, m, a, omega, UP;
            xm=xm, N=N, tol=tol, sfe=sfe, lfe=lfe)
        return boundary_condition == DOWN ?
            _combine_isem_down(Rin, Rup) : _combine_isem_out(Rin, Rup)
    end
    route = _isem_route(s, l, m, a, omega, boundary_condition; xm=xm, N=N, tol=tol, sfe=sfe, lfe=lfe)
    return ISEM.teukolsky_radial_function(route)
end

@doc raw"""
    GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order, method="auto", data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, rsmp=nothing, lambda=nothing, xm=nothing, N=nothing, sfe=nothing, lfe=nothing)

Compute the GSN function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency [which *can be complex*]) 
and boundary condition specified by `boundary_condition`, which can be either

    - `IN` for purely-ingoing at the horizon,
    - `UP` for purely-outgoing at infinity,
    - `OUT` for purely-outgoing at the horizon,
    - `DOWN` for purely-ingoing at infinity.

Note that the `OUT` and `DOWN` solutions are constructed by linearly combining the `IN` and `UP` solutions.

The GSN function is numerically solved, for real values of `omega`, in the interval of *tortoise coordinates* $r_{*} \in$ `[rsin, rsout]` using the ODE solver (from `DifferentialEquations.jl`) specified by `ODE_algorithm` (default: `AutoVern9(Rosenbrock23(autodiff=false))`) 
with tolerance specified by `tolerance` (default: `1e-12`); `tol` is an alias that takes precedence over `tolerance`.
A precomputed spin-weighted spheroidal eigenvalue can be passed as `lambda`; by default it is computed from `s`, `l`, `m` and `a*omega`. The solution method is determined by the keyword `method` (default: `auto`).
The default `auto` tries `ISEM` (iterative series expansion matching of the GSN equation) first and falls back to the legacy automatic radial route only if ISEM fails and `xm` was not supplied. With explicit `xm`, the failure is reported instead of switching to a solver that cannot use that junction. Explicit `method="ISEM"` is strict; the ODE paths (`"Riccati"`, `"linear"`) can be selected explicitly. At the superradiant threshold `omega = m * Omega_H` the `UP` solution is always built with ISEM.
Omit `xm` for automatic matching; an explicit `xm` can introduce a cancellation error in the incidence amplitude that the matching tolerance does not bound.
By default the data type used is `ComplexF64` (i.e. double-precision floating-point number) but it can be changed by
specifying `data_type` (e.g. `Complex{BigFloat}` for complex arbitrary precision number).

For the ODE methods (`"linear"`, `"Riccati"`) and complex `omega`, the GSN function is solved along a rotated path on the complex plane of $r_*$, 
where the path consists of two broken line segments parametrized by a real variable/a new coordinate $\rho$.
The angle with which the path is rotated is determined by the frequency $\omega$, the Kerr spin parameter $a$ and the azimuthal index $m$, 
such that the GSN function still behaves like a plane wave near the horizon and spatial infinity. 
At $\rho = 0$, the path intersects with the real axis at $r_{*} = r_{*}^{\rm{mp}}$ (`rsmp`, default to be nothing, which means it will be determined automatically).
Both the numerical and semi-analytical GSN solutions are evaluated as functions of $\rho$ instead of the now-complex $r_{*}$.

While the numerical GSN solution is only accurate in the range `[rsin, rsout]`, 
the full GSN solution is constructed by smoothly attaching the asymptotic solutions near horizon (up to `horizon_expansion_order`-th order) 
and infinity (up to `infinity_expansion_order`-th order).

Note, however, when `omega = 0`, the exact GSN function expressed using Gauss hypergeometric functions will be returned (i.e., instead of being solved numerically). Small nonzero frequencies retain their finite-frequency boundary conditions.
In this case, only `s`, `l`, `m`, `a`, `omega`, `boundary_condition` will be parsed.

Return a `GSNRadialFunction` object which contains all the information about the GSN solution.
"""
function GSN_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order,
    method="auto", data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, rsmp=nothing,
    xm=nothing, sfe=nothing, lfe=nothing, N=nothing, lambda=nothing
)
    _check_radial_options(method, xm)
    requested_tolerance = tol === nothing ? tolerance : tol
    if _is_auto_method(method) && !_is_static_frequency(omega) &&
            !ExtremalRadial.is_exact_extremal_spin(a)
        try
            return _gsn_from_isem(
                s, l, m, a, omega, boundary_condition;
                xm=xm, N=N, tol=requested_tolerance, sfe=sfe, lfe=lfe)
        catch err
            # The legacy route cannot honor a requested ISEM junction.
            xm === nothing || rethrow()
            _warn_low_frequency_fallback(omega)
            @warn "ISEM failed under method=\"auto\"; falling back to the legacy automatic radial route." s l m a omega boundary_condition failure=sprint(showerror, err)
        end
    end
    # Compile the legacy solver only when it is actually requested or needed.
    return Base.inferencebarrier(_gsn_radial_legacy)(
        s, l, m, a, omega, boundary_condition, rsin, rsout;
        horizon_expansion_order, infinity_expansion_order, method, data_type,
        ODE_algorithm, tolerance, tol, rsmp, xm, sfe, lfe,
        N, lambda)
end

function _gsn_radial_legacy(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order,
    infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order,
    method="auto", data_type=Solutions._DEFAULTDATATYPE,
    ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, rsmp=nothing,
    xm=nothing, sfe=nothing, lfe=nothing,
    N=nothing, lambda=nothing
)
    requested_tolerance = tol === nothing ? tolerance : tol
    tolerance = requested_tolerance === nothing ?
        Solutions._DEFAULTTOLERANCE : requested_tolerance
    if _is_static_frequency(omega)
        return GSN_radial(s, l, m, a, zero(omega), boundary_condition)
    elseif ExtremalRadial.is_exact_extremal_spin(a)
        return _gsn_from_extremal(
            s, l, m, a, omega, boundary_condition, rsin, rsout;
            horizon_expansion_order, infinity_expansion_order,
            method, data_type, ODE_algorithm, tolerance, lambda)
    elseif _is_isem_method(method) ||
            (_is_horizon_superradiance_frequency(a, m, omega) && boundary_condition == UP)
        # The ODE routes have no outgoing solution at the superradiant threshold; ISEM does.
        return _gsn_from_isem(s, l, m, a, omega, boundary_condition;
            xm=xm, N=N, tol=requested_tolerance, sfe=sfe, lfe=lfe)
    else
        sawwarn = Ref(false)
        logger = EarlyFilteredLogger(
            (log)->(log.level == Logging.Warn && (sawwarn[] = true); log.level ≥ Logging.Error),
            current_logger())
        # Compute the SWSH eigenvalue
        lambda = lambda === nothing ? _swsh_eigenvalue(s, l, m, a * omega) :
            lambda
        # Fill in the mode information
        mode = Mode(s, l, m, a, omega, lambda)
        horizon_scale, horizon_expansion_order = _adaptive_horizon_expansion_order(a, horizon_expansion_order)
        infinity_scale, infinity_expansion_order = _adaptive_infinity_expansion_order(omega, infinity_expansion_order)
        rsin = min(horizon_scale * rsin, rsin)
        rsout = max(10 * infinity_scale * rsout, rsout)
        if boundary_condition == IN
            # Solve for Xin
            if isa(omega, Real)
                if method == "auto"
                    # For real frequencies, we first try the Riccati form
                    sol = with_logger(logger) do
                        GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="Riccati", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
                    end

                    if sawwarn[]
                        return GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="linear", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
                    else 
                        return sol
                    end
                end

                if method == "Riccati"
                    Phiinsoln = Solutions.solve_Phiin(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=horizon_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
                    if !SciMLBase.successful_retcode(Phiinsoln)
                        _warn_solve_failure(method, Phiinsoln.retcode)
                    end
                    # Then convert to Xin
                    Xinsoln = Solutions.Xsoln_from_Phisoln(Phiinsoln)
                elseif method == "linear"
                    Phiinsoln = nothing
                    Xinsoln = Solutions.solve_Xin(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=horizon_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
                    if !SciMLBase.successful_retcode(Xinsoln)
                        _warn_solve_failure(method, Xinsoln.retcode)
                    end
                else
                    error("Method must be 'Riccati', 'linear', or 'ISEM'")
                end

                # Extract the incidence and reflection amplitudes (NOTE: transmisson amplitude is *always* 1)
                Bref_SN, Binc_SN = Solutions.BrefBinc_SN_from_Xin(s, m, a, omega, lambda, Xinsoln, rsout; order=infinity_expansion_order)

                # Construct the full, 'semi-analytical' GSN solution
                semianalytical_Xinsoln(rs) = Solutions.semianalytical_Xin(s, m, a, omega, lambda, Xinsoln, rsin, rsout, horizon_expansion_order, infinity_expansion_order, rs)

                return GSNRadialFunction(
                    mode,
                    IN,
                    rsin,
                    rsout,
                    rsmp,
                    horizon_expansion_order,
                    infinity_expansion_order,
                    data_type(1),
                    Binc_SN,
                    Bref_SN,
                    Xinsoln,
                    Phiinsoln,
                    semianalytical_Xinsoln,
                    UNIT_GSN_TRANS,
                    method
                )
            else
                if method == "auto"
                    # For complex frequencies, we first try the linear form
                    sol = with_logger(logger) do
                        GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="linear", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
                    end
                    if sawwarn[]
                        return GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="Riccati", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
                    else 
                        return sol
                    end
                end

                p = omega - m*Kerr.omega_horizon(a)

                # First solve for r in terms of rho,
                # the distance along the rotated path on the complex plane
                rho_min = rsin
                rho_max = rsout

                if isnothing(rsmp)
                    r_from_rho, rsmp = ComplexFrequencies.solve_r_from_rho(
                        s, m, a, omega, lambda,
                        -angle(p), -angle(omega),
                        rho_min, rho_max; sign_pos=ComplexFrequencies.determine_sign(omega), sign_neg=ComplexFrequencies.determine_sign(p)
                    )
                    if isnan(rsmp)
                        rsmp = 0 # Optimization has failed actually
                    end
                else
                    r_from_rho = ComplexFrequencies.solve_r_from_rho(
                        a, -angle(p), -angle(omega),
                        rsmp, rho_min, rho_max; sign_pos=ComplexFrequencies.determine_sign(omega), sign_neg=ComplexFrequencies.determine_sign(p)
                    )
                end

                if method == "Riccati"
                    Phiinsoln, Phiinsoln_pos, Phiinsoln_neg = ComplexFrequencies.solve_Phiin(
                        s, m, a, -angle(omega), -angle(p),
                        omega, lambda, r_from_rho,
                        rsmp, rho_min, rho_max;
                        initialconditions_order=horizon_expansion_order, dtype=data_type,
                        odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance
                    )
                    _warn_complex_solve_failures(
                        method, (Phiinsoln_pos, Phiinsoln_neg))
                    # Then convert to Xin
                    Xinsoln = Solutions.Xsoln_from_Phisoln(Phiinsoln)
                elseif method == "linear"
                    Phiinsoln = nothing
                    Xinsoln, Xinsoln_pos, Xinsoln_neg = ComplexFrequencies.solve_Xin(
                        s, m, a, -angle(omega), -angle(p),
                        omega, lambda, r_from_rho,
                        rsmp, rho_min, rho_max;
                        initialconditions_order=horizon_expansion_order, dtype=data_type,
                        odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance
                    )
                    _warn_complex_solve_failures(
                        method, (Xinsoln_pos, Xinsoln_neg))
                else
                    error("Method must be 'Riccati', 'linear', or 'ISEM'")
                end

                # Extract the incidence and reflection amplitudes (NOTE: transmisson amplitude is *always* 1)
                Bref_SN, Binc_SN = ComplexFrequencies.BrefBinc_SN_from_Xin(
                    s, m, a, -angle(omega), omega, lambda, Xinsoln, r_from_rho, rsmp, rho_max; order=infinity_expansion_order
                )

                # Construct the full, 'semi-analytical' GSN solution *in rho*
                semianalytical_Xinsoln_rho(rho) = ComplexFrequencies.semianalytical_Xin(s, m, a, -angle(omega), -angle(p), omega, lambda, Xinsoln, r_from_rho, rsmp, rho_min, rho_max, horizon_expansion_order, infinity_expansion_order, rho)

                return GSNRadialFunction(
                    mode,
                    IN,
                    rsin,
                    rsout,
                    rsmp,
                    horizon_expansion_order,
                    infinity_expansion_order,
                    data_type(1),
                    Binc_SN,
                    Bref_SN,
                    Xinsoln,
                    Phiinsoln,
                    semianalytical_Xinsoln_rho,
                    UNIT_GSN_TRANS,
                    method
                )
            end
        elseif boundary_condition == UP
            # Solve for Xup
            if isa(omega, Real)
                if method == "auto"
                    # For real frequencies, we first try the Riccati form
                    sol = with_logger(logger) do
                        GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="Riccati", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
                    end
                    if sawwarn[]
                        return GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="linear", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
                    else 
                        return sol
                    end
                end

                if method == "Riccati"
                    # Actually solve for Phiup first
                    Phiupsoln = Solutions.solve_Phiup(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=infinity_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
                    if !SciMLBase.successful_retcode(Phiupsoln)
                        _warn_solve_failure(method, Phiupsoln.retcode)
                    end
                    # Then convert to Xup
                    Xupsoln = Solutions.Xsoln_from_Phisoln(Phiupsoln)
                elseif method == "linear"
                    Phiupsoln = nothing
                    Xupsoln = Solutions.solve_Xup(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=infinity_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
                    if !SciMLBase.successful_retcode(Xupsoln)
                        _warn_solve_failure(method, Xupsoln.retcode)
                    end
                else
                    error("Method must be 'Riccati', 'linear', or 'ISEM'")
                end

                # Extract the incidence and reflection amplitudes (NOTE: transmisson amplitude is *always* 1)
                Cref_SN, Cinc_SN = Solutions.CrefCinc_SN_from_Xup(s, m, a, omega, lambda, Xupsoln, rsin; order=horizon_expansion_order)

                # Construct the full, 'semi-analytical' GSN solution
                semianalytical_Xupsoln(rs) = Solutions.semianalytical_Xup(s, m, a, omega, lambda, Xupsoln, rsin, rsout, horizon_expansion_order, infinity_expansion_order, rs)

                return GSNRadialFunction(
                    mode,
                    UP,
                    rsin,
                    rsout,
                    rsmp,
                    horizon_expansion_order,
                    infinity_expansion_order,
                    data_type(1),
                    Cinc_SN,
                    Cref_SN,
                    Xupsoln,
                    Phiupsoln,
                    semianalytical_Xupsoln,
                    UNIT_GSN_TRANS,
                    method
                )
            else
                if method == "auto"
                    # For complex frequencies, we first try the linear form
                    sol = with_logger(logger) do
                        GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="linear", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
                    end
                    if sawwarn[]
                        return GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="Riccati", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
                    else 
                        return sol
                    end
                end

                p = omega - m*Kerr.omega_horizon(a)

                # First solve for r in terms of rho,
                # the distance along the rotated path on the complex plane
                rho_min = rsin
                rho_max = rsout

                if isnothing(rsmp)
                    r_from_rho, rsmp = ComplexFrequencies.solve_r_from_rho(
                        s, m, a, omega, lambda,
                        -angle(p), -angle(omega),
                        rho_min, rho_max; sign_pos=ComplexFrequencies.determine_sign(omega), sign_neg=ComplexFrequencies.determine_sign(p)
                    )
                    if isnan(rsmp)
                        rsmp = 0 # Optimization has failed actually
                    end
                else
                    r_from_rho = ComplexFrequencies.solve_r_from_rho(
                        a, -angle(p), -angle(omega),
                        rsmp, rho_min, rho_max; sign_pos=ComplexFrequencies.determine_sign(omega), sign_neg=ComplexFrequencies.determine_sign(p)
                    )
                end

                if method == "Riccati"
                    Phiupsoln, Phiupsoln_pos, Phiupsoln_neg = ComplexFrequencies.solve_Phiup(
                        s, m, a, -angle(omega), -angle(p),
                        omega, lambda, r_from_rho,
                        rsmp, rho_min, rho_max;
                        initialconditions_order=infinity_expansion_order, dtype=data_type,
                        odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance
                    )
                    _warn_complex_solve_failures(
                        method, (Phiupsoln_pos, Phiupsoln_neg))
                    # Then convert to Xup
                    Xupsoln = Solutions.Xsoln_from_Phisoln(Phiupsoln)
                elseif method == "linear"
                    Phiupsoln = nothing
                    Xupsoln, Xupsoln_pos, Xupsoln_neg = ComplexFrequencies.solve_Xup(
                        s, m, a, -angle(omega), -angle(p),
                        omega, lambda, r_from_rho,
                        rsmp, rho_min, rho_max;
                        initialconditions_order=infinity_expansion_order, dtype=data_type,
                        odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance
                    )
                    _warn_complex_solve_failures(
                        method, (Xupsoln_pos, Xupsoln_neg))
                else
                    error("Method must be 'Riccati', 'linear', or 'ISEM'")
                end

                # Extract the incidence and reflection amplitudes (NOTE: transmisson amplitude is *always* 1)
                Cref_SN, Cinc_SN = ComplexFrequencies.CrefCinc_SN_from_Xup(
                    s, m, a, -angle(p), omega, lambda, Xupsoln, r_from_rho, rsmp, rho_min; order=horizon_expansion_order
                )

                # Construct the full, 'semi-analytical' GSN solution *in rho*
                semianalytical_Xupsoln_rho(rho) = ComplexFrequencies.semianalytical_Xup(s, m, a, -angle(omega), -angle(p), omega, lambda, Xupsoln, r_from_rho, rsmp, rho_min, rho_max, horizon_expansion_order, infinity_expansion_order, rho)

                return GSNRadialFunction(
                    mode,
                    UP,
                    rsin,
                    rsout,
                    rsmp,
                    horizon_expansion_order,
                    infinity_expansion_order,
                    data_type(1),
                    Cinc_SN,
                    Cref_SN,
                    Xupsoln,
                    Phiupsoln,
                    semianalytical_Xupsoln_rho,
                    UNIT_GSN_TRANS,
                    method
                )
            end
        elseif boundary_condition == DOWN
            # Construct Xdown from Xin and Xup, instead of solving the ODE numerically with the boundary condition

            # Solve for Xin *and* Xup first
            Xin = GSN_radial(s, l, m, a, omega, IN, rsin, rsout, horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
            Xup = GSN_radial(s, l, m, a, omega, UP, rsin, rsout, horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)

            # Xdown is a linear combination of Xin and Xup with the following coefficients
            Btrans = Xin.transmission_amplitude # Should really be just 1
            Binc = Xin.incidence_amplitude
            Bref = Xin.reflection_amplitude
            Ctrans = Xup.transmission_amplitude # Should really be just 1
            Cinc = Xup.incidence_amplitude
            Cref = Xup.reflection_amplitude

            _full_Xdown_solution(rs) = Binc^-1 .* (Xin.GSN_solution(rs) .- Bref/Ctrans .* Xup.GSN_solution(rs))

            # These solutions are "normalized" in the sense that Xdown -> exp(-i*omega*rs) near infinity
            # NOTE The definition of the "incidence" and "reflection" amplitudes follow 2101.04592, Eq. (93)
            return GSNRadialFunction(
                Xin.mode,
                DOWN,
                rsin,
                rsout,
                rsmp,
                horizon_expansion_order,
                infinity_expansion_order,
                data_type(1),
                Btrans/Binc - (Bref*Cref)/(Binc*Ctrans),
                (-Bref*Cinc)/(Binc*Ctrans),
                missing,
                missing,
                _full_Xdown_solution,
                UNIT_GSN_TRANS,
                method
            )
        elseif boundary_condition == OUT
            # Construct Xout from Xin and Xup, instead of solving the ODE numerically with the boundary condition

            # Solve for Xin *and* Xup first
            Xin = GSN_radial(s, l, m, a, omega, IN, rsin, rsout, horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)
            Xup = GSN_radial(s, l, m, a, omega, UP, rsin, rsout, horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)

            # Xout is a linear combination of Xin and Xup with the following coefficients
            Btrans = Xin.transmission_amplitude # Should really be just 1
            Binc = Xin.incidence_amplitude
            Bref = Xin.reflection_amplitude
            Ctrans = Xup.transmission_amplitude # Should really be just 1
            Cinc = Xup.incidence_amplitude
            Cref = Xup.reflection_amplitude

            _full_Xout_solution(rs) = Cinc^-1 .* (Xup.GSN_solution(rs) .- Cref/Btrans .* Xin.GSN_solution(rs))

            # These solutions are "normalized" in the sense that Xout -> exp(i*p*rs) near the horizon
            # NOTE The definition of the "incidence" and "reflection" amplitudes follow 2101.04592, Eq. (93)
            return GSNRadialFunction(
                Xin.mode,
                OUT,
                rsin,
                rsout,
                rsmp,
                horizon_expansion_order,
                infinity_expansion_order,
                data_type(1),
                Ctrans/Cinc - (Bref*Cref)/(Btrans*Cinc),
                (-Binc*Cref)/(Btrans*Cinc),
                missing,
                missing,
                _full_Xout_solution,
                UNIT_GSN_TRANS,
                method
            )
        else
            error("Boundary condition must be IN or UP")
        end
    end
end

@doc raw"""
    GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition; method="auto", tol=nothing, xm=nothing, N=nothing, sfe=nothing, lfe=nothing)

Compute the GSN function for the mode `(s, l, m, a, omega)` with boundary condition `IN`, `UP`, `OUT` or `DOWN`,
choosing the numerical domain automatically (see the method with explicit `rsin, rsout` for the meaning of `method`
and the other keywords). Negative `a` is handled through the `(m, a) -> (-m, -a)` symmetry, and `omega = 0` returns
the exact static `IN` or `UP` solution.

Return a `GSNRadialFunction`.
"""
function GSN_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition;
    method="auto", tolerance=nothing, tol=nothing,
    xm=nothing, sfe=nothing, lfe=nothing, N=nothing
)
    _check_radial_options(method, xm)
    if a < zero(a)
        positive = GSN_radial(
            s, l, -m, -a, omega, boundary_condition;
            method, tolerance, tol, xm, sfe, lfe, N)
        return ISEM._spin_reflect_gsn_radial_function(positive, m, a)
    end
    requested_tolerance = tol === nothing ? tolerance : tol
    if !_is_static_frequency(omega)
        rsin, rsout = ExtremalRadial.is_exact_extremal_spin(a) ?
            _extremal_default_bounds(m, a, omega) : (_DEFAULT_rsin, _DEFAULT_rsout)
        return GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout;
            method=method, tolerance=requested_tolerance, xm=xm, sfe=sfe, lfe=lfe, N=N)
    else
        omega0 = zero(omega)
        teuk_func = Teukolsky_radial(s, l, m, a, omega0, boundary_condition)
        GSN_solution = Solutions.Sasaki_Nakamura_function_from_Teukolsky_radial_function(s, m, a, omega0, teuk_func.mode.lambda, teuk_func.Teukolsky_solution)

        return GSNRadialFunction(
                teuk_func.mode,
                boundary_condition,
                missing,
                missing,
                missing,
                missing,
                missing,
                missing,
                missing,
                missing,
                missing,
                missing,
                GSN_solution,
                UNIT_TEUKOLSKY_TRANS,
                "static"
            )
    end
end

@doc raw"""
    GSN_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, method="auto", xm=nothing, N=nothing, sfe=nothing, lfe=nothing)

Compute the GSN function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency)
with the purely-ingoing boundary condition at the horizon (`IN`) and the purely-outgoing boundary condition at infinity (`UP`).

For the ODE methods and _real_ frequencies, the numerical inner boundary (rsin) and outer boundary (rsout) are set to `-50` and `1000`,
while the order of the asymptotic expansion at the horizon and infinity are determined automatically.
As for _complex_ frequencies, the numerical inner and the outer boundaries are determined automatically,
while the asymptotic expansions at the horizon and infinity are taken to order `25`.
These domain and order choices apply to the ODE routes; with the default `method="auto"` both solutions are first constructed with `ISEM`.

Return the tuple `(Xin, Xup)` of `GSNRadialFunction`s.
"""
function GSN_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, method="auto", xm=nothing, sfe=nothing, lfe=nothing, N=nothing)
    _check_radial_options(method, xm)
    if a < zero(a)
        positive = GSN_radial(
            s, l, -m, -a, omega;
            data_type, ODE_algorithm, tolerance, tol, method, xm,
            sfe, lfe, N)
        return map(f -> ISEM._spin_reflect_gsn_radial_function(f, m, a),
            positive)
    end
    requested_tolerance = tol === nothing ? tolerance : tol
    tolerance = requested_tolerance === nothing ?
        Solutions._DEFAULTTOLERANCE : requested_tolerance
    if _is_static_frequency(omega)
        omega0 = zero(omega)
        Xin = GSN_radial(s, l, m, a, omega0, IN)
        Xup = GSN_radial(s, l, m, a, omega0, UP)
    elseif ExtremalRadial.is_exact_extremal_spin(a)
        rsin, rsout = _extremal_default_bounds(m, a, omega)
        Xin = GSN_radial(s, l, m, a, omega, IN, rsin, rsout;
            method, data_type, ODE_algorithm, tolerance)
        Xup = GSN_radial(s, l, m, a, omega, UP, rsin, rsout;
            method, data_type, ODE_algorithm, tolerance)
    elseif _is_horizon_superradiance_frequency(a, m, omega) || _is_auto_method(method)
        Xin = GSN_radial(s, l, m, a, omega, IN;
            method=method, tolerance=requested_tolerance, xm=xm,
            sfe=sfe, lfe=lfe, N=N)
        Xup = GSN_radial(s, l, m, a, omega, UP;
            method=method, tolerance=requested_tolerance, xm=xm,
            sfe=sfe, lfe=lfe, N=N)
    elseif _is_isem_method(method)
        Xin = _gsn_from_isem(s, l, m, a, omega, IN;
            xm=xm, N=N, tol=requested_tolerance, sfe=sfe, lfe=lfe)
        Xup = _gsn_from_isem(s, l, m, a, omega, UP;
            xm=xm, N=N, tol=requested_tolerance, sfe=sfe, lfe=lfe)
    elseif isa(omega, Real)
        # Solve for Xin and Xup using the default settings
        Xin = GSN_radial(s, l, m, a, omega, IN, _DEFAULT_rsin, _DEFAULT_rsout,
            horizon_expansion_order=_DEFAULT_horizon_expansion_order,
            infinity_expansion_order=_DEFAULT_infinity_expansion_order,
            data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=nothing, method=method
        )
        Xup = GSN_radial(s, l, m, a, omega, UP, _DEFAULT_rsin, _DEFAULT_rsout,
            horizon_expansion_order=_DEFAULT_horizon_expansion_order,
            infinity_expansion_order=_DEFAULT_infinity_expansion_order,
            data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=nothing, method=method
        )

        # For real frequencies, bump up the expansion order until the solution is "sane"
        # The maximum expansion order to use
        _MAX_horizon_expansion_order = 50
        _MAX_infinity_expansion_order = 50
        # Step size when increasing the expansion order
        _STEP_horizon_expansion_order = 5
        _STEP_infinity_expansion_order = 5
        while(!Solutions.check_XinXup_sanity(Xin, Xup))
            new_horizon_expansion_order = Xin.horizon_expansion_order + _STEP_horizon_expansion_order
            new_infinity_expansion_order = Xup.infinity_expansion_order + _STEP_infinity_expansion_order

            if new_horizon_expansion_order > _MAX_horizon_expansion_order || new_infinity_expansion_order > _MAX_infinity_expansion_order
                @warn "Failed to solve the equation within a reasonable range of expansion orders. Consider widening the integration domain."
                continue
            end

            # Re-solve Xin and Xup using the updated settings
            Xin = GSN_radial(s, l, m, a, omega, IN, rsin, rsout,
                horizon_expansion_order=new_horizon_expansion_order,
                infinity_expansion_order=new_infinity_expansion_order,
                data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=nothing, method=method
            )
            Xup = GSN_radial(s, l, m, a, omega, UP, rsin, rsout,
                horizon_expansion_order=new_horizon_expansion_order,
                infinity_expansion_order=new_infinity_expansion_order,
                data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=nothing, method=method
            )
        end
    else
        # Use a different strategy for complex frequencies, widen the integration domain until the solution is "sane"
        _MIN_absrho = 5000
        _MAX_absrho = 50000
        _STEP_absrho = 5000

        rhoin = -_MIN_absrho
        rhoout = _MIN_absrho

        p = omega - m*Kerr.omega_horizon(a)
        lambda = _swsh_eigenvalue(s, l, m, a * omega)
        _, rsmp = ComplexFrequencies.solve_r_from_rho(
            s, m, a, omega, lambda,
            -angle(p), -angle(omega),
            rhoin, rhoout; sign_pos=ComplexFrequencies.determine_sign(omega), sign_neg=ComplexFrequencies.determine_sign(p)
        )

        while(isnan(rsmp))
            rhoin -= _STEP_absrho
            rhoout += _STEP_absrho

            if rhoout >= _MAX_absrho
                # Auto determination of rsmp has failed
                rsmp = 0
                continue
            end

            _, rsmp = ComplexFrequencies.solve_r_from_rho(
                s, m, a, omega, lambda,
                -angle(p), -angle(omega),
                rhoin, rhoout; sign_pos=ComplexFrequencies.determine_sign(omega), sign_neg=ComplexFrequencies.determine_sign(p)
            )
        end
        # Solve Xin and Xup using the updated settings
        Xin = GSN_radial(s, l, m, a, omega, IN, rhoin, rhoout,
            horizon_expansion_order=_DEFAULT_horizon_expansion_order_for_cplx_freq,
            infinity_expansion_order=_DEFAULT_infinity_expansion_order_for_cplx_freq,
            data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp, method=method
        )
        Xup = GSN_radial(s, l, m, a, omega, UP, rhoin, rhoout,
            horizon_expansion_order=_DEFAULT_horizon_expansion_order_for_cplx_freq,
            infinity_expansion_order=_DEFAULT_infinity_expansion_order_for_cplx_freq,
            data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp, method=method
        )
    end

    return (Xin, Xup)
end

# The power of multiple dispatch
(gsn_func::GSNRadialFunction)(rs) = gsn_func.GSN_solution(rs)[1] # Only return X(rs), discarding the first derivative


@doc raw"""
    Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order, method="auto", data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, rsmp=nothing, xm=nothing, N=nothing, sfe=nothing, lfe=nothing)

Compute the Teukolsky function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency [which *can be complex*]) 
and boundary condition specified by `boundary_condition`, which can be either

    - `IN` for purely-ingoing at the horizon,
    - `UP` for purely-outgoing at infinity,
    - `OUT` for purely-outgoing at the horizon,
    - `DOWN` for purely-ingoing at infinity.

Note that the `OUT` and `DOWN` solutions are constructed by linearly combining the `IN` and `UP` solutions.

The full GSN solution is converted to the corresponding Teukolsky solution $(R(r), dR/dr)$ and 
the incidence, reflection and transmission amplitude are converted from the GSN formalism to the Teukolsky formalism 
with the normalization convention that the transmission amplitude is normalized to 1 (i.e. `normalization_convention=UNIT_TEUKOLSKY_TRANS`).
When `method` is `"auto"`, the homogeneous Teukolsky solution tries
`ISEM` first and falls back to the legacy automatic radial route only if
ISEM fails. Explicit `method="ISEM"` is strict.

Note, however, when `omega = 0`, the exact Teukolsky function expressed using Gauss hypergeometric functions will be returned (i.e., instead of using the GSN formalism). 
In this case, only `s`, `l`, `m`, `a`, `omega`, `boundary_condition` will be parsed.

For the ODE methods (`"linear"`, `"Riccati"`) and complex `omega`, the Teukolsky function is evaluated as a function of $r$,
where the value at the corresponding location $\rho = r_{*}(r) \in \mathbb{R}$ along the rotated integration path on the complex plane is returned.
"""
function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order,
    method="auto", data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, rsmp=nothing,
    xm=nothing, sfe=nothing, lfe=nothing, N=nothing
)
    _check_radial_options(method, xm)
    requested_tolerance = tol === nothing ? tolerance : tol
    tolerance = requested_tolerance === nothing ?
        Solutions._DEFAULTTOLERANCE : requested_tolerance
    if _is_static_frequency(omega)
        return Teukolsky_radial(s, l, m, a, zero(omega), boundary_condition; method="static")
    elseif ExtremalRadial.is_exact_extremal_spin(a)
        # Continue below through the common GSN-to-Teukolsky conversion.
        nothing
    elseif _is_isem_method(method) ||
            (_is_horizon_superradiance_frequency(a, m, omega) && boundary_condition == UP)
        return _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, N=N, tol=requested_tolerance, sfe=sfe, lfe=lfe)
    end

    # Solve for the GSN solution
    gsn_func = GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout;
        horizon_expansion_order, infinity_expansion_order, method, data_type,
        ODE_algorithm, tolerance=requested_tolerance, rsmp, xm, sfe, lfe, N)

    return _teukolsky_from_gsn(gsn_func,s,m,a,omega,data_type,requested_tolerance)
end

function _teukolsky_from_gsn(gsn_func,s,m,a,omega,data_type,requested_tolerance)

    # Convert asymptotic amplitudes from GSN to Teukolsky formalism
    if gsn_func.boundary_condition == IN
        transmission_amplitude_conv_factor = ConversionFactors.Btrans(s, m, a, omega, gsn_func.mode.lambda)
        incidence_amplitude = ConversionFactors.Binc(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.incidence_amplitude / transmission_amplitude_conv_factor
        reflection_amplitude = ConversionFactors.Bref(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.reflection_amplitude / transmission_amplitude_conv_factor
    elseif gsn_func.boundary_condition == UP
        transmission_amplitude_conv_factor = ConversionFactors.Ctrans(s, m, a, omega, gsn_func.mode.lambda)
        incidence_amplitude = ConversionFactors.Cinc(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.incidence_amplitude / transmission_amplitude_conv_factor
        reflection_amplitude = ConversionFactors.Cref(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.reflection_amplitude / transmission_amplitude_conv_factor
    elseif gsn_func.boundary_condition == DOWN
        # The "transmission amplitude" transforms like Binc
        transmission_amplitude_conv_factor = ConversionFactors.Binc(s, m, a, omega, gsn_func.mode.lambda)
        # The "incidence amplitude" transforms like Btrans
        incidence_amplitude = ConversionFactors.Btrans(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.incidence_amplitude / transmission_amplitude_conv_factor
        # The "reflection amplitude" transforms like Cinc
        reflection_amplitude = ConversionFactors.Cinc(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.reflection_amplitude / transmission_amplitude_conv_factor
    elseif gsn_func.boundary_condition == OUT
        # The "tranmission amplitude" transforms like Cinc
        transmission_amplitude_conv_factor = ConversionFactors.Cinc(s, m, a, omega, gsn_func.mode.lambda)
        # The "incidence amplitude" transforms like Bref
        incidence_amplitude = ConversionFactors.Bref(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.incidence_amplitude / transmission_amplitude_conv_factor
        # The "reflection amplitude" transformslike Binc
        reflection_amplitude = ConversionFactors.Binc(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.reflection_amplitude / transmission_amplitude_conv_factor
    else
        error("Does not understand the boundary condition applied to the solution")
    end

    # Convert the GSN solution to the Teukolsky solution
    # For real-frequency ISEM UP solutions, R(r) at large r is taken from the Teukolsky asymptotic series, which
    # returns `nothing` outside its convergence radius; the GSN-to-Teukolsky
    # transformation is used there instead.
    infinity_series = gsn_func.boundary_condition == UP && isreal(omega) &&
            gsn_func.method == "ISEM" &&
            !ExtremalRadial.is_exact_extremal_spin(a) ?
        ISEM.RadialFunctions.teukolsky_infinity_series(
            gsn_func.mode,max(gsn_func.infinity_expansion_order,gsn_func.numerical_GSN_solution.ordinary_order),
            requested_tolerance === nothing ? ISEM.Iteration._ENDPOINT_SEED_TARGET :
                max(requested_tolerance,ISEM.Iteration._ENDPOINT_SEED_TARGET)) : nothing
    converted_solution = Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function(
        gsn_func.mode.s,
        gsn_func.mode.m,
        gsn_func.mode.a,
        gsn_func.mode.omega,
        gsn_func.mode.lambda,
        gsn_func.GSN_solution)
    teuk_func = r -> begin
        if infinity_series !== nothing
            endpoint = infinity_series(r)
            endpoint !== nothing && return Solutions.SA[endpoint[1],endpoint[2]]
        end
        return converted_solution(r)/transmission_amplitude_conv_factor
    end

    return TeukolskyRadialFunction(
        gsn_func.mode,
        gsn_func.boundary_condition,
        data_type(1),
        incidence_amplitude,
        reflection_amplitude,
        missing,
        gsn_func,
        teuk_func,
        UNIT_TEUKOLSKY_TRANS
    )
end

@doc raw"""
    Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition; method="auto", tol=nothing, xm=nothing, N=nothing, sfe=nothing, lfe=nothing)

Compute the Teukolsky function for the mode `(s, l, m, a, omega)` with boundary condition `IN`, `UP`, `OUT` or `DOWN`,
choosing the numerical domain automatically (see the method with explicit `rsin, rsout` for the meaning of `method`
and the other keywords). Negative `a` is handled through the `(m, a) -> (-m, -a)` symmetry, and `omega = 0` returns
the exact static `IN` or `UP` solution.

Return a `TeukolskyRadialFunction`.
"""
function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition; method="auto", xm=nothing, sfe=nothing, lfe=nothing, N=nothing, tol=nothing
)
    _check_radial_options(method, xm)
    if a < zero(a)
        positive = Teukolsky_radial(
            s, l, -m, -a, omega, boundary_condition;
            method, xm, sfe, lfe, N, tol)
        return ISEM._spin_reflect_teukolsky_radial_function(positive, m, a)
    end
    if ExtremalRadial.is_exact_extremal_spin(a) && !_is_static_frequency(omega)
        rsin, rsout = _extremal_default_bounds(m, a, omega)
        return Teukolsky_radial(
            s, l, m, a, omega, boundary_condition, rsin, rsout;
            method, tolerance=tol === nothing ? Solutions._DEFAULTTOLERANCE : tol)
    elseif !_is_static_frequency(omega) && (_is_isem_method(method) ||
            (_is_horizon_superradiance_frequency(a, m, omega) && boundary_condition == UP))
        return _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, N=N, tol=tol, sfe=sfe, lfe=lfe)
    elseif !_is_static_frequency(omega)
        return Teukolsky_radial(s, l, m, a, omega, boundary_condition,
            _DEFAULT_rsin, _DEFAULT_rsout;
            method, xm, sfe, lfe, N, tol)
    else
        # Compute the SWSH eigenvalue
        lambda = spin_weighted_spherical_eigenvalue(s, l, m)
        # Fill in the mode information
        mode = Mode(s, l, m, a, zero(omega), lambda)
        if boundary_condition == IN
            teuk_func = Solutions.solve_static_Rin(s, l, m, a)
        elseif boundary_condition == UP
            teuk_func = Solutions.solve_static_Rup(s, l, m, a)
        else
            error("Does not understand the boundary condition applied to the solution")
        end

        return TeukolskyRadialFunction(
            mode,
            boundary_condition,
            1,
            missing,
            missing,
            missing,
            missing,
            teuk_func,
            UNIT_TEUKOLSKY_TRANS
        )
    end
end

@doc raw"""
    Teukolsky_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, method="auto", xm=nothing, N=nothing, sfe=nothing, lfe=nothing)

Compute the Teukolsky function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency)
with the purely-ingoing boundary condition at the horizon (`IN`) and the purely-outgoing boundary condition at infinity (`UP`).
When `method` is `"auto"`, this overload tries `ISEM` first and falls back to the legacy automatic radial route only if ISEM fails. Explicit `method="ISEM"` is strict.

For the ODE methods and _real_ frequencies, the numerical inner boundary (rsin) and outer boundary (rsout) are set to `-50` and `1000`,
while the order of the asymptotic expansion at the horizon and infinity are determined automatically.
As for _complex_ frequencies, the numerical inner and the outer boundaries are determined automatically,
while the asymptotic expansions at the horizon and infinity are taken to order `25`.

Return the tuple `(Rin, Rup)` of `TeukolskyRadialFunction`s.
"""
function Teukolsky_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=nothing, tol=nothing, method="auto", xm=nothing, sfe=nothing, lfe=nothing, N=nothing)
    _check_radial_options(method, xm)
    if a < zero(a)
        positive = Teukolsky_radial(
            s, l, -m, -a, omega;
            data_type, ODE_algorithm, tolerance, tol, method, xm,
            sfe, lfe, N)
        return map(
            f -> ISEM._spin_reflect_teukolsky_radial_function(f, m, a),
            positive)
    end
    requested_tolerance = tol === nothing ? tolerance : tol
    tolerance = requested_tolerance === nothing ?
        Solutions._DEFAULTTOLERANCE : requested_tolerance
    if _is_static_frequency(omega) || _is_horizon_superradiance_frequency(a, m, omega) || _use_isem_method(method)
        Rin = Teukolsky_radial(s, l, m, a, omega, IN; method=method, xm=xm, sfe=sfe, lfe=lfe, N=N, tol=requested_tolerance)
        Rup = Teukolsky_radial(s, l, m, a, omega, UP; method=method, xm=xm, sfe=sfe, lfe=lfe, N=N, tol=requested_tolerance)
        return (Rin, Rup)
    end

    # NOTE This is not the most efficient implementation but ensures self-consistency
    Xin, Xup = GSN_radial(s, l, m, a, omega; data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, method=method)
    Rin = Teukolsky_radial(s, l, m, a, omega, IN, Xin.rsin, Xin.rsout;
        horizon_expansion_order=Xin.horizon_expansion_order,
        infinity_expansion_order=Xin.infinity_expansion_order,
        data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=Xin.rsmp, method=method)
    Rup = Teukolsky_radial(s, l, m, a, omega, UP, Xup.rsin, Xup.rsout;
        horizon_expansion_order=Xup.horizon_expansion_order,
        infinity_expansion_order=Xup.infinity_expansion_order,
        data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=Xup.rsmp, method=method)
    return (Rin, Rup)
end

# The power of multiple dispatch
_teukolsky_radial_value(sol, r) = sol(r)[1]
_teukolsky_radial_value(sol::ISEM.RadialFunctions.TeukolskySolution, r) =
    sol.radial_value(r)

function (teuk_func::TeukolskyRadialFunction)(r)
    sol = teuk_func.Teukolsky_solution
    if sol isa ISEM.RadialFunctions.TeukolskySolution
        return sol.radial_value(r)
    end
    return _teukolsky_radial_value(sol, r)
end


include("Inhomogeneous/AsymptoticExpansionCoefficientsY.jl")
include("Inhomogeneous/InitialConditionsY.jl")
include("Inhomogeneous/SolutionsY.jl")
include("Inhomogeneous/GridSampling.jl")
include("Inhomogeneous/ConvolutionIntegrals.jl")
include("Inhomogeneous/ModeSummation.jl")


function Base.show(io::IO, ::MIME"text/plain", ysol::SolutionsY.YSolutionResult)
    println(io, "YSolution(")
    print(io, "    basis = "); show(io, ysol.basis_type); println(io, ",")
    print(io, "    mode = "); show(io, "text/plain", ysol.mode); println(io, ",")
    println(io, "    rin = $(ysol.asymptotic.rin),")
    println(io, "    rout = $(ysol.asymptotic.rout),")
    print(io, ")")
end

struct PointParticleMode
    s::Int # spin weight
    l::Int # harmonic index
    m::Int # azimuthal index
    n::Int # radial index
    k::Int # polar index
    a # Kerr spin parameter
    omega::Union{Real, Complex, Nothing} # absent harmonics have no evaluated frequency
    lambda # SWSH eigenvalue
end

function Base.show(io::IO, ::MIME"text/plain", mode::PointParticleMode)
    print(io, "Mode(s = $(mode.s), l = $(mode.l), m = $(mode.m), n = $(mode.n), k = $(mode.k), a = $(mode.a), omega = $(mode.omega), lambda = $(mode.lambda))")
end

struct GSNPointParticleMode
    mode::PointParticleMode # Information about the mode, where the frequency will be computed from orbital parameters
    amplitude # In GSN formalism
    energy_flux # Identical in both formalisms
    angular_momentum_flux # Identical in both formalisms
    Carter_const_flux # Identical in both formalisms
    trajectory
    Y_solution # The Y solution used to compute the fluxes
    SWSH # The spin-weighted spheroidal harmonics used to compute the fluxes
    method # Method used in computing the convolution integral
end

struct TeukolskyPointParticleMode
    mode::PointParticleMode # Information about the mode, where the frequency will be computed from orbital parameters
    amplitude # In Teukolsky formalism
    energy_flux # Identical in both formalisms
    angular_momentum_flux # Identical in both formalisms
    Carter_const_flux # Identical in both formalisms
    trajectory
    Y_solution # The Y solution used to compute the fluxes
    SWSH # The spin-weighted spheroidal harmonics used to compute the fluxes
    method # Method used in computing the convolution integral
end

struct TeukolskyPointParticleFlux
    a
    p
    e
    x
    orbit_type::Symbol
    infinity_energy_flux
    infinity_angular_momentum_flux
    infinity_carter_constant_flux
    horizon_energy_flux
    horizon_angular_momentum_flux
    horizon_carter_constant_flux
    total_modes::Int
    tolerance
    truncation_floor
    reached::NamedTuple
    cost
    result
end

function _flux_last_shell_energy(result, key::Symbol)
    shell_list = _flux_get(result, key)
    if shell_list === nothing || isempty(shell_list)
        return nothing
    end
    return shell_list[end]
end

function _flux_reached_print_value(flux::TeukolskyPointParticleFlux, key::Symbol, value)
    value !== nothing && return value
    nmax = _flux_get(flux.result, :n_max)
    nmax === nothing && return value

    shell_key = if key == :n_reached_inf
        :infinity_energy_flux_list
    elseif key == :n_reached_hor
        :horizon_energy_flux_list
    else
        return value
    end

    last_shell_energy = _flux_last_shell_energy(flux.result, shell_key)
    if last_shell_energy === nothing
        return "$nmax (hit nmax = $nmax)"
    end
    return "$nmax (hit nmax = $nmax; last shell energy = $last_shell_energy)"
end

function _flux_hit_nmax_relative_energy(flux::TeukolskyPointParticleFlux, reached_key::Symbol, shell_key::Symbol, total_energy)
    haskey(flux.reached, reached_key) || return nothing
    flux.reached[reached_key] === nothing || return nothing
    nmax = _flux_get(flux.result, :n_max)
    nmax === nothing && return nothing

    last_shell_energy = _flux_last_shell_energy(flux.result, shell_key)
    last_shell_energy === nothing && return nothing
    denom = max(abs(total_energy), eps(Float64))
    return abs(last_shell_energy) / denom
end

function _flux_tolerance_print_value(flux::TeukolskyPointParticleFlux)
    nmax = _flux_get(flux.result, :n_max)
    nmax === nothing && return flux.tolerance

    parts = String[]
    rel_inf = _flux_hit_nmax_relative_energy(
        flux,
        :n_reached_inf,
        :infinity_energy_flux_list,
        flux.infinity_energy_flux,
    )
    rel_hor = _flux_hit_nmax_relative_energy(
        flux,
        :n_reached_hor,
        :horizon_energy_flux_list,
        flux.horizon_energy_flux,
    )

    rel_inf !== nothing && push!(parts, "infinity: $rel_inf")
    rel_hor !== nothing && push!(parts, "horizon: $rel_hor")
    isempty(parts) && return flux.tolerance

    suggested_nmax = max(nmax + 1, ceil(Int, 1.6 * nmax))
    ratios = join(parts, "; ")
    return "$(flux.tolerance) (hit nmax; last-shell relative energy = $ratios; increase nmax manually, e.g. nmax = $suggested_nmax)"
end

_flux_truncation_floor_print_value(flux::TeukolskyPointParticleFlux) = flux.truncation_floor

function Base.show(io::IO, ::MIME"text/plain", teuk_mode::TeukolskyPointParticleMode)
    s = teuk_mode.mode.s
    
    if s == -2
        suffix = "_inf" 
    elseif s == 2
        suffix = "_hor" 
    else
        error("Invalid spin weight s = $s. Teukolsky formalism requires s ∈ {-2, 2}.")
    end

    # Print structured output with frame-specific suffixes
    println(io, "TeukolskyPointParticleMode(")
    print(io, "    mode = "); show(io, "text/plain", teuk_mode.mode); println(io, ",")
    println(io, "    amplitude$suffix = $(teuk_mode.amplitude),")
    println(io, "    energy_flux$suffix = $(teuk_mode.energy_flux),")
    println(io, "    angular_momentum_flux$suffix = $(teuk_mode.angular_momentum_flux),")
    println(io, "    Carter_const_flux$suffix = $(teuk_mode.Carter_const_flux),")
    println(io, "    method = $(teuk_mode.method),")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", gsn_mode::GSNPointParticleMode)
    s = gsn_mode.mode.s
    
    if s == -2
        suffix = "_inf" 
    elseif s == 2
        suffix = "_hor" 
    else
        error("Invalid spin weight s = $s. GSN formalism requires s ∈ {-2, 2}.")
    end

    println(io, "GSNPointParticleMode(")
    print(io, "    mode = "); show(io, "text/plain", gsn_mode.mode); println(io, ",")
    println(io, "    amplitude$suffix = $(gsn_mode.amplitude),")
    println(io, "    energy_flux$suffix = $(gsn_mode.energy_flux),")
    println(io, "    angular_momentum_flux$suffix = $(gsn_mode.angular_momentum_flux),")
    println(io, "    Carter_const_flux$suffix = $(gsn_mode.Carter_const_flux),")
    println(io, "    method = $(gsn_mode.method),")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", flux::TeukolskyPointParticleFlux)
    println(io, "TeukolskyPointParticleFlux(")
    println(io, "    orbital_parameters(a = $(flux.a), p = $(flux.p), e = $(flux.e), x = $(flux.x)),")
    println(io, "    orbit_type = $(flux.orbit_type),")
    println(io, "    infinity_energy_flux = $(flux.infinity_energy_flux),")
    println(io, "    infinity_angular_momentum_flux = $(flux.infinity_angular_momentum_flux),")
    println(io, "    infinity_carter_constant_flux = $(flux.infinity_carter_constant_flux),")
    println(io, "    horizon_energy_flux = $(flux.horizon_energy_flux),")
    println(io, "    horizon_angular_momentum_flux = $(flux.horizon_angular_momentum_flux),")
    println(io, "    horizon_carter_constant_flux = $(flux.horizon_carter_constant_flux),")
    println(io, "    total_modes = $(flux.total_modes),")
    if haskey(flux.reached, :n_reached_inf) && haskey(flux.reached, :n_reached_hor)
        inf_reached = _flux_reached_print_value(flux, :n_reached_inf, flux.reached[:n_reached_inf])
        hor_reached = _flux_reached_print_value(flux, :n_reached_hor, flux.reached[:n_reached_hor])
        println(io, "    n_reached = (infinity = $inf_reached, horizon = $hor_reached),")
        if haskey(flux.reached, :convolution_integral)
            println(io, "    convolution_integral = $(flux.reached[:convolution_integral]),")
        end
    elseif haskey(flux.reached, :l_reached_inf) && haskey(flux.reached, :l_reached_hor)
        println(io, "    l_reached = (infinity = $(flux.reached[:l_reached_inf]), horizon = $(flux.reached[:l_reached_hor])),")
        if haskey(flux.reached, :convolution_integral)
            println(io, "    convolution_integral = $(flux.reached[:convolution_integral]),")
        end
    else
        for key in keys(flux.reached)
            println(io, "    $key = $(_flux_reached_print_value(flux, key, flux.reached[key])),")
        end
    end
    println(io, "    tolerance = $(_flux_tolerance_print_value(flux)),")
    println(io, "    truncation_floor = $(_flux_truncation_floor_print_value(flux)),")
    println(io, "    cost = $(flux.cost) seconds,")
    print(io, ")")
end

_flux_get(result, key::Symbol, default=nothing) = hasproperty(result, key) ? getproperty(result, key) : default

function _teukolsky_flux_orbit_type(e, x)
    if e == 0.0 && (x == 1.0 || x == -1.0)
        return :circular
    elseif x == 1.0 || x == -1.0
        return :eccentric
    elseif e == 0.0
        return :inclined
    else
        return :generic
    end
end

function _teukolsky_flux_geodesic_params(a, p, e, x)
    if x == -1.0
        return (-a, p, e, 1.0)
    end
    return (a, p, e, x)
end

function _teukolsky_flux_bound_orbit(a, p, e, x)
    ag, pg, eg, xg = _teukolsky_flux_geodesic_params(a, p, e, x)
    orbit = kerr_geo_orbit_type_metadata(ag, pg, eg, xg)
    return orbit.family == "Stable" ||
        (orbit.family == "Critical" && orbit.shape == "Circular" &&
         orbit.energy_regime == "Elliptic")
end

function _flux_tail_levin_branch_strategy(label::AbstractString, enabled, start_n)
    if enabled === true
        if start_n === nothing
            return "$label: ISEM adaptive trapezoidal for all computed n; tail ISEM adaptive Levin enabled but not triggered"
        end
        return "$label: ISEM adaptive trapezoidal for n < $start_n; tail ISEM adaptive Levin for n >= $start_n"
    end
    return "$label: ISEM adaptive trapezoidal for all n; tail ISEM adaptive Levin inactive"
end

function _flux_tail_levin_strategy(enabled, start_n)
    if enabled === true
        if start_n === nothing
            return "ISEM adaptive trapezoidal for all computed n; tail ISEM adaptive Levin enabled but not triggered"
        end
        return "ISEM adaptive trapezoidal for n < $start_n; tail ISEM adaptive Levin for n >= $start_n"
    end
    return "ISEM adaptive trapezoidal for all n; tail ISEM adaptive Levin inactive"
end

function _flux_tail_levin_same(result)
    return _flux_get(result, :tail_levin_infinity_enabled) == _flux_get(result, :tail_levin_horizon_enabled) &&
        _flux_get(result, :tail_levin_start_inf) == _flux_get(result, :tail_levin_start_hor)
end

function _teukolsky_flux_convolution_integral_strategy(result, orbit_type::Symbol)
    if orbit_type == :circular
        return "no convolution integral is needed"
    elseif orbit_type == :inclined
        return "ISEM adaptive trapezoidal for all modes"
    end
    if _flux_tail_levin_same(result)
        return _flux_tail_levin_strategy(
            _flux_get(result, :tail_levin_infinity_enabled),
            _flux_get(result, :tail_levin_start_inf),
        )
    end
    inf = _flux_tail_levin_branch_strategy(
        "infinity",
        _flux_get(result, :tail_levin_infinity_enabled),
        _flux_get(result, :tail_levin_start_inf),
    )
    hor = _flux_tail_levin_branch_strategy(
        "horizon",
        _flux_get(result, :tail_levin_horizon_enabled),
        _flux_get(result, :tail_levin_start_hor),
    )
    return "$inf; $hor"
end

function _teukolsky_flux_convolution_integral_info(result, orbit_type::Symbol)
    strategy = _teukolsky_flux_convolution_integral_strategy(result, orbit_type)
    if orbit_type == :circular || orbit_type == :inclined
        return (strategy = strategy,)
    end
    if _flux_tail_levin_same(result)
        return (
            strategy = strategy,
            tail_levin_nmin = _flux_get(result, :tail_levin_nmin),
            tail_levin_local_n = _flux_get(result, :tail_levin_local_n),
            tail_levin_max_depth = _flux_get(result, :tail_levin_max_depth),
        )
    end
    return (
        strategy = strategy,
        infinity = (
            tail_levin_enabled = _flux_get(result, :tail_levin_infinity_enabled),
            tail_levin_start_n = _flux_get(result, :tail_levin_start_inf),
        ),
        horizon = (
            tail_levin_enabled = _flux_get(result, :tail_levin_horizon_enabled),
            tail_levin_start_n = _flux_get(result, :tail_levin_start_hor),
        ),
        tail_levin_nmin = _flux_get(result, :tail_levin_nmin),
        tail_levin_local_n = _flux_get(result, :tail_levin_local_n),
        tail_levin_max_depth = _flux_get(result, :tail_levin_max_depth),
    )
end

function _teukolsky_flux_reached(result, orbit_type::Symbol)
    if orbit_type == :circular
        return (
            l_reached = _flux_get(result, :l_reached),
            l_reached_inf = _flux_get(result, :l_reached_inf),
            l_reached_hor = _flux_get(result, :l_reached_hor),
            convolution_integral = _teukolsky_flux_convolution_integral_info(result, orbit_type),
        )
    elseif orbit_type == :inclined
        return (
            k_reached_inf = _flux_get(result, :k_reached_inf),
            k_reached_hor = _flux_get(result, :k_reached_hor),
            convolution_integral = _teukolsky_flux_convolution_integral_info(result, orbit_type),
        )
    else
        return (
            n_reached_inf = _flux_get(result, :n_reached_inf),
            n_reached_hor = _flux_get(result, :n_reached_hor),
            convolution_integral = _teukolsky_flux_convolution_integral_info(result, orbit_type),
        )
    end
end


_pointparticle_float(a, p, e, x) = (Float64(a), Float64(p), Float64(e), Float64(x))

"""
    with_pointparticle_submission(f, a, p, e, x; Nmax=2^14, Kmax=2^12, stop_requested=nothing, lifecycle_observer=nothing)

Create a submission for the orbit `(a, p, e, x)` and the sampling-grid caps `Nmax`, `Kmax`, call `f(submission)`
and return its result. Mode and flux calls inside `f` that pass `submission=submission` share one read-only
trajectory and sampling geometry, while every admitted mode keeps its own private cache.
`stop_requested` is an optional zero-argument callback; when it returns `true` no further modes are admitted.
On normal return, error or cancellation, all admitted modes and registered tasks are drained before the internal
caches are released; returned results remain valid.
"""
function with_pointparticle_submission(f::Function,a,p,e,x;Nmax::Int=2^14,Kmax::Int=2^12,
        stop_requested=nothing,lifecycle_observer=nothing)
    # The point-particle routes run in Float64; any real input (e.g. an integer or BigFloat 0) is converted.
    a, p, e, x = _pointparticle_float(a, p, e, x)
    ModeSummation.with_public_submission(f,a,p,e,x;Nmax,Kmax,stop_requested,lifecycle_observer)
end
export with_pointparticle_submission

"""
    cancel_pointparticle_submission!(submission; reason="cooperative cancellation")

Stop admitting new modes and tasks to `submission`; modes already running finish normally, and later admission
attempts throw `PublicSubmissionCancelled` with `reason`. Return `nothing`.
"""
function cancel_pointparticle_submission!(scope;reason="cooperative cancellation")
    ModeSummation.request_submission_cancel!(scope;reason)
end
export cancel_pointparticle_submission!

"""
    spawn_pointparticle_submission!(f, submission)

Start `f()` on a new task registered with `submission` and return the `Task`. The task is joined before the
submission is released; spawning after admission has closed throws `PublicSubmissionCancelled`.
"""
function spawn_pointparticle_submission!(f::Function,scope)
    ModeSummation.spawn_submission_task!(f,scope)
end
export spawn_pointparticle_submission!

@doc raw"""
    Teukolsky_pointparticle_flux(a, p, e, x; tol=1e-8, lmax=30, nmax=500, kmax=20)

Compute total point-particle fluxes by automatically selecting the mode-summation strategy from the orbit type.
For eccentric and generic tails, `tail_levin=nothing` switches to adaptive Levin
for `n >= levin_nmin` once the uniform trapezoidal grid would be much more
expensive than the adaptive one, `true` forces adaptive Levin from `levin_nmin`,
and `false` keeps adaptive trapezoidal sampling. `levin_local_n` controls the local adaptive-Levin
interval count when that path is selected.
`minimum_consecutive` is the number of consecutive shells below the tolerance
required to stop a summation; the default `nothing` uses 5 for circular
equatorial orbits and 2 otherwise.
For generic orbits, `generic_beta=1e-4` raises the per-mode cutoff of the
n != 0 shells to at least `generic_beta * tol * abs(branch_total)`, and
`generic_deep=1e-5` lets an outward l walk stop as soon as one mode is at most
that fraction of its cutoff. Peak-search walks use the plain shell cutoff and
confirmation count, without the tail budget or the deep stop. Setting either
value to zero disables that rule.
The generic l walks of k != 0 branches start at the expected l peak, an m
shell that is still growing does not count towards the stop, and the positive-m
walk only stops more than 5 past the m with the largest shell.
Pass `submission` (see `with_pointparticle_submission`) to reuse the orbit's
sampling geometry across calls.
The keyword `lifecycle_observer` is a diagnostic hook and is not guaranteed to be stable.

Return a `TeukolskyPointParticleFlux` with the total energy, angular-momentum and
Carter-constant fluxes at infinity and at the horizon, or `nothing` (with a
warning) when `(a, p, e, x)` is not a bound orbit.
"""
function Teukolsky_pointparticle_flux(a, p, e, x; submission=nothing, tol = 1e-8, lmax = 30, nmax = 500, kmax = 20, minimum_consecutive = nothing, N = 64, N0 = N, K = 16, K0 = K, Nmax = (submission === nothing ? 2^14 : submission.Nmax), Kmax = (submission === nothing ? 2^12 : submission.Kmax), sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, truncation_floor = 1e-16, mode_abs_floor = truncation_floor, zero_low_flux = false, threaded_sampling = false, neg_branch_scale = 0.1, generic_beta = ModeSummation.DEFAULT_GENERIC_BETA, generic_deep = ModeSummation.DEFAULT_GENERIC_DEEP, tail_levin = nothing, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = mode_abs_floor, levin_local_n::Int = ConvolutionIntegrals.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N, levin_max_depth::Int = 8, info::Bool = false, stop_requested=nothing, lifecycle_observer=nothing)
    # The point-particle routes run in Float64; any real input (e.g. an integer or BigFloat 0) is converted.
    a, p, e, x = _pointparticle_float(a, p, e, x)
    if !_teukolsky_flux_bound_orbit(a, p, e, x)
        @warn "The specified parameters do not correspond to a bound orbit." a p e x
        return nothing
    end

    orbit_type = _teukolsky_flux_orbit_type(e, x)
    t0 = time()
    result = ConvolutionIntegrals.with_y_radial_info(info) do
        if orbit_type == :circular
            ModeSummation.circular_mode_summation(a, p; x, tol = tol, lmax = lmax, min_consecutive = something(minimum_consecutive, 5), submission=submission, stop_requested=stop_requested, lifecycle_observer=lifecycle_observer)
        elseif orbit_type == :eccentric
            ModeSummation.eccentric_mode_summation(a, p, e; x, N = N, N0 = N0, Nmax = Nmax, tol = tol, lmax = lmax, nmax = nmax, minimum_consecutive = something(minimum_consecutive, 2), sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, tail_levin = tail_levin === nothing ? :auto : tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_local_n = levin_local_n, levin_max_depth = levin_max_depth, submission=submission, stop_requested=stop_requested, lifecycle_observer=lifecycle_observer)
        elseif orbit_type == :inclined
            ModeSummation.inclined_mode_summation(a, p, x; K = K, K0 = K0, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, minimum_consecutive = something(minimum_consecutive, 2), sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, submission=submission, stop_requested=stop_requested, lifecycle_observer=lifecycle_observer)
        else
            ModeSummation.generic_mode_summation(a, p, e, x; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, nmax = nmax, minimum_consecutive = something(minimum_consecutive, 2), sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, neg_branch_scale = neg_branch_scale, beta=generic_beta, deep=generic_deep, tail_levin = tail_levin === nothing ? :auto : tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_local_n = levin_local_n, levin_max_depth = levin_max_depth, submission=submission, stop_requested=stop_requested, lifecycle_observer=lifecycle_observer)
        end
    end
    cost = time() - t0

    return TeukolskyPointParticleFlux(
        a,
        p,
        e,
        x,
        orbit_type,
        result.infinity_energy_flux,
        result.infinity_angular_momentum_flux,
        result.infinity_carter_constant_flux,
        result.horizon_energy_flux,
        result.horizon_angular_momentum_flux,
        result.horizon_carter_constant_flux,
        result.total_modes,
        _flux_get(result, :tolerance, tol),
        (
            infinity = _flux_get(result, :infinity_truncation_floor, mode_abs_floor),
            horizon = _flux_get(result, :horizon_truncation_floor, mode_abs_floor),
        ),
        _teukolsky_flux_reached(result, orbit_type),
        cost,
        result,
    )
end

# Orbit family of the mode cache for a single-mode call: ISEM methods classify the
# orbit with a 1e-12 tolerance on e and |x|, as their convolution routes do; the
# legacy methods use the submission's family.
function _pointparticle_cache_family(method, e, x, fallback::Symbol)
    method in ("auto", "isem_trapezoidal", "isem_levin") || return fallback
    radial_zero = isapprox(e, 0.0; atol = 1e-12)
    equatorial = isapprox(abs(x), 1.0; atol = 1e-12)
    return radial_zero ? (equatorial ? :circular : :inclined) :
        (equatorial ? :eccentric : :generic)
end

# For ISEM methods with 0 < |e| <= 1e-12 the convolution routes build their own
# e = 0 orbit, so they must not share the submission geometry built for the raw e.
_pointparticle_share_geometry(method, e) =
    !(method in ("auto", "isem_trapezoidal", "isem_levin") &&
        !isequal(e, 0.0) && isapprox(e, 0.0; atol = 1e-12))

@doc raw"""
    Teukolsky_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int, K::Int)

Compute the amplitude of the inhomogeneous Teukolsky solution _at infinity_ (for s = - 2) and _at horizon_ (for s = + 2) due to a point particle on a generic timelike bound orbit around a Kerr black hole,
with a spin parameter of `a`, for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `n` the radial index and `k` the polar index).
The orbit is specified by `p` the semi-latus rectum, `e` the eccentricity and `x` the inclination parameter ($x \equiv \cos \theta_{\mathrm{inc}}$).

In addition, we compute also the energy, angular momentum and Carter constant flux at infinity (for s = - 2) and the horizon (for s = + 2). Note that these values are formalism-independent.

The numerical method to compute the convolution integral is specified by `method` (default: `auto`).
Use `method = "isem_trapezoidal"` for the ISEM trapezoidal path and `method = "isem_levin"` for the ISEM Levin path.
Legacy non-ISEM paths remain available as `method = "trapezoidal"` and `method = "levin"`.
Both ISEM convolution paths use `ISEM` for their homogeneous radial input.
The ISEM trapezoidal path starts from `N = 256`, `K = 64` and doubles the grid until the mode converges to `tol`, up to `Nmax = 2^12` and `Kmax = 2^9`.
For adaptive ISEM Levin, `levin_max_depth` controls the maximum bisection depth; `Nmax` and `Kmax` are the grid caps of the trapezoidal and non-adaptive paths.
Pass `submission` (see `with_pointparticle_submission`) to reuse the orbit's sampling geometry across calls.
The keyword `lifecycle_observer` is a diagnostic hook and is not guaranteed to be stable.

Return a `TeukolskyPointParticleMode` holding the mode information, the Teukolsky amplitude, the three fluxes,
the trajectory, the Y solution, the spin-weighted spheroidal harmonic and the method settings used.
"""
function Teukolsky_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", submission=nothing, N::Int=-1, K::Int=-1, Nmax::Int = (submission === nothing ? 2^12 : submission.Nmax), Kmax::Int = (submission === nothing ? 2^9 : submission.Kmax), tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, truncation_floor::Float64 = 1e-16, mode_abs_floor::Float64 = truncation_floor, zero_low_flux::Bool = false, threaded_sampling::Bool = false, levin_max_depth::Int = 8, levin_adaptive::Bool = true, stop_requested=nothing, lifecycle_observer=nothing)
    # The point-particle routes run in Float64; any real input (e.g. an integer or BigFloat 0) is converted.
    a, p, e, x = _pointparticle_float(a, p, e, x)
    method = lowercase(String(method))
    if method == "auto"
        method = "isem_trapezoidal"
    end

    if method == "isem_trapezoidal"
        N = N <= 0 ? 256 : N
        K = K <= 0 ? 64 : K
    elseif method == "isem_levin"
        N = N <= 0 ? 256 : N
        K = K <= 0 ? 32 : K
    elseif method == "trapezoidal"
        N = N <= 0 ? 256 : N
        K = K <= 0 ? 64 : K
    elseif method == "levin"
        N = N <= 0 ? 256 : N
        K = K <= 0 ? 32 : K
    else
        error("Currently supported methods are \"auto\", \"isem_trapezoidal\", \"isem_levin\", \"trapezoidal\", and \"levin\".")
    end

    scope, owned = ModeSummation.ensure_public_submission(submission, a, p, e, x;
        Nmax=Nmax, Kmax=Kmax, stop_requested=stop_requested, lifecycle_observer=lifecycle_observer)
    return ModeSummation._with_public_scope(scope, owned) do
      output = nothing
      Y_mode = nothing
      try
        output = ModeSummation.with_submission_operation(scope) do
          cache_family = _pointparticle_cache_family(method, e, x, scope.family)
          share_geometry = _pointparticle_share_geometry(method, e)
          ModeSummation.with_submission_mode(scope; mode=(s,l,m,n,k), family=cache_family, share_geometry=share_geometry) do cache
        if method == "isem_trapezoidal"
            output = ConvolutionIntegrals.convolution_integral_trapezoidal_isem(a, p, e, x, s, l, m, n, k; N=N, K=K, Nmax=Nmax, Kmax=Kmax, tol=tol, sample_tol=sample_tol, max_flux=max_flux, mode_abs_floor=mode_abs_floor, zero_low_flux=zero_low_flux, threaded_sampling=threaded_sampling, cache=cache, trajectory=ModeSummation.submission_orbit!(scope))
        elseif method == "isem_levin"
            output = ConvolutionIntegrals.convolution_integral_levin_isem(a, p, e, x, s, l, m, n, k; N=N, K=K, Nmax=Nmax, Kmax=Kmax, tol=tol, sample_tol=sample_tol, max_flux=max_flux, mode_abs_floor=mode_abs_floor, zero_low_flux=zero_low_flux, threaded_sampling=threaded_sampling, adaptive_max_depth=levin_max_depth, adaptive=levin_adaptive, cache=cache, trajectory=ModeSummation.submission_orbit!(scope))
        elseif method == "trapezoidal"
            output = ConvolutionIntegrals.convolution_integral_trapezoidal(a, p, e, x, s, l, m, n, k; N=N, K=K)
        elseif method == "levin"
            output = ConvolutionIntegrals.convolution_integral_levin(a, p, e, x, s, l, m, n, k; N=N, K=K)
        end
            output
          end
        end

    Y_mode = _pointparticle_output_mode(output["YSolution"], s, l, m, a, output["omega"])
    return TeukolskyPointParticleMode(
            PointParticleMode(s, l, m, n, k, a, output["omega"], Y_mode.lambda),
            output["Amplitude"],
            output["EnergyFlux"],
            output["AngularMomentumFlux"],
            output["CarterConstantFlux"],
            output["Trajectory"],
            output["YSolution"],
            output["SWSH"],
            (method=method, radial_method="ISEM",
                radial_sfe=output["omega"] === nothing ? nothing :
                    ConvolutionIntegrals._inhomogeneous_small_frequency_expansion_control(l, a, output["omega"]), N=N, K=K,
                truncation_floor=mode_abs_floor)
        )
      finally
        # Inner do-block assignments can box this partial result. Release only
        # these local references on return/error; rich API objects stay intact.
        output = nothing
        Y_mode = nothing
      end
    end
end

function _pointparticle_output_mode(Y_solution, s, l, m, a, omega)
    if hasproperty(Y_solution, :mode)
        return Y_solution.mode
    elseif Y_solution isa AbstractDict
        if haskey(Y_solution, :mode)
            return Y_solution[:mode]
        elseif haskey(Y_solution, "mode")
            return Y_solution["mode"]
        elseif haskey(Y_solution, :params)
            return Y_solution[:params]
        elseif haskey(Y_solution, "params")
            return Y_solution["params"]
        end
    elseif Y_solution === nothing
        lambda = omega === nothing ? missing :
            _swsh_eigenvalue(s, l, m, a * omega)
        return (s = s, l = l, m = m, a = a, omega = omega, lambda = lambda)
    end
    error("Could not extract mode from YSolution output.")
end

@doc raw"""
    GSN_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int, K::Int)

Compute the amplitude of the inhomogeneous GSN solution _at infinity_ (for s = - 2) and _at horizon_ (for s = + 2) due to a point particle on a generic timelike bound orbit around a Kerr black hole,
with a spin parameter of `a`, for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `n` the radial index and `k` the polar index).
The orbit is specified by `p` the semi-latus rectum, `e` the eccentricity and `x` the inclination parameter ($x \equiv \cos \theta_{\mathrm{inc}}$).

In addition, we compute also the energy, angular momentum and Carter constant flux at infinity (for s = - 2) and the horizon (for s = + 2). Note that these values are formalism-independent.

The numerical method to compute the convolution integral is specified by `method` (default: `auto`).
Use `method = "isem_trapezoidal"` for the ISEM trapezoidal path and `method = "isem_levin"` for the ISEM Levin path.
Legacy non-ISEM paths remain available as `method = "trapezoidal"` and `method = "levin"`.
Both ISEM convolution paths use `ISEM` for their homogeneous radial input.
We sample the trajectory over a grid of size N x K, where N and K are the number of Chebyshev nodes in the radial and the polar direction, respectively.
Note that they must be powers of 2.
The remaining keywords are those of `Teukolsky_pointparticle_mode`.

Return a `GSNPointParticleMode`, whose amplitude is the Teukolsky amplitude divided by the GSN-to-Teukolsky
transmission conversion factor (`Btrans` for s = + 2, `Ctrans` for s = - 2).
"""
function GSN_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", submission=nothing, N::Int=-1, K::Int=-1, Nmax::Int = (submission === nothing ? 2^12 : submission.Nmax), Kmax::Int = (submission === nothing ? 2^9 : submission.Kmax), tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, truncation_floor::Float64 = 1e-16, mode_abs_floor::Float64 = truncation_floor, zero_low_flux::Bool = false, threaded_sampling::Bool = false, levin_max_depth::Int = 8, levin_adaptive::Bool = true, stop_requested=nothing, lifecycle_observer=nothing)
    Teukolsky_mode = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x; method=method, N=N, K=K, Nmax=Nmax, Kmax=Kmax, tol=tol, sample_tol=sample_tol, max_flux=max_flux, truncation_floor=truncation_floor, mode_abs_floor=mode_abs_floor, zero_low_flux=zero_low_flux, threaded_sampling=threaded_sampling, levin_max_depth=levin_max_depth, levin_adaptive=levin_adaptive, submission=submission, stop_requested=stop_requested, lifecycle_observer=lifecycle_observer)
    if iszero(Teukolsky_mode.amplitude)
        T_to_SN = one(Teukolsky_mode.amplitude)
    elseif s == 2
        T_to_SN = ConversionFactors.Btrans(s, m, a, Teukolsky_mode.mode.omega, Teukolsky_mode.mode.lambda)
    elseif s == -2
        T_to_SN = ConversionFactors.Ctrans(s, m, a, Teukolsky_mode.mode.omega, Teukolsky_mode.mode.lambda)
    end
    return GSNPointParticleMode(
        Teukolsky_mode.mode,
        Teukolsky_mode.amplitude / T_to_SN,
        Teukolsky_mode.energy_flux,
        Teukolsky_mode.angular_momentum_flux,
        Teukolsky_mode.Carter_const_flux,
        Teukolsky_mode.trajectory,
        Teukolsky_mode.Y_solution,
        Teukolsky_mode.SWSH,
        Teukolsky_mode.method
    )
end

include("QNM/QNM.jl")
using .QNM
export QNMMode, QNMBranch, QNMResult, QNMEstimate, QNMFailure,
    QNMEndpointResult,
    QNMPairResult
export ordinary, mirror, Ordinary, Mirror
export LeaverResult, ISEMValidationResult, ExcitationFactorResult
export qnm_frequency, qnm_sequence, validate_qnm_with_isem
export qnm_excitation_factor, qnm, qnm_pair
export angular_A_to_lambda, lambda_to_angular_A

end
