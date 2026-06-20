module GeneralizedSasakiNakamura

include("Homogeneous/Kerr.jl")
include("Homogeneous/Coordinates.jl")
include("Homogeneous/AsymptoticExpansionCoefficients.jl")
include("Homogeneous/InitialConditions.jl")
include("Homogeneous/ConversionFactors.jl")
include("Homogeneous/Transformation.jl")
include("Homogeneous/Potentials.jl")
include("Homogeneous/Solutions.jl")
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

function _adaptive_horizon_expansion_order(a, order::Int)
    spin = abs(float(real(a)))
    spin >= 1 && return max(order, 50)
    scale = log10(inv(max(1 - spin, eps(Float64))))
    return scale, max(order, _DEFAULT_horizon_expansion_order, ceil(Int, 10 * max(scale, 0) - 1e-12))
end

function _adaptive_infinity_expansion_order(omega, order::Int)
    freq = abs(omega)
    freq == 0 && return order
    scale = -log10(float(freq))
    return scale, max(order, _DEFAULT_infinity_expansion_order, ceil(Int, 10 * max(scale, 0) - 1e-12))
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

struct GSNRadialFunction
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
    GSN_solution # Store the *full* GSN solution where asymptotic solutions are smoothly attached
    normalization_convention::NormalizationConvention # The normalization convention used for the *stored* GSN solution
    method # The method used to solve the GSN equation
end

_display_bool_flag(x) = (x === missing || x === nothing) ? x : Bool(x)
_display_isem_N(x) = (x === missing || x === nothing) ? "adaptive" : string(x)

# Implement pretty-printing for GSNRadialFunction
# Mostly to suppress the printing of the numerical solution
function Base.show(io::IO, ::MIME"text/plain", gsn_func::GSNRadialFunction)
    println(io, "GSNRadialFunction(")
    print(io, "    mode = "); show(io, "text/plain", gsn_func.mode); println(io, ",")
    println(io, "    boundary_condition = $(gsn_func.boundary_condition),")
    if gsn_func.method == "ISEM" && gsn_func.numerical_GSN_solution isa NamedTuple
        params = gsn_func.numerical_GSN_solution
        N_display = hasproperty(params, :N) ? _display_isem_N(params.N) : "adaptive"
        println(io, "    matching_point(xm = $(params.xm), rhom = $(params.rhom)),")
        println(io, "    N = $(N_display),")
        println(io, "    tolerance = $(params.tol),")
        println(io, "    small/large_frequency_expansion(sfe = $(_display_bool_flag(params.sfe)), lfe = $(_display_bool_flag(params.lfe))),")
        println(io, "    Teukolsky-Starobinsky_identity(InInf = $(_display_bool_flag(params.TSinInf)), OutInf = $(_display_bool_flag(params.TSoutInf)), InHor = $(_display_bool_flag(params.TSinHor)), OutHor = $(_display_bool_flag(params.TSoutHor))),")
    else
        println(io, "    rsin = $(gsn_func.rsin),")
        println(io, "    rsout = $(gsn_func.rsout),")
        println(io, "    rsmp = $(gsn_func.rsmp),")
        println(io, "    horizon_expansion_order = $(gsn_func.horizon_expansion_order),")
        println(io, "    infinity_expansion_order = $(gsn_func.infinity_expansion_order),")
    end
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

struct TeukolskyRadialFunction
    mode::Mode # Information about the mode
    boundary_condition::BoundaryCondition # The boundary condition that this radial function statisfies
    transmission_amplitude # In Teukolsky formalism
    incidence_amplitude # In Teukolsky formalism
    reflection_amplitude # In Teukolsky formalism
    P_solution::Union{Function, Missing} # Store the full P solution if available
    GSN_solution::Union{GSNRadialFunction, Missing} # Store the full GSN solution if available
    Teukolsky_solution # Store the full Teukolsky solution
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
_use_isem_method(method) = method == "auto" || method == "ISEM"
const _STATIC_OMEGA_TOL = 1e-12
_is_static_frequency(omega) = abs(omega) < _STATIC_OMEGA_TOL
_is_horizon_superradiance_frequency(a, m, omega) = isreal(omega) && !_is_static_frequency(omega) && abs(omega - m * Kerr.omega_horizon(a)) < _STATIC_OMEGA_TOL

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

function _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    if boundary_condition == IN || boundary_condition == UP
        return ISEM.Teukolsky_radial(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    elseif boundary_condition == DOWN
        Rin, Rup = ISEM.Teukolsky_radial(s, l, m, a, omega; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        return _combine_isem_down(Rin, Rup)
    elseif boundary_condition == OUT
        Rin, Rup = ISEM.Teukolsky_radial(s, l, m, a, omega; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
        return _combine_isem_out(Rin, Rup)
    else
        error("Does not understand the boundary condition applied to the solution")
    end
end

function _gsn_from_isem(s, l, m, a, omega, boundary_condition; xm=nothing, rhom=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, use_gsn_asymptotic_patches=true, gsn_horizon_delta_r_max=ISEM._GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min=ISEM._GSN_INFINITY_PHASE_MIN)
    return ISEM.GSN_radial(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
end

function _try_legacy_riccati_then_linear(context, riccati_build, linear_build)
    saw_warn = Ref(false)
    logger = EarlyFilteredLogger(
        log -> begin
            if log.level == Logging.Warn
                saw_warn[] = true
            end
            return log.level >= Logging.Error
        end,
        current_logger(),
    )

    try
        sol = with_logger(logger) do
            riccati_build()
        end
        if saw_warn[]
            @info "method = \"auto\" tried Riccati, received a warning, and switched to method = \"linear\"." context=context
            return linear_build()
        end
        return sol
    catch err
        @info "method = \"auto\" tried Riccati, received an error, and switched to method = \"linear\"." context=context error=sprint(showerror, err)
        return linear_build()
    end
end

function _auto_try_isem_then_legacy(context, isem_build, riccati_build, linear_build)
    saw_warn = Ref(false)
    logger = EarlyFilteredLogger(
        log -> begin
            if log.level == Logging.Warn
                saw_warn[] = true
            end
            return log.level >= Logging.Error
        end,
        current_logger(),
    )

    try
        sol = with_logger(logger) do
            isem_build()
        end
        if saw_warn[]
            @info "method = \"auto\" tried ISEM, received an ISEM warning, and switched to legacy auto (Riccati, then linear). Use method = \"ISEM\" to force ISEM." context=context
            return _try_legacy_riccati_then_linear(context, riccati_build, linear_build)
        end
        return sol
    catch err
        @info "method = \"auto\" tried ISEM, received an ISEM error, and switched to legacy auto (Riccati, then linear). Use method = \"ISEM\" to force ISEM." context=context error=sprint(showerror, err)
        return _try_legacy_riccati_then_linear(context, riccati_build, linear_build)
    end
end

@doc raw"""
    GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order, method="auto", data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, rsmp=nothing)

Compute the GSN function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency [which *can be complex*]) 
and boundary condition specified by `boundary_condition`, which can be either

    - `IN` for purely-ingoing at the horizon,
    - `UP` for purely-outgoing at infinity,
    - `OUT` for purely-outgoing at the horizon,
    - `DOWN` for purely-ingoing at infinity.

Note that the `OUT` and `DOWN` solutions are constructed by linearly combining the `IN` and `UP` solutions, respectively.

The GSN function is numerically solved, for real values of `omega`, in the interval of *tortoise coordinates* $r_{*} \in$ `[rsin, rsout]` using the ODE solver (from `DifferentialEquations.jl`) specified by `ODE_algorithm` (default: `Vern9()`) 
with tolerance specified by `tolerance` (default: `1e-12`). The solution method is determined by the keyword `method` (default: `auto`).
The default `auto` uses `ISEM`; the legacy ODE paths remain available as `linear` and `Riccati`.
By default the data type used is `ComplexF64` (i.e. double-precision floating-point number) but it can be changed by
specifying `data_type` (e.g. `Complex{BigFloat}` for complex arbitrary precision number).

With complex values of `omega`, the GSN function is solved along a rotated path on the complex plane of $r_*$, 
where the path consists of two broken line segments parametrized by a real variable/a new coordinate $\rho$.
The angle with which the path is rotated is determined by the frequency $\omega$, the Kerr spin parameter $a$ and the azimuthal index $m$, 
such that the GSN function still behaves like a plane wave near the horizon and spatial infinity. 
At $\rho = 0$, the path intersects with the real axis at $r_{*} = r_{*}^{\rm{mp}}$ (`rsmp`, default to be nothing, which means it will be determined automatically).
Both the numerical and semi-analytical GSN solutions are evaluated as functions of $\rho$ instead of the now-complex $r_{*}$.

While the numerical GSN solution is only accurate in the range `[rsin, rsout]`, 
the full GSN solution is constructed by smoothly attaching the asymptotic solutions near horizon (up to `horizon_expansion_order`-th order) 
and infinity (up to `infinity_expansion_order`-th order). Therefore, the now-semi-analytical GSN solution is *accurate everywhere*.

Note, however, when `omega = 0`, the exact GSN function expressed using Gauss hypergeometric functions will be returned (i.e., instead of being solved numerically). 
In this case, only `s`, `l`, `m`, `a`, `omega`, `boundary_condition` will be parsed.

Return a `GSNRadialFunction` object which contains all the information about the GSN solution.
"""
function GSN_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order,
    method="auto", data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, tol=nothing, rsmp=nothing,
    xm=nothing, rhom=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, N=nothing, use_gsn_asymptotic_patches=true, gsn_horizon_delta_r_max=ISEM._GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min=ISEM._GSN_INFINITY_PHASE_MIN
)
    tolerance = tol === nothing ? tolerance : tol
    if _is_static_frequency(omega)
        return GSN_radial(s, l, m, a, zero(omega), boundary_condition)
    elseif _is_horizon_superradiance_frequency(a, m, omega) && boundary_condition == UP
        return _gsn_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tolerance, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    elseif _is_auto_method(method)
        context = (function_name = "GSN_radial", s = s, l = l, m = m, a = a, omega = omega, boundary_condition = boundary_condition)
        return _auto_try_isem_then_legacy(
            context,
            () -> _gsn_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tolerance, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min),
            () -> GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="Riccati", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp),
            () -> GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="linear", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp),
        )
    elseif method == "ISEM"
        return _gsn_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tolerance, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    else
        sawwarn = Ref(false)
        logger = EarlyFilteredLogger(
            (log)->(log.level == Logging.Warn && (sawwarn[] = true); log.level ≥ Logging.Error),
            current_logger())
        # Compute the SWSH eigenvalue
        lambda = spin_weighted_spheroidal_eigenvalue(s, l, m, a*omega)
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
                        @warn "There might be an issue with solving the Riccati form of the GSN equation. Consider using method=\"linear\" instead."
                    end
                    # Then convert to Xin
                    Xinsoln = Solutions.Xsoln_from_Phisoln(Phiinsoln)
                elseif method == "linear"
                    Phiinsoln = nothing
                    Xinsoln = Solutions.solve_Xin(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=horizon_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
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
                    Phiinsoln, _, _ = ComplexFrequencies.solve_Phiin(
                        s, m, a, -angle(omega), -angle(p),
                        omega, lambda, r_from_rho,
                        rsmp, rho_min, rho_max;
                        initialconditions_order=horizon_expansion_order, dtype=data_type,
                        odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance
                    )
                    # Then convert to Xin
                    Xinsoln = Solutions.Xsoln_from_Phisoln(Phiinsoln)
                elseif method == "linear"
                    Phiinsoln = nothing
                    Xinsoln, _, _ = ComplexFrequencies.solve_Xin(
                        s, m, a, -angle(omega), -angle(p),
                        omega, lambda, r_from_rho,
                        rsmp, rho_min, rho_max;
                        initialconditions_order=horizon_expansion_order, dtype=data_type,
                        odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance
                    )
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
                        @warn "There might be an issue with solving the Riccati form of the GSN equation. Consider using method=\"linear\" or method=\"ISEM\" instead."
                    end
                    # Then convert to Xup
                    Xupsoln = Solutions.Xsoln_from_Phisoln(Phiupsoln)
                elseif method == "linear"
                    Phiupsoln = nothing
                    Xupsoln = Solutions.solve_Xup(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=infinity_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
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
                    Phiupsoln, _, _ = ComplexFrequencies.solve_Phiup(
                        s, m, a, -angle(omega), -angle(p),
                        omega, lambda, r_from_rho,
                        rsmp, rho_min, rho_max;
                        initialconditions_order=infinity_expansion_order, dtype=data_type,
                        odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance
                    )
                    # Then convert to Xup
                    Xupsoln = Solutions.Xsoln_from_Phisoln(Phiupsoln)
                elseif method == "linear"
                    Phiupsoln = nothing
                    Xupsoln, _, _ = ComplexFrequencies.solve_Xup(
                        s, m, a, -angle(omega), -angle(p),
                        omega, lambda, r_from_rho,
                        rsmp, rho_min, rho_max;
                        initialconditions_order=infinity_expansion_order, dtype=data_type,
                        odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance
                    )
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

function GSN_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition;
    method="auto", tolerance=Solutions._DEFAULTTOLERANCE, tol=nothing,
    xm=nothing, rhom=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, N=nothing, use_gsn_asymptotic_patches=true, gsn_horizon_delta_r_max=ISEM._GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min=ISEM._GSN_INFINITY_PHASE_MIN
)
    tolerance = tol === nothing ? tolerance : tol
    if !_is_static_frequency(omega)
        return GSN_radial(s, l, m, a, omega, boundary_condition, _DEFAULT_rsin, _DEFAULT_rsout;
            method=method, tolerance=tolerance, xm=xm, rhom=rhom, sfe=sfe, lfe=lfe,
            TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, N=N, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
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
    GSN_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, method="auto")

Compute the GSN function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency)
with the purely-ingoing boundary condition at the horizon (`IN`) and the purely-outgoing boundary condition at infinity (`UP`).

Note that for _real_ frequencies, the numerical inner boundary (rsin) and outer boundary (rsout) are set to the default values `_DEFAULT_rsin` and `_DEFAULT_rsout`, respectively,
while the order of the asymptotic expansion at the horizon and infinity are determined automatically.
As for _complex_ frequencies, the numerical inner and the outer boundaries are determined automatically,
while the order of the asymptotic expansion at the horizon and infinity are set to `_DEFAULT_horizon_expansion_order_for_cplx_freq` and `_DEFAULT_infinity_expansion_order_for_cplx_freq`, respectively.
"""
function GSN_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, tol=nothing, method="auto", xm=nothing, rhom=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, N=nothing, use_gsn_asymptotic_patches=true, gsn_horizon_delta_r_max=ISEM._GSN_HORIZON_DELTA_R_MAX, gsn_infinity_phase_min=ISEM._GSN_INFINITY_PHASE_MIN)
    tolerance = tol === nothing ? tolerance : tol
    if _is_static_frequency(omega)
        omega0 = zero(omega)
        Xin = GSN_radial(s, l, m, a, omega0, IN)
        Xup = GSN_radial(s, l, m, a, omega0, UP)
    elseif _is_horizon_superradiance_frequency(a, m, omega)
        Xin = GSN_radial(s, l, m, a, omega, IN; method=method, tolerance=tolerance, xm=xm, rhom=rhom, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, N=N, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
        Xup = GSN_radial(s, l, m, a, omega, UP; method=method, tolerance=tolerance, xm=xm, rhom=rhom, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, N=N, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    elseif _is_auto_method(method)
        Xin = GSN_radial(s, l, m, a, omega, IN; method=method, tolerance=tolerance, xm=xm, rhom=rhom, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, N=N, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
        Xup = GSN_radial(s, l, m, a, omega, UP; method=method, tolerance=tolerance, xm=xm, rhom=rhom, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, N=N, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
    elseif method == "ISEM"
        Xin = _gsn_from_isem(s, l, m, a, omega, IN; xm=xm, rhom=rhom, N=N, tol=tolerance, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
        Xup = _gsn_from_isem(s, l, m, a, omega, UP; xm=xm, rhom=rhom, N=N, tol=tolerance, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, use_gsn_asymptotic_patches=use_gsn_asymptotic_patches, gsn_horizon_delta_r_max=gsn_horizon_delta_r_max, gsn_infinity_phase_min=gsn_infinity_phase_min)
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
        lambda = spin_weighted_spheroidal_eigenvalue(s, l, m, a*omega)
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
    Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order, method="auto", data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, rsmp=nothing)

Compute the Teukolsky function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency [which *can be complex*]) 
and boundary condition specified by `boundary_condition`, which can be either

    - `IN` for purely-ingoing at the horizon,
    - `UP` for purely-outgoing at infinity,
    - `OUT` for purely-outgoing at the horizon,
    - `DOWN` for purely-ingoing at infinity.

Note that the `OUT` and `DOWN` solutions are constructed by linearly combining the `IN` and `UP` solutions, respectively.

The full GSN solution is converted to the corresponding Teukolsky solution $(R(r), dR/dr)$ and 
the incidence, reflection and transmission amplitude are converted from the GSN formalism to the Teukolsky formalism 
with the normalization convention that the transmission amplitude is normalized to 1 (i.e. `normalization_convention=UNIT_TEUKOLSKY_TRANS`).
When `method` is `"auto"` or `"ISEM"`, the homogeneous Teukolsky solution is obtained directly from `ISEM`.

Note, however, when `omega = 0`, the exact Teukolsky function expressed using Gauss hypergeometric functions will be returned (i.e., instead of using the GSN formalism). 
In this case, only `s`, `l`, `m`, `a`, `omega`, `boundary_condition` will be parsed.

With complex values of `omega`, the Teukolsky function is evaluated as a function of $r$,
where the value at the corresponding location $\rho = r_{*}(r) \in \mathbb{R}$ along the rotated integration path on the complex plane is returned.
"""
function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order,
    method="auto", data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, tol=nothing, rsmp=nothing,
    xm=nothing, rhom=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, N=nothing
)
    tolerance = tol === nothing ? tolerance : tol
    if _is_static_frequency(omega)
        return Teukolsky_radial(s, l, m, a, zero(omega), boundary_condition; method="static")
    elseif _is_horizon_superradiance_frequency(a, m, omega) && boundary_condition == UP
        return _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tolerance, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    elseif _is_auto_method(method)
        context = (function_name = "Teukolsky_radial", s = s, l = l, m = m, a = a, omega = omega, boundary_condition = boundary_condition)
        return _auto_try_isem_then_legacy(
            context,
            () -> _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tolerance, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor),
            () -> Teukolsky_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="Riccati", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp),
            () -> Teukolsky_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method="linear", data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp),
        )
    elseif method == "ISEM"
        return _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tolerance, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    end

    # Solve for the GSN solution
    gsn_func = GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, method=method, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance, rsmp=rsmp)

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
    teuk_func(r) = Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function(
        gsn_func.mode.s,
        gsn_func.mode.m,
        gsn_func.mode.a,
        gsn_func.mode.omega,
        gsn_func.mode.lambda,
        gsn_func.GSN_solution
    )(r) / transmission_amplitude_conv_factor

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

function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition; method="auto", xm=nothing, rhom=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing, N=nothing, tol=nothing
)
    if _is_horizon_superradiance_frequency(a, m, omega) && boundary_condition == UP
        return _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    elseif !_is_static_frequency(omega) && _is_auto_method(method)
        context = (function_name = "Teukolsky_radial", s = s, l = l, m = m, a = a, omega = omega, boundary_condition = boundary_condition)
        return _auto_try_isem_then_legacy(
            context,
            () -> _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor),
            () -> Teukolsky_radial(s, l, m, a, omega, boundary_condition, _DEFAULT_rsin, _DEFAULT_rsout; method="Riccati", tolerance = tol === nothing ? Solutions._DEFAULTTOLERANCE : tol),
            () -> Teukolsky_radial(s, l, m, a, omega, boundary_condition, _DEFAULT_rsin, _DEFAULT_rsout; method="linear", tolerance = tol === nothing ? Solutions._DEFAULTTOLERANCE : tol),
        )
    elseif !_is_static_frequency(omega) && method == "ISEM"
        return _teukolsky_from_isem(s, l, m, a, omega, boundary_condition; xm=xm, rhom=rhom, N=N, tol=tol, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor)
    elseif !_is_static_frequency(omega)
        return Teukolsky_radial(s, l, m, a, omega, boundary_condition, _DEFAULT_rsin, _DEFAULT_rsout; method=method)
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
    Teukolsky_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, method="auto")

Compute the Teukolsky function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency)
with the purely-ingoing boundary condition at the horizon (`IN`) and the purely-outgoing boundary condition at infinity (`UP`).
When `method` is `"auto"` or `"ISEM"`, this overload uses `ISEM` for the homogeneous solutions.

Note that for _real_ frequencies, the numerical inner boundary (rsin) and outer boundary (rsout) are set to the default values `_DEFAULT_rsin` and `_DEFAULT_rsout`, respectively,
while the order of the asymptotic expansion at the horizon and infinity are determined automatically.
As for _complex_ frequencies, the numerical inner and the outer boundaries are determined automatically,
while the order of the asymptotic expansion at the horizon and infinity are set to `_DEFAULT_horizon_expansion_order_for_cplx_freq` and `_DEFAULT_infinity_expansion_order_for_cplx_freq`, respectively.
"""
function Teukolsky_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, tol=nothing, method="auto", xm=nothing, rhom=nothing, sfe=nothing, lfe=nothing, TSinInf=nothing, TSoutInf=nothing, TSinHor=nothing, TSoutHor=nothing)
    tolerance = tol === nothing ? tolerance : tol
    if _is_static_frequency(omega) || _is_horizon_superradiance_frequency(a, m, omega) || _use_isem_method(method)
        Rin = Teukolsky_radial(s, l, m, a, omega, IN; method=method, xm=xm, rhom=rhom, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, tol=tolerance)
        Rup = Teukolsky_radial(s, l, m, a, omega, UP; method=method, xm=xm, rhom=rhom, sfe=sfe, lfe=lfe, TSinInf=TSinInf, TSoutInf=TSoutInf, TSinHor=TSinHor, TSoutHor=TSoutHor, tol=tolerance)
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
(teuk_func::TeukolskyRadialFunction)(r) = teuk_func.Teukolsky_solution(r)[1] # Only return R(r), discarding the first derivative


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
    omega::Union{Real, Complex} # frequency
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
    for key in keys(flux.reached)
        println(io, "    $key = $(_flux_reached_print_value(flux, key, flux.reached[key])),")
    end
    println(io, "    tolerance = $(_flux_tolerance_print_value(flux)),")
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
    try
        KG = kerr_geo_orbit(ag, pg, eg, xg)
        return typeof(KG) != Vector{String}
    catch
        return false
    end
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

function _teukolsky_flux_convolution_integral_strategy(result, orbit_type::Symbol)
    if orbit_type == :circular || orbit_type == :inclined
        return "ISEM adaptive trapezoidal for all modes"
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
        tail_levin_max_depth = _flux_get(result, :tail_levin_max_depth),
    )
end

function _teukolsky_flux_reached(result, orbit_type::Symbol)
    if orbit_type == :circular
        return (
            l_reached = _flux_get(result, :l_reached),
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

@doc raw"""
    Teukolsky_pointparticle_flux(a, p, e, x; tol=1e-8, lmax=30, nmax=500, kmax=20)

Compute total point-particle fluxes by automatically selecting the mode-summation strategy from the orbit type.
"""
function Teukolsky_pointparticle_flux(a, p, e, x; tol = 1e-8, lmax = 30, nmax = 500, kmax = 20, minimum_consecutive = 2, N = 64, N0 = N, K = 16, K0 = K, Nmax = 2^14, Kmax = 2^12, sample_tol = 1e-3, record::Bool = false, record_path = nothing, fast = true, mode_abs_floor = 1e-16, zero_low_flux = false, threaded_sampling = false, neg_branch_scale = 0.1, tail_levin = nothing, tail_levin_infinity = nothing, tail_levin_horizon = nothing, levin_nmin = 50, levin_mode_abs_floor = 1e-16, levin_max_depth::Int = 8)
    if !_teukolsky_flux_bound_orbit(a, p, e, x)
        @warn "The specified parameters do not correspond to a bound orbit." a p e x
        return nothing
    end

    orbit_type = _teukolsky_flux_orbit_type(e, x)
    t0 = time()
    result = if orbit_type == :circular
        ModeSummation.circular_mode_summation(x == -1.0 ? -a : a, p; tol = tol, lmax = lmax, min_consecutive = minimum_consecutive)
    elseif orbit_type == :eccentric
        ModeSummation.eccentric_mode_summation(x == -1.0 ? -a : a, p, e; N = N, N0 = N0, Nmax = Nmax, tol = tol, lmax = lmax, nmax = nmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, tail_levin = tail_levin === nothing ? true : tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_max_depth = levin_max_depth)
    elseif orbit_type == :inclined
        ModeSummation.inclined_mode_summation(a, p, x; K = K, K0 = K0, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling)
    else
        ModeSummation.generic_mode_summation(a, p, e, x; N0 = N0, K0 = K0, Nmax = Nmax, Kmax = Kmax, tol = tol, lmax = lmax, kmax = kmax, nmax = nmax, minimum_consecutive = minimum_consecutive, sample_tol = sample_tol, record = record, record_path = record_path, fast = fast, mode_abs_floor = mode_abs_floor, zero_low_flux = zero_low_flux, threaded_sampling = threaded_sampling, neg_branch_scale = neg_branch_scale, tail_levin = tail_levin === nothing ? true : tail_levin, tail_levin_infinity = tail_levin_infinity, tail_levin_horizon = tail_levin_horizon, levin_nmin = levin_nmin, levin_mode_abs_floor = levin_mode_abs_floor, levin_max_depth = levin_max_depth)
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
        _teukolsky_flux_reached(result, orbit_type),
        cost,
        result,
    )
end

@doc raw"""
    Teukolsky_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int, K::Int)

Compute the amplitude of the inhomogeneous Teukolsky solution _at infinity_ (for s = - 2) and _at horizon_ (for s = + 2) due to a point particle on a generic timelike bound orbit around a Kerr black hole,
with a spin parameter of `a`, for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `n` the radial index and `k` the polar index).
The orbit is specified by `p` the semi-latus rectum, `e` the eccentricity and `x` the inclination parameter ($x \equiv \cos \theta_{\mathrm{inc}}$).

In addition, we compute also the energy, angular momentum and Carter constant flux at infinity (for s = - 2) and the horizon (for s = + 2). Note that these values are formalism-independent.

The numerical method to compute the convolution integral is specified by `method` (default: `auto`).
Use `method = "isem_trapezoidal"` for the ISEM trapezoidal path and `method = "isem_levin"` for the ISEM Levin path.
Legacy non-ISEM paths remain available as `method = "trapezoidal"` and `method = "levin"`.
For adaptive ISEM Levin, `levin_max_depth` controls the maximum bisection depth; `Nmax` and `Kmax` are fixed-grid caps used by trapezoidal and non-adaptive paths.
"""
function Teukolsky_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int=-1, K::Int=-1, Nmax::Int = 2^12, Kmax::Int = 2^9, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false, levin_max_depth::Int = 8)
    method = lowercase(String(method))
    if method == "auto"
        # For now, choose "isem_trapezoidal"
        method = "isem_trapezoidal"
    end

    # Fill-in default values for N/K if none is given
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
    end

    if method == "isem_trapezoidal"
        output = ConvolutionIntegrals.convolution_integral_trapezoidal_isem(a, p, e, x, s, l, m, n, k; N=N, K=K)
    elseif method == "isem_levin"
        output = ConvolutionIntegrals.convolution_integral_levin_isem(a, p, e, x, s, l, m, n, k; N=N, K=K, Nmax=Nmax, Kmax=Kmax, tol=tol, sample_tol=sample_tol, max_flux=max_flux, mode_abs_floor=mode_abs_floor, zero_low_flux=zero_low_flux, threaded_sampling=threaded_sampling, adaptive_max_depth=levin_max_depth)
    elseif method == "trapezoidal"
        output = ConvolutionIntegrals.convolution_integral_trapezoidal(a, p, e, x, s, l, m, n, k; N=N, K=K)
    elseif method == "levin"
        output = ConvolutionIntegrals.convolution_integral_levin(a, p, e, x, s, l, m, n, k; N=N, K=K)
    else
        error("Currently supported methods are \"auto\", \"isem_trapezoidal\", \"isem_levin\", \"trapezoidal\", and \"levin\".")
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
            (method=method, N=N, K=K)
        )
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
        lambda = omega === nothing ? missing : spin_weighted_spheroidal_eigenvalue(s, l, m, a * omega)
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
We sample the trajectory over a grid of size N x K, where N and K are the number of Chebyshev nodes in the radial and the polar direction, respectively.
Note that they must be powers of 2.
"""
function GSN_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int=-1, K::Int=-1, Nmax::Int = 2^12, Kmax::Int = 2^9, tol = 1e-8, sample_tol::Float64 = 1e-3, max_flux = 1.0, mode_abs_floor::Float64 = 1e-16, zero_low_flux::Bool = false, threaded_sampling::Bool = false, levin_max_depth::Int = 8)
    Teukolsky_mode = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x; method=method, N=N, K=K, Nmax=Nmax, Kmax=Kmax, tol=tol, sample_tol=sample_tol, max_flux=max_flux, mode_abs_floor=mode_abs_floor, zero_low_flux=zero_low_flux, threaded_sampling=threaded_sampling, levin_max_depth=levin_max_depth)
    if s == 2
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

end
