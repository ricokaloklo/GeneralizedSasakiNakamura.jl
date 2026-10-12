# Iterative series expansion matching (ISEM) for the generalized Sasaki-Nakamura (GSN) equation.
#
# The GSN equation is written in a compact coordinate on (r_+, infinity). Local series solutions
# at the horizon and at infinity (or, for large frequencies, their eikonal and large-frequency
# forms; for small frequencies the MST series) are continued by ordinary-point series and matched
# at an interior point, which gives the IN and UP solutions with their asymptotic amplitudes.
# Complex frequencies use the rotated-contour rational route of ComplexFrequency, cross-checked
# where needed against the P-equation solver in PEquation. Y_radial builds the auxiliary Y function
# used by the point-particle source integrals.
module ISEM

using ..GeneralizedSasakiNakamura
using ..Coordinates
using ..ConversionFactors
using ..Solutions

include("TeukolskyTransformation.jl")
using .TeukolskyTransformation

include("Parameters.jl")
using .Parameters
include("CoefficientTables.jl")
using .CoefficientTables
include("GeneratedBackend.jl")
using .GeneratedBackend
include("OrdinaryPointExpansion.jl")
using .OrdinaryPointExpansion
include("LargeFrequencyExpansion.jl")
using .LargeFrequencyExpansion
include("NearExtreme.jl")
using .NearExtreme
include("Eikonal.jl")
using .Eikonal
include("LocalSolutionAtZero.jl")
using .LocalSolutionAtZero
include("HorizonLargeFrequencyExpansion.jl")
using .HorizonLargeFrequencyExpansion
include("LocalSolutionAtInfinity.jl")
using .LocalSolutionAtInfinity
include("InfinityLargeFrequencyExpansion.jl")
using .InfinityLargeFrequencyExpansion
include("Iteration.jl")
using .Iteration
include("ManoSuzukiTakasugiInfinity.jl")
using .ManoSuzukiTakasugiInfinity
include("Matching.jl")
using .Matching
include("AsymptoticAmplitudes.jl")
using .AsymptoticAmplitudes
include("RadialFunctions.jl")
using .RadialFunctions
include("ComplexRational.jl")
using .ComplexRational
include("ComplexFrequency.jl")
using .ComplexFrequency
include("NegativeImaginaryAxisAmplitudes.jl")
using .NegativeImaginaryAxisAmplitudes

register_coefficient_backend!(GeneratedBackend.generated_coefficients)

include("PEquation/PEquation.jl")
using .PEquation: PEquation

export BoundaryCondition, IN, UP, OUT, DOWN
export NormalizationConvention, UNIT_GSN_TRANS, UNIT_TEUKOLSKY_TRANS
export Mode, GSNRadialFunction, TeukolskyRadialFunction
export YRadialFunction, Y_radial

# Frequencies within this distance of zero count as zero when choosing the mirror mode.
const _STATIC_OMEGA_TOL = 1e-12

_is_static_frequency(omega) = iszero(omega)
_is_exact_extremal_spin(a) = isreal(a) && isone(abs(real(a)))

# Y is defined for s = -2 with IN and s = +2 with UP only.
function _validate_y_branch(s, boundary_condition)
    if !(s == -2 && boundary_condition == IN || s == 2 && boundary_condition == UP)
        error("Y_radial currently supports only s = -2 with IN and s = +2 with UP.")
    end
end

# r-derivative of the Y phase, K / Delta with K = (r^2 + a^2) omega - a m.
function _y_phase_derivative(r, a, m, omega)
    return ((r^2 + a^2) * omega - a * m) / (r^2 - 2 * r + a^2)
end

# Y phase at |a| = 1 (Delta = (r - 1)^2): omega r* + a m / (r - 1), with
# r* = r + 2 log(r - 1) - 2 / (r - 1) - log 4.
function _extremal_y_phase(r, a, m, omega)
    h = r - one(r)
    h > zero(h) || throw(DomainError(r,
        "the exact-extremal source-adapted phase requires r > r_+ = 1"))
    detuning = 2omega - a * m
    return omega * r - omega * log(oftype(r, 4)) +
        2omega * log(h) - detuning / h
end

# Teukolsky transmission amplitude of the Y branch: B_trans for s = -2 (IN), C_trans for s = +2 (UP).
function _y_branch_transmission_amplitude(s, m, a, omega, lambda)
    if s == -2
        return ConversionFactors.Btrans(s, m, a, omega, lambda)
    elseif s == 2
        return ConversionFactors.Ctrans(s, m, a, omega, lambda)
    else
        error("Unsupported Y branch.")
    end
end

"""
    YRadialFunction

Callable auxiliary ISEM radial solution returned by `Y_radial`: `Y(r)` evaluates the solution at radius `r`.
Stores the mode, boundary condition, asymptotic amplitudes, the evaluators of the underlying P, Teukolsky and X
solutions where available, and the normalization convention.
"""
struct YRadialFunction
    mode::Mode
    boundary_condition::BoundaryCondition
    transmission_amplitude
    incidence_amplitude
    reflection_amplitude
    P_solution::Union{Function, Missing}
    Teukolsky_solution::Union{Function, Missing}
    X_solution::Union{Function, Missing}
    Y_scalar_solution::Function
    Y_solution::Function
    normalization_convention::NormalizationConvention
end

function Base.show(io::IO, ::MIME"text/plain", y_func::YRadialFunction)
    println(io, "YRadialFunction(")
    print(io, "    mode = "); show(io, "text/plain", y_func.mode); println(io, ",")
    println(io, "    boundary_condition = $(y_func.boundary_condition),")
    println(io, "    transmission_amplitude = $(y_func.transmission_amplitude),")
    println(io, "    incidence_amplitude = $(y_func.incidence_amplitude),")
    println(io, "    reflection_amplitude = $(y_func.reflection_amplitude),")
    println(io, "    normalization_convention = $(y_func.normalization_convention)")
    print(io, ")")
end

function Base.show(io::IO, y_func::YRadialFunction)
    print(io, "YRadialFunction(mode = Mode(s = $(y_func.mode.s), l = $(y_func.mode.l), m = $(y_func.mode.m), a = $(y_func.mode.a), omega = $(y_func.mode.omega), lambda = $(y_func.mode.lambda)), boundary_condition = $(y_func.boundary_condition))")
end

# Frequencies with Re(omega) < 0 that are real, or nearly real with |Re(omega)| <= 1, are built from
# the mirror mode (m -> -m, omega -> -conj(omega)) by complex conjugation.
function _is_reflected_lhp_frequency(omega)
    real_part = real(omega)
    real_part < -_STATIC_OMEGA_TOL || return false
    isreal(omega) && return true
    return abs(real_part) <= 1.0 &&
        abs(imag(omega)) <= 1.0e-3 * abs(real_part)
end
_conjugate_or_missing(x) = x === missing ? missing : conj(x)
_conjugate_solution_value(x) = x isa Tuple ? map(_conjugate_solution_value, x) : _conjugate_or_missing(x)
_conjugate_solution_function(f) = f === missing ? missing : (args...; kwargs...) -> _conjugate_solution_value(f(args...; kwargs...))

struct ConjugatedTeukolskySolution{F,P,R}
    full_value::F
    pair_value::P
    radial_value::R
end

(solution::ConjugatedTeukolskySolution)(args...; kwargs...) =
    solution.full_value(args...; kwargs...)

function _conjugate_teukolsky_solution(solution)
    solution === missing && return missing
    full_value = _conjugate_solution_function(solution)
    hasproperty(solution, :pair_value) && hasproperty(solution, :radial_value) ||
        return full_value
    return ConjugatedTeukolskySolution(
        full_value,
        _conjugate_solution_function(solution.pair_value),
        _conjugate_solution_function(solution.radial_value),
    )
end

function _conjugate_teukolsky_radial_function(f::TeukolskyRadialFunction, m::Int, omega)
    mode = Mode(f.mode.s, f.mode.l, m, f.mode.a, omega, conj(f.mode.lambda))
    return TeukolskyRadialFunction(
        mode,
        f.boundary_condition,
        _conjugate_or_missing(f.transmission_amplitude),
        _conjugate_or_missing(f.incidence_amplitude),
        _conjugate_or_missing(f.reflection_amplitude),
        _conjugate_solution_function(f.P_solution),
        f.GSN_solution === missing ? missing : _conjugate_gsn_radial_function(f.GSN_solution, m, omega),
        _conjugate_teukolsky_solution(f.Teukolsky_solution),
        f.normalization_convention,
    )
end

function _conjugate_gsn_radial_function(f::GSNRadialFunction, m::Int, omega)
    mode = Mode(f.mode.s, f.mode.l, m, f.mode.a, omega, conj(f.mode.lambda))
    return GSNRadialFunction(
        mode,
        f.boundary_condition,
        f.rsin,
        f.rsout,
        f.rsmp,
        f.horizon_expansion_order,
        f.infinity_expansion_order,
        _conjugate_or_missing(f.transmission_amplitude),
        _conjugate_or_missing(f.incidence_amplitude),
        _conjugate_or_missing(f.reflection_amplitude),
        f.numerical_GSN_solution,
        _conjugate_solution_function(f.numerical_Riccati_solution),
        _conjugate_solution_function(f.GSN_solution),
        f.normalization_convention,
        f.method,
    )
end

function _conjugate_y_radial_function(f::YRadialFunction, m::Int, omega)
    mode = Mode(f.mode.s, f.mode.l, m, f.mode.a, omega, conj(f.mode.lambda))
    return YRadialFunction(
        mode,
        f.boundary_condition,
        _conjugate_or_missing(f.transmission_amplitude),
        _conjugate_or_missing(f.incidence_amplitude),
        _conjugate_or_missing(f.reflection_amplitude),
        _conjugate_solution_function(f.P_solution),
        _conjugate_teukolsky_solution(f.Teukolsky_solution),
        _conjugate_solution_function(f.X_solution),
        _conjugate_solution_function(f.Y_scalar_solution),
        _conjugate_solution_function(f.Y_solution),
        f.normalization_convention,
    )
end

function _spin_reflection_metadata(metadata, source_a, source_m, target_a, target_m)
    metadata isa NamedTuple || return metadata
    return merge(metadata, (
        spin_reflection_applied=true,
        spin_reflection_map=:a_m_to_minus_a_minus_m,
        spin_reflection_source_a=source_a,
        spin_reflection_source_m=source_m,
        spin_reflection_target_a=target_a,
        spin_reflection_target_m=target_m,
    ))
end

function _spin_reflect_teukolsky_radial_function(
        f::TeukolskyRadialFunction, m::Int, a)
    mode = Mode(f.mode.s, f.mode.l, m, a, f.mode.omega, f.mode.lambda)
    return TeukolskyRadialFunction(
        mode,
        f.boundary_condition,
        f.transmission_amplitude,
        f.incidence_amplitude,
        f.reflection_amplitude,
        f.P_solution,
        f.GSN_solution,
        f.Teukolsky_solution,
        f.normalization_convention,
    )
end


function _spin_reflect_gsn_radial_function(f::GSNRadialFunction, m::Int, a)
    mode = Mode(f.mode.s, f.mode.l, m, a, f.mode.omega, f.mode.lambda)
    metadata = _spin_reflection_metadata(
        f.numerical_GSN_solution, f.mode.a, f.mode.m, a, m)
    return GSNRadialFunction(
        mode,
        f.boundary_condition,
        f.rsin,
        f.rsout,
        f.rsmp,
        f.horizon_expansion_order,
        f.infinity_expansion_order,
        f.transmission_amplitude,
        f.incidence_amplitude,
        f.reflection_amplitude,
        metadata,
        f.numerical_Riccati_solution,
        f.GSN_solution,
        f.normalization_convention,
        f.method,
    )
end


function _spin_reflect_y_radial_function(f::YRadialFunction, m::Int, a)
    mode = Mode(f.mode.s, f.mode.l, m, a, f.mode.omega, f.mode.lambda)
    return YRadialFunction(
        mode,
        f.boundary_condition,
        f.transmission_amplitude,
        f.incidence_amplitude,
        f.reflection_amplitude,
        f.P_solution,
        f.Teukolsky_solution,
        f.X_solution,
        f.Y_scalar_solution,
        f.Y_solution,
        f.normalization_convention,
    )
end

_kerr_delta(r, a) = r^2 + a^2 - 2 * r

function _y_fallback_horizon_expansion_order(a)
    denom = max(1 - abs(float(a)), eps(Float64))
    return max(20, ceil(Int, 5 * max(log10(inv(denom)), 0)))
end

function _y_fallback_infinity_expansion_order(omega)
    ω = abs(float(real(omega)))
    if ω <= 0
        return 30
    end
    return max(20, ceil(Int, 10 * max(-log10(ω), 0) + 10))
end

# Y from a GSN solution: the exact-extremal GSN route when |a| = 1, otherwise the Riccati or
# linear integration of the GSN equation.
function _build_y_solution_from_gsn(s, l, m, a, omega, boundary_condition;
        method="Riccati", tol=nothing, allow_complex::Bool=false)
    _validate_y_branch(s, boundary_condition)
    iszero(omega) && error("Y_radial currently requires nonzero omega.")
    (allow_complex || isreal(omega)) ||
        error("Y_radial currently supports complex omega only on the exact-extremal ISEM GSN route.")

    parent = parentmodule(@__MODULE__)
    X_func = if _is_exact_extremal_spin(a)
        getfield(parent, :GSN_radial)(
            s, l, m, a, omega, boundary_condition;
            method=method,
            tolerance=tol === nothing ? Solutions._DEFAULTTOLERANCE : tol,
        )
    else
        getfield(parent, :GSN_radial)(
            s,
            l,
            m,
            a,
            omega,
            boundary_condition,
            getfield(parent, :_DEFAULT_rsin),
            getfield(parent, :_DEFAULT_rsout);
            method=method,
            tolerance = tol === nothing ? Solutions._DEFAULTTOLERANCE : tol,
            horizon_expansion_order = _y_fallback_horizon_expansion_order(a),
            infinity_expansion_order = _y_fallback_infinity_expansion_order(omega),
        )
    end
    lambda = X_func.mode.lambda
    branch_transmission_amplitude = _y_branch_transmission_amplitude(s, m, a, omega, lambda)
    teukolsky_from_gsn_matrix = r -> Solutions.Teukolsky_radial_function_from_Sasaki_Nakamura_function_conversion_matrix(s, m, a, omega, lambda, r)

    if _is_exact_extremal_spin(a)
        function teukolsky_state(r)
            rs = rstar_from_r(a, r)
            return teukolsky_from_gsn_matrix(r) * X_func.GSN_solution(rs)
        end
        Xsoln = r -> X_func.GSN_solution(rstar_from_r(a, r))[1]
        function Ysoln(r)
            radial = teukolsky_state(r)
            R, dRdr = radial
            phase = _extremal_y_phase(r, a, m, omega)
            phase_derivative = _y_phase_derivative(r, a, m, omega)
            if s == -2
                prefactor = branch_transmission_amplitude * exp(-im * phase)
                Y = prefactor * R / r^2
                dYdr = prefactor * (dRdr / r^2 - 2R / r^3 -
                    im * phase_derivative * R / r^2)
            else
                delta = _kerr_delta(r, a)
                delta_derivative = 2r - 2
                prefactor = branch_transmission_amplitude * exp(im * phase)
                Y = prefactor * delta^2 * R / r^2
                dYdr = prefactor * (
                    delta^2 * dRdr / r^2 +
                    (2delta * delta_derivative / r^2 -
                     2delta^2 / r^3 +
                     im * phase_derivative * delta^2 / r^2) * R)
            end
            return (Y, dYdr, Xsoln(r), zero(abs(Y)))
        end
        Yscalar = r -> Ysoln(r)[1]
        conv_incidence = s == -2 ?
            ConversionFactors.Binc(s, m, a, omega, lambda) :
            ConversionFactors.Cinc(s, m, a, omega, lambda)
        conv_reflection = s == -2 ?
            ConversionFactors.Bref(s, m, a, omega, lambda) :
            ConversionFactors.Cref(s, m, a, omega, lambda)
        incidence_amplitude = X_func.incidence_amplitude === missing ?
            missing : conv_incidence * X_func.incidence_amplitude
        reflection_amplitude = X_func.reflection_amplitude === missing ?
            missing : conv_reflection * X_func.reflection_amplitude
        mode = Mode(s, l, m, a, omega, lambda)
        return YRadialFunction(
            mode, boundary_condition, branch_transmission_amplitude,
            incidence_amplitude, reflection_amplitude, missing,
            teukolsky_state, Xsoln, Yscalar, Ysoln,
            UNIT_TEUKOLSKY_TRANS)
    end

    Ysoln = TeukolskyTransformation.GSN_to_Y_solution_from_matrix(
        r -> rstar_from_r(a, r),
        teukolsky_from_gsn_matrix,
        X_func.GSN_solution,
        s,
        m,
        a,
        omega,
        lambda,
        branch_transmission_amplitude,
    )

    conv_incidence = s == -2 ?
        ConversionFactors.Binc(s, m, a, omega, lambda) :
        ConversionFactors.Cinc(s, m, a, omega, lambda)
    conv_reflection = s == -2 ?
        ConversionFactors.Bref(s, m, a, omega, lambda) :
        ConversionFactors.Cref(s, m, a, omega, lambda)
    incidence_amplitude = X_func.incidence_amplitude === missing ? missing : conv_incidence * X_func.incidence_amplitude
    reflection_amplitude = X_func.reflection_amplitude === missing ? missing : conv_reflection * X_func.reflection_amplitude
    Xsoln = r -> X_func.GSN_solution(rstar_from_r(a, r))[1]
    Yscalar = r -> Ysoln(r)[1]

    mode = Mode(s, l, m, a, omega, lambda)
    return YRadialFunction(mode, boundary_condition, branch_transmission_amplitude, incidence_amplitude, reflection_amplitude, missing, missing, Xsoln, Yscalar, Ysoln, UNIT_TEUKOLSKY_TRANS)
end

"""
    Y_radial(s, l, m, a, omega, boundary_condition; method="auto", xm=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing)

Construct the auxiliary radial function Y used by the point-particle convolution integrals for the mode
`(s, l, m, a, omega)` and return a `YRadialFunction`; only `s = -2` with `IN` and `s = +2` with `UP` are
supported, and `omega` must be nonzero. `method = "auto"` or `"ISEM"` builds Y with ISEM; `"Riccati"` or
`"linear"` build it from a numerically integrated GSN solution. Negative `a` and frequencies with negative
real part are obtained from the mirror modes. `xm`, `N`, `tol` and the expansion switches are passed to ISEM
as in `GSN_radial`.
"""
function Y_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition::BoundaryCondition; method="auto", xm=nothing, N=nothing, tol=nothing, sfe=nothing, lfe=nothing)
    GeneralizedSasakiNakamura._check_radial_options(method, xm)
    if a < zero(a)
        positive = Y_radial(
            s, l, -m, -a, omega, boundary_condition;
            method, xm, N, tol, sfe, lfe)
        return _spin_reflect_y_radial_function(positive, m, a)
    end
    if _is_static_frequency(omega)
        _validate_y_branch(s, boundary_condition)
        error("Y_radial currently requires nonzero omega.")
    elseif _is_exact_extremal_spin(a)
        return _build_y_solution_from_gsn(
            s, l, m, a, omega, boundary_condition;
            method=method == "auto" ? "ISEM" : method, tol=tol, allow_complex=true)
    elseif _is_reflected_lhp_frequency(omega)
        positive_func = Y_radial(s, l, -m, a, -conj(omega), boundary_condition; method=method, xm=xm, N=N, tol=tol, sfe=sfe, lfe=lfe)
        return _conjugate_y_radial_function(positive_func, m, omega)
    end

    if method == "auto" || method == "ISEM"
        _validate_y_branch(s, boundary_condition)
        branch = boundary_condition == IN ? :IN : :UP
        route = radial_route(s, l, m, a, omega, branch;
            N=N, xm=xm, tol=tol,
            small_frequency_expansion=sfe === nothing ? :auto : sfe,
            large_frequency_expansion=lfe === nothing ? :auto : lfe)
        return y_radial_function(route)
    elseif method == "Riccati" || method == "linear"
        return _build_y_solution_from_gsn(s, l, m, a, omega, boundary_condition; method=method, tol=tol)
    else
        error("Method must be 'auto', 'ISEM', 'Riccati', or 'linear'")
    end
end

(y_func::YRadialFunction)(coord) = y_func.Y_scalar_solution(coord)

end
