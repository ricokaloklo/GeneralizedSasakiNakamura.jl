module GeneralizedSasakiNakamura

include("Kerr.jl")
include("Coordinates.jl")
include("AsymptoticExpansionCoefficients.jl")
include("InitialConditions.jl")
include("ConversionFactors.jl")
include("Potentials.jl")
include("Transformation.jl")
include("Solutions.jl")

using .Coordinates
export r_from_rstar, rstar_from_r # Useful to be exposed
using .Solutions

using SpinWeightedSpheroidalHarmonics
using DifferentialEquations # Should have been compiled by now

export GSN_radial, Teukolsky_radial

# IN for purely-ingoing at the horizon and UP for purely-outgoing at infinity
@enum BoundaryCondition begin
    IN = 1
    UP = 2
end
export IN, UP # Use these to specify the BC

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
    omega # frequency
    lambda # SWSH eigenvalue
end

# Implement pretty-printing for Mode type
# REPL
function Base.show(io::IO, ::MIME"text/plain", mode::Mode)
    print(io, "Mode(s=$(mode.s), l=$(mode.l), m=$(mode.m), a=$(mode.a), omega=$(mode.omega), lambda=$(mode.lambda))")
end


struct GSNRadialFunction
    mode::Mode # Information about the mode
    boundary_condition::BoundaryCondition # The boundary condition that this radial function statisfies
    rsin # The numerical inner boundary where the GSN equation is numerically evolved
    rsout # The numerical outer boundary where the GSN equation is numerically evolved
    horizon_expansion_order::Int # The order of the asymptotic expansion at the horizon
    infinity_expansion_order::Int # The order of the asymptotic expansion at infinity
    transmission_amplitude # In GSN formalism
    incidence_amplitude # In GSN formalism
    reflection_amplitude # In GSN formalism
    numerical_GSN_solution::Union{ODESolution, Nothing} # Store the numerical solution to the GSN equation in [rsin, rsout]
    numerical_Riccati_solution::Union{ODESolution, Nothing} # Store the numerical solution to the GSN equation in the Riccati form if applicable
    GSN_solution # Store the *full* GSN solution where asymptotic solutions are smoothly attached
    normalization_convention::NormalizationConvention # The normalization convention used for the *stored* GSN solution
end

# Implement pretty-printing for GSNRadialFunction
# Mostly to suppress the printing of the numerical solution
function Base.show(io::IO, ::MIME"text/plain", gsn_func::GSNRadialFunction)
    println(io, "GSNRadialFunction(")
    print(io, "    mode="); show(io, "text/plain", gsn_func.mode); println(io, ",")
    println(io, "    boundary_condition=$(gsn_func.boundary_condition),")
    println(io, "    rsin=$(gsn_func.rsin),")
    println(io, "    rsout=$(gsn_func.rsout),")
    println(io, "    horizon_expansion_order=$(gsn_func.horizon_expansion_order),")
    println(io, "    infinity_expansion_order=$(gsn_func.infinity_expansion_order),")
    println(io, "    transmission_amplitude=$(gsn_func.transmission_amplitude),")
    println(io, "    incidence_amplitude=$(gsn_func.incidence_amplitude),")
    println(io, "    reflection_amplitude=$(gsn_func.reflection_amplitude),")
    println(io, "    normalization_convention=$(gsn_func.normalization_convention)")
    print(io, ")")
end

struct TeukolskyRadialFunction
    mode::Mode # Information about the mode
    boundary_condition::BoundaryCondition # The boundary condition that this radial function statisfies
    transmission_amplitude # In Teukolsky formalism
    incidence_amplitude # In Teukolsky formalism
    reflection_amplitude # In Teukolsky formalism
    GSN_solution::GSNRadialFunction # Store the full GSN solution
    Teukolsky_solution # Store the full Teukolsky solution
    normalization_convention::NormalizationConvention # The normalization convention used for the *stored* Teukolsky solution
end

# Implement pretty-printing for TeukolskyRadialFunction
# Mostly to suppress the printing of the numerical solution
function Base.show(io::IO, ::MIME"text/plain", teuk_func::TeukolskyRadialFunction)
    println(io, "TeukolskyRadialFunction(")
    print(io, "    mode="); show(io, "text/plain", teuk_func.mode); println(io, ",")
    println(io, "    boundary_condition=$(teuk_func.boundary_condition),")
    println(io, "    transmission_amplitude=$(teuk_func.transmission_amplitude),")
    println(io, "    incidence_amplitude=$(teuk_func.incidence_amplitude),")
    println(io, "    reflection_amplitude=$(teuk_func.reflection_amplitude),")
    println(io, "    normalization_convention=$(teuk_func.normalization_convention)")
    print(io, ")")
end

@doc raw"""
    GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=3, infinity_expansion_order::Int=6, data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE)

Compute the GSN function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency) 
and boundary condition (specified by `boundary_condition` which can be either `IN` for purely-ingoing at the horizon or `UP` for purely-outgoing at infinity).

The GSN function is numerically solved in the interval `[rsin, rsout]` using the ODE solver (from `DifferentialEquations.jl`) specified by `ODE_algorithm` (default: `Vern9()`) 
with tolerance specified by `tolerance` (default: `1e-12`). By default the data type used is `ComplexF64` (i.e. double-precision floating-point number) but it can be changed by 
specifying `data_type` (e.g. `Complex{BigFloat}` for complex arbitrary precision number).

While the numerical GSN solution is only accurate in the range `[rsin, rsout]`, 
the full GSN solution is constructed by smoothly attaching the asymptotic solutions near horizon (up to `horizon_expansion_order`-th order) 
and infinity (up to `infinity_expansion_order`-th order). Therefore, the now-semi-analytical GSN solution is *accurate everywhere*.
"""
function GSN_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=3, infinity_expansion_order::Int=6,
    data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE
)
    # Compute the SWSH eigenvalue
    lambda = spin_weighted_spheroidal_eigenvalue(s, l, m, a*omega)
    # Fill in the mode information
    mode = Mode(s, l, m, a, omega, lambda)
    if boundary_condition == IN
        # Solve for Xin
        # NOTE For now we do *not* implement intelligent switching between the Riccati and the GSN form

        # Actually solve for Phiin first
        Phiinsoln = Solutions.solve_Phiin(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=horizon_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
        # Then convert to Xin
        Xinsoln = Solutions.Xsoln_from_Phisoln(Phiinsoln)

        # Extract the incidence and reflection amplitudes (NOTE: transmisson amplitude is *always* 1)
        Bref_SN, Binc_SN = Solutions.BrefBinc_SN_from_Xin(s, m, a, omega, lambda, Xinsoln, rsout; order=infinity_expansion_order)

        # Construct the full, 'semi-analytical' GSN solution
        semianalytical_Xinsoln(rs) = Solutions.semianalytical_Xin(s, m, a, omega, lambda, Xinsoln, rsin, rsout, horizon_expansion_order, infinity_expansion_order, rs)

        return GSNRadialFunction(
            mode,
            IN,
            rsin,
            rsout,
            horizon_expansion_order,
            infinity_expansion_order,
            data_type(1),
            Binc_SN,
            Bref_SN,
            nothing,
            Phiinsoln,
            semianalytical_Xinsoln,
            UNIT_GSN_TRANS
        )
    elseif boundary_condition == UP
        # Solve for Xup
        # NOTE For now we do *not* implement intelligent switching between the Riccati and the GSN form

        # Actually solve for Phiup first
        Phiupsoln = Solutions.solve_Phiup(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=infinity_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
        # Then convert to Xup
        Xupsoln = Solutions.Xsoln_from_Phisoln(Phiupsoln)

        # Extract the incidence and reflection amplitudes (NOTE: transmisson amplitude is *always* 1)
        Cref_SN, Cinc_SN = Solutions.CrefCinc_SN_from_Xup(s, m, a, omega, lambda, Xupsoln, rsin; order=horizon_expansion_order)

        # Construct the full, 'semi-analytical' GSN solution
        semianalytical_Xupsoln(rs) = Solutions.semianalytical_Xup(s, m, a, omega, lambda, Xupsoln, rsin, rsout, horizon_expansion_order, infinity_expansion_order, rs)

        return GSNRadialFunction(
            mode,
            UP,
            rsin,
            rsout,
            horizon_expansion_order,
            infinity_expansion_order,
            data_type(1),
            Cinc_SN,
            Cref_SN,
            nothing,
            Phiupsoln,
            semianalytical_Xupsoln,
            UNIT_GSN_TRANS
        )
    else
        error("Boundary condition must be IN or UP")
    end
end

# The power of multiple dispatch
(gsn_func::GSNRadialFunction)(rs) = gsn_func.GSN_solution(rs)[1] # Only return X(rs), discarding the first derivative

function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=3, infinity_expansion_order::Int=6,
    data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE
)
    # Solve for the GSN solution
    gsn_func = GSN_radial(s, l, m, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance)

    # Convert asymptotic amplitudes from GSN to Teukolsky formalism
    if gsn_func.boundary_condition == IN
        transmission_amplitude_conv_factor = ConversionFactors.Btrans(s, m, a, omega, gsn_func.mode.lambda)
        incidence_amplitude = ConversionFactors.Binc(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.incidence_amplitude / transmission_amplitude_conv_factor
        reflection_amplitude = ConversionFactors.Bref(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.reflection_amplitude / transmission_amplitude_conv_factor
    elseif gsn_func.boundary_condition == UP
        transmission_amplitude_conv_factor = ConversionFactors.Ctrans(s, m, a, omega, gsn_func.mode.lambda)
        incidence_amplitude = ConversionFactors.Cinc(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.incidence_amplitude / transmission_amplitude_conv_factor
        reflection_amplitude = ConversionFactors.Cref(s, m, a, omega, gsn_func.mode.lambda) * gsn_func.reflection_amplitude / transmission_amplitude_conv_factor
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
        gsn_func,
        teuk_func,
        UNIT_TEUKOLSKY_TRANS
    )
end

# The power of multiple dispatch
(teuk_func::TeukolskyRadialFunction)(r) = teuk_func.Teukolsky_solution(r)[1] # Only return R(r), discarding the first derivative

end
