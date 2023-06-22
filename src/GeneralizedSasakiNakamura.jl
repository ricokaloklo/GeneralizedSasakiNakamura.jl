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

export GSN_radial

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

struct GSN_radial_function
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

struct Teukolsky_radial_function
    mode::Mode # Information about the mode
    boundary_condition::BoundaryCondition # The boundary condition that this radial function statisfies
    transmission_amplitude # In Teukolsky formalism
    incidence_amplitude # In Teukolsky formalism
    reflection_amplitude # In Teukolsky formalism
    GSN_solution::GSN_radial_function # Store the full GSN solution
    Teukolsky_solution # Store the full Teukolsky solution
    normalization_convention::NormalizationConvention # The normalization convention used for the *stored* Teukolsky solution
end

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

        return GSN_radial_function(
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

        return GSN_radial_function(
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

end
