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

# Default values
_DEFAULT_rsin = -50
_DEFAULT_rsout = 1000
_DEFAULT_horizon_expansion_order = 3
_DEFAULT_infinity_expansion_order = 6

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
    numerical_GSN_solution::Union{ODESolution, Missing} # Store the numerical solution to the GSN equation in [rsin, rsout]
    numerical_Riccati_solution::Union{ODESolution, Missing} # Store the numerical solution to the GSN equation in the Riccati form if applicable
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
    GSN_solution::Union{GSNRadialFunction, Missing} # Store the full GSN solution
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
    GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order, data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE)

Compute the GSN function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency) 
and boundary condition specified by `boundary_condition`, which can be either

    - `IN` for purely-ingoing at the horizon,
    - `UP` for purely-outgoing at infinity,
    - `OUT` for purely-outgoing at the horizon,
    - `DOWN` for purely-ingoing at infinity.

Note that the `OUT` and `DOWN` solutions are constructed by linearly combining the `IN` and `UP` solutions, respectively.

The GSN function is numerically solved in the interval of *tortoise coordinates* $r_{*} \in$ `[rsin, rsout]` using the ODE solver (from `DifferentialEquations.jl`) specified by `ODE_algorithm` (default: `Vern9()`) 
with tolerance specified by `tolerance` (default: `1e-12`). By default the data type used is `ComplexF64` (i.e. double-precision floating-point number) but it can be changed by 
specifying `data_type` (e.g. `Complex{BigFloat}` for complex arbitrary precision number).

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
    data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE
)
    if omega == 0
        return GSN_radial(s, l, m, a, omega, boundary_condition)
    else
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
                missing,
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
                missing,
                Phiupsoln,
                semianalytical_Xupsoln,
                UNIT_GSN_TRANS
            )
        elseif boundary_condition == DOWN
            # Construct Xdown from Xin and Xup, instead of solving the ODE numerically with the boundary condition

            # Solve for Xin *and* Xup first
            Xin = GSN_radial(s, l, m, a, omega, IN, rsin, rsout, horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance)
            Xup = GSN_radial(s, l, m, a, omega, UP, rsin, rsout, horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance)

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
                horizon_expansion_order,
                infinity_expansion_order,
                data_type(1),
                Btrans/Binc - (Bref*Cref)/(Binc*Ctrans),
                (-Bref*Cinc)/(Binc*Ctrans),
                missing,
                missing,
                _full_Xdown_solution,
                UNIT_GSN_TRANS
            )
        elseif boundary_condition == OUT
            # Construct Xout from Xin and Xup, instead of solving the ODE numerically with the boundary condition

            # Solve for Xin *and* Xup first
            Xin = GSN_radial(s, l, m, a, omega, IN, rsin, rsout, horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance)
            Xup = GSN_radial(s, l, m, a, omega, UP, rsin, rsout, horizon_expansion_order=horizon_expansion_order, infinity_expansion_order=infinity_expansion_order, data_type=data_type, ODE_algorithm=ODE_algorithm, tolerance=tolerance)

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
                horizon_expansion_order,
                infinity_expansion_order,
                data_type(1),
                Ctrans/Cinc - (Bref*Cref)/(Btrans*Cinc),
                (-Binc*Cref)/(Btrans*Cinc),
                missing,
                missing,
                _full_Xout_solution,
                UNIT_GSN_TRANS
            )
        else
            error("Boundary condition must be IN or UP")
        end
    end
end

@doc raw"""
    GSN_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition)

Compute the exact *static* (`omega = 0`) GSN function using Gauss hypergeometric functions 
for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency) 
and boundary condition (specified by `boundary_condition` which can be either `IN` for purely-ingoing at the horizon or `UP` for purely-outgoing at infinity).

Note that the GSN function is transformed from the exact Teukolsky function using the GSN transformation.
"""
function GSN_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition
)
    if omega != 0
        error("Cannot compute the GSN function for a nonstatic (omega != 0) case without specifying rsin and rsout")
    else
        teuk_func = Teukolsky_radial(s, l, m, a, omega, boundary_condition)
        GSN_solution = Solutions.Sasaki_Nakamura_function_from_Teukolsky_radial_function(s, m, a, omega, teuk_func.mode.lambda, teuk_func.Teukolsky_solution)

        return GSNRadialFunction(
                teuk_func.mode,
                boundary_condition,
                missing,
                missing,
                -1,
                -1,
                missing,
                missing,
                missing,
                missing,
                missing,
                GSN_solution,
                UNIT_TEUKOLSKY_TRANS
            )
    end
end

# The power of multiple dispatch
(gsn_func::GSNRadialFunction)(rs) = gsn_func.GSN_solution(rs)[1] # Only return X(rs), discarding the first derivative

@doc raw"""
    Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout; horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order, data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE)

Compute the Teukolsky function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency) 
and boundary condition specified by `boundary_condition`, which can be either

    - `IN` for purely-ingoing at the horizon,
    - `UP` for purely-outgoing at infinity,
    - `OUT` for purely-outgoing at the horizon,
    - `DOWN` for purely-ingoing at infinity.

Note that the `OUT` and `DOWN` solutions are constructed by linearly combining the `IN` and `UP` solutions, respectively.

The full GSN solution is converted to the corresponding Teukolsky solution $(R(r), dR/dr)$ and 
the incidence, reflection and transmission amplitude are converted from the GSN formalism to the Teukolsky formalism 
with the normalization convention that the transmission amplitude is normalized to 1 (i.e. `normalization_convention=UNIT_TEUKOLSKY_TRANS`).

Note, however, when `omega = 0`, the exact Teukolsky function expressed using Gauss hypergeometric functions will be returned (i.e., instead of using the GSN formalism). 
In this case, only `s`, `l`, `m`, `a`, `omega`, `boundary_condition` will be parsed.
"""
function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order,
    data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE
)
    if omega == 0
        return Teukolsky_radial(s, l, m, a, 0, boundary_condition)
    else

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
            gsn_func,
            teuk_func,
            UNIT_TEUKOLSKY_TRANS
        )
    end
end

@doc raw"""
    Teukolsky_radial(s::Int, l::Int, m::Int, a, omega, boundary_condition)

Compute the exact *static* (`omega = 0`) Teukolsky function using Gauss hypergeometric functions 
for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency) 
and boundary condition (specified by `boundary_condition` which can be either `IN` for purely-ingoing at the horizon or `UP` for purely-outgoing at infinity).
"""
function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition
)
    if omega != 0
        error("Cannot compute the Teukolsky function for a nonstatic (omega != 0) case without specifying rsin and rsout")
    else
        # Compute the SWSH eigenvalue
        lambda = spin_weighted_spherical_eigenvalue(s, l, m)
        # Fill in the mode information
        mode = Mode(s, l, m, a, omega, lambda)
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
            teuk_func,
            UNIT_TEUKOLSKY_TRANS
        )
    end
end

# The power of multiple dispatch
(teuk_func::TeukolskyRadialFunction)(r) = teuk_func.Teukolsky_solution(r)[1] # Only return R(r), discarding the first derivative

end
