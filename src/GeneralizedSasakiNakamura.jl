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

export GSN_radial, Teukolsky_radial # Homogeneous solutions
export GSN_pointparticle_mode, Teukolsky_pointparticle_mode # Inhomogeneous solutions

# Default values
_DEFAULT_rsin = -50
_DEFAULT_rsout = 1000
_DEFAULT_horizon_expansion_order = 3
_DEFAULT_infinity_expansion_order = 6

_DEFAULT_rhoin = -5000
_DEFAULT_rhoout = 5000
_DEFAULT_horizon_expansion_order_for_cplx_freq = 25
_DEFAULT_infinity_expansion_order_for_cplx_freq = 25

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
    print(io, "Mode(s=$(mode.s), l=$(mode.l), m=$(mode.m), a=$(mode.a), omega=$(mode.omega), lambda=$(mode.lambda))")
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

# Implement pretty-printing for GSNRadialFunction
# Mostly to suppress the printing of the numerical solution
function Base.show(io::IO, ::MIME"text/plain", gsn_func::GSNRadialFunction)
    println(io, "GSNRadialFunction(")
    print(io, "    mode="); show(io, "text/plain", gsn_func.mode); println(io, ",")
    println(io, "    boundary_condition=$(gsn_func.boundary_condition),")
    println(io, "    rsin=$(gsn_func.rsin),")
    println(io, "    rsout=$(gsn_func.rsout),")
    println(io, "    rsmp=$(gsn_func.rsmp),")
    println(io, "    horizon_expansion_order=$(gsn_func.horizon_expansion_order),")
    println(io, "    infinity_expansion_order=$(gsn_func.infinity_expansion_order),")
    println(io, "    transmission_amplitude=$(gsn_func.transmission_amplitude),")
    println(io, "    incidence_amplitude=$(gsn_func.incidence_amplitude),")
    println(io, "    reflection_amplitude=$(gsn_func.reflection_amplitude),")
    println(io, "    normalization_convention=$(gsn_func.normalization_convention)")
    println(io, "    method=$(gsn_func.method)")
    print(io, ")")
end

function Base.show(io::IO, gsn_func::GSNRadialFunction)
    print(io, "GSNRadialFunction(mode=Mode(s=$(gsn_func.mode.s), l=$(gsn_func.mode.l), m=$(gsn_func.mode.m), a=$(gsn_func.mode.a), omega=$(gsn_func.mode.omega), lambda=$(gsn_func.mode.lambda)), boundary_condition=$(gsn_func.boundary_condition))")
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

function Base.show(io::IO, teuk_func::TeukolskyRadialFunction)
    print(io, "TeukolskyRadialFunction(mode=Mode(s=$(teuk_func.mode.s), l=$(teuk_func.mode.l), m=$(teuk_func.mode.m), a=$(teuk_func.mode.a), omega=$(teuk_func.mode.omega), lambda=$(teuk_func.mode.lambda)), boundary_condition=$(teuk_func.boundary_condition))")
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
with tolerance specified by `tolerance` (default: `1e-12`). The ODE to be solved is determined by the keyword `method` (default: `auto`),
which can be either `linear` (solving the GSN equation in a linear form) or `Riccati` (solving the GSN equation in a nonlinear form).
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
    method="auto", data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, rsmp=nothing
)
    if omega == 0
        return GSN_radial(s, l, m, a, omega, boundary_condition)
    else
        sawwarn = Ref(false)
        logger = EarlyFilteredLogger(
            (log)->(log.level == Logging.Warn && (sawwarn[] = true); log.level ≥ Logging.Error),
            current_logger())
        # Compute the SWSH eigenvalue
        lambda = spin_weighted_spheroidal_eigenvalue(s, l, m, a*omega)
        # Fill in the mode information
        mode = Mode(s, l, m, a, omega, lambda)
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
                    error("Method must be either 'Riccati' or 'linear'")
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
                    error("Method must be either 'Riccati' or 'linear'")
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
                        @warn "There might be an issue with solving the Riccati form of the GSN equation. Consider using method=\"linear\" instead."
                    end
                    # Then convert to Xup
                    Xupsoln = Solutions.Xsoln_from_Phisoln(Phiupsoln)
                elseif method == "linear"
                    Phiupsoln = nothing
                    Xupsoln = Solutions.solve_Xup(s, m, a, omega, lambda, rsin, rsout; initialconditions_order=infinity_expansion_order, dtype=data_type, odealgo=ODE_algorithm, abstol=tolerance, reltol=tolerance)
                else
                    error("Method must be either 'Riccati' or 'linear'")
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
                    error("Method must be either 'Riccati' or 'linear'")
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
    s::Int, l::Int, m::Int, a, omega, boundary_condition
)
    if omega != 0
        return GSN_radial(s, l, m, a, omega, boundary_condition, _DEFAULT_rsin, _DEFAULT_rsout)
    else
        teuk_func = Teukolsky_radial(s, l, m, a, omega, boundary_condition)
        GSN_solution = Solutions.Sasaki_Nakamura_function_from_Teukolsky_radial_function(s, m, a, omega, teuk_func.mode.lambda, teuk_func.Teukolsky_solution)

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
                method
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
function GSN_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, method="auto")
    if omega == 0
        Xin = GSN_radial(s, l, m, a, omega, IN)
        Xup = GSN_radial(s, l, m, a, omega, UP)
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

Note, however, when `omega = 0`, the exact Teukolsky function expressed using Gauss hypergeometric functions will be returned (i.e., instead of using the GSN formalism). 
In this case, only `s`, `l`, `m`, `a`, `omega`, `boundary_condition` will be parsed.

With complex values of `omega`, the Teukolsky function is evaluated as a function of $r$,
where the value at the corresponding location $\rho = r_{*}(r) \in \mathbb{R}$ along the rotated integration path on the complex plane is returned.
"""
function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition, rsin, rsout;
    horizon_expansion_order::Int=_DEFAULT_horizon_expansion_order, infinity_expansion_order::Int=_DEFAULT_infinity_expansion_order,
    method="auto", data_type=Solutions._DEFAULTDATATYPE,  ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, rsmp=nothing
)
    if omega == 0
        return Teukolsky_radial(s, l, m, a, 0, boundary_condition)
    else

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
            gsn_func,
            teuk_func,
            UNIT_TEUKOLSKY_TRANS
        )
    end
end

function Teukolsky_radial(
    s::Int, l::Int, m::Int, a, omega, boundary_condition
)
    if omega != 0
        return Teukolsky_radial(s, l, m, a, omega, boundary_condition, _DEFAULT_rsin, _DEFAULT_rsout)
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

@doc raw"""
    Teukolsky_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, method="auto")

Compute the Teukolsky function for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `a` the Kerr spin parameter, and `omega` the frequency)
with the purely-ingoing boundary condition at the horizon (`IN`) and the purely-outgoing boundary condition at infinity (`UP`).

Note that for _real_ frequencies, the numerical inner boundary (rsin) and outer boundary (rsout) are set to the default values `_DEFAULT_rsin` and `_DEFAULT_rsout`, respectively,
while the order of the asymptotic expansion at the horizon and infinity are determined automatically.
As for _complex_ frequencies, the numerical inner and the outer boundaries are determined automatically,
while the order of the asymptotic expansion at the horizon and infinity are set to `_DEFAULT_horizon_expansion_order_for_cplx_freq` and `_DEFAULT_infinity_expansion_order_for_cplx_freq`, respectively.
"""
function Teukolsky_radial(s::Int, l::Int, m::Int, a, omega; data_type=Solutions._DEFAULTDATATYPE, ODE_algorithm=Solutions._DEFAULTSOLVER, tolerance=Solutions._DEFAULTTOLERANCE, method="auto")
    if omega == 0
        Rin = Teukolsky_radial(s, l, m, a, omega, IN)
        Rup = Teukolsky_radial(s, l, m, a, omega, UP)
    else
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
    end

    return (Rin, Rup)
end

# The power of multiple dispatch
(teuk_func::TeukolskyRadialFunction)(r) = teuk_func.Teukolsky_solution(r)[1] # Only return R(r), discarding the first derivative


include("Inhomogeneous/AsymptoticExpansionCoefficientsY.jl")
include("Inhomogeneous/InitialConditionsY.jl")
include("Inhomogeneous/SolutionsY.jl")
include("Inhomogeneous/GridSampling.jl")
include("Inhomogeneous/ConvolutionIntegrals.jl")


function Base.show(io::IO, ::MIME"text/plain", ysol::SolutionsY.YSolutionResult)
    println(io, "YSolution(")
    print(io, "    basis="); show(io, ysol.basis_type); println(io, ",")
    print(io, "    mode="); show(io, "text/plain", ysol.mode); println(io, ",")
    println(io, "    rin=$(ysol.asymptotic.rin),")
    println(io, "    rout=$(ysol.asymptotic.rout),")
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
    print(io, "Mode(s=$(mode.s), l=$(mode.l), m=$(mode.m), n=$(mode.n), k=$(mode.k), a=$(mode.a), omega=$(mode.omega), lambda=$(mode.lambda))")
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
    print(io, "    mode="); show(io, "text/plain", teuk_mode.mode); println(io, ",")
    println(io, "    amplitude$suffix=$(teuk_mode.amplitude),")
    println(io, "    energy_flux$suffix=$(teuk_mode.energy_flux),")
    println(io, "    angular_momentum_flux$suffix=$(teuk_mode.angular_momentum_flux),")
    println(io, "    Carter_const_flux$suffix=$(teuk_mode.Carter_const_flux),")
    println(io, "    method=$(teuk_mode.method),")
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
    print(io, "    mode="); show(io, "text/plain", gsn_mode.mode); println(io, ",")
    println(io, "    amplitude$suffix=$(gsn_mode.amplitude),")
    println(io, "    energy_flux$suffix=$(gsn_mode.energy_flux),")
    println(io, "    angular_momentum_flux$suffix=$(gsn_mode.angular_momentum_flux),")
    println(io, "    Carter_const_flux$suffix=$(gsn_mode.Carter_const_flux),")
    println(io, "    method=$(gsn_mode.method),")
    print(io, ")")
end

@doc raw"""
    Teukolsky_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int, K::Int)

Compute the amplitude of the inhomogeneous Teukolsky solution _at infinity_ (for s = - 2) and _at horizon_ (for s = + 2) due to a point particle on a generic timelike bound orbit around a Kerr black hole,
with a spin parameter of `a`, for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `n` the radial index and `k` the polar index).
The orbit is specified by `p` the semi-latus rectum, `e` the eccentricity and `x` the inclination parameter ($x \equiv \cos \theta_{\mathrm{inc}}$).

In addition, we compute also the energy, angular momentum and Carter constant flux at infinity (for s = - 2) and the horizon (for s = + 2). Note that these values are formalism-independent.

The numerical method to compute the convolution integral is specified by `method` (default: `auto`), which can either be `trapezoidal` or `levin`.
We sample the trajectory over a grid of size N x K, where N and K are the number of Chebyshev nodes in the radial and the polar direction, respectively.
Note that they must be powers of 2.
"""
function Teukolsky_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int=-1, K::Int=-1)
    if method == "auto"
        # For now, choose "trapezoidal"
        method = "trapezoidal"
    end

    # Fill-in default values for N/K if none is given
    if method == "trapezoidal"
        N = N <= 0 ? 256 : N
        K = K <= 0 ? 64 : K
    elseif method == "levin"
        N = N <= 0 ? 256 : N
        K = K <= 0 ? 32 : K
    end

    if method == "trapezoidal"
        output = ConvolutionIntegrals.convolution_integral_trapezoidal(a, p, e, x, s, l, m, n, k; N=N, K=K)
    elseif method == "levin"
        output = ConvolutionIntegrals.convolution_integral_levin(a, p, e, x, s, l, m, n, k; N=N, K=K)
    else
        error("Currently only support method = \"trapezoidal\" or \"levin\"")
    end

    return TeukolskyPointParticleMode(
            PointParticleMode(s, l, m, n, k, a, output["omega"], output["YSolution"].mode.lambda),
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

@doc raw"""
    GSN_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int, K::Int)

Compute the amplitude of the inhomogeneous GSN solution _at infinity_ (for s = - 2) and _at horizon_ (for s = + 2) due to a point particle on a generic timelike bound orbit around a Kerr black hole,
with a spin parameter of `a`, for a given mode (specified by `s` the spin weight, `l` the harmonic index, `m` the azimuthal index, `n` the radial index and `k` the polar index).
The orbit is specified by `p` the semi-latus rectum, `e` the eccentricity and `x` the inclination parameter ($x \equiv \cos \theta_{\mathrm{inc}}$).

In addition, we compute also the energy, angular momentum and Carter constant flux at infinity (for s = - 2) and the horizon (for s = + 2). Note that these values are formalism-independent.

The numerical method to compute the convolution integral is specified by `method` (default: `auto`), which can either be `trapezoidal` or `levin`.
We sample the trajectory over a grid of size N x K, where N and K are the number of Chebyshev nodes in the radial and the polar direction, respectively.
Note that they must be powers of 2.
"""
function GSN_pointparticle_mode(s::Int, l::Int, m::Int, n::Int, k::Int, a, p, e, x; method="auto", N::Int=-1, K::Int=-1)
    Teukolsky_mode = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x; method=method, N=N, K=K)
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