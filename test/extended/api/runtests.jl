# Public API: exports, enums, struct layouts, show methods, argument validation of the radial interfaces.
# Field lists are copied from src/GeneralizedSasakiNakamura.jl, src/ISEM/ISEM.jl and src/QNM/Types.jl;
# a change to them is an API change and must be deliberate.

if !@isdefined(GSN)
    include(joinpath(@__DIR__, "..", "load_package.jl"))
end
include(joinpath(@__DIR__, "..", "common", "helpers.jl"))

@testset "api" begin

@testset "exports" begin
    exported = (
        :r_from_rstar, :rstar_from_r, :GSN_radial, :Teukolsky_radial, :Y_radial,
        :GSN_pointparticle_mode, :Teukolsky_pointparticle_mode, :Teukolsky_pointparticle_flux,
        :BoundaryCondition, :NormalizationConvention, :Mode, :GSNRadialFunction, :TeukolskyRadialFunction,
        :YRadialFunction, :ISEM, :IN, :UP, :OUT, :DOWN, :UNIT_GSN_TRANS, :UNIT_TEUKOLSKY_TRANS,
        :QNMMode, :QNMBranch, :QNMResult, :QNMEstimate, :QNMFailure, :QNMEndpointResult, :QNMPairResult,
        :ordinary, :mirror, :Ordinary, :Mirror, :LeaverResult, :ISEMValidationResult, :ExcitationFactorResult,
        :qnm_frequency, :qnm_sequence, :validate_qnm_with_isem, :qnm_excitation_factor, :qnm, :qnm_pair,
        :angular_A_to_lambda, :lambda_to_angular_A,
    )
    public_names = names(GSN)
    for name in exported
        @test name in public_names
        @test isdefined(GSN, name)
    end
    # Internal types stay unexported.
    for name in (:PointParticleMode, :TeukolskyPointParticleMode, :GSNPointParticleMode, :TeukolskyPointParticleFlux)
        @test isdefined(GSN, name)
        @test !(name in public_names)
    end
    @test GSN.QNM isa Module && GSN.ISEM isa Module && GSN.ISEM.DirectGSN isa Module
end

@testset "enums" begin
    @test Int(IN) == 1 && Int(UP) == 2 && Int(OUT) == 3 && Int(DOWN) == 4
    @test IN isa BoundaryCondition
    @test Int(UNIT_GSN_TRANS) == 1 && Int(UNIT_TEUKOLSKY_TRANS) == 2
    @test UNIT_GSN_TRANS isa NormalizationConvention
    @test ordinary isa QNMBranch && mirror isa QNMBranch
    @test Ordinary === ordinary && Mirror === mirror
end

@testset "struct layouts" begin
    @test fieldnames(Mode) == (:s, :l, :m, :a, :omega, :lambda)
    @test fieldnames(GSNRadialFunction) == (:mode, :boundary_condition, :rsin, :rsout, :rsmp,
        :horizon_expansion_order, :infinity_expansion_order, :transmission_amplitude, :incidence_amplitude,
        :reflection_amplitude, :numerical_GSN_solution, :numerical_Riccati_solution, :GSN_solution,
        :normalization_convention, :method)
    @test fieldnames(TeukolskyRadialFunction) == (:mode, :boundary_condition, :transmission_amplitude,
        :incidence_amplitude, :reflection_amplitude, :P_solution, :GSN_solution, :Teukolsky_solution,
        :normalization_convention)
    @test fieldnames(YRadialFunction) == (:mode, :boundary_condition, :transmission_amplitude,
        :incidence_amplitude, :reflection_amplitude, :P_solution, :Teukolsky_solution, :X_solution,
        :Y_scalar_solution, :Y_solution, :normalization_convention)
    @test fieldnames(GSN.PointParticleMode) == (:s, :l, :m, :n, :k, :a, :omega, :lambda)
    for T in (GSN.TeukolskyPointParticleMode, GSN.GSNPointParticleMode)
        @test fieldnames(T) == (:mode, :amplitude, :energy_flux, :angular_momentum_flux, :Carter_const_flux,
            :trajectory, :Y_solution, :SWSH, :method)
    end
    @test fieldnames(GSN.TeukolskyPointParticleFlux) == (:a, :p, :e, :x, :orbit_type,
        :infinity_energy_flux, :infinity_angular_momentum_flux, :infinity_carter_constant_flux,
        :horizon_energy_flux, :horizon_angular_momentum_flux, :horizon_carter_constant_flux,
        :total_modes, :tolerance, :truncation_floor, :reached, :cost, :result)
    @test fieldnames(QNMMode) == (:s, :l, :m, :n, :branch)
end

@testset "radial return types and show" begin
    X = GSN_radial(-2, 2, 2, 0.68, 0.3, IN)
    R = Teukolsky_radial(-2, 2, 2, 0.68, 0.3, IN)
    Y = Y_radial(-2, 2, 2, 0.68, 0.3, IN)
    @test X isa GSNRadialFunction && R isa TeukolskyRadialFunction && Y isa YRadialFunction
    @test X.mode isa Mode && X.mode.s == -2 && X.mode.omega == 0.3
    @test X.transmission_amplitude isa Complex
    @test isfinite(X(0.0)) && isfinite(R(10.0)) && isfinite(Y(10.0))
    @test length(X.GSN_solution(0.0)) >= 2
    @test length(R.Teukolsky_solution(10.0)) >= 2
    pair = GSN_radial(-2, 2, 2, 0.68, 0.3)
    @test pair isa Tuple{GSNRadialFunction, GSNRadialFunction}
    @test pair[1].boundary_condition == IN && pair[2].boundary_condition == UP
    tpair = Teukolsky_radial(-2, 2, 2, 0.68, 0.3)
    @test tpair isa Tuple{TeukolskyRadialFunction, TeukolskyRadialFunction}
    @test tpair[1].incidence_amplitude ≈ R.incidence_amplitude rtol=1e-14
    # Negative spin through the pair interface maps to (-a, -m).
    npair = Teukolsky_radial(-2, 2, -2, -0.68, 0.3)
    @test npair[1].mode.a == -0.68 && npair[1].mode.m == -2
    @test npair[1].incidence_amplitude == tpair[1].incidence_amplitude
    # show
    @test startswith(sprint(show, MIME("text/plain"), X.mode), "Mode(s = -2, l = 2, m = 2")
    textX = sprint(show, MIME("text/plain"), X)
    @test startswith(textX, "GSNRadialFunction(\n") && endswith(textX, ")")
    @test occursin("method = GSN-ISEM", textX)
    @test startswith(sprint(show, X), "GSNRadialFunction(mode = Mode(s = -2")
    textR = sprint(show, MIME("text/plain"), R)
    @test startswith(textR, "TeukolskyRadialFunction(\n")
    @test occursin("normalization_convention = UNIT_TEUKOLSKY_TRANS", textR)
    @test startswith(sprint(show, R), "TeukolskyRadialFunction(mode = Mode(")
    @test startswith(sprint(show, MIME("text/plain"), pair), "(\n")
    @test startswith(sprint(show, MIME("text/plain"), tpair), "(\n")
    @test startswith(sprint(show, MIME("text/plain"), Y), "YRadialFunction(\n")
    @test startswith(sprint(show, Y), "YRadialFunction(mode = Mode(")
end

@testset "radial argument validation" begin
    for radial in (GSN_radial, Teukolsky_radial)
        @test_throws ArgumentError radial(-2, 2, 2, 0.68, 0.3, IN; rhom=1.0)
        @test_throws ArgumentError radial(-2, 2, 2, 0.68, 0.3, IN; xm=-1.0)
        @test_throws ArgumentError radial(-2, 2, 2, 0.68, 0.3, IN; xm=1.0)
        @test_throws ArgumentError radial(-2, 2, 2, 0.68, 0.3, IN; TSinHor=true)
        @test_throws DomainError radial(-2, 2, 2, 1.0, 1.0, IN; method="GSN-ISEM")
    end
    @test_throws ArgumentError GSN_radial(-2, 2, 2, 1.0, 0.3, IN; method="bogus")   # exact-extremal method label
    @test_throws ErrorException quietly(() -> GSN_radial(-2, 2, 2, 0.5, 0.3, IN; method="bogus"))
    # Y_radial supports s = -2 with IN and s = +2 with UP only.
    @test_throws ErrorException Y_radial(-2, 2, 2, 0.68, 0.3, UP)
    @test_throws ErrorException Y_radial(2, 2, 2, 0.68, 0.3, IN)
    @test_throws ErrorException Y_radial(-2, 2, 2, 0.68, 0.3, IN; method="bogus")
    Yup = Y_radial(2, 2, 2, 0.68, 0.3, UP)
    @test isfinite(Yup(10.0))
    # tol and tolerance are aliases on the GSN-ISEM route.
    a1 = GSN_radial(-2, 2, 2, 0.5, 0.4, IN; tol=1e-10)
    a2 = GSN_radial(-2, 2, 2, 0.5, 0.4, IN; tolerance=1e-10)
    @test a1.incidence_amplitude == a2.incidence_amplitude
end

@testset "point-particle return types and show" begin
    mode = Teukolsky_pointparticle_mode(-2, 2, 2, 0, 0, 0.5, 8.0, 0.0, 1.0)
    @test mode isa GSN.TeukolskyPointParticleMode
    @test mode.mode isa GSN.PointParticleMode
    @test mode.method.radial_method == "GSN-ISEM"
    @test mode.method.method == "isem_trapezoidal"
    @test haskey(mode.method, :truncation_floor)
    text = sprint(show, MIME("text/plain"), mode)
    @test startswith(text, "TeukolskyPointParticleMode(\n")
    @test occursin("energy_flux_inf", text)
    @test startswith(sprint(show, MIME("text/plain"), mode.mode), "Mode(s = -2, l = 2, m = 2, n = 0, k = 0")
    hor = Teukolsky_pointparticle_mode(2, 2, 2, 0, 0, 0.5, 8.0, 0.0, 1.0)
    @test occursin("energy_flux_hor", sprint(show, MIME("text/plain"), hor))
    gsn = GSN_pointparticle_mode(-2, 2, 2, 0, 0, 0.5, 8.0, 0.0, 1.0)
    @test gsn isa GSN.GSNPointParticleMode
    @test startswith(sprint(show, MIME("text/plain"), gsn), "GSNPointParticleMode(\n")
    F = Teukolsky_pointparticle_flux(0.5, 8.0, 0.0, 1.0; tol=1e-6)
    @test F isa GSN.TeukolskyPointParticleFlux
    @test F.orbit_type == :circular
    @test F.truncation_floor isa NamedTuple && haskey(F.truncation_floor, :infinity)
    @test F.cost >= 0
    @test occursin("l_reached", sprint(show, MIME("text/plain"), F))
end

@testset "QNM result types and show" begin
    failure = qnm(0.0, -2, 2, 2, 8)
    @test failure isa QNMFailure
    @test failure.status == :failed
    @test startswith(sprint(show, failure), "QNMFailure(mode=")
    @test startswith(sprint(show, MIME("text/plain"), failure), "QNMFailure(\n")
    @test :stage in propertynames(failure)
    endpoint = qnm(1.0, -2, 2, 2, 0)
    @test endpoint isa QNMEndpointResult
    @test startswith(sprint(show, endpoint), "QNMEndpointResult(mode=")
    root = qnm_frequency(QNMMode(-2, 2, 2, 0), 0.3)
    @test root isa LeaverResult
    @test root.convention == :overtone
    @test root.overtone_index == 0
    @test root.status == :accepted
    @test root.precision_bits >= 53
end

end # api
