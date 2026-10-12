using GeneralizedSasakiNakamura
using Test

@testset "GeneralizedSasakiNakamura.jl" begin
    @testset "Public radial interfaces" begin
        X = GSN_radial(-2, 2, 2, 0.68, 0.3, IN)
        R = Teukolsky_radial(-2, 2, 2, 0.68, 0.3, IN)
        Y = Y_radial(-2, 2, 2, 0.68, 0.3, IN)
        r = 10.0
        rstar = rstar_from_r(0.68, r)

        @test X.method == "ISEM"
        @test R.GSN_solution.method == "ISEM"
        @test isfinite(X(rstar))
        @test isfinite(R(r))
        @test isfinite(Y(r))
        @test r_from_rstar(0.68, rstar) ≈ r
        for radial in (GSN_radial, Teukolsky_radial)
            @test_throws MethodError radial(-2, 2, 2, 0.68, 0.3, IN; rhom=1.0)
            @test_throws ArgumentError radial(-2, 2, 2, 0.68, 0.3, IN; xm=-1.0)
        end
    end

    @testset "Default Teukolsky tolerance" begin
        # Independent complex-path and MST reference amplitudes.
        cases = (
            (omega=0.3 - 3.0im,
             reflection=0.000034049211003554543018 - 0.000017510919791042320907im,
             incidence=5.0661768763345224235 + 0.25720529286268218959im),
            (omega=0.37467168441804183 - 0.0889623156889357im,
             reflection=-0.79350698537260631 + 1.0219268994327333im,
             incidence=-0.06972551142642598 - 0.053077386482722241im),
        )
        for c in cases
            _, Rup = Teukolsky_radial(-2, 2, 2, 0.0, c.omega)
            @test Rup.reflection_amplitude / Rup.incidence_amplitude ≈
                c.reflection / c.incidence rtol=1e-10
            route = GeneralizedSasakiNakamura.ISEM.radial_route(
                -2, 2, 2, 0.0, c.omega, :UP)
            @test Rup.reflection_amplitude / Rup.transmission_amplitude ≈
                route.teukolsky_reflection / route.teukolsky_transmission rtol=1e-14
        end
    end

    @testset "Near-positive-imaginary-axis amplitudes" begin
        # Independent high-precision integration of the Teukolsky equation.
        for (boundary, incidence, reflection) in (
            (IN, 3846.386902042131 + 4.542548119603213e-6im,
             45.67541400622709 + 41.544821068135516im),
            (UP, 360.59877206644977 + 3.507391420322909e-7im,
             214.59326356735199 + 7.213446066030415e-7im),
        )
            R = Teukolsky_radial(-2, 2, 2, 0.0, 1e-10 + 0.3im, boundary)
            @test R.incidence_amplitude / R.transmission_amplitude ≈ incidence rtol=1e-10
            @test R.reflection_amplitude / R.transmission_amplitude ≈ reflection rtol=1e-10
        end
    end

    @testset "Negative-frequency amplitudes" begin
        # Independent Teukolsky references reflected by (m, omega) -> (-m, -conj(omega)).
        for (l, m, omega, incidence, reflection) in (
            (2, -2, -1.0 - 0.05im,
             1.2856521453922971 - 0.44862664931401232im,
             -2.9559869085975771e-4 - 1.4693045725576587e-3im),
            (3, -3, -1.0,
             0.5043157658448172 + 0.24588840385554175im,
             0.00012198157685068585 + 0.18138339086652683im),
        )
            R = Teukolsky_radial(-2, l, m, 0.9, omega, IN)
            @test R.incidence_amplitude / R.transmission_amplitude ≈ incidence rtol=1e-9
            @test R.reflection_amplitude / R.transmission_amplitude ≈ reflection rtol=1e-9
        end
    end

    @testset "Infinity-normalized radial value" begin
        # Independent 256-bit integration of the Teukolsky equation.
        R = Teukolsky_radial(-2, 2, 2, 0.5, 0.19660151220340238, IN)
        @test R(6.0) / R.incidence_amplitude ≈
            0.05226426135520719944 - 0.06822458072753078286im rtol=1e-12
    end

    @testset "Explicit junction failure" begin
        for radial in (GSN_radial, Teukolsky_radial)
            @test_throws ErrorException radial(
                -2, 2, 2, 0.7, 0.4198385710926635 - 2.09355961067422im, UP;
                xm=0.6)
        end
    end

    @testset "Reflected amplitude uncertainty" begin
        direct = GeneralizedSasakiNakamura.ISEM
        positive = direct.radial_route(-2, 4, 4, 0.7, 0.1 - 0.05im, :IN)
        negative = direct.radial_route(-2, 4, -4, 0.7, -0.1 - 0.05im, :IN)
        errors = direct.Matching.amplitude_errors(positive)
        @test all(isfinite, errors)
        @test errors == direct.Matching.amplitude_errors(negative)
    end

    @testset "Real-frequency MST branch" begin
        direct = GeneralizedSasakiNakamura.ISEM
        mst = direct.MSTInfinity
        parameters = direct.gsn_parameters(-2, 3, 3, 0.9, 1.0)
        nu = mst._mst_data(parameters).params.nu
        amplitudes = mst.mst_principal_amplitudes(parameters, :IN)
        @test abs(imag(nu)) ≈ 1.2619741165113905 atol=1e-10 rtol=0
        @test amplitudes.teuk[1] ≈
            0.5043157658448172 - 0.24588840385554175im rtol=1e-10
    end

    @testset "Exact-extremal radial interface" begin
        coordinates = GeneralizedSasakiNakamura.Coordinates
        extremal = GeneralizedSasakiNakamura.ExtremalRadial
        for a in (-1.0, 1.0)
            @test coordinates._is_exact_extremal_spin(a)
            @test extremal.is_exact_extremal_spin(a)
            @test GeneralizedSasakiNakamura.ISEM._is_exact_extremal_spin(a)
            @test !coordinates._is_exact_extremal_spin(sign(a) * prevfloat(1.0))
            @test !extremal.is_exact_extremal_spin(sign(a) * prevfloat(1.0))
            @test !GeneralizedSasakiNakamura.ISEM._is_exact_extremal_spin(
                sign(a) * prevfloat(1.0))
            @test_throws DomainError extremal._require_nonsynchronous(a, 2, a)
            @test extremal._require_nonsynchronous(a, 2, nextfloat(a)) ==
                2nextfloat(a) - 2a
            @test extremal._require_nonsynchronous(a, 0, 1e-13) == 2e-13
        end
        X = GSN_radial(
            -2, 2, 2, 1.0, 0.3, IN, -40.0, 60.0;
            method="ISEM", tolerance=1.0e-10)
        @test X.method == "ISEM"
        @test all(isfinite, X.GSN_solution(0.0))
        @test r_from_rstar(-1.0, rstar_from_r(-1.0, 2.0)) ≈ 2.0
        @test_throws DomainError GSN_radial(
            -2, 2, 2, 1.0, 1.0, IN; method="ISEM")
    end

    @testset "QNM interface" begin
        ordinary_mode = qnm(0.68, -2, 2, 2, 0)
        mirror_mode = qnm(0.68, -2, 2, 2, 0, mirror)

        @test ordinary_mode isa QNMResult
        @test mirror_mode isa QNMResult
        @test ordinary_mode.status == mirror_mode.status == :accepted
        @test real(ordinary_mode.omega) > 0
        @test real(mirror_mode.omega) < 0
        @test imag(ordinary_mode.omega) < 0
        @test imag(mirror_mode.omega) < 0
        @test ordinary_mode.excitation_factor ≈
            ordinary_mode.reflection_amplitude /
            (2 * ordinary_mode.omega * ordinary_mode.incidence_derivative)

        display_text = sprint(show, MIME("text/plain"), ordinary_mode)
        @test length(split(display_text, '\n')) == 9
        @test startswith(display_text, "QuasiNormalMode(\n")
        @test endswith(display_text, "    formalism = GSN)")

        detailed = qnm(
            0.68, -2, 2, 2, 0;
            primary_guess=ordinary_mode.omega, detailed=true)
        @test detailed.X isa GSNRadialFunction
        @test detailed.Y isa YRadialFunction
        @test detailed.R isa TeukolskyRadialFunction

        overtone7 = qnm_frequency(QNMMode(-2, 2, 2, 7), 0.0)
        overtone10 = qnm_frequency(QNMMode(-2, 2, 2, 10), 0.0)
        @test overtone7.status == overtone10.status == :accepted
        @test abs(overtone10.omega - overtone7.omega) > 0.1

        # Overtone label in the Python qnm / Cook convention.
        overtone29 = qnm_frequency(QNMMode(-2, 2, 2, 29), 0.5)
        @test overtone29.status == :accepted
        @test overtone29.omega ≈ 0.2725174928809329 - 6.734330683370764im atol=1e-10 rtol=0

        for (branch, omega8) in (
                (ordinary, 0.30101833878061446 - 1.6659917852234734im),
                (mirror, -0.09835632569909068 - 2.078506281482729im))
            overtone8 = qnm(0.5, -2, 2, 2, 8, branch)
            @test overtone8 isa QNMResult
            @test overtone8.omega ≈ omega8 atol=1e-9 rtol=0
        end

        special = qnm_frequency(QNMMode(-2, 2, 2, 8), 0.0)
        endpoint = qnm(1.0, -2, 2, 2, 0)
        @test special.omega == -2.0im
        @test endpoint isa QNMEndpointResult
        @test ismissing(endpoint.excitation_factor)
    end

    @testset "Schwarzschild flux rotation" begin
        equatorial = Teukolsky_pointparticle_flux(0.0, 10.0, 0.0, 1.0)
        for x in (-0.5, 0.5)
            tilted = Teukolsky_pointparticle_flux(0.0, 10.0, 0.0, x)
            @test tilted.infinity_energy_flux ≈ equatorial.infinity_energy_flux rtol=1e-8
            @test tilted.horizon_energy_flux ≈ equatorial.horizon_energy_flux rtol=1e-8
            @test tilted.infinity_angular_momentum_flux ≈
                x * equatorial.infinity_angular_momentum_flux rtol=1e-8
            @test tilted.horizon_angular_momentum_flux ≈
                x * equatorial.horizon_angular_momentum_flux rtol=1e-8
        end
    end

    @testset "Point-particle mode" begin
        mode = Teukolsky_pointparticle_mode(
            -2, 2, 2, 0, 0, 0.5, 8.0, 0.0, 1.0)
        @test mode.method.radial_method == "ISEM"
        @test isfinite(mode.amplitude)
        @test isfinite(mode.energy_flux)
        for constructor in (Teukolsky_pointparticle_mode, GSN_pointparticle_mode)
            absent = constructor(-2, 2, 2, 1, 0, 0.5, 8.0, 0.0, 1.0)
            @test iszero(absent.amplitude) && iszero(absent.energy_flux)
            @test absent.mode.omega === nothing
        end
        for e in (0.0, 0.1), x in (1.0, -1.0, 1.0)
            n = iszero(e) ? 0 : 1
            signed_mode = Teukolsky_pointparticle_mode(
                -2, 2, 2, n, 0, 0.5, 11.0, e, x; N=64, Nmax=128)
            orbit = GeneralizedSasakiNakamura.KerrGeodesics.kerr_geo_orbit(0.5, 11.0, e, x)
            frequencies = orbit["Frequencies"]
            @test signed_mode.mode.omega ≈
                (2frequencies["ϒϕ"] + n*frequencies["ϒr"]) / frequencies["ϒt"]
        end
        @test_throws ArgumentError Teukolsky_pointparticle_mode(
            -2, 2, 2, 1, 1, 0.0, 20.0, 1.1, 0.7)
    end
end
