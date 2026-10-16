# Single point-particle modes: Teukolsky_pointparticle_mode / GSN_pointparticle_mode.
# References: data/mode_references.jl and the 40-digit circular references copied below (provenance in
# ../README.md, section "Point-particle modes").

if !@isdefined(GSN)
    include(joinpath(@__DIR__, "..", "load_package.jl"))
end
include(joinpath(@__DIR__, "..", "common", "helpers.jl"))
include(data_file("mode_references.jl"))

const CF = GSN.ConversionFactors
const MS = GSN.ModeSummation
const KG = GSN.KerrGeodesics

mode_inf(l, m, n, k, a, p, e, x; kw...) = Teukolsky_pointparticle_mode(-2, l, m, n, k, a, p, e, x; kw...)
mode_hor(l, m, n, k, a, p, e, x; kw...) = Teukolsky_pointparticle_mode(2, l, m, n, k, a, p, e, x; kw...)

@testset "pointparticle_mode" begin

# ---------------------------------------------------------------------------------------------------
# Schwarzschild circular orbit, (l, m) = (2, +-2): package-independent references from three methods
# (Zerilli-Moncrief, MST, direct Teukolsky integration; 40-50 agreeing digits) at p = 10; direct Teukolsky
# integration only (50 digits) at p = 1000 and 10000. p = 10 measured -4.3e-12 (auto); p = 1000 and 10000
# probe the near-static regime (omega = 6.3e-5 and 6.3e-7) and use looser tolerances.
# ---------------------------------------------------------------------------------------------------
CIRCULAR_A0 = (
    (p=10.0, Einf=2.684397739551056800957597197786816721611e-5, Ehor=5.654138734536942973248997741647276606086e-9,
     rtol_inf=1e-10, rtol_hor=1e-9, level=:quick),
    (p=1000.0, Einf=3.1849724968834935461165736672274896085871713296801e-15,
     Ehor=3.209661437135926463982881220433809266204747992753e-27, rtol_inf=1e-7, rtol_hor=1e-6, level=:full),
    (p=10000.0, Einf=3.1984098320174248359757748087428211703184129871002e-20,
     Ehor=3.200960607383136877018940040561288621664276635493e-36, rtol_inf=1e-5, rtol_hor=1e-4, level=:full),
)

@testset "Schwarzschild circular (2,2): independent references" begin
    for c in CIRCULAR_A0
        at_level(c.level) || continue
        for m in (2, -2)
            @testset "p=$(c.p) m=$m" begin
                inf = mode_inf(2, m, 0, 0, 0.0, c.p, 0.0, 1.0)
                hor = mode_hor(2, m, 0, 0, 0.0, c.p, 0.0, 1.0)
                @test relerr(inf.energy_flux, c.Einf) <= c.rtol_inf
                @test relerr(hor.energy_flux, c.Ehor) <= c.rtol_hor
                # Circular: dL/dt = dE/dt / Omega_phi with Omega_phi = p^(-3/2) (prograde), Carter flux 0.
                @test inf.angular_momentum_flux ≈ inf.energy_flux * c.p^1.5 rtol=1e-10
                @test hor.angular_momentum_flux ≈ hor.energy_flux * c.p^1.5 rtol=1e-10
                @test inf.Carter_const_flux == 0 || abs(inf.Carter_const_flux) <= 1e-12 * inf.angular_momentum_flux
                @test inf.mode.omega ≈ m * c.p^-1.5 rtol=1e-12
                @test inf.method.radial_method == "ISEM"
            end
        end
    end
    # All convolution methods agree at p = 10 (measured: auto -4.3e-12, legacy trapezoidal +2.8e-11).
    for method in ("isem_trapezoidal", "isem_levin", "trapezoidal", "levin")
        mode = quietly(() -> mode_inf(2, 2, 0, 0, 0.0, 10.0, 0.0, 1.0; method=method))
        @test mode.method.method == method
        @test relerr(mode.energy_flux, CIRCULAR_A0[1].Einf) <= 1e-9
    end
end

# ---------------------------------------------------------------------------------------------------
# 56 spherical-orbit modes (7 orbits x 8 modes, n = 0). Infinity: vs published GSN 0.9.0 (<= 4e-13 from
# independent mpmath values) with rtol 3e-10 (frozen self-values differ by up to 1.37e-10). Horizon: vs the
# frozen self-values (rtol 1e-9) and published (rtol 1e-7; published horizon <= 6e-8 off). The a = 0
# (2,2,k=1) mode vanishes by rotation symmetry (both return ~1e-34 noise): tested as zero.
# ---------------------------------------------------------------------------------------------------
@testset "56 spherical modes" begin
    rows = level_subset(MODES_56, 7)
    for r in rows
        @testset "a=$(r.a) p=$(r.p) x=$(r.x) (l,m,k)=($(r.l),$(r.m),$(r.k))" begin
            inf = mode_inf(r.l, r.m, 0, r.k, r.a, r.p, 0.0, r.x)
            hor = mode_hor(r.l, r.m, 0, r.k, r.a, r.p, 0.0, r.x)
            @test inf.mode.omega ≈ r.omega rtol=1e-12
            @test hor.mode.omega ≈ r.omega rtol=1e-12
            vanishing = abs(r.pub_Einf_pp) < 1e-25
            if vanishing
                @test abs(inf.energy_flux) <= 1e-25
                @test abs(hor.energy_flux) <= 1e-25
            else
                @test relerr(inf.energy_flux, r.pub_Einf_pp) <= 3e-10
                @test relerr(hor.energy_flux, r.dev_Ehor_pp) <= 1e-9
                @test relerr(hor.energy_flux, r.pub_Ehor_pp) <= 1e-7
            end
            if at_level(:full) && !vanishing
                flux = MS.inclined_mode_flux(r.a, r.p, r.x, r.l, r.m, r.k)
                @test relerr(flux.infinity.energy_flux, r.pub_Einf) <= 3e-10
                @test relerr(flux.infinity.angular_momentum_flux, r.pub_Lzinf) <= 3e-10
                @test relerr(flux.infinity.carter_constant_flux, r.pub_Qinf) <= 3e-10
                @test relerr(flux.horizon.energy_flux, r.dev_Ehor) <= 1e-9
                # Mode-flux helper and single-mode API are the same computation (<= 1.6e-13 measured).
                @test flux.infinity.energy_flux ≈ inf.energy_flux rtol=1e-11
            end
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Cross-path consistency with the tol = 1e-10 summation tables (frozen self-values from the summation's own
# convolution path; the single-mode default path samples independently): rtol 1e-6 infinity, 1e-5 horizon.
# ---------------------------------------------------------------------------------------------------
@testset "summation-table modes (eccentric, generic, spherical)" begin
    for r in level_subset(TABLE_MODES, 8)
        @testset "$(basename(r.src)) (l,m,k,n)=($(r.l),$(r.m),$(r.k),$(r.n))" begin
            inf = mode_inf(r.l, r.m, r.n, r.k, r.a, r.p, r.e, r.x)
            @test inf.mode.omega ≈ r.omega rtol=1e-12
            @test relerr(inf.energy_flux, r.Einf) <= 1e-6
            if at_level(:full) && isfinite(r.Ehor) && abs(r.Ehor) > 1e-14 * abs(r.Einf)
                hor = mode_hor(r.l, r.m, r.n, r.k, r.a, r.p, r.e, r.x)
                @test relerr(hor.energy_flux, r.Ehor) <= 1e-5
            end
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Harmonic frequencies omega = (m Y_phi + n Y_r + k Y_theta) / Y_t from kerr_geo_orbit (registered-test
# relation), for all orbit types including retrograde x = -1.
# ---------------------------------------------------------------------------------------------------
@testset "harmonic frequencies" begin
    orbits = ((0.5, 11.0, 0.0, 1.0), (0.5, 11.0, 0.0, -1.0), (0.5, 11.0, 0.1, 1.0), (0.5, 11.0, 0.1, -1.0),
              (0.9, 8.0, 0.0, 0.4), (0.9, 8.0, 0.3, 0.4), (0.0, 12.0, 0.2, -0.6), (0.7, 9.0, 0.4, -0.3))
    harmonics = ((2, 2, 0, 0), (2, 1, 1, 0), (3, -2, 2, 1), (2, 2, 1, 1), (4, 3, -1, 2))
    for (a, p, e, x) in level_subset(collect(orbits), 3), (l, m, n, k) in harmonics
        equatorial = abs(x) == 1
        circular = e == 0
        absent = (circular && n != 0) || (equatorial && k != 0)
        mode = Teukolsky_pointparticle_mode(-2, l, m, n, k, a, p, e, x; N=64, Nmax=128)
        if absent
            @test mode.mode.omega === nothing
            @test iszero(mode.amplitude) && iszero(mode.energy_flux)
        else
            @test mode.mode.omega ≈ harmonic_frequency(a, p, e, x, m, n, k) rtol=1e-12
            @test isfinite(mode.amplitude) && isfinite(mode.energy_flux)
            @test mode.energy_flux >= 0
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Absent harmonics and static modes: zero amplitude without a radial solve; absent -> omega === nothing,
# static (omega = 0 at infinity, omega = m Omega_H at the horizon) -> omega kept, fluxes exactly zero
# (src/Inhomogeneous/ConvolutionIntegrals.jl _zero_radiative_mode).
# ---------------------------------------------------------------------------------------------------
@testset "absent and static modes" begin
    for constructor in (Teukolsky_pointparticle_mode, GSN_pointparticle_mode)
        absent = constructor(-2, 2, 2, 1, 0, 0.5, 8.0, 0.0, 1.0)
        @test iszero(absent.amplitude) && iszero(absent.energy_flux)
        @test absent.mode.omega === nothing
        absent_k = constructor(-2, 2, 2, 0, 1, 0.5, 8.0, 0.3, 1.0)
        @test absent_k.mode.omega === nothing && iszero(absent_k.amplitude)
        absent_n = constructor(2, 3, 1, 2, 1, 0.5, 8.0, 0.0, 0.5)
        @test absent_n.mode.omega === nothing && iszero(absent_n.energy_flux)
    end
    for (a, p, e, x) in ((0.5, 8.0, 0.0, 1.0), (0.9, 10.0, 0.0, 0.5), (0.0, 10.0, 0.0, 1.0))
        static_inf = Teukolsky_pointparticle_mode(-2, 2, 0, 0, 0, a, p, e, x)
        @test static_inf.mode.omega == 0
        @test iszero(static_inf.amplitude)
        @test static_inf.energy_flux == 0 && static_inf.angular_momentum_flux == 0 && static_inf.Carter_const_flux == 0
        static_hor = Teukolsky_pointparticle_mode(2, 3, 0, 0, 0, a, p, e, x)
        @test static_hor.mode.omega == 0
        @test static_hor.energy_flux == 0
    end
end

# ---------------------------------------------------------------------------------------------------
# Symmetry: the (l, -m, -n, -k) harmonic has omega -> -omega and equal fluxes (data: (2,-2,0) rows of MODES_56
# equal the (2,2,0) rows to 1e-15). Amplitudes have equal modulus.
# ---------------------------------------------------------------------------------------------------
@testset "(-m, -n, -k) symmetry" begin
    cases = ((0.9, 10.0, 0.0, 0.3, 2, 2, 0, 1), (0.5, 7.5, 0.0, -0.5, 3, 2, 0, 1), (0.9, 10.0, 0.3, 1.0, 2, 1, 2, 0),
             (0.7, 9.0, 0.2, 0.5, 3, 3, -1, 2), (0.0, 10.0, 0.2, 0.5, 2, 1, 1, 1), (0.9, 6.0, 0.0, 0.8, 4, 4, 0, 0))
    for (a, p, e, x, l, m, n, k) in level_subset(collect(cases), 2), s in (-2, 2)
        plus = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x)
        minus = Teukolsky_pointparticle_mode(s, l, -m, -n, -k, a, p, e, x)
        @test minus.mode.omega ≈ -plus.mode.omega rtol=1e-13
        @test minus.energy_flux ≈ plus.energy_flux rtol=1e-10
        @test minus.angular_momentum_flux ≈ plus.angular_momentum_flux rtol=1e-10
        @test isapprox(minus.Carter_const_flux, plus.Carter_const_flux; rtol=1e-10, atol=1e-30)
        @test abs(minus.amplitude) ≈ abs(plus.amplitude) rtol=1e-10
    end
end

# ---------------------------------------------------------------------------------------------------
# GSN_pointparticle_mode = Teukolsky amplitude / C_trans (s = -2) or B_trans (s = +2), same fluxes
# (src/GeneralizedSasakiNakamura.jl GSN_pointparticle_mode).
# ---------------------------------------------------------------------------------------------------
@testset "GSN vs Teukolsky amplitudes" begin
    for (s, l, m, n, k, a, p, e, x) in ((-2, 2, 2, 0, 0, 0.5, 8.0, 0.0, 1.0), (2, 2, 2, 0, 0, 0.5, 8.0, 0.0, 1.0),
                                         (-2, 3, 2, 1, 0, 0.9, 10.0, 0.3, 1.0), (2, 2, 1, 0, 1, 0.9, 10.0, 0.0, 0.5))
        T = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x)
        G = GSN_pointparticle_mode(s, l, m, n, k, a, p, e, x)
        factor = s == -2 ? CF.Ctrans(s, m, a, T.mode.omega, T.mode.lambda) : CF.Btrans(s, m, a, T.mode.omega, T.mode.lambda)
        @test G.amplitude * factor ≈ T.amplitude rtol=1e-12
        @test isequal(G.energy_flux, T.energy_flux)
        @test isequal(G.angular_momentum_flux, T.angular_momentum_flux)
        @test isequal(G.Carter_const_flux, T.Carter_const_flux)
        @test G.mode.omega == T.mode.omega
    end
end

# ---------------------------------------------------------------------------------------------------
# Convolution methods agree on eccentric / generic modes (isem_levin vs isem_trapezoidal default).
# ---------------------------------------------------------------------------------------------------
if at_level(:full)
    @testset "convolution methods agree" begin
        for (s, l, m, n, k, a, p, e, x) in ((-2, 2, 2, 2, 0, 0.9, 10.0, 0.5, 1.0), (-2, 2, 1, 1, 1, 0.9, 10.0, 0.2, 0.3),
                                             (2, 2, 2, 1, 0, 0.0, 10.0, 0.2, 1.0))
            trap = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x; method="isem_trapezoidal")
            levin = Teukolsky_pointparticle_mode(s, l, m, n, k, a, p, e, x; method="isem_levin")
            @test levin.energy_flux ≈ trap.energy_flux rtol=1e-8
            @test levin.mode.omega == trap.mode.omega
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Argument validation.
# ---------------------------------------------------------------------------------------------------
@testset "argument validation" begin
    @test_throws ArgumentError Teukolsky_pointparticle_mode(-2, 2, 2, 1, 1, 0.0, 20.0, 1.1, 0.7)
    @test_throws ErrorException Teukolsky_pointparticle_mode(-2, 2, 2, 0, 0, 0.5, 8.0, 0.0, 1.0; method="bogus")
    @test_throws ErrorException GSN_pointparticle_mode(-2, 2, 2, 0, 0, 0.5, 8.0, 0.0, 1.0; method="bogus")
end

# ---------------------------------------------------------------------------------------------------
# Orbit constants through the KerrGeodesics adapter. Schwarzschild closed forms (e = 0):
# E = (p-2)/sqrt(p(p-3)), L = p/sqrt(p-3), Lz = x L, Q = L^2 (1-x^2) >= 0 (KerrGeodesics 0.3 returned
# Q = -10.714 at (0, 10, 0, 0.5)), Mino frequencies
# Y_theta = L, |Y_phi| = L with sign(x), Y_r = sqrt(p(p-6)/(p-3)), Y_theta/Y_t = p^(-3/2).
# ---------------------------------------------------------------------------------------------------
@testset "Schwarzschild orbit constants" begin
    for p in (6.5, 8.0, 10.0, 20.0, 50.0), x in (1.0, 0.9, 0.5, 0.1, -0.1, -0.5, -1.0)
        orbit = KG.kerr_geo_orbit(0.0, p, 0.0, x)
        L = p / sqrt(p - 3)
        @test orbit["Energy"] ≈ (p - 2) / sqrt(p * (p - 3)) rtol=1e-12
        @test orbit["AngularMomentum"] ≈ x * L rtol=1e-12
        @test isapprox(orbit["CarterConstant"], L^2 * (1 - x^2); rtol=1e-12, atol=1e-12)
        @test orbit["CarterConstant"] >= -1e-12
        f = orbit["Frequencies"]
        @test f["ϒθ"] ≈ L rtol=1e-12
        @test abs(f["ϒϕ"]) ≈ L rtol=1e-12
        @test sign(f["ϒϕ"]) == sign(x)
        @test f["ϒr"] ≈ sqrt(p * (p - 6) / (p - 3)) rtol=1e-12
        @test f["ϒθ"] / f["ϒt"] ≈ p^-1.5 rtol=1e-12
    end
    @test KG.kerr_geo_orbit(0.0, 10.0, 0.0, 0.5)["CarterConstant"] ≈ 75 / 7 rtol=1e-14
end

# Frozen KerrGeodesics 0.4 values for the orbits of the 56-mode table.
ORBITS_DEV = (
    (a=0.9, p=10.0, x=0.8, E=0.952833026330822, L=2.8059200349681057, Q=4.455526873681355),
    (a=0.9, p=10.0, x=0.3, E=0.9546077474302431, L=1.0959001169531528, Q=12.208813278044707),
    (a=0.9, p=10.0, x=-0.5, E=0.9585945373643364, L=-1.9781035223892252, Q=11.787946767474146),
    (a=0.9, p=6.0, x=0.3, E=0.9328154653197315, L=0.9477933611239344, Q=9.17865137400366),
    (a=0.5, p=6.0, x=0.8, E=0.9316386339345527, L=2.493914493279256, Q=3.5104147946493836),
    (a=0.5, p=7.5, x=-0.5, E=0.9501836723940091, L=-1.844499283266101, Q=10.224748628291698),
)
@testset "Kerr orbit constants (frozen KerrGeodesics 0.4 values)" begin
    for o in ORBITS_DEV
        orbit = KG.kerr_geo_orbit(o.a, o.p, 0.0, o.x)
        @test orbit["Energy"] ≈ o.E rtol=1e-12
        @test orbit["AngularMomentum"] ≈ o.L rtol=1e-12
        @test orbit["CarterConstant"] ≈ o.Q rtol=1e-12
    end
end

include("generic_trapezoid_regression.jl")
include("generic_stop_rule_regression.jl")
include("generic_interpolation_regression.jl")
include("generic_levin_regression.jl")
include("eccentric_interpolation_regression.jl")

end # pointparticle_mode
