# Teukolsky_pointparticle_flux totals.
#
# Tolerances: at tol = 1e-8 (full) infinity totals within 2e-8 relative; at
# tol = 1e-6 (quick) within 2e-6. Horizon totals are judged in units of the matching infinity total
# (|dEhor| <= rel |Einf|, likewise Lz, Q): the summation's stop rule is relative to the total flux, so a
# small horizon flux is not resolved to rel of itself.
#

if !@isdefined(GSN)
    include(joinpath(@__DIR__, "..", "load_package.jl"))
end
include(joinpath(@__DIR__, "..", "common", "helpers.jl"))
include(data_file("flux_references.jl"))

const FLUX_TOL = GSN_TEST_LEVEL === :quick ? 1e-6 : 1e-8
const REL = 2 * FLUX_TOL

flux(a, p, e, x; kwargs...) = Teukolsky_pointparticle_flux(a, p, e, x; tol=FLUX_TOL, kwargs...)

"Compare the six totals of `F` with the NamedTuple `want` (fields Einf, Lzinf, Qinf, Ehor, Lzhor, Qhor)."
function check_totals(F, want; rel=REL)
    @test F isa GSN.TeukolskyPointParticleFlux
    pairs_inf = ((F.infinity_energy_flux, want.Einf), (F.infinity_angular_momentum_flux, want.Lzinf),
                 (F.infinity_carter_constant_flux, want.Qinf))
    for (got, ref) in pairs_inf
        if iszero(ref)
            @test abs(got) <= rel * abs(want.Lzinf)
        else
            @test relerr(got, ref) <= rel
        end
    end
    scale_q = max(abs(want.Qinf), abs(want.Lzinf))
    @test abs(F.horizon_energy_flux - want.Ehor) <= rel * abs(want.Einf)
    @test abs(F.horizon_angular_momentum_flux - want.Lzhor) <= rel * abs(want.Lzinf)
    @test abs(F.horizon_carter_constant_flux - want.Qhor) <= rel * scale_q
    return nothing
end

find_orbit(table, a, p, e, x) = only(r for r in table if r.a == a && r.p == p && r.e == e && r.x == x)

"a = 0 rotation identities: E = E_eq, Lz = x Lz_eq, Q = 2 L (1-x^2) Lz_eq, with L = p / sqrt(p - 3 - e^2)."
function rotated_equatorial(p, e, x)
    eq = find_orbit(EQUATORIAL_TOL10, 0.0, p, e, 1.0)
    L = p / sqrt(p - 3 - e^2)
    return (Einf=eq.Einf, Lzinf=x * eq.Lzinf, Qinf=2L * (1 - x^2) * eq.Lzinf,
            Ehor=eq.Ehor, Lzhor=x * eq.Lzhor, Qhor=2L * (1 - x^2) * eq.Lzhor)
end

@testset "flux (tol = $FLUX_TOL)" begin

# ---------------------------------------------------------------------------------------------------
# Circular equatorial, prograde and retrograde (x = -1); references: tol = 1e-10 runs of the circular path
# (frozen self-values, EQUATORIAL_TOL10). a = 0: x = -1 mirrors x = +1.
# ---------------------------------------------------------------------------------------------------
@testset "circular equatorial" begin
    for (a, x) in ((0.9, 1.0), (0.9, -1.0), (0.0, 1.0), (0.0, -1.0))
        @testset "a=$a x=$x" begin
            F = flux(a, 10.0, 0.0, x)
            @test F.orbit_type == :circular
            check_totals(F, find_orbit(EQUATORIAL_TOL10, a, 10.0, 0.0, x))
            @test F.total_modes > 0
            @test F.tolerance == FLUX_TOL
            @test haskey(F.reached, :l_reached_inf)
            text = sprint(show, MIME("text/plain"), F)
            @test startswith(text, "TeukolskyPointParticleFlux(")
            @test occursin("orbit_type = circular", text)
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Eccentric equatorial: (0.9, 12, 0.5, +-1), (0, 10, 0.2, +-1) and the lower-sideband regression
# (0.9, 10, 0.5, 1): the n-shell stop must also cover the omega = 0 side (lower sideband).
# ---------------------------------------------------------------------------------------------------
@testset "eccentric equatorial" begin
    orbits = ((0.9, 12.0, 0.5, 1.0), (0.9, 12.0, 0.5, -1.0), (0.0, 10.0, 0.2, 1.0), (0.0, 10.0, 0.2, -1.0),
              (0.9, 10.0, 0.5, 1.0))
    orbits = GSN_TEST_LEVEL === :quick ? orbits[1:1] : orbits
    for (a, p, e, x) in orbits
        @testset "($a, $p, $e, $x)" begin
            F = flux(a, p, e, x)
            @test F.orbit_type == :eccentric
            check_totals(F, find_orbit(EQUATORIAL_TOL10, a, p, e, x))
            @test haskey(F.reached, :n_reached_inf)
        end
    end
    # a = 0: the x = -1 orbit is the x = +1 orbit rotated; E equal, Lz opposite.
    if at_level(:full)
        Fp = flux(0.0, 10.0, 0.2, 1.0)
        Fm = flux(0.0, 10.0, 0.2, -1.0)
        @test relerr(Fm.infinity_energy_flux, Fp.infinity_energy_flux) <= REL
        @test relerr(-Fm.infinity_angular_momentum_flux, Fp.infinity_angular_momentum_flux) <= REL
    end
end

# ---------------------------------------------------------------------------------------------------
# a = 0 rotation identities: spherical x in (+-0.5, +-0.1) and generic (0, 10, 0.2, +-0.5) against the
# equatorial references.
# ---------------------------------------------------------------------------------------------------
@testset "Schwarzschild rotation identities" begin
    spherical = GSN_TEST_LEVEL === :quick ? (-0.5,) : (0.5, -0.5, 0.1, -0.1)
    for x in spherical
        @testset "spherical x=$x" begin
            F = flux(0.0, 10.0, 0.0, x)
            @test F.orbit_type == :inclined
            check_totals(F, rotated_equatorial(10.0, 0.0, x))
        end
    end
    generic = at_level(:extended) ? (0.5, -0.5) : at_level(:full) ? (0.5,) : ()
    for x in generic
        @testset "generic (0, 10, 0.2, $x)" begin
            F = flux(0.0, 10.0, 0.2, x)
            @test F.orbit_type == :generic
            check_totals(F, rotated_equatorial(10.0, 0.2, x))
        end
    end
    if at_level(:extended)
        # The retrograde high-n m spectrum can dip below the cutoff before rising.
        for x in (-0.5, 0.5)
            F = flux(0.0, 10.0, 0.5, x)
            check_totals(F, rotated_equatorial(10.0, 0.5, x))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Kerr spherical orbits. Reference: independent full-grid sum (SPHERICAL_FULLGRID) when present, otherwise the
# tol = 1e-10 candidate run (CANDIDATE_TOL10, frozen self-reference; agrees with the full grid to <= 3e-11 at
# the orbits where both exist). Rows without an independent full-grid reference are marked as skipped.
# ---------------------------------------------------------------------------------------------------
SPHERICAL_ORBITS = [(r.a, r.p, r.x) for r in CANDIDATE_TOL10 if r.e == 0]
@testset "Kerr spherical" begin
    orbits = GSN_TEST_LEVEL === :quick ? [(0.9, 10.0, 0.3)] : SPHERICAL_ORBITS
    for (a, p, x) in orbits
        @testset "($a, $p, 0, $x)" begin
            fullgrid = [r for r in SPHERICAL_FULLGRID if r.a == a && r.p == p && r.x == x]
            reference = isempty(fullgrid) ? find_orbit(CANDIDATE_TOL10, a, p, 0.0, x) : only(fullgrid)
            if isempty(fullgrid)
                @test_skip "independent full-grid reference for ($a, $p, 0, $x)" == ""
            end
            F = flux(a, p, 0.0, x)
            @test F.orbit_type == :inclined
            check_totals(F, reference)
            @test haskey(F.reached, :k_reached_inf)
            # Sign of the angular-momentum flux follows the orbit's direction.
            @test sign(F.infinity_angular_momentum_flux) == sign(x)
            @test F.infinity_carter_constant_flux > 0
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Kerr generic (0.9, 10, 0.2, +-0.3), extended level only (~200-250 s each at tol 1e-8). Reference: independent
# full-grid sums (GENERIC_FULLGRID, no stop rule; edge layers <= 7.4e-11).
# ---------------------------------------------------------------------------------------------------
if at_level(:extended)
    @testset "Kerr generic" begin
        for x in (0.3, -0.3)
            F = flux(0.9, 10.0, 0.2, x)
            @test F.orbit_type == :generic
            check_totals(F, find_orbit(GENERIC_FULLGRID, 0.9, 10.0, 0.2, x))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Unbound / unstable orbits return nothing with a warning (src/GeneralizedSasakiNakamura.jl
# _teukolsky_flux_bound_orbit). (0.9, 10, 0.5, -1) is not bound.
# ---------------------------------------------------------------------------------------------------
@testset "unbound orbits" begin
    result = @test_logs (:warn, r"bound orbit") match_mode=:any Teukolsky_pointparticle_flux(0.9, 10.0, 0.5, -1.0)
    @test result === nothing
    # Unstable circular orbits inside the ISCO are accepted by design (_teukolsky_flux_bound_orbit: Critical,
    # Circular, Elliptic), e.g. Schwarzschild p = 5.
    @test Teukolsky_pointparticle_flux(0.0, 5.0, 0.0, 1.0; tol=1e-6) !== nothing
end

# Leading quadrupole limit, not a finite-frequency accuracy certificate.
@testset "large-p circular quadrupole limit (items 1/12)" begin
    previous = Inf
    for p in (1e5, 1e6, 1e7, 2e7, 3e7, 1.6e8, 1e9)
        mode = Teukolsky_pointparticle_mode(-2, 2, 2, 0, 0, 0.5, p, 0.0, 1.0)
        v = (mode.mode.omega / 2)^(1 / 3)
        correction = abs(mode.energy_flux / ((16 / 5) * v^10) - 1)
        @test isfinite(correction)
        @test correction < previous
        @test correction <= 6 / p
        previous = correction
    end
    total = Teukolsky_pointparticle_flux(0.5, 2e7, 0.0, 1.0)
    @test abs(total.infinity_energy_flux / ((32 / 5) * (2e7)^(-5)) - 1) <= 6 / 2e7
    @test total.total_modes > 0
end

# Fixed l=2:7, n=-30:30 sums, independently refined at relative mode tolerance 1e-11.
# The outer l/n layers contribute less than 7e-14; no summation stopping rule is used.
@testset "large-p eccentric full-box totals" begin
    for (p, infinity, horizon) in (
        (1e3, 7.069635134875905e-15, -6.639717219612091e-23),
        (1e4, 7.091345568411108e-20, -2.092300164964513e-30),
        (3e4, 2.9189641258868893e-22, -5.5217849374485675e-34),
    )
        total = Teukolsky_pointparticle_flux(0.5, p, 0.3, 1.0; tol=1e-8)
        @test abs(total.infinity_energy_flux / infinity - 1) <= 1e-8
        @test abs(total.horizon_energy_flux / horizon - 1) <= 1e-8
    end
end

# ---------------------------------------------------------------------------------------------------
# Argument validation (ModeSummation.jl): thrown before any mode is computed.
# ---------------------------------------------------------------------------------------------------
@testset "argument validation" begin
    ecc = (0.9, 10.0, 0.3, 1.0)
    inc = (0.9, 10.0, 0.0, 0.5)
    gen = (0.9, 10.0, 0.3, 0.5)
    for orbit in (ecc, inc, gen)
        @test_throws ArgumentError Teukolsky_pointparticle_flux(orbit...; tol=0.0)
        @test_throws ArgumentError Teukolsky_pointparticle_flux(orbit...; tol=-1e-8)
        @test_throws ArgumentError Teukolsky_pointparticle_flux(orbit...; lmax=1)
        @test_throws ArgumentError Teukolsky_pointparticle_flux(orbit...; minimum_consecutive=0)
    end
    @test_throws ArgumentError Teukolsky_pointparticle_flux(0.9, 10.0, 0.0, 1.0; minimum_consecutive=0)
    for orbit in (ecc, gen)
        @test_throws ArgumentError Teukolsky_pointparticle_flux(orbit...; nmax=-1)
        @test_throws ArgumentError Teukolsky_pointparticle_flux(orbit...; levin_local_n=0)
        @test_throws ArgumentError Teukolsky_pointparticle_flux(orbit...; levin_max_depth=-1)
    end
    # Sampling-grid sizes (checked at the start of eccentric_mode_summation / inclined_mode_summation).
    @test_throws ArgumentError Teukolsky_pointparticle_flux(ecc...; N0=48)
    @test_throws ArgumentError Teukolsky_pointparticle_flux(ecc...; Nmax=1000)
    @test_throws ArgumentError Teukolsky_pointparticle_flux(ecc...; N0=2^15)          # N0 > Nmax
    @test_throws ArgumentError Teukolsky_pointparticle_flux(inc...; K0=12)
    @test_throws ArgumentError Teukolsky_pointparticle_flux(inc...; Kmax=1000)
    @test_throws ArgumentError Teukolsky_pointparticle_flux(inc...; K0=2^13)          # K0 > Kmax
    for orbit in (inc, gen)
        @test_throws ArgumentError Teukolsky_pointparticle_flux(orbit...; kmax=-1)
    end
    @test_throws ArgumentError Teukolsky_pointparticle_flux(ecc...; tail_levin=:sometimes)
end

end # flux
