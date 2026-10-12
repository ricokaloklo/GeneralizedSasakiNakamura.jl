# One testset per previously observed failure or known limitation. Cases whose main test lives in another
# group are only re-pinned here when a short, direct reproducer exists.

if !@isdefined(GSN)
    include(joinpath(@__DIR__, "..", "load_package.jl"))
end
include(joinpath(@__DIR__, "..", "common", "helpers.jl"))
include(data_file("radial_references.jl"))

const D = GSN.ISEM
dominant(x, y) = max(abs(x), abs(y))
point(id) = only(p for p in STRONG_DAMPING_POINTS if p.id == id)
grid_point(id) = only(g for g in GRID_POINTS if g.id == id)

@testset "regression" begin

# Snapshot of earlier registered tests; current Pkg.test remains separate.
include(joinpath(@__DIR__, "registered_suite.jl"))

# a = 0.7, l = m = 2, n = 10 root, UP C_up buried at the junction (1.78e-6 in C_ref units).
@testset "buried C_up at the a=0.7 n10 root" begin
    q = point("Q_a07_m2_n10_at_f64")
    R = Teukolsky_radial(-2, 2, 2, q.a, q.omega, UP)
    cup, cref = teukolsky_ratios(R)
    scale = dominant(q.cup, q.cref)
    @test scaled_err(cref, q.cref, scale) <= 1e-9
    @test scaled_err(cup, q.cup, scale) <= 1e-5
    @test_broken scaled_err(cup, q.cup, scale) <= 1e-10   # buried horizon coefficient
end

# a = 0 n12 / n20 roots, UP C_up (1.0e-7 / 1.2e-7 in C_ref units).
@testset "a=0 high-overtone UP C_up" begin
    for (id, floor) in (("Q_a0_n12_at_f64", 1e-6), ("Q_a0_n20_at_f64", 1e-6), ("Q_a0_n12_off_f64", 5e-9))
        q = point(id)
        R = Teukolsky_radial(-2, 2, 2, q.a, q.omega, UP)
        cup, cref = teukolsky_ratios(R)
        scale = dominant(q.cup, q.cref)
        @test scaled_err(cup, q.cup, scale) <= floor
        @test scaled_err(cref, q.cref, scale) <= 1e-9
    end
end

# Suppressed reflection: compare the per-coefficient relative error with the independent reference.
@testset "reflection error coverage" begin
    for id in (32, 62)
        g = grid_point(id)
        route = D.radial_route(g.s, g.l, g.m, g.a, g.omega, :IN)
        errors = D.Matching.amplitude_errors(route)
        ref = only(r for r in GRID_REFERENCES if r.id == id && r.branch == "IN")
        _, reflected = teukolsky_ratios(D.teukolsky_radial_function(route))
        @test relerr(reflected, ref.reflection) <= errors[2]
        @test isfinite(errors[1])
    end
end

# This checks the error-slot interface only; strong-damping coverage remains unresolved.
@testset "R(r) error slot" begin
    route = D.radial_route(-2, 2, 2, 0.7, 0.5 - 0.1im, :IN)
    R = D.teukolsky_radial_function(route)
    state = R.Teukolsky_solution(6.0)
    @test length(state) == 4
    @test all(isfinite, state[1:3])
    # The error slot carries a positive estimate.
    @test state[4] > 0
end

# R(10) at real frequency (grid id45 UP, id70 IN, id111 UP): normalization floor.
@testset "R(r) floor at real frequency" begin
    refs = Dict((r.id, r.branch) => r for r in GRID_REFERENCES)
    for (id, branch) in ((45, "UP"), (70, "IN"), (111, "UP"))
        g = grid_point(id)
        ref = refs[(id, branch)]
        R = Teukolsky_radial(g.s, g.l, g.m, g.a, g.omega, branch == "IN" ? IN : UP)
        err = relerr(R(10.0) / R.transmission_amplitude, ref.R10)
        @test err <= 1e-6
        @test err <= 1e-12
    end
end

# Error estimates short by x1.12-x1.32 (two-ray step budget, ys tables, MST eps floor).
@testset "under-reports stay within x2" begin
    for (id, branch, j) in (("D_l2_a0.9_w2.0-0.1i", :UP, 1), ("Q_a07_m2_n10_at_f64", :UP, 2),
                            ("Q_a0_n9_off_f64", :UP, 1), ("Q_a0_n20_off_f64", :UP, 1), ("Q_a09_m2_n0_off_f64", :IN, 1))
        q = point(id)
        route = D.radial_route(q.s, q.l, q.m, q.a, q.omega, branch)
        values = teukolsky_ratios(D.teukolsky_radial_function(route))
        refs = branch == :IN ? (q.binc, q.bref) : (q.cup, q.cref)
        reported = D.Matching.amplitude_errors(route)[j]
        actual = relerr(values[j], refs[j])
        @test actual <= 2 * reported
        @test_broken actual <= reported   # amplitude error estimate
    end
end

# Real-frequency suppressed amplitudes, in dominant units on the real-frequency (R_) validation points.
@testset "suppressed real-frequency amplitudes" begin
    for id in ("R_l3_a0.7_w3.0", "R_l2_a0.9_w2.0", "R_l2_a0.9_w3.0", "R_l4_a0.9_w3.0")
        q = point(id)
        Rin = Teukolsky_radial(-2, q.l, q.m, q.a, q.omega, IN)
        binc, bref = teukolsky_ratios(Rin)
        @test scaled_err(bref, q.bref, dominant(q.binc, q.bref)) <= 1e-9
        @test relerr(binc, q.binc) <= 1e-9
    end
end

# Open: mano_suzuki_takasugi_nu_complex rejects nu at grid id32/id62 by Float64 evaluation noise. Not observable through public
# values (MST condition 1e10-1e11 there); needs an internal reproducer with a BigFloat nu reference.
@testset "nu gate (open)" begin
    @test_skip false
end

# B_inc at grid 24/34/54/64 IN (omega = -0.1-0.05i, m flipped) is the exact conjugate of the Re omega > 0
# partner.
@testset "Re omega < 0 two-ray cells" begin
    for id in (24, 34, 54, 64)
        g = grid_point(id)
        P = Teukolsky_radial(g.s, g.l, g.m, g.a, g.omega, IN)
        N = Teukolsky_radial(g.s, g.l, -g.m, g.a, -conj(g.omega), IN)
        pi_, pr = teukolsky_ratios(P)
        ni, nr = teukolsky_ratios(N)
        @test abs(ni - conj(pi_)) <= 1e-11 * dominant(pi_, pr)
        @test abs(nr - conj(pr)) <= 1e-11 * dominant(pi_, pr)
    end
end

# Open: explicit-control calls at Re omega < 0 are not reflected at the route level, so
# the two signs are separate two-ray computations. Their disagreement must stay within the two reports.
@testset "explicit controls at Re omega < 0" begin
    g = grid_point(24)
    positive = D.radial_route(g.s, g.l, g.m, g.a, g.omega, :IN; N=48)
    negative = D.radial_route(g.s, g.l, -g.m, g.a, -conj(g.omega), :IN; N=48)
    rp = teukolsky_ratios(D.teukolsky_radial_function(positive))
    rn = teukolsky_ratios(D.teukolsky_radial_function(negative))
    ep = D.Matching.amplitude_errors(positive)
    en = D.Matching.amplitude_errors(negative)
    for j in 1:2
        @test abs(rn[j] - conj(rp[j])) <= en[j] * abs(rn[j]) + ep[j] * abs(rp[j]) + 1e-14 * dominant(rp...)
    end
    @test_broken rn[1] == conj(rp[1])   # explicit-control reflection
end

# Explicit xm with method = "auto" rethrows instead of silently using the legacy ODE solver (which returns
# C_up 4.52e4 off in C_ref units at every xm >= 0.5).
@testset "no silent legacy fallback with explicit xm" begin
    q = point("Q_a07_m2_n10_at_f64")
    @test_throws ErrorException GSN_radial(-2, 2, 2, q.a, q.omega, UP; xm=0.7)
    @test_throws ErrorException Teukolsky_radial(-2, 2, 2, q.a, q.omega, UP; xm=0.8)
    X = GSN_radial(-2, 2, 2, q.a, q.omega, UP)
    @test X.method == "ISEM"
end

# The p-consensus route reports |p - two-ray| + E(two-ray), which must be finite.
@testset "finite reports on p-consensus routes" begin
    for id in (24, 34, 54, 64, 84, 94, 114, 124, 154, 180)
        g = grid_point(id)
        route = D.radial_route(g.s, g.l, g.m, g.a, g.omega, :IN)
        @test all(isfinite, D.Matching.amplitude_errors(route))
    end
end

# Default Teukolsky_radial and GSN_radial use the same ISEM route and tolerance and agree to rounding
# (an ODE tolerance on this route gives C_ref/C_up off by 2.0e-9 at a = 0, w = 0.3-3i).
@testset "default tolerance passthrough" begin
    for (s, l, m, a, w) in ((-2, 2, 2, 0.0, 0.3 - 3.0im), (-2, 2, 2, 0.7, 0.5 - 0.1im), (2, 3, 1, 0.9, 1.0 - 0.5im))
        for bc in (IN, UP)
            X = GSN_radial(s, l, m, a, w, bc)
            Xs = GSN_radial(s, l, m, a, w, bc; method="ISEM")
            @test X.incidence_amplitude == Xs.incidence_amplitude
            @test X.reflection_amplitude == Xs.reflection_amplitude
            # Teukolsky: auto converts the GSN amplitudes, strict reads the route's Teukolsky amplitudes.
            R = Teukolsky_radial(s, l, m, a, w, bc)
            Rs = Teukolsky_radial(s, l, m, a, w, bc; method="ISEM")
            ri, rr = teukolsky_ratios(R)
            si, sr = teukolsky_ratios(Rs)
            @test abs(ri - si) <= 1e-12 * dominant(si, sr)
            @test abs(rr - sr) <= 1e-12 * dominant(si, sr)
        end
    end
end

# Detailed QNM normalization bridge: residual and Richardson drifts within their gates.
@testset "detailed QNM normalization bridge" begin
    detailed = qnm(0.68, -2, 2, 2, 0; detailed=true)
    @test detailed.status == :accepted
    @test detailed.excitation.bridge_residual <= 1e-10
    @test detailed.excitation.direction_drift <= 1e-7
    @test detailed.excitation.step_drift <= 1e-7
end

# Weak-damping UP state: R(5), R(10) at grid 64/134/38/164/8 UP (omega = 0.1-0.05i) against the independent
# 896-bit grid references (GRID_REFERENCES).
@testset "weak-damping UP state" begin
    refs = Dict((r.id, r.branch) => r for r in GRID_REFERENCES)
    for id in (64, 134, 38, 164, 8)
        g = grid_point(id)
        ref = refs[(id, "UP")]
        R = Teukolsky_radial(g.s, g.l, g.m, g.a, g.omega, UP)
        @test relerr(R(10.0) / R.transmission_amplitude, ref.R10) <= 1e-8
        @test relerr(R(5.0) / R.transmission_amplitude, ref.R5) <= 1e-8
    end
end

# KerrGeodesics 0.3 gave Q = -10.714 at (a, p, e, x) = (0, 10, 0, 0.5); the correct value is +75/7.
@testset "KerrGeodesics Carter-constant sign" begin
    for x in (0.5, -0.5, 0.2)
        Q = GSN.KerrGeodesics.kerr_geo_orbit(0.0, 10.0, 0.0, x)["CarterConstant"]
        @test Q ≈ 100 / 7 * (1 - x^2) rtol=1e-13
    end
end

# Any stop rule must keep a retrograde spherical orbit consistent with its a = 0 rotation identity
# (an l-streak stop rule fails this check).
@testset "retrograde spherical stop rule (a = 0)" begin
    F = Teukolsky_pointparticle_flux(0.0, 10.0, 0.0, -0.5; tol=1e-6)
    E_eq = 6.150372548996145e-5          # a = 0, p = 10 circular prograde total (EQUATORIAL_TOL10)
    Lz_eq = 0.001944918571340336
    @test relerr(F.infinity_energy_flux, E_eq) <= 2e-6
    @test relerr(F.infinity_angular_momentum_flux, -0.5 * Lz_eq) <= 2e-6
end

end # regression

# GSN_pointparticle_mode and Teukolsky_pointparticle_mode must use the same grid defaults; otherwise the two
# public routes disagree (Carter flux 1.2e-2 for this generic mode, 3e-8 for spherical k = 1).
@testset "GSN vs Teukolsky point-particle mode defaults" begin
    for c in ((-2, 2, 2, 1, 1, 0.9, 10.0, 0.2, 0.3), (2, 2, 1, 0, 1, 0.9, 10.0, 0.0, 0.5))
        T = Teukolsky_pointparticle_mode(c...)
        G = GSN_pointparticle_mode(c...)
        @test G.energy_flux ≈ T.energy_flux rtol=1e-10
        @test isequal(G.Carter_const_flux, T.Carter_const_flux)
    end
end
