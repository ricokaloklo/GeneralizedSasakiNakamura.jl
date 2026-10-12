# ISEM paths: complex contours, endpoints, reflections, junctions and reported errors.

if !@isdefined(GSN)
    include(joinpath(@__DIR__, "..", "load_package.jl"))
end
include(joinpath(@__DIR__, "..", "common", "helpers.jl"))
include(data_file("radial_references.jl"))
include(data_file("positive_axis_references.jl"))

const D = GSN.ISEM
const MST = D.MSTInfinity

dominant(x, y) = max(abs(x), abs(y))
point(id) = only(p for p in STRONG_DAMPING_POINTS if p.id == id)
grid_point(id) = only(g for g in GRID_POINTS if g.id == id)

@testset "isem" begin

# ---------------------------------------------------------------------------------------------------
# Near the positive imaginary axis, a = 0, omega = eps + 0.3i: a junction that puts the infinity ray through
# the r* image of r = 0 returns O(1)-wrong amplitudes while reporting ~1e-12. Independent 384-bit references
# (POSITIVE_AXIS_ROWS, POSITIVE_AXIS_POINTS); measured errors IN <= 3.02e-11, UP <= 8.65e-15. Tolerance 1e-9.
# ---------------------------------------------------------------------------------------------------
@testset "near the positive imaginary axis" begin
    for row in POSITIVE_AXIS_ROWS
        bc = row.branch == "IN" ? IN : UP
        @testset "$(row.id) $(row.branch)" begin
            R = Teukolsky_radial(-2, 2, 2, row.a, row.omega, bc)
            inc, refl = teukolsky_ratios(R)
            @test relerr(inc, row.incidence) <= 1e-9
            @test relerr(refl, row.reflection) <= 1e-9
            @test isapprox(R.mode.lambda, row.lambda; rtol=1e-12)
            route = D.radial_route(-2, 2, 2, row.a, row.omega, Symbol(row.branch))
            errors = D.Matching.amplitude_errors(route)
            dv = teukolsky_ratios(D.teukolsky_radial_function(route))
            # The error estimate is of the right order (e.g. reported 4.3e-11 vs actual 4.07e-11).
            @test relerr(dv[1], row.incidence) <= max(2 * errors[1], 1e-13)
            @test relerr(dv[2], row.reflection) <= max(2 * errors[2], 1e-13)
            # At a = 0 the solution does not depend on m, so the exact reflection gives the mirrored point:
            # R_{m}(-conj(omega)) = conj R_{-m}(omega) = conj R_{m}(omega).
            if row.a == 0
                M = Teukolsky_radial(-2, 2, 2, 0.0, -conj(row.omega), bc)
                mi, mr = teukolsky_ratios(M)
                @test relerr(mi, conj(row.incidence)) <= 1e-9
                @test relerr(mr, conj(row.reflection)) <= 1e-9
            end
        end
    end
    @testset "on the axis (Re omega = +0.0) and near it" begin
        for row in POSITIVE_AXIS_POINTS
            bc = row.branch == "IN" ? IN : UP
            R = Teukolsky_radial(-2, 2, 2, row.a, row.omega, bc)
            @test relerr(first(teukolsky_ratios(R)), row.incidence) <= 1e-9
        end
    end
    # GSN_radial uses the same direct route and must agree after amplitude conversion.
    for (a, w) in ((0.0, 1e-10 + 0.3im), (0.0, 1e-8 + 0.3im), (0.7, 1e-10 + 0.3im))
        for bc in (IN, UP)
            X = GSN_radial(-2, 2, 2, a, w, bc)
            R = Teukolsky_radial(-2, 2, 2, a, w, bc)
            CFm = GSN.ConversionFactors
            lam = X.mode.lambda
            f = bc == IN ? CFm.Binc(-2, 2, a, w, lam) / CFm.Btrans(-2, 2, a, w, lam) :
                CFm.Cinc(-2, 2, a, w, lam) / CFm.Ctrans(-2, 2, a, w, lam)
            @test isapprox(f * X.incidence_amplitude / X.transmission_amplitude,
                first(teukolsky_ratios(R)); rtol=1e-10)
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Negative frequencies through the exact reflection R_{s,l,m}(omega) = conj R_{s,l,-m}(-conj(omega)).
# References: independent validation rows (real and weakly damped, s = -2) carried to (-m, -conj(omega)).
# ---------------------------------------------------------------------------------------------------
@testset "negative frequencies through the exact reflection" begin
    rows = [p for p in STRONG_DAMPING_POINTS if startswith(p.id, "R_") || startswith(p.id, "D_")]
    for p in level_subset(rows, 5)
        @testset "$(p.id) reflected" begin
            w = -conj(p.omega)
            Rin = Teukolsky_radial(p.s, p.l, -p.m, p.a, w, IN)
            Rup = Teukolsky_radial(p.s, p.l, -p.m, p.a, w, UP)
            values = (teukolsky_ratios(Rin)..., teukolsky_ratios(Rup)...)
            refs = conj.((p.binc, p.bref, p.cup, p.cref))
            scales = (dominant(p.binc, p.bref), dominant(p.binc, p.bref), dominant(p.cup, p.cref), dominant(p.cup, p.cref))
            for j in 1:4
                @test scaled_err(values[j], refs[j], scales[j]) <= p.tol[j]
            end
            lam_pos = p.lam !== nothing ? p.lam :
                p.E - 2 * p.a * p.m * p.omega + p.a^2 * p.omega^2 - p.s * (p.s + 1)
            @test abs(Rin.mode.lambda - conj(lam_pos)) <= 1e-12 * max(1, abs(lam_pos))
        end
    end
    # Registered-suite cases (test/runtests.jl "Negative-frequency amplitudes").
    for (l, m, omega, incidence, reflection) in (
        (2, -2, -1.0 - 0.05im, 1.2856521453922971 - 0.44862664931401232im,
         -2.9559869085975771e-4 - 1.4693045725576587e-3im),
        (3, -3, -1.0, 0.5043157658448172 + 0.24588840385554175im,
         0.00012198157685068585 + 0.18138339086652683im))
        R = Teukolsky_radial(-2, l, m, 0.9, omega, IN)
        @test R.incidence_amplitude / R.transmission_amplitude ≈ incidence rtol=1e-9
        @test R.reflection_amplitude / R.transmission_amplitude ≈ reflection rtol=1e-9
    end
    # Mirror qnm(0.68) n0 frequency: B_inc vanishes.
    X = GSN_radial(-2, 2, 2, 0.68, -0.3111628576854571 - 0.0887546633199609im, IN)
    @test abs(X.incidence_amplitude) / max(1, abs(X.reflection_amplitude), abs(X.transmission_amplitude)) < 1e-10
    # Real omega < 0 with Im = +0.0 and -0.0 give identical values.
    for bc in (IN, UP)
        A = Teukolsky_radial(-2, 3, 3, 0.9, complex(-1.0, 0.0), bc)
        B = Teukolsky_radial(-2, 3, 3, 0.9, complex(-1.0, -0.0), bc)
        @test A.incidence_amplitude == B.incidence_amplitude
        @test A.reflection_amplitude == B.reflection_amplitude
    end
    # Pairs straddling Re omega = 0 (+-1e-3, +-1e-6) are exact conjugates.
    for (a, im_part) in ((0.7, -2.0), (0.0, -1.0)), re in (1e-3, 1e-6), bc in (IN, UP)
        P = Teukolsky_radial(-2, 2, 2, a, complex(re, im_part), bc)
        N = Teukolsky_radial(-2, 2, -2, a, complex(-re, im_part), bc)
        @test N.incidence_amplitude / N.transmission_amplitude ≈ conj(P.incidence_amplitude / P.transmission_amplitude) rtol=1e-12
        @test N.reflection_amplitude / N.transmission_amplitude ≈ conj(P.reflection_amplitude / P.transmission_amplitude) rtol=1e-12
    end
end

# ---------------------------------------------------------------------------------------------------
# Reflected amplitude uncertainty: the weak-damping p-consensus IN routes (grid 24/34/54/64/84/94/114/124/
# 154/180) report finite errors, identical across the reflection (test/runtests.jl "Reflected amplitude
# uncertainty", extended to all ten cells).
# ---------------------------------------------------------------------------------------------------
@testset "reported errors across the reflection" begin
    for id in (24, 34, 54, 64, 84, 94, 114, 124, 154, 180)
        g = grid_point(id)
        positive = D.radial_route(g.s, g.l, g.m, g.a, g.omega, :IN)
        negative = D.radial_route(g.s, g.l, -g.m, g.a, -conj(g.omega), :IN)
        ep = D.Matching.amplitude_errors(positive)
        en = D.Matching.amplitude_errors(negative)
        @test all(isfinite, ep)
        @test ep == en
        rp = teukolsky_ratios(D.teukolsky_radial_function(positive))
        rn = teukolsky_ratios(D.teukolsky_radial_function(negative))
        @test rn[1] == conj(rp[1]) || isapprox(rn[1], conj(rp[1]); rtol=1e-14)
        @test rn[2] == conj(rp[2]) || isapprox(rn[2], conj(rp[2]); rtol=1e-14)
    end
    # Reports are relative to the individual reflection coefficient, not the dominant amplitude.
    for id in (32, 62)
        g = grid_point(id)
        route = D.radial_route(g.s, g.l, g.m, g.a, g.omega, :IN)
        ref = only(r for r in GRID_REFERENCES if r.id == id && r.branch == "IN")
        _, reflected = teukolsky_ratios(D.teukolsky_radial_function(route))
        @test relerr(reflected, ref.reflection) <= D.Matching.amplitude_errors(route)[2]
    end
end

# ---------------------------------------------------------------------------------------------------
# Explicit junction at the a = 0.7, l = m = 2, n = 10 root (UP): xm >= 0.5 cannot be certified (a real-axis
# horizon path gives C_up off by 13..3e9 in C_ref units while reporting 8e-8..2e-2). The direct route raises;
# the public auto method rethrows instead of silently using the legacy solver (main-module
# `xm === nothing || rethrow()`).
# ---------------------------------------------------------------------------------------------------
@testset "explicit junction" begin
    q = point("Q_a07_m2_n10_at_f64")
    for xm in (0.5, 0.6, 0.7, 0.8)
        @test_throws ErrorException D.radial_route(-2, 2, 2, q.a, q.omega, :UP; xm=xm)
    end
    for radial in (GSN_radial, Teukolsky_radial)
        @test_throws ErrorException radial(-2, 2, 2, 0.7, 0.4198385710926635 - 2.09355961067422im, UP; xm=0.6)
        @test_throws ErrorException radial(-2, 2, 2, q.a, q.omega, UP; xm=0.6, method="ISEM")
    end
    # Junctions that certify still give the reference C_ref; C_up stays within the buried-C_up tolerance class.
    scale = dominant(q.cup, q.cref)
    for xm in (0.3, 0.4)
        R = Teukolsky_radial(-2, 2, 2, q.a, q.omega, UP; xm=xm)
        cup, cref = teukolsky_ratios(R)
        @test scaled_err(cref, q.cref, scale) <= 1e-9
        # At xm=0.4 the error is 1.8e-5 in units of the small C_ref, about 5e-10 in C_trans units.
        err_cup = scaled_err(cup, q.cup, scale)
        if xm == 0.4
            @test err_cup <= 1.816e-5
            @test_broken err_cup <= 1e-5   # explicit-junction C_up limit
        else
            @test err_cup <= 1e-5
        end
    end
    # xm must lie in (0, 1); keywords of the 0.9 ISEM (rhom, TS switches, asymptotic patches) no longer exist.
    for radial in (GSN_radial, Teukolsky_radial), xm in (-1.0, 0.0, 1.0, 1.5)
        @test_throws ArgumentError radial(-2, 2, 2, 0.68, 0.3, IN; xm=xm)
        @test_throws ArgumentError radial(-2, 2, 2, 0.68, 0.3, IN; xm=xm, method="ISEM")
    end
    for radial in (GSN_radial, Teukolsky_radial)
        @test_throws MethodError radial(-2, 2, 2, 0.68, 0.3, IN; rhom=1.0)
        @test_throws MethodError radial(-2, 2, 2, 0.68, 0.3, IN; TSinInf=true)
        @test_throws MethodError radial(-2, 2, 2, 0.68, 0.3, IN; TSoutHor=false, method="ISEM")
    end
    @test_throws MethodError GSN_radial(-2, 2, 2, 0.68, 0.3, IN; use_gsn_asymptotic_patches=true)
end

# ---------------------------------------------------------------------------------------------------
# Known limitation: pinned UP C_up at the Q_a09_m2_n5 root with N = 96, tol = 1e-14: error 4.85e-8 in
# C_ref units, reported 8.81e-9 (x5.5 short).
# ---------------------------------------------------------------------------------------------------
@testset "pinned UP at the a=0.9 n5 root" begin
    q = point("Q_a09_m2_n5_at_f64")
    route = D.radial_route(-2, 2, 2, q.a, q.omega, :UP; N=96, tol=1e-14)
    cup, cref = teukolsky_ratios(D.teukolsky_radial_function(route))
    errors = D.Matching.amplitude_errors(route)
    scale = dominant(q.cup, q.cref)
    @test scaled_err(cref, q.cref, scale) <= 1e-9
    @test scaled_err(cup, q.cup, scale) <= 1e-6                 # regression bound (4.85e-8 measured)
    @test_broken scaled_err(cup, q.cup, scale) <= 2 * errors[1] * abs(cup) / scale   # honest error report
    # The public default is better through the MST overlay (8.65e-12 measured).
    R = Teukolsky_radial(-2, 2, 2, q.a, q.omega, UP)
    @test scaled_err(first(teukolsky_ratios(R)), q.cup, scale) <= 1e-9
end

# ---------------------------------------------------------------------------------------------------
# Real-frequency MST branch: nu is complex outside |cos 2 pi nu| <= 1 (monodromy trichotomy).
# Reference: independent validation rows (27 real points); worst measured error 2.8e-10 relative.
# ---------------------------------------------------------------------------------------------------
@testset "real-frequency MST principal amplitudes" begin
    rows = [p for p in STRONG_DAMPING_POINTS if startswith(p.id, "R_")]
    for p in level_subset(rows, 5)
        @testset "$(p.id)" begin
            params = D.gsn_parameters(p.s, p.l, p.m, p.a, real(p.omega))
            nu = MST._mst_data(params).params.nu
            @test isfinite(nu)
            amplitudes = MST.mst_principal_amplitudes(params, :IN)
            @test relerr(amplitudes.teuk[1], p.binc) <= 1e-9
            @test scaled_err(amplitudes.teuk[2], p.bref, dominant(p.binc, p.bref)) <= 1e-9
        end
    end
    # Registered case: R_l3_a0.9_w1 with nu = 2.5 + 1.262i (mod integers, nu -> -nu-1).
    params = D.gsn_parameters(-2, 3, 3, 0.9, 1.0)
    nu = MST._mst_data(params).params.nu
    @test abs(imag(nu)) ≈ 1.2619741165113905 atol=1e-10 rtol=0
    @test MST.mst_principal_amplitudes(params, :IN).teuk[1] ≈ 0.5043157658448172 - 0.24588840385554175im rtol=1e-10
end

# ---------------------------------------------------------------------------------------------------
# Angular eigenvalue input error. The SWSH eigenvalue carries a few ulp (8e-16 absolute at a = 0.7, l = 2,
# omega = 1); the tests fix that floor so that a regression in the angular input is visible.
# ---------------------------------------------------------------------------------------------------
@testset "angular eigenvalue input" begin
    for p in STRONG_DAMPING_POINTS
        p.lam === nothing && continue
        c = isreal(p.omega) ? p.a * real(p.omega) : p.a * p.omega
        lam = GSN._swsh_eigenvalue(p.s, p.l, p.m, c)
        tol = isreal(p.omega) ? 1e-13 : 1e-12
        @test abs(lam - p.lam) <= tol * max(1, abs(p.lam))
        if isreal(p.omega)
            @test isapprox(D.swsh_eigenvalue(p.s, p.l, p.m, p.a * real(p.omega)), p.lam; rtol=1e-13)
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Default tolerance of Teukolsky_radial: the default call must use the ISEM tolerance, not the ODE
# tolerance.
# ---------------------------------------------------------------------------------------------------
@testset "default Teukolsky tolerance" begin
    for c in ((w=0.3 - 3.0im, cref=0.000034049211003554543018 - 0.000017510919791042320907im,
               cup=5.0661768763345224235 + 0.25720529286268218959im),
              (w=0.37467168441804183 - 0.0889623156889357im, cref=-0.79350698537260631 + 1.0219268994327333im,
               cup=-0.06972551142642598 - 0.053077386482722241im))
        Rin, Rup = Teukolsky_radial(-2, 2, 2, 0.0, c.w)
        @test Rup.reflection_amplitude / Rup.incidence_amplitude ≈ c.cref / c.cup rtol=1e-10
        route = D.radial_route(-2, 2, 2, 0.0, c.w, :UP)
        @test Rup.reflection_amplitude / Rup.transmission_amplitude ≈
            route.teukolsky_reflection / route.teukolsky_transmission rtol = 1e-14
    end
end

# ---------------------------------------------------------------------------------------------------
# Method selection: "auto" == strict "ISEM" for default calls; the legacy ODE routes stay available
# and agree at easy points.
# ---------------------------------------------------------------------------------------------------
@testset "method selection and legacy agreement" begin
    for (s, l, m, a, w) in ((-2, 2, 2, 0.5, 0.5), (0, 2, 1, 0.3, 0.8), (2, 3, -2, 0.7, 0.4 - 0.1im), (-1, 2, 0, 0.9, 1.3))
        for bc in (IN, UP)
            auto = GSN_radial(s, l, m, a, w, bc)
            strict = GSN_radial(s, l, m, a, w, bc; method="ISEM")
            @test auto.method == strict.method == "ISEM"
            @test auto.incidence_amplitude == strict.incidence_amplitude
            @test auto.reflection_amplitude == strict.reflection_amplitude
            scale = max(abs(auto.incidence_amplitude), abs(auto.reflection_amplitude)) / abs(auto.transmission_amplitude)
            # Legacy ODE routes only at real omega (their complex-path accuracy is not a tested contract).
            legacy_methods = isreal(w) ? ("linear", "Riccati") : ()
            for method in legacy_methods
                legacy = quietly(() -> GSN_radial(s, l, m, a, w, bc; method=method))
                @test legacy.method == method
                ratio_auto = auto.incidence_amplitude / auto.transmission_amplitude
                ratio_legacy = legacy.incidence_amplitude / legacy.transmission_amplitude
                @test abs(ratio_auto - ratio_legacy) <= 1e-6 * scale
            end
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# The full-state MST backend on the qnm(.68) stencil is rejected, i.e. it fails closed.
# complex_nia_jump is not verified (no independent reference).
# ---------------------------------------------------------------------------------------------------
@testset "fail-closed and unverified paths" begin
    w = 0.5239751042900836 - 0.08151262363120163im
    @test_throws Exception D.radial_route(-2, 2, 2, 0.68, w, :IN; backend=:mst)
    @test_skip false   # needs independent a = 0 and a = 0.7 references at -0.3i for the NIA jump
    for w in (1e-6, 1e-5), bc in (IN, UP)
        P = Teukolsky_radial(-2, 2, 2, 0.7, w, bc)
        N = Teukolsky_radial(-2, 2, -2, 0.7, -w, bc)
        @test all(isfinite, teukolsky_ratios(P))
        @test teukolsky_ratios(N)[1] ≈ conj(teukolsky_ratios(P)[1]) rtol=1e-12
    end
end

@testset "near-static finite-frequency references" begin
    for module_ in (GSN, GSN.ISEM), w in (-1e-13, 1e-13)
        @test !module_._is_static_frequency(w)
    end
    for w in (-1e-13, 1e-13)
        @test !GSN._is_horizon_superradiance_frequency(.5, 0, w)
    end
    for a in (.5, .9), m in (1, 2, 3)
        w = m * GSN.Kerr.omega_horizon(a)
        @test GSN._is_horizon_superradiance_frequency(a, m, w)
        @test !GSN._is_horizon_superradiance_frequency(a, m, nextfloat(w))
        @test !GSN._is_horizon_superradiance_frequency(a, m, prevfloat(w))
    end
    rows = split.(readlines(data_file("near_static.tsv"))[2:end], '\t')
    ultra(row) = startswith(row[1], "EM_") || startswith(row[1], "M0_")
    selected = level_subset(filter(row -> !ultra(row), rows), 24)
    append!(selected, filter(ultra, rows))
    for row in selected
        s, l, m = parse.(Int, row[2:4])
        a, w = parse.(Float64, row[5:6])
        bc = row[7] == "IN" ? IN : UP
        reference = complex.(parse.(Float64, row[8:2:14]), parse.(Float64, row[9:2:15]))
        R = Teukolsky_radial(s, l, m, a, w, bc)
        inc, refl = teukolsky_ratios(R)
        # A shared amplitude scale does not certify a suppressed reflection coefficient.
        scale = max(abs(reference[1]), abs(reference[2]))
        @test scaled_err(inc, reference[1], scale) <= 1e-10
        @test scaled_err(refl, reference[2], scale) <= 1e-10
        @test R(5.0) / R.transmission_amplitude ≈ reference[3] rtol=5e-12
        @test R(10.0) / R.transmission_amplitude ≈ reference[4] rtol=5e-12
    end
    # Independent Teukolsky ODE integrations at 40/65 and 45/70 decimal digits.
    for (r, w, value, derivative) in (
        (1e6, 1.999999999e-9, 9.113446900378605e23+4.1163399025270184e23im,
            3.644832859220209e18+1.647753514567056e18im),
        (2e7, 2.23606797748729e-11, 1.4588376487548126e29+6.571085527056393e28im,
            2.917577459427722e22+1.314434673055514e22im),
        (1.6e8, 1e-12, 5.975915229095591e32+2.690373189724911e32im,
            1.493960878980218e25+6.726331430625239e24im))
        R = Teukolsky_radial(-2, 2, 2, .5, w, IN)
        @test R.mode.omega == w
        @test R(r) ≈ value rtol=5e-12
        @test R.Teukolsky_solution(r)[2] ≈ derivative rtol=5e-12
    end
    R = Teukolsky_radial(-2, 2, 2, .5, 1e-12, IN)
    @test R.incidence_amplitude ≈ -1.5390791994649324e60+3.419610389751223e60im rtol=5e-12
    R = Teukolsky_radial(2, 11, -1, .5, -2.6e-10, UP)
    @test R(7692.) ≈ 2.680381593764424e35-2.5626945604972966e40im rtol=5e-12
    @test R.incidence_amplitude ≈ 2.382837501785559e100+1.854735148541830e100im rtol=5e-12
    for w in (-1e-13, 1e-13)
        R = Teukolsky_radial(-2, 2, 2, .5, w, IN)
        @test R.mode.omega == w
        @test R(10.) ≈ Teukolsky_radial(-2, 2, 2, .5, 0., IN)(10.) rtol=1e-11
    end
    for w in (-1e-13, 1e-13), s in (-2, -1, 0, 1, 2)
        X = GSN_radial(s, 2, 0, .5, w, UP)
        R = Teukolsky_radial(s, 2, 0, .5, w, UP)
        @test X.mode.omega == R.mode.omega == w
        @test X.method == "ISEM"
        @test all(isfinite, (R(5.), R.incidence_amplitude, X.incidence_amplitude))
        if s == 2
            Y = Y_radial(s, 2, 0, .5, w, UP)
            @test Y.mode.omega == w
            @test isfinite(Y(5.))
        end
    end
end

@testset "high-l near-zone UP boundary solution" begin
    # Direct Teukolsky integration at 60/90 decimal digits, with different
    # outer radii and Taylor step sizes; not an MST or static reference.
    for (l,m,w,r,value,derivative) in (
        (16,16,5.679891948853166e-6,3e4,
            -5939.007601126625-127.28195072097409im,
            3.761404827668372+.07638759120253598im),
        (17,16,5.679891948853166e-6,3e4,
            20768.486995678705-1.0291431344800493e6im,
            -13.15676825250353+686.1024241020361im),
        (18,16,5.679891948853166e-6,3e4,
            1.902560369736828e8+3.6278670806021436e6im,
            -133180.7617156858-2419.23827802867im),
        (23,23,1.8576810288076106e-6,1e3,
            2.887454640048551e65+3.9011845440541655e68im,
            -7.952530760069528e63-1.015323585022585e67im),
        (30,30,.00044046087194570704,1e3,
            1.1374023082001997e33+4.262780137633929e31im,
            -3.756941824312484e31-1.373480525496288e30im))
        for sign in (-1,1)
            R=Teukolsky_radial(2,l,sign*m,.5,sign*w,UP)
            expected=sign==1 ? (value,derivative) : conj.((value,derivative))
            state=R.Teukolsky_solution(r)
            @test state[1] ≈ expected[1] rtol=5e-12
            @test state[2] ≈ expected[2] rtol=5e-12
            @test R(r) ≈ expected[1] rtol=5e-12
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# High-spin lower-half-plane IN fallback (main module _use_lhp_highspin_p_fallback): Re omega <= -10,
# |Im omega| <= 1e-3 |Re omega|, |a| >= 0.99 uses the conjugated p-equation route of the positive partner.
# ---------------------------------------------------------------------------------------------------
@testset "high-spin LHP IN fallback" begin
    for (a, w) in ((0.99, -12.0 - 0.001im), (0.995, -15.0 - 0.01im))
        X = GSN_radial(-2, 2, 2, a, w, IN)
        @test X.method == "ISEM"
        @test all(isfinite, (X.incidence_amplitude, X.reflection_amplitude, X.transmission_amplitude))
        P = GSN_radial(-2, 2, -2, a, -conj(w), IN)
        scale = max(abs(P.incidence_amplitude), abs(P.reflection_amplitude)) / abs(P.transmission_amplitude)
        @test abs(X.incidence_amplitude / X.transmission_amplitude -
            conj(P.incidence_amplitude / P.transmission_amplitude)) <= 1e-6 * scale
    end
end

@testset "real-axis normalization reference" begin
    p = only(p for p in GRID_POINTS if p.id == 70)
    reference = only(p for p in GRID_REFERENCES if p.id == 70 && p.branch == "IN")
    solution = Teukolsky_radial(p.s,p.l,p.m,p.a,p.omega,IN)
    for (r,value) in ((5.0,reference.R5),(10.0,reference.R10))
        @test solution(r)/solution.transmission_amplitude ≈ value rtol=2e-12
    end
    route = GSN.ISEM.radial_route(p.s,p.l,p.m,p.a,p.omega,:IN)
    @test route.metadata.real_axis_state_anchor_status == :uncertified_reference
end

end # isem
