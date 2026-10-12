# Homogeneous radial solutions: GSN_radial / Teukolsky_radial against package-independent references,
# plus invariants on parameter grids. Reference provenance: ../README.md, "Radial amplitudes".
#
# Conventions used below (package source, verified against the references in this file):
# * Teukolsky_radial(..., IN/UP) is normalized to unit Teukolsky transmission (UNIT_TEUKOLSKY_TRANS);
#   incidence/reflection are B_inc, B_ref (IN) and C_up, C_ref (UP), as ratios to the transmission
#   amplitude.
# * Exact conjugation symmetry of the Teukolsky equation on real r (src/ISEM/ISEM.jl
#   _conjugate_teukolsky_radial_function; DirectComplexFrequency.direct_gsn_radial comment):
#       R_{s,l,m}(omega; r) = conj R_{s,l,-m}(-conj(omega); r),  lambda -> conj(lambda),
#   so the amplitude ratios map to their complex conjugates.
# * Spin reflection (src/GeneralizedSasakiNakamura.jl, `a < zero(a)` branches):
#       Teukolsky_radial(s,l,m,-a,omega) returns the (s,l,-m,a,omega) solution with mode relabelled,
#   so amplitudes are identical.

if !@isdefined(GSN)
    include(joinpath(@__DIR__, "..", "load_package.jl"))
end
include(joinpath(@__DIR__, "..", "common", "helpers.jl"))
include(data_file("radial_references.jl"))

const D = GSN.ISEM.DirectGSN
const CF = GSN.ConversionFactors

dominant(x, y) = max(abs(x), abs(y))
branch_symbol(bc) = bc == IN ? :IN : :UP

"Teukolsky amplitude conversion (incidence factor, reflection factor) over the transmission factor."
function conversion_ratios(s, m, a, omega, lambda, bc)
    if bc == IN
        t = CF.Btrans(s, m, a, omega, lambda)
        return CF.Binc(s, m, a, omega, lambda) / t, CF.Bref(s, m, a, omega, lambda) / t
    else
        t = CF.Ctrans(s, m, a, omega, lambda)
        return CF.Cinc(s, m, a, omega, lambda) / t, CF.Cref(s, m, a, omega, lambda) / t
    end
end

const SUPERRADIANT_TOL = 1e-12
is_threshold(a, m, omega) = isreal(omega) && abs(omega - m * a / (2 * kerr_rplus(a))) < SUPERRADIANT_TOL

@testset "homogeneous" begin

# ---------------------------------------------------------------------------------------------------
# 1. Independent references: G1-G10 (strong damping), 24 QNM-root/off-root Float64 points, 39 real and
#    weakly damped validation points (s = -2). 4 amplitudes, lambda and the amplitude Wronskian per point.
#    Tolerances: dominant units, 1e-9 unless the row records a wider class (data/radial_references.jl).
# ---------------------------------------------------------------------------------------------------
@testset "independent references: public Teukolsky_radial" begin
    for p in level_subset(STRONG_DAMPING_POINTS, 6)
        @testset "$(p.id)" begin
            Rin = Teukolsky_radial(p.s, p.l, p.m, p.a, p.omega, IN)
            Rup = Teukolsky_radial(p.s, p.l, p.m, p.a, p.omega, UP)
            values = (teukolsky_ratios(Rin)..., teukolsky_ratios(Rup)...)
            refs = (p.binc, p.bref, p.cup, p.cref)
            scales = (dominant(p.binc, p.bref), dominant(p.binc, p.bref),
                dominant(p.cup, p.cref), dominant(p.cup, p.cref))
            for j in 1:4
                @test all(isfinite, (values[j],))
                @test scaled_err(values[j], refs[j], scales[j]) <= p.tol[j]
                # Non-suppressed amplitudes must also be right relative to themselves.
                if abs(refs[j]) >= 0.1 * scales[j]
                    @test relerr(values[j], refs[j]) <= 10 * p.tol[j]
                end
            end
            @test Rin.normalization_convention == UNIT_TEUKOLSKY_TRANS
            @test Rup.normalization_convention == UNIT_TEUKOLSKY_TRANS
            @test Rin.boundary_condition == IN && Rup.boundary_condition == UP
            # Angular eigenvalue: independent Leaver/spectral value (input lambda error is a few ulp).
            lam_ref = p.lam !== nothing ? p.lam :
                p.E - 2 * p.a * p.m * p.omega + p.a^2 * p.omega^2 - p.s * (p.s + 1)
            @test abs(Rin.mode.lambda - lam_ref) <= 1e-12 * max(1, abs(lam_ref))
            @test isapprox(Rin.mode.lambda, Rup.mode.lambda; rtol=1e-14)
            # Amplitude Wronskian C_up (4ikr+ + 2s kappa) = 2i omega B_inc: first on the references (checks the
            # convention used by these tests), then on the package values.
            if abs(p.binc) > 1e-6 * scales[1] && abs(p.cup) > 1e-6 * scales[3]
                lr, rr = amplitude_wronskian_sides(p.s, p.m, p.a, p.omega, p.binc, p.cup)
                @test abs(lr - rr) <= 1e-12 * max(abs(lr), abs(rr))
                lv, rv = amplitude_wronskian_sides(p.s, p.m, p.a, p.omega, values[1], values[3])
                @test abs(lv - rv) <= 1e-8 * max(abs(lv), abs(rv))
            end
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# 2. The same points through the direct route, with its own error report (direct_amplitude_errors).
#    Coverage rule (from measured errors at these points):
#    * complex omega, non-suppressed amplitude: actual relative error <= 2 x reported
#      (known under-reports x1.12-1.32);
#    * real omega: <= 100 x reported (measured under-reports up to x51 at omega = 1);
#    * incidence at a QNM root (`_at` points, value ~0): reported >= 0.5, i.e. "no digits" is reported.
# ---------------------------------------------------------------------------------------------------
@testset "independent references: direct route and reported errors" begin
    for p in level_subset(STRONG_DAMPING_POINTS, 6), (bi, branch) in ((1, :IN), (2, :UP))
        @testset "$(p.id) $(branch)" begin
            route = D.direct_gsn_radial(p.s, p.l, p.m, p.a, p.omega, branch)
            R = D.direct_teukolsky_radial_function(route)
            values = teukolsky_ratios(R)
            errors = D.DirectMatching.direct_amplitude_errors(route)
            refs = branch == :IN ? (p.binc, p.bref) : (p.cup, p.cref)
            scale = dominant(refs...)
            @test all(isfinite, errors)   # no Inf reports
            for j in 1:2
                jj = 2 * (bi - 1) + j
                @test scaled_err(values[j], refs[j], scale) <= p.tol[jj]
                actual = relerr(values[j], refs[j])
                suppressed = abs(refs[j]) < 1e-6 * scale
                if endswith(p.id, "_at_f64") && j == 1
                    @test errors[j] >= 0.5
                elseif !suppressed
                    factor = isreal(p.omega) ? 100.0 : 2.0
                    @test actual <= factor * errors[j] + 1e-14
                end
            end
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# 3. Independent (id, branch) Teukolsky references on the s = -2..2 grid: amplitudes in dominant units and
#    R(5), R(10) relative. R tolerance class from the measured package error (meas_R5, meas_R10):
#    <= 1e-8 -> 1e-6;  <= 1e-5 -> 1e-4 (R(r) floor at real/weak damping);
#    larger -> @test_broken (R(r) past the dominant/subdominant crossover at strong damping).
# ---------------------------------------------------------------------------------------------------
GRID_BY_ID = Dict(g.id => g for g in GRID_POINTS)
# Amplitude cells with measured errors above 2e-10 (strong damping, a >= 0.9): wider tolerances.
GRID_AMPLITUDE_EXCEPTIONS = Dict((142, "IN") => 2e-8, (178, "UP") => 2e-8, (159, "UP") => 2e-8,
    (72, "IN") => 5e-9, (72, "UP") => 5e-9, (157, "UP") => 5e-9)

@testset "independent grid references (896-bit)" begin
    for ref in level_subset(GRID_REFERENCES, 7)
        g = GRID_BY_ID[ref.id]
        bc = ref.branch == "IN" ? IN : UP
        @testset "grid$(ref.id) $(ref.branch) (s=$(g.s), l=$(g.l), m=$(g.m), a=$(g.a), omega=$(g.omega))" begin
            R = Teukolsky_radial(g.s, g.l, g.m, g.a, g.omega, bc)
            inc, refl = teukolsky_ratios(R)
            scale = dominant(ref.incidence, ref.reflection)
            tol = get(GRID_AMPLITUDE_EXCEPTIONS, (ref.id, ref.branch), 1e-9)
            @test scaled_err(inc, ref.incidence, scale) <= tol
            @test scaled_err(refl, ref.reflection, scale) <= tol
            for (rr, refR, meas) in ((5.0, ref.R5, ref.meas_R5), (10.0, ref.R10, ref.meas_R10))
                value = R(rr) / R.transmission_amplitude
                err = relerr(value, refR)
                if isnan(meas) && abs(imag(g.omega)) >= 1
                    # Not measured, and in the strong-damping crossover regime: no class to assert.
                    @test_skip err <= 1e-6
                elseif isnan(meas) || meas <= 1e-8
                    @test err <= 1e-6
                elseif meas <= 1e-5
                    @test err <= 1e-4
                else
                    @test_broken err <= 1e-6   # strong-damping radial crossover
                end
            end
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# 4. Real-axis cells at strong damping: independent lambda (Leaver CF + spectral, dps 160),
#    B_inc / B_ref, and R(r)/B_trans on r in {1.5..20}. R classes as in 3.
# ---------------------------------------------------------------------------------------------------
@testset "independent real-axis cells (cloud)" begin
    for c in level_subset(CLOUD_CELLS, 4)
        @testset "cell$(c.id) $(c.branch)" begin
            bc = c.branch == "IN" ? IN : UP
            R = Teukolsky_radial(c.s, c.l, c.m, c.a, c.omega, bc)
            @test abs(R.mode.lambda - c.lambda) <= 1e-11 * max(1, abs(c.lambda))
            inc, refl = teukolsky_ratios(R)
            scale = dominant(c.binc, c.bref)
            @test scaled_err(inc, c.binc, scale) <= 1e-9
            @test scaled_err(refl, c.bref, scale) <= 1e-9
            for point in c.R
                r_plus = kerr_rplus(c.a)
                point.r > r_plus || continue
                err = relerr(R(point.r) / R.transmission_amplitude, point.R)
                if isnan(point.meas)
                    @test_skip err <= 1e-6   # not measured (strong-damping crossover regime)
                elseif point.meas <= 1e-8
                    @test err <= 1e-6
                elseif point.meas <= 1e-5
                    @test err <= 1e-4
                else
                    @test_broken err <= 1e-6   # strong-damping radial crossover
                end
            end
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# 5. Invariants on the 180-point review grid (s = -2..2, l = 2..4, a in {0, 1e-4, .3, .7, .9, .99},
#    real omega 0.05..3, complex omega down to Im = -5, both signs of Re omega).
# ---------------------------------------------------------------------------------------------------
function check_radial_invariants(s, l, m, a, omega; wronskian_radii=(4.0, 6.0, 10.0, 20.0))
    threshold = is_threshold(a, m, omega)
    Rin = Teukolsky_radial(s, l, m, a, omega, IN)
    Rup = Teukolsky_radial(s, l, m, a, omega, UP)
    for R in (Rin, Rup)
        @test R.normalization_convention == UNIT_TEUKOLSKY_TRANS
        @test isapprox(R.transmission_amplitude, 1; atol=1e-12)
        @test all(isfinite, (R.incidence_amplitude, R.reflection_amplitude))
        @test isapprox(R.mode.lambda, Rin.mode.lambda; rtol=1e-14)
    end
    binc, bref = teukolsky_ratios(Rin)
    cup, cref = teukolsky_ratios(Rup)
    # GSN <-> Teukolsky amplitude conversion (ConversionFactors) on the separately built GSN function.
    for (R, bc) in ((Rin, IN), (Rup, UP))
        X = GSN_radial(s, l, m, a, omega, bc)
        @test X.method == "GSN-ISEM"
        @test X.normalization_convention == UNIT_GSN_TRANS
        @test X.boundary_condition == bc
        @test isapprox(X.mode.lambda, R.mode.lambda; rtol=1e-14)
        fi, fr = conversion_ratios(s, m, a, omega, X.mode.lambda, bc)
        ti, tr = teukolsky_ratios(R)
        scale = dominant(ti, tr)
        @test abs(fi * X.incidence_amplitude / X.transmission_amplitude - ti) <= 1e-10 * scale
        @test abs(fr * X.reflection_amplitude / X.transmission_amplitude - tr) <= 1e-10 * scale
        for rs in (-20.0, 0.0, 30.0)
            @test isfinite(X(rs))
        end
    end
    # Amplitude Wronskian identity (see helpers.jl).
    if !threshold && abs(binc) > 1e-6 * dominant(binc, bref) && abs(cup) > 1e-6 * dominant(cup, cref)
        lhs, rhs = amplitude_wronskian_sides(s, m, a, omega, binc, cup)
        @test abs(lhs - rhs) <= 1e-8 * max(abs(lhs), abs(rhs))
    end
    # Radial Wronskian Delta^(s+1)(Rin Rup' - Rup Rin') = 2 i omega B_inc (unit transmissions), constant in r.
    # Only real or weakly damped omega: R(r) at strong damping is not certified. Tolerance 1e-5: R(r) floor
    # up to 4.5e-7 at real omega (grid id7).
    if !threshold && abs(imag(omega)) <= 0.1
        expected = 2im * omega * binc
        for r in wronskian_radii
            W = radial_wronskian(Rin, Rup, s, a, r) / (Rin.transmission_amplitude * Rup.transmission_amplitude)
            @test abs(W - expected) <= 1e-5 * abs(expected)
        end
    end
    # Conjugation symmetry (m, omega) -> (-m, -conj(omega)): exact for default calls.
    if !threshold && !is_threshold(a, -m, -conj(omega))
        for (R, bc) in ((Rin, IN), (Rup, UP))
            P = Teukolsky_radial(s, l, -m, a, -conj(omega), bc)
            pi_, pr = teukolsky_ratios(P)
            ti, tr = teukolsky_ratios(R)
            scale = dominant(ti, tr)
            # Exact for the direct route; 1e-11 allows the separately evaluated ConversionFactors.
            @test abs(pi_ - conj(ti)) <= 1e-11 * scale
            @test abs(pr - conj(tr)) <= 1e-11 * scale
            @test isapprox(P.mode.lambda, conj(R.mode.lambda); rtol=1e-13)
            if abs(imag(omega)) <= 0.1
                v1 = R(6.0) / R.transmission_amplitude
                v2 = P(6.0) / P.transmission_amplitude
                @test abs(v2 - conj(v1)) <= 1e-10 * abs(v1)
            end
        end
    end
    # Spin reflection (a, m) -> (-a, -m): the same solution, relabelled.
    if a > 0
        for (R, bc) in ((Rin, IN), (Rup, UP))
            N = Teukolsky_radial(s, l, -m, -a, omega, bc)
            @test N.mode.a == -a && N.mode.m == -m
            @test N.incidence_amplitude == R.incidence_amplitude
            @test N.reflection_amplitude == R.reflection_amplitude
        end
    end
    return nothing
end

@testset "invariants on the 180-point grid" begin
    for g in level_subset(GRID_POINTS, 9)
        @testset "grid$(g.id) (s=$(g.s), l=$(g.l), m=$(g.m), a=$(g.a), omega=$(g.omega))" begin
            check_radial_invariants(g.s, g.l, g.m, g.a, g.omega)
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# 6. Cartesian sweep (full level): s = -2..2 x (l, m) x a x omega, invariants only.
# ---------------------------------------------------------------------------------------------------
SWEEP_LM = ((2, 2), (2, -1), (3, 0), (3, 3), (4, -2))
SWEEP_A = (0.0, 0.5, 0.9, 0.99)
SWEEP_OMEGA = (0.1, 0.8, 2.0, 0.6 - 0.2im, 1.5 - 0.8im, -0.4 - 0.3im)

@testset "invariant sweep" begin
    # 600 cases; ~10 radial constructions each. quick: every 29th, full: every 4th, extended: all.
    cases = [(s, l, m, a, w) for s in -2:2 for (l, m) in SWEEP_LM if l >= abs(s)
             for a in SWEEP_A for w in SWEEP_OMEGA]
    stride = GSN_TEST_LEVEL === :quick ? 29 : GSN_TEST_LEVEL === :full ? 4 : 1
    cases = cases[1:stride:end]
    for (s, l, m, a, w) in cases
        @testset "sweep (s=$s, l=$l, m=$m, a=$a, omega=$w)" begin
            check_radial_invariants(s, l, m, a, w; wronskian_radii=(5.0, 12.0))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# 7. Near-horizon normalization of X_in: X_in(r*) e^{i p r*} -> X_trans as r* -> -inf, p = omega - m Omega_H.
#    At r* = -40, r - r_+ ~ exp(-40 (r_+ - r_-) / (2 r_+)) <= 6e-8 for a <= 0.7; tolerance 1e-4.
# ---------------------------------------------------------------------------------------------------
@testset "IN normalization at the horizon" begin
    for s in (-2, 0, 2), (l, m) in ((2, 2), (3, -1)), a in (0.0, 0.7), w in (0.3, 1.2)
        l >= abs(s) || continue
        X = GSN_radial(s, l, m, a, w, IN)
        p = w - m * a / (2 * kerr_rplus(a))
        rs = -40.0
        @test isapprox(X(rs) * exp(1im * p * rs), X.transmission_amplitude; rtol=1e-4)
    end
end

# ---------------------------------------------------------------------------------------------------
# 8. OUT / DOWN combinations (main-module combiners): amplitudes from the IN/UP formulas of
#    _combine_direct_gsn_down/out and finite solutions.
# ---------------------------------------------------------------------------------------------------
@testset "OUT and DOWN boundary conditions" begin
    for (s, l, m, a, w) in ((-2, 2, 2, 0.7, 0.5), (2, 3, -1, 0.3, 1.1 - 0.1im), (0, 2, 0, 0.9, 0.25))
        Xin = GSN_radial(s, l, m, a, w, IN)
        Xup = GSN_radial(s, l, m, a, w, UP)
        Xdown = GSN_radial(s, l, m, a, w, DOWN)
        Xout = GSN_radial(s, l, m, a, w, OUT)
        Bt, Bi, Br = Xin.transmission_amplitude, Xin.incidence_amplitude, Xin.reflection_amplitude
        Ct, Ci, Cr = Xup.transmission_amplitude, Xup.incidence_amplitude, Xup.reflection_amplitude
        @test Xdown.boundary_condition == DOWN && Xout.boundary_condition == OUT
        @test Xdown.transmission_amplitude == 1 && Xout.transmission_amplitude == 1
        @test isapprox(Xdown.incidence_amplitude, Bt / Bi - Br * Cr / (Bi * Ct); rtol=1e-12)
        @test isapprox(Xdown.reflection_amplitude, -Br * Ci / (Bi * Ct); rtol=1e-12)
        @test isapprox(Xout.incidence_amplitude, Ct / Ci - Br * Cr / (Bt * Ci); rtol=1e-12)
        @test isapprox(Xout.reflection_amplitude, -Bi * Cr / (Bt * Ci); rtol=1e-12)
        for rs in (-10.0, 5.0, 40.0)
            @test isfinite(Xdown(rs)) && isfinite(Xout(rs))
        end
        Rdown = Teukolsky_radial(s, l, m, a, w, DOWN)
        Rout = Teukolsky_radial(s, l, m, a, w, OUT)
        @test Rdown.boundary_condition == DOWN && Rout.boundary_condition == OUT
        @test isfinite(Rdown(8.0)) && isfinite(Rout(8.0))
    end
end

# ---------------------------------------------------------------------------------------------------
# 8b. Superradiant threshold omega = m Omega_H (non-extremal): UP takes the original-ISEM threshold builder
#     (main module `_is_horizon_superradiance_frequency`); results must be finite and normalized.
# ---------------------------------------------------------------------------------------------------
@testset "superradiant threshold" begin
    for (s, l, m, a) in ((-2, 2, 2, 0.7), (-2, 3, 1, 0.9), (2, 2, 2, 0.5))
        w = m * a / (2 * kerr_rplus(a))
        for bc in (IN, UP)
            R = Teukolsky_radial(s, l, m, a, w, bc)
            @test R.normalization_convention == UNIT_TEUKOLSKY_TRANS
            # The threshold UP builder carries no incidence/reflection amplitude (missing); IN carries all three.
            amps = (R.transmission_amplitude, R.incidence_amplitude, R.reflection_amplitude)
            @test all(v -> (bc == UP && v === missing) || isfinite(v), amps)
            @test isfinite(R(8.0))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# 9. Exactly static frequency: exact hypergeometric solutions, method "static".
# ---------------------------------------------------------------------------------------------------
@testset "static frequency" begin
    for s in (-2, -1, 0, 1, 2), (l, m) in ((2, 0), (2, 2), (3, -1)), a in (0.0, 0.5, 0.9)
        l >= abs(s) || continue
        for bc in (IN, UP), w in (0.0, -0.0)
            R = Teukolsky_radial(s, l, m, a, w, bc)
            @test R.mode.omega == 0
            @test R.mode.lambda == GSN.SpinWeightedSpheroidalHarmonics.spin_weighted_spherical_eigenvalue(s, l, m)
            X = GSN_radial(s, l, m, a, w, bc)
            @test X.method == "static"
            @test isfinite(R(6.0))
            @test isfinite(X(5.0))
        end
    end
    # Terminating regularized hypergeometric polynomials, including their derivatives.
    for (s, l, coefficient, power) in ((1, 2, 1.5, 3), (2, 3, 3.75, 5)),
        a in (0.0, 0.5, 0.9), r in (3.0, 6.0, 20.0)
        gamma = sqrt(1-a^2)
        R = Teukolsky_radial(s, l, 0, a, 0.0, IN)
        value, derivative = R.Teukolsky_solution(r)
        @test value ≈ coefficient*(r-1)/gamma^power rtol=2e-13
        @test derivative ≈ coefficient/gamma^power rtol=2e-13
    end
end

# ---------------------------------------------------------------------------------------------------
# 10. Exact-extremal Kerr (a = +-1): finite GSN-ISEM solutions; synchronous omega = m Omega_H = m/2 is a
#     branch point and must raise DomainError (src/Homogeneous/ExtremalRadial.jl _require_nonsynchronous).
# ---------------------------------------------------------------------------------------------------
@testset "exact extremal" begin
    X = GSN_radial(-2, 2, 2, 1.0, 0.3, IN, -40.0, 60.0; method="GSN-ISEM", tolerance=1.0e-10)
    @test X.method == "GSN-ISEM"
    @test all(isfinite, X.GSN_solution(0.0))
    for (s, l, m, w) in ((-2, 2, 2, 0.3), (-2, 2, 1, 0.8), (0, 2, 0, 0.5), (2, 3, 3, 1.0))
        for bc in (IN, UP)
            Xe = GSN_radial(s, l, m, 1.0, w, bc)
            @test Xe.method == "GSN-ISEM"
            @test all(isfinite, (Xe.transmission_amplitude, Xe.incidence_amplitude))
            @test isfinite(Xe(0.0))
            Re = Teukolsky_radial(s, l, m, 1.0, w, bc)
            @test isfinite(Re(5.0))
        end
    end
    for (m, w) in ((1, 0.5), (2, 1.0), (3, 1.5))
        @test_throws DomainError GSN_radial(-2, 3, m, 1.0, w, IN; method="GSN-ISEM")
        @test_throws DomainError GSN_radial(-2, 3, m, 1.0, w, UP)
    end
    @test r_from_rstar(-1.0, rstar_from_r(-1.0, 2.0)) ≈ 2.0
    @test r_from_rstar(1.0, rstar_from_r(1.0, 3.5)) ≈ 3.5
end

# ---------------------------------------------------------------------------------------------------
# 11. Coordinates: r_from_rstar inverts rstar_from_r; r* is increasing in r.
# ---------------------------------------------------------------------------------------------------
@testset "tortoise coordinate" begin
    for a in (0.0, 0.3, 0.7, 0.9, 0.99, 0.999), r in (1.0001, 1.01, 1.5, 2.0, 3.0, 6.0, 20.0, 1e3, 1e5)
        rp = kerr_rplus(a)
        rr = rp * r + (r > 10 ? r : 0.0)   # points outside the horizon
        rr <= rp && continue
        rs = rstar_from_r(a, rr)
        @test isapprox(r_from_rstar(a, rs), rr; rtol=1e-10)
        @test rstar_from_r(a, rr * 1.001) > rs
    end
end

end # homogeneous
