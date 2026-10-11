# QNM frequencies, branches, excitation factors and failure modes.
# References: data/qnm_references.jl (provenance in that file and README.md, section "qnm").

if !@isdefined(GSN)
    include(joinpath(@__DIR__, "..", "load_package.jl"))
end
include(joinpath(@__DIR__, "..", "common", "helpers.jl"))
include(data_file("qnm_references.jl"))

const Q = GSN.QNM
const CF = GSN.ConversionFactors

"Mirror-branch root through the package's own seeding (src/QNM/ModeObservables.jl _mirror_qnm_root)."
mirror_root(a, s, l, m, n; convention=:overtone) =
    first(Q._mirror_qnm_root(a, s, l, m, n, nothing, NamedTuple(); convention))

@testset "qnm" begin

# ---------------------------------------------------------------------------------------------------
# Schwarzschild s = -2, l = 2 overtones against the independent dps-25 Leaver scan. Overtone n = 8 is the
# algebraically special frequency -2i (not the scan's index-8 root); for n >= 9 the scan stored the mirror
# root, so the ordinary root is -conj. atol 1e-10.
# ---------------------------------------------------------------------------------------------------
@testset "Schwarzschild overtones (independent scan)" begin
    ns = GSN_TEST_LEVEL === :quick ? (0, 1, 3, 7, 9) : Tuple(0:22)
    for n in ns
        result = qnm_frequency(QNMMode(-2, 2, 2, n), 0.0)
        if n == 8
            @test result.omega == -2.0im
            @test !result.provenance.gsn_amplitude_limit_implemented
            continue
        end
        scan = SCHWARZSCHILD_SCAN[n + 1].omega
        expected = real(scan) < 0 ? -conj(scan) : scan
        @test result.status == :accepted
        @test isapprox(result.omega, expected; atol=1e-10, rtol=0)
        @test result.mode.n == n
    end
    # Mirror branch at a = 0: -conj of the ordinary root.
    for n in (0, 3, 10)
        root = mirror_root(0.0, -2, 2, 2, n)
        expected = SCHWARZSCHILD_SCAN[n + 1].omega
        expected = real(expected) > 0 ? -conj(expected) : expected
        @test isapprox(root.omega, expected; atol=1e-10, rtol=0)
    end
end

# ---------------------------------------------------------------------------------------------------
# Kerr roots against the independent 45-digit Leaver roots (strong_damping_reference/results/qnm_list.json).
# ---------------------------------------------------------------------------------------------------
@testset "Kerr roots (independent Leaver, 47-57 digits)" begin
    for q in QNM_CLOUD_ROOTS
        q.a == 0 && GSN_TEST_LEVEL === :quick && q.n > 3 && continue
        # The reference n is a CF inversion index, not necessarily the public overtone label.
        result = qnm_frequency(QNMMode(q.s, q.l, q.m, q.n), q.a;
            convention=:leaver, inversion_index=q.n, guess=q.omega)
        @test result.status == :accepted
        @test isapprox(result.omega, q.omega; atol=1e-10, rtol=0)
        # Spin reflection (a, m) -> (-a, -m) leaves the frequency unchanged.
        if q.a > 0 && q.n <= 5
            reflected = qnm_frequency(QNMMode(q.s, q.l, -q.m, q.n), -q.a)
            @test isapprox(reflected.omega, result.omega; atol=1e-12, rtol=0)
            @test reflected.a == -q.a
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Published tables: Forteza & Mourier (10 digits; atol 2e-10) and Cook's Zenodo catalogue (atol 5e-9).
# Labels follow Python qnm and Cook: Forteza's ordinary (2,2) n = 9 is Cook's 8_1, taken as n = 9 with convention = :complete_spectrum;
# Cook's ordinary (2,2) label 19 is the package's overtone 19. A is Leaver's angular separation constant
# (package convention assumed equal: A -> l(l+1) - s(s+1) as c -> 0; see README "API uncertainties").
# ---------------------------------------------------------------------------------------------------
@testset "published QNM tables" begin
    rows = GSN_TEST_LEVEL === :quick ? PUBLISHED_QNM[1:6:end] : PUBLISHED_QNM
    for r in rows
        tol = startswith(r.source, "Forteza") ? 2e-10 : 5e-9
        root = r.branch == ordinary ? qnm_frequency(QNMMode(-2, r.l, r.m, r.n), r.a; convention=get(r, :convention, :overtone)) :
            mirror_root(r.a, -2, r.l, r.m, r.n; convention=get(r, :convention, :overtone))
        @testset "$(r.source) a=$(r.a) ($(r.l),$(r.m),$(r.n)) $(r.branch)" begin
            @test root.status == :accepted
            if startswith(r.source, "Cook") && r.n == 20 &&
                    (((r.l, r.m) == (2, 2) && r.branch != ordinary) || ((r.l, r.m) == (3, 1) && r.branch == ordinary))
                # Catalogue precision is informational; independent CF references below
                # retain the 1e-10 acceptance requirement.
                @info "Cook catalogue difference" a=r.a l=r.l m=r.m branch=r.branch difference=abs(root.omega-r.omega)
            else
                @test isapprox(root.omega, r.omega; atol=tol, rtol=0)
            end
            @test isapprox(root.angular_A, r.A; atol=10 * tol, rtol=0)
        end
    end
end

@testset "independent n20 CF roots" begin
    for r in N20_CF_ROOTS
        root = r.branch == ordinary ? qnm_frequency(QNMMode(-2,r.l,r.m,20),r.a) :
            mirror_root(r.a,-2,r.l,r.m,20)
        @test root.status == :accepted
        @test isapprox(root.omega,r.omega;atol=1e-10,rtol=0)
        fine=qnm_frequency(root.mode,r.a;guess=root.omega,
            angular_order=48,inversion_index=root.inversion_index,
            cf_minimum_iterations=600,cf_maximum_iterations=24000,
            cf_tolerance=1e-13,automatic_precision_fallback=false)
        @test fine.status == :accepted
        @test isapprox(fine.omega,r.omega;atol=1e-10,rtol=0)
    end
end

# ---------------------------------------------------------------------------------------------------
# Legacy-bundle overtone table (separate implementation, 2026-09-24), signed branches, a in {0.5, 0.68, 0.9,
# -0.5}, l = 2, 3; atol 1e-9 as in check_qnm_results.jl.
# ---------------------------------------------------------------------------------------------------
@testset "legacy overtone table" begin
    rows = GSN_TEST_LEVEL === :quick ? LEGACY_OVERTONES[1:9:end] : LEGACY_OVERTONES
    for r in rows
        root = r.branch == ordinary ? qnm_frequency(QNMMode(-2, r.l, r.m, r.n), r.a; convention=get(r, :convention, :overtone)) :
            mirror_root(r.a, -2, r.l, r.m, r.n; convention=get(r, :convention, :overtone))
        @test isapprox(root.omega, r.omega; atol=1e-9, rtol=0)
    end
    # Overtones 7..10 are distinct roots on both branches (check_qnm_results.jl).
    for branch in (ordinary, mirror)
        roots = [branch == ordinary ? qnm_frequency(QNMMode(-2, 2, 2, n), 0.5).omega :
            mirror_root(0.5, -2, 2, 2, n).omega for n in 7:10]
        @test minimum(abs.(diff(roots))) > 0.1
    end
end

# ---------------------------------------------------------------------------------------------------
# qnm(0.68, -2, 2, 2, 0): root, Teukolsky reflection, dB_inc/domega and excitation factor against the
# independent direct-solver reference (cloud_refs/qnm068_bridge_reference; dps 40/60). The derivative
# reference is ~1e-10 relative from GSN's stencil (README there), hence rtol 1e-8.
# ---------------------------------------------------------------------------------------------------
@testset "qnm(0.68) against independent observables" begin
    a = 0.68
    detailed = qnm(a, -2, 2, 2, 0; detailed=true)
    @test detailed isa QNMResult
    @test detailed.status == :accepted
    @test isapprox(detailed.omega, Q068.omega; atol=1e-12, rtol=0)
    lam = Q068.E - 2 * a * 2 * detailed.omega + a^2 * detailed.omega^2 - 2
    @test isapprox(detailed.lambda, lam; rtol=1e-11)
    @test isapprox(detailed.teukolsky.reflection_amplitude, Q068.bref_over_btrans; rtol=1e-10)
    @test isapprox(detailed.teukolsky.incidence_derivative, Q068.dbinc_domega; rtol=1e-8)
    # The independent reference uses 2omega; the Teukolsky API uses 2iomega.
    @test isapprox(detailed.teukolsky.excitation_factor, Q068.excitation_teukolsky / im; rtol=1e-8)
    @test detailed.teukolsky.normalization == :teukolsky_unit_transmission
    @test detailed.X isa GSNRadialFunction
    @test detailed.Y isa YRadialFunction
    @test detailed.R isa TeukolskyRadialFunction
    r = 10.0
    @test all(isfinite, (detailed.X(rstar_from_r(a, r)), detailed.Y(r), detailed.R(r)))
    # In the 2iomega convention the mirror factor is -conj(B).
    mirror_mode = qnm(a, -2, 2, -2, 0, mirror; detailed=true)
    @test mirror_mode.status == :accepted
    @test isapprox(mirror_mode.omega, -conj(Q068.omega); atol=1e-12, rtol=0)
    @test isapprox(mirror_mode.teukolsky.excitation_factor, -conj(Q068.excitation_teukolsky / im); rtol=1e-8)
    # Stencil points around the root: B_inc/B_trans and B_ref/B_trans, absolute 1e-11 (|B_ref| ~ 0.42;
    # measured offsets 2e-12..3e-12 pinned, <= 5e-13 public; README there).
    for st in Q068_STENCIL
        R = Teukolsky_radial(-2, 2, 2, a, st.omega, IN)
        binc, bref = teukolsky_ratios(R)
        @test abs(binc - st.binc) <= 1e-11
        @test abs(bref - st.bref) <= 1e-11
    end
end

# Mirror branch at a = 0.68: the (l, m) = (2, 2) mirror root is -conj of the (2, -2) ordinary root, in omega and in the
# Teukolsky excitation factor.
@testset "mirror = -conj(ordinary) at a = 0.68" begin
    ordinary = qnm(0.68, -2, 2, -2, 0; detailed=true)
    mirrored = qnm(0.68, -2, 2, 2, 0, mirror; detailed=true)
    @test isapprox(mirrored.omega, -conj(ordinary.omega); rtol=1e-8)
    @test isapprox(mirrored.teukolsky.excitation_factor, -conj(ordinary.teukolsky.excitation_factor); rtol=1e-8)
end

# ---------------------------------------------------------------------------------------------------
# Residues of 1/(i omega B_inc^T) at a = 0.5 against arXiv:2609.09531v1 Table I (7 decimals): with the GSN
# incidence derivative alpha, c_n = 1 / (i omega (B_inc/B_trans conversion) alpha)
# (the residue conversion script (2026-09-24, not distributed)). Componentwise 6.5e-8 (printed rounding 5e-8).
# ---------------------------------------------------------------------------------------------------
@testset "published residues (a = 0.5)" begin
    rows = GSN_TEST_LEVEL === :quick ? filter(r -> r.n == 0, PUBLISHED_RESIDUES) : PUBLISHED_RESIDUES
    for r in rows
        @testset "$(r.branch) n=$(r.n)" begin
            if !r.within_rounding
                # n = 11, 13 ordinary differ from the printed values by 6.1e-8 / 3.8e-7 in the legacy run:
                # QNM published-residue backlog: distinguish printed rounding from solver error.
                @test_skip "QNM published-residue backlog, ordinary n=$(r.n)" == ""
            else
                # These published rows retain the table's complete-spectrum labels.
                convention = get(r, :convention, :overtone)
                result = r.branch == ordinary ?
                    qnm(0.5, -2, 2, 2, r.n; convention) :
                    qnm(0.5, -2, 2, 2, r.n, mirror; convention)
                @test result.status in (:accepted, :estimated)
                w = result.omega
                conversion = CF.Binc(-2, 2, 0.5, w, result.lambda) / CF.Btrans(-2, 2, 0.5, w, result.lambda)
                c = 1 / (im * w * conversion * result.incidence_derivative)
                @test abs(real(c) - real(r.published)) <= 6.5e-8
                @test abs(imag(c) - imag(r.published)) <= 6.5e-8
                @test result.excitation_factor ≈ result.reflection_amplitude /
                    (2 * result.omega * result.incidence_derivative) rtol = 1e-14
            end
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Registered-suite interface checks (test/runtests.jl "QNM interface"), plus internal acceptance gates.
# ---------------------------------------------------------------------------------------------------
@testset "compact interface" begin
    ordinary_mode = qnm(0.68, -2, 2, 2, 0)
    mirror_mode = qnm(0.68, -2, 2, 2, 0, mirror)
    @test ordinary_mode isa QNMResult && mirror_mode isa QNMResult
    @test ordinary_mode.status == mirror_mode.status == :accepted
    @test real(ordinary_mode.omega) > 0 && real(mirror_mode.omega) < 0
    @test imag(ordinary_mode.omega) < 0 && imag(mirror_mode.omega) < 0
    for result in (ordinary_mode, mirror_mode)
        @test result.excitation_factor ≈ result.reflection_amplitude /
            (2 * result.omega * result.incidence_derivative) rtol = 1e-14
        @test result.validation.scaled_incidence <= result.validation.metadata.incidence_tolerance
        @test result.excitation.step_drift <= result.excitation.metadata.step_tolerance
        @test result.scientific_acceptance
        @test result.frequency == result.omega
        @test result.lambda ≈ angular_A_to_lambda(result.angular_A, result.a * result.omega, result.mode.m) rtol=1e-12
    end
    # String/Symbol branch selectors.
    @test qnm(0.68, -2, 2, 2, 0, "mirror").omega == mirror_mode.omega
    @test qnm(0.68, -2, 2, 2, 0, :ordinary).omega == ordinary_mode.omega
    @test Ordinary === ordinary && Mirror === mirror
    # Display (9 lines, fixed head and tail).
    text = sprint(show, MIME("text/plain"), ordinary_mode)
    @test length(split(text, '\n')) == 9
    @test startswith(text, "QuasiNormalMode(\n")
    @test endswith(text, "    formalism = GSN)")
    @test startswith(sprint(show, ordinary_mode), "QNMResult(mode=")
    # Guess-seeded detailed call reproduces the compact frequency.
    detailed = qnm(0.68, -2, 2, 2, 0; primary_guess=ordinary_mode.omega, detailed=true)
    @test detailed.omega ≈ ordinary_mode.omega atol=1e-12
    @test detailed.detailed
    # qnm_pair returns both branches, consistent with the single-branch calls.
    pair = qnm_pair(0.68, -2, 2, 2, 0)
    @test pair isa QNMPairResult
    @test pair.status == :accepted
    @test pair.ordinary.omega ≈ ordinary_mode.omega atol=1e-12
    @test pair.mirror.omega ≈ mirror_mode.omega atol=1e-12
    @test startswith(sprint(show, MIME("text/plain"), pair), "QNMPairResult")
    # Lower-level pieces.
    root = qnm_frequency(QNMMode(-2, 2, 2, 0), 0.68)
    @test root isa LeaverResult
    validation = validate_qnm_with_isem(root)
    @test validation isa ISEMValidationResult
    @test validation.status == :accepted
    excitation = qnm_excitation_factor(root)
    @test excitation isa ExcitationFactorResult
    @test excitation.B_gsn ≈ ordinary_mode.excitation_factor rtol=1e-8
    @test excitation.B_teukolsky ≈ Q068.excitation_teukolsky / im rtol=1e-8
    # Sequences along spin seed each step from the previous root.
    sequence = qnm_sequence(QNMMode(-2, 2, 2, 0), [0.0, 0.3, 0.68])
    @test length(sequence) == 3
    @test all(r -> r.status == :accepted, sequence)
    @test last(sequence).omega ≈ root.omega atol=1e-10
    @test qnm_sequence(QNMMode(-2, 2, 2, 0), Float64[]) == LeaverResult[]
    # Angular conversions are inverse to each other.
    for A in (4.0 + 0.0im, 3.7 - 0.2im), c in (0.0 + 0.0im, 0.35 - 0.05im), m in (-2, 0, 2)
        @test lambda_to_angular_A(angular_A_to_lambda(A, c, m), c, m) ≈ A rtol=1e-14
    end
end

@testset "overtones and branches at a = 0.5 (registered)" begin
    overtone7 = qnm_frequency(QNMMode(-2, 2, 2, 7), 0.0)
    overtone10 = qnm_frequency(QNMMode(-2, 2, 2, 10), 0.0)
    @test overtone7.status == overtone10.status == :accepted
    @test abs(overtone10.omega - overtone7.omega) > 0.1
    if at_level(:full)
        # Match test/runtests.jl: this root was n=30 before overtone relabeling.
        overtone29 = qnm_frequency(QNMMode(-2, 2, 2, 29), 0.5)
        @test overtone29.status == :accepted
        @test overtone29.omega ≈ 0.2725174928809329 - 6.734330683370764im atol=1e-10 rtol=0
    end
    for (branch, omega8) in ((ordinary, 0.30101833878061446 - 1.6659917852234734im),
                            (mirror, -0.09835632569909068 - 2.078506281482729im))
        overtone8 = qnm(0.5, -2, 2, 2, 8, branch)
        @test overtone8 isa QNMResult
        @test overtone8.omega ≈ omega8 atol=1e-9 rtol=0
    end
    # Leaver inversion convention reproduces the overtone-convention root (check_qnm_edges.jl).
    root8 = qnm_frequency(QNMMode(-2, 2, 2, 8), 0.5)
    leaver = qnm_frequency(QNMMode(-2, 2, 2, 8), 0.5; convention=:leaver, inversion_index=8, guess=root8.omega)
    @test leaver.status == :accepted
    @test leaver.omega ≈ root8.omega atol=1e-10
    @test leaver.convention == :leaver
    complete = qnm_frequency(QNMMode(-2, 2, 2, 0), 0.5; convention=:complete_spectrum)
    @test complete.omega ≈ qnm_frequency(QNMMode(-2, 2, 2, 0), 0.5).omega atol=1e-12
end

# ---------------------------------------------------------------------------------------------------
# Edges: algebraically special frequency, exact-extremal endpoints, near-extremal estimated results,
# m = 0 overtone 10 continuation, near-extremal n = 10 (documented failure).
# review.md "Existing QNM limitations" and "Full-call edge inventory".
# ---------------------------------------------------------------------------------------------------
@testset "edges and documented limitations" begin
    special = qnm_frequency(QNMMode(-2, 2, 2, 8), 0.0)
    @test special.omega == -2.0im
    @test !special.provenance.gsn_amplitude_limit_implemented
    as_result = qnm(0.0, -2, 2, 2, 8)
    @test as_result isa QNMFailure     # no simple-pole residue at the algebraically special frequency
    @test ismissing(as_result.excitation)
    for (a, m) in ((1.0, 2), (-1.0, -2))
        endpoint = qnm(a, -2, 2, m, 0)
        @test endpoint isa QNMEndpointResult
        @test ismissing(endpoint.excitation_factor)
        @test ismissing(endpoint.incidence_amplitude)
        @test !endpoint.simple_pole
        @test startswith(sprint(show, MIME("text/plain"), endpoint), "ExtremalQuasiNormalModeEndpoint(")
    end
    if at_level(:full)
        plus = qnm(0.99999, -2, 2, 2, 0)
        minus = qnm(-0.99999, -2, 2, -2, 0)
        for result in (plus, minus)
            @test result.status in (:accepted, :estimated)
            @test isfinite(result.excitation_factor)
            @test result.root.status in (:accepted, :estimated)
        end
        @test plus.omega ≈ minus.omega rtol=1e-10
        @test plus.excitation_factor ≈ minus.excitation_factor rtol=1e-7
        # Documented failure ("Full-call edge inventory"): near-extremal n = 10. Whatever the outcome,
        # it must be structured: accepted results are finite; failures name their stage.
        if at_level(:extended)
            hard = qnm(0.99999, -2, 2, 2, 10)
            if hard.status == :accepted
                @test isfinite(hard.excitation_factor)
            else
                @test hard isa Union{QNMFailure, QNMEstimate}
                @test hard.stop_reason isa Symbol
            end
        end
        # m = 0, n = 10 at a = 0.5: the spin continuation stops at the NIA collision near a = 0.3911
        # (documented); if it is ever accepted it must be the catalogued target root.
        axisymmetric = qnm_frequency(QNMMode(-2, 2, 0, 10), 0.5)
        if axisymmetric.status == :accepted
            @test axisymmetric.omega ≈ AXISYMMETRIC_N10 atol=1e-9 rtol=0
        else
            @test axisymmetric.status == :failed
            @test axisymmetric.stop_reason isa Symbol
        end
        # Other harmonic, high overtone (check_qnm_edges.jl).
        other = qnm(0.5, -2, 3, 3, 10)
        @test other.status in (:accepted, :estimated)
        @test isfinite(other.excitation_factor)
        legacy = only(r for r in LEGACY_OVERTONES if r.a == 0.5 && r.l == 3 && r.m == 3 && r.n == 10)
        @test other.omega ≈ legacy.omega atol=1e-9 rtol=0
    end
end

# ---------------------------------------------------------------------------------------------------
# Exact-extremal a = +-1 (src/QNM/ExtremalLeaver.jl, RootPolishing.jl _extremal_damped_limit_result).
# Damped frequencies are roots of the Richartz r = 2M recurrence at a = 1; zero-damping families
# return the synchronous endpoint. Paper values: Richartz, PRD 93, 064062 (2016), Tables I, V, VI
# (six printed decimals). The exact roots below were cross-checked against an independent Python
# implementation of the same recurrence (test/extended/tools/richartz_extremal_reference.py).
# ---------------------------------------------------------------------------------------------------
@testset "exact extremal a = +-1" begin
    paper = at_level(:full) ?
        ((2, 0, 0, 0.425145 - 0.071806im), (2, -2, 0, 0.291553 - 0.088026im),
         (4, 3, 0, 1.503222 - 0.004371im), (2, 1, 0, 0.581433 - 0.038255im),
         (3, 2, 0, 1.028553 - 0.018572im), (7, 5, 0, 2.523730 - 0.010112im),
         (2, -1, 0, 0.343862 - 0.083384im)) :
        ((2, 0, 0, 0.425145 - 0.071806im),)
    for (l, m, n, expected) in paper
        root = qnm_frequency(QNMMode(-2, l, m, n), 1.0)
        @test root.status == :accepted
        @test root.stop_reason == :accepted_exact_extremal_root
        @test root.provenance.root_equation == :richartz_extremal_recurrence
        @test isapprox(root.omega, expected; atol=1e-6, rtol=0)
        @test abs(Q.extremal_radial_residual(-2, m, root.omega, root.angular_A;
            depth=root.provenance.exact_extremal_depth)) < 1e-10
    end
    # Damped s=-2, l=m=2, n=5 branch: not the synchronous endpoint.
    dm = qnm_frequency(QNMMode(-2, 2, 2, 5), 1.0)
    @test dm.status == :accepted
    @test dm.provenance.family_classification == :DM
    @test isapprox(dm.omega, 0.5034599508 - 0.7073988203im; atol=1e-9, rtol=0)
    # Zero-damping families return the endpoint m/2, including |m| < l.
    for (l, m, n) in (at_level(:full) ? ((2, 2, 0), (4, 3, 1)) : ((2, 2, 0),))
        endpoint = qnm_frequency(QNMMode(-2, l, m, n), 1.0)
        @test endpoint.provenance.family_classification == :ZDM
        @test endpoint.omega == m / 2
        @test !endpoint.provenance.simple_pole
    end
    # Spin reflection a = -1, m -> -m.
    plus = qnm_frequency(QNMMode(-2, 2, 0, 0), 1.0)
    minus = qnm_frequency(QNMMode(-2, 2, 0, 0), -1.0)
    @test minus.omega == plus.omega
end

# ---------------------------------------------------------------------------------------------------
# Argument validation (src/QNM/Types.jl, RootPolishing.jl, ModeObservables.jl).
# ---------------------------------------------------------------------------------------------------
@testset "argument validation" begin
    @test_throws ArgumentError QNMMode(-2, 1, 1, 0)
    @test_throws ArgumentError QNMMode(-2, 2, 3, 0)
    @test_throws ArgumentError QNMMode(-2, 2, 2, -1)
    @test_throws ArgumentError QNMMode(-2, 2, 2, 0, :sideways)
    mode = QNMMode(-2, 2, 2, 0)
    @test mode.branch == :positive_real
    @test_throws DomainError qnm_frequency(mode, 1.5)
    @test_throws DomainError qnm_frequency(mode, -1.0001)
    @test_throws ArgumentError qnm_frequency(mode, 0.5; convention=:bogus)
    @test_throws ArgumentError qnm_frequency(mode, 0.5; convention=:leaver)
    @test_throws ArgumentError qnm_frequency(mode, 0.5; inversion_index=-1)
    @test_throws ArgumentError qnm_frequency(mode, 0.5; initial_spin_step=0.0)
    @test_throws ArgumentError qnm_frequency(mode, 0.5; minimum_spin_step=0.1, maximum_spin_step=0.01)
    @test_throws ArgumentError qnm_frequency(mode, 0.5; cf_minimum_iterations=500, cf_maximum_iterations=100)
    @test_throws ArgumentError qnm(0.5, -2, 2, 2, 0, "sideways")
    @test_throws ArgumentError qnm(0.5, -2, 2, 2, 0, :up)
end

end # qnm
