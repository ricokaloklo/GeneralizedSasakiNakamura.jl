@testset "basis coverage failure is retried, then recorded" begin
    CI = GSN.ConvolutionIntegrals
    calls = Ref(0)
    covered_when_ordinary = () -> begin
        calls[] += 1
        get(task_local_storage(), CI._FORCE_ORDINARY_KEY, false) ? Dict("EnergyFlux" => 1.0) :
            error("ISEM GSN basis does not cover x=0.93")
    end
    @test CI._with_coverage_retry(covered_when_ordinary, 2, 2, 2, 5, 0)["EnergyFlux"] == 1.0
    @test calls[] == 2
    failed = @test_logs (:warn,) CI._with_coverage_retry(() -> error("ISEM GSN basis does not cover x=0.93"), 2, 2, 2, 5, 0)
    @test failed["FailedMode"] && failed["EnergyFlux"] == 0.0
    @test_throws ErrorException CI._with_coverage_retry(() -> error("unrelated"), 2, 2, 2, 5, 0)
    @test_throws InterruptException CI._with_coverage_retry(() -> throw(InterruptException()), 2, 2, 2, 5, 0)
end

@testset "polar orbits with x = 0 are finite" begin
    f0 = GSN.Teukolsky_pointparticle_flux(0.9, 8.0, 0.0, 0.0; tol = 1e-6)
    fm = GSN.Teukolsky_pointparticle_flux(0.9, 8.0, 0.0, -0.0; tol = 1e-6)
    @test isfinite(f0.infinity_energy_flux) && isfinite(f0.horizon_energy_flux)
    @test f0.infinity_energy_flux == fm.infinity_energy_flux
    # the x -> 0 limit: the mean of x = +-1e-3 is 1.8694714e-4 (an O(x^2) = 1.5e-7 shift)
    @test isapprox(f0.infinity_energy_flux, 1.869471077e-4; rtol = 1e-7)
    @test haskey(f0.reached, :small_x_interpolation)
end

@testset "nearly polar orbits are interpolated through x = 0" begin
    # at a = 0 the flux does not depend on x; |x| = 1e-6 used to be off by 1e-5 and 1e-8 by a factor 6
    ref = GSN.Teukolsky_pointparticle_flux(0.0, 25.0, 0.0, 0.3; tol = 1e-8).infinity_energy_flux
    for x in (1e-8, -1e-6, 0.0)
        @test isapprox(GSN.Teukolsky_pointparticle_flux(0.0, 25.0, 0.0, x; tol = 1e-8).infinity_energy_flux, ref; rtol = 1e-8)
    end
end

@testset "retrograde equatorial fallbacks keep the sign of L" begin
    f = GSN.Teukolsky_pointparticle_flux(0.9, 10.0, 0.0, -1.0; tol = 1e-6)
    g = MS.inclined_mode_summation(0.9, 10.0, -1.0; tol = 1e-6)
    @test f.infinity_angular_momentum_flux < 0
    @test g.infinity_angular_momentum_flux ≈ f.infinity_angular_momentum_flux rtol = 1e-6  # minimum_consecutive 2 vs 5
end

@testset "near-static modes that fail in the radial solver are recorded as zero" begin
    CI = GSN.ConvolutionIntegrals
    MSTE = GSN.ISEM.ManoSuzukiTakasugiInfinity.ManoSuzukiTakasugiCertificateError
    Ωφ = 1e-10
    for err in (ErrorException("nonfinite ISEM GSN matching scale"), ErrorException("ISEM GSN log-scaled route lacks a physical MST handoff state"),
                ErrorException("Non-finite flux for mode (s, l, m, n, k) = (-2, 7, -6, 0, 6)"), MSTE(:physical_in_residual, "MST physical-in residual rejected"))
        @test CI._near_static_overflow(err, 1e-20, Ωφ)
        @test !CI._near_static_overflow(err, 1e-15, Ωφ)          # not near-static
    end
    @test !CI._near_static_overflow(ErrorException("unrelated"), 1e-20, Ωφ)
    @test !CI._near_static_overflow(ArgumentError("ISEM"), 1e-20, Ωφ)
    zeroed = CI._with_coverage_retry(() -> error("ISEM GSN infinity recurrence lost Float64 credibility"), -2, 3, -3, 0, 3; omega = 1e-20, Ωφ = Ωφ)
    @test zeroed["NearStaticOverflow"] && zeroed["EnergyFlux"] == 0.0    # a high-l, low-frequency mode whose matching overflows is negligible: (8.1e-5 * 167)^(2 * 67) ~ 1e-250
    overflow = ErrorException("nonfinite ISEM GSN matching numerator scale")
    @test CI._negligible_overflow(overflow, 8.1e-5, 68, 167.0)
    @test !CI._negligible_overflow(overflow, 8.1e-5, 4, 167.0)
    @test !CI._negligible_overflow(ErrorException("unrelated"), 8.1e-5, 68, 167.0)
    low = CI._with_coverage_retry(() -> throw(overflow), -2, 68, -67, 84, 0; omega = 8.1e-5, Ωφ = 1.5e-2, rmax = 167.0)
    @test low["EnergyFlux"] == 0.0
end

@testset "eccentric m walks reach the peak of each n shell" begin
    # References: Nasipak's FEWv2 Kerr eccentric equatorial fluxes (pybhpt), rows 167369 and 83 of regions A + B.
    # Retrograde (a < 0 there, x = -1 here): the negative-m branch used to stop at m = -2 for n = 22, 25, 27.. (-7.6e-4).
    f = GSN.Teukolsky_pointparticle_flux(0.5247987322667358, 9.517944054644722, 0.5312200481750099, -1.0; tol = 1e-8)
    @test isapprox(f.infinity_energy_flux, 2.0813683218179032e-4; rtol = 1e-7)
    @test isapprox(f.infinity_angular_momentum_flux, -3.7099764781389254e-3; rtol = 1e-7)
    # Near the ISSO with e = 0.009 the n = 6 shell peaks at m = l ~ 20; the walk used to stop at m = 2 (-4e-6).
    g = GSN.Teukolsky_pointparticle_flux(0.9320585935700156, 2.09532934584992, 0.00914066143877242, 1.0; tol = 1e-8)
    @test isapprox(g.infinity_energy_flux, 0.049183278297751945; rtol = 2e-7)
end

@testset "geometric tails of the l and n sums" begin
    # a geometric tail r = 0.9 after a steady fall: factor r / (1 - r)
    shells = [0.9^k for k in 0:20]
    factor, r = MS._geometric_shell_tail(shells, 3)
    @test r ≈ 0.9 rtol = 1e-12
    @test factor ≈ 9.0 rtol = 1e-10
    @test MS._geometric_shell_tail([1.0, 0.5, 0.8, 0.4, 0.2], 3) == (0.0, 0.0)   # not a steady fall
    @test MS._geometric_shell_tail([1.0, -0.5, 0.25, -0.12, 0.06], 3) == (0.0, 0.0)  # sign changes
    @test MS._geometric_shell_tail([1.0, 0.99, 0.98, 0.97, 0.96], 3) == (0.0, 0.0)   # r above N_TAIL_MAX_RATIO
    # the n tail only at high eccentricity and long n walks
    @test MS._geometric_n_tail_factor(shells, 0.8, 150, 3)[1] ≈ 9.0 rtol = 1e-10
    @test MS._geometric_n_tail_factor(shells, 0.5, 150, 3) == (0.0, 0.0)
    @test MS._geometric_n_tail_factor(shells, 0.8, 50, 3) == (0.0, 0.0)
    @test MS._geometric_n_tail_factor(shells, 0.8, nothing, 3) == (0.0, 0.0)
    # the zoom-whirl n stop spans one period Ups_phi / Ups_r of the n shells
    @test MS._n_stop_run(2, 0.5) == 3
    @test MS._n_stop_run(2, 1 / 29.4) == 30
end

@testset "circular l tail and near-ISSO m walks against FEWv2" begin
    # References: Nasipak's FEWv2 Kerr equatorial fluxes, rows 50712 (circular) and 80542 of regions A + B.
    # Circular, a = 0.973, p = 1.757: the l sum stopped 1.4e-9 low; with the l tail +2.6e-10.
    c = GSN.Teukolsky_pointparticle_flux(0.9731817507537753, 1.7566650042516612, 0.0, 1.0; tol = 1e-8)
    @test isapprox(c.infinity_energy_flux, 0.07253744207060168; rtol = 6e-10)
    @test c.result.l_tail_infinity.energy > 0
    # Nearly circular near the ISSO: the n = 1 shell dips at |m| = 25..29 and carries 1.5e-8 at |m| = 31..38, which
    # the m walk dropped (-2.0e-7) before it had to pass the |m| reached by the previous shell.
    f = GSN.Teukolsky_pointparticle_flux(0.9909069337322672, 1.5914023608769856, 0.029324604931803146, 1.0; tol = 1e-8)
    @test isapprox(f.infinity_energy_flux, 0.08365952830274999; rtol = 5e-8)
end

@testset "eccentric refinement floor (horizon)" begin
    # FEWv2 row 569994 (a = 0.439, p = 22.7, e = 0.844): thousands of horizon modes near the branch floor kept
    # coarse values, and E_hor at tol 1e-8 was 4.5e-7 off the tol 1e-10 value -1.0555654212e-10 (now 6e-8).
    f = GSN.Teukolsky_pointparticle_flux(0.4388530721544722, 22.723645351550264, 0.84375, 1.0; tol = 1e-8)
    @test isapprox(f.horizon_energy_flux, -1.0555654212e-10; rtol = 2e-7)
    @test f.result.n_tail_infinity.energy > 0
end

@testset "ISEM small-frequency route at high angular compression" begin
    CI = GSN.ConvolutionIntegrals
    a = 0.9974711493723011
    # (l + 1/2)/kappa = 400: the small-frequency route returned UP solutions off by up to 1e96; auto now takes the
    # ordinary route, which agrees with the linear ODE started at r_out = 2e5 to 3e-10.
    for s in (-2, 2)
        auto = GSN.Teukolsky_radial(s, 28, 28, a, 0.072, GSN.UP; method = "ISEM")
        ordinary = GSN.Teukolsky_radial(s, 28, 28, a, 0.072, GSN.UP; method = "ISEM", sfe = false)
        @test all(r -> isapprox(auto(r), ordinary(r); rtol = 1e-12), (5.0, 50.0, 198.0))
    end
    # the Wronskian self-check: constant for the solutions used, a large spread for the small-frequency route there
    @test CI._isem_wronskian_check(2, 28, 28, a, 0.072, (16.8, 57.7, 198.0)) <= CI.ISEM_WRONSKIAN_CHECK_TOL
    sfe_in = GSN.Teukolsky_radial(2, 28, 28, a, 0.072, GSN.IN; method = "ISEM", sfe = true)
    sfe_up = GSN.Teukolsky_radial(2, 28, 28, a, 0.072, GSN.UP; method = "ISEM", sfe = true)
    W(r) = (pin = sfe_in.Teukolsky_solution.pair_value(r); pup = sfe_up.Teukolsky_solution.pair_value(r);
        (r^2 - 2r + a^2)^3 * (pin[1] * pup[2] - pin[2] * pup[1]))
    @test abs(W(198.0) / W(16.8) - 1) > 1e3
    @test CI._is_basis_coverage_error(ErrorException("ISEM radial self-check failed: the Wronskian varies"))
    @test CI._orbit_check_radii(Dict("r" => [4.0, 9.0, 16.0])) == (4.0, 8.0, 16.0)
    @test CI._orbit_check_radii(Dict("p" => 10.0)) == (8.0, 10.0, 12.5)
end

@testset "failed negligible modes are zeroed quietly" begin
    CI = GSN.ConvolutionIntegrals
    fail = () -> error("ISEM GSN basis does not cover x=0.93")
    quiet = @test_logs CI._with_coverage_retry(fail, -2, 9, -9, 0, 4; omega = 5e-7, Ωφ = 1.4e-7, negligible_if_failed = true)
    @test quiet["NearStaticOverflow"] && quiet["EnergyFlux"] == 0.0 && !get(quiet, "FailedMode", false)
    # (l, m, k) = (9, -9, 4) at (a, p, x) = (0.998, 3.7e4, 0.004): |omega| r ~ 0.03, near-zone bound ~ 1e-35 of l = 2
    @test CI._negligible_low_frequency_mode(-2, 9, 5e-7, 1.4e-7, 0.998, 3.7e4, 0.0; bound = CI.NEGLIGIBLE_FAILED_MODE_BOUND)
    @test !CI._negligible_low_frequency_mode(-2, 9, 5e-7, 1.4e-7, 0.998, 3.7e4, 0.0)
    @test !CI._negligible_low_frequency_mode(-2, 2, 5e-7, 1.4e-7, 0.998, 3.7e4, 0.0; bound = CI.NEGLIGIBLE_FAILED_MODE_BOUND)
end

@testset "generic n tail at high eccentricity" begin
    # (a, p, e, x) = (0.5, 20, 0.85, 0.8) with lmax = 3, kmax = 2 (cheap shells): the n walk stops by itself at
    # n = 145 at tol 1e-6; the geometric tail (r = 0.886) brings E_inf from 5.8e-6 to 1.2e-7 of the tol 1e-9 value
    # 1.0528034712962e-6 of the same caps.
    f = GSN.Teukolsky_pointparticle_flux(0.5, 20.0, 0.85, 0.8; tol = 1e-6, lmax = 3, kmax = 2)
    tail = f.result.n_tail_infinity
    @test f.result.n_reached_inf !== nothing && f.result.n_reached_inf >= 100
    @test tail.energy > 0 && 0 < tail.ratio < 1
    @test isapprox(f.infinity_energy_flux, 1.0528034712962e-6; rtol = 5e-7)
    @test !isapprox(f.infinity_energy_flux - tail.energy, 1.0528034712962e-6; rtol = 3e-6)
end

@testset "eccentric adaptive Levin at e = 0.95" begin
    # Deep-refinement references (two-wave and single-phase Levin, 11 levels, agree to 2e-10; trapezoid agrees).
    # (12, -11, 208) once stopped at level 2 where two under-resolved levels agreed on 2.1e-13; (2, 2, 400) once
    # stopped at relative energy change 1e-6 relaxed to 1e-2 with depth, 1.7e-5 off.
    CI = GSN.ConvolutionIntegrals
    scope, owned = MS.ensure_public_submission(nothing, 0.9, 8.0, 0.95, 1.0)
    MS._with_public_scope(scope, owned) do
        MS.with_submission_operation(scope) do
            KG = MS.submission_master!(scope)
            levin(l, m, n) = MS.with_submission_mode(scope; mode = (-2, l, m, n, 0)) do cache
                CI.convolution_integral_eccentric_adaptive_levin_isem(KG, -2, l, m, n; tol = 1e-8,
                    mode_abs_floor = 1e-16, energy_abs_tol = 2.85e-15, cache = cache)["EnergyFlux"]
            end
            @test levin(12, -11, 208) < 1e-20
            ref = levin(2, 2, 400)
            deep = MS.with_submission_mode(scope; mode = (-2, 2, 2, 400, 0)) do cache
                CI.convolution_integral_eccentric_adaptive_levin_isem(KG, -2, 2, 2, 400; tol = 1e-8, tol0 = 1e-12,
                    max_depth = 11, mode_abs_floor = 1e-16, cache = cache)["EnergyFlux"]
            end
            @test ref ≈ deep rtol = 1e-8
        end
    end
end
