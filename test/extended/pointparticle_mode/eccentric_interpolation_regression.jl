@testset "eccentric checked radial interpolation" begin
    ci = GSN.ConvolutionIntegrals
    gs = GSN.GridSampling
    # The last three cases pin the three non-default branches of the table build (found by scanning the
    # V3 high-n grid): deferred (no refinement large enough to pay for the table), fallback (no 512-point
    # table passes the self-check; exact samples are used), and a table accepted after one doubling.
    cases = (
        ((0.9, 6.0, 0.7), (-2, 2, 2, 10), :any),
        ((0.9, 6.0, 0.7), (-2, 9, 9, 100), :any),
        ((0.9, 6.0, 0.7), (-2, 12, 11, 150), :any),
        ((0.9, 6.0, 0.7), (2, 4, 3, 60), :any),
        ((0.5, 8.0, 0.9), (-2, 6, 4, 200), :any),
        ((0.5, 8.0, 0.9), (2, 3, 2, 120), :any),
        ((0.95, 7.0, 0.85), (2, 12, 12, 200), :deferred),
        ((0.5, 8.0, 0.9), (2, 2, 2, 200), :fallback),
        ((0.5, 8.0, 0.9), (2, 5, 5, 150), :doubled),
    )
    for (orbit, mode, branch) in cases
        @testset "$orbit $mode" begin
            s, l, m, n = mode
            kernel = ci.convolution_integral_eccentric_trapezoidal_isem(orbit..., s, l, m, n, 64)
            sample = kernel["Trajectory"]
            omega = kernel["omega"]
            harmonic = ci.spin_weighted_spheroidal_harmonic(s, l, m, orbit[1] * omega; method = "jacobi")
            solution = ci._isem_y_solution(s, l, m, orbit[1], omega; lambda = harmonic.lambda)
            exact = s == -2 ? gs.y_sample_m2_isem(solution, sample) : gs.y_sample_p2_isem(solution, sample)
            for component in ("Y", "Yp", "X")
                delta = maximum(abs.(kernel["YSolution"][component] .- exact[component]))
                @test delta <= gs.YINTERP_CHECK_RTOL * maximum(abs.(exact[component]))
            end
            @test kernel["YInterpolationPoints"] <= gs.YINTERP_MAX_POINTS
            branch == :deferred && @test kernel["YInterpolationDeferred"]
            branch == :fallback && @test kernel["YInterpolationFallback"]
            if branch == :doubled
                @test !kernel["YInterpolationDeferred"] && !kernel["YInterpolationFallback"]
                @test kernel["YInterpolationDoublings"] >= 1
            end
            if kernel["YInterpolationDeferred"]
                @test kernel["YInterpolationPoints"] == 0
                @test isnan(kernel["YInterpolationSelfCheckError"])
                @test !kernel["YInterpolationFallback"]
                @test all(isequal(kernel["YSolution"][q], exact[q]) for q in ("Y", "Yp", "X"))
            elseif kernel["YInterpolationFallback"]
                @test kernel["YInterpolationPoints"] == gs.YINTERP_MAX_POINTS
                @test all(isequal(kernel["YSolution"][q], exact[q]) for q in ("Y", "Yp", "X"))
            else
                @test kernel["YInterpolationSelfCheckError"] <= gs.YINTERP_CHECK_RTOL
            end
            reference = ci._eccentric_flux_from_sample(sample, exact, harmonic, s, orbit[1], omega, m, n)
            # Interpolated radii agree with exact samples to the self-check tolerance, but cancellation in
            # small high-n modes amplifies that difference in the mode flux; e.g. (0.5,8,0.9) n=200 has
            # E≈6e-15 and a relative change of 2e-7. Modes below 1e-12 are therefore checked absolutely.
            @test abs(kernel["EnergyFlux"] - reference["EnergyFlux"]) <=
                  1e-8 * max(abs(reference["EnergyFlux"]), 1e-12)
        end
    end
end
