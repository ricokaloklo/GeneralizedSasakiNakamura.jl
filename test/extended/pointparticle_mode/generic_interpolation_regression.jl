@testset "generic checked interpolation and nested trapezoid" begin
    ci = GSN.ConvolutionIntegrals
    gs = GSN.GridSampling
    cases = (
        ((0.9, 6.0, 0.7, 0.5), (-2, 8, 6, 50, 2)),
        ((0.9, 6.0, 0.7, 0.5), (2, 20, 18, 84, 3)),
        ((0.9, 10.0, 0.2, 0.3), (-2, 12, 9, 100, 3)),
        ((0.9, 10.0, 0.2, 0.3), (2, 12, 9, 150, 3)),
        ((0.9, 12.0, 0.8, 0.5), (-2, 30, 30, 150, 0)),
        ((0.9, 12.0, 0.8, 0.5), (2, 21, 16, 100, 3)),
        ((0.5, 8.0, 0.9, 0.3), (-2, 25, 3, 150, 15)),
        ((0.5, 8.0, 0.9, 0.3), (2, 12, 9, 50, 3)),
    )
    for (orbit, mode) in cases
        @testset "$orbit $mode" begin
            MS.with_public_submission(orbit...) do scope
                MS.with_submission_operation(scope) do
                    master = MS.submission_master!(scope)
                    MS.with_submission_mode(scope; mode = mode) do cache
                        s, l, m, n, k = mode
                        kernel = ci.generic_mode_flux_from_master_cached!(cache, master, mode...)
                        sample = kernel["Trajectory"]
                        omega = kernel["omega"]
                        harmonic = ci.spin_weighted_spheroidal_harmonic(s, l, m, orbit[1]*omega; method="jacobi")
                        solution = ci._generic_isem_y_solution(s, l, m, orbit[1], omega; lambda=harmonic.lambda)
                        exact = s == -2 ? gs.y_sample_m2_isem(solution, sample) : gs.y_sample_p2_isem(solution, sample)
                        for component in ("Y", "Yp", "X")
                            delta = maximum(abs.(kernel["YSolution"][component] .- exact[component]))
                            scale = maximum(abs.(exact[component]))
                            @test delta <= gs.YINTERP_CHECK_RTOL * scale
                        end
                        @test kernel["YInterpolationPoints"] <= gs.YINTERP_MAX_POINTS
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
                        reference = ci._generic_flux_from_sample(sample, exact, kernel["SWSH"],
                            mode..., orbit[1], omega, master["Frequencies"]["ϒθ"])
                        # Relative mode-flux checks apply to the resolved, non-tiny contribution.
                        if abs(reference["EnergyFlux"]) > 1e-15
                            @test kernel["EnergyFlux"] ≈ reference["EnergyFlux"] rtol=1e-8
                            @test kernel["CarterConstantFlux"] ≈ reference["CarterConstantFlux"] rtol=1e-8
                        end
                    end
                end
            end
        end
    end
end