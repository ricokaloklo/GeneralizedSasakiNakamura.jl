@testset "generic optional Levin correctness" begin
    ci = GSN.ConvolutionIntegrals
    @test ci._generic_adaptive_effective_theta_intervals(16, 2, 8) == 32
    @test ci._generic_adaptive_effective_theta_intervals(16, 14, 20) == 128
    @test ci._generic_adaptive_effective_theta_intervals(256, 14, 20) == 256
    nr, nc = 33, 5
    x = pi .* ([cos(pi*j/(nr-1)) for j in 0:nr-1] .+ 1.0) ./ 2.0
    phase = 17.0 .* x .+ 0.2 .* x.^2
    derivative = 17.0 .+ 0.4 .* x
    factor = ci._levin_1d_local_phase_factor(derivative, phase; local_order=nr)
    values = ComplexF64[exp(0.1im*j*x[i])*(1+x[i]^2/j) for i in 1:nr, j in 1:nc]
    batched = ci._levin_1d_columns_factored(values, factor)
    scalar = [ci._levin_1d_integral_local_phase_factored(collect(values[:,j]), factor) for j in 1:nc]
    @test batched ≈ scalar rtol=1e-12
    weights = [1.0/j for j in 1:nc]
    theta_phase = [0.3*j for j in 1:nc]
    first, second = ci._generic_adaptive_radial_levin_theta_cc(values, 2.0 .* values,
        weights, theta_phase, factor, factor)
    reference = sum(weights .* exp.(1im .* theta_phase) .* scalar)
    @test first ≈ reference rtol=1e-12
    @test second ≈ 2.0*reference rtol=1e-12
    orbit = (0.9, 6.0, 0.7, 0.5)
    mode = (-2, 8, 6, 50, 2)
    MS.with_public_submission(orbit...) do scope
        MS.with_submission_operation(scope) do
            master = MS.submission_master!(scope)
            frequencies = master["Frequencies"]
            s,l,m,n,k = mode
            omega = (m*frequencies["ϒϕ"]+n*frequencies["ϒr"]+k*frequencies["ϒθ"])/frequencies["ϒt"]
            MS.with_submission_mode(scope; mode=mode) do cache
                key = ci.GenericAdaptiveLevin2DKey(3, 0, 6, 0)
                phase_nr = ci.DEFAULT_ADAPTIVE_LEVIN_LOCAL_N + 1
                sample = ci._generic_adaptive_segment_sample!(cache, master, key, phase_nr, 33)
                plan = ci._levin_1d_local_plan(phase_nr, phase_nr)
                for direction in (1,-1)
                    basis = ci._generic_adaptive_phase_basis_build(sample, s, orbit[1], direction)
                    phase_values = real.(ci._glevin_radial(basis, m, n, k))
                    numerical = plan.segments[1].D * phase_values[plan.perm]
                    analytic = ci._generic_adaptive_radial_gprime(sample, key, s, orbit[1], omega,
                        m,n,direction)[plan.perm]
                    @test maximum(abs.(numerical.-analytic)) <= 1e-12 * maximum(abs.(analytic))
                end
                trapezoid = ci.generic_mode_flux_from_master_cached!(cache, master, mode...)
                levin = ci.generic_mode_flux_from_master_cached_adaptive_levin!(cache, master, mode...;
                    sample_tol=1e-3, tol=1e-8, mode_abs_floor=0.0)
                @test levin["AdaptiveLevinTol0"] == min(1e-3,3sqrt(1e-8))
                @test levin["AdaptiveLevinEffectiveThetaIntervals"] >= 32
                @test levin["EnergyFlux"] ≈ trapezoid["EnergyFlux"] rtol=3sqrt(1e-8)
                @test levin["CarterConstantFlux"] ≈ trapezoid["CarterConstantFlux"] rtol=3sqrt(1e-8)
            end
        end
    end
end