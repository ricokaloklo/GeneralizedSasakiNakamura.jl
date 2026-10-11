@testset "generic trapezoid radial aliasing and Carter average" begin
    ci = GSN.ConvolutionIntegrals
    orbit = (0.9, 6.0, 0.7, 0.5)
    mode = (-2, 20, 7, 85, 14)
    public_mode = Teukolsky_pointparticle_mode(mode..., orbit...)
    MS.with_public_submission(orbit...) do scope
        MS.with_submission_operation(scope) do
            master = MS.submission_master!(scope)
            MS.with_submission_mode(scope; mode = mode) do cache
                cached = ci.generic_mode_flux_from_master_cached!(cache, master, mode...)
                uncached = ci.generic_mode_flux_from_master(master, mode...)
                @test cached["EnergyFlux"] ≈ public_mode.energy_flux rtol=1e-8
                @test uncached["EnergyFlux"] ≈ public_mode.energy_flux rtol=1e-8
                @test 1e-18 < cached["EnergyFlux"] < 1e-15
                @test cached["N_sample"] >= nextpow(2, 4 * abs(mode[4]))
                @test uncached["N_sample"] >= nextpow(2, 4 * abs(mode[4]))
                low = ci._cached_generic_sample!(cache, master, 64, 16)
                high = ci._cached_generic_sample!(cache, master, 64, 64)
                for (m, n, k) in ((2, 10, 0), (2, 1, -2), (11, 80, 0))
                    f = master["Frequencies"]
                    omega = (m*f["ϒϕ"] + n*f["ϒr"] + k*f["ϒθ"]) / f["ϒt"]
                    @test ci._generic_carter_factor_min_grid(low, m, omega) ≈
                          ci._generic_carter_factor_min_grid(high, m, omega) rtol=1e-12
                end
            end
        end
    end
    mktempdir() do directory
        record_path = joinpath(directory, "capped_generic.h5")
        MS.generic_mode_summation(orbit...; lmax=2, kmax=0, nmax=3,
            tail_levin=false, record=true, record_path=record_path)
        MS.h5open(record_path, "r") do file
            observed = 0
            for nkey in keys(file)
                startswith(nkey, "n_") || continue
                for kkey in keys(file[nkey]), mkey in keys(file[nkey][kkey]),
                    lkey in keys(file[nkey][kkey][mkey]),
                    boundary in keys(file[nkey][kkey][mkey][lkey])
                    entry = file[nkey][kkey][mkey][lkey][boundary]
                    n = read(entry["n"])
                    @test read(entry["N_sample"]) >= nextpow(2, max(1, 4 * abs(n)))
                    observed += 1
                end
            end
            @test observed > 0
        end
    end
end
