@testset "generic deep stop retains branch traversal" begin
    values = Dict(2=>0.5, 3=>0.5, 4=>8.0, 5=>10.0, 6=>10.0,
                  7=>1e-6, 8=>4e-6, 9=>0.5, 10=>0.1)
    function walked(deep)
        seen = Int[]
        MS._branch_l_walk(2, 6, 10, 1.0, 2, true; deep=deep) do l
            push!(seen,l)
            values[l]
        end
        seen
    end
    @test walked(0.0) == [2,3,5,4,6,7,8]
    @test walked(MS.DEFAULT_GENERIC_DEEP) == [2,3,5,4,6,7]
end

@testset "generic stopping keywords reach the public flux route" begin
    orbit = (0.9, 10.0, 0.2, 0.3)
    public_flux = GSN.Teukolsky_pointparticle_flux(orbit...;
        lmax=2, kmax=0, nmax=3, tail_levin=false, generic_beta=0.0, generic_deep=0.0)
    direct_flux = MS.generic_mode_summation(orbit...;
        lmax=2, kmax=0, nmax=3, tail_levin=false, beta=0.0, deep=0.0)
    for field in (:infinity_energy_flux, :infinity_angular_momentum_flux,
                  :infinity_carter_constant_flux, :horizon_energy_flux,
                  :horizon_angular_momentum_flux, :horizon_carter_constant_flux)
        @test isequal(getproperty(public_flux, field), getproperty(direct_flux, field))
    end
end

@testset "generic tail budget preserves horizon peak discovery" begin
    # Horizon branch at (0.9,6,0.7,0.5), n=56, m=2, k=7.
    # Raising both search cutoffs leaves l=4:6 unseen, including the peak at l=5.
    values = Dict(2=>1.8694408464497535e-18, 3=>3.2505517211626002e-16,
        4=>3.483988658771645e-14, 5=>6.442512267728932e-13,
        6=>1.961873128265081e-14, 7=>6.756120858072641e-16,
        8=>2.3245606848495408e-17, 9=>4.888890535213275e-19,
        10=>2.3516093606890773e-21)
    budget_cutoff = MS.DEFAULT_GENERIC_BETA * 1e-6 * 9.2e-6
    function peak_search(search_cutoff)
        seen=Int[]
        MS._branch_l_walk(2,9,10,budget_cutoff,2,true;
            deep=MS.DEFAULT_GENERIC_DEEP, search_threshold=search_cutoff) do l
            push!(seen,l)
            values[l]
        end
        seen
    end
    raised_search=peak_search(budget_cutoff)
    original_search=peak_search(1e-16)
    @test raised_search == [2,3,8,7,9,10]
    @test original_search == collect(2:9)
    @test 5 in original_search && !(5 in raised_search)
    @test length(unique(original_search)) == length(original_search)
end