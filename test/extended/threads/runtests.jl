using LinearAlgebra

if Threads.nthreads() == 1
    println("Thread checks skipped: run with --threads=4.")
else
    BLAS.set_num_threads(1)
    radial_values(sol) = (sol.transmission_amplitude,sol.incidence_amplitude,
        sol.reflection_amplitude,sol(5.),sol(10.))
    mode_values(sol) = (sol.amplitude,sol.energy_flux,
        sol.angular_momentum_flux,sol.Carter_const_flux)
    jobs = Function[]
    for s in (-2,0,2), w in (.3,.5-.1im), bc in (IN,UP)
        push!(jobs,()->radial_values(Teukolsky_radial(s,2,2,.7,w,bc)))
    end
    for (a,w) in ((0.,.01+.3im),(.9,.1-2im),(1.,.7),(-1.,.7)), bc in (IN,UP)
        push!(jobs,()->radial_values(GSN_radial(-2,2,2,a,w,bc)))
    end
    for branch in (ordinary,mirror)
        push!(jobs,()->begin
            sol=qnm(.68,-2,2,2,0,branch;detailed=true)
            @test sol.status == :accepted
            (sol.omega,sol.lambda,sol.reflection_amplitude,
                sol.incidence_derivative,sol.excitation_factor,
                sol.X(5.),sol.R(5.),sol.Y(5.))
        end)
    end
    for n in (12,20)
        push!(jobs,()->begin
            sol=qnm_frequency(QNMMode(-2,2,2,n),.5)
            @test sol.status == :accepted
            (sol.omega,sol.lambda,sol.root_residual,sol.cf_error)
        end)
    end
    for args in ((-2,2,2,0,0,.7,8.,0.,1.),
                 (-2,2,2,1,0,.9,10.,.3,1.),
                 (-2,2,1,0,1,.9,10.,0.,.5),
                 (-2,2,2,1,1,.9,10.,.2,.3))
        push!(jobs,()->mode_values(Teukolsky_pointparticle_mode(args...)))
    end
    @testset "serial and threaded public calls" begin
        expected=[job() for job in jobs]
        for order in (eachindex(jobs),reverse(eachindex(jobs)),
                circshift(collect(eachindex(jobs)),7))
            actual=Vector{Any}(undef,length(jobs))
            @sync for i in order
                Threads.@spawn actual[i]=jobs[i]()
            end
            @test isequal(actual,expected)
        end
    end
    @testset "concurrent evaluation of one result" begin
        radii=[3.,5.,10.,50.]
        for w in (.05,.5-.1im), bc in (IN,UP)
            reference=Teukolsky_radial(-2,2,2,.5,w,bc)
            expected=[reference.Teukolsky_solution(r) for r in radii]
            shared=Teukolsky_radial(-2,2,2,.5,w,bc)
            actual=Vector{Any}(undef,length(radii))
            @sync for i in reverse(eachindex(radii))
                Threads.@spawn actual[i]=shared.Teukolsky_solution(radii[i])
            end
            @test isequal(actual,expected)
        end
    end
    @testset "nested threaded sampling" begin
        args = ((-2,2,2,1,0,.9,10.,.3,1.),
                (-2,3,1,2,1,.9,10.,.2,.3))
        expected = [mode_values(Teukolsky_pointparticle_mode(job...;
            threaded_sampling=false)) for job in args]
        for _ in 1:3
            actual = Vector{Any}(undef,length(args))
            @sync for i in eachindex(args)
                Threads.@spawn actual[i] = mode_values(
                    Teukolsky_pointparticle_mode(args[i]...;threaded_sampling=true))
            end
            @test isequal(actual,expected)
        end
    end
end
