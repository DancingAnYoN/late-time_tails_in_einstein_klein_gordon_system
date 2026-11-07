using Test
using MMIEF

@testset "Legacy implementation equivalence" begin
    # Keep the frozen monolith as a *reference*, not as a runtime dependency.
    if !isdefined(Main,:MMIEFCore)
        include(joinpath(@__DIR__,"..","..","legacy","MMIEFCore.jl"))
    end
    for amp in (0.0,1e-4)
        kwargs=(rmax=80.0,dr=0.5,cfl=0.2,tmax=2.0,
                pulse_amplitude=amp,output_dt=0.5,snapshot_dt=0.5,
                profile_stride=2)
        before=Main.MMIEFCore.run_evolution(Main.MMIEFCore.Config(;kwargs...);
                                            verbose=false)
        after=MMIEF.run_evolution(MMIEF.Config(;kwargs...);verbose=false)
        @test before.times==after.times
        @test before.shifts ≈ after.shifts atol=1e-10
        @test before.masses_h ≈ after.masses_h atol=1e-10
        @test maximum(abs.(before.final_state.-after.final_state))<1e-9
        for j in eachindex(before.signals)
            @test before.signals[j] ≈ after.signals[j] atol=1e-10
        end
    end
end
