using Test
using MMIEF

@testset "Versioned output and round-trip" begin
    p=Config(rmax=80.0,dr=0.5,tmax=1.0,pulse_amplitude=0.0,
             output_dt=0.5,snapshot_dt=0.5,profile_stride=2)
    data=run_evolution(p;verbose=false)
    mktempdir() do dir
        output=save_run(data,joinpath(dir,"case"))
        @test isfile(joinpath(output,"metadata.toml"))
        @test isfile(joinpath(output,"config.toml"))
        @test isfile(joinpath(output,"diagnostics.csv"))
        @test isfile(joinpath(output,"observers.csv"))
        @test isfile(joinpath(output,"snapshots.csv"))
        recovered=load_run(output)
        @test recovered.times==data.times
        @test recovered.shifts==data.shifts
        @test recovered.masses_h==data.masses_h
        @test recovered.signals==data.signals
        @test recovered.frame_times==data.frame_times
        @test recovered.frame_profiles==data.frame_profiles
        @test recovered.final_state==data.final_state
        @test recovered.settings.observers==data.settings.observers
        @test isfile(save_csv(data,joinpath(dir,"legacy.csv")))
    end
end
