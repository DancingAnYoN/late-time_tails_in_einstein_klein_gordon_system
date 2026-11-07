# Optional GLMakie smoke test; run from root:
#   julia --project=viz test/visual_smoke.jl
using Test
using MMIEF
include(joinpath(@__DIR__,"..","viz","Visualization.jl"))
using .MMIEFVisualization

@testset "GLMakie snapshot rendering after loading saved data" begin
    p=Config(rmax=70.0,dr=0.5,tmax=2.0,pulse_amplitude=0.0,
             output_dt=0.5,snapshot_dt=0.5,profile_stride=2)
    result=run_evolution(p;verbose=false)
    mktempdir() do dir
        save_run(result,dir)
        loaded=load_run(dir)
        make_dashboard(loaded;outdir=dir,movie=false,show_window=false)
        @test isfile(joinpath(dir,"dashboard.png"))
        make_pulse_movie(loaded;outdir=dir,fps=4,show_window=false)
        @test isfile(joinpath(dir,"scalar_pulse_absorption.mp4"))
    end
end
