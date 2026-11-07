using Test
using MMIEF

@testset "Stationary Schwarzschild evolution" begin
    p=Config(rmax=80.0,dr=0.5,cfl=0.2,tmax=2.0,
             pulse_amplitude=0.0,output_dt=0.5,snapshot_dt=0.5)
    data=run_evolution(p;verbose=false)
    @test all(==(0.0),data.shifts)
    @test abs(data.masses_h[end]-p.M0)<1e-12
    @test maximum(abs,data.horizon_error)<1e-9
    @test maximum(abs.(data.final_state[:,1] .-
        [schwarzschild_fields(r,p.M0).a for r in data.r])) < 1e-8
end

@testset "Weak coupled field short evolution" begin
    p=Config(rmax=80.0,dr=0.5,cfl=0.2,tmax=4.0,
             output_dt=0.5,snapshot_dt=1.0)
    data=run_evolution(p;verbose=false)
    @test all(isfinite,data.final_state)
    @test data.masses_h[end] >= p.M0-1e-10
    @test data.masses_h[end] < p.M0+0.01
    @test abs(data.horizon_error[end])<1e-8
    @test all(length(y)==length(data.times) for y in data.signals)
    @test maximum(abs,data.final_state[:,5])>0
end
