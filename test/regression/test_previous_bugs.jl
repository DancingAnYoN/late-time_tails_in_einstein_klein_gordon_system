using Test
using MMIEF

@testset "Regression: no forbidden GLMakie core dependency" begin
    # Headless solver package does not import visualization.
    @test !isdefined(MMIEF,:GLMakie)
    @test isdefined(MMIEF,:Evolution)
    @test isdefined(MMIEF,:Integrators)
end

@testset "Regression: vacuum horizon shift" begin
    p=Config(rmax=80.0,dr=0.5,cfl=0.2,tmax=0.1,
             pulse_amplitude=0.0,output_dt=0.1,snapshot_dt=0.1)
    r,u=initial_data(p)
    wk=MMIEF.Evolution.EvolutionWorkspace(u,r)
    MMIEF.Evolution.setup_vacuum_balance!(p,wk.rhswork)
    MMIEF.Constraints.project_hamiltonian!(u,0.0,p,r)
    f,drift=MMIEF.Evolution.evolve_one_step!(u,0.0,p.dr*p.cfl,p,wk)
    @test f == 0.0
    @test isfinite(drift)
    @test isapprox((2*p.M0+f)/2,p.M0;atol=1e-12)
end
