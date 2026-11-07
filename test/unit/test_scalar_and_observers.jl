using Test
using MMIEF

@testset "Scalar pulse initial data" begin
    p=Config(rmax=80.0,dr=0.25,tmax=1.0)
    @test MMIEF.ScalarField.gaussian(2.0,p)==(0.0,0.0)
    g,gp=MMIEF.ScalarField.gaussian(20.0,p)
    @test isapprox(g,p.pulse_amplitude;atol=1e-16)
    @test abs(gp) < 1e-14
    @test MMIEF.ScalarField.gaussian(8.0,p)[1] > 0
end

@testset "Fixed areal radius observers" begin
    p=Config(rmax=80.0,dr=0.25,tmax=1.0)
    r=MMIEF.radial_grid(p)
    u=zeros(length(r),5)
    # Polynomial profile: 4-point Lagrange interpolation is exact.
    for i in eachindex(r)
        u[i,5]=r[i]^3 - 2*r[i]^2 + 5*r[i] + 1
        u[i,1]=1.0
    end
    f=0.17; J_obs=15.4; x=J_obs-f
    expected=x^3-2*x^2+5*x+1
    @test MMIEF.Observers.sample_fixed_J(u,r,f,J_obs,p;time=2.0) ≈ expected atol=1e-9
    e=try
        MMIEF.Observers.sample_fixed_J(u,r,9.0,10.0,p;time=25.0)
        nothing
    catch err
        err
    end
    @test e isa ErrorException
    @test occursin("t/M0=25",sprint(showerror,e))
end
