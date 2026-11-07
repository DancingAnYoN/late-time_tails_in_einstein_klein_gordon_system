using Test
using MMIEF

@testset "Outer scalar boundary projection" begin
    p=Config(rmax=80.0,dr=0.5,tmax=1.0,pulse_amplitude=0.0)
    r,u=initial_data(p)
    u[end,3]=0.003
    u[end,5]=-0.002
    MMIEF.BoundaryConditions.outer_projection!(u,0.0,p,r)
    fd=MMIEF.Geometry.fdot_horizon(u,0.0,p)
    a=u[end,1]; B=u[end,2]; J=r[end]
    beta=(fd+a*J*B)/(1+a*J*B)
    c=1-2*beta
    residual=(1-beta)*(u[end,4]+u[end,3])+(fd+c)*u[end,5]/J
    @test abs(residual)<1e-14
end

@testset "Synthetic QNM and tail recovery" begin
    t=collect(0.0:0.5:125.0)
    gamma=0.094
    omega=0.111
    y=exp.(-gamma.*t).*cos.(omega.*t.+0.4)
    fit=fit_qnm(t,y;tlo=45.0,thi=100.0,mass=1.0)
    @test fit !== nothing
    if fit !== nothing
        @test abs(fit.omegaM-omega)<0.015
        @test abs(fit.gammaM-gamma)<0.015
    end
    t2=collect(100.0:1.0:500.0)
    y2=t2.^(-3)
    pidx=local_power_index(t2,y2;peak_halfwidth=0.0,
                            fit_halfwidth=30.0,tstart=120.0)
    kept=[x for x in pidx if isfinite(x)]
    @test length(kept)>50
    @test abs(kept[end]-3.0)<0.05
end
