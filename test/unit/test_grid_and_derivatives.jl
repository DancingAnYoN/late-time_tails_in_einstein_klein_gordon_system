using Test
using MMIEF
const FD=MMIEF.FiniteDifferences

@testset "Grid and finite differences" begin
    p=Config(rmax=80.0,dr=0.25,tmax=1.0)
    g=MMIEF.Grid.RadialGrid(p)
    @test g.r[1] == 2*p.M0
    @test g.r[end] ≈ p.rmax
    @test all(isapprox.(diff(g.r),p.dr;atol=1e-14))
    @test_throws ErrorException MMIEF.Grid.radial_grid(Config(dr=-0.1))
    function derivative_error(h)
        r=collect(2.0:h:12.0)
        v=sin.(r)
        dv=similar(v)
        FD.deriv4!(dv,v,h)
        return maximum(abs.(dv[3:end-2].-cos.(r[3:end-2])))
    end
    e1=derivative_error(0.2)
    e2=derivative_error(0.1)
    @test 12 < e1/e2 < 20
    r=collect(0.0:0.1:5.0)
    v=1.0 .+ r .+ 2 .*r.^2
    dv=similar(v)
    FD.deriv4!(dv,v,0.1)
    @test maximum(abs.(dv .- (1 .+ 4 .*r))) < 1e-9
    u=fill(1.0,length(r),5)
    du=zeros(size(u))
    FD.ko6!(du,u,0.1,0.06)
    @test maximum(abs,du)<1e-12
end
