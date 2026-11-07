using Test
using MMIEF
const GEO=MMIEF.Geometry
const CNS=MMIEF.Constraints

@testset "Exact geometrical identities" begin
    for J in (2.0,3.0,5.0,10.0,80.0)
        s=GEO.schwarzschild_fields(J,1.0)
        @test GEO.misner_sharp_mass(J,s.a,s.B) ≈ 1.0 atol=1e-12
        gauge=GEO.lapse_shift(s.a,J,s.B,0.0)
        @test gauge.beta ≈ s.beta atol=1e-12
        @test gauge.alpha ≈ s.a*(1-s.beta) atol=1e-12
        @test gauge.alpha > 0
    end
    @test maximum(abs,MMIEF.vacuum_test(dr=0.25,rmax=80.0)) < 1e-10
end

@testset "Constrained Gaussian initial geometry" begin
    p=Config(rmax=80.0,dr=0.25,tmax=1.0)
    r,u=MMIEF.initial_data(p)
    d=MMIEF.diagnostics(u,0.0,r,p)
    @test all(isfinite,u)
    @test d.M_out >= p.M0
    @test abs(d.horizon_error)<1e-12
    @test d.ham_l2 < 1e-3
    @test d.scalar_constraint_l2 < 1e-6
    u2=copy(u)
    CNS.project_hamiltonian!(u2,0.0,p,r)
    @test maximum(abs.(u2[:,1].-u[:,1])) < 1e-3
end
