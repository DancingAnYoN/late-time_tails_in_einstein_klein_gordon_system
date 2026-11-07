using Test
using MMIEF
const INT=MMIEF.Integrators

@testset "Generic RK4 and scalar shift" begin
    function solve(n)
        u=reshape(Float64[1.0],1,1)
        f=1.0
        w=INT.RK4Workspace(u)
        rhs! = (du,x,y)->(du .= -x; -2*y)
        h=1.0/n
        for _ in 1:n
            f,_=INT.rk4_step!(u,f,h,w;rhs! = rhs!)
        end
        return abs(u[1]-exp(-1)), abs(f-exp(-2))
    end
    e1,f1=solve(10)
    e2,f2=solve(20)
    @test e1/e2 > 12
    @test f1/f2 > 12

    # Earlier `2f2 + 2f3` was parsed as 200+2000, not RK4 weights.
    u=ones(2,5)
    w=INT.RK4Workspace(u)
    static_rhs! = (du,x,shift)->(fill!(du,0.0); 0.0)
    next_shift,_=INT.rk4_step!(u,0.0,0.1,w;rhs! = static_rhs!)
    @test next_shift == 0.0
    @test all(u .== 1.0)

    # Optional stage projections are called on all three intermediate stages
    # plus the accepted state; the pre-projection residual is measured first.
    calls=Ref(0)
    projection! = (x,shift)->(calls[]+=1; nothing)
    INT.rk4_step!(u,0.0,0.1,w;rhs! = static_rhs!,project! = projection!)
    @test calls[] == 4
end
