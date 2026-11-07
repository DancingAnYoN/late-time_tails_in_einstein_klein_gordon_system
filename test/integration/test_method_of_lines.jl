using Test
using MMIEF

@testset "Method of Lines: periodic one-dimensional wave equation" begin
    # First-order 1D wave system p_t=q_x, q_t=p_x with p=q=sin(x+t).
    N=64
    dx=2π/N
    x=collect(0:N-1).*dx
    u=zeros(Float64,N,2)
    u[:,1].=sin.(x)
    u[:,2].=sin.(x)
    work=MMIEF.Integrators.RK4Workspace(u)
    function wave_rhs!(du, state, shift)
        @inbounds for i in 1:N
            prev2=mod1(i-2,N); prev1=mod1(i-1,N)
            next1=mod1(i+1,N); next2=mod1(i+2,N)
            du[i,1]=(state[prev2,2]-8*state[prev1,2]+
                      8*state[next1,2]-state[next2,2])/(12*dx)
            du[i,2]=(state[prev2,1]-8*state[prev1,1]+
                      8*state[next1,1]-state[next2,1])/(12*dx)
        end
        return 0.0
    end
    noop=(state,shift)->nothing
    model=MMIEF.MethodOfLines.MOLSystem(wave_rhs!,noop,noop,
                                        (state,shift)->0.0)
    t=0.0; f=0.0
    while t<1-1e-12
        dt=MMIEF.MethodOfLines.cfl_dt(dx,0.2,1.0,1-t)
        f,_=MMIEF.MethodOfLines.mol_step!(u,f,dt,work,model)
        t+=dt
    end
    expected=sin.(x.+t)
    @test maximum(abs.(u[:,1].-expected))<1e-4
    @test maximum(abs.(u[:,2].-expected))<1e-4
    @test f==0.0
end
