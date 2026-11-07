module Diagnostics
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..FiniteDifferences: deriv4!
using ..Geometry: misner_sharp_mass
using ..Constraints: hamiltonian_residual
export diagnostics
function diagnostics(u::Matrix{Float64}, f::Float64, r::Vector{Float64},p::Config)
    n = length(r)
    da=zeros(n); dB=zeros(n); dphi=zeros(n)
    @views begin
        deriv4!(da,u[:,IA],p.dr)
        deriv4!(dB,u[:,IB],p.dr)
        deriv4!(dphi,u[:,IFIELD],p.dr)
    end
    sumH=0.0; sumDef=0.0; count=0
    @inbounds for i in 2:n-2
        J=r[i]+f
        a,B,Phi,Pi=u[i,IA],u[i,IB],u[i,IPHI],u[i,IPI]
        A=B+J*dB[i]-4π*J*Phi*Pi/a
        H=da[i]+(a^3-a)/(2J)+(a^3*J/2)*B*(2A+B)-2π*J*a*(Phi^2+Pi^2)
        sumH += H^2
        sumDef += (Phi-dphi[i])^2
        # Momentum residual identically zero by algebraic reconstruction.
        count += 1
    end
    jh=r[1]+f
    ah=u[1,IA]; bh=u[1,IB]
    outerJ=r[end]+f
    aend=u[end,IA]; bend=u[end,IB]
    madm=(outerJ/2)*(1-1/aend^2+(outerJ*bend)^2)
    horizon_error=ah*jh*bh-1
    return (ham_l2=sqrt(sumH/max(count,1)),
            scalar_constraint_l2=sqrt(sumDef/max(count,1)),
            horizon_error=horizon_error, M_h=jh/2, M_out=madm)
end

end # module Diagnostics
