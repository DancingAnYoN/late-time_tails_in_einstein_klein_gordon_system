module EinsteinEquations
using ..Types: IA, IB
export gravity_rhs!
"""The two evolved Einstein equations in the MMIEF gauge."""
function gravity_rhs!(du,u,J,A,beta,dB,dBeta,dABeta)
    @inbounds for i in eachindex(J)
        ai=u[i,IA]; Bi=u[i,IB]; Ai=A[i]; b=beta[i]; Ji=J[i]
        du[i,IA] = -ai^2*(1-b)*Ai + dABeta[i]
        du[i,IB] = b*dB[i] + ai*(1-b)*Bi*(Ai+2*Bi) +
                   (1-b)*(ai-1/ai)/Ji^2 + dBeta[i]/(ai*Ji)
    end
    return nothing
end
end # module EinsteinEquations
