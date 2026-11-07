module BoundaryConditions
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..Geometry: fdot_horizon
export outer_projection!, outer_rhs!
"""Project only the incoming scalar characteristic at the distant outer boundary.
Uses the asymptotic outgoing spherical-wave condition, Eq. (4.11).
"""
function outer_projection!(u::Matrix{Float64}, f::Float64, p::Config, r::Vector{Float64})
    i = length(r)
    fd = fdot_horizon(u,f,p)
    q = u[i,IA]*(r[i]+f)*u[i,IB]
    den = 1+q
    den > 0 || error("Outer gauge denominator <= 0")
    beta = (fd+q)/den
    omb = 1-beta
    omb > 0 || error("Outer lapse factor <= 0")
    c = 1-2beta
    u[i,IPI] = -u[i,IPHI] - (fd+c)*u[i,IFIELD]/((r[i]+f)*omb)
    return nothing
end

"""Outer radial RHS following the original code's approximate outgoing conditions.
Precomputed derivatives are passed explicitly; this routine is not the source
of their discretization.
"""
function outer_rhs!(du,u,J,beta,fd,dAmetric,dB,dPhi,dBeta)
    i=size(u,1)
    Ji=J[i]; b=beta[i]; c=1-2*b
    du[i,IA] = -c*dAmetric[i]-(fd+c)*(u[i,IA]-1)/Ji
    du[i,IB] = -c*dB[i]-2*(fd+c)*u[i,IB]/Ji
    du[i,IPHI] = -c*dPhi[i]-(fd+c-2*Ji*dBeta[i])*u[i,IPHI]/Ji +
                  (fd+c+2*Ji*dBeta[i])*u[i,IFIELD]/Ji^2
    du[i,IPI] = -du[i,IPHI] # projection enforces actual Pi value
    return nothing
end
end # module BoundaryConditions
