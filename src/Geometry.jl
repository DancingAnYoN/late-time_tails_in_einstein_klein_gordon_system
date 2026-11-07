module Geometry
using ..Types: Config, IA, IB, IPHI, IPI
export schwarzschild_fields, schwarzschild_Bprime, fdot_horizon,
       lapse_shift, misner_sharp_mass, horizon_radius
"""Vacuum Schwarzschild values for f=0, useful as an exact stationary benchmark."""
function schwarzschild_fields(J::Float64, M::Float64)
    a = sqrt(1 + 2M/J)
    den = (J*(J+2M))^1.5
    B = 2M*(J+2M)/den
    A = -2M*(J+M)/den
    beta = 2M/(J+2M)
    return (a=a, B=B, A=A, beta=beta)
end

"""Analytic radial derivative of the Schwarzschild B component."""
@inline schwarzschild_Bprime(r::Float64,M::Float64) =
    -2M*(2r+3M)/(r^(5/2)*(r+2M)^(3/2))

"""Gauge reconstruction using J=r+f and incoming null slicing."""
function lapse_shift(a::Real,J::Real,B::Real,fdot::Real)
    q = a*J*B
    denominator = 1+q
    denominator>0 || error("Invalid gauge denominator (1+aJB)<=0")
    beta=(fdot+q)/denominator
    alpha=a*(1-beta)
    alpha>0 || error("Nonpositive lapse")
    return (alpha=alpha,beta=beta)
end
@inline misner_sharp_mass(J,a,B) = J/2*(1-1/a^2+(J*B)^2)
@inline horizon_radius(p::Config,f::Real)=2*p.M0+f
"""Horizon-locked coordinate growth rate, Eq. (2.20)."""
@inline function fdot_horizon(u::Matrix{Float64}, f::Float64, p::Config)
    J = 2*p.M0+f
    return 4π*(J/u[1,IA])^2*(u[1,IPHI]+u[1,IPI])^2
end

end # module Geometry
