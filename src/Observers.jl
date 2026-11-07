module Observers
using Printf
using ..Types: Config, IFIELD
using ..Geometry: fdot_horizon
export sample_fixed_J
"""Cubic Lagrange interpolation at fixed physical areal radius J_obs.

A fixed-J observer must remain outside the horizon and within the shifted
physical grid. If not, throw an error with the actual geometry and time;
silently clamping the interpolation point would produce a false waveform.
"""
function sample_fixed_J(u::Matrix{Float64}, r::Vector{Float64}, f::Float64,
                        J_obs::Float64, p::Config; time::Float64=NaN)
    x = J_obs-f
    if !(isfinite(x) && r[1] <= x <= r[end])
        fd = isfinite(f) ? fdot_horizon(u,f,p) : NaN
        error(@sprintf("""
        Fixed-J observer is outside the computational grid at t/M0=%.9g.
        Requested J_obs/M0=%.9g; f/M0=%.9g; observer coordinate r_obs/M0=%.9g.
        Coordinate grid r/M0=[%.9g, %.9g]; physical J/M0=[%.9g, %.9g].
        Horizon mass M_h/M0=%.9g; fdot=%.9g; initial pulse amplitude=%.9g.
        If J_obs <= J_h, the observer is inside/at the tracked horizon.
        For a weak pulse this suggests numerical instability; inspect the
        mass-budget diagnostic, lower the CFL number, and test --vacuum.
        """, time/p.M0, J_obs/p.M0, f/p.M0, x/p.M0,
        r[1]/p.M0, r[end]/p.M0, (r[1]+f)/p.M0, (r[end]+f)/p.M0,
        (r[1]+f)/(2*p.M0), fd, p.pulse_amplitude))
    end
    n = length(r)
    i0 = clamp(floor(Int,(x-r[1])/p.dr),1,n-3)
    i0 = max(1,i0)  # 4-point stencil i0:i0+3
    val = 0.0
    @inbounds for j in i0:i0+3
        w = 1.0
        for k in i0:i0+3
            if k != j
                w *= (x-r[k])/(r[j]-r[k])
            end
        end
        val += w*u[j,IFIELD]
    end
    return val
end

end # module Observers
