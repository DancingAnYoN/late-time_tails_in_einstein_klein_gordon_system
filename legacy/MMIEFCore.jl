module MMIEFCore

using Printf

export Config, RunData, initial_data, diagnostics, run_evolution, save_csv,
       fit_qnm, local_power_index, schwarzschild_fields, vacuum_test,
       deriv4!, radial_grid, gaussian

# EKG in minimally modified ingoing Eddington-Finkelstein (MMIEF) coordinates.
# Metric: ds² = -α²dt² + a²(dr+βdt)² + J²dΩ², J=r+f(t).
# Units G=c=1, signature (-+++), K_ij = -1/2 L_n γ_ij.
# Physical field variables Φ=∂r φ and Π=a/α (∂t φ - β Φ).

Base.@kwdef struct Config
    M0::Float64 = 1.0
    rmax::Float64 = 600.0
    dr::Float64 = 0.25
    cfl::Float64 = 0.25
    tmax::Float64 = 400.0
    pulse_r0::Float64 = 20.0
    pulse_sigma::Float64 = 2.0
    pulse_amplitude::Float64 = 1e-4
    cutoff_start::Float64 = 6.0
    cutoff_end::Float64 = 9.0
    ko_strength::Float64 = 0.06
    output_dt::Float64 = 0.5
    snapshot_dt::Float64 = 2.0
    profile_stride::Int = 4
    observers::Vector{Float64} = [10.0, 30.0, 50.0]
    # Fail early when the horizon mass exceeds the initial outer mass by
    # far more than can be attributed to truncation error. A weak, isolated
    # ingoing pulse cannot increase the BH mass by many solar/geometric units.
    check_mass_budget::Bool = true
    mass_budget_fraction::Float64 = 0.01
end

const IA, IB, IPHI, IPI, IFIELD = 1, 2, 3, 4, 5

function check_config(p::Config)
    rh = 2*p.M0
    p.M0 > 0 || error("M0 must be positive")
    p.dr > 0 && p.rmax > rh + 12*p.dr || error("Invalid radial grid")
    nsteps = (p.rmax-rh)/p.dr
    abs(nsteps-round(nsteps)) < 1e-8 ||
        error("(rmax-2M0)/dr must be an integer to place the outer boundary exactly")
    0 < p.cfl <= 0.5 || error("Choose 0 < cfl <= 0.5 for this implementation")
    p.tmax > 0 && p.output_dt > 0 && p.snapshot_dt > 0 || error("Invalid time settings")
    p.pulse_sigma > 0 && p.pulse_amplitude >= 0 || error("Invalid pulse")
    rh < p.cutoff_start < p.cutoff_end < p.pulse_r0 || error("Cutoff must lie between horizon and pulse")
    all(J -> rh < J < p.rmax, p.observers) || error("Observers must lie in domain")
    0 <= p.ko_strength <= 0.15 || error("Use a small KO strength")
    p.profile_stride >= 1 || error("profile_stride must be positive")
    p.mass_budget_fraction >= 0 || error("mass_budget_fraction must be nonnegative")
end

function radial_grid(p::Config)
    check_config(p)
    N = floor(Int, (p.rmax - 2*p.M0) / p.dr) + 1
    return [2*p.M0 + (i - 1) * p.dr for i in 1:N]
end

"""Fourth-order accurate 5-point derivative, including 4th-order one-sided endpoints."""
function deriv4!(d::AbstractVector{Float64}, y::AbstractVector{Float64}, h::Float64)
    n = length(y)
    n >= 7 || error("At least 7 grid points are required")
    z = 1.0 / (12h)
    @inbounds begin
        d[1] = z*(-25y[1] + 48y[2] - 36y[3] + 16y[4] - 3y[5])
        d[2] = z*(-3y[1] - 10y[2] + 18y[3] - 6y[4] + y[5])
        for i in 3:n-2
            d[i] = z*(y[i-2] - 8y[i-1] + 8y[i+1] - y[i+2])
        end
        d[n-1] = z*(3y[n] + 10y[n-1] - 18y[n-2] + 6y[n-3] - y[n-4])
        d[n] = z*(25y[n] - 48y[n-1] + 36y[n-2] - 16y[n-3] + 3y[n-4])
    end
    return d
end

"""Sixth-derivative KO filter: negative semidefinite for Fourier modes."""
function ko6!(du::Matrix{Float64}, u::Matrix{Float64}, h::Float64, strength::Float64)
    strength == 0 && return
    n = size(u, 1)
    weight = strength / (64h)
    @inbounds for col in 1:size(u, 2), i in 4:n-3
        du[i,col] += weight*(u[i-3,col] - 6u[i-2,col] + 15u[i-1,col]
                            - 20u[i,col] + 15u[i+1,col] - 6u[i+2,col] + u[i+3,col])
    end
end

"""Gaussian pulse multiplied by a C-infinity cutoff; return (phi, dphi/dr)."""
function gaussian(r::Float64, p::Config)
    q = (r - p.pulse_r0) / p.pulse_sigma
    g = p.pulse_amplitude * exp(-q*q)
    gp = -2*q/p.pulse_sigma * g
    if r <= p.cutoff_start
        return 0.0, 0.0
    elseif r >= p.cutoff_end
        return g, gp
    end
    w = p.cutoff_end-p.cutoff_start
    x = (r-p.cutoff_start)/w
    # Standard smooth bump: chi=exp(-1/x)/(exp(-1/x)+exp(-1/(1-x))).
    e1 = exp(-1/x)
    e2 = exp(-1/(1-x))
    chi = e1 / (e1+e2)
    chip = chi*(1-chi)*(1/x^2 + 1/(1-x)^2)/w
    return chi*g, chi*gp+chip*g
end

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

"""
Initial constraint ODEs for delta_m=m-M0 and delta_B=B-B_Schw(M0),
with the Marsa--Choptuik Schwarzschild-inspired A(m,r) ansatz. Evolving
only differences exactly preserves the vacuum initial-data background.
Ingoing scalar data use ∂t(r φ)=∂r(r φ), f(0)=fdot(0)=0.
"""
function initial_ode(r::Float64, dm::Float64, dB::Float64, p::Config)
    m = p.M0 + dm
    B = schwarzschild_fields(r,p.M0).B + dB
    inva2 = 1 - 2m/r + (r*B)^2
    (isfinite(inva2) && inva2 > 0) || error("Non-positive initial inverse a² at r=$r")
    a = inv(sqrt(inva2))
    A = -2m*(r+m)/(r*(r+2m))^1.5
    phi, Phi = gaussian(r,p)
    q = a*r*B
    (1+q > 0) || error("Initial gauge denominator failed")
    # Pi = Phi + phi/(r(1-beta)); 1/(1-beta) = 1+q at fdot=0.
    Pi = Phi + (1+q)*phi/r
    mp = 4π*r*r*((Phi*Phi+Pi*Pi)/(2a*a) + (r*B/a)*Phi*Pi)
    Bp = (A-B)/r + 4π*Phi*Pi/a - schwarzschild_Bprime(r,p.M0)
    return mp,Bp
end

function initial_data(p::Config)
    r = radial_grid(p)
    n = length(r)
    u = Matrix{Float64}(undef, n, 5)
    dm = Vector{Float64}(undef, n)
    dB = Vector{Float64}(undef, n)
    dm[1] = 0.0
    dB[1] = 0.0 # a_h J_h B_h=1 for the Schwarzschild reference
    h = p.dr
    @inbounds for i in 1:n-1
        x = r[i]
        m0, b0 = dm[i], dB[i]
        k1m,k1b = initial_ode(x,m0,b0,p)
        k2m,k2b = initial_ode(x+h/2,m0+h*k1m/2,b0+h*k1b/2,p)
        k3m,k3b = initial_ode(x+h/2,m0+h*k2m/2,b0+h*k2b/2,p)
        k4m,k4b = initial_ode(x+h,m0+h*k3m,b0+h*k3b,p)
        dm[i+1] = m0 + h*(k1m+2k2m+2k3m+k4m)/6
        dB[i+1] = b0 + h*(k1b+2k2b+2k3b+k4b)/6
    end
    @inbounds for i in 1:n
        J = r[i]
        mi = p.M0 + dm[i]
        Bi = schwarzschild_fields(J,p.M0).B + dB[i]
        a = inv(sqrt(1 - 2mi/J + (J*Bi)^2))
        phi,Phi = gaussian(J,p)
        Pi = Phi + (1+a*J*Bi)*phi/J
        u[i,IA] = a
        u[i,IB] = Bi
        u[i,IPHI] = Phi
        u[i,IPI] = Pi
        u[i,IFIELD] = phi
    end
    return r,u
end

"""Radial Misner--Sharp mass equation for constraint projection.
This uses J'=1; 1/a^2=1-2m/J+(JB)^2.
"""
@inline function mass_rhs(J::Float64, m::Float64, B::Float64,
                          Phi::Float64, Pi::Float64)
    ainvsq = 1 - 2m/J + (J*B)^2
    (isfinite(ainvsq) && ainvsq > 0) || error("Hamiltonian projection: nonpositive 1/a² at J=$J")
    ainv = sqrt(ainvsq)
    return 4π*J^2 * (0.5*(Phi^2+Pi^2)*ainvsq + J*B*ainv*Phi*Pi)
end

"""
Project the Hamiltonian constraint on each RK4 stage by integrating the
Misner--Sharp mass outward from m(J_h)=J_h/2, then reconstructing a.
4th-order radial RK with cubic midpoint interpolation, in vacuum exact.
This is a constrained-RK4 modification of the paper's original method.
"""
function project_hamiltonian!(u::Matrix{Float64}, f::Float64,
                              p::Config, r::Vector{Float64})
    n=length(r)
    m=(r[1]+f)/2
    h=p.dr
    @inbounds for i in 1:n
        J=r[i]+f
        if i > 1
            j=i-1
            J0=r[j]+f
            B0, P0, Pi0 = u[j,IB], u[j,IPHI], u[j,IPI]
            B1, P1, Pi1 = u[i,IB], u[i,IPHI], u[i,IPI]
            if 2 <= j <= n-2
                # Fourth-order midpoint interpolation of stage fields.
                Bmid=(-u[j-1,IB]+9B0+9B1-u[j+2,IB])/16
                Pmid=(-u[j-1,IPHI]+9P0+9P1-u[j+2,IPHI])/16
                Pimid=(-u[j-1,IPI]+9Pi0+9Pi1-u[j+2,IPI])/16
            else
                Bmid=(B0+B1)/2
                Pmid=(P0+P1)/2
                Pimid=(Pi0+Pi1)/2
            end
            k1=mass_rhs(J0,m,B0,P0,Pi0)
            k2=mass_rhs(J0+h/2,m+h*k1/2,Bmid,Pmid,Pimid)
            k3=mass_rhs(J0+h/2,m+h*k2/2,Bmid,Pmid,Pimid)
            k4=mass_rhs(J,m+h*k3,B1,P1,Pi1)
            m+=h*(k1+2k2+2k3+k4)/6
        end
        inva2=1-2m/J+(J*u[i,IB])^2
        (isfinite(inva2) && inva2 > 0) || error("Hamiltonian projection failed at J=$J")
        u[i,IA]=inv(sqrt(inva2))
    end
    return nothing
end

mutable struct Workspace
    r::Vector{Float64}
    J::Vector{Float64}
    A::Vector{Float64}
    beta::Vector{Float64}
    dAmetric::Vector{Float64}
    dB::Vector{Float64}
    dPhi::Vector{Float64}
    dBeta::Vector{Float64}
    dABeta::Vector{Float64}
    dFFlux::Vector{Float64}
    dPFlux::Vector{Float64}
    aBeta::Vector{Float64}
    fFlux::Vector{Float64}
    pFlux::Vector{Float64}
    vacuum_a_rhs::Vector{Float64}
    vacuum_B_rhs::Vector{Float64}
end

function Workspace(r::Vector{Float64})
    a = [zeros(length(r)) for _ in 1:15]
    return Workspace(r, a...)
end

"""Precompute the exact Schwarzschild discrete RHS as a well-balanced
truncation-error correction, giving stationary vacuum to roundoff.
This subtracts only O(dr^4) terms (plus boundary/dissipation errors).
"""
function setup_vacuum_balance!(p::Config, w::Workspace)
    n=length(w.r)
    uv=zeros(n,5); dv=similar(uv)
    for i in 1:n
        s=schwarzschild_fields(w.r[i],p.M0)
        uv[i,IA]=s.a
        uv[i,IB]=s.B
    end
    rhs!(dv,uv,0.0,p,w;subtract_background=false)
    @views w.vacuum_a_rhs .= dv[:,IA]
    @views w.vacuum_B_rhs .= dv[:,IB]
    return nothing
end

"""Horizon-locked coordinate growth rate, Eq. (2.20)."""
@inline function fdot_horizon(u::Matrix{Float64}, f::Float64, p::Config)
    J = 2*p.M0+f
    return 4π*(J/u[1,IA])^2*(u[1,IPHI]+u[1,IPI])^2
end

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

"""Evaluate the fully coupled five-field RHS plus scalar reconstruction.
The 5th column is phi; the sixth evolution variable f is scalar.
"""
function rhs!(du::Matrix{Float64}, u::Matrix{Float64}, f::Float64,
              p::Config, w::Workspace; subtract_background::Bool=true)
    n = length(w.r)
    h = p.dr
    @views begin
        a=u[:,IA]; B=u[:,IB]; Phi=u[:,IPHI]; Pi=u[:,IPI]; phi=u[:,IFIELD]
        deriv4!(w.dB,B,h)
        deriv4!(w.dAmetric,a,h)
        deriv4!(w.dPhi,Phi,h)
        fd = fdot_horizon(u,f,p)
        fd < 1 || error(@sprintf("MMIEF lapse collapse (fdot=%.4g >= 1)",fd))
        @inbounds for i in 1:n
            J = w.r[i]+f
            w.J[i] = J
            ai = a[i]; Bi = B[i]
            (isfinite(ai) && ai > 0 && J > 0) || error("Invalid geometry at index $i")
            w.A[i] = Bi+J*w.dB[i]-4π*J*Phi[i]*Pi[i]/ai
            q = ai*J*Bi
            den = 1+q
            den > 0 || error("Gauge denominator <= 0 at grid index $i")
            beta = (fd+q)/den
            w.beta[i] = beta
            w.aBeta[i] = ai*beta
            w.fFlux[i] = beta*Phi[i]+(1-beta)*Pi[i]
            w.pFlux[i] = J*J*(beta*Pi[i]+(1-beta)*Phi[i])
        end
        deriv4!(w.dBeta,w.beta,h)
        deriv4!(w.dABeta,w.aBeta,h)
        deriv4!(w.dFFlux,w.fFlux,h)
        deriv4!(w.dPFlux,w.pFlux,h)
        @inbounds for i in 1:n
            J=w.J[i]; ai=a[i]; Bi=B[i]; Ai=w.A[i]; beta=w.beta[i]
            du[i,IA] = -ai*ai*(1-beta)*Ai+w.dABeta[i]
            du[i,IB] = beta*w.dB[i]+ai*(1-beta)*Bi*(Ai+2Bi) +
                       (1-beta)*(ai-1/ai)/(J*J)+w.dBeta[i]/(ai*J)
            du[i,IPHI] = w.dFFlux[i]
            du[i,IPI] = w.dPFlux[i]/(J*J)-2fd*Pi[i]/J
            du[i,IFIELD] = w.fFlux[i]
        end
        ko6!(du,u,h,p.ko_strength)
        # Far-boundary approximate outgoing relations (Marsa--Choptuik
        # Eqs. 4.10, 4.22, 4.23). Source-free and imperfect at finite r.
        i = n
        J=w.J[i]; beta=w.beta[i]; c=1-2beta
        du[i,IA] = -c*w.dAmetric[i]-(fd+c)*(a[i]-1)/J
        du[i,IB] = -c*w.dB[i]-2(fd+c)*B[i]/J
        du[i,IPHI] = -c*w.dPhi[i] - (fd+c-2J*w.dBeta[i])*Phi[i]/J +
                      (fd+c+2J*w.dBeta[i])*phi[i]/(J*J)
        # Pi at the last point is algebraically projected after each RK stage.
        du[i,IPI] = -du[i,IPHI]
        if subtract_background
            @inbounds for j in 1:n
                du[j,IA] -= w.vacuum_a_rhs[j]
                du[j,IB] -= w.vacuum_B_rhs[j]
            end
        end
        return fd
    end
end

struct RKWorkspace
    k1::Matrix{Float64}
    k2::Matrix{Float64}
    k3::Matrix{Float64}
    k4::Matrix{Float64}
    temp::Matrix{Float64}
    rhswork::Workspace
end
RKWorkspace(u::Matrix{Float64},r::Vector{Float64}) =
    RKWorkspace(similar(u),similar(u),similar(u),similar(u),similar(u),Workspace(r))

function rk4_step!(u::Matrix{Float64}, f::Float64, dt::Float64,
                   p::Config, wk::RKWorkspace)
    k1,k2,k3,k4,tmp = wk.k1,wk.k2,wk.k3,wk.k4,wk.temp
    f1 = rhs!(k1,u,f,p,wk.rhswork)
    @inbounds for i in eachindex(u)
        tmp[i] = u[i]+dt*k1[i]/2
    end
    f2arg = f+dt*f1/2
    project_hamiltonian!(tmp,f2arg,p,wk.rhswork.r)
    outer_projection!(tmp,f2arg,p,wk.rhswork.r)
    f2 = rhs!(k2,tmp,f2arg,p,wk.rhswork)
    @inbounds for i in eachindex(u)
        tmp[i] = u[i]+dt*k2[i]/2
    end
    f3arg = f+dt*f2/2
    project_hamiltonian!(tmp,f3arg,p,wk.rhswork.r)
    outer_projection!(tmp,f3arg,p,wk.rhswork.r)
    f3 = rhs!(k3,tmp,f3arg,p,wk.rhswork)
    @inbounds for i in eachindex(u)
        tmp[i] = u[i]+dt*k3[i]
    end
    f4arg = f+dt*f3
    project_hamiltonian!(tmp,f4arg,p,wk.rhswork.r)
    outer_projection!(tmp,f4arg,p,wk.rhswork.r)
    f4 = rhs!(k4,tmp,f4arg,p,wk.rhswork)
    @inbounds for i in eachindex(u)
        u[i] += dt*(k1[i]+2k2[i]+2k3[i]+k4[i])/6
    end
    # IMPORTANT: write multiplication explicitly! In Julia, `2f2` and `2f3`
    # are Float32 literals (200 and 2000), NOT 2*f2 and 2*f3.
    fnew = f + dt*(f1 + 2*f2 + 2*f3 + f4)/6
    # Pre-projection horizon drift is an independent numerical diagnostic;
    # post-projection residual is zero by construction.
    preprojection_error = u[1,IA]*(wk.rhswork.r[1]+fnew)*u[1,IB]-1
    project_hamiltonian!(u,fnew,p,wk.rhswork.r)
    outer_projection!(u,fnew,p,wk.rhswork.r)
    return fnew, preprojection_error
end

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

mutable struct RunData
    settings::Config
    r::Vector{Float64}
    times::Vector{Float64}
    signals::Vector{Vector{Float64}}
    shifts::Vector{Float64}
    masses_h::Vector{Float64}
    masses_out::Vector{Float64}
    ham_l2::Vector{Float64}
    scalar_constraint_l2::Vector{Float64}
    horizon_error::Vector{Float64}
    horizon_preprojection::Vector{Float64}
    frame_times::Vector{Float64}
    frame_shift::Vector{Float64}
    frame_profiles::Vector{Vector{Float32}}
    final_state::Matrix{Float64}
end

function RunData(p::Config,r::Vector{Float64},u::Matrix{Float64})
    return RunData(p,r,Float64[],[Float64[] for _ in p.observers],Float64[],
                   Float64[],Float64[],Float64[],Float64[],Float64[],
                   Float64[],Float64[],Float64[],Vector{Float32}[],u)
end

function store_sample!(data::RunData,t::Float64,u::Matrix{Float64},f::Float64,
                       preprojection_error::Float64=0.0)
    # Sample all observers before modifying RunData, so a failed interpolation
    # cannot leave its arrays with inconsistent lengths.
    wave = [sample_fixed_J(u,data.r,f,J,data.settings;time=t)
            for J in data.settings.observers]
    d=diagnostics(u,f,data.r,data.settings)
    push!(data.times,t)
    push!(data.shifts,f)
    push!(data.masses_h,d.M_h)
    push!(data.masses_out,d.M_out)
    push!(data.ham_l2,d.ham_l2)
    push!(data.scalar_constraint_l2,d.scalar_constraint_l2)
    push!(data.horizon_error,d.horizon_error)
    push!(data.horizon_preprojection,preprojection_error)
    for k in eachindex(data.settings.observers)
        push!(data.signals[k],wave[k])
    end
end

function store_frame!(data::RunData,t::Float64,u::Matrix{Float64},f::Float64)
    push!(data.frame_times,t)
    push!(data.frame_shift,f)
    stride=data.settings.profile_stride
    push!(data.frame_profiles, Float32[
        (data.r[i]+f)*u[i,IFIELD] for i in 1:stride:length(data.r)])
end

"""Full coupled MMIEF evolution with constraint projection and well-balanced vacuum RHS.
A numerical research prototype; not yet validated as an accurate QNM / Price-law solver.
"""
function run_evolution(p::Config; verbose::Bool=true)
    r,u=initial_data(p)
    f=0.0
    outer_projection!(u,f,p,r)
    wk=RKWorkspace(u,r)
    setup_vacuum_balance!(p,wk.rhswork)
    # Start on the Hamiltonian constraint manifold.
    project_hamiltonian!(u,f,p,r)
    outer_projection!(u,f,p,r)
    initial_outer_mass=diagnostics(u,f,r,p).M_out
    # The tolerance here is deliberately generous: it is a fail-fast check
    # against catastrophic unphysical growth, NOT a precision error estimate.
    mass_budget_tol=max(p.mass_budget_fraction*p.M0,
                        0.1*max(0.0,initial_outer_mass-p.M0))
    data=RunData(p,r,u)
    store_sample!(data,0.0,u,f)
    store_frame!(data,0.0,u,f)
    if verbose
        d=diagnostics(u,f,r,p)
        @printf("Initial: N=%d, M_outer=%.10f, ||H||2=%.3e, |aJB-1|_h=%.3e\n",
                length(r),d.M_out,d.ham_l2,abs(d.horizon_error))
    end
    t=0.0; step=0
    next_sample=p.output_dt
    next_frame=p.snapshot_dt
    next_report=25*p.M0
    while t < p.tmax-1e-13
        # The radial scalar characteristic speeds are -1 and (1-2beta).
        # Recompute the maximum characteristic speed in the current geometry.
        # The gravitational principal part also needs separate stability checks.
        fd=fdot_horizon(u,f,p)
        maxspeed=1.0
        @inbounds for i in eachindex(r)
            q=u[i,IA]*(r[i]+f)*u[i,IB]
            b=(fd+q)/(1+q)
            maxspeed=max(maxspeed,abs(1-2b))
        end
        dt=min(p.cfl*p.dr/maxspeed, p.tmax-t)
        f, preproj_error=rk4_step!(u,f,dt,p,wk)
        t+=dt
        step+=1
        # Catch a runaway horizon *before* fixed-J interpolation fails.
        # The total horizon mass should not exceed the initial mass at the
        # outer boundary in this isolated, outgoing-boundary setup.
        if p.check_mass_budget
            Mh=(r[1]+f)/2
            if !(isfinite(Mh) && Mh <= initial_outer_mass+mass_budget_tol)
                error(@sprintf("""
                Unphysical horizon-mass growth detected at t/M0=%.9g (step %d).
                M_h/M0=%.9g, initial outer mass/M0=%.9g, f/M0=%.9g,
                fdot=%.9g, tolerance/M0=%.9g.
                The weak pulse cannot physically grow the black hole this much.
                Investigate numerical stability, constraints, and boundaries;
                do not fix this by clamping observer coordinates.
                """, t/p.M0,step,Mh/p.M0,initial_outer_mass/p.M0,
                f/p.M0,fdot_horizon(u,f,p),mass_budget_tol/p.M0))
            end
        end
        # Every observer must remain within the physical domain.
        # At f >= J_obs-r_h, that observer is behind the inner boundary.
        for J_obs in p.observers
            sample_x = J_obs-f
            if !(isfinite(sample_x) && r[1] <= sample_x <= r[end])
                # Use the same informative message as the sampling function.
                sample_fixed_J(u,r,f,J_obs,p;time=t)
            end
        end
        if step%20==0
            @views av=u[:,IA]
            if !all(isfinite,u) || !isfinite(f) || minimum(av)<=0
                error("NaN / Inf / nonpositive a at t=$t; reduce amplitude or CFL")
            end
            fd=fdot_horizon(u,f,p)
            fd < 1 || error("MMIEF gauge has collapsed: fdot >= 1")
        end
        if t+1e-12>=next_sample || t>=p.tmax-1e-12
            store_sample!(data,t,u,f,preproj_error)
            next_sample+=p.output_dt
        end
        if t+1e-12>=next_frame || t>=p.tmax-1e-12
            store_frame!(data,t,u,f)
            next_frame+=p.snapshot_dt
        end
        if verbose && t+1e-12>=next_report
            d=diagnostics(u,f,r,p)
            @printf("t/M0=%7.1f  M_h=%1.8f  f/M0=%9.3e  fdot=%9.2e  H_L2=%8.2e  horizon_err=%8.2e\n",
                    t/p.M0,d.M_h,f/p.M0,fdot_horizon(u,f,p),d.ham_l2,abs(d.horizon_error))
            next_report+=25*p.M0
        end
    end
    data.final_state=copy(u)
    return data
end

function save_csv(data::RunData,filename::AbstractString)
    open(filename,"w") do io
        cols=["t","T_asymptotic","f","M_h","M_outer","hamiltonian_L2","scalar_definition_L2","horizon_error",
              "horizon_preprojection_error"]
        for J in data.settings.observers
            push!(cols,@sprintf("phi_J%.2f",J))
        end
        println(io,join(cols,","))
        for i in eachindex(data.times)
            vals=[data.times[i],data.times[i]-data.shifts[i],data.shifts[i],data.masses_h[i],data.masses_out[i],
                  data.ham_l2[i],data.scalar_constraint_l2[i],
                  data.horizon_error[i],data.horizon_preprojection[i]]
            append!(vals,[sig[i] for sig in data.signals])
            println(io,join([@sprintf("%.16e",x) for x in vals],","))
        end
    end
    return filename
end

"""Exploratory damped-sinusoid fit by a (gamma,omega) grid search.
Not a robust multi-overtone or Prony/QNM estimator: inspect window dependence.
"""
function fit_qnm(t::Vector{Float64}, y::Vector{Float64};
                 tlo::Float64=45.0, thi::Float64=100.0,
                 mass::Float64=1.0)
    inds=findall(i -> tlo<=t[i]<=thi, eachindex(t))
    length(inds)>=25 || return nothing
    tt=[t[i]-tlo for i in inds]
    yy=y[inds]
    best=(Inf, 0.0,0.0,0.0,0.0)
    # Dimensionless search in M omega; similar scalar l=0 Schwarzschild domain.
    for (gmin,gmax,wmin,wmax,nsteps) in ((0.045,0.16,0.055,0.18,30),
                                          (NaN,NaN,NaN,NaN,23))
        if isnan(gmin)
            _,g0,w0,_,_=best
            dg=0.115/(30-1)*2
            dw=0.125/(30-1)*2
            gmin=max(0.02,g0-2dg); gmax=g0+2dg
            wmin=max(0.02,w0-2dw); wmax=w0+2dw
        end
        for gammaM in range(gmin,gmax,length=nsteps), omegaM in range(wmin,wmax,length=nsteps)
            gamma=gammaM/mass; omega=omegaM/mass
            scc=0.0; sss=0.0; scs=0.0; syc=0.0; sys=0.0; syy=0.0
            for k in eachindex(tt)
                decay=exp(-gamma*tt[k])
                C=decay*cos(omega*tt[k]); S=decay*sin(omega*tt[k]); yk=yy[k]
                scc+=C*C; sss+=S*S; scs+=C*S
                syc+=yk*C; sys+=yk*S; syy+=yk*yk
            end
            det=scc*sss-scs*scs
            det<=1e-25 && continue
            ampC=(syc*sss-sys*scs)/det
            ampS=(sys*scc-syc*scs)/det
            residual=max(0.0,syy-(ampC*syc+ampS*sys))
            if residual<best[1]
                best=(residual,gammaM,omegaM,ampC,ampS)
            end
        end
    end
    isfinite(best[1]) || return nothing
    res,gM,wM,C,S=best
    return (gammaM=gM, omegaM=wM, C=C, S=S, tlo=tlo, thi=thi,
            rel_error=sqrt(res/max(sum(abs2,yy),eps(Float64))), mass=mass)
end

"""Estimate local -dlog(envelope)/dlog(t) with a running peak envelope.
Window parameters are in coordinate time units. Avoid zero crossings by using
maximum |phi| over ±peak_halfwidth, then regress over ±fit_halfwidth.
"""
function local_power_index(t::Vector{Float64},y::Vector{Float64};
                           peak_halfwidth::Float64=15.0,
                           fit_halfwidth::Float64=25.0,
                           tstart::Float64=100.0)
    n=length(t)
    env=zeros(Float64,n)
    for i in 1:n
        l=searchsortedfirst(t,t[i]-peak_halfwidth)
        h=searchsortedlast(t,t[i]+peak_halfwidth)
        vmax=0.0
        @inbounds for j in l:h
            vmax=max(vmax,abs(y[j]))
        end
        env[i]=vmax
    end
    out=fill(NaN,n)
    flooramp=max(maximum(env)*1e-12,1e-29)
    for i in 1:n
        if t[i] < tstart || env[i] < flooramp
            continue
        end
        lo=searchsortedfirst(t,max(tstart,t[i]-fit_halfwidth))
        hi=searchsortedlast(t,t[i]+fit_halfwidth)
        hi-lo<12 && continue
        sx=0.0; sy=0.0; sxx=0.0; sxy=0.0; count=0
        for j in lo:hi
            if t[j]>0 && env[j]>=flooramp
                x=log(t[j]); z=log(env[j])
                sx+=x; sy+=z; sxx+=x*x; sxy+=x*z; count+=1
            end
        end
        den=count*sxx-sx*sx
        if count>12 && den>1e-12
            out[i]=-(count*sxy-sx*sy)/den
        end
    end
    return out
end

"""Zero-field Schwarzschild invariance test: reports discrete RHS residuals."""
function vacuum_test(; dr::Float64=0.25,rmax::Float64=100.0)
    p=Config(dr=dr,rmax=rmax,tmax=1.0,pulse_amplitude=0.0)
    r=radial_grid(p)
    u=zeros(Float64,length(r),5)
    for i in eachindex(r)
        s=schwarzschild_fields(r[i],p.M0)
        u[i,IA]=s.a
        u[i,IB]=s.B
    end
    du=similar(u)
    w=Workspace(r)
    setup_vacuum_balance!(p,w)
    rhs!(du,u,0.0,p,w)
    return [maximum(abs, @view(du[:,c])) for c in 1:5]
end

end # module
