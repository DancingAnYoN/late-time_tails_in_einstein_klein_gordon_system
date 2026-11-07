module Evolution
using Printf
using ..Types: Config, IA, IB, IPHI, IPI, IFIELD
using ..Grid: radial_grid
using ..FiniteDifferences: deriv4!, ko6!
using ..Geometry: schwarzschild_fields, fdot_horizon
using ..Constraints: project_hamiltonian!, reconstruct_A!
using ..ScalarField: scalar_rhs!
using ..EinsteinEquations: gravity_rhs!
using ..InitialData: initial_data
using ..BoundaryConditions: outer_projection!, outer_rhs!
using ..Integrators: RK4Workspace
using ..MethodOfLines: MOLSystem, mol_step!, cfl_dt
using ..Observers: sample_fixed_J
using ..Diagnostics: diagnostics
using ..DataOutput: RunData, store_sample!, store_frame!
export run_evolution, vacuum_test, rhs!, Workspace, setup_vacuum_balance!,
       evolve_one_step!, EvolutionWorkspace
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

"""Coupled method-of-lines RHS. Physics formulas are implemented in their
own modules; only derivative assembly and coupling happen here.
"""
function rhs!(du::Matrix{Float64}, u::Matrix{Float64}, f::Float64,
              p::Config, w::Workspace; subtract_background::Bool=true)
    n=length(w.r); h=p.dr
    @views begin
        a=u[:,IA]; B=u[:,IB]; Phi=u[:,IPHI]; Pi=u[:,IPI]
        deriv4!(w.dB,B,h)
        deriv4!(w.dAmetric,a,h)
        deriv4!(w.dPhi,Phi,h)
        fd=fdot_horizon(u,f,p)
        fd < 1 || error(@sprintf("MMIEF lapse collapse (fdot=%.4g >= 1)",fd))
        @inbounds for i in 1:n
            Ji=w.r[i]+f; ai=a[i]; Bi=B[i]
            (isfinite(ai) && ai>0 && Ji>0) || error("Invalid geometry at index $i")
            w.J[i]=Ji
        end
        reconstruct_A!(w.A,u,w.J,w.dB)
        @inbounds for i in 1:n
            Ji=w.J[i]; ai=a[i]; Bi=B[i]
            q=ai*Ji*Bi; den=1+q
            den>0 || error("Gauge denominator <=0 at index $i")
            beta=(fd+q)/den
            w.beta[i]=beta
            w.aBeta[i]=ai*beta
            w.fFlux[i]=beta*Phi[i]+(1-beta)*Pi[i]
            w.pFlux[i]=Ji^2*(beta*Pi[i]+(1-beta)*Phi[i])
        end
        deriv4!(w.dBeta,w.beta,h)
        deriv4!(w.dABeta,w.aBeta,h)
        deriv4!(w.dFFlux,w.fFlux,h)
        deriv4!(w.dPFlux,w.pFlux,h)
        gravity_rhs!(du,u,w.J,w.A,w.beta,w.dB,w.dBeta,w.dABeta)
        scalar_rhs!(du,u,w.J,fd,w.dFFlux,w.dPFlux,w.fFlux)
        ko6!(du,u,h,p.ko_strength)
        outer_rhs!(du,u,w.J,w.beta,fd,w.dAmetric,w.dB,w.dPhi,w.dBeta)
        if subtract_background
            @inbounds for i in 1:n
                du[i,IA]-=w.vacuum_a_rhs[i]
                du[i,IB]-=w.vacuum_B_rhs[i]
            end
        end
        return fd
    end
end

struct EvolutionWorkspace
    rk4::RK4Workspace
    rhswork::Workspace
end
EvolutionWorkspace(u::Matrix{Float64},r::Vector{Float64}) =
    EvolutionWorkspace(RK4Workspace(u),Workspace(r))

"""Advance an RK4 stage with constraint/boundary projection callbacks."""
function evolve_one_step!(u::Matrix{Float64},f::Float64,dt::Float64,
                          p::Config,w::EvolutionWorkspace)
    rhs_eval! = (du,x,shift)->rhs!(du,x,shift,p,w.rhswork)
    project! = (x,shift)->project_hamiltonian!(x,shift,p,w.rhswork.r)
    boundary! = (x,shift)->outer_projection!(x,shift,p,w.rhswork.r)
    preprojection = (x,shift)->x[1,IA]*(w.rhswork.r[1]+shift)*x[1,IB]-1
    model=MOLSystem(rhs_eval!,project!,boundary!,preprojection)
    return mol_step!(u,f,dt,w.rk4,model)
end
"""Full coupled MMIEF evolution with constraint projection and well-balanced vacuum RHS.
A numerical research prototype; not yet validated as an accurate QNM / Price-law solver.
"""
function run_evolution(p::Config; verbose::Bool=true)
    r,u=initial_data(p)
    f=0.0
    outer_projection!(u,f,p,r)
    wk=EvolutionWorkspace(u,r)
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
        dt=cfl_dt(p.dr,p.cfl,maxspeed,p.tmax-t)
        f, preproj_error=evolve_one_step!(u,f,dt,p,wk)
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

end # module Evolution
