module Analysis
export fit_qnm, local_power_index
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

end # module Analysis
