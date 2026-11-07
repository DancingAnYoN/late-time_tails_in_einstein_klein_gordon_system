module Dashboard
using GLMakie
using Printf
using MMIEF: RunData, fit_qnm, local_power_index

export make_dashboard

"""
GLMakie visualization for the MMIEF evolution.

The simulation is run beforehand; animation uses stored snapshots only.
The QNM fit and tail indicators are exploratory, not convergence-certified.
"""

@inline p2(x::Real, y::Real) = Point2f(Float32(x), Float32(y))

function trace_points(t, y, n, M)
    return Point2f[p2(t[i] / M, y[i]) for i in 1:n]
end

function log_points(t, y, n, M)
    points = Point2f[]
    for i in 1:n
        if t[i] > 0.0 && abs(y[i]) > 1e-28
            push!(points, p2(log10(t[i] / M), log10(abs(y[i]))))
        end
    end
    return points
end

function residual_points(t, y, n, M)
    return Point2f[p2(t[i] / M, log10(max(abs(y[i]), 1e-18))) for i in 1:n]
end

function lpi_points(t, y, n, M)
    points = Point2f[]
    for i in 1:n
        if isfinite(y[i]) && -4.0 <= y[i] <= 12.0
            push!(points, p2(t[i] / M, y[i]))
        end
    end
    return points
end

"""
    make_dashboard(data; outdir="results", movie=true, fps=20, show_window=false)

Save `dashboard.png` and optionally `dashboard.mp4` in outdir.
Set show_window=true to open an interactive GLMakie window after rendering.
GLMakie requires a working OpenGL context, even when show_window=false.
"""
function make_dashboard(data::RunData;
                        outdir::AbstractString="results", movie::Bool=true,
                        fps::Int=20, show_window::Bool=false)
    isempty(data.frame_times) && error("No simulation snapshots to plot")
    mkpath(outdir)
    GLMakie.activate!()
    p = data.settings
    M = p.M0
    times = data.times
    # T = t - f(t) is asymptotically normalized in these gauge conventions.
    T = times .- data.shifts
    npts = length(T)
    nframes = length(data.frame_times)
    kref = argmin(abs.(p.observers .- 30*M))
    refwave = data.signals[kref]

    # An illustrative, single-mode fit (not a precision QNM measurement).
    fit = nothing
    if T[end] >= 95*M && maximum(abs, refwave) > 1e-24
        fit = fit_qnm(T, refwave; tlo=45*M, thi=min(110*M, T[end]), mass=data.masses_h[end])
    end
    lpis = [local_power_index(T, signal;
                peak_halfwidth=12*M, fit_halfwidth=24*M, tstart=100*M)
            for signal in data.signals]

    fig = Figure(size=(1550, 1040), fontsize=17)
    Label(fig[0, 1:2], "Einstein–Klein–Gordon | MMIEF | Method of Lines + RK4",
          fontsize=23, font=:bold)

    ax_pulse = Axis(fig[1,1], title="Scalar pulse and apparent horizon",
                    xlabel="J / M₀", ylabel="J φ")
    ax_mass = Axis(fig[1,2], title="Horizon and outer mass",
                   xlabel="t / M₀", ylabel="Mass / M₀")
    ax_wave = Axis(fig[2,1], title="Scalar field at fixed areal radii",
                   xlabel="T / M₀", ylabel="φ(T, Jₒᵦₛ)")
    ax_tail = Axis(fig[2,2], title="Power-law tail (log–log)",
                   xlabel="log₁₀(T / M₀)", ylabel="log₁₀ |φ|")
    ax_lpi = Axis(fig[3,1], title="Local decay exponent (peak envelope)",
                  xlabel="T / M₀", ylabel="−d ln(envelope) / d ln T")
    ax_constraints = Axis(fig[3,2], title="Constraint diagnostics",
                          xlabel="t / M₀", ylabel="log₁₀ |residual|")

    colors = [:dodgerblue, :darkorange, :seagreen, :purple, :firebrick]
    pulse_obs = Observable(Point2f[])
    horizon_obs = Observable(Float32[2.0])
    lines!(ax_pulse, pulse_obs; color=:royalblue, linewidth=2.5)
    vlines!(ax_pulse, horizon_obs; color=:firebrick, linestyle=:dash,
            linewidth=2, label="Apparent horizon")
    xlims!(ax_pulse, 2.0, min(p.rmax/M, 90.0))
    ylim = max(1e-9, maximum(maximum(abs, profile)
                            for profile in data.frame_profiles))
    ylims!(ax_pulse, -1.2*ylim, 1.2*ylim)

    horizon_mass = Observable(Point2f[])
    outer_mass = Observable(Point2f[])
    lines!(ax_mass, horizon_mass; color=:firebrick, linewidth=2.3,
           label="Mₕ = Jₕ/2")
    lines!(ax_mass, outer_mass; color=:slategray, linewidth=1.8,
           linestyle=:dash, label="m at outer boundary")
    axislegend(ax_mass, position=:rb)

    obs_wave = [Observable(Point2f[]) for _ in p.observers]
    obs_tail = [Observable(Point2f[]) for _ in p.observers]
    obs_lpi = [Observable(Point2f[]) for _ in p.observers]
    for k in eachindex(p.observers)
        col = colors[mod1(k, length(colors))]
        label = @sprintf("J = %.0f M₀", p.observers[k]/M)
        lines!(ax_wave, obs_wave[k]; color=col, linewidth=2,
               label=label)
        lines!(ax_tail, obs_tail[k]; color=col, linewidth=2,
               label=label)
        lines!(ax_lpi, obs_lpi[k]; color=col, linewidth=2,
               label=label)
    end
    axislegend(ax_wave, position=:rt)
    axislegend(ax_tail, position=:lb)
    hlines!(ax_lpi, [3.0]; color=:black, linestyle=:dash, linewidth=2)
    xlims!(ax_lpi, 0, max(T[end]/M,1.0))
    ylims!(ax_lpi, -1, 8)

    if fit !== nothing
        overlay = Point2f[]
        for t in T
            if fit.tlo <= t <= fit.thi
                tau = t-fit.tlo
                envelope = exp(-fit.gammaM / fit.mass * tau)
                value = envelope * (fit.C*cos(fit.omegaM / fit.mass * tau) +
                                    fit.S*sin(fit.omegaM / fit.mass * tau))
                push!(overlay, p2(t/M, value))
            end
        end
        lines!(ax_wave, overlay; color=:black, linestyle=:dash,
               linewidth=2, label="Exploratory QNM fit")
    end

    if T[end] > 120*M
        tref = 150*M <= T[end] ? 150*M : 0.8*T[end]
        ref_idx = clamp(searchsortedfirst(T, tref), 1, npts)
        ilo = max(1,ref_idx-10)
        ihi = min(npts,ref_idx+10)
        reference_amplitude = maximum(abs, view(refwave,ilo:ihi))
        if reference_amplitude > 1e-25
            tail_ref = Point2f[]
            for t in range(100*M, T[end], length=100)
                amp = reference_amplitude * (tref/t)^3
                push!(tail_ref, p2(log10(t/M),log10(amp)))
            end
            lines!(ax_tail, tail_ref; color=:black, linestyle=:dash,
                   linewidth=2, label="Reference slope −3")
        end
    end

    ham_obs = Observable(Point2f[])
    def_obs = Observable(Point2f[])
    horizon_drift_obs = Observable(Point2f[])
    lines!(ax_constraints, ham_obs; color=:dodgerblue, label="Hamiltonian L₂")
    lines!(ax_constraints, def_obs; color=:seagreen, label="Φ − ∂ᵣφ L₂")
    lines!(ax_constraints, horizon_drift_obs; color=:firebrick,
           label="Preprojection horizon drift")
    axislegend(ax_constraints, position=:rb)

    r_frame = data.r[1:p.profile_stride:end]
    function update_frame!(k::Int)
        k = clamp(k, 1, nframes)
        t = data.frame_times[k]
        f = data.frame_shift[k]
        profile = data.frame_profiles[k]
        pulse_obs[] = Point2f[p2((r_frame[i]+f)/M, profile[i])
                               for i in eachindex(profile)]
        horizon_obs[] = Float32[(2*p.M0+f)/M]
        idx = clamp(searchsortedlast(times, t+1e-10), 1, npts)
        horizon_mass[] = trace_points(times, data.masses_h, idx, M)
        outer_mass[] = trace_points(times, data.masses_out, idx, M)
        for j in eachindex(p.observers)
            obs_wave[j][] = trace_points(T, data.signals[j], idx, M)
            obs_tail[j][] = log_points(T, data.signals[j], idx, M)
            obs_lpi[j][] = lpi_points(T, lpis[j], idx, M)
        end
        ham_obs[] = residual_points(times, data.ham_l2, idx, M)
        def_obs[] = residual_points(times, data.scalar_constraint_l2, idx, M)
        horizon_drift_obs[] = residual_points(times,
                                   data.horizon_preprojection, idx, M)
        ax_pulse.title = @sprintf("Scalar profile at t/M₀ = %.1f", t/M)
    end

    update_frame!(nframes)
    png = joinpath(outdir, "dashboard.png")
    save(png, fig)
    println("Saved GLMakie dashboard: ", png)
    if movie
        mp4 = joinpath(outdir, "dashboard.mp4")
        record(fig, mp4, 1:nframes; framerate=fps) do i
            update_frame!(i)
        end
        println("Saved GLMakie animation: ", mp4)
        update_frame!(nframes)
    end
    if show_window
        display(fig)
    end
    return (figure=fig, qnm_fit=fit, local_power_indices=lpis)
end

end # module Dashboard
