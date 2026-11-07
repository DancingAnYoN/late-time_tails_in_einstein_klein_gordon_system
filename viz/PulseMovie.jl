module PulseMovie
using GLMakie
using Printf
using MMIEF: RunData

export make_pulse_movie

"""
    make_pulse_movie(data::RunData; outdir="results/pulse_movie", fps=24,
                     every=1, radial_max=75.0, show_window=false)

Render a focused MP4 of the *computed* MMIEF scalar field using saved solver
snapshots. The vertical marker is the apparent horizon at physical areal
radius J_h=2M0+f(t). The left plot uses physical areal radius J; the bottom
right spacetime diagram uses the fixed numerical coordinate r and therefore
is exactly aligned across snapshot times.

No artificial propagation or absorption is injected into the frames.
"""
function make_pulse_movie(data::RunData;
                          outdir::AbstractString=joinpath("results", "pulse_movie"),
                          fps::Int=24, every::Int=1,
                          radial_max::Real=75.0,
                          show_window::Bool=false)
    GLMakie.activate!()
    mkpath(outdir)
    p=data.settings
    M=p.M0
    nframes=length(data.frame_times)
    nframes >= 2 || error("At least two saved snapshots are needed to make a movie")
    every >= 1 || error("every must be >= 1")
    fps >= 1 || error("fps must be >= 1")
    frames=collect(1:every:nframes)
    frames[end] != nframes && push!(frames, nframes)
    rvec=data.r[1:p.profile_stride:end]
    nr=length(rvec)
    all(length(profile)==nr for profile in data.frame_profiles) ||
        error("Inconsistent profile lengths in RunData")

    # The normalized plotted field is psi/(M0*epsilon), with psi=J*phi.
    # For a vacuum test the normalization is 1, producing an identically zero line.
    scale=p.pulse_amplitude > 0 ? M*p.pulse_amplitude : 1.0
    normalized=[Float64.(profile)./scale for profile in data.frame_profiles]
    amplitude=max(maximum(maximum(abs, profile) for profile in normalized), 1e-9)
    ymax=1.18*amplitude
    jmax=min(Float64(radial_max), (data.r[end]+minimum(data.frame_shift))/M)
    jmax>3 || error("radial_max must extend beyond the horizon")

    fig=Figure(size=(1480,900),fontsize=18,backgroundcolor=RGBf(0.98,0.985,0.995))
    Label(fig[0, 1:2], "Scalar pulse falling into a Schwarzschild black hole",
          fontsize=29, font=:bold, color=RGBf(0.09,0.16,0.27))
    time_label=Observable("t / M₀ = 0.0")
    Label(fig[1, 1:2], time_label, fontsize=20,
          color=RGBf(0.20,0.36,0.55))

    ax_pulse=Axis(fig[2:3,1],
        title="Radial scalar pulse | physical areal radius",
        xlabel="J / M₀", ylabel="J φ / (M₀ ε)",
        backgroundcolor=:white)
    ax_mass=Axis(fig[2,2],
        title="Apparent-horizon growth",
        xlabel="t / M₀", ylabel="10⁶ × (Mₕ/M₀ − 1)",
        backgroundcolor=:white)
    ax_map=Axis(fig[3,2],
        title="Spacetime history | fixed numerical grid",
        xlabel="r / M₀", ylabel="t / M₀",
        backgroundcolor=:white)
    colsize!(fig.layout, 1, Relative(0.60))

    # Dynamic horizon shading and a dashed horizon line in the main plot.
    shade=Observable(Point2f[])
    poly!(ax_pulse,shade;color=(:black,0.095),strokewidth=0)
    hJ=Observable(Float32[2.0])
    vlines!(ax_pulse,hJ;color=:firebrick,linewidth=3,linestyle=:dash,
            label="Apparent horizon")
    hlines!(ax_pulse,[0.0];color=(:black,0.25),linewidth=1)
    pulse=Observable(Point2f[])
    lines!(ax_pulse,pulse;color=:royalblue,linewidth=3.6,
           label="Scalar wave Jφ")
    xlims!(ax_pulse,0,jmax)
    ylims!(ax_pulse,-ymax,ymax)
    axislegend(ax_pulse,position=:rt,labelsize=15)

    # Horizon growth is small for epsilon=1e-4 but not zero in the coupled run.
    total_t=max(data.frame_times[end]/M,0.5)
    mh_initial=data.masses_h[1]
    mchanges=1e6 .* ((data.masses_h .- mh_initial)./M)
    mmax=max(maximum(mchanges),1e-5)
    mass_line=Observable(Point2f[])
    lines!(ax_mass,mass_line;color=:firebrick,linewidth=3)
    xlims!(ax_mass,0,total_t)
    ylims!(ax_mass,-0.05*mmax,1.18*mmax)

    # r is fixed across all frames. Thus this heatmap does not falsely imply
    # a stationary *physical* radius in the presence of the moving f(t).
    shown_times=Float64.(data.frame_times)./M
    rscaled=Float64.(rvec)./M
    map_data=fill(Float32(NaN),nr,nframes)
    map_obs=Observable(copy(map_data))
    heatmap!(ax_map,rscaled,shown_times,map_obs;
             colormap=:balance,colorrange=(-amplitude,amplitude),nan_color=:white)
    vlines!(ax_map,[2.0];color=:black,linewidth=2.4,linestyle=:dash)
    now_t=Observable(Float32[0.0])
    hlines!(ax_map,now_t;color=:firebrick,linewidth=2)
    xlims!(ax_map,2,min(jmax,rscaled[end]))
    ylims!(ax_map,0,total_t)

    Label(fig[4,1:2],
        "MMIEF Einstein–Klein–Gordon | Method of Lines + RK4 | " *
        "Constrained evolution | Numerics require convergence tests",
        fontsize=14,color=RGBf(0.35,0.42,0.49))

    last_filled=Ref(0)
    function update!(k::Int)
        t=data.frame_times[k]
        f=data.frame_shift[k]
        prof=normalized[k]
        Jh=(2M+f)/M
        pulse[]=Point2f[Point2f(Float32((rvec[i]+f)/M),Float32(prof[i]))
                         for i in eachindex(rvec)]
        hJ[]=Float32[Jh]
        shade[]=Point2f[
            Point2f(0,-ymax), Point2f(Float32(Jh),-ymax),
            Point2f(Float32(Jh),ymax), Point2f(0,ymax)]
        idx=clamp(searchsortedlast(data.times,t+1e-10),1,length(data.times))
        mass_line[]=Point2f[Point2f(Float32(data.times[j]/M),
            Float32(1e6*(data.masses_h[j]-mh_initial)/M)) for j in 1:idx]
        # Write only the newly revealed columns; keep future columns blank.
        for j in (last_filled[]+1):k
            map_data[:,j].=Float32.(normalized[j])
        end
        last_filled[]=k
        map_obs[]=copy(map_data)
        now_t[]=Float32[t/M]
        time_label[]=@sprintf("t / M₀ = %5.1f    |    Jₕ / M₀ = %.6f",t/M,Jh)
    end

    # Save a clean final frame as well as the MP4.
    path=joinpath(outdir,"scalar_pulse_absorption.mp4")
    record(fig,path,frames;framerate=fps) do k
        update!(k)
    end
    save(joinpath(outdir,"scalar_pulse_final.png"),fig)
    @printf("Saved scalar-pulse GLMakie movie: %s (%d frames, %d fps)\n",
            path,length(frames),fps)
    if show_window
        display(fig)
    end
    return path
end
end # module
