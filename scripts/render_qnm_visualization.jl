using Printf
using GLMakie

# -------------------------
# Lightweight CSV utilities
# -------------------------
function _split_csv_line(line::String)
    # Simple CSV splitter; sufficient for our generated numeric/text files
    return String.(strip.(split(chomp(line), ',')))
end

function read_simple_csv(path::AbstractString)
    isfile(path) || error("CSV file not found: $path")
    lines = readlines(path)
    isempty(lines) && error("Empty CSV file: $path")
    header = _split_csv_line(lines[1])
    cols = Dict(h => String[] for h in header)
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = _split_csv_line(line)
        if length(vals) < length(header)
            append!(vals, fill("", length(header) - length(vals)))
        elseif length(vals) > length(header)
            vals = vals[1:length(header)]
        end
        for (h, v) in zip(header, vals)
            push!(cols[h], v)
        end
    end
    return cols
end

function find_col(cols::AbstractDict{<:AbstractString,<:AbstractVector{<:AbstractString}}, aliases::Vector{String})
    keys_lower = Dict(lowercase(k) => k for k in keys(cols))
    for alias in aliases
        lk = lowercase(alias)
        if haskey(keys_lower, lk)
            return keys_lower[lk]
        end
    end
    return nothing
end

tofloat(s::AbstractString) = try parse(Float64, s) catch; NaN end

function get_num_col(cols::AbstractDict{<:AbstractString,<:AbstractVector{<:AbstractString}}, aliases::Vector{String})
    key = find_col(cols, aliases)
    key === nothing && return nothing, Float64[]
    return key, [tofloat(v) for v in cols[key]]
end

# -------------------------
# Data extraction
# -------------------------
function run_mass_scales(run_dir::AbstractString)
    config_path = joinpath(run_dir, "config.toml")
    isfile(config_path) || error("Missing config.toml: $config_path")
    # Read this single numeric key without adding dependencies to the viz environment.
    # save_run() emits an ordinary numeric TOML assignment: M0 = 1.0.
    M0 = NaN
    for line in eachline(config_path)
        pair = split(split(line, '#'; limit=2)[1], '='; limit=2)
        if length(pair) == 2 && strip(pair[1]) == "M0"
            M0 = parse(Float64, strip(pair[2]))
            break
        end
    end
    isfinite(M0) && M0 > 0 || error("Invalid initial mass M0")

    Mf = M0
    diag_path = joinpath(run_dir, "diagnostics.csv")
    if isfile(diag_path)
        diagnostics = read_simple_csv(diag_path)
        _, mh = get_num_col(diagnostics, ["M_h", "M_horizon"])
        valid = filter(isfinite, mh)
        if !isempty(valid)
            Mf = valid[end]
        end
    else
        @warn "diagnostics.csv not found; using Mf=M0 for theoretical comparison"
    end
    isfinite(Mf) && Mf > 0 || error("Invalid final horizon mass Mf")
    return M0, Mf
end

function load_observer_data(run_dir::AbstractString, M0::Float64)
    path = joinpath(run_dir, "observers.csv")
    cols = read_simple_csv(path)

    _, T = get_num_col(cols, ["T_asymptotic", "time_asymptotic", "T", "time", "t"])
    _, J = get_num_col(cols, ["J_obs", "Jobs", "J", "J_over_M0"])
    _, phi = get_num_col(cols, ["phi", "varphi", "field", "scalar"])

    if isempty(T) || isempty(J) || isempty(phi)
        error("Could not identify columns T/J_obs/phi in observers.csv")
    end

    # Store both the time coordinate and observer radius in units of M0.
    # This matches the convention of qnm_comparison.csv (J_over_M0 and Tlo_over_M0).
    Jnorm = J ./ M0
    Tnorm = T ./ M0
    groups = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
    for j in sort(unique(filter(isfinite, Jnorm)))
        idx = findall(x -> isfinite(x) && isapprox(x, j; atol=1e-8, rtol=0), Jnorm)
        good = [i for i in idx if isfinite(Tnorm[i]) && isfinite(phi[i])]
        order = sortperm(Tnorm[good])
        groups[j] = (Tnorm[good][order], phi[good][order])
    end
    return groups
end

function load_qnm_comparison(run_dir::AbstractString, M0::Float64)
    path = joinpath(run_dir, "qnm_comparison.csv")
    isfile(path) || return nothing
    cols = read_simple_csv(path)

    # The QNM checker produces these normalized fields exactly:
    # J_over_M0,Tlo_over_M0,Thi_over_M0,Mf_omegaR,Mf_gamma,
    # omegaR_error_pct,gamma_error_pct,relative_fit_residual.
    jkey, J = get_num_col(cols, ["J_over_M0", "J_obs", "Jobs", "J"])
    _, wR = get_num_col(cols, ["Mf_omegaR", "Mf_omegaR_num", "MfomegaR", "omegaR_Mf", "omegaR"])
    _, gamma = get_num_col(cols, ["Mf_gamma", "Mf_Gamma", "Mf_Gamma_num", "MfGamma", "Gamma_Mf", "gamma"])
    _, residual = get_num_col(cols, ["relative_fit_residual", "rel_residual", "relative_residual", "residual", "fit_residual"])
    t1key, T1 = get_num_col(cols, ["Tlo_over_M0", "T_start", "t_start", "fit_start", "T1", "t1"])
    t2key, T2 = get_num_col(cols, ["Thi_over_M0", "T_end", "t_end", "fit_end", "T2", "t2"])

    isempty(J) && error("qnm_comparison.csv is missing its observer-radius column. Headers: $(join(keys(cols), ", "))")
    isempty(wR) && error("qnm_comparison.csv is missing Mf_omegaR. Headers: $(join(keys(cols), ", "))")
    isempty(gamma) && error("qnm_comparison.csv is missing Mf_gamma. Headers: $(join(keys(cols), ", "))")

    if jkey != "J_over_M0"
        J ./= M0
    end
    if t1key !== nothing && t1key != "Tlo_over_M0"
        T1 ./= M0
    end
    if t2key !== nothing && t2key != "Thi_over_M0"
        T2 ./= M0
    end
    n = maximum((length(J), length(wR), length(gamma), length(residual), length(T1), length(T2)))
    rows = NamedTuple[]
    for i in 1:n
        push!(rows, (
            J_obs = i <= length(J) ? J[i] : NaN,
            Mf_omegaR = i <= length(wR) ? wR[i] : NaN,
            Mf_Gamma = i <= length(gamma) ? gamma[i] : NaN,
            residual = i <= length(residual) ? residual[i] : NaN,
            T_start = i <= length(T1) ? T1[i] : NaN,
            T_end = i <= length(T2) ? T2[i] : NaN,
        ))
    end
    return rows
end

function best_fit_rows(rows)
    rows === nothing && return Dict{Float64,NamedTuple}()
    byJ = Dict{Float64,Vector{NamedTuple}}()
    for r in rows
        isfinite(r.J_obs) || continue
        get!(byJ, r.J_obs, NamedTuple[])
        push!(byJ[r.J_obs], r)
    end
    out = Dict{Float64,NamedTuple}()
    for (J, rs) in byJ
        if any(isfinite.(getfield.(rs, :residual)))
            idx = argmin([isfinite(r.residual) ? r.residual : Inf for r in rs])
            out[J] = rs[idx]
        else
            out[J] = rs[1]
        end
    end
    return out
end

# -------------------------
# Ringdown model helper
# -------------------------
function build_ringdown_model(t::Vector{Float64}, y::Vector{Float64}, omega::Float64, gamma::Float64)
    isempty(t) && return Float64[], NaN, NaN, NaN
    dt = t .- t[1]
    basis_c = exp.(-gamma .* dt) .* cos.(omega .* dt)
    basis_s = exp.(-gamma .* dt) .* sin.(omega .* dt)
    cc = sum(abs2, basis_c)
    ss = sum(abs2, basis_s)
    cs = sum(basis_c .* basis_s)
    yc = sum(y .* basis_c)
    ys = sum(y .* basis_s)
    determinant = cc*ss - cs*cs
    if !isfinite(determinant) || determinant <= 1e-14*max(cc*ss, eps(Float64))
        return fill(NaN, length(t)), NaN, NaN, Inf
    end
    coeff_c = (yc*ss - ys*cs) / determinant
    coeff_s = (ys*cc - yc*cs) / determinant
    yfit = coeff_c .* basis_c .+ coeff_s .* basis_s
    relative_residual = sqrt(sum(abs2, yfit .- y)) / max(sqrt(sum(abs2, y)), eps(Float64))
    return yfit, coeff_c, coeff_s, relative_residual
end

function nearest_key(dict::Dict{Float64,<:Any}, target::Float64)
    ks = collect(keys(dict))
    isempty(ks) && return nothing
    return ks[argmin(abs.(ks .- target))]
end

# -------------------------
# Plotting
# -------------------------
function make_qnm_overview(run_dir::AbstractString; show_window::Bool=false)
    theory_wR = 0.110455
    theory_g  = 0.104896

    M0, Mf = run_mass_scales(run_dir)
    obs = load_observer_data(run_dir, M0)
    rawfits = load_qnm_comparison(run_dir, M0)
    fits = best_fit_rows(rawfits)
    rawfits === nothing && @warn "qnm_comparison.csv not found; showing waveforms only. Run qnm_check first for fit overlays."

    outdir = abspath(run_dir)
    isempty(obs) && error("No valid observer waveforms in observers.csv")

    # Overview figure
    fig = Figure(size=(1700, 1100))
    ax1 = Axis(fig[1,1], xlabel="T / M0", ylabel="ϕ", title="Observer waveforms")
    ax2 = Axis(fig[1,2], xlabel="T / M0", ylabel="log10 |ϕ|", title="Amplitude envelope")
    ax3 = Axis(fig[2,1], xlabel="T / M0", ylabel="ϕ", title="Ringdown windows and fitted models")
    ax4 = Axis(fig[2,2], xlabel="observer index", ylabel="dimensionless value", title="QNM frequencies vs theory")

    Jvals = sort(collect(keys(obs)))
    xticks_vals = 1:length(Jvals)
    xticks_labels = [@sprintf("J=%.0f", j) for j in Jvals]
    ax4.xticks = (xticks_vals, xticks_labels)

    fit_indices = Int[]
    ωR_vals = Float64[]
    Γ_vals  = Float64[]
    has_fit_lines = false

    for (k, J) in enumerate(Jvals)
        t, y = obs[J]
        lines!(ax1, t, y, label=@sprintf("J=%.0f", J))
        lines!(ax2, t, log10.(abs.(y) .+ 1e-300), label=@sprintf("J=%.0f", J))

        if haskey(fits, J)
            fit = fits[J]
            if isfinite(fit.Mf_omegaR) && isfinite(fit.Mf_Gamma)
                push!(fit_indices, k)
                push!(ωR_vals, fit.Mf_omegaR)
                push!(Γ_vals, fit.Mf_Gamma)
            end

            if isfinite(fit.T_start) && isfinite(fit.T_end)
                idx = findall(tt -> tt >= fit.T_start && tt <= fit.T_end, t)
                if length(idx) >= 5 && isfinite(fit.Mf_omegaR) && isfinite(fit.Mf_Gamma)
                    tw = t[idx]
                    yw = y[idx]
                    # t is in units M0; input fitted rates are in units Mf.
                    rate_scale = M0 / Mf
                    yfit, _, _, relres = build_ringdown_model(tw, yw,
                        fit.Mf_omegaR * rate_scale, fit.Mf_Gamma * rate_scale)
                    lines!(ax3, tw, yw, label=@sprintf("data J=%.0f", J))
                    lines!(ax3, tw, yfit, linestyle=:dash, label=@sprintf("fit J=%.0f (res=%.3g)", J, relres))
                    has_fit_lines = true
                end
            end
        end
    end

    axislegend(ax1, position=:rb)
    axislegend(ax2, position=:rb)
    has_fit_lines && axislegend(ax3, position=:rb)

    # Frequency comparison (plot only observers with a usable numerical fit)
    if !isempty(fit_indices)
        scatter!(ax4, Float64.(fit_indices) .- 0.08, ωR_vals, markersize=14, label="Mf·ωR (fit)")
        scatter!(ax4, Float64.(fit_indices) .+ 0.08, Γ_vals, markersize=14, label="Mf·Γ (fit)")
    end
    lines!(ax4, [0.5, length(Jvals)+0.5], [theory_wR, theory_wR], linestyle=:dash, label="theory Mf·ωR")
    lines!(ax4, [0.5, length(Jvals)+0.5], [theory_g, theory_g], linestyle=:dot, label="theory Mf·Γ")
    axislegend(ax4, position=:rt)

    # Add text summary
    summary_lines = String[]
    push!(summary_lines, @sprintf("Theory: Mf·ωR = %.6f, Mf·Γ = %.6f", theory_wR, theory_g))
    for J in Jvals
        if haskey(fits, J)
            fit = fits[J]
            δω = 100*(fit.Mf_omegaR/theory_wR - 1)
            δγ = 100*(fit.Mf_Gamma/theory_g - 1)
            push!(summary_lines,
                @sprintf("J=%.0f: fit Mf·ωR=%.5f (%.1f%%), Mf·Γ=%.5f (%.1f%%)",
                        J, fit.Mf_omegaR, δω, fit.Mf_Gamma, δγ))
        end
    end
    Label(fig[3,1:2], join(summary_lines, "\n"), tellwidth=false, justification=:left)

    overview_path = joinpath(outdir, "qnm_overview.png")
    save(overview_path, fig)

    # Per-observer fit figures
    fit_paths = String[]
    for J in Jvals
        haskey(fits, J) || continue
        t, y = obs[J]
        fit = fits[J]
        if !(isfinite(fit.T_start) && isfinite(fit.T_end) && isfinite(fit.Mf_omegaR) && isfinite(fit.Mf_Gamma))
            continue
        end
        idx = findall(tt -> tt >= fit.T_start && tt <= fit.T_end, t)
        length(idx) >= 5 || continue
        tw = t[idx]
        yw = y[idx]
        rate_scale = M0 / Mf
        yfit, c1, c2, relres = build_ringdown_model(tw, yw,
            fit.Mf_omegaR * rate_scale, fit.Mf_Gamma * rate_scale)

        figJ = Figure(size=(1400, 800))
        axJ1 = Axis(figJ[1,1], xlabel="T / M0", ylabel="ϕ", title=@sprintf("Ringdown fit at J=%.0f", J))
        axJ2 = Axis(figJ[1,2], xlabel="T / M0", ylabel="residual", title="Fit residual")
        lines!(axJ1, tw, yw, label="data")
        lines!(axJ1, tw, yfit, linestyle=:dash, label="damped-sinusoid fit")
        axislegend(axJ1, position=:rb)
        lines!(axJ2, tw, yw .- yfit)

        δω = 100*(fit.Mf_omegaR/theory_wR - 1)
        δγ = 100*(fit.Mf_Gamma/theory_g - 1)
        text = @sprintf("Theory: Mf·ωR = %.6f, Mf·Γ = %.6f\nFit:    Mf·ωR = %.6f, Mf·Γ = %.6f\nErrors: ΔωR = %.1f%%, ΔΓ = %.1f%%\nWindow: [%.2f, %.2f]\nRelative residual: %.3g\nModel coefficients: C=%.3e, S=%.3e",
                        theory_wR, theory_g,
                        fit.Mf_omegaR, fit.Mf_Gamma,
                        δω, δγ,
                        fit.T_start, fit.T_end,
                        relres, c1, c2)
        Label(figJ[2,1:2], text, tellwidth=false, justification=:left)
        pathJ = joinpath(outdir, @sprintf("qnm_fit_J%.0f.png", J))
        save(pathJ, figJ)
        push!(fit_paths, pathJ)
    end

    if show_window
        display(fig)
    end

    println("Saved QNM overview: $overview_path")
    for p in fit_paths
        println("Saved QNM fit figure: $p")
    end

    return (overview_path=overview_path, fit_paths=fit_paths)
end

function qnm_viz_main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_qnm_visualization.jl <run_directory> [--show]")
    run_dir = abspath(args[1])
    allowed = Set(["--show"])
    all(x -> x in allowed, args[2:end]) || error("Unknown render option")
    return make_qnm_overview(run_dir; show_window=("--show" in args))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    qnm_viz_main(ARGS)
end
