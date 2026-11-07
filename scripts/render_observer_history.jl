using GLMakie
using Printf

# -------------------------
# Lightweight CSV utilities
# -------------------------
function _split_csv_line(line::String)
    return String.(strip.(split(chomp(line), ',')))
end

function read_simple_csv(path::AbstractString)
    isfile(path) || error("CSV file not found: $path")
    lines = readlines(path)
    isempty(lines) && error("Empty CSV file: $path")
    header = _split_csv_line(lines[1])
    cols = Dict{String,Vector{String}}(h => String[] for h in header)
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

function find_col(cols::Dict{String,Vector{String}}, aliases::Vector{String})
    lower_map = Dict(lowercase(k) => k for k in keys(cols))
    for a in aliases
        key = get(lower_map, lowercase(a), nothing)
        key === nothing || return key
    end
    return nothing
end

tofloat(s::AbstractString) = try parse(Float64, s) catch; NaN end

function get_num_col(cols::Dict{String,Vector{String}}, aliases::Vector{String})
    key = find_col(cols, aliases)
    key === nothing && return nothing, Float64[]
    return key, [tofloat(v) for v in cols[key]]
end

# -------------------------
# Observer data loading
# -------------------------
function load_observers(run_dir::AbstractString)
    path = joinpath(run_dir, "observers.csv")
    cols = read_simple_csv(path)

    _, T   = get_num_col(cols, ["T", "t_asymptotic", "time_asymptotic", "time", "t"])
    _, J   = get_num_col(cols, ["J_obs", "Jobs", "J", "J_over_M0"])
    _, phi = get_num_col(cols, ["phi", "varphi", "field", "scalar"])

    if isempty(T) || isempty(J) || isempty(phi)
        error("Could not identify T/J_obs/phi columns in observers.csv")
    end

    data = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
    for target in [10.0, 30.0, 50.0]
        idx = findall(k -> isfinite(J[k]) && abs(J[k] - target) < 1e-8, eachindex(J))
        if !isempty(idx)
            data[target] = (T[idx], phi[idx])
        end
    end
    return data
end

# -------------------------
# Plotting
# -------------------------
function make_observer_history(run_dir::AbstractString; show_window::Bool=false)
    data = load_observers(run_dir)
    isempty(data) && error("No observer data found for J = 10, 30, 50 in $(joinpath(run_dir, "observers.csv"))")

    fig = Figure(size=(1600, 1200))
    ax1 = Axis(fig[1,1], xlabel="T / M0", ylabel="ϕ", title="Entire waveform history at fixed areal radii")
    ax2 = Axis(fig[1,2], xlabel="T / M0", ylabel="log10 |ϕ|", title="Amplitude history: log10 |ϕ(T)|")
    ax3 = Axis(fig[2,1], xlabel="T / M0", ylabel="|ϕ|", title="Absolute amplitude history")
    ax4 = Axis(fig[2,2], xlabel="T / M0", ylabel="cumulative max |ϕ|", title="Running envelope max{|ϕ|}")

    for J in sort(collect(keys(data)))
        t, phi = data[J]
        label = @sprintf("J = %.0f M0", J)
        lines!(ax1, t, phi, label=label)
        lines!(ax2, t, log10.(abs.(phi) .+ 1e-300), label=label)
        lines!(ax3, t, abs.(phi), label=label)
        env = accumulate(max, abs.(phi))
        lines!(ax4, t, env, label=label)
    end

    axislegend(ax1, position=:rb)
    axislegend(ax2, position=:rb)
    axislegend(ax3, position=:rb)
    axislegend(ax4, position=:rb)

    Label(fig[3,1:2],
          "This figure shows the full stored time history of the scalar amplitude at the fixed-areal-radius observers J = 10M₀, 30M₀, and 50M₀.\n" *
          "Top left: raw waveform. Top right: logarithmic amplitude. Bottom left: absolute amplitude. Bottom right: running envelope maximum.",
          tellwidth=false, justification=:left)

    outpath = joinpath(abspath(run_dir), "observer_history_J10_J30_J50.png")
    save(outpath, fig)
    println("Saved observer-history figure: $outpath")

    if show_window
        display(fig)
    end
    return outpath
end

function observer_history_main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_observer_history.jl <run_directory> [--show]")
    run_dir = abspath(args[1])
    allowed = Set(["--show"])
    all(x -> x in allowed, args[2:end]) || error("Unknown render option")
    return make_observer_history(run_dir; show_window=("--show" in args))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    observer_history_main(ARGS)
end
