# Standalone interpretation of the MMIEF scalar l=0 QNM comparison table.
# Requires only Printf and TOML, both declared by the MMIEF core project.
# Run from MMIEF_Modular root:
#   include("scripts/interpret_qnm_results.jl")
#   QNMInterpretation.interpret_qnm("results/standard")
# Or from PowerShell:
#   julia --project=. scripts/interpret_qnm_results.jl results/standard

module QNMInterpretation
using Printf
using TOML

export interpret_qnm

# Fundamental scalar, spherical (l=0), Schwarzschild QNM for exp(-i*omega*t).
# Gamma > 0 is defined by omega = omega_R - i*Gamma.
const THEORY_OMEGA_R_M = 0.110455
const THEORY_GAMMA_M   = 0.104896

function _read_numeric_csv(path::AbstractString, required::Vector{String})
    isfile(path) || error("Missing file: $path. Run qnm_check(\"results/standard\") first.")
    rows = readlines(path)
    isempty(rows) && error("CSV is empty: $path")
    names = String.(strip.(split(strip(rows[1]), ',')))
    all(in(names), required) || error("Unexpected headers in $(basename(path)). Missing: $(setdiff(required, names))")
    records = Vector{Dict{String,Float64}}()
    for (j, line) in enumerate(rows[2:end])
        line_number = j + 1
        isempty(strip(line)) && continue
        cells = strip.(split(line, ','))
        length(cells) == length(names) || error("Column count mismatch in $(basename(path)):$line_number")
        vals = tryparse.(Float64, cells)
        any(isnothing, vals) && error("Non-numeric value in $(basename(path)):$line_number")
        record = Dict(names[i] => vals[i]::Float64 for i in eachindex(names))
        push!(records, record)
    end
    return records
end

function _median(values::Vector{Float64})
    isempty(values) && return NaN
    ordered = sort(values)
    n = length(ordered)
    return isodd(n) ? ordered[(n+1)÷2] : (ordered[n÷2]+ordered[n÷2+1])/2
end

_pct_error(x, ref) = 100.0 * (x / ref - 1.0)
_spread_pct(x, ref) = length(x) < 2 ? NaN : 100.0 * (maximum(x)-minimum(x))/ref

function _linear_sample(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    length(xs) == length(ys) || error("Diagnostic time/mass length mismatch")
    isempty(xs) && return nothing
    (x < xs[1] || x > xs[end]) && return nothing
    i = searchsortedlast(xs, x)
    i == length(xs) && return ys[end]
    xs[i+1] > xs[i] || error("Repeated diagnostic times")
    w = (x-xs[i])/(xs[i+1]-xs[i])
    return (1-w)*ys[i] + w*ys[i+1]
end

function _mass_history(run_dir::AbstractString)
    diagnostic_path = joinpath(run_dir, "diagnostics.csv")
    config_path = joinpath(run_dir, "config.toml")
    if !isfile(diagnostic_path) || !isfile(config_path)
        return nothing
    end
    config = TOML.parsefile(config_path)
    haskey(config,"M0") || return nothing
    M0 = Float64(config["M0"])
    M0 > 0 || error("M0 must be positive")
    rows = _read_numeric_csv(diagnostic_path, ["T_asymptotic","M_h"])
    isempty(rows) && return nothing
    ts = [r["T_asymptotic"] for r in rows]
    mh = [r["M_h"] for r in rows]
    issorted(ts) && all(diff(ts).>0) || error("Diagnostic asymptotic times are not strictly increasing")
    all(isfinite, mh) && all(>(0),mh) || error("Invalid horizon mass in diagnostics.csv")
    return (M0=M0, Mf=mh[end], ts=ts, mass=mh)
end

_fmt(x;digits=4) = isfinite(x) ? @sprintf("%.*f", digits, x) : "n/a"

"""
    interpret_qnm(run_directory="results/standard"; freq_tolerance_pct=5.0,
                  spread_tolerance_pct=5.0, max_fit_residual=0.20,
                  max_horizon_mass_drift_pct=0.5, min_windows=3)

Interpret *exploratory* MMIEF scalar l=0 fits saved in qnm_comparison.csv.
Writes qnm_interpretation.md and qnm_interpretation_summary.csv.
A positive verdict is only preliminary compatibility, never QNM confirmation.
"""
function interpret_qnm(run_directory::AbstractString="results/standard";
                       freq_tolerance_pct::Float64=5.0,
                       spread_tolerance_pct::Float64=5.0,
                       max_fit_residual::Float64=0.20,
                       max_horizon_mass_drift_pct::Float64=0.5,
                       min_windows::Int=3)
    all(x->x>0, (freq_tolerance_pct,spread_tolerance_pct,max_fit_residual,
                 max_horizon_mass_drift_pct)) || error("All tolerances must be positive")
    min_windows >= 2 || error("min_windows must be >= 2")
    run_dir = abspath(run_directory)
    required = ["J_over_M0","Tlo_over_M0","Thi_over_M0","Mf_omegaR",
                "Mf_gamma","relative_fit_residual"]
    rows = _read_numeric_csv(joinpath(run_dir,"qnm_comparison.csv"), required)
    rows = [row for row in rows if all(isfinite, values(row)) &&
             row["Mf_omegaR"]>0 && row["Mf_gamma"]>0 &&
             row["relative_fit_residual"]>=0 &&
             row["Thi_over_M0"]>row["Tlo_over_M0"]]
    history = _mass_history(run_dir)

    locations = sort(unique(row["J_over_M0"] for row in rows))
    report_path = joinpath(run_dir, "qnm_interpretation.md")
    summary_path = joinpath(run_dir, "qnm_interpretation_summary.csv")
    summary_rows = NamedTuple[]
    reasons_by_observer = Vector{String}()

    for J_over_M0 in locations
        selected = [r for r in rows if r["J_over_M0"] == J_over_M0]
        w = [r["Mf_omegaR"] for r in selected]
        g = [r["Mf_gamma"] for r in selected]
        residuals = [r["relative_fit_residual"] for r in selected]
        wmed, gmed = _median(w), _median(g)
        eR, eI = _pct_error(wmed, THEORY_OMEGA_R_M), _pct_error(gmed, THEORY_GAMMA_M)
        sR, sI = _spread_pct(w, THEORY_OMEGA_R_M), _spread_pct(g, THEORY_GAMMA_M)
        rmed = _median(residuals)
        drifts = Float64[]
        if history !== nothing
            for r in selected
                tlo = r["Tlo_over_M0"]*history.M0
                thi = r["Thi_over_M0"]*history.M0
                mlo = _linear_sample(history.ts, history.mass, tlo)
                mhi = _linear_sample(history.ts, history.mass, thi)
                if mlo !== nothing && mhi !== nothing
                    push!(drifts, 100*abs(mhi-mlo)/history.Mf)
                end
            end
        end
        mass_drift = length(drifts) == length(selected) ? maximum(drifts) : NaN
        period = 2pi * (history===nothing ? 1.0 : history.Mf) / THEORY_OMEGA_R_M
        M0 = history===nothing ? 1.0 : history.M0
        cycles = _median([M0*(r["Thi_over_M0"]-r["Tlo_over_M0"])/period for r in selected])
        flags = String[]
        length(selected) >= min_windows || push!(flags,"fewer than $min_windows fitted windows")
        abs(eR) <= freq_tolerance_pct || push!(flags,"real frequency differs by $( _fmt(abs(eR);digits=1) )%")
        abs(eI) <= freq_tolerance_pct || push!(flags,"damping rate differs by $( _fmt(abs(eI);digits=1) )%")
        isfinite(sR) && sR <= spread_tolerance_pct || push!(flags,"real frequency changes substantially across windows")
        isfinite(sI) && sI <= spread_tolerance_pct || push!(flags,"damping rate changes substantially across windows")
        rmed <= max_fit_residual || push!(flags,"median relative fit residual is too high")
        isfinite(mass_drift) || push!(flags,"horizon-mass stability could not be checked")
        isfinite(mass_drift) && mass_drift <= max_horizon_mass_drift_pct ||
            (isfinite(mass_drift) && push!(flags,"horizon mass changes during fit windows"))
        cycles >= 1.0 || push!(flags,"typical fit window contains < 1 predicted oscillation")
        status = isempty(flags) ? "Plausible candidate (not confirmed)" : "Inconclusive / needs inspection"
        push!(reasons_by_observer,"J/M0=$(_fmt(J_over_M0;digits=1)): $status" *
             (isempty(flags) ? "" : " — " * join(flags,"; ")))
        push!(summary_rows,(J=J_over_M0,n=length(selected),omegaR=wmed,gamma=gmed,
                            errorR=eR,errorI=eI,spreadR=sR,spreadI=sI,
                            fit_error=rmed,mass_drift=mass_drift,cycles=cycles,status=status))
    end

    open(summary_path,"w") do io
        println(io,"J_over_M0,windows,Mf_omegaR_median,Mf_gamma_median,omegaR_error_pct,gamma_error_pct,omegaR_window_spread_pct,gamma_window_spread_pct,median_relative_fit_residual,max_horizon_mass_drift_pct,median_predicted_cycles,status")
        for x in summary_rows
            vals = (@sprintf("%.8g",x.J),string(x.n),@sprintf("%.10g",x.omegaR),
                    @sprintf("%.10g",x.gamma),@sprintf("%.8g",x.errorR),
                    @sprintf("%.8g",x.errorI),@sprintf("%.8g",x.spreadR),
                    @sprintf("%.8g",x.spreadI),@sprintf("%.8g",x.fit_error),
                    @sprintf("%.8g",x.mass_drift),@sprintf("%.8g",x.cycles),x.status)
            println(io,join(vals,","))
        end
    end

    open(report_path,"w") do io
        println(io,"# Scalar l=0 QNM interpretation — MMIEF")
        println(io)
        println(io,"Run: `$(run_dir)`")
        println(io)
        println(io,"## Theoretical Schwarzschild benchmark")
        println(io)
        println(io,"For the fundamental massless scalar mode (l=0) with time dependence exp(-i omega T):")
        println(io)
        println(io,"- M_f omega_R = **$(THEORY_OMEGA_R_M)**")
        println(io,"- M_f Gamma = **$(THEORY_GAMMA_M)** (omega = omega_R - i Gamma)")
        println(io,"- Period = **$(_fmt(2pi/THEORY_OMEGA_R_M;digits=3)) M_f**; e-folding time = **$(_fmt(1/THEORY_GAMMA_M;digits=3)) M_f**")
        if history !== nothing
            println(io,"- Saved final horizon mass M_f = **$(_fmt(history.Mf;digits=9))** (initial M_0=$(_fmt(history.M0;digits=9)))")
        else
            println(io,"- Mass history unavailable. The fit CSV was already normalized by M_f, but stationarity and time-window duration cannot be verified.")
        end
        println(io)
        println(io,"## Numerical results (median across available fit windows)")
        println(io)
        println(io,"| J/M0 | Windows | M_f omega_R | M_f Gamma | error R | error damping | spread R | spread damping | fit residual | mass drift |")
        println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for x in summary_rows
            println(io,"| $(_fmt(x.J;digits=1)) | $(x.n) | $(_fmt(x.omegaR;digits=6)) | $(_fmt(x.gamma;digits=6)) | $(_fmt(x.errorR;digits=2))% | $(_fmt(x.errorI;digits=2))% | $(_fmt(x.spreadR;digits=2))% | $(_fmt(x.spreadI;digits=2))% | $(_fmt(x.fit_error;digits=3)) | $(_fmt(x.mass_drift;digits=3))% |")
        end
        if isempty(summary_rows)
            println(io)
            println(io,"**No valid QNM fits found.** Verify that the run is long enough, observer data are saved, and the fitter produced rows.")
        end
        println(io)
        println(io,"## Interpretation")
        println(io)
        for line in reasons_by_observer
            println(io,"- ",line)
        end
        if isempty(reasons_by_observer)
            println(io,"- No conclusion: the comparison file has no usable fits.")
        end
        if length(summary_rows)>=2
            cross_R = _spread_pct([x.omegaR for x in summary_rows],THEORY_OMEGA_R_M)
            cross_I = _spread_pct([x.gamma for x in summary_rows],THEORY_GAMMA_M)
            println(io)
            println(io,"Across-observer spread: omega_R=$(_fmt(cross_R;digits=2))% of theory; Gamma=$(_fmt(cross_I;digits=2))% of theory.")
            if cross_R > spread_tolerance_pct || cross_I > spread_tolerance_pct
                println(io,"**Warning:** fitted frequencies are not consistent across observer radii at the chosen $(_fmt(spread_tolerance_pct;digits=1))% spread threshold.")
            end
        end
        println(io)
        println(io,"Comparison criteria (heuristic, **not** statistical confidence intervals):")
        println(io,"- Absolute error <= $(_fmt(freq_tolerance_pct;digits=1))% for each dimensionless frequency component")
        println(io,"- Full range of window-to-window estimates <= $(_fmt(spread_tolerance_pct;digits=1))% of the theoretical component")
        println(io,"- Median normalized RMS residual <= $(_fmt(max_fit_residual;digits=2))")
        println(io,"- Maximum horizon-mass drift <= $(_fmt(max_horizon_mass_drift_pct;digits=2))% of M_f over the fit windows")
        println(io,"- At least $min_windows windows and preferably more than one expected oscillation per window")
        println(io)
        println(io,"**Caution:** The scalar l=0 QNM damps rapidly (~9.53 M_f), so single-mode fits over short and overlapping windows can be degenerate. Even a close match does not demonstrate a genuine QNM without grid refinement, boundary independence, fit-window stability, and consistent results at multiple observer radii. Different observer distances introduce travel-time shifts; the checker uses a Schwarzschild tortoise-coordinate approximation.")
    end

    @printf("\nScalar l=0 Schwarzschild prediction: Mf*omegaR=%.6f, Mf*Gamma=%.6f\n",
            THEORY_OMEGA_R_M,THEORY_GAMMA_M)
    if isempty(reasons_by_observer)
        println("No valid fitted frequencies were found; inspect qnm_comparison.csv.")
    else
        for line in reasons_by_observer
            println(line)
        end
    end
    if length(summary_rows)>=2
        rvalues=[x.omegaR for x in summary_rows]
        ivalues=[x.gamma for x in summary_rows]
        @printf("Across-observer median-frequency ranges: omegaR %.2f%%, Gamma %.2f%% of theory\n",
                _spread_pct(rvalues,THEORY_OMEGA_R_M),_spread_pct(ivalues,THEORY_GAMMA_M))
    else
        println("Fewer than two observers with valid fits: no radius-consistency check possible.")
    end
    println("Written: ",report_path)
    println("Written: ",summary_path)
    return (report=report_path,summary=summary_path,observers=length(summary_rows))
end
end # module

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) <= 1 || error("Usage: julia --project=. scripts/interpret_qnm_results.jl [results/standard]")
    QNMInterpretation.interpret_qnm(isempty(ARGS) ? "results/standard" : ARGS[1])
end
