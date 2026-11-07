# Compare MMIEF scalar l=0 ringdown estimates with Schwarzschild theory.
# Run at project root:
#   julia --project=. scripts/check_qnm_vs_theory.jl results/standard
# Or in VS Code Julia REPL:
#   include("scripts/check_qnm_vs_theory.jl")
#   qnm_check("results/standard")
#
# No GLMakie required. Results are exploratory, not a QNM validation.

using MMIEF
using Printf

# Load interpretation code at top-level to avoid Julia 1.12 world-age issues.
include(joinpath(@__DIR__, "interpret_qnm_results.jl"))

const SCALAR_L0_OMEGA_R = 0.110455
const SCALAR_L0_GAMMA   = 0.104896

# Schwarzschild tortoise radius, used only to shift proposed fit windows
# between observers. This is a heuristic delay, not a QNM fit parameter.
function schwarzschild_tortoise(J::Float64, M::Float64)
    J > 2M || throw(ArgumentError("Observer must lie outside final Schwarzschild radius"))
    return J + 2M * log(J / (2M) - 1)
end

function qnm_check(run_directory::AbstractString="results/standard")
    run_dir = abspath(run_directory)
    data = MMIEF.load_run(run_dir)
    p = data.settings
    M0 = p.M0
    Mf = data.masses_h[end]
    isfinite(Mf) && Mf > 0 || error("Final horizon mass must be positive and finite")
    T = data.times .- data.shifts  # Asymptotically normalized time in MMIEF gauge.
    issorted(T) || error("Asymptotic times are not monotone")

    @printf("\nMMIEF scalar l=0 fundamental QNM comparison\n")
    @printf("Run: %s\n", run_dir)
    @printf("Initial mass M0=%.9g, final horizon mass Mf=%.9g\n", M0, Mf)
    @printf("Reference (Schwarzschild): Mf*omega_R=%.6f, Mf*gamma=%.6f\n",
            SCALAR_L0_OMEGA_R, SCALAR_L0_GAMMA)
    @printf("Reference physical frequencies: omega_R=%.9f, gamma=%.9f\n",
            SCALAR_L0_OMEGA_R/Mf, SCALAR_L0_GAMMA/Mf)
    @printf("T span [%.3f, %.3f] M0\n", T[1]/M0, T[end]/M0)
    @printf("NOTE: These single-mode fits must be checked for window and resolution stability.\n")

    # These windows are proposed at J=30 M0; shift by outgoing null travel time
    # to compare approximately corresponding waveform phases at other observers.
    baseline_windows = ((35.0,80.0), (40.0,90.0), (45.0,100.0),
                        (50.0,110.0), (55.0,120.0))
    J_ref = 30*M0
    rstar_ref = schwarzschild_tortoise(J_ref, Mf)
    out_path = joinpath(run_dir, "qnm_comparison.csv")
    open(out_path, "w") do io
        println(io, "J_over_M0,Tlo_over_M0,Thi_over_M0,Mf_omegaR,Mf_gamma,omegaR_error_pct,gamma_error_pct,relative_fit_residual")
        for (j, J) in enumerate(p.observers)
            phi = data.signals[j]
            rstar_shift = schwarzschild_tortoise(J, Mf) - rstar_ref
            @printf("\nObserver J/M0 = %.1f (estimated time shift %.2f M0)\n",
                    J/M0, rstar_shift/M0)
            @printf("  %12s %12s %13s %13s %12s %12s %12s\n",
                    "Tlo/M0", "Thi/M0", "Mf*omegaR", "Mf*gamma", "deltaR %", "deltaI %", "fit err")
            for (left, right) in baseline_windows
                tlo = M0*left + rstar_shift
                thi = M0*right + rstar_shift
                if tlo < T[1] || thi > T[end]
                    @printf("  %12.2f %12.2f  (window outside saved data)\n", tlo/M0, thi/M0)
                    continue
                end
                fit = MMIEF.fit_qnm(T, phi; tlo=tlo, thi=thi, mass=Mf)
                if fit === nothing
                    @printf("  %12.2f %12.2f  (insufficient data / fit failed)\n",tlo/M0,thi/M0)
                    continue
                end
                dR = 100*(fit.omegaM/SCALAR_L0_OMEGA_R - 1)
                dI = 100*(fit.gammaM/SCALAR_L0_GAMMA - 1)
                @printf("  %12.2f %12.2f %13.6f %13.6f %+12.2f %+12.2f %12.3g\n",
                        tlo/M0, thi/M0, fit.omegaM, fit.gammaM, dR, dI, fit.rel_error)
                @printf(io, "%.8g,%.8g,%.8g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                        J/M0, tlo/M0, thi/M0, fit.omegaM, fit.gammaM, dR, dI, fit.rel_error)
            end
        end
    end
    println("\nSaved comparison table: ",out_path)
    println("Do not identify a QNM from numerical proximity alone; require stability across windows, radii and resolutions.")
    QNMInterpretation.interpret_qnm(run_dir)
    return out_path
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) <= 1 || error("Usage: julia --project=. scripts/check_qnm_vs_theory.jl [results/standard]")
    qnm_check(isempty(ARGS) ? "results/standard" : ARGS[1])
end
