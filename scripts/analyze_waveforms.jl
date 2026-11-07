using MMIEF
using Printf

function analyze_main(args::Vector{String}=ARGS)
    length(args)==1 || error("Usage: analyze_waveforms.jl <run_directory>")
    data=load_run(abspath(args[1]))
    times=data.times.-data.shifts
    M=data.masses_h[end]
    println("Exploratory analysis only — validate convergence and fit windows.")
    for (j,J) in enumerate(data.settings.observers)
        phi=data.signals[j]
        @printf("J/M0 = %.2f; max |phi|=%.5e\n",J/data.settings.M0,maximum(abs,phi))
        fit=fit_qnm(times,phi;tlo=45*data.settings.M0,
                    thi=min(110*data.settings.M0,times[end]),mass=M)
        if fit===nothing
            println("    QNM fit unavailable (short run).")
        else
            @printf("    exploratory M omega_R=%.6f, M gamma=%.6f, residual=%.3e\n",
                    fit.omegaM,fit.gammaM,fit.rel_error)
        end
    end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    analyze_main(ARGS)
end
