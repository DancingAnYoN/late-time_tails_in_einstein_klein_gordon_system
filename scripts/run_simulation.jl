# Julia project runner: NO GLMakie dependency.
# Execute from project root:
#   julia --project=. scripts/run_simulation.jl --smoke --vacuum
# Or REPL: include("scripts/run_simulation.jl"); main(["--smoke"])
using MMIEF
using Printf

function configuration(mode::String, vacuum::Bool)
    if mode=="--smoke"
        return Config(rmax=70.0,dr=0.5,cfl=0.2,tmax=2.0,output_dt=0.5,
                      snapshot_dt=0.5,profile_stride=2,ko_strength=0.045,
                      pulse_amplitude=vacuum ? 0.0 : 1e-4)
    elseif mode=="--quick"
        return Config(rmax=160.0,dr=0.5,cfl=0.2,tmax=65.0,output_dt=0.5,
                      snapshot_dt=1.0,profile_stride=2,ko_strength=0.045,
                      pulse_amplitude=vacuum ? 0.0 : 1e-4)
    elseif mode=="--standard"
        return Config(rmax=700.0,dr=0.25,cfl=0.2,tmax=350.0,output_dt=0.5,
                      snapshot_dt=2.0,profile_stride=4,ko_strength=0.045,
                      pulse_amplitude=vacuum ? 0.0 : 1e-4)
    elseif mode=="--long"
        return Config(rmax=1700.0,dr=0.25,cfl=0.2,tmax=800.0,output_dt=0.5,
                      snapshot_dt=4.0,profile_stride=8,ko_strength=0.045,
                      pulse_amplitude=vacuum ? 0.0 : 1e-4)
    else
        error("Unrecognized mode $mode")
    end
end

function main(args::Vector{String}=ARGS)
    if "--help" in args || "-h" in args
        println("""
        MMIEF Einstein–Klein–Gordon numerical simulation
        Usage: julia --project=. scripts/run_simulation.jl [options]
        --smoke | --quick | --standard | --long   (default --quick)
        --vacuum    Schwarzschild vacuum test
        --help      This help
        Output: results/<mode>[_vacuum]/, containing versioned CSV + TOML.
        Run visualization separately via --project=viz.
        """)
        return nothing
    end
    allowed=Set(["--smoke","--quick","--standard","--long","--vacuum"])
    all(x->x in allowed,args) || error("Unrecognized flag. See --help")
    modes=[x for x in args if x in ("--smoke","--quick","--standard","--long")]
    length(modes)<=1 || error("Specify only one simulation size")
    mode=isempty(modes) ? "--quick" : only(modes)
    p=configuration(mode,"--vacuum" in args)
    runname=mode[3:end]*("--vacuum" in args ? "_vacuum" : "")
    target=joinpath(@__DIR__,"..","results",runname)
    println("Running $runname: dr=$(p.dr), tmax=$(p.tmax), rmax=$(p.rmax)")
    data=run_evolution(p)
    output=save_run(data,target)
    @printf("Saved standardized output: %s\nFinal horizon mass: %.12g\n",
            output,data.masses_h[end])
    return data
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
