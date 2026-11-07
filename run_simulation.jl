# VS Code / Julia REPL convenience entry point.
# Load once: include("run_simulation.jl")
# Then: main(["--smoke", "--vacuum"])
include(joinpath(@__DIR__,"scripts","run_simulation.jl"))
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main(ARGS)
end
