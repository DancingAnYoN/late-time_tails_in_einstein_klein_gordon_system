# Requires visualization environment; see README.
#   julia --project=viz scripts/render_movie.jl results/quick
import MMIEF
include(joinpath(@__DIR__,"..","viz","Visualization.jl"))
import .MMIEFVisualization

function movie_main(args::Vector{String}=ARGS)
    isempty(args) && error("Usage: render_movie.jl <run_directory> [--show]")
    run_dir=abspath(args[1])
    all(x->x=="--show",args[2:end]) || error("Unknown render option")
    run=MMIEF.load_run(run_dir)
    return MMIEFVisualization.make_pulse_movie(run;outdir=run_dir,
                                               show_window=("--show" in args))
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    movie_main(ARGS)
end
