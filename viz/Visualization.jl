module MMIEFVisualization
using MMIEF
using GLMakie
include("Dashboard.jl")
include("PulseMovie.jl")
using .Dashboard: make_dashboard
using .PulseMovie: make_pulse_movie
export make_dashboard, make_pulse_movie
end # module MMIEFVisualization
