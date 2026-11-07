using Test
using MMIEF
include("unit/test_grid_and_derivatives.jl")
include("unit/test_rk4.jl")
include("integration/test_method_of_lines.jl")
include("unit/test_geometry_and_constraints.jl")
include("unit/test_scalar_and_observers.jl")
include("unit/test_boundaries_and_analysis.jl")
include("integration/test_evolution.jl")
include("integration/test_data_output.jl")
include("regression/test_previous_bugs.jl")
include("regression/test_legacy_equivalence.jl")
println("Core package tests complete. Run test/visual_smoke.jl with --project=viz separately.")
