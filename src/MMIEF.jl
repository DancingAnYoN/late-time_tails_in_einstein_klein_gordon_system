module MMIEF
# Include every submodule once, at package load time (Julia 1.12 world-age safe).
include("Types.jl")
include("Grid.jl")
include("FiniteDifferences.jl")
include("Geometry.jl")
include("ScalarField.jl")
include("EinsteinEquations.jl")
include("Constraints.jl")
include("InitialData.jl")
include("BoundaryConditions.jl")
include("Integrators.jl")
include("MethodOfLines.jl")
include("Observers.jl")
include("Diagnostics.jl")
include("DataOutput.jl")
include("Evolution.jl")
include("Analysis.jl")

using .Types: Config, SimulationState
using .Grid: radial_grid, RadialGrid
using .Geometry: schwarzschild_fields
using .InitialData: initial_data
using .Diagnostics: diagnostics
using .DataOutput: RunData, save_run, load_run, save_csv
using .Evolution: run_evolution, vacuum_test
using .Analysis: fit_qnm, local_power_index

export Config, SimulationState, RadialGrid, radial_grid, schwarzschild_fields,
       initial_data, diagnostics, run_evolution, vacuum_test, RunData,
       save_run, load_run, save_csv, fit_qnm, local_power_index
end # module MMIEF
