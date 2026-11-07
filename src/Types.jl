module Types
export Config, IA, IB, IPHI, IPI, IFIELD, SimulationState, fieldview
Base.@kwdef struct Config
    M0::Float64 = 1.0
    rmax::Float64 = 600.0
    dr::Float64 = 0.25
    cfl::Float64 = 0.25
    tmax::Float64 = 400.0
    pulse_r0::Float64 = 20.0
    pulse_sigma::Float64 = 2.0
    pulse_amplitude::Float64 = 1e-4
    cutoff_start::Float64 = 6.0
    cutoff_end::Float64 = 9.0
    ko_strength::Float64 = 0.06
    output_dt::Float64 = 0.5
    snapshot_dt::Float64 = 2.0
    profile_stride::Int = 4
    observers::Vector{Float64} = [10.0, 30.0, 50.0]
    # Fail early when the horizon mass exceeds the initial outer mass by
    # far more than can be attributed to truncation error. A weak, isolated
    # ingoing pulse cannot increase the BH mass by many solar/geometric units.
    check_mass_budget::Bool = true
    mass_budget_fraction::Float64 = 0.01
end
const IA, IB, IPHI, IPI, IFIELD = 1, 2, 3, 4, 5

"""Full method-of-lines state: five radial fields and scalar horizon shift."""
mutable struct SimulationState
    fields::Matrix{Float64}
    f::Float64
end
"""Read-only named field-column accessor (without copying)."""
function fieldview(state::SimulationState, name::Symbol)
    index = Dict(:a=>IA, :B=>IB, :Phi=>IPHI, :Pi=>IPI, :phi=>IFIELD)
    haskey(index,name) || throw(ArgumentError("Unknown field $name"))
    return view(state.fields,:,index[name])
end
end # module Types
