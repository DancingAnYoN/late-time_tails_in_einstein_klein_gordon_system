module Grid
using ..Types: Config
export RadialGrid, check_config, radial_grid

struct RadialGrid
    r::Vector{Float64}
    dr::Float64
end
function check_config(p::Config)
    rh = 2*p.M0
    p.M0 > 0 || error("M0 must be positive")
    p.dr > 0 && p.rmax > rh + 12*p.dr || error("Invalid radial grid")
    nsteps = (p.rmax-rh)/p.dr
    abs(nsteps-round(nsteps)) < 1e-8 ||
        error("(rmax-2M0)/dr must be an integer to place the outer boundary exactly")
    0 < p.cfl <= 0.5 || error("Choose 0 < cfl <= 0.5 for this implementation")
    p.tmax > 0 && p.output_dt > 0 && p.snapshot_dt > 0 || error("Invalid time settings")
    p.pulse_sigma > 0 && p.pulse_amplitude >= 0 || error("Invalid pulse")
    rh < p.cutoff_start < p.cutoff_end < p.pulse_r0 || error("Cutoff must lie between horizon and pulse")
    all(J -> rh < J < p.rmax, p.observers) || error("Observers must lie in domain")
    0 <= p.ko_strength <= 0.15 || error("Use a small KO strength")
    p.profile_stride >= 1 || error("profile_stride must be positive")
    p.mass_budget_fraction >= 0 || error("mass_budget_fraction must be nonnegative")
end

function radial_grid(p::Config)
    check_config(p)
    N = floor(Int, (p.rmax - 2*p.M0) / p.dr) + 1
    return [2*p.M0 + (i - 1) * p.dr for i in 1:N]
end

RadialGrid(p::Config) = RadialGrid(radial_grid(p), p.dr)
end # module Grid
