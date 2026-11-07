# Refactor validation and scientific limitations

## Software validation checklist

1. `Pkg.test()` — core unit, integration, and regression tests.
2. `julia --project=. scripts/run_simulation.jl --smoke --vacuum`
   — `f` and `fdot` identically zero; horizon mass exactly `M0`.
3. `julia --project=. scripts/run_simulation.jl --smoke`
   — all finite; versioned CSV+TOML output created.
4. `julia --project=viz test/visual_smoke.jl`
   — PNG and MP4 created after reading saved data.
5. Compare `--quick` waveforms/constraint norms against the old monolithic
   implementation at **identical** time steps and parameters.

## Previous defects and how this refactor guards against them

- Invalid GLMakie UUID -> exact registry UUID in `viz/Project.toml`.
- `2f2` interpreted as a Float32 literal -> generic RK4 uses `2*f2`.
- World-age failure after `include` inside main -> packages loaded at top level.
- Fixed-J observer leaving grid -> explicit guard with diagnostic geometry.
- Unphysical horizon growth -> fail-fast budget diagnostic in Evolution.
- Solver and visualization coupling -> separate optional visualization project.
- Data not reusable -> versioned portable CSV and TOML, `load_run`.

## Scientific caveats

The solver is a **partially constrained** MMIEF implementation. `A` comes from
momentum constraint and `a` is projected onto Hamiltonian constraint after RK
stages. Schwarzschild residual subtraction is a well-balancing modification.
This should be validated separately from unprojected ADM evolution. Boundary
stencils and constraint projection may limit the global order of accuracy;
interior fourth-order derivative and RK4 unit tests do not prove overall
fourth-order PDE convergence. The outgoing boundary is approximate and will
produce reflections. The scalar fit routines are exploratory. Validate mass,
constraint behavior, radius-dependent waveforms, resolution convergence and
outer-boundary dependence before interpreting quasinormal ringing and tails.

Julia runtime and GLMakie were unavailable to the build environment; **all
runtime tests must be executed on the user's installation**.
