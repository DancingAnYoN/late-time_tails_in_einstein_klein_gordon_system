# MMIEF — modular Julia Einstein–Klein–Gordon solver

This is a refactor of the existing MMIEF partially constrained research
prototype, **not** a new verified high-order solver. Geometry and physics
follow Marsa–Choptuik; the time integrator is Method of Lines + RK4 with
Hamiltonian projection at intermediate stages and Schwarzschild discrete
residual subtraction. Ingoing null slicing: `alpha=a*(1-beta)`; shifted
areal radius: `J=r+f(t)`.

## 1. Install and run the numerical solver (no GLMakie needed)

Open the folder `MMIEF_Modular` in VS Code, start a Julia REPL, and type:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using MMIEF
include("run_simulation.jl")
main(["--smoke", "--vacuum"])
main(["--quick"])
```

Or from PowerShell inside the project root:

```powershell
julia --project=. scripts/run_simulation.jl --smoke --vacuum
julia --project=. scripts/run_simulation.jl --quick
```

Available simulation sizes are `--smoke`, `--quick`, `--standard`, `--long`;
add `--vacuum` for Schwarzschild vacuum. Results go into
`results/smoke_vacuum/`, `results/quick/`, etc.

**Do not** use `--no-movie` or `--no-gl` with the numerical runner: visualization
is now a separate process. This is intentional.

## 2. Install the separate GLMakie visualization environment

From the root folder in the Julia REPL:

```julia
using Pkg
Pkg.activate("viz")
Pkg.develop(path=abspath("."))
Pkg.instantiate()
```

If you use a fresh Julia REPL, the terminal equivalents are:

```powershell
julia --project=viz scripts/render_dashboard.jl results/quick
julia --project=viz scripts/render_movie.jl results/quick
```

These load **saved numerical results**; they do not rerun the PDE simulation.
Outputs include `dashboard.png` and `scalar_pulse_absorption.mp4`.
Add `--movie` to `render_dashboard.jl` to render the six-panel dashboard MP4.

When manually switching between the core and viz environments, activate the
right one with `Pkg.activate(".")` or `Pkg.activate("viz")`.

## 3. Tests

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

Tests are grouped as `test/unit`, `test/integration`, and `test/regression`.
Core tests do not load GLMakie. For a separate OpenGL-dependent smoke test:

```powershell
julia --project=viz test/visual_smoke.jl
```

## 4. Data format

Every run directory contains:

- `config.toml`: parameters used in the run.
- `metadata.toml`: schema v1, coordinate and unit conventions.
- `diagnostics.csv`: horizon mass, outer mass and constraint norms.
- `observers.csv`: long-format `(t,T_asymptotic,J_obs,phi)` time series.
- `snapshots.csv`: long-format `(t,f,r,J,Jphi)` movie input data.
- `final_state.csv`: complete five-column evolved state at final time.

You can read it from any Julia process with `MMIEF.load_run("results/quick")`.
All CSV values are plain text and intended to be portable. Large `--long`
runs may create sizeable snapshot files; adjust `snapshot_dt` and
`profile_stride` if disk use becomes too high. **There is not yet a restart
checkpoint format**: final_state.csv alone cannot resume integration.

## 5. Source layout

- `src/Types.jl`, `Grid.jl`: configuration, field indexes, grid.
- `src/FiniteDifferences.jl`, `Integrators.jl`, `MethodOfLines.jl`: independent MoL methods.
- `src/Geometry.jl`, `ScalarField.jl`, `EinsteinEquations.jl`:
  continuum formulas in semi-discrete RHS form.
- `src/Constraints.jl`, `InitialData.jl`, `BoundaryConditions.jl`:
  constraint solve/projection, Gaussian, excision & outgoing boundaries.
- `src/Evolution.jl`: coupled RHS and run driver.
- `src/Diagnostics.jl`, `Observers.jl`, `DataOutput.jl`:
  monitoring, fixed-areal-radius sampling, versioned CSV/TOML.
- `src/Analysis.jl`: exploratory QNM and tail estimators.
- `viz/Dashboard.jl`, `PulseMovie.jl`: GLMakie (optional environment).
- `scripts/`: standalone runner, renderer and analysis commands.
- `legacy/MMIEFCore.jl`: frozen monolithic solver used for equivalence tests.

## Validation status

Static source refactoring and independent Python checks were possible in the
creation environment; Julia and GLMakie were not installed there. The included
Julia tests **have not been executed here**. Run `Pkg.test()` and the optional
GLMakie smoke test locally before treating simulation output as reliable.
The original research prototype is also not convergence-validated for late
time Price tails or quasinormal modes. See `VALIDATION.md`.
