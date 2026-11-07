# Architecture and test contracts

## Module boundaries

| Module | Inputs | Outputs | Test | Dependency policy |
|---|---|---|---|---|
| Types | Config parameters | Config, state constants | grid/type tests | Base only |
| Grid | Config | radial grid | unit/grid | Types |
| FiniteDifferences | arrays + dr | derivatives, KO terms | unit/FD | Base only |
| Integrators | semi-discrete RHS callback + state | RK4 accepted step | unit/RK4 | Base only |
| MethodOfLines | RHS + projection + boundary hooks | MOL step + CFL choice | integration/wave | Integrators |
| Geometry | geometry + shift | lapse, beta, mass, fdot | unit/geometry | Types |
| ScalarField | profile, fluxes | Gaussian and wave RHS | unit/scalar | Types |
| EinsteinEquations | geometry + differentiated fields | gravitational RHS | vacuum and integration | Types |
| Constraints | stage fields | A, Hamiltonian projection | unit/constraints | Types |
| InitialData | Config | fields on initial slice | constraints tests | Grid + Physics |
| BoundaryConditions | fields at outer end | outgoing scalar and metric RHS | unit/boundary | Geometry |
| Evolution | Config and all numerical physics modules | RunData | evolution integration | All core modules |
| Observers | field, J_obs, f | phi at fixed J | unit/interpolation | Geometry |
| Diagnostics | evolved fields | mass/constraint diagnostics | unit/geometry + evolution | Geometry |
| DataOutput | RunData | schema v1 CSV/TOML | integration/round-trip | No GLMakie |
| Analysis | waveform arrays | exploratory fits | unit/synthetic signals | Base only |
| Visualization | saved RunData | PNG/MP4 | visual_smoke | GLMakie (optional) |

## Design principles

- Include Julia modules **once** at package initialization. Never `include` a
  new module inside the running `main` function.
- Do not put `using GLMakie` inside `src/`. Visualization is a separate Julia
  environment under `viz/`.
- Keep boundary and projection callbacks explicit in the numerical step; do
  not secretly enforce constraints inside the generic RK4 routine.
- Store physical coordinates `J=r+f` with snapshot profiles and document
  the coordinate time `T=t-f` alongside `t`.
- Disallow silent observer clamping and unphysical mass growth.
- Preserve the original monolithic `legacy/MMIEFCore.jl` strictly as a
  test oracle. Remove only after independent resolution convergence.

## Validation levels

1. **Syntax, dependency, and package checks:** `using MMIEF`, `Pkg.test()`.
2. **Unit tests:** stencils, RK4, geometry, gauge, interpolation, analysis.
3. **Functional tests:** periodic 1D wave equation, stationary Schwarzschild,
   short weak-field Einstein–Klein–Gordon, data output round trip.
4. **Regression tests:** artificial horizon growth; old/new solver comparison.
5. **OpenGL smoke:** render PNG and MP4 from saved data in isolated viz env.
6. **Physics acceptance:** radial/time convergence; horizon/outer-mass
   conservation; boundary independence; QNM and Price-tail fit stability.

Levels 1–5 are software checks, not by themselves sufficient for level 6.
