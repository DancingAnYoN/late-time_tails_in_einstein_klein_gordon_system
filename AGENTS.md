# Agent instructions — MMIEF Einstein–Klein–Gordon project

These repository-level instructions apply to code changes, numerical experiments,
and scientific explanations involving this repository. Use the skill
`.agents/skills/mmief-ekg/SKILL.md` when working specifically on this implementation;
use the two more general skills for reusable derivations or algorithms.

## Scope and scientific status

- The repository implements a **partially constrained**, spherically symmetric,
  massless Einstein–Klein–Gordon (EKG) evolution in modified ingoing
  Eddington–Finkelstein (MMIEF) coordinates, inspired by Marsa–Choptuik (1996).
- This is a **research prototype**, not a convergence-validated QNM/Price-tail solver.
- Never represent a visual match, a grid-search fit, or a passing unit test as
  physical convergence or independent verification. Distinguish *derived*,
  *implemented*, *tested*, and *not yet validated* claims.
- Read `PROJECT_CONTEXT.md` for the current project-specific snapshot and
  `EQUATIONS.md`, `VALIDATION.md`, `ARCHITECTURE.md` for code-aligned documentation.

## Mandatory mathematical conventions

- Signature `(-,+,+,+)`, units `G=c=1`, MTW curvature conventions.
- Extrinsic curvature `K_ij = -(1/2) L_n gamma_ij`.
- Metric `ds²=-alpha² dt²+a²(dr+beta dt)²+J² dOmega²`.
- Gauge: `J=r+f(t)` and `alpha=a(1-beta)`; `r` is the fixed numerical grid.
- `A=K^r_r`, `B=K^theta_theta`, `Phi=partial_r(varphi)`,
  `Pi=(a/alpha)(partial_t(varphi)-beta Phi)`.
- The scalar is `varphi`, not a generic wave variable `u` when discussing physics.
- Courant number `C=Delta t/Delta r` (reserve capital C for it).
- Observer time `T=t-f(t)`; a fixed areal radius `J_obs` corresponds to
  **moving grid coordinate** `r_obs(t)=J_obs-f(t)`.

## Code architecture and commands

- Source: Julia package `src/MMIEF.jl`; CLI entry: `scripts/run_simulation.jl`.
- MATLAB/NumPy-style reasoning is not a substitute for verifying Julia semantics.
- Core environment: `julia --project=. scripts/run_simulation.jl --smoke --vacuum`.
- Standard run: `julia --project=. scripts/run_simulation.jl --standard`.
- Core tests: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- GLMakie is in **separate** `viz/` environment. First-time setup, from root:
  `julia --project=viz -e 'using Pkg; Pkg.develop(path=abspath(".")); Pkg.instantiate()'`.
- Movie: `julia --project=viz scripts/render_movie.jl results/standard`.
- Observer history: `julia --project=viz scripts/render_observer_history.jl results/standard`.
- QNM comparison: `julia --project=. scripts/check_qnm_vs_theory.jl results/standard`.
- QNM plots: `julia --project=viz scripts/render_qnm_visualization.jl results/standard`.
- A new clone does not contain run data; generate it first. Generated outputs live
  under `results/` and should not be committed by default.

## Safe modification workflow

1. Identify the relevant equations in `EQUATIONS.md`, and inspect the corresponding
   `src/` implementation before changing either.
2. Preserve physics conventions and the boundary/horizon contracts; state any
   changed assumptions explicitly and update the reference documentation.
3. Keep generic MoL/RK4 numerics independent from physical constraint projection.
4. Do not import GLMakie in `src/`; retain a separate visualization environment.
5. Add/adjust focused tests. For physics-affecting changes, run vacuum,
   constraints, convergence and outer-boundary checks as appropriate.
6. Do not edit `legacy/MMIEFCore.jl` casually; it is a regression oracle.
7. Avoid Julia `2f2` and `2f3` in RK4 sums (Julia may parse them as Float32
   literals); use `2*f2`, `2*f3`.
8. Avoid including freshly defined visualization modules inside a function and
   immediately invoking them (Julia world-age issue). Prefer top-level imports.
9. Use `String.(strip.(split(...)))` or generic `AbstractString` interfaces when
   loading CSV headings; `split` returns `SubString` values.
10. If you cannot execute Julia/GLMakie, say so and report only the checks done.

## QNM/tail discipline

- The scalar spherical fundamental Schwarzschild reference is
  `M_f omega = 0.110455 - 0.104896 i` using `exp(-i omega T)`.
  This is **not** the gravitational `l=2` reference.
- Check final-mass normalization, time coordinate, ringdown interval, fit
  residual, parameter-window sensitivity, multiple observers and grid resolutions.
- The MMIEF stationary outgoing null travel time is
  `Delta T = Delta r + 4M ln[(r2-2M)/(r1-2M)]`; do not silently substitute
  the Schwarzschild tortoise delay with coefficient `2M`.
- Generic massless scalar `l=0` fixed-radius Price tail: `varphi ~ T^(-3)`;
  inspect finite-radius/finite-domain reflections and local power index drift.
- See `.agents/skills/black-hole-perturbations/references/`.

## Output/communication requirements

- Include reproducible terminal and REPL commands when providing a workflow.
- Treat simulations and fits as calculations with uncertainties and assumptions.
- Prefer narrow diffs, deterministic input/output schema, and versioned metadata.
- Do not include `.git/`, simulation results, local paths or credential files in
  a GitHub-ready source archive.
