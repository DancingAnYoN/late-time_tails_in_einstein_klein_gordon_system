---
name: mmief-ekg
description: Work with the Julia MMIEF solver for the coupled spherically symmetric massless Einstein–Klein–Gordon system. Use when editing its evolution equations, constraints, horizon-locking gauge, RK4 integration, fixed-areal-radius observers, simulations, GLMakie plots, QNM analysis or validation.
---

# MMIEF Einstein–Klein–Gordon workflow

This skill is **project-specific**. At the repository root read `AGENTS.md`,
`PROJECT_CONTEXT.md` and (when relevant) `EQUATIONS.md` before modifying source.
If those files are not available, use the self-contained references below.

## Choose the reference by task

- [Field definitions, gauge and full equations](references/equations.md).
- [Derivation from the ADM 3+1 equations](references/derivations.md).
- [Julia modules, numerical algorithm and CLI commands](references/numerical-methods.md).
- [Scientific validation requirements and known defects](references/validation.md).

## Procedure

1. Identify whether the request is about continuum mathematics, numerical
   implementation, running a simulation, visualization, or fitting observations.
2. Check units, time coordinate and whether observer radii are **areal** (`J`)
   or computational (`r`). Set `J=r+f(t)`; for fixed `J_obs`, sample at
   `r_obs=J_obs-f(t)`.
3. Before changing a field equation, derive it with the specified extrinsic-
   curvature sign and check the Schwarzschild vacuum substitution.
4. Keep **momentum reconstruction** of `A` and **mass projection** of `a`
   distinct from their formal continuum evolution equations.
5. If changing the solver, add focused tests and run vacuum/core checks first.
6. If analyzing QNMs or Price tails, inspect fit windows, mass/time
   normalization, boundary effects, and grid convergence. Never infer a QNM
   merely from apparent oscillations or one fitted frequency.
7. Give precise reproducible Julia commands, output paths, and a statement
   of what has and has not been tested.

## Core conventions

- `(-,+,+,+)`, `G=c=1`, `K_ij=-(1/2) L_n gamma_ij`, MTW curvature.
- `ds²=-alpha² dt²+a²(dr+beta dt)²+J²dOmega²`, `J=r+f`,
  `alpha=a(1-beta)`.
- `A=K^r_r`, `B=K^theta_theta`, `Phi=d_r varphi`,
  `Pi=(a/alpha)(d_t varphi-beta Phi)`.
- Scalar QNM target: `M_f omega_l0=0.110455-0.104896i` in the weak-field
  Schwarzschild limit; fixed-radius generic scalar tail `~T^(-3)`.

## Avoid known mistakes

- Do **not** accidentally use `2f2` in Julia RK4 weighted sums.
- Do **not** mix `T=t-f(t)` with Schwarzschild Killing time for propagation
  corrections without deriving the coordinate transformation.
- Do **not** load GLMakie into the core `src/` package.
- Do **not** call a script-defined REPL function before `include`ing it.
- Do **not** claim physical validation merely because `Pkg.test()` passes.
