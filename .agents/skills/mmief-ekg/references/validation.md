# MMIEF test hierarchy and scientific validation

## Distinguish verification types

- **Static review:** conventions and source interfaces appear consistent.
- **Software unit/integration checks:** algorithms run with expected outputs.
- **Physics validation:** independent benchmarks and convergent observations.
- **Scientific interpretation:** QNM/tail parameter inference with uncertainty.

One level does **not** imply the next. Julia and GLMakie execution may be
unavailable in an AI sandbox; make that explicit.

## Practical checklist

1. `Pkg.test()` with isolated core Julia environment.
2. `--smoke --vacuum`: verify stationary Schwarzschild, `f=0`,
   `fdot=0`, `M_h=M0`, and small constraint residuals.
3. Short weak-field run with finite fields, sensible horizon mass increase
   bounded by outer mass, and no observer leaving the domain.
4. Method-of-lines/RK4 test on an independent known 1D wave solution;
   compare error slopes under refinement.
5. Fixed-background scalar-on-Schwarzschild run for QNM and tail tests,
   independent of backreaction.
6. Convergence tests at multiple `Delta r`, with compatible time steps,
   for scalar observers, mass, and constraints.
7. Outer-boundary sensitivity: vary `rmax` and boundary scheme.
8. QNM extraction: test synthetic decaying sinusoids and mixtures with
   prompt transient and tail, scan fitting intervals and observer radii.
9. Price-law extraction: use several observers, consistent late-time windows,
   local power index, noise floor and finite-time extrapolation checks.

## Known physical and computational caveats

- Boundary closures / constraint projection may reduce global PDE order.
- Schwarzschild well balancing deliberately modifies numerical RHS.
- Outer Sommerfeld conditions are approximate and generate reflections.
- Fitted scalar l=0 QNM damping is difficult: amplitude attenuates by
  roughly `exp(-0.104896*(2pi/0.110455)) ~= 0.00256` per cycle.
- Prior exploratory QNM reports differed by ~15–18% in real frequencies
  and ~53–61% in damping; windows were unstable. This is an **inconclusive**
  numerical finding, not an established continuum discrepancy.
- `scripts/check_qnm_vs_theory.jl` in the baseline archive contains a
  tortoise-based observer delay; stationary MMIEF null travel requires
  `Delta T = Delta r+4M ln[(r2-2M)/(r1-2M)]`.
- The final horizon mass may differ from initial `M0`, affecting QNM units.
- The compiled standalone ADM LaTeX notes are useful mathematical
  documentation; reproducing their formulas in software is a separate test.

## Reporting template

For each scientific result record:

- Source revision / git commit and run configuration.
- `Delta r`, time step/CFL, `rmax`, dissipation, boundary condition.
- Initial pulse amplitude, center, width and scalar cutoff.
- Observer radius `J_obs`, sampled time coordinate, and mass normalization.
- Fitting interval, model, residual, window scan and uncertainty.
- Constraint norms, apparent-horizon residual, mass budget.
- Cross-resolution or cross-domain comparison.
- Explicit verdict: verified / consistent but provisional / inconclusive.
