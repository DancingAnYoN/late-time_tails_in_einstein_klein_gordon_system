# Validation hierarchy for spherical numerical relativity

1. **Symbolic checks:** metric inverse, tensor signs, spherical reductions,
   vacuum substitution, dimensional consistency.
2. **Discrete operator tests:** centered/one-sided differences, RK4 convergence
   on a known ODE, CFL and artificial dissipation tests.
3. **Analytic baseline:** Minkowski, Schwarzschild in the chosen slicing,
   vacuum constraints and zero scalar solution.
4. **Independent PDE baseline:** scalar 1D wave / fixed-background scattering,
   with known analytic or high-accuracy reference.
5. **Constraints and budgets:** Hamiltonian, momentum, scalar definition,
   horizon expansion, Misner–Sharp/ADM mass.
6. **Multiple resolutions:** at least three grid scales, controlled timestep
   scaling, matching physical observers and times; report measured convergence.
7. **Boundaries:** vary outer radius and boundary prescription; identify
   potential reflection contamination via causal travel times.
8. **Gauge robustness:** compare geometric invariant outputs, not merely
   coordinate locations; check lapse positivity/coordinate pathology.
9. **Scientific observables:** compare with literature only after 1–8.

For a quantity `Q_h`, fourth-order Richardson behavior predicts
`Q_h-Q_{h/2} ~= 16 (Q_{h/2}-Q_{h/4})` when errors are asymptotic and
other error sources are controlled. This is a *test*, not a guaranteed
property of a particular implementation.

When reporting, state which levels were actually performed.
