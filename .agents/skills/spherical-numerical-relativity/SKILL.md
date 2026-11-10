---
name: spherical-numerical-relativity
description: Derive, implement, explain, or validate 3+1 general relativity equations in spherical symmetry, including ADM constraints, lapse/shift gauges, scalar and fluid matter sources, horizons, Misner–Sharp mass, finite differences and Method of Lines/RK4 schemes. Use for mathematical derivations or numerical PDE projects beyond the MMIEF repository.
---

# Spherical numerical relativity

This is a **portable** general-purpose skill, not tied to the MMIEF codebase.
Always establish metric signature, curvature convention, extrinsic-curvature
sign and gauge before writing equations.

## Reference selection

- [ADM and matter projections](references/adm-foundations.md) for `3+1` setup.
- [Spherical reductions and constraints](references/spherical-geometry.md)
  for explicit spatial curvature, `A,B` and mass aspect.
- [Method of Lines and stability](references/numerical-methods.md)
  for discrete evolutions, boundaries and projection.
- [Validation and diagnostics](references/validation.md) for scientific checks.

## Procedure for an equation derivation

1. State the metric, normal, induced metric, `K_ij` convention and all
   coordinate conditions; never silently assume areal gauge `J=r`.
2. Calculate inverse metric, determinant, spatial connection and curvature.
3. Project stress-energy to `rho`, `j_i`, `S_ij` with stated signs.
4. Write the full ADM constraints and evolution, then impose spherical
   symmetry, *then* gauge conditions. Show intermediate contractions.
5. Classify each resulting relation as evolution, constraint, gauge,
   algebraic reconstruction or boundary data.
6. Check vacuum/known exact solutions and dimensional consistency.

## Procedure for a numerical implementation

1. Choose primary state variables and constraints to evolve or project.
2. Express the spatial differential operator and characteristic structure.
3. Keep numerics (grid, derivatives, RK4) separate from physics (RHS,
   constraint projection, boundary conditions).
4. Verify discrete derivatives and RK4 on known analytic solutions.
5. Test constraint propagation, energy/mass balance, horizon motion,
   convergence and outer-boundary contamination.
6. Clearly distinguish formal interior order from measured global convergence.

When a specific project supplies different sign conventions or gauges,
re-derive formulas rather than blindly copying a formula from this skill.
