# Perturbation/QNM/tail validation matrix

| Question | Minimum evidence |
|---|---|
| Does the code solve the scalar wave equation? | Analytic 1D wave test and fixed-background Schwarzschild benchmark |
| Is apparent ringdown a QNM? | Stable complex frequency, low residual, window scans, multiple observers and resolution tests |
| Is damping physically accurate? | Fit data above noise, stable `Gamma`, mass normalization and model sensitivity |
| Is a late-time tail resolved? | Extended power-law regime, local power index plateau, outer-boundary and resolution tests |
| Is coupled EKG backreaction relevant? | Compare fixed vs self-gravitating runs at several pulse amplitudes, monitor mass and constraints |

## Common reference values

- Fundamental scalar Schwarzschild `l=0`: `M omega=0.110455-0.104896 i`.
- Generic scalar fixed-areal-radius `l=0` Price tail: `T^-3`.

Both are *benchmark values under appropriate assumptions*, not code
outputs. A synthetic signal agreeing with a target cannot validate a PDE
solver; it validates only parts of an estimator.

## Reproducibility record

For each figure/result include: source revision, background mass/final mass,
initial data, evolution and observer time coordinates, extraction radii,
numerical domain/time step, boundary scheme, dissipation, chosen fit windows,
fit model and normalized residual, resolution and `rmax` studies, final
constraint/mass diagnostics and clear physical caveats.
