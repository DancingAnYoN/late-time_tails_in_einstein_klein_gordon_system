---
name: black-hole-perturbations
description: Study or numerically extract black-hole quasinormal modes, ringdown, perturbation potentials, wave scattering and Price-law late-time tails. Use for scalar waves on Schwarzschild, Regge–Wheeler/Zerilli comparisons, horizon/infinity boundary conditions, fit-window selection, spectral analysis and QNM convergence tests.
---

# Black-hole perturbations and waveform analysis

This skill is portable: derive the perturbation equation from the model
before applying quoted Schwarzschild benchmark frequencies or tail exponents.

## Read the relevant reference

- [Schwarzschild scalar equation and effective potential](references/schwarzschild-scalar.md).
- [QNM definitions, extraction and uncertainty](references/qnm-extraction.md).
- [Price-law tails and local power index](references/late-time-tails.md).
- [Numerical validation and observer/time choices](references/validation.md).

## Workflow

1. Identify the perturbation field, spin, angular multipole `l`, spacetime
   dimension, asymptotics, time coordinate and background parameters.
2. Write the governing 1+1 wave equation and effective potential; verify
   the horizon and asymptotic boundary conditions and field normalization.
3. Evolve a fixed-background reference before analyzing nonlinear coupled data.
4. Plot signed waveforms and `log|waveform|` at multiple physical observers;
   separate prompt response, QNM-like ringdown and late-time tails.
5. Scan multiple fit intervals and initial data, record residuals and
   conditioning, and compare frequency estimates across observers.
6. Report both `omega_R` and `Gamma=-Im(omega)>0` in units of final `M`, and
   uncertainty/systematics, not only one best-fit number.
7. Check resolution, outer-boundary sensitivity, late-time noise floor and
   finite-time effects before announcing agreement with theory.

## Cautions

- Scalar `l=0` Schwarzschild frequency differs substantially from
  gravitational `l=2`; never conflate the two.
- `exp(-i omega t)` means damped modes have `Im(omega)<0`.
- Early-time transients may resemble an oscillation; fitting a short
  fraction of a cycle can be ill-conditioned.
- A logarithmic plot near zero crossings can mislead; inspect signed data
  and an envelope separately.
- Price-law powers depend on spin/multipole, data class, observer location
  (fixed radius vs null infinity), background and nonlinearities.
