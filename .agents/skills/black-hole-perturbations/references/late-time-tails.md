# Price-law late-time tails

In the asymptotically flat Schwarzschild exterior, for generic compact
massless scalar initial data and fixed areal radius, the multipole
components exhibit leading decay

$$\varphi_l(T,J_{obs})\sim T^{-(2l+3)}$$

under standard non-exceptional assumptions. In particular `l=0` gives
`T^-3`. Exceptional data classes, different observers (e.g. future null
infinity), non-asymptotically-flat backgrounds, different dimensions and
nonlinearities can alter this exponent. Do not present this as a theorem
for all data without specifying assumptions.

## Local power index

If `|varphi|~T^-p`, then

$$p(T)=-\frac{d\log|\varphi(T)|}{d\log T}
\longrightarrow p.$$

Raw `log|varphi|` becomes singular at zero crossings, so use a defensible
peak envelope or a sliding-window fit with an explicit noise threshold.
For numerical `T_i, y_i`, local regression of `log envelope` against
`log T` gives a more stable diagnostic than two-point differencing.
Avoid a **running maximum from initial time** as a tail envelope:
it is monotone and remembers the early large pulse, so it cannot track
late-time decay. Instead use a localized max or peak sequence.

## Validation

- Show scalar magnitude and sign separately.
- Choose late-time fit interval by observed power-law regime, not a fixed
  window before checking waveform transitions.
- Vary observer radii at fixed physical `J` and test consistency.
- Increase `rmax` and evolve long enough that boundary reflections
  do not precede the intended tail window.
- Vary resolution, KO dissipation and extraction method; estimate the floor.
- Compare against a fixed-background result with known expected exponent.

In a dynamical EKG spacetime, separate linear scattering tails from
backreaction-induced nonlinear contributions whenever possible.
