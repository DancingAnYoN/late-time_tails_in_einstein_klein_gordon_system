# Extracting and interpreting quasinormal ringing

## Conventions

For time dependence `exp(-i omega T)`, write

$$\omega=\omega_R-i\Gamma,\quad\Gamma>0,\qquad
\varphi(T)\approx A e^{-\Gamma(T-T_0)}
\cos[\omega_R(T-T_0)+\delta].$$

At fixed areal radius the weak-field Schwarzschild scalar `l=0` target is

$$M_f\omega_R=0.110455,\qquad M_f\Gamma=0.104896.$$

Thus the period is `P/Mf=2pi/0.110455~=56.89` and damping time is
`tau/Mf=1/0.104896~=9.53`; the amplitude shrinks by about `390`
per complete cycle. A single 50–60 `M` interval may contain almost no
usable complete cycle after prompt contamination.

## Suggested extraction procedure

1. Load **signed** waveform and observer metadata; preserve sampling times.
2. Determine the physical observer radius, time convention and final-mass
   normalization. Check whether final horizon mass has settled in ringdown.
3. Plot prompt transient, candidate ringdown and later tail on both linear
   and `log|field|` axes, but avoid fitting across zero crossings blindly.
4. Use multiple time intervals; avoid extrapolating from fewer than a
   useful number of oscillations unless using specialized, uncertainty-aware
   complex-exponential methods and verifying conditioning.
5. Fit `e^(-Gamma t)(C cos(omega_R t)+S sin(omega_R t))` with amplitudes
   solved linearly for each candidate `(Gamma,omega_R)`. Record parameter
   boundaries, grid search resolution and relative residual.
6. If necessary, compare nonlinear least-squares, matrix-pencil/Prony
   and multi-overtone models on controlled synthetic data.
7. Scan window start/end, observer radius, grid resolution and outer boundary.
   Report spread/systematics rather than only best-fit point estimates.
8. Cross-check fixed Schwarzschild scalar runs *before* coupled EKG.

## Recommended diagnostics

- `delta_R=100*(Mf*omega_R_num/0.110455 - 1)%`.
- `delta_I=100*(Mf*Gamma_num/0.104896 - 1)%`.
- Relative residual `||fit-data||_2/||data||_2` in the chosen window.
- Fit-window stability (`omega_R`, `Gamma` vs start/end).
- Horizon mass drift during fit.
- Multiple-observer consistency and allowed travel-time shift.
- Absolute amplitude and noise floor at ringdown times.

## Failure modes

- Differentiate `M0` from `Mf`; using wrong mass scales both components.
- Window contains the prompt signal or a polynomial tail.
- Damping is fitted against numerical noise; single-overtone model inadequate.
- Search optimum lies at the edge of `(Gamma,omega_R)` domain.
- Observer time coordinate is mistaken for Schwarzschild Killing time.
- A QNM plot appears convincing but fails time-window or resolution tests.

The QNM asymptotic boundary-value problem and the inverse problem of
estimating QNMs from a finite waveform are distinct; treat them separately.
