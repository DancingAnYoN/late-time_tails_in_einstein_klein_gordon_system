# Project context — late-time tails in the Einstein–Klein–Gordon system

> Portable, model-independent context snapshot for a new AI conversation or
> coding-agent session. **Source baseline:** the user's
> `late-time_tails_in_einstein_klein_gordon_system.zip`, repackaged as the
> GitHub source archive. **Snapshot date:** 2026-11-09.
>
> This is a guide to the current repository, not a certification that its
> numerical predictions agree with theoretical black-hole perturbation theory.

## 1. Research goal

Study a weak, ingoing, spherically symmetric, massless scalar pulse interacting
with a Schwarzschild black hole **with gravitational backreaction**. Follow
absorption, horizon growth, scalar scattering, QNM ringdown and late-time tails.
Use a partially constrained 3+1 evolution in modified ingoing
Eddington–Finkelstein (MMIEF) coordinates, implemented in Julia with Method of
Lines, fourth-order Runge–Kutta and GLMakie visualizations.

The primary reference is R. L. Marsa and M. W. Choptuik,
*Black-hole–scalar-field interactions in spherical symmetry*, Phys. Rev. D **54**
4929 (1996), https://arxiv.org/abs/gr-qc/9607034.
The repository is a subsequent modular implementation, not an assertion of
exact identity with the paper's discretization.

## 2. Mathematical conventions

- Natural units `G=c=1`; signature `(-,+,+,+)`; MTW curvature conventions.
- Foliation normal `n^mu=(1/alpha,-beta/alpha,0,0)`.
- Extrinsic curvature `K_ij=-(1/2) L_n gamma_ij`.
- `t` = evolution coordinate; `r` = fixed uniform numerical coordinate.
- `J=r+f(t)` = **physical areal radius**; `T=t-f(t)` = exported observer time.
- Fixed-`J` observers are interpolated at moving `r=J_obs-f(t)`.

## 3. Geometric and matter variables

The ADM metric is

$$
 ds^2=-\alpha^2dt^2+a^2(dr+\beta\,dt)^2+J^2d\Omega^2,
 \qquad J(t,r)=r+f(t).
$$

Use `A=K^r{}_r`, `B=K^theta{}_theta=K^phi{}_phi`, and `K=A+2B`.
Scalar reduction:

$$
 \Phi=\partial_r\varphi, \qquad
 \Pi=\frac{a}{\alpha}(\partial_t\varphi-\beta\Phi).
$$

The time-evolved state is `(a,B,Phi,Pi,varphi,f)`; auxiliary fields
`A,alpha,beta` are reconstructed. The physical equations are
`G_{mu nu}=8 pi T_{mu nu}`, `square_g varphi=0`, with
`T_{mu nu}=d_mu varphi d_nu varphi - (1/2)g_{mu nu}(d varphi)^2`.

## 4. Gauge, reconstruction and horizon locking

$$
\alpha=a(1-\beta),\quad J'=1,\quad\dot J=\dot f,
\quad v=t+r\text{ null};
\quad c_{\rm in}=-1,\ c_{\rm out}=1-2\beta.
$$

$$
\begin{aligned}
A &= B+JB'-\frac{4\pi J}{a}\Phi\Pi,\\
\dot f&=\left.\frac{4\pi J^2}{a^2}(\Phi+\Pi)^2\right|_{r_h},\\
\beta&=\frac{\dot f+aJB}{1+aJB},\\
\alpha&=\frac{a(1-\dot f)}{1+aJB}.
\end{aligned}
$$

The apparent-horizon condition is `a_h J_h B_h=1`, with fixed numerical
`r_h=2M0`, physical radius `J_h=2M0+f`, and horizon mass `M_h=J_h/2`.

## 5. Evolution equations (dot = d/dt, prime = d/dr)

$$
\begin{aligned}
\dot a &= -a^2(1-\beta) A+(a\beta)',\\
\dot B &= \beta B'+a(1-\beta)B(A+2B)
  +\frac{1-\beta}{J^2}\left(a-\frac1a\right)+\frac{\beta'}{aJ},\\
\dot\Phi &= [\beta\Phi+(1-\beta)\Pi]',\\
\dot\Pi &= \frac{1}{J^2}[J^2(\beta\Pi+(1-\beta)\Phi)]'
                   -\frac{2\dot f}{J}\Pi,\\
\dot\varphi&=\beta\Phi+(1-\beta)\Pi,\\
\dot f &= \left.4\pi J^2(\Phi+\Pi)^2/a^2\right|_{r_h}.
\end{aligned}
$$

## 6. Constraints and mass

Hamiltonian and radial momentum constraints:

$$
\begin{aligned}
0&=a'+\frac{a^3-a}{2J}+\frac{a^3 J}{2} B(2A+B)
                -2\pi Ja(\Phi^2+\Pi^2),\\
0&=B'+\frac{B-A}{J}-\frac{4\pi}{a}\Phi\Pi,\\
0&=\Phi-\partial_r\varphi.
\end{aligned}
$$

Misner–Sharp mass:

$$
m=\frac J2\left(1-\frac1{a^2}+J^2B^2\right),\qquad
m'=4\pi J^2\left[\frac{\Phi^2+\Pi^2}{2a^2}
+\frac{JB}{a}\Phi\Pi\right].
$$

Reconstruct `a=[1-2m/J+J^2 B^2]^(-1/2)` from the integrated mass
constraint at RK4 stages. `A` is reconstructed from the momentum constraint.
The evolution is **partially constrained**, not unconstrained ADM.

## 7. Initial/boundary data and reference vacuum

- Set `M0=1`, `f(0)=0`, grid begins at `r_h=2M0`.
- Scalar `varphi_0=epsilon chi(r) exp[-(r-r0)^2/sigma^2]`,
  `Phi_0=varphi_0'`, `Pi_0=Phi_0+varphi_0/[r(1-beta_0)]` (approximately ingoing).
- Initial geometry solves coupled radial constraints; do not assume Schwarzschild
  geometry is exact once the scalar amplitude is nonzero.
- Vacuum Schwarzschild: `a=sqrt(1+2M/r)`, `beta=2M/(r+2M)`,
  `B=2M(r+2M)/[r(r+2M)]^(3/2)`,
  `A=-2M(r+M)/[r(r+2M)]^(3/2)`.
- Outer boundary uses approximate outgoing/Sommerfeld conditions, including
  `(1-beta)(Pi+Phi)+(1-2beta+fdot) varphi/J ~= 0`.
- Check boundary reflections and constraint projection effects; no rigorous
  continuum transparent boundary has been implemented.

## 8. Numerical method and code map

- Julia 1.10+ compatible package (`Project.toml`), typically used in Julia 1.12
  under Windows and VS Code; core package `MMIEF`.
- Uniform coordinate grid, fourth-order interior differences, explicit RK4,
  Courant factor `C=Delta t/Delta r`, sixth-order KO dissipation, mass projection
  at stages, and Schwarzschild discrete-residual subtraction for vacuum balance.
- `src/Types.jl` and `src/Grid.jl`: state/configuration/grid.
- `src/FiniteDifferences.jl`, `Integrators.jl`, `MethodOfLines.jl`: numerics.
- `src/Geometry.jl`, `EinsteinEquations.jl`, `ScalarField.jl`: geometry/RHS.
- `src/Constraints.jl`, `InitialData.jl`, `BoundaryConditions.jl`: constraints/BCs.
- `src/Evolution.jl`, `Diagnostics.jl`, `Observers.jl`, `DataOutput.jl`: driver,
  monitoring, interpolated fixed-areal-radius sampling, CSV/TOML.
- `src/Analysis.jl`: exploratory single-mode QNM fit, local power index.
- `viz/`: optional GLMakie dashboard/pulse movie environment.
- `scripts/`: CLI simulation, waveform analysis, QNM check/interpretation,
  observer history and QNM graphics.
- `legacy/MMIEFCore.jl`: frozen regression oracle.

Run data schema: `results/<mode>/config.toml`, `metadata.toml`,
`diagnostics.csv`, `observers.csv`, `snapshots.csv`, `final_state.csv`.
The source-only GitHub archive excludes generated result data.

## 9. How to run (project root)

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=. scripts/run_simulation.jl --smoke --vacuum
julia --project=. scripts/run_simulation.jl --quick
julia --project=. scripts/run_simulation.jl --standard
julia --project=. scripts/check_qnm_vs_theory.jl results/standard
```

One-time GLMakie setup, from project root:

```julia
using Pkg
Pkg.activate("viz")
Pkg.develop(path=abspath("."))
Pkg.instantiate()
```

Then:

```powershell
julia --project=viz scripts/render_dashboard.jl results/standard
julia --project=viz scripts/render_movie.jl results/standard
julia --project=viz scripts/render_observer_history.jl results/standard
julia --project=viz scripts/render_qnm_visualization.jl results/standard
```

REPL simulation entry point: `include("run_simulation.jl"); main(["--standard"])`.
REPL movie entry: `include("scripts/render_movie.jl"); movie_main(["results/standard"])`.

## 10. QNM/tail diagnostics and what is currently known

For a weakly perturbed final Schwarzschild black hole, the fundamental
**massless scalar, spherical `l=0`** QNM has
`M_f omega = 0.110455 - 0.104896 i` with `exp(-i omega T)`.
Period `P/M_f ~= 56.89`, e-fold damping time `tau/M_f ~= 9.53`;
signal drops by about a factor of 390 per period.
For generic massless scalar compact initial data the fixed-radius Price tail is
`varphi ~ T^(-3)` (non-exceptional leading contribution).

The user reported an earlier *exploratory* QNM interpretation with median
frequency deviations approximately `15.1%` at `J/M0=10` and `17.6%` at
`J/M0=30`, and damping-rate deviations `53.3%`, `61.2%` respectively;
fits varied substantially across windows. This is evidence of an inconclusive
extraction, **not** verified physical error bars. Source-only archive contains
no raw waveforms, so these numbers cannot be independently recomputed here.

Potential issue to audit: `scripts/check_qnm_vs_theory.jl` in this source
baseline shifts candidate windows by the Schwarzschild tortoise coordinate
`r_* = r+2M ln(r/(2M)-1)`. But at fixed `f=0`, MMIEF `T=t` gives outgoing
travel time `Delta T=Delta r+4M ln[(r2-2M)/(r1-2M)]`.
This coordinate issue may bias comparisons; correct it and test window choices
before interpreting observer-dependent fits. General dynamical spacetime
requires greater care than either stationary delay approximation.

Other risks: QNM frequency/decay search-grid boundaries; few visible cycles;
ringdown-window contamination from prompt response or power-law tails;
resolution, outer boundary and observer interpolation errors.

## 11. Validation priorities, in order

1. `--smoke --vacuum`: horizon shift `f=0`, mass stationary, small constraints.
2. Synthetic QNM fitter test with known frequency/damping and added transients.
3. Fixed Schwarzschild scalar wave test with the same analysis pipeline.
4. Refine grid and vary time step; check constraints, mass budget, waveform
   convergence and moving fixed-`J` observers.
5. Increase `rmax`; check that outgoing boundary reflection does not pollute
   selected QNM and tail intervals.
6. Compare fit across windows and observer locations, after coordinate and
   mass normalization corrections.
7. Only after these checks claim agreement or disagreement with QNM/Price laws.

See `VALIDATION.md`, `.agents/skills/mmief-ekg/references/validation.md` and
`.agents/skills/black-hole-perturbations/references/`.

## 12. References

- Marsa & Choptuik, PRD 54 (1996): https://arxiv.org/abs/gr-qc/9607034
- E. Gourgoulhon, *3+1 Formalism and Bases of Numerical Relativity*:
  https://arxiv.org/abs/gr-qc/0703035
- Existing repository files: `EQUATIONS.md`, `ARCHITECTURE.md`,
  `VALIDATION.md`, `README.md`.
- Extended human-readable derivation: `mmief_ekg_3plus1_derivation.tex`
  from the separate LaTeX notes package (not required by these skills).
