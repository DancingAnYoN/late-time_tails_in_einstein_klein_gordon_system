# Equations implemented in MMIEF_RK4_GL

Reference: R. L. Marsa and M. W. Choptuik, *Black-hole–scalar-field
interactions in spherical symmetry*, Phys. Rev. D **54** 4929 (1996),
https://arxiv.org/abs/gr-qc/9607034.

## Conventions and gauge

We set G=c=1, use metric signature (-+++), and the sign convention
`K_ij = -1/2 L_n gamma_ij`.

\[
 ds^2=-\alpha^2dt^2+a^2(dr+\beta\,dt)^2+J^2d\Omega^2,\qquad
 J(t,r)=r+f(t),\qquad \alpha=a(1-\beta).
\]

A dot denotes \(\partial_t\), a prime denotes \(\partial_r\).
The advanced null coordinate is \(v=t+r\); radial null speeds are
\(c_{\rm in}=-1\), \(c_{\rm out}=1-2\beta\).

\[
 A=K^r{}_r,\qquad B=K^\theta{}_\theta,\qquad
 \Phi=\partial_r\varphi,\qquad
 \Pi=\frac{a}{\alpha}(\partial_t\varphi-\beta\Phi).
\]

The lapse, shift and radial curvature component are reconstructed:

\[
\begin{aligned}
 A&=B+JB'-\frac{4\pi J}{a}\Phi\Pi,\\
 \dot f&=\left.\frac{4\pi J^2}{a^2}(\Phi+\Pi)^2\right|_{r=r_h},\\
 \beta&=\frac{\dot f+aJB}{1+aJB},\\
 \alpha&=\frac{a(1-\dot f)}{1+aJB}.
\end{aligned}
\]

The apparent horizon is kept at \(r_h=2M_0\), with
\(a(t,r_h)J(t,r_h)B(t,r_h)=1\) and physical horizon mass \(M_h=J_h/2\).

## Evolution equations

The fields integrated in time are \(a,B,\Phi,\Pi,\varphi,f\):

\[
\boxed{\begin{aligned}
 \dot a&=-a^2(1-\beta)A+(a\beta)',\\
 \dot B&=\beta B'+a(1-\beta)B(A+2B)
     +\frac{1-\beta}{J^2}\left(a-\frac1a\right)
     +\frac{\beta'}{aJ},\\
 \dot\Phi&=[\beta\Phi+(1-\beta)\Pi]',\\
 \dot\Pi&=\frac1{J^2}\left[J^2(\beta\Pi+(1-\beta)\Phi)\right]'
      -\frac{2\dot f}{J}\Pi,\\
 \dot\varphi&=\beta\Phi+(1-\beta)\Pi,\\
 \dot f&=\left.\frac{4\pi J^2}{a^2}(\Phi+\Pi)^2\right|_{r_h}.
\end{aligned}}
\]

The numerical implementation reconstructs \(A\) from the momentum
constraint at each RK4 stage, and optionally monitors the constraints.

## Constraints

\[
\boxed{\begin{aligned}
 0 &= a'+\frac{a^3-a}{2J}+\frac{a^3J}{2}B(2A+B)
             -2\pi Ja(\Phi^2+\Pi^2),\\
 0 &= B'+\frac{B-A}{J}-\frac{4\pi}{a}\Phi\Pi,\\
 0 &= \Phi-\partial_r\varphi.
\end{aligned}}
\]

The spherical Misner–Sharp mass obeys

\[
 m=\frac J2\left(1-\frac1{a^2}+J^2B^2\right),\qquad
 m'=4\pi J^2\left[\frac{\Phi^2+\Pi^2}{2a^2}
                     +\frac{JB}{a}\Phi\Pi\right].
\]

We use the latter mass ODE for Hamiltonian projection at every RK4 stage.
This procedure enforces the constraint algebraically and is **not**
identical to pure unconstrained ADM evolution.

## Initial data

Set \(M_0=1\), \(f(0)=0\), \(r_h=2M_0\).
Start from a Gaussian profile with a smooth cutoff near the horizon:

\[
 \varphi_0=\varepsilon\,\chi(r)
       \exp[-(r-r_0)^2/\sigma^2],\qquad
 \Phi_0=\varphi'_0,\qquad
 \Pi_0=\Phi_0+\frac{\varphi_0}{r(1-\beta_0)}.
\]

The initial metric is not kept exactly Schwarzschild when \(\varepsilon\ne0\).
We solve the two constraints radially with the initial extrinsic-curvature
ansatz

\[
 A_0(r)=-\frac{2m(r)(r+m(r))}{[r(r+2m(r))]^{3/2}},
 \qquad a_0(r)=[1-2m(r)/r+(rB_0(r))^2]^{-1/2}.
\]

At the horizon, in the scalar-free cutoff region,
\(m(r_h)=M_0\), \(B(r_h)=1/(2\sqrt2M_0)\).

The unperturbed Schwarzschild fields are

\[
 a_S=\sqrt{1+2M/r},\quad
 B_S=\frac{2M(r+2M)}{[r(r+2M)]^{3/2}},\quad
 A_S=-\frac{2M(r+M)}{[r(r+2M)]^{3/2}},\quad
 \beta_S=\frac{2M}{r+2M}.
\]

## Outer and inner boundaries

- Inner boundary: fixed apparent-horizon coordinate \(r=r_h\), with
  one-sided spatial derivatives and no imposed incoming physical scalar
  characteristic from the excised interior.
- Outer boundary: approximate outgoing waves. The scalar condition is

\[
 (1-\beta)(\Pi+\Phi)
 +\frac{1-2\beta+\dot f}{J}\varphi\simeq0.
\]

For the geometric fields the leading large-J falloffs
\(a-1\sim J^{-1}\), \(B\sim J^{-2}\) guide outgoing conditions.
All finite-boundary conditions are approximate; boundary reflections must
be checked by increasing \(r_{\max}\).

## Numerical discretization

Uniform fixed grid \(r_i=r_h+i\Delta r\), with
\(J_i(t)=r_i+f(t)\), and Courant factor \(C=\Delta t/\Delta r\).
Interior finite differences are fourth-order:

\[
 (D_ru)_i=\frac{u_{i-2}-8u_{i-1}+8u_{i+1}-u_{i+2}}{12\Delta r}.
\]

The PDE system becomes ODEs by Method of Lines and is advanced with
classical RK4. The effective timestep is reduced automatically if the
computed radial scalar characteristic speed exceeds one in magnitude.

A sixth-order KO filter and Schwarzschild well-balancing are used as
numerical stabilizers. These are implementation choices and must be tested
under refinement.

## Waveforms and diagnostics

Fixed-areal-radius observers have moving numerical locations
\(r_{obs}(t)=J_{obs}-f(t)\). We sample at \(J_{obs}=10M_0,30M_0,50M_0\).
Also export the asymptotically normalized EF time \(T=t-f(t)\).

For weak fields, we expect \(\ell=0\) Schwarzschild scalar QNM ringing,
and ultimately a generic fixed-radius Price tail \(\varphi\sim T^{-3}\).
The GLMakie QNM and local-power estimators are exploratory until
convergence and boundary tests are completed.
