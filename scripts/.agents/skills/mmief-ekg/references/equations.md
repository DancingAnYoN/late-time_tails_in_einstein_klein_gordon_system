# MMIEF–EKG equations (project reference)

Definitions: signature `(-,+,+,+)`, `G=c=1`, `K_ij=-1/2 L_n gamma_ij`,
primes `d/dr`, dots `d/dt`.

## Field equations

$$G_{\mu\nu}=8\pi T_{\mu\nu},\qquad \square_g\varphi=0,$$

$$T_{\mu\nu}=\partial_\mu\varphi\partial_\nu\varphi
-\tfrac12g_{\mu\nu}(\nabla\varphi)^2.$$

Metric, gauge and reduction:

$$ds^2=-\alpha^2dt^2+a^2(dr+\beta dt)^2+J^2d\Omega^2,$$

$$J=r+f(t),\quad \alpha=a(1-\beta),\quad
A=K^r{}_r,\quad B=K^\theta{}_\theta,\quad K=A+2B,$$

$$\Phi=\partial_r\varphi,\qquad
\Pi=\frac a\alpha(\partial_t\varphi-\beta\Phi).$$

## Evolution (field state a, B, Phi, Pi, varphi, f)

$$\begin{aligned}
\dot a&=-a^2(1-\beta)A+(a\beta)',\\
\dot B&=\beta B'+a(1-\beta)B(A+2B)
 +\frac{1-\beta}{J^2}\left(a-\frac1a\right)+\frac{\beta'}{aJ},\\
\dot\Phi&=[\beta\Phi+(1-\beta)\Pi]',\\
\dot\Pi&=\frac{1}{J^2}[J^2(\beta\Pi+(1-\beta)\Phi)]'
-\frac{2\dot f}{J}\Pi,\\
\dot\varphi&=\beta\Phi+(1-\beta)\Pi,\\
\dot f&=\left.\frac{4\pi J^2}{a^2}(\Phi+\Pi)^2\right|_{r=r_h}.
\end{aligned}$$

## Constraints

$$\begin{aligned}
0&=a'+\frac{a^3-a}{2J}+\frac{a^3J}{2}B(2A+B)
-2\pi Ja(\Phi^2+\Pi^2),\\
0&=B'+\frac{B-A}{J}-\frac{4\pi}{a}\Phi\Pi,\\
0&=\Phi-\partial_r\varphi.
\end{aligned}$$

Reconstruction:

$$A=B+JB'-\frac{4\pi J}{a}\Phi\Pi,\quad
\beta=\frac{\dot f+aJB}{1+aJB},\quad
\alpha=\frac{a(1-\dot f)}{1+aJB}.$$

## Horizon and quasi-local mass

$$\theta_+=\frac2J\left(\frac1a-JB\right),\quad
\theta_-=-\frac2J\left(\frac1a+JB\right),$$

$$a_hJ_hB_h=1,\quad r_h=2M_0,\quad J_h=2M_0+f(t),\quad
M_h=\frac{J_h}{2}.$$

$$m=\frac J2\left(1-\frac1{a^2}+J^2B^2\right),$$

$$m'=4\pi J^2\left[\frac{\Phi^2+\Pi^2}{2a^2}
+\frac{JB}{a}\Phi\Pi\right],\quad
a=\left[1-\frac{2m}{J}+J^2B^2\right]^{-1/2}.$$

## Initial and boundary prescriptions

Scalar initial pulse:

$$\varphi_0=\epsilon\chi(r)e^{-(r-r_0)^2/\sigma^2},\quad
\Phi_0=\varphi_0',\quad
\Pi_0=\Phi_0+\frac{\varphi_0}{r(1-\beta_0)}.$$

Schwarzschild stationary test (at `f=0`):

$$a_S=\sqrt{1+2M/r},\quad\beta_S=\frac{2M}{r+2M},$$

$$B_S=\frac{2M(r+2M)}{[r(r+2M)]^{3/2}},\quad
A_S=-\frac{2M(r+M)}{[r(r+2M)]^{3/2}}.$$

Approximate scalar outgoing boundary:

$$(1-\beta)(\Pi+\Phi)
+\frac{1-2\beta+\dot f}{J}\varphi\simeq0.$$

Do not treat the finite-radius outer boundary as perfectly transparent.
The code solves the constraints on the initial slice; with nonzero scalar
amplitude the geometry is not exactly stationary Schwarzschild.

## Domain and numerical observers

The grid is uniform in `r` (not `J`); observer interpolation uses
`r_obs(t)=J_obs-f(t)` and exported time `T=t-f(t)`.
Physical waves are analyzed at `J_obs=(10,30,50)M0` in standard runs.
