# Derivation map: ADM 3+1 to MMIEF scalar Einstein equations

The complete narrative LaTeX derivation exists separately as
`mmief_ekg_3plus1_derivation.tex` in the notes package. This condensed
reference shows each essential intermediate identity so an agent can check
signs and reproduce the implemented equations.

## 1. Spherical 3-geometry and ADM equations

$$\gamma_{ij}dx^idx^j=a^2dr^2+J^2d\Omega^2,\quad
n^\mu=\alpha^{-1}(1,-\beta,0,0),\quad
K^i{}_j=\operatorname{diag}(A,B,B).$$

The ADM evolution and constraint equations are

$$\partial_t\gamma_{ij}=-2\alpha K_{ij}+\mathcal L_\beta\gamma_{ij},$$

$$\begin{aligned}
(\partial_t-\mathcal L_\beta)K_{ij}
&=-D_iD_j\alpha+\alpha({}^{(3)}R_{ij}+KK_{ij}-2K_{ik}K^k{}_j)\\
&\quad-8\pi\alpha\left(S_{ij}-\frac12\gamma_{ij}(S-\rho)\right),\\
{}^{(3)}R+K^2-K_{ij}K^{ij}&=16\pi\rho,\\
D_j(K^j{}_i-\delta^j{}_i K)&=8\pi j_i.
\end{aligned}$$

For the massless scalar,

$$\rho=\frac{\Phi^2+\Pi^2}{2a^2},\quad
j_r=-\frac{\Phi\Pi}{a},\quad
S^r{}_r=\rho,\quad S^\theta{}_\theta=S^\phi{}_\phi
=\frac{\Pi^2-\Phi^2}{2a^2}.$$

**Important correction:** The angular stress is
`S^theta_theta=(Pi²-Phi²)/(2a²)`, *not* always `-rho`;
when reducing the angular ADM equation, the explicit angular stress term
cancels because `S^theta_theta - (S-rho)/2 = 0`.

## 2. Spatial Christoffel and Ricci terms

At fixed time the nonzero radial-angular connections include

$$\Gamma^r{}_{rr}=\frac{a'}a,\quad
\Gamma^r{}_{\theta\theta}=-\frac{JJ'}{a^2},\quad
\Gamma^\theta{}_{r\theta}=\frac{J'}J.$$

The scalar spatial curvature for general `J(t,r)` is

$$ {}^{(3)}R=\frac{2}{J^2}
\left(1-\frac{J'^2}{a^2}-\frac{2JJ''}{a^2}
+\frac{2JJ'a'}{a^3}\right). $$

Setting the MMIEF spatial gauge `J'=1, J''=0` yields

$$ {}^{(3)}R=\frac{2}{J^2}\left(1-\frac1{a^2}\right)
+\frac{4a'}{a^3J}. $$

The angular spatial Ricci component and Hessian are

$$ {}^{(3)}R^\theta{}_\theta
=\frac{1-1/a^2+Ja'/a^3}{J^2},\qquad
D^\theta D_\theta\alpha=\frac{\alpha'}{a^2J}. $$

## 3. Hamiltonian and momentum constraints

Use `K²-K_ij K^ij = (A+2B)²-(A²+2B²)=2B(2A+B)`.
Insert the curvature and `rho` into the Hamiltonian constraint, and
multiply by `a³ J/4`:

$$0=a'+\frac{a^3-a}{2J}
+\frac{a^3 J}{2}B(2A+B)-2\pi Ja(\Phi^2+\Pi^2).$$

For the radial momentum constraint,

$$D_j K^j{}_r-D_rK=2\frac{A-B}{J}-2B'=8\pi j_r,$$

so

$$ B'+\frac{B-A}{J}-\frac{4\pi\Phi\Pi}{a}=0,\qquad
A=B+JB'-\frac{4\pi J}{a}\Phi\Pi.$$

## 4. Evolution of a and B

The `rr` metric evolution gives

$$\partial_t(a^2)=-2\alpha a^2 A+2a^2\beta'+2aa'\beta,$$

hence

$$\dot a=-\alpha a A+(a\beta)'=-a^2(1-\beta)A+(a\beta)'.$$

The mixed angular extrinsic-curvature ADM evolution reduces to

$$ (\partial_t-\beta\partial_r)B
=-\frac{\alpha'}{a^2J}
+\alpha\left[\frac{1-1/a^2+Ja'/a^3}{J^2}
+(A+2B)B\right]. $$

Insert `alpha=a(1-beta)` and
`alpha'=(1-beta)a'-a beta'`; the two `a'` terms cancel:

$$\dot B=\beta B'+a(1-\beta)B(A+2B)
+\frac{1-\beta}{J^2}\left(a-\frac1a\right)
+\frac{\beta'}{aJ}.$$

The angular ADM matter term cancels **for a massless minimally coupled
scalar** by the precise stress identities above. Do not assume this for
arbitrary matter or a scalar potential.

## 5. Gauge and lapse/shift reconstruction

From the angular metric evolution,

$$\partial_t(J^2)=-2\alpha J^2 B+2\beta JJ'.$$

With `J'=1`, `Jdot=fdot`, one obtains

$$\dot f=\beta-\alpha JB.$$

Combined with `alpha=a(1-beta)`:

$$\beta=\frac{\dot f+aJB}{1+aJB},\qquad
\alpha=\frac{a(1-\dot f)}{1+aJB}.$$

Radial null rays have `dr/dt=-beta±alpha/a`:

$$ c_{in}=-1,\qquad c_{out}=1-2\beta. $$

## 6. Scalar evolution from covariant conservation law

The coordinate volume element is
`√(-g)=alpha*a*J²*sin(theta)`, and
`Phi=varphi'`, `Pi=(a/alpha)(dot(varphi)-beta Phi)` imply

$$\partial_t\varphi=\beta\Phi+\frac{\alpha}{a}\Pi
=\beta\Phi+(1-\beta)\Pi.$$

The wave equation `partial_mu(sqrt(-g) g^{mu nu} partial_nu varphi)=0`
in flux form becomes

$$\partial_t(J^2\Pi)
=\partial_r\left[J^2\left(\beta\Pi+\frac\alpha a\Phi\right)\right].$$

Since `dot J=dot f`:

$$\dot\Pi=\frac1{J^2}
\left[J^2(\beta\Pi+(1-\beta)\Phi)\right]'
-\frac{2\dot f}{J}\Pi.$$

Finally `dot Phi=(dot varphi)'`, which gives the stated `Phi` equation.

## 7. Trapped surface and horizon locking

For radial outward/inward null vectors `l_±=n±s`, with
`s^i=(1/a,0,0)`, the areal expansions are

$$\theta_\pm=\frac{2}{J}\left(\pm\frac1a-JB\right).$$

The **outer** expansion `theta_+=0` gives `aJB=1`. Set the horizon
coordinate `r_h=2M0` (fixed) and differentiate `aJB=1` along `r_h`.
Substitute the `a` and `B` evolution equations plus Hamiltonian and
momentum constraints; solve the resulting equation for the gauge motion.
The resulting (ingoing-branch) horizon-locking ODE is

$$\dot f=\left.\frac{4\pi J^2}{a^2}
 (\Phi+\Pi)^2\right|_{r_h}.$$

This is a *gauge condition obtained by preserving the horizon condition*,
not an independent Einstein field equation. It ensures the area radius of
the horizon changes as the black hole accretes positive ingoing scalar flux.

## 8. Mass aspect

Define the invariant Misner–Sharp mass

$$1-\frac{2m}{J}=\nabla_\mu J\nabla^\mu J.$$

In ADM coordinates and with `J'=1`,
`n^mu partial_mu J = (dot f-beta)/alpha=-JB`; hence

$$\nabla J\cdot\nabla J = \frac1{a^2}-J^2 B^2.$$

Therefore

$$m=\frac J2\left(1-\frac1{a^2}+J^2 B^2\right).$$

Differentiate radially, substitute the Hamiltonian constraint for `a'`
and the momentum constraint for `B'`, then simplify:

$$m'=4\pi J^2\left[
\frac{\Phi^2+\Pi^2}{2a^2}+\frac{JB}{a}\Phi\Pi\right].$$

In the program this ODE is integrated radially for the Hamiltonian
projection, after which `a` is recovered algebraically. The projection
changes the numerical scheme; it is not the same as solving unconstrained
ADM everywhere.

For source-aligned definitions see the repository's `EQUATIONS.md`,
`src/Geometry.jl`, `src/EinsteinEquations.jl`, `src/ScalarField.jl`,
`src/Constraints.jl`.
