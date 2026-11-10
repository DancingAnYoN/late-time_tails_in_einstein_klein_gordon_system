# Spherical 3+1 geometry and exact checks

General spherically symmetric ADM ansatz:

$$ds^2=-\alpha^2dt^2+a^2(dr+\beta dt)^2+J^2d\Omega^2,$$

where `J(t,r)` is the areal radius; only in an *areal gauge* is `J=r`.
The 3-metric is `gamma_ij=diag(a²,J²,J² sin²theta)`.
Let `A=K^r{}_r`, `B=K^theta{}_theta=K^phi{}_phi`, `K=A+2B`.

## Extrinsic curvature kinematics

The angular metric evolution yields

$$ \dot J=\beta J'-\alpha JB,\qquad
B=\frac{\beta J'-\dot J}{\alpha J}. $$

The radial metric evolution yields

$$\dot a=-\alpha aA+(a\beta)' .$$

The extrinsic contraction is

$$K^2-K_{ij}K^{ij}=2B(2A+B).$$

## Spatial Ricci quantities

For general `J(t,r)`:

$$ {}^{(3)}R=
\frac2{J^2}\left(1-\frac{J'^2}{a^2}
-\frac{2JJ''}{a^2}+\frac{2JJ'a'}{a^3}\right),$$

$$ {}^{(3)}R^r{}_r=-\frac{2J''}{a^2J}+
\frac{2J'a'}{a^3J}, $$

$$ {}^{(3)}R^\theta{}_\theta=
\frac{1-J'^2/a^2-JJ''/a^2+JJ'a'/a^3}{J^2}. $$

For `J'=1` and `J''=0`:

$$ {}^{(3)}R=\frac2{J^2}(1-1/a^2)+\frac{4a'}{a^3J}. $$

## Matter in the massless radial scalar reduction

Set `Phi=d_r varphi` and `Pi=a n^mu d_mu varphi`:

$$ \rho=\frac{\Phi^2+\Pi^2}{2a^2},\quad
j_r=-\frac{\Phi\Pi}{a},\quad
S^r{}_r=\rho,\quad
S^\theta{}_\theta=\frac{\Pi^2-\Phi^2}{2a^2}. $$

The correct angular stress is important when reducing ADM angular
extrinsic curvature evolution; the term
`S^theta_theta - (S-rho)/2` vanishes in the massless scalar case.

## Radial constraints for J'=1

$$0=a'+\frac{a^3-a}{2J}
+\frac{a^3 J}{2}B(2A+B)-2\pi Ja(\Phi^2+\Pi^2),$$

$$0=B'+\frac{B-A}{J}-\frac{4\pi\Phi\Pi}{a}.$$

For general `J'` the constraints contain additional `J'`/`J''` terms;
do not transplant these simplified formulas into a different gauge.

## Spherical mass and horizon

The invariant Misner–Sharp mass is

$$m=\frac J2(1-g^{\mu\nu}\partial_\mu J\partial_\nu J).$$

For general `J(t,r)` in ADM form,

$$g^{\mu\nu}\partial_\mu J\partial_\nu J
=\frac{J'^2}{a^2}-\frac{(\dot J-\beta J')^2}{\alpha^2}
=\frac{J'^2}{a^2}-J^2B^2.$$

Thus

$$m=\frac J2\left(1-\frac{J'^2}{a^2}+J^2B^2\right).$$

The outward/inward null expansions of symmetry spheres (for `J'>0`) are

$$\theta_+=\frac2J\left(\frac{J'}a-JB\right),\qquad
\theta_-=\frac2J\left(-\frac{J'}a-JB\right).$$

A marginally outer trapped sphere has `J'/a-JB=0`; in `J'=1`
this simplifies to `aJB=1`. At this surface, `2m=J`.

## Check cases

- Minkowski: `a=1`, `alpha=1`, `beta=0`, `J=r`, `A=B=0`, `m=0`.
- Schwarzschild in an ingoing, horizon-penetrating gauge:
  `a=sqrt(1+2M/r)`, `beta=2M/(r+2M)`, `alpha=1/a`, and
  `B=2M(r+2M)/[r(r+2M)]^(3/2)` with `J=r`.
- Derive formulas anew for other gauge conditions, negative `J'`,
  scalar potentials, higher dimensions or exotic stress tensors.
