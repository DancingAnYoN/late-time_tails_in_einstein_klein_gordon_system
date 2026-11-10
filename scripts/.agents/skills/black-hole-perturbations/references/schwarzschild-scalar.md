# Massless scalar perturbations on Schwarzschild

Schwarzschild metric in usual Killing time `t_S`:

$$ds^2=-F(r)dt_S^2+F(r)^{-1}dr^2+r^2d\Omega^2,\qquad
F(r)=1-\frac{2M}{r}.$$

For a scalar `square_g varphi=0`, decompose

$$\varphi(t_S,r,\theta,\phi)=\sum_{lm}
\frac{\psi_{lm}(t_S,r)}{r}Y_{lm}(\theta,\phi).$$

Define tortoise coordinate by `dr_*/dr=1/F`:

$$r_*=r+2M\ln\left(\frac{r}{2M}-1\right).$$

Then each multipole obeys

$$\left[-\partial_{t_S}^2+\partial_{r_*}^2-V_l(r)\right]
\psi_{lm}=0,$$

$$V_l(r)=F(r)\left[\frac{l(l+1)}{r^2}+\frac{2M}{r^3}\right].$$

For `l=0`:

$$V_0(r)=\frac{2M}{r^3}\left(1-\frac{2M}{r}\right).$$

With `exp(-i omega t_S)`, QNM boundary conditions are ingoing at the
future horizon, `psi~exp(-i omega r_*)` for `r_*->-infinity`, and
outgoing at spatial infinity, `psi~exp(+i omega r_*)` for
`r_*->+infinity`. Both are understood by analytic continuation for
complex `omega`; naive growing spatial exponentials are not normalizable
modes in the ordinary bound-state sense.

## Time coordinate caveat: MMIEF

The **stationary** MMIEF metric used by the Julia project has
`beta=2M/(r+2M)`, `alpha/a=1-beta` and `f=0`.
Outgoing radial null rays satisfy

$$ \frac{dr}{dT}=1-2\beta=\frac{r-2M}{r+2M}, $$

so

$$ \Delta T=(r_2-r_1)
+4M\ln\frac{r_2-2M}{r_1-2M}. $$

This is **not** the tortoise-coordinate travel time, which has coefficient
`2M` in the logarithmic term. One cannot transfer a fitting-window shift
between time coordinates without doing the coordinate transformation.
In a genuinely dynamical EKG metric the stationary expression is only
an approximation.

## Scalar vs gravitational modes

The reference scalar `l=0` QNM is approximately

$$M\omega_{00}^{(s=0)}=0.110455-0.104896\,i.$$

A commonly quoted gravitational axial/polar `l=2` fundamental mode is

$$M\omega_{20}^{(s=2)}\simeq0.373672-0.088962\,i.$$

Do not substitute one for the other. Metric perturbations in spherical
symmetry have no independent radiative `l>=2` gravitational degrees of
freedom; coupled spherical scalar dynamics can nevertheless display
spherical scalar QNM-like ringing in the weak-field regime.

Background reference: Chandrasekhar, *The Mathematical Theory of Black
Holes*; Berti–Cardoso–Starinets review https://arxiv.org/abs/0905.2975.
