# ADM 3+1 reference: conventions, evolution and sources

Adopt signature `(-,+,+,+)` and
`K_ij=-(1/2) L_n gamma_ij`. A spacetime foliation `Sigma_t` has

$$ds^2=-\alpha^2dt^2+\gamma_{ij}(dx^i+\beta^i dt)
(dx^j+\beta^j dt),$$

$$n_\mu=(-\alpha,0,0,0),\qquad
n^\mu=\alpha^{-1}(1,-\beta^i),\qquad
\gamma_{\mu\nu}=g_{\mu\nu}+n_\mu n_\nu.$$

Projection of stress-energy:

$$\rho=T_{\mu\nu}n^\mu n^\nu,\quad
j_i=-\gamma_i{}^\mu n^\nu T_{\mu\nu},\quad
S_{ij}=\gamma_i{}^\mu\gamma_j{}^\nu T_{\mu\nu},\quad
S=\gamma^{ij}S_{ij}.$$

ADM constraint equations:

$$ {}^{(3)}R+K^2-K_{ij}K^{ij}=16\pi\rho, $$

$$ D_j(K^j{}_i-\delta^j{}_i K)=8\pi j_i. $$

ADM evolution equations:

$$ \partial_t\gamma_{ij}=-2\alpha K_{ij}+\mathcal L_\beta\gamma_{ij}, $$

$$\begin{aligned}
(\partial_t-\mathcal L_\beta)K_{ij}
={}&-D_iD_j\alpha+\alpha\left({}^{(3)}R_{ij}
+K K_{ij}-2K_{ik}K^k{}_j\right)\\
&-8\pi\alpha\left[S_{ij}-\tfrac12\gamma_{ij}(S-\rho)\right].
\end{aligned}$$

For a real massless scalar,

$$T_{\mu\nu}=\nabla_\mu\varphi\nabla_\nu\varphi
-\frac12g_{\mu\nu}(\nabla\varphi)^2.$$

Define `Pi_n=n^mu d_mu varphi` and `Phi_i=D_i varphi`. Then

$$\rho=\frac12(\Pi_n^2+\Phi_i\Phi^i),\quad
j_i=-\Pi_n\Phi_i,\quad
S_{ij}=\Phi_i\Phi_j
+\frac12\gamma_{ij}(\Pi_n^2-\Phi_k\Phi^k).$$

The trace and scalar pressure depend on spatial gradients. If instead
`Pi=a Pi_n` in spherical radial gauge, the factors of `a` change;
see the spherical reference. A scalar potential `V(varphi)` adds
`+V` to `rho`, `-V gamma_ij` to `S_ij`, and a source to `square varphi`;
never reuse massless cancellation identities without recalculating.
