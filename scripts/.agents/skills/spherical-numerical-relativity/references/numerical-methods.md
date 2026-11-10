# Method of Lines in spherical relativity

## Workflow

1. Choose primary variables and gauge; express `partial_t U = F(U, d_rU, ...)`.
2. Choose a fixed grid coordinate (`r`) and distinguish it from physical
   areal radius (`J(t,r)`).
3. Identify characteristic speeds and whether boundaries are inflow/outflow.
4. Implement centered/biased finite differences and explicit boundary closures.
5. Convert the spatial PDE to a finite ODE system (`dY/dt=F_h(Y)`).
6. Evolve with RK4 while applying declared constraint/boundary operations.
7. Inspect mass, constraints, horizon behavior and continuum convergence.

Interior fourth-order first derivative:

$$ (D_r\varphi)_i=\frac{\varphi_{i-2}-8\varphi_{i-1}
+8\varphi_{i+1}-\varphi_{i+2}}{12\Delta r}. $$

Classical RK4 for the *coupled* state vector `Y`:

$$\begin{aligned}
k_1&=F_h(Y_n),\\
k_2&=F_h(Y_n+\tfrac12\Delta t k_1),\\
k_3&=F_h(Y_n+\tfrac12\Delta t k_2),\\
k_4&=F_h(Y_n+\Delta t k_3),\\
Y_{n+1}&=Y_n+\tfrac{\Delta t}{6}(k_1+2k_2+2k_3+k_4).
\end{aligned}$$

Courant factor `C=Delta t/Delta r`; for unit wave speed and a given
stencil/RK4, confirm stability empirically and via linear stability analysis.
Never assume `C<1` alone proves stability of every discrete operator.

## Projection/constrained evolution

A constraint projection changes the time-discrete map:
`Y_{n+1} = Projection(RK4(Y_n))`, possibly with projection at each
intermediate stage. Do not describe this as pure unconstrained evolution.
Keep this projection distinct from the generic integrator API. Note that
projection and boundary closures can reduce the global convergence order.

## Spherical numerics considerations

- Use parity regularity at `r=0` if simulating regular-center problems.
  Do **not** apply regular-center formulas at an excised black-hole horizon.
- Beware `1/J` and `1/J²` coefficients; analytic regularization may be needed.
- In horizon-penetrating coordinates, count incoming characteristics at each
  computational boundary, rather than imposing all fields arbitrarily.
- Outer Sommerfeld conditions are approximate for curved backgrounds;
  use `rmax` variation and causal arrival-time estimates.
- Dissipation: specify operator, coefficient, location and refinement scaling.
- Long-time tails can fall below the spatial truncation error long before the
  evolution fails a basic stability check.
- For a scalar wave on a fixed potential background,
  `partial_t varphi=Pi`, `partial_t Pi=D_xx varphi-V(x,t)varphi`
  is a useful independent MoL/RK4 test.
