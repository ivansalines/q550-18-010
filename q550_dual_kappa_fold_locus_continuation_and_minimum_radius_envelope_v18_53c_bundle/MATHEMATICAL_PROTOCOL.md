# Mathematical protocol — Q550 v18.53c

## Fixed-charge relative equilibrium

Work at the historical model parameters, prescribed *constant* signed external drive `s`, and fixed
`J = r1²ω1 − r2²ω2`. On the Δ=π sheet, `x=1/2` and `ω2=(r1²ω1−J)/r2²`.
The three stationary equations are the archived v18.53b residuals:

\[
 F_1=\omega_1^2-\frac{U'(r_1)}{r_1}+\kappa_0 x\frac{r_2}{r_1},
 \quad F_2=\omega_2^2-\frac{U'(r_2)}{r_2}+\kappa_0 x\frac{r_1}{r_2},
 \quad F_3=\omega_1+\omega_2-sx.
\]

\[
 U'(r)=2r-2.4r^2+0.8r^3-\frac{e_{\mathrm{reserve}}^2}{r^3},\quad e_{\mathrm{reserve}}=0.01.
\]

Let `a=(log r1,log r2,ω1)`. A simple fold for fixed `J` satisfies

\[
 F(a;J,\kappa_0,s)=0, \qquad \det F_a(a;J,\kappa_0,s)=0.
\]

The augmented Newton solver computes `F_a` analytically from the archived equations and normalizes its determinant by the product of the two largest singular values, improving near-null conditioning. The exact rank and root residual are checked independently.

For SVD left/right null vectors `w,v`, test the local conditions

\[
 w^{\mathsf T}F_s\neq0,\qquad
 w^{\mathsf T}F_{aa}[v,v]\neq0,
\]

with three curvature stencils. The reduced 8D Jacobian must have a real near-zero eigenmode and seven remaining modes with negative real part at the fold. Authentic 11D RHS and 8D Jacobian validation is mandatory for the source-certified status.

## Fold curve, implicit derivative and radius envelope

Write `z=(a,s)` and `G(z,κ)=(F,normalized-det(F_a))=0`. The fold locus has tangent

\[
 \frac{dz}{d\kappa}=-G_z^{-1}G_\kappa.
\]

In particular,

\[
 \frac{dr_2}{d\kappa}=r_2\frac{d(\log r_2)}{d\kappa},
\]

which is cross-checked against two independent centered finite differences of *freshly solved folds*. The 102-point narrow scan tests for sign changes and component switches. These are **sampled diagnostics, not a rigorous interval enclosure**. In particular, a continuous interior stationary point or disconnected branch cannot be mathematically excluded from a grid and local Newton solvers alone.

The two endpoint radii are evaluated at different fold drives:
`(κ_original,s_fold(κ_original))` and `(κ_autonomous,s_fold(κ_autonomous))`. The separate reference comparison at `s=-1` correctly holds the drive fixed for both κ prescriptions.

## Conservation, topology and physical limits

- This is `min(r1,r2)` **of a phase-locked critical equilibrium**; not `min_t r_i(t)` over physical orbits.
- A **fold in a stationary constant-drive parameter family** is not an observed bifurcation during a historical forward–return ramp.
- A near-zero eigenvalue and nondegenerate normal-form coefficients support a **local numerical saddle-node candidate**, not an interval-proof theorem.
- `J` is an exactly conserved continuous counterrotation charge of the documented angular sector; it does not imply conservation of total energy or quantized topological charge.
