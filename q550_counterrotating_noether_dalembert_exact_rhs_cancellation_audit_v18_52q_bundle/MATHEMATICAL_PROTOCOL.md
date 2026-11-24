# Analytical derivation: exact angular transfer balance

Let `l_i=log(r_i)`, `v_i=dr_i/(r_i dt)`, `theta_i`, `omega_i=dtheta_i/dt`.
The archived Q550 angular equations (physical model time `t`) are

\[
\dot\omega_1=\kappa_0x\frac{r_2}{r_1}\sin\Delta-2v_1\omega_1,
\qquad
\dot\omega_2=\kappa_0x\frac{r_1}{r_2}\sin\Delta-2v_2\omega_2,
\quad \Delta=\theta_1+\theta_2-\phi_0-\Theta(t).
\]

With \(J_i=r_i^2\omega_i\), the product rule gives

\[
\dot J_1=2r_1^2v_1\omega_1+r_1^2\dot\omega_1
=\kappa_0 x r_1r_2\sin\Delta
=\dot J_2,
\qquad
\boxed{\dot C=0,\;C=J_1-J_2.}
\]

The original support coordinate is `s=Theta/(2*pi)` with positive factor
`dt/ds=2*pi/x` (where `0<x<1`). Multiplication of both identical torque
terms by the same positive clock factor preserves cancellation. Any smooth
monotone material-time reparameterization preserves it as well.

## Angular symmetry and generalized Noether balance

Under `theta_1 -> theta_1 + epsilon`, `theta_2 -> theta_2 - epsilon` the
argument `Delta` remains **unchanged**. The angular kinetic term
`T_angle=0.5*(r1²*omega1²+r2²*omega2²)` and interaction
`V_int=kappa0*x*r1*r2*cos(Delta)` are invariant.

Consequently `p_theta_i=dL/domega_i=r_i²*omega_i`, and the infinitesimal
generator `X=∂_{theta1}-∂_{theta2}` gives the Noether–d’Alembert balance
`d(p_theta1-p_theta2)/dt=Q_theta1-Q_theta2`. In the archived angular equations
there are **no extra angular generalized forces**, so the RHS is zero. Radial
damping and driven auxiliary `x` dynamics are present and prevent claiming a
purely conservative action for the entire 11D forced system unless all those
forces are modeled explicitly. A compatible Lagrange–d’Alembert representation
is obtained by putting their discrepancies in the radial and `x` generalized
forces, which annihilate `X`. Thus the `C` conservation is exactly symmetry-
protected in the documented forced model, not just a coordinate identity.

## Counterfactual breaker

Add an angular forcing `f1` to `omega1_dot` and `f2` to `omega2_dot`.
Then `dC/dt=r1²*f1-r2²*f2`. This will be nonzero unless the angular
forces are exactly balanced. A deliberate synthetic break must FAIL the charge
test, preventing tautological certification.

## Conditional noncollapse

Cauchy–Schwarz gives `C² <= (r1²+r2²) * (r1²*omega1² + r2²*omega2²)`;
writing `R²=r1²+r2²` and `T_ang=0.5*(r1²*omega1²+r2²*omega2²)`:

\[
T_{ang}\ge\frac{C^2}{2R^2}.
\]

If `C≠0` and `T_ang(t) <= K<∞` *for the entire future interval* then
`R(t)>=|C|/sqrt(2K)>0`. The observed finite-horizon upper bound is NOT a
guaranteed all-time energy bound. No individual radius or 5D relative
separation noncollapse is inferred. The continuous real charge is not an
integer winding number.
