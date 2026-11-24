# Q550 v18.52w — Mathematical protocol

## Historical drive-zero condition

Source `Design('v18_50p','beta:2:2',0.24)` defines a signed cubic clock `s(t)=cubic_clock_sign(t,des)`, with `s(0.12)=0` and opposite signs immediately before/after. The return phase is `Theta(t)=Theta_0-omega*(cubic_phase_excursion(t,des)-Phi_0)` and the authenticated cubic source satisfies `Phi'(t)=-s(t)`, hence `Theta'(t)=omega s(t)` in the source support coordinate `t=event+xi`.

We locate the zero using the original source callable, not by fitting a polynomial to the archived phase. Original and autonomous seed contexts have slightly different `event`, so `xi_zero=0.12-event(seed)` is seed-dependent. The v18.51r midpoint trajectories must bracket both events.

## Genuine material motion

The original signed RHS receives `F_11(Y,Theta,kappa,s)`. We record `F_11[10]=omega s` and the independent material tangent `F_11[:10]`. The Cartesian 5D embedding is

`Q_5=(r1 cos(theta1), r1 sin(theta1), r2 cos(theta2), r2 sin(theta2), x)`, `r_i=exp(log r_i)`, `x=sigmoid(z)`, with `Q'_5=Jacobian(Q_5)*F_11[:10]`.

We check `||Q'_5||` at `s=0`, the direction angles at +/-0.005, 0.010, 0.020 support units, and `||F_11[:10]||`. A material full 10D tangent norm includes coordinates of mixed physical units and must not be interpreted as a calibrated energy magnitude. The geometrical 5D tangent is an independent complementary measurement. The old v18.51r tabulated `Y` is interpolated with PCHIP; its derivative is compared with the authentic source RHS, and the mismatch is reported.

**Work balance:** along original return trajectories, `dE/dtau=-D_rho-D_x+W_x+W_Theta`, `W_Theta=kappa*x*r1*r2*sin(Delta)*Theta_dot_tau`, and `Theta_dot_tau=x*s(t)` with `dtau/dxi=omega/x`. At `s=0`, the phase work power vanishes, but dissipation, x-source work and material velocities need not vanish. Integrals across small symmetric windows compare the measured energy change against source work and dissipation using the actual historical states. No smooth extrapolation creates a genuine forward trajectory.

**Noether-d'Alembert angular difference:** `J1=r1^2*omega1`, `J2=r2^2*omega2`, `Jdiff=J1-J2`. Conservation of this angular-sector quantity is independent of drive orientation but is NOT full energy conservation.

## One external-phase-turn recurrence scout

The v18.52l long return states occupy support `[12.8,25.6]`. From the historical phase function find pairs `(xi_a,xi_b)` for `Theta(xi_b)=Theta(xi_a)-2pi`; the source phase is shared across perturbations. Compare the entire physical 10D state modulo the action `(theta1,theta2)->(theta1+epsilon,theta2-epsilon)`. A complete reduced coordinate set is

`(log r1, v1, omega1, log r2, v2, omega2, z, vz, [theta1+theta2] mod 2pi)`.

For each pair calculate nonangular differences and `wrap((theta1_b+theta2_b)-(theta1_a+theta2_a))` instead of fitting an alignment angle. The drive phase is matched modulo `2pi` explicitly. The descriptive scaled residual divides mixed-unit norm by `max(1, ||nonangular_reference||)` and is NOT a physically justified metric. The raw componentwise residuals, gauge shift and charge drift are also exported.

**Certification hierarchy:** drive-zero clock rate zero vs material tangent nonzero (local statement); one external phase turn near recurrence (finite-window numerical statement); repeated full section return modulo symmetry (requires more observed periods); attracting limit cycle (requires stability/Floquet or contracting Poincare return and autonomous/periodic forcing); conservative full flow (requires independent energy/structure proof, not Jdiff alone). No stage is promoted automatically to a later stage.
