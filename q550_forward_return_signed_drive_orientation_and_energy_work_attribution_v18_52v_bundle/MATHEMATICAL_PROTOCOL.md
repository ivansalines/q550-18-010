# Mathematical protocol — Q550 v18.52v

## Two physically distinct operations
- Source-defined forward: `Theta_xi=+omega`, `clock_sign=+1`.
- Source-defined return: `Theta_xi=-omega*cubic_clock_sign(xi,design)`, while the original source passes `cubic_clock_sign(xi,design)` into the physical RHS. The actual source is interrogated without assuming that the argument is the phase derivative's own sign.
- A **passive coordinate relabel** `Theta_new=-Theta_old` also changes the coupling `Delta=theta1+theta2-phi0-Theta_old`; it is not the return protocol merely by notation.

## Energy work identity (material time)
Define `Delta=theta1+theta2-phi0-Theta`, `J_i=r_i^2 omega_i`, `x=sigmoid(z)` and `kappa_eff=kappa_0*x`.

`dE/dtau = -gamma_rho*(r1dot^2+r2dot^2) - gamma_x*xdot^2 + r1*r2*(mu_x*sin(Delta)+kappa_0*cos(Delta))*xdot + kappa_0*x*r1*r2*sin(Delta)*Theta_dot`.

At the **same instantaneous 10D state and Theta**, changing only `Theta_dot -> -Theta_dot` flips the last term exactly. For trajectories evolved under different signs, their states and Delta will generally differ. A sign flip therefore does NOT imply equal and opposite integrated work in an actual forward/return experiment.

For `c(tau)=Theta_dot/x`, the clock-adapted frame `H_c=E - c*x*(J1+J2)/2` has

`H_c'=-D_rho-D_x+W_x - [c_dot*x+c*x_dot]*(J1+J2)/2`

provided `Theta_dot=c*x` and `Jsum_dot=2*kappa_0*x*r1*r2*sin(Delta)`.

The `c_dot` term is essential during a smooth reversal ramp. Neither this identity nor the return ledger alone gives a uniform all-time bound.

## Source clock and material time
In the historical source `xi` is support time, `d tau/d xi = (2*pi*engine.slope)/x`. Thus `Theta_dot_tau = (x/(2*pi*slope))*Theta_prime_xi`. The original return source has `Theta_prime_xi=-engine.omega*cubic_clock_sign`, where `engine.omega=2*pi*slope`. Therefore `Theta_dot_tau=-x*cubic_clock_sign` when both are evaluated for the same trajectory. For kicks driven by the *shared* reference clock, use the actual saved tau for every kick, not `x_kick` as a surrogate.

## No non-existent forward trajectory
The v18.52l native data contain an extended **return** segment. v18.52v quantifies the authentic return work and matched-state sign-control sensitivity. It does not fabricate a second forward trajectory, and does not prove a physical work loop or true time reversal.
