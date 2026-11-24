# Q550 v18.52t — mathematical protocol

## Question and preserved model

Determine the actual orientation of the external phase clock used to evaluate the
v18.52l 10D states, separately from the native 11D structural generator. Test the
ordinary energy budget and rotating-frame energies under both possible signs.
Do not adjust the clock sign to minimize the ledger residual; measure the phase
orientation independently first.

The native v18.52l array consists of 10 physical dynamical variables plus an
**elapsed material-time column**. It is NOT an 11D original model state; in
particular column 11 is not `Theta`. External `Theta` is authenticated from the
v18.51s source and `ctx['theta'](ctx['event']+xi)` per seed. A single externally
prescribed phase is shared by all kick magnitudes of a seed. Thus the reference
may obey a precise signed clock relation while kicks have a small mismatch of
`x_kick` relative to `x_reference`.

Set `r_i=exp(y_i)`, `x=sigmoid(z)`, `Delta=theta1+theta2-phi0-Theta`,
`k_eff=kappa0*x`, and `J_i=r_i^2*omega_i`. Use the seven SHA-authenticated
Q550 model parameters exactly as in `inputs_config.json`.

The ordinary physical energy is

`E = 1/2 sum_i (rdot_i^2 + r_i^2 omega_i^2) + m_x xdot^2/2
    + sum_i U(r_i) + V_B(x) + kappa0*x*r1*r2*cos(Delta)`.

For **arbitrary prescribed external phase velocity** `Theta_dot=w(t)`, the
mechanical ten equations yield the exact identity

`E_dot = -D_rho - D_x + W_x + kappa0*x*r1*r2*sin(Delta)*w`,

where

`D_rho=gamma_rho*(rdot1^2+rdot2^2)`, `D_x=gamma_x*xdot^2`, and
`W_x=r1*r2*(mu_x*sin(Delta)+kappa0*cos(Delta))*xdot`.

The angular equations additionally give

`J_sum_dot=2*kappa0*x*r1*r2*sin(Delta)`, `J_diff_dot=0`.

For `sigma in {-1,+1}`, define

`H_sigma=E - sigma*x*(J1+J2)/2`.

Exactly when the phase clock satisfies `Theta_dot=sigma*x`,

`H_sigma_dot= -D_rho-D_x+W_x - sigma*xdot*(J1+J2)/2`.

If the prescribed clock is **not** `sigma*x`, the right-hand side acquires a
residual work term

`kappa0*x*r1*r2*sin(Delta)*(Theta_dot-sigma*x)`.

No equality `H_+=H_-` is presumed. The initial phase anchors must remain
unchanged when comparing signs: a post hoc phase change is not a signed-clock
comparison. A physical orientation-reversal transformation would require
transforming the entire dynamical system and angular variables consistently.

## Measurements

1. Authenticate `v18.51s` and the 10D v18.52l archive by exact SHA256.
2. Reconstruct external phase per seed, differentiate on the native material-time
   clock, and compare `Theta_delta - sigma*integral(x dt)` for both signs.
3. Probe the original `structural_jet` off orbit. Normalize its support clock by
   the documented `x/(2*pi*engine.slope)` factor and compare **all 11** components
   to the historical 11D generator, including the source clock component.
4. Evaluate energy and both signed frame balances on each of the 20 authentic
   series without integrating a new trajectory.
5. Include maximum cumulative energy mismatch and integrated absolute work
   mismatch, not merely endpoint closure.
6. State explicitly if the 11th generator component and externally applied clock
   obey opposite orientations. Conservation in the ten mechanical coordinates
   is not a certificate of full eleven-dimensional clock equivalence.
7. Separate finite-horizon energy budgets from a uniform global bound. No
   topological winding or 5D relative noncollapse claim follows.

## Algebraic verification

`scripts/symbolic_signed.py` performs a SymPy exact check when available and
`scripts/phase_energy.py` plus `tests/test_synthetic.py` independently test the
full arbitrary-phase-rate identity off orbit. Neither substitutes for the
SHA-authenticated generator checks.
