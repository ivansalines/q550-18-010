# v18.52v-r1 — Mathematical protocol (English)

## 1. Source-certified sign

The actual historical source defines `s(xi)=cubic_clock_sign(xi)` and `Phi(xi)=cubic_phase_excursion(xi)` with `Phi'=-s` for xi>0 (piecewise regular parts, excluding boundary nondifferentiability). The original *return* phase is `Theta=Theta0-omega*(Phi-Phi0)`; therefore

`Theta_xi_prime=+omega*s`.

The original signed RHS computes `F[10]=s*(2*pi*slope)`; `omega=2*pi*slope` is independently checked. The forward protocol uses `Theta_xi_prime=+omega`. Since the native material time obeys `tau_xi_prime=omega/x_native`, along each native state (including kicked trajectories):

`Theta_tau_prime=s*x_native`.

The external Theta is shared as a function of xi, not as an independently-integrated oscillator per kicked orbit. Differentiation with respect to each orbit's own tau therefore carries its own x-native factor. A generic passive `Theta -> -Theta` transformation is not asserted.

## 2. Native fixed model parameter

The historical seed kappa values are `k_original=1/6`, `k_autonomous≈0.16659961041764054`, but the **native 52l trajectories** evolved under the shared `k_fixed=(k_original+k_autonomous)/2≈0.1666331385421536`. The *physical energy* must use this `k_fixed`, never the individual historical seed kappa. This is an authenticated genealogy invariant, not a fitted parameter.

## 3. Energy and work

`Delta = theta1+theta2-phi0-Theta`.

`E_tau=-gamma_rho*(rdot1²+rdot2²)-gamma_x*xdot²+r1*r2*(mu_x*sinDelta+k_fixed*cosDelta)*xdot + k_fixed*x*r1*r2*sinDelta*Theta_tau`.

`W_drive=integral(k_fixed*x*r1*r2*sinDelta * Theta_tau, d tau)`.

Frozen-state sign-switched comparison: `Theta_tau -> -Theta_tau` at the **same recorded state** reverses instantaneous drive power. It is not a new physically integrated forward branch and does not establish forward-return hysteresis work.

The time-varying frame `H=E - (c*x/2)*(J1+J2)` with `c=Theta_tau/x=s` has source power `H_tau= -D_rho-D_x+W_x - (x*c_tau+c*x_tau)*Jsum/2 + k*x*r1*r2*sinDelta*(Theta_tau-c*x)`. The last factor vanishes for source-consistent c. Include c_tau if reversal ramp is included.

The counterrotational first integral `Jdiff=J1-J2` remains an algebraic identity for source-consistent sector angular equations regardless of external phase profile. Neither global uniform energy bounds nor physical separation noncollapse are concluded.

## 4. Quadrature convergence

For each native series, evaluate at `stride=[1,2,4,8,16,32]` (inherited sampling). Compare maximum absolute integrated energy residual and numerical phase-rate estimate to the exact source-driven phase rate. The intentionally wrong seed-specific kappa check is presented as a falsifier and must not be used as a calibration.
