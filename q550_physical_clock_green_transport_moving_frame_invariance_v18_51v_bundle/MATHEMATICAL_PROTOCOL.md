# Mathematical protocol — Q550 v18.51v

## 1. Causal source-response transport in xi

The authenticated v18.51u audit produces the model-generated first variation at the historical midpoint kappa:

`S_j(xi) = integral_0^xi Phi_xi(xi,u) b_j(u) du`.

For the natural second-plane radial/azimuthal observable `L(xi)`, the two-time density is

`K_xi(xi,u) = L(xi) Phi_xi(xi,u) b_j(u)`.

Only `u <= xi` is causally admissible. The source split is `v2`, `omega2`, `other`, and `total`; the fourth channel is the exact signed sum of the first three.

## 2. Physical clock as a trajectory-dependent map

Let `tau(xi; kappa)` be the already integrated physical clock and `c(xi;kappa)=d tau/d xi >0` the clock rate, both saved by v18.51r. Along one **fixed reference history**,

`d tau_source = c(u) du` and `K_tau(tau(xi),tau(u)) = K_xi(xi,u)/c(u)`.

Thus, for a frozen fixed-xi tangent equation,

`integral_0^tau K_tau(tau_obs,tau_source) d tau_source = integral_0^xi K_xi(xi,u) du`.

The positivity of `c`, trajectory-integrated tau and independent nonuniform-tau composite Simpson quadrature are checked. Note that this is **not** the statement that a parameter derivative at fixed tau is unchanged.

## 3. The required clock correction for distinct kappa trajectories

For any physical Cartesian configuration `Q(xi;kappa)` and positive clock `tau(xi;kappa)`, the chain rule gives

`(partial Q / partial kappa)|tau = (partial Q / partial kappa)|xi - (dQ/dxi)/c * (partial tau/partial kappa)|xi`.

This is the crucial difference between a relabeled frozen-clock kernel and actual parameter sensitivity measured at equal physical time. The historical original and autonomous kappas are distinct. The audit uses their nonlinear v18.51r orbits to compute central secants and compares the chain-rule prediction against **independent cubic interpolation of both trajectories on the same tau grid**. The midpoint v18.51t tangent sensitivity provides a second comparison. The extra clock-shift contribution is kept as a **separate response channel**, rather than being spuriously assigned to `v2` or `omega2`.

Secants and cubic interpolation have finite errors, especially when the kappa-dependent clock difference is very small. No interval error certificate is claimed.

## 4. Passive moving physical-plane frames

At observation xi, let `E_theta(xi)` be the plane rotation from the natural radial/azimuthal frame into fixed Cartesian lab coordinates. For any declared smooth observer rotation `psi(xi)`, write

`F(xi)=E_theta(xi) R(psi(xi))`,

`K_psi(xi,u)=R(-psi(xi)) K_radial(xi,u)` and `I_psi(xi)=R(-psi(xi)) I_radial(xi)`.

Then `F I_psi=E_theta I_radial` and `F K_psi=E_theta K_radial`; 2D Euclidean vector norms and vector-valued channel recomposition are unchanged. Individual components need not be invariant, and neither does a signed component-based cancellation fraction.

For `w(xi)` expressed in the rotating frame with angular coordinate `alpha=theta2+psi`, its physical (Cartesian) derivative is

`d v_lab/d xi = R(alpha) [w' + alpha' J w]`,

with `J=[[0,-1],[1,0]]`. The connection contribution is **essential**; dropping it produces a false derivative discrepancy under observer rotation. Two unrelated observer schedules are tested:

- `psi_1(xi)=0.38 sin(2 pi z) +0.27 z`
- `psi_2(xi)=0.31 cos(3 pi z) +0.17 z +0.08`

where `z=(xi-xi_0)/(xi_end-xi_0)`.

A smooth `O(2)` passive change is a test of representational consistency, **not** evidence for invariance under arbitrary active physical perturbations or under different state-space metrics.

## 5. Response-weighted source ages

The origin and observation physical times are reported for the same events. Response-weighted source-age quantiles are transported from xi to tau using the authenticated positive clock. Source ages are not intrinsic relaxation times. Signed and absolute channel contributions remain tied to the chosen observation and source decomposition.

## 6. Quality controls

- SHA256-matched nested genealogy U -> T -> S -> R; authenticated cached numeric arrays.
- A fast analytic synthetic test for positive time reparameterization, causal support, moving-frame connection and a parameter-dependent physical clock.
- `d tau/dxi` vs the archived positive physical `clock`.
- Reparameterized two-time kernel density and independent nonuniform-tau quadrature against U's reference integrals.
- Four-source-channel additivity, causality and observer-frame Cartesian reconstruction.
- Connection-covariant derivative reconstruction under two separate passive frame schedules.
- Fixed-tau matched nonlinear kappa secant vs the clock-corrected fixed-xi sensitivity.
- Independent v18.51t tangent sensitivity against the distinct-kappa nonlinear secant.
- Strictly numerical classification; no interval proof or material memory claim.
