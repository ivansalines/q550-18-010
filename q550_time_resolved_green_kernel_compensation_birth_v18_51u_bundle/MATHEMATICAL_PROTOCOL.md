# Q550 v18.51u — Mathematical protocol

## 1. Authenticated history and explicit chart

The verified v18.51t server sendback contains the exact physical v18.51s→r→q ancestry, analytic model-generated Jacobian samples, parameter-source samples, channel-forced states and independently reconstructed Green responses. The source domain is the same event-relative interval `xi ∈ [0,0.2]`, at 513 uniformly spaced observation nodes, for two original physical seed families at the historical parameter midpoint.

The explicit historical state ordering is

`Y=(logrho1,v1,theta1,omega1,logrho2,v2,theta2,omega2,logitx,vz)`.

Source splitting uses the **measured** vector

\[
 b(u)=\partial_{\kappa_0}f(Y_m(u),u;\kappa_m),\qquad
 b_{v_2}=e_5 b_5,\quad b_{\omega_2}=e_7 b_7,\quad
 b_{\mathrm{other}}=b-b_{v_2}-b_{\omega_2}.
\]

The third channel must not be dropped; it includes the first oscillator's parameter-source response. Source labels are coordinate-conditional and cannot be presumed independent experimental force actuators.

## 2. Two-time Green kernel

The original analytic state Jacobian samples `J(\xi)` are cubic-spline interpolated and the 10×10 fundamental solution is integrated on the 513-point observation grid:

\[
\partial_\xi \Phi(\xi,0)=J(\xi)\Phi(\xi,0),\qquad \Phi(0,0)=I.
\]

The propagation from an earlier source time `u` to an observation `xi` is

\[
\Phi(\xi,u)=\Phi(\xi,0)\Phi(u,0)^{-1}, \qquad 0\leq u\leq\xi.
\]

Projection is explicitly made at the **observation**, not the source. With the Cartesian immersion of the second oscillator, use

\[
L_r(\xi)=e_r(\theta_2(\xi))^{\mathsf T}D_Y Q_{34}(Y_m(\xi)),\quad
L_a(\xi)=e_a(\theta_2(\xi))^{\mathsf T}D_Y Q_{34}(Y_m(\xi)),
\]

where `Q34 = (rho2 cos theta2, rho2 sin theta2)`, `e_r=(cos theta2,sin theta2)` and `e_a=(-sin theta2,cos theta2)`. For this chart, `L_r=rho2 d(logrho2)` and `L_a=rho2 d(theta2)`, which are checked against the immersion rather than assumed from the start.

For source channel `j` and output `p ∈ {r,a}`,

\[
K_{p,j}(\xi,u)=
\begin{cases}
 L_p(\xi)\Phi(\xi,u)b_j(u), &u\leq\xi,\\
 0,&u>\xi,
\end{cases}
\qquad
S_{p,j}(\xi)=\int_0^\xi K_{p,j}(\xi,u)\,du.
\]

The full triangular arrays are kept at float64 precision; *no physical future-to-past response is represented*. The source value at `u=xi` has zero instantaneous positional feedthrough in this particular historical model because its direct κ source occupies velocity slots.

## 3. Independent tests and quadrature

- Compare `S_r,S_a` for all four channels and all 513 observations to the independently integrated forced sensitivities archived by v18.51t.
- Require `K_total=K_v2+K_omega2+K_other` and the analogous cumulative response recomposition.
- Independently restart a 10D tangent IVP at each of `u=0,0.05,0.10,0.15` to `xi=0.20`, initialized with the actual `b_j(u)`; compare its radial and azimuthal endpoint projections against the corresponding Green-kernel columns.
- Refine the historical source sampling interpolation between 129 and 257 nodes, while preserving the same physical reference midpoint trajectory.
- Compare the endpoint response under 257-source-node composite Simpson quadrature and 513-source-node composite Simpson quadrature.
- Monitor the condition number of `Phi(xi,0)` and abort if it becomes numerically unsuitable.

At even grid indices the cumulative integral is exact composite Simpson summation of the *sampled* kernel. Odd grid indices use composite Simpson up to `n-3` with a final Simpson 3/8 panel, except the first interval which uses trapezoid. This avoids conflating a last-panel first-order quadrature artifact with a physical phase effect.

## 4. Compensation birth is an operational threshold

At every observation, define the signed radial component sums

\[
R_j(\xi)=\int_0^\xi K_{r,j}(\xi,u)\,du,\qquad
C_r(\xi)=1-\frac{|\sum_j R_j(\xi)|}{\sum_j|R_j(\xi)|}.
\]

The `xi=0` ratio is undefined; mark it as missing. Record the **first sampled observation** where `C_r` reaches thresholds of 1%, 5% and 10%, with an additional small-amplitude guard for the first observed direct `v2`/`omega2` radial opposition. The marker depends on the thresholds, chart, grid and time horizon; it is **not** a singular dynamical creation event.

At the final observation, measure response-weighted source ages using positive Simpson weights applied to `|K(xi,u)|`: median source time, mean `xi-u`, and share of absolute kernel weight in the first half of source history. Positive and negative signed kernel masses are separately reported. An absolute-weighted mean lag is **not** a transit delay or memory lifetime.

## 5. Fixed source-window decomposition

At observation sections `xi=0.05,0.10,0.15,0.20`, integrate the source kernel separately across fixed consecutive windows of width `0.025`. Every window has 64 fine grid steps, so its composite Simpson integral reassembles the full integral without a panel mismatch. This answers whether a late-observed cancellation is dominated by early or recent source contributions without silently erasing a sign change within a channel.

## 6. Interpretational boundaries

These are numerical linear-response diagnostics along the midpoint orbit. The causal order `u≤xi` is a property of the initial-value Green operator and is **not** by itself a certificate of intrinsic material memory. No observer-independent radial/azimuthal causal decomposition, full 10D return, phase-space material holonomy or κ₀ shadowing for arbitrary long times is claimed. The verified horizon is `[0,0.2]` only.
