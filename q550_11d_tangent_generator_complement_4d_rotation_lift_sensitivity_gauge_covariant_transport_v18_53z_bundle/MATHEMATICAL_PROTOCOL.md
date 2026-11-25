# Q550 v18.53z — mathematical protocol

## Authenticated model

The original Q550 11D trace is `z(xi)`; its last variable is the material time accumulator. The first ten variables are log-radius, angular and logit coordinates and velocities for the two coupled oscillators and the third degree of freedom. In the historical xi clock,

\[ dz_{0:10}/d\xi = (\omega/x) f_{\tau,10}(z_{0:10},\phi_{ext}(\xi),\kappa),\quad d\tau/d\xi=\omega/x,\quad x=\operatorname{expit}(z_8). \]

The function `f_{tau,10}` is the **exact historical material RHS** implemented in the SHA-pinned `torque_math_v18_53f_authenticated.py`. The phase is inferred by solving the model's linear `cos Delta, sin Delta` decomposition against the **authenticated** Frobenius-isometric 10D source. The overdetermined 10D source residual and the `cos²+sin²=1` defect are mandatory checks. This is source-conditioned phase recovery and is not an independent external-phase observation.

At fixed external phase, define the **11D tangent generator**

\[A_{11}(\xi) = D_z f_\xi(z(\xi),\xi).\]

Its final material-time column is zero. Central-difference Jacobians at `h` and `h/2` are compared in real Q550 states; a synthetic test uses known model terms. The derivative of the archived 11D cubic-spline trace provides an additional diagnostic, **not** a second physical solution.

## Nonunique 5D horizontal lifts

Let `q(z)` be the authenticated mass-normalized physical Cartesian embedding and `J_q=D_z q` (5×11, rank 5). A right inverse `L` obeys `J_q L=I_5`. There is no unique lift because the 11D state has six additional directions.

We test precisely:

1. `fixed_generalized_velocity`: perturb only the five configuration coordinates, holding all generalized velocity coordinates and the material-time accumulator fixed.
2. `fixed_cartesian_velocity`: choose the velocity-coordinate perturbations so that the **instantaneous Cartesian velocity** is also unchanged, `D_z V L=0`; material time remains fixed.

These are different, model-defined probes; neither is declared *the* physical material transport.

## Finite-lag tangent map and complement

For the declared finite xi-horizon `h=0.07`, define the midpoint approximation

\[ T_h(\xi)=J_q(z(\xi+h))\exp\!\{h A_{11}(\xi+h/2)\}\,L(z(\xi)). \]

This is a **finite-lag approximation** of the projected variational map and must not be misrepresented as an exact fundamental matrix. Let `U_h` denote its orthogonal polar factor. Define a skew rate by

\[K_h=(U_h-U_h^T)/(2h),\quad e=V/\|V\|,\quad P_e=I-ee^T,\quad B_h=P_e K_h P_e.\]

`B_h` is skew and annihilates `e`. The model-conditioned connection is

\[ \Omega_z = \Omega_{\rm min}+B_h,\quad \Omega_{\rm min}=e'e^T-ee'^T. \]

The additional rotation is integrated in the y minimum-turn co-moving frame using a Lie-group step. The transported 5D frame `R_z` acts through `Lambda²(R_z^T)` on the **same** original 10D physical residual paths as y.

## Comparisons and non-additivity

There are 204 original paths (two seeds × three kappa rules × two epoch metric calibrations × 17 lambda values), three temporal anchor assignments and two horizontal lifts, producing 1,224 rows. The polygon in `Lambda² R^10` is a 45-component oriented diagnostic, not a physical 45D configuration. Cross-polygon angles and arc lengths are not additively causal contributions.

Gauge tests impose an artificial `Q(ξ) in SO(5)` and verify return to the common frame plus covariance of the induced exterior-square generator. The alternative finite horizon `h/2` is compared on the authentic paths near lambda=4. Only the local spin field is horizon-refined, not a new nonlinear physical orbit.

## Scientific boundary

The generator is authentic to the source-conditioned model; the **complement projection is horizontal-lift dependent** and its discretization is finite lag. Without a physical principle selecting the lift, and without an independent material return section in the formal clock extension, one must report `NO_UNIQUE_MATERIAL_CONNECTION` and `NO_MATERIAL_HOLONOMY` regardless of numerical covariance success.
