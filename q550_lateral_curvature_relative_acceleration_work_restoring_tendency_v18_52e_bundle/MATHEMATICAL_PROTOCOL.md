# Q550 v18.52e — Mathematical protocol

## Separate three distinct effects

For each fixed-parameter physical velocity kick, reconstruct the authentic nonlinear 10D trajectories on a common elapsed material clock `tau`. Let `Q,V,a=dV/dtau` be the kinetic Cartesian embedding and its first two physical-time derivatives. Relative variables use perturbation minus its **own** unperturbed baseline at the same tau:

`dQ=Q_kick-Q_base`, `dV=V_kick-V_base`, `da=a_kick-a_base`.

The following **identities** hold at matched material tau:

- Relative velocity norm: `K_rel = 1/2 |dV|^2`; `dK_rel/dtau = dV·da = P_rel`.
- Relative position separation: `D_rel = 1/2 |dQ|^2`; `dD_rel/dtau = dQ·dV = G_rel`.
- Second derivative of position separation: `dG_rel/dtau = |dV|^2 + dQ·da`.

`P_rel<0` indicates instantaneous decrease of **squared relative speed**, NOT dissipated physical energy or whole-state recovery. `G_rel<0` indicates instantaneous contraction of squared relative position. `dQ·da<0` is only a partial inward-acceleration projection; it does NOT by itself imply that `G_rel<0` because the positive `|dV|^2` term contributes.

## Lateral curvature

Reference tangent `t=V/|V|`. Geometric acceleration `a_perp=a-(a·t)t`, curvature magnitude `k_g=|a_perp|/|V|^2`. This is geometric curvature in the explicitly specified five-dimensional kinetic Cartesian chart, not a coordinate-free curved-manifold intrinsic without a declared metric. Compare reference/kick curvature and their signed difference. Also decompose `da` relative to **baseline tangent**. Distinguish that from the difference between two trajectories' own lateral accelerations.

## Clock and interpolation

All vector comparisons are at common elapsed material time `tau`, using cubic splines of original and kicked trajectories, each with its own monotonically increasing material clock. External drive phases may differ at fixed tau; record the phase offset. Audit independently the spline derivative of `1/2|dV|²` and of `1/2|dQ|²` against the model vector identities, with precision-aware tolerances. Do not confuse derivatives on an interpolation mesh with validated interval bounds.

## Recovery classifications

- `RELATIVE_VELOCITY_CONTRACTION_WITH_POSITION_DIVERGENCE`: predominantly negative `P_rel`, predominantly positive `G_rel`, and terminal `|dV|` smaller than initial while terminal `|dQ|` larger.
- `POSITION_RETURN_CANDIDATE`: a nontrivial interval in which `G_rel<0` and positional distance decreases after a previous maximum.
- `NO_FINITE_WINDOW_POSITION_RETURN`: no such interval in the sampled window.

No finite-window test certifies asymptotic return, dynamical stability, dissipative energy loss, or constraint-preserving global perturbations. Evaluate both signs and both kick amplitudes for both physical seeds.

## Scientific falsification

Numerical PASS means interpolation identities, source genealogy, fixed-tau comparisons, and real nonlinear trajectories are consistent. Scientific status is separately reported, even if `P_rel` becomes positive or some perturbation fails to recover. All source code and comments are in English.
