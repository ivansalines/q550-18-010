# Q550 — v18.53o: Full-10D Residual Valley Oriented Direction & Boundary-Source Transport Audit

## Objective

Follow the **unprojected** 10D bivector residual through the numerical near-return valley found by v18.53n-fix1. Separately measure (i) scalar norm, (ii) unit direction and its 10D turning, and (iii) the physically reconstructed integrand at both temporal endpoints. The orientation analysis uses no SVD moving-plane sign choice.

### Mathematical contract

For origin `o`, three-window width `w`, and the accumulated bivector `C(xi)`:

`N(o,w) = C(o+3*w)-C(o)`, `r=||N||`, `u=N/r`.

With `g=C'`, `J_o=g_end-g_start`, `J_w=3*g_end`.

`d_o r = u·J_o`, `d_w r = u·J_w` and `d_o u=(J_o-u*d_o r)/r`.

At a true unconstrained two-parameter stationary point with `r>0`, `N·g_start=N·g_end=0`. These are geometric endpoint orthogonality conditions, **not zero residual and not a material return**. A nonzero `d_o u` at that point is compatible with a locally stationary norm.

## Data provenance

- The exact authenticated v18.53n-fix1 *physical server sendback* is embedded, unchanged, including the full m/l/k/j/... executable ancestry.
- Requires the user's original v18.51s parent source and native Q550 NPZ (hard SHA-256 pins, not bundled). No synthetic stand-in is allowed during server execution.
- All new quantitative results use six genuine 11D Q550 integrations and four tighter IVP controls, rather than re-labeling archived CSV rows as a new physical run.
- A new `derived_primitives.npz` stores the six derived 10D cumulative curves, plus their sampled physical sources and the common xi grid, to support independent lightweight replay later.

## Main measurements

1. **Fixed-width curve:** origin spanning ±0.35 around each parent n-fix1 minimum, with w held at the parent optimal value; records 10D unit direction, scalar norm, angular turning rate, directed 2-form of orientation and contributions from both endpoints.
2. **Local valley-following curve:** bounded scalar optimization of w for each origin, with analytic chain-rule including w'(o); distinguishes rotation along the valley from rotation in a fixed-width slice. Boundary width candidates are marked.
3. **Endpoint orthogonality:** independently test dot products of N with g_start/g_end, and the gradient J^T N, at the parent optimum.
4. **Numerical controls:** compare primary, decimated quadrature, PCHIP primitive interpolation, and four tighter physical trajectories at the same parent minimum. Report angle and norm changes, without pretending PCHIP is another dynamical model.
5. **Orientation gauge:** global ambient orthogonal-coordinate changes leave all scalar norms, angular speeds and cosines invariant. No transported-frame/geometric holonomy is claimed.

## Server run

```bash
cd /root/elab1/hysteresis
B=q550_full_10d_residual_valley_oriented_direction_boundary_source_transport_v18_53o_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
  bash "$B/scripts/run_on_server.sh" --require-parent --require-native --require-dynamic \
  2>&1 | tee v18_53o_server.log
```

Creates `q550_full_10d_residual_valley_oriented_direction_boundary_source_transport_v18_53o_bundle_server_sendback.zip` beside the root directory.

## Local-only checks (no physical Q550 IVPs)

```bash
python3.13 scripts/verify_manifest.py
python3.13 tests/test_synthetic.py
python3.13 scripts/run_audit.py --preflight
```

## Scientific limitations

- The trajectory after `xi=25.6` is the **formal continuation** of the archived numerical model.
- 10D bivector-return metrics are not a physical closure of all 11D state variables.
- Derivative identities for piecewise-linear cumulative data are exact **within interpolation cells**, not certified physical derivatives at knots.
- The selected valley is local, not a proof of global uniqueness; its width optimization must report any active boundaries.
- No inferred physical memory, phase-locking of the coupled oscillators, or actual mechanical spin is established by an oriented `u ∧ u'` in bivector-coordinate space.
