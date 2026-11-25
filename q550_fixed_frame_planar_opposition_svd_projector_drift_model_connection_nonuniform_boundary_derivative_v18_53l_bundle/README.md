# Q550 — v18.53l
## Fixed-Frame Planar Opposition, SVD-Projector Drift, Model-Connection & Nonuniform Boundary-Derivative Audit

This bundle tests the observation raised by the independently reviewed v18.53k experiment: around three equal windows of width 12.8, the moving-SVD planar projected residual at widths 12.775 and 12.825 is almost oppositely directed, whereas the third singular contribution remains similarly oriented. The central falsification question is whether this apparent inversion **survives the refusal to rotate the observer plane to chase each new width**.

### Execution on the Q550 research server

```bash
cd /root/elab1/hysteresis
B=q550_fixed_frame_planar_opposition_svd_projector_drift_model_connection_nonuniform_boundary_derivative_v18_53l_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
bash "$B/scripts/run_on_server.sh" \
  --require-parent --require-native --require-dynamic \
  2>&1 | tee v18_53l_server.log
```

If the two large parent files are stored elsewhere, use `--parent-v1851s /path/to/authentic_parent.zip --native-npz /path/to/authentic_10D.npz`. These **two original large historical artifacts are not packaged** but their SHA-256 values are pinned in `inputs_config.json`. The same authenticated discovery mechanism used by v18.53k is retained.

A **light local verification** without physical integration:

```bash
bash run_on_server.sh --preflight
```

The `--preflight` mode explicitly reports `PREFLIGHT_PASS__NO_NEW_Q550_IVPS`. If a mandatory server parent is absent or has a different digest, physical execution exits nonzero and labels the resulting report as incomplete. It never silently substitutes synthetic data.

### Pinned genealogy

- Full original v18.53k authenticated server sendback, embedded unmodified in `inputs/v18_53k_authenticated_server_sendback.zip`; v18.53k source/output manifests and all 756 original rows checked at runtime.
- Embedded v18.53k includes the v18.53j sendback, the earlier code ancestry, the native loader, material 5D embedding, torque decomposition, report, and regression checks.
- 10 Q550 IVPs are recomputed on the user's server: six unique seed/κ₀ trajectories and four refined repeats. Widths **do not** correspond to additional independent dynamics: they postprocess the same six trajectories.
- Two historical seeds (`original`, `autonomous`), three κ₀ rules (`shared`, `original`, `autonomous`), six origins and 21 widths exactly as v18.53k: 756 triplets, plus 36 paired-width comparisons.

### Three distinctly defined 10D bivector-space observers

**Fixed:** top-two rank-2 projector of the three-residual SVD evaluated at width 12.8, held constant along the whole width scan for each seed/κ₀/origin.

**Moving SVD:** top-two right singular projector of the actual three interval-integrated 5×5 skew bivectors, rebuilt at each width. Depends explicitly on the chosen windows.

**Model-induced transported:** apply a declared minimal-radial SO(5) kinetic-frame lift derived from the integrated 5D position Q(xi), and its induced exterior-square SO(10) rotation, to the anchor projector. Reference and target are evaluated at the right endpoint xi=origin+3*width. This is one reproducible model-based *kinematic* connection. It is **not** evidence of a unique model-canonical material connection, or of physical holonomy.

All observers share the same Frobenius-isometric ambient 10D coordinates. Rotations used for the third observer are explicitly exported by the `frame_projector_geometry.py` generator, not guessed from the SVD result.

### Separation of changes

For two consecutive widths `a,b`, compute the **exact finite secant identity**:

`P_b N_b - P_a N_a = [(P_a+P_b)/2] (N_b-N_a) + (P_b-P_a) [(N_a+N_b)/2]`.

The first term is the true integrated net variation projected by an average operator; the second is variation induced by changing the projector. The average operator is NOT claimed to be an orthogonal projector. For the fixed observer, the projector drift term is exactly zero.

For net bivector `N(w)`, the boundary derivative of the three-window *sum* is the exact FTC identity

`dN/dw = 3 B_tau(origin+3w) * d tau/d xi (origin+3w)`.

Test both three-node and five-node **nonuniform** finite-difference weights, with an independent endpoint source interpolated from genuine Q550 traces. At adjacent widths additionally reconcile the net difference to the accumulated integral between the corresponding endpoints. Do not interpret a coarse nonuniform-stencil error as physical motion.

### Report and sendback

The server writes:

- `results/three_plane_comparison.csv` — all 756 comparisons, 10D net and three projected vectors, rank gaps and frame conditioning.
- `results/fixed_vs_moving_vs_model_opposition.csv` — 36 paired-width (12.775/12.825) cosines, model-pullback comparisons.
- `results/exact_secant_boundary_vs_projector_drift.csv` — 720 exact plane/projector-change decompositions × 3 observers = 2,160 rows.
- `results/nonuniform_endpoint_derivative.csv` — 3-node and 5-node FD vs exact FTC endpoint source.
- `results/genuine_integrated_ivps.csv`, `results/refined_integration.csv` — solver genealogy.
- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, diagnostic plots, and source/output SHA-256 manifests.

Expected sendback filename:

`q550_fixed_frame_planar_opposition_svd_projector_drift_model_connection_nonuniform_boundary_derivative_v18_53l_bundle_server_sendback.zip`

**Scientific limits.** For xi>25.6 the external clock is model continuation rather than independent historical observation. No claim of a physical flip, exact three-window closure, new spatial axis, unique dynamical connection, or material holonomy follows automatically from the present study. See `SCIENTIFIC_DECISION_RULES.md`.
