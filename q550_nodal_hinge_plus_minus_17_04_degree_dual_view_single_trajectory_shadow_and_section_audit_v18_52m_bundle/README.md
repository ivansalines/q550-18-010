# Q550 — v18.52m

**Nodal Hinge ±17.04° — Dual-View Single-Trajectory Shadow & Section Audit**

## Scientific question

Could the apparent twin lobes of Q550 be **two observation shadows of a single physical 5D relative trajectory**, revealed by tilting a sheet around the projected nodal cut like a flag, rather than by rotating a graph inside its plane?

This is a **postprocessing audit only**. There is no fresh nonlinear integration or change to forcing, kappa, global phase, physical metric, kinetic embedding, material clock, or 10D state.

## Hinge geometry

For each authentic physical 5D relative vector `q(t) = Q_kick(t) - Q_reference(t)` on the same matched material time:

- Prespecified candidate hinge: physical `Q1` axis, `e1=(1,0,0,0,0)`.
- Original visible transverse axis: `e2=(0,1,0,0,0)`.
- Hidden normalized direction `h` in physical span `(Q3,Q4,Q5)`.
- For `theta = 17.04 degrees`, observer axes are
  `e2_plus = cos(theta)*e2 + sin(theta)*h`,
  `e2_minus = cos(theta)*e2 - sin(theta)*h`.
- Both observer planes share exactly the same fixed `Q1` hinge.
- Each reported pair of 2D curves is evaluated at the **same sample times of the same single 5D physical trajectory**.

The original horizontal `Q2=0` cut is a **4-dimensional hyperplane in physical R5**, not a unique physical 1D hinge. Choosing the `Q1` line as hinge is a *hypothesis to audit*, not a certified dynamical axis. The 17.04-degree candidate is predeclared and not fitted.

## Independent controls and exclusions

1. At fixed theta, compare the hidden axes `Q3`, `Q4`, `Q5`, and an explicit PCA hidden direction.
2. The PCA direction is trained **only on the first half** of the `original_negative_1.0e-04` trajectory and frozen for all other trajectories; report the held-out half's principal-axis angle.
3. Test theta = 0, 10, 17.04, 25, 35 degrees, with both signs. Zero tilt is a degeneracy control.
4. Compare cross-view crossing pairs AND consecutive crossings within each fixed view using full 5D posture, radial D ratios, and antipodal residuals; never claim return from apparent overlap.
5. Monitor observer-dependent projected-phase turns, intrinsic 5D hinge-perpendicular radius, axial drift, and opposing-view gate windows. These are **diagnostics**, not a topological invariant or cone/tube proof.
6. A planar synthetic trajectory must produce identical +/- shadows. A nonplanar synthetic trace must separate them.

**Important:** The identities `w_+ + w_- = 2*cos(theta)*Q2` and `w_+ - w_- = 2*sin(theta)*(h·Q_hidden)` are guaranteed by construction and are not physical symmetry evidence.

## Input provenance

The *actual* v18.52l server output is expected in the sibling directory:

`../q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_5d_vectors_and_traces.npz`

The source must match SHA-256 from the exact historical `OUTPUT_MANIFEST.json` pinned as `reference/PARENT_v18_52l_OUTPUT_MANIFEST.json`:

`d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`

The input is **not copied** into this small bundle: the full 5D file remains with the user's authenticated v18.52l server sendback. If the sibling path differs, override `--source`.

## Server commands

From `/root/elab1/hysteresis`:

```bash
unzip q550_nodal_hinge_plus_minus_17_04_degree_dual_view_single_trajectory_shadow_and_section_audit_v18_52m_bundle.zip
cd q550_nodal_hinge_plus_minus_17_04_degree_dual_view_single_trajectory_shadow_and_section_audit_v18_52m_bundle

OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh \
  2>&1 | tee ../v18_52m_server.log
```

For a different source location:

```bash
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh --source /path/to/extended_horizon_5d_vectors_and_traces.npz
```

## Outputs

- `results/report.json` and `results/SCIENTIFIC_REPORT_SERVER.md`;
- `results/paired_hinge_views_single_trajectory.png`;
- `results/single_trajectory_reduced_3d_and_hinge_radius.png`;
- `results/paired_observer_gate_windows.png`;
- `results/angle_and_hidden_direction_sensitivity.csv`;
- `results/tilted_hinge_section_crossings_5d.csv`;
- `results/cross_view_sequential_posture_pairs_5d.csv`;
- `results/single_observer_consecutive_section_return_map_5d.csv`;
- `results/single_trajectory_hinge_observer_metrics.csv`;
- `SOURCE_MANIFEST.json` and `OUTPUT_MANIFEST.json`, plus compact ZIP sendback retaining scripts and reports.

The script never asserts that two lobes, a double helix, autosimilar manifolds, or a topological invariant have been certified. Such interpretations remain conditional on coordinate robustness and future dynamical/topological diagnostics.
