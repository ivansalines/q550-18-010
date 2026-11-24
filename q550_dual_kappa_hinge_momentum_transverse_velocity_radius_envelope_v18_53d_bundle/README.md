# Q550 — v18.53d: Dual-κ₀ Hinge-Momentum, Transverse-Velocity & Radius-Envelope Audit

**Purpose:** test the suspicion that an inclined observer can see a projected return while the underlying movement advances in an unseen direction, and distinguish that effect from the two real rotor angular momenta along the authenticated original/autonomous κ₀ fold locus.

This ZIP is self-contained for the **v18.53c fold postprocessing**: authenticated report, selected states, complete CSV, output manifest and historical source scripts are included. The large SHA-256-pinned v18.52l physical relative 5D trace is an **external input**, not copied into the ZIP. This is a light, no-IVP postprocessing audit; it does not re-integrate Q550.

## On the server

```bash
cd /root/elab1/hysteresis
unzip q550_dual_kappa_hinge_momentum_transverse_velocity_radius_envelope_v18_53d_bundle.zip
cd q550_dual_kappa_hinge_momentum_transverse_velocity_radius_envelope_v18_53d_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
 bash run_on_server.sh --require-trace 2>&1 | tee ../v18_53d_server.log
```

The script auto-finds:

`/root/elab1/hysteresis/q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_5d_vectors_and_traces.npz`

Alternatively set the exact file explicitly using `--trace-npz /absolute/path/to/extended_horizon_5d_vectors_and_traces.npz`. The hash **must** match the pinned one.

If the file is missing, the script outputs an explicitly **partial** fold-only report and still creates a sendback. `--require-trace` exits nonzero on this incomplete result. Output is at the parent level: `q550_dual_kappa_hinge_momentum_transverse_velocity_radius_envelope_v18_53d_bundle_server_sendback.zip`.

## Files generated

- `results/fold_rotor_kinematics.csv`: all 166 authenticated fold points (102 in the dual-κ interval), including radii, angular and tangential velocities, individual rotor moments, conserved difference.
- `results/fold_rotor_report.json`: endpoints, differences, ranges, distinct signed-drive warning, and the inherited same-drive radius comparison.
- `results/physical_5d_relative_projection_shadow_hidden_depth.csv`: eight authentic 5D relative traces, five fixed observer inclinations, closest separated early/late screen pairs, corresponding full 5D gaps, unseen depth and velocity mismatch, and **kinematic** wedge values.
- `results/relative_5d_hinge_wedge_summary.json`: ranges of observer-plane radii and kinematic area rates; not canonical momentum.
- `results/dual_kappa_rotor_radius_speed_momentum.png` and `results/tilted_observer_shadow_vs_hidden_5d_gap.png`: diagnostic images.
- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, `OUTPUT_MANIFEST.json`: decision and provenance.

## Interpretation safety

The fold data are absolute relative-equilibrium rotor radii and angular rates. The 5D `dQ,dV` data are perturbation-minus-reference **separations** at a shared historical κ₀. They are **not** absolute 5D trajectories and cannot yield a canonical hinge charge without absolute coordinates, a canonical momentum map, and calibrated hinge geometry. The observed four-point slope has not been fit. Do not promote a near 2D screen return to a true 5D return or a global lower-radius bound.

Documentation and source comments are English as requested.
