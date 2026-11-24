# Q550 — v18.52o

**Full-Path Five-Dimensional Homothety & Gate-Conditioned Scaling Audit**

## Purpose

The authentic v18.52n server audit found approximately 3:1 physical 5D scaling between the inherited `1e-4` and `3e-4` perturbation amplitudes at matched gate crossings. This audit answers the *next independent question*: **Does the 3:1 relation hold throughout the entire material-time path, or does its accuracy improve systematically at the nodal-hinge gate crossings?**

No new Q550 nonlinear integration, 5D trajectory creation, or numerical optimization is performed. The scripts analyze pinned v18.52l physical vectors and reuse observational conventions from the genuine v18.52m/v18.52n sendbacks. The exact four coordinates that originally suggested slope `-0.308` are still not supplied, and are not invented here.

## Genuine source and genealogical pins

The required file is the same authenticated v18.52l source used in v18.52m/n:

`../q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_5d_vectors_and_traces.npz`

Expected source SHA-256:

`d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`

This source has 8 labeled traces with physical 5D `dQ,dV`, the material clock `tau`, and `D,G` diagnostics. Each inherited trace contains 65,537 samples. There are **four matched 1:3 perturbation pairs**, not eight independent 1:3 comparisons. Files in `reference/` pin the genuine v18.52l, m and n metadata. The bundle deliberately does NOT embed the ~111 MB parent source; it reads it from the server's existing genealogy.

## Run the scientific audit on the server

```bash
cd /root/elab1/hysteresis
unzip q550_full_path_five_dimensional_homothety_and_gate_conditioned_scaling_audit_v18_52o_bundle.zip
cd q550_full_path_five_dimensional_homothety_and_gate_conditioned_scaling_audit_v18_52o_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh 2>&1 | tee ../v18_52o_server.log
```

If the source path differs, pass `--source /absolute/path/to/extended_horizon_5d_vectors_and_traces.npz` to `bash run_on_server.sh`. The server script runs lineage verification, synthetic tests, authenticated source preflight, full-resolution analysis, a SHA-256 manifest, and then packages a sendback ZIP with the **same complete root name and descriptive ZIP filename**:

`q550_full_path_five_dimensional_homothety_and_gate_conditioned_scaling_audit_v18_52o_bundle_server_sendback.zip`

**Do not use `--sample-stride` for scientific analysis.** `--sample-stride 128` exists only for quick code testing on a reduced subset and is forcibly labeled `DEVELOPMENT_SMOKE_NOT_SCIENTIFIC`. The full run uses `--sample-stride 1` by default. `--output-stride 16` only controls the size of an exported chart-support CSV; metrics, events, gate windows, and all summary statistics are computed on the complete grid.

## Key scientific questions

1. At identical material time `tau`, does the physical vector satisfy `q(3epsilon,tau) ≈ 3 q(epsilon,tau)` throughout the horizon? Report the 5D vector residual, not just the ratio of 5D norms or a rotated 2D shadow.
2. Does the gate-conditioned residual become systematically **smaller** than an equal-phase midpoint baseline, non-gate samples, and matched low-amplitude-norm quartiles? For the midpoint-minus-gate signed differences, **positive** means gate improvement.
3. Are observer gate results robust to changing the hinge slope, transverse tilt, and hidden-axis direction? Gate locations themselves depend on the observer; the underlying 5D residual does not.
4. Do slight deviations of the material clocks or the choice of temporal interpolation explain an apparent scaling advantage? Explicit negative/control comparisons are included.
5. How large are the even (positive/negative kick symmetry) components relative to odd components? This is a diagnostic of deviations from an ideal odd first-order response, not an independent theorem.

## Output files

The `results/` directory includes:

- `SCIENTIFIC_REPORT_SERVER.md`, `report.json`: numerical summaries and limitations.
- `full_path_5d_homothety_and_gate_alignment.png`: 5D norm ratio, relative vector residual, global-normalized vector residual over the entire horizon.
- `gate_vs_midphase_relative_residual.png`: all signed gate-vs-midphase comparisons.
- `full_path_decimated_5d_scaling.csv`: time curves sampled for practical plotting, not used for metrics.
- `gate_event_5d_residuals.csv`: registered small-kick gate event residuals and true 5D orientations.
- `midphase_controls.csv`: matched midpoint of consecutive same-view gates.
- `gate_windows_vs_full_path.csv`: fixed gate-neighborhood widths and observer opposite-side window.
- `norm_quartile_conditioned_controls.csv`: controls for changing norms/relative denominators.
- `observer_slope_tilt_hidden_axis_sensitivity.csv`: slope, angle, and hidden-axis negative controls.
- `amplitude_dependent_gate_time_shift.csv`: large-minus-small timing drift of gate crossings.
- `perturbation_even_odd_diagnostics.csv`: even/odd positive/negative kick diagnostic.

All actual generating Python scripts, the source lineage, command recipes, SHA-256 hashes, dependency notes, and the authentic results are included in the sendback ZIP.

**Interpretation boundary:** Since the 3:1 amplitude ratio was specified at initialization, an accurately maintained ratio is evidence for *persistence of approximately linear response*, not a spontaneously generated self-similar topological invariant. The earlier four observational points behind the -0.308 slope remain a separately open provenance question. No double helix, tube, cone or hidden physical gate theorem is presumed.
