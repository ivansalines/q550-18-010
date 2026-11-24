# Q550 — v18.52n

**Four-Point Nodal-Axis Reconstruction & Paired-Hinge Observer Audit**

## Why this audit exists

The observational proposal was refined: an apparent nodal seam may have a **nonzero slope**

`m = Q2/Q1 approximately -0.308`, corresponding to `atan(m) approximately -17.119 degrees`.

The observer is additionally proposed to be tilted **around that sloped seam** by `+17.04` and `-17.04` degrees (a physical sheet hinge, NOT an in-plane rotation of an existing 2D plot). These are two independent degrees of freedom. For reference, `tan(-17.04 degrees) ≈ -0.30649`.

**A crucial provenance limit:** the actual four coordinate pairs used by the observer to infer `-0.308` have NOT been supplied. In particular, the four genuine `Q2=0` section crossings from v18.52l are *not* evidence for the nonzero slope: their y-coordinate is zero by construction. The default audit therefore labels `m=-0.308` as a *predeclared observer hypothesis*, not a reconstructed four-point fit. It independently produces the historical four-crossing negative control.

## Actual physical data (no new integration)

Parent v18.52l authenticated physical 5D relative traces, eight kick labels, 65,537 material-time samples per trajectory, `dQ` and `dV`, with inherited `D` and `G` identity checks.

Default source relative to this bundle:

`../q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_5d_vectors_and_traces.npz`

Expected historical SHA-256:

`d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`

The parent v18.52m authentic server report and v18.52l/m output manifests are pinned in `reference/`. The v18.52m learned hidden PCA axis is recomputed and compared to the historical result.

## Full server execution

```bash
cd /root/elab1/hysteresis
unzip q550_four_point_nodal_axis_reconstruction_and_paired_hinge_observer_audit_v18_52n_bundle.zip
cd q550_four_point_nodal_axis_reconstruction_and_paired_hinge_observer_audit_v18_52n_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_52n_server.log
```

If the parent file is located elsewhere, append

`--source /absolute/path/to/extended_horizon_5d_vectors_and_traces.npz`

to the `bash run_on_server.sh ...` command.

The bundle is **small** and does not contain the large v18.52l physical data or re-integrate the original model. Running the whole analysis on 65,537 samples per kick is the authoritative server audit. `--sample-stride 128` is **development smoke only** and writes an explicit non-scientific stage marker. Do not mistake smoke results for server evidence.

## Optional exact four-point reconstruction

Populate a separate file (not the header-only template) containing **exactly four physical 2D coordinate pairs**:

```csv
point_id,Q1,Q2
P1,<coordinate1>,<coordinate2>
P2,<coordinate1>,<coordinate2>
P3,<coordinate1>,<coordinate2>
P4,<coordinate1>,<coordinate2>
```

Replace the `<...>` placeholders with floating-point numbers; this example is **not** an executable four-point file. The program will report an explicit affine OLS fit, constrained-through-origin OLS fit, total least squares orientation, four leave-one-out affine slopes, residuals, and axis angle difference. It does *not* silently change the fixed observer axis:

```bash
OPENBLAS_NUM_THREADS=1 bash run_on_server.sh \
  --four-points /absolute/path/to/observed_four_points.csv
```

If and only if the user wants the geometry to use the **fitted affine line** `Q2=m*Q1+b`, add `--adopt-four-point-fit` explicitly. `--intercept` allows a specified fixed intercept without adopting a fit. Four points are observational coordinates, not an automatic proof of a material invariant.

## Geometry and controls

See `MATHEMATICAL_PROTOCOL.md`. The audit compares these independent hinge axes:

- original v18.52m `m=0` Q1-axis observer;
- user proposed `m=-0.308` observer;
- exact slope `tan(-17.04°)` for comparison;
- opposite-signed `m=+0.308` control.

Each is evaluated at transverse tilt `0°, 10°, 17.04°, 25°, 35°`, with hidden axes Q3/Q4/Q5 and a fixed PCA direction trained only on the first half of a single historical primary trace. It measures interval-resolved opposite-side gate windows, crossings and actual 5D posture, serial section return, intrinsic 5D radius about the candidate line, apparent ellipse axis ratios, and cone-fit diagnostic. None of these by itself proves a double helix, tube, cone, topological invariant, or autosimilar manifold.

## Outputs

`results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, `results/four_point_reconstruction.json`, `results/historical_four_Q2_zero_negative_control.json`, 6 CSV diagnostics, two plots (plus optional four-point fit plot), `results/PREFLIGHT.json`; `OUTPUT_MANIFEST.json` hashes reproducible sources and results. The generated sendback ZIP includes all scripts, README/protocol, reference provenance and outputs, using the complete matching root/filename.

No Python script needed to reproduce the geometry has been omitted from this bundle.
