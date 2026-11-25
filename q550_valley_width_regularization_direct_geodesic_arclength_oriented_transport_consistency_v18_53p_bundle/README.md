# Q550 — v18.53p: Optimized-Valley Width Regularity, Direct Geodesic Arclength & Oriented Transport Consistency Audit

## Scope

This is a **derived numerical convergence audit**, using the complete, SHA-256 pinned, physical v18.53o server sendback as its sole scientific input. **No new Q550 physical IVPs are performed or claimed.** The v18.53o sendback contains six prior physical cases and four refined IVPs, together with the full executable ancestry reaching back to historical model solvers. The input is retained byte-for-byte in `inputs/`.

The motivating discrepancy in v18.53o is that direct adjacent unit-vector arcs yield approximately **0.0098–0.0099 rad more angular path length** than a trapezoidal integral of the tangent speed along locally reoptimized-width valley curves, while fixed-width curves have substantially better consistency. This difference is not evidence of micro-holonomy.

## Scientific tests

1. Reconstruct each full 10D residual `N(o,w)=C(o+3w)-C(o)` from authenticated cumulative primitives.
2. Recompute the local width minimizer `w(o)` at 141, 281 and 561 origin samples (optional 1121 with `--deep`). For piecewise-linear `C`, search **every intersected endpoint segment analytically** for a global bounded minimum. For PCHIP `C`, use bracketed optimization over all sampled candidate valleys.
3. Compare the sum of **direct great-circle angles** of consecutive normalized 10D residuals against the trapezoidal integral of (a) chain-rule angular speed with finite-difference `w'`, (b) the same with interpolated PCHIP `w'`, and (c) finite-difference unit-vector speed.
4. Compare 45-component oriented two-forms: a geodesic segment integral vs the instantaneous `u ∧ u'` trapezoidal integral. Inspect convergence and any change in orientation sign/coherence.
5. Measure width total variation, derivative jumps, sampled clock-knot contacts, radial fundamental-theorem consistency, and stability under alternative interpolants.
6. Replay and check the original v18.53o 141-point curves and report any change caused by using the **global** piecewise-linear minimizer rather than the parent scalar optimizer.

## Execute on the user's server

```bash
cd /root/elab1/hysteresis
B=q550_valley_width_regularization_direct_geodesic_arclength_oriented_transport_consistency_v18_53p_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
  bash "$B/scripts/run_on_server.sh" 2>&1 | tee v18_53p_server.log
```

Optional deeper (derived) resolution:

```bash
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
  bash "$B/scripts/run_on_server.sh" --deep 2>&1 | tee v18_53p_deep_server.log
```

The script automatically generates a full-name ZIP ending in `_bundle_server_sendback.zip` **next to** the extracted bundle. A prior run's output will be replaced on a second run. The corrected packager uses `Path(root.name) / relative_path` rather than concatenating strings and `PosixPath` objects.

No original 11D parent NPZ is required for this **post-processing only** run because the exact physical-derived 10D primitives are already pinned within the v18.53o sendback. NumPy, SciPy and Python 3.13 are required; Matplotlib is optional for plots.

## Local validation without real Q550 scans

```bash
python3.13 scripts/verify_manifest.py
python3.13 tests/test_synthetic.py
python3.13 tests/test_authentic_smoke.py
python3.13 scripts/run_audit.py --preflight
```

Synthetic tests exercise geometry; the smoke test performs seven small objective optimizations on **already authenticated data**, no physical ODE integration.

## Results

- `arclength_convergence_summary.csv`: 6 cases × 2 primitive interpolants × 3 grid resolutions × 2 curves = **72 summaries**, or 96 with optional deep resolution.
- `width_and_direction_samples.csv`: origin, optimized width, 10D norm, angular speed, direct-arc increments and clock-knot clearance.
- `resolution_pair_comparison.csv`: changes across grid resolutions.
- `historical_v18_53o_replay.csv`: original discrepancies and read-only replay controls.
- `unit_direction_and_oriented_area.npz`: all unit vectors, origin grids, widths, and 45D integrated oriented forms for independent reconstruction.
- `report.json`, `SCIENTIFIC_REPORT_SERVER.md`, two PNG convergence plots when Matplotlib is available, and `OUTPUT_MANIFEST.json` in sendback.

## Interpretation limits

This audit is not a new physical evolution and does not establish a real-material periodic return, micro-holonomy, an internal spin, phase locking with the separate coupled-oscillator model, or closure in the full 11D state. The high-ξ endpoints use formal continuation beyond ξ=25.6. Interpolation changes are observer/numerical-method changes, not independent physical seed histories. All scientific reports and code comments remain in English.
