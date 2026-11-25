# Q550 v18.53y — 11D Velocity-Generated Minimal-Turn Connection, Induced 10D Transport and Gauge Covariance

A self-contained, genealogy-preserving audit that constructs a *kinematic* SO(5) minimum-turn frame from an authenticated physical 11D trajectory, induces Lambda² frame changes on integrated 10D bivector residuals, and tests their dependence on temporal anchoring and imposed gauge. It does **not** claim a unique material connection or holonomy.

## Server execution

```bash
cd /root/elab1/hysteresis
B=q550_11d_velocity_generated_minimal_turn_bivector_connection_gauge_covariant_transport_v18_53y_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 bash "$B/scripts/run_on_server.sh" --require-parent 2>&1 | tee v18_53y_server.log
```

The server run authenticates inputs, runs synthetic tests, audits 204 metric continuation paths (612 start/midpoint/end lifts), performs gauge and refinement checks, creates CSV/NPZ/PNG/JSON output and packages the complete reproducible `..._server_sendback.zip` including Python sources, inputs, hashes and results. No physical 11D/21D IVP is repeated.

`--preflight-only` authenticates the lineage without loading the full 11D trajectory. `--smoke-only` completes a single three-lambda authentic case without writing final scientific outputs.

## Main outputs

- `results/dynamic_frame_oriented_transport.csv`: full 204×3 comparisons.
- `results/moving_frame_gauge_covariance.csv`: imposed time-dependent frame covariance checks.
- `results/frame_resolution_convergence.csv`: dense and stride-4 frame comparisons.
- `results/selected_dynamic_10d_45d_transport.npz`: selected authentic transported vectors and polygons.
- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`: scientific decision and limitations.
- `results/dynamic_connection_oriented_transport_original.png` and `_autonomous.png`: figures.

## Dependencies

Python 3.13; NumPy, SciPy, Matplotlib. The full v18.53x server sendback is included byte-for-byte with its nested authenticated ancestry. No external native Q550 source NPZ is required to run this analysis.
