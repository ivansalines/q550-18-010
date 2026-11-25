# Q550 — v18.53x
**Metric-Identification Gauge Dependence & Fixed-Physical-Frame Oriented Transport Audit**

This bundle studies whether the ~6.56-degree parent v18.53w change of the 45D oriented polygon is intrinsic, or depends on the orthonormal identification of different kinetic-sensitivity metrics. The mathematical comparison is rooted in the model's historical mass-normalized Cartesian 5D kinetic coordinates.

## Inputs and genealogy
- `inputs/v18_53w_authenticated_server_sendback.zip`: exact original v18.53w sendback, with all 204 parent nodes and the nested v/u/t physics ancestry; pinned by SHA-256.
- `inputs/u_metric_geometry_authenticated.py`: exact historic metric helper; SHA-256 pinned.
- `SOURCE_MANIFEST.json`: exact SHA-256 of all distributed scripts, config, documents and input archives.

## Run on the server
```bash
cd /root/elab1/hysteresis
B=q550_metric_identification_gauge_dependence_fixed_physical_frame_oriented_transport_v18_53x_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 \
  bash "$B/scripts/run_on_server.sh" --require-parent \
  2>&1 | tee v18_53x_server.log
```

`--preflight-only` runs source, geometry and ancestry tests only. `--smoke-only` checks three authentic parent nodes and the moving-frame logic. Full mode uses the authenticated 21D source **without new physical IVPs** and produces a sendback ZIP next to the extracted bundle.

## Outputs
- `metric_identification_gauge_scan.csv`: 816 comparisons (204 parent nodes x 4 SO(5)-induced gauges); identity, intraplane and two mixed gauges.
- `fixed_physical_frame_path.csv`: 204 45D polygons in one canonical mass-normalized physical frame.
- `moving_center_fixed_physical_frame.csv`: 60 selected moving-center curves reconstructed from the authentic primitive and parent width minimizer.
- `moving_gauge_connection_checks.csv`: 24 explicit pure-gauge connection checks (12 families x 2 metrics).
- `selected_fixed_frame_vectors.npz`: selected normalized 10D paths and 45D polygons in the fixed/gauge frames.
- `gauge_fixed_frame_comparison_original.png`, `gauge_fixed_frame_comparison_autonomous.png`: visual comparison.
- `report.json`, `SCIENTIFIC_REPORT_SERVER.md`, `OUTPUT_MANIFEST.json`: classifications, precise errors and SHA-256 hashes.

All scientific content, formulas, documentation and script comments are in English; run instructions can be used as written from the server workspace.

## Interpretation discipline
A cross-metric angle is defined only after a choice of identification between the orthonormal metric frames. A large change of that angle under a pure-gauge change is *not* a physical rotation of the 11D state. Fixed-frame path geometry is a more direct diagnostic of motion, but still not an invariant material transport theorem.
