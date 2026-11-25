# Q550 — v18.53k
## Residual-Bivector Transverse SVD Direction Transport, Fine Window-Width Scan & Boundary Attribution Audit

**Purpose.** Test whether the nearly 2D equiangular three-bivector configuration at Δξ=12.8, and the small residual along the third 10D bivector-coordinate SVD direction, survives narrow changes of equal window width, as well as time-origin shifts. Distinguish rank-3 projection, actual integrated torque, and purely algebraic endpoint effects. This is a *falsification* audit, not a certification of a physical 3D motion axis, genuine holonomy, or material recurrence.

**Genealogy and engines.** Contains the entire authenticated v18.53j server sendback as an immutable pinned ancestor; the original v18.53j reference scripts; h/f/g/i reference code and source audit inputs; and the independent j SVD-review report. The source manifest authenticates all bundled scripts/inputs, and the internal j source/output manifests are recomputed at run time. Real Q550 integration still depends on these two LARGE historical artifacts already on the research server (SHA-256 values pinned in `inputs_config.json`):

- v18.51s original source parent ZIP (SHA `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`)
- v18.52l native 10D archive NPZ (SHA `92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12`)

The runner automatically looks for them using the same loader as v18.53j and supports `--parent-v1851s` and `--native-npz` for explicit paths. It refuses mismatched ancestry.

### Run on server

```bash
cd /root/elab1/hysteresis
B=q550_residual_bivector_transverse_svd_direction_transport_fine_window_width_scan_boundary_attribution_v18_53k_bundle
unzip -q "$B.zip"
set -o pipefail
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 bash "$B/scripts/run_on_server.sh" \
  --require-parent --require-native --require-dynamic \
  2>&1 | tee v18_53k_server.log
```

If file discovery fails, explicitly specify locations at the end of the command, for example `--parent-v1851s /path/to/verified-parent.zip --native-npz /path/to/verified-native.npz`. Run `bash "$B/run_on_server.sh" --preflight` for a light ancestry/configuration/synthetic-data check. **Preflight does not integrate the model** and is not a physical result.

### Predeclared scope

- Six unique Q550 IVPs (two historical seeds × three κ₀ choices), four tightened integration repeats, same genuine dynamical model and material clock as j.
- Window origins ξ = 12.8, 13.6, 14.0, 14.4, 14.8, 15.2.
- Widths 11.8, 12, 12.2, 12.4, 12.6, 12.65, 12.7, 12.725, 12.75, 12.775, **12.8**, 12.825, 12.85, 12.875, 12.9, 12.95, 13.0, 13.1, 13.2, 13.4, 13.6.
- 6 × 3 × 6 × 21 = **756** triplets; j baseline coordinates and 72 parent-j coarse-scan data points reconciled (four original origins × three matching widths × six seed/κ₀ configurations).
- Full interval ends at ξ ≤ 56.0; past ξ=25.6 is an externally continued historical clock and model trajectory, not new independent historical measurement.

### Expected sendback

`q550_residual_bivector_transverse_svd_direction_transport_fine_window_width_scan_boundary_attribution_v18_53k_bundle_server_sendback.zip`

Contains `results/report.json`, `SCIENTIFIC_REPORT_SERVER.md`, all 756 triplet rows in `fine_width_svd_geometry.csv`, candidate signed changes, sampled minima, coarse j reconciliation, IVP logs, integrator refinement, quadrature decimation, offgrid checks, plot PNGs, scripts, self-contained scientific protocol, and SHA256 ledgers. Outputs from a missing-input run are *partial*, labeled honestly; the script exits nonzero when required inputs are absent.

**Known limitation.** An O(10) SVD direction v₃ of three bivector coordinates is an observation-space line. It is not a third local physical Cartesian coordinate. Near-equal leading singular values do not invalidate the *rank-2 plane* when the σ₂−σ₃ gap is stable; near-vanishing σ₃ does invalidate the oriented v₃ coefficient. Frame continuity checks explicitly record such cases.

To replay the predecessor independent-SVD analysis using the included source script:

```bash
python3.13 reference/v18_53j_independent_svd_review.py --sendback inputs/v18_53j_authenticated_server_sendback.zip --out results/recomputed_j_svd
```
