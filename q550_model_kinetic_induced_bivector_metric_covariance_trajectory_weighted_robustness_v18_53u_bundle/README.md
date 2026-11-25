# Q550 — v18.53u

## Model Kinetic-Induced Bivector Metric, Covariance & Trajectory-Weighted Valley Robustness Audit

**Goal.** Answer whether the full-10D residual valley and oriented transport remain recognizable when their inner product is derived consistently from the physical 5D kinetic geometry, rather than naively relabeling the historical Euclidean norm as arbitrary.

**Critical fact from the historical source.** The material position q already contains `sqrt(m_x)*x`. Its canonical kinetic metric is the Euclidean identity in the q coordinates, so the 10D bivector metric is also identity in the historical `sqrt(2)*B_ij` convention. A log-polar coordinate pullback MUST give the same measurement when the bivector is transformed covariantly. A nontrivial alternative is explicitly labeled as a *model-conditioned kinetic-occupation sensitivity*, not a model-imposed canonical metric.

## Run on the Q550 server

```bash
cd /root/elab1/hysteresis
B=q550_model_kinetic_induced_bivector_metric_covariance_trajectory_weighted_robustness_v18_53u_bundle
unzip -q "$B.zip"
PYTHON_BIN=python3.13 OPENBLAS_NUM_THREADS=1 bash "$B/scripts/run_on_server.sh" --require-parent 2>&1 | tee v18_53u_server.log
```

No new 11D or 21D IVPs are needed. This audit processes the authentic 8961-node on-the-fly v18.53t states, cumulative bivectors and physical source. It generates the normal `..._bundle_server_sendback.zip` with CSV, plots, NPZ, README, scripts, pinned original t sendback, and verified manifests. Do not rerun historic IVPs.

## Case design

- 2 historical seeds × 3 historical κ0 rules = 6 authentic physical cases.
- Canonical physical kinetic metric `I5` → induced `I10`.
- 2 frozen kinetic-occupation epochs (`xi<=25.6`, entire formal trace), each at λ=1 and λ=3; five metric families per case, 30 case/metric rows.
- 141 origin points per case/metric curve (identical `xi0` interval in each historical case). Re-optimize width in the historical ±0.05 band.
- Pullback identities in the log-polar/material coordinates using exact historical 5D embedding.
- Model-velocity identity verified against byte-identical, SHA-pinned `torque_math.py` extracted from the authenticated v18.53k ancestry within v18.53t.

## Outputs

`results/metric_valley_summary.csv`, `results/metric_valley_curves.csv`, `results/metric_comparison.csv`, `results/metric_oriented_10d_data.npz`, `results/kinetic_metric_valley_<seed>.png`, `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, `OUTPUT_MANIFEST.json`.

`--preflight` verifies the outer authentic v18.53t parent CRC, manifest SHA entries, expected decision and physical embedding matching without new physical IVPs. `tests/test_metric_geometry.py` checks mathematical identities on synthetic data. Source code and comments intentionally stay in English.

**Limitations:** This uses the formal clock continuation beyond xi=25.6, and it is not a certificate of physical holonomy, full-state recurrence or a fundamentally new kinetic Lagrangian.
