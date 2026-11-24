# Q550 — v18.51y

**Model-Constrained Source Subspaces & Physical Attribution Admissibility Audit**

## Research question

Which assignments of parameter-forcing contributions preserve the historical model's velocity-component source supports, and which changes are only coordinate relabellings? This version continues v18.51x after v18.51w, v18.51t and v18.51s.

## Reproduction

```bash
cd /root/elab1/hysteresis
unzip q550_model_constrained_source_subspaces_physical_attribution_admissibility_v18_51y_bundle.zip
cd q550_model_constrained_source_subspaces_physical_attribution_admissibility_v18_51y_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_51y_server.log
```

The full package embeds the verbatim X server sendback. The lightweight package contains SHA-256-locked T/W/X source and result snapshots plus the parent scripts and reports, without the 140 MB X ancestry. Both generate the same numerical analysis; the full package also independently verifies the X parent's nested CRC, source/output hashes and provenance. To add full genealogy to a lightweight unpacked bundle, put the original `q550_physical_source_channel_recombination_geometric_attribution_robustness_v18_51x_bundle_server_sendback.zip` in `inputs/` and rerun.

## Output

- `results/report.json`: scientific classification and per-seed measurements.
- `results/source_support_admissibility_ledger.csv`: support-leak tests at bookkeeping-mix levels.
- `results/coordinate_covariance_seed_summary.csv`: chart-covariant residuals and negative controls.
- `results/model_constrained_source_covariance_arrays.npz`: 257-sample traces and 513-sample original timeline.
- `results/physical_source_admissibility_and_covariance.png`: comparisons between source support, dynamic coupling and chart artifacts.
- `results/SCIENTIFIC_REPORT_SERVER.md`: concise generated scientific report.
- `SHA256SUMS_OUTPUT.json`: independent output checksums.
- `_server_sendback.zip`: reproducible output bundle, matching root folder.

## Version validity

Synthetic tests and archival results checks are separate. Expected scientific status is **CONDITIONAL_FIXED_CHART_SOURCE_SUPPORT_AND_COORDINATE_COVARIANCE_DIAGNOSTIC**, not a proof of independent physical actuators or material-memory law. All source and document text is in English.
