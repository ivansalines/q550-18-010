# Q550 — v18.52v Forward–Return Drive Reversal & Signed-Work Attribution Audit

No heavy orbit integrations. Uses the authenticated v18.51s ancestor and the v18.52l native 10D states; compares the **source-defined** return `theta(t)`, `cubic_clock_sign(t, Design)` and `rhs_from_original`, measures the actual archived return-drive energy work, and computes an **instantaneous frozen-state counterfactual**, not a synthetic forward orbit.

```bash
cd /root/elab1/hysteresis
unzip q550_forward_return_signed_drive_orientation_and_energy_work_attribution_v18_52v_bundle.zip
cd q550_forward_return_signed_drive_orientation_and_energy_work_attribution_v18_52v_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_52v_server.log
```

Dependencies: numpy, scipy (staged ancestor may use it), matplotlib optional. `run_on_server.sh` verifies all sources and synthetic identities, then stages the **SHA-256 authenticated** v18.51s genealogy and native v18.52l states. Optional `--parent-v1851s /absolute/file.zip`, `--states /absolute/file.npz`, `--stride 4`, `--require-parent`, `--require-native`. A missing ancestor produces a partial result; no substitute phase is invented. On errors, wrapper still packages any available diagnostics into its sendback. For full server certification, use `--require-parent --require-native` if archives are known present.

Files:
- `authenticated_control_sources.json`: exact source excerpts and ancestor SHA-256.
- `historical_forward_return_clock_protocol.csv`: historical beta-ramp sign and phase-excursion derivative.
- `original_source_signed_rhs_power_probes.csv`: source RHS and power under ± sign at the same state.
- `native_signed_energy_work_ledger.csv`: measured archived physical energy work and frozen-state counterfactual.
- `reference_cumulative_work_curves.csv`, plots, `report.json`, `SCIENTIFIC_REPORT_SERVER.md`.

**Scientific limitation:** native v18.52l is a long **return-branch** record; it contains no matched authentic forward physical trajectory. Do not report the counterfactual as an observed forward/return hysteretic work difference. No claims of overall time reversal, global energy bound, topology or collapse theorem.
