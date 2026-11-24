# Q550 — v18.52u Phase-Clock Orientation Genealogy & Full-State Conjugacy Audit

Goal: locate the **actual origin** of `dTheta_external/dtau=-x_reference`, compare it with the original 11D `dTheta_internal/dtau=+x`, and test whether the discrepancy is mere relabeling or a distinct prescribed drive.

## Server
From `/root/elab1/hysteresis`:

```bash
unzip q550_phase_clock_orientation_genealogy_and_full_state_conjugacy_v18_52u_bundle.zip
cd q550_phase_clock_orientation_genealogy_and_full_state_conjugacy_v18_52u_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_52u_server.log
```

Default searches original v18.51s sendback in workspace root and native v18.52l extracted results. Optional `--parent-v1851s /absolute/path` and `--states /absolute/path`, or `--require-parent` / `--require-native` to fail if missing. The wrapper **still packages the sendback on errors**.

You should return `..._v18_52u_bundle_server_sendback.zip`. Results: `source_genealogy.json`, `parent_clock_source_excerpt.txt`, `clock_force_discrimination.csv`, `full_field_map_residuals.csv`, `reference_clock_phase_traces.csv`, `report.json`, `SCIENTIFIC_REPORT_SERVER.md`, optional plots.

Development without v18.51s produces a clearly **partial** result; it does not forge the original theta function.

`run_on_server.sh` verifies source checksums and synthetic falsification controls, then reads saved states and ancestor functions; **no IVP integrations**.
