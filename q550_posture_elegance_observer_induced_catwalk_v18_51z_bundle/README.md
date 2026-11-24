# Q550 — v18.51z: Posture, Elegance & Observer-Induced Catwalk Audit

A reproducible scientific audit of **signed angular transport of parameter-sensitivity response vectors** in the second Cartesian plane. The title is playful; the mathematical objective is to separate source-response geometry from passive observer effects.

## Run

```bash
cd /root/elab1/hysteresis
unzip q550_posture_elegance_observer_induced_catwalk_v18_51z_bundle.zip
cd q550_posture_elegance_observer_induced_catwalk_v18_51z_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_51z_server.log
```

The runner verifies sources, runs synthetic controls, replays both seed families using the y sendback, runs independent output checks, records SHA256 hashes, and creates a sendback ZIP next to the bundle folder.

## Included genealogy

- `inputs/q550_model_constrained_source_subspaces_physical_attribution_admissibility_v18_51y_bundle_server_sendback.zip` — unmodified v18.51y server result; contains authenticated T, W, X inputs, original scripts, and the Y report.
- `SOURCE_MANIFEST.json` — SHA256 of all inputs, scripts, tests and interpretive documents.
- `scripts/run_audit.py` — actual replay, tables, figures and scientific checks.
- `tests/test_synthetic.py` — independent synthetic identities and negative controls.
- `tests/test_results.py` — independent checks of generated files.
- `results/*` — local replay results for reference.
- `SHA256SUMS_OUTPUT.json` — output hashes; recreated by the server.

No expensive nonlinear solver or Green integral is run in v18.51z; parent model-derived arrays are replayed. The local result is **not** a replacement for independent server reproducibility.

## Primary outputs

- `results/posture_observer_turning_arrays.npz`
- `results/turning_budget_sections.csv`
- `results/report.json`
- `results/SCIENTIFIC_REPORT_SERVER.md`
- `results/channel_posture_and_observer_catwalk.png`
- `results/physical_clock_posture_turning.png`

A positive-clock reparametrization and a common passive observer rotation are transformations of the **same model-derived response**, not new physics. The origin has zero response, hence undefined angle.
