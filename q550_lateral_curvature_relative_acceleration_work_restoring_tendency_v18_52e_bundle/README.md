# Q550 — v18.52e: Lateral Curvature, Relative-Acceleration Work & Restoring-Tendency Audit

The parent v18.52d showed 98.47% lateral acceleration response to a permanent change of kappa0, while eight finite initial velocity kicks moved away in combined phase-space distance. The new experiment resolves the **signed vector work** and **position separation rate** at fixed material time for the same historical nonlinear model; it is deliberately possible for relative velocity to contract while position separation expands.

## Recommended server commands

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_lateral_curvature_relative_acceleration_work_restoring_tendency_v18_52e_bundle_source_only_bundle.zip
cd q550_lateral_curvature_relative_acceleration_work_restoring_tendency_v18_52e_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_52e_server.log
```

The lightweight bundle requires your previously run **exact** `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` at `/root/elab1/hysteresis/` (or in `inputs/`). Do not substitute a reconstructed file with different bytes. The full bundle includes that exact 84 MB archive.

Optional preflight (no nonlinear integrations): `python3.13 scripts/run_audit.py --mode preflight`.

Server stage: rebuild authentic engine from the v18.51s lineage, check each reference history, integrate eight velocity-kick IVPs, evaluate vectors and their acceleration from the exact nonlinear right hand side, interpolate to the same material time, compute the signed kinematic budgets, compare independent derivatives, and package a scientific sendback. The eight IVPs are intentional and belong on your server; no heavy solver run was performed in the chat environment.

## Outputs

`results/report.json`; `results/relative_work_curvature_ledger.csv`; `results/kick_classifications.csv`; `results/relative_work_vs_position_separation.png`; `results/lateral_curvature_and_restoring_acceleration.png`; `results/physical_relative_vector_diagnostics.npz`; `results/SCIENTIFIC_REPORT_SERVER.md`; SHA-256 manifests; ZIP sendback. Preflight outputs are labeled separately.

## Environment

Python 3.13, NumPy, SciPy, Matplotlib. Absolute time/drive phase is explicitly distinguished from elapsed material time. The source code and documents use English only.
