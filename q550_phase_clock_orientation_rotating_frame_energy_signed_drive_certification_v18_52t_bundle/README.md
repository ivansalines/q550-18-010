# Q550 — v18.52t: Phase-Clock Orientation, Rotating-Frame Energy & Signed-Drive Certification Audit

**Goal:** settle the measured phase sign without changing native 10D states or historical parameters. Inspect an apparent conflict between the v18.52l externally prescribed phase (`dTheta/dtau` near `-x`) and the historical full 11D generator (`dTheta/dt=+x`). Test both signed ordinary-energy and rotating-frame budgets, including the original engine's 11th-component clock law.

**Run on the original Linux server** (Python 3.13 preferred):

```bash
cd /root/elab1/hysteresis
unzip q550_phase_clock_orientation_rotating_frame_energy_signed_drive_certification_v18_52t_bundle.zip
cd q550_phase_clock_orientation_rotating_frame_energy_signed_drive_certification_v18_52t_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_52t_server.log
```

The script finds the native v18.52l results and the authenticated v18.51s
sendback in `/root/elab1/hysteresis`. Optional overrides:

```bash
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  --states /absolute/path/to/extended_horizon_10d_states.npz \
  --parent-v1851s /absolute/path/to/v18_51s_server_sendback.zip
```

`--stride 1` is the default (full saved resolution); use `--stride 16` for a
preflight only. If an ancestor is unavailable, the run produces an **incomplete**
report instead of inventing that phase. A mismatched hash stops execution.

The output contains `results/report.json`, `SCIENTIFIC_REPORT_SERVER.md`,
`symbolic_signed_identity.json`, signed energy ledgers per native series, and
figures for the original and autonomous reference series. A self-contained
`*_v18_52t_bundle_server_sendback.zip` is packaged at the parent directory.

Python dependencies: numpy, scipy, matplotlib (plots), sympy (optional exact
symbolic verification). Original `v18.51s` lineage loader may have additional
requirements already present from prior Q550 audits.

**Do not interpret a successful ten-state energy ledger as proof that the full
11D original model uses the same clock orientation as the externally supplied
v18.52l phase.** See `SCIENTIFIC_DECISION_RULES.md`.
