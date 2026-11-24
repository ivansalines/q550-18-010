# Q550 — v18.52d
## Perturbation-to-Lateral-Acceleration Transfer & Recovery Audit

An executable, two-layer, source-grounded audit of historical Q550 nonlinear 10D dynamics.

**Layer 1: Lightweight preflight** uses the exact v18.51r 513-section Cartesian states, velocities and integrated material clocks for both seed families and two historical kappa0 choices. It derives a kinematic finite-parameter acceleration difference by differentiating archived physical velocities with respect to tau. This step does not integrate the nonlinear model.

**Layer 2: Server physics** reconstructs the historical 19-ancestor 10D generator from the SHA256-verified v18.51s sendback, computes physical accelerations directly from the nonlinear RHS, and executes paired signed transverse initial-velocity kicks at a fixed kappa0. These integrations are intentionally reserved for the user's server.

Inputs: `inputs/historical_v18_51r_physical_states.npz` (byte-identical archive member); the actual v18.52c server sendback (24 physical parameter probes); source files `reference/kinetic_metric_v18_49h.py` and `reference/parent_v18_51t_run_audit.py`, taken from authenticated ancestors.

## Local preflight (no heavy solver)

```bash
PYTHON_BIN=python3.13
$PYTHON_BIN scripts/verify_source_manifest.py
$PYTHON_BIN tests/test_synthetic.py
$PYTHON_BIN scripts/run_audit.py --mode preflight
```

## Full run on server

Put the *exact* original v18.51s server sendback in `/root/elab1/hysteresis/` or `inputs/`. Its name is:

`q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`

```bash
cd /root/elab1/hysteresis/q550_perturbation_to_lateral_acceleration_transfer_recovery_v18_52d_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
   2>&1 | tee ../v18_52d_server.log
```

The runner authenticates sources and the v18.52c result, tests synthetic projection/return cases, recovers the historical generator, evaluates physical acceleration on all six archived nonlinear branches, independently integrates two baselines and 8 signed transverse perturbations, compares material-time-matched trajectories, writes results, verifies the scientific numerical checks, hashes outputs, then creates the server sendback ZIP.

Options: `python3.13 scripts/run_audit.py --mode server --skip-kicks` is a debug mode (not the default scientific experiment), and `--kick-relative-amplitudes 1e-4 3e-4` sets the signed amplitude families.

## Expected results

- `LOCAL_PREFLIGHT.json` — no-generator interpolation diagnostic
- `fixed_tau_lateral_acceleration_ledger.csv` — 2 seeds x 513 material-time comparisons
- `transverse_kick_recovery_ledger.csv` — signed paired transient tests
- `physical_kick_trajectories_and_diagnostics.npz` — numerical time series from kick runs
- `lateral_acceleration_vs_recovery.png` — acceleration and separation chart
- `report.json`, `SCIENTIFIC_REPORT_SERVER.md`, `SHA256SUMS_OUTPUT.json` — provenance, limits and checks

All code and scientific document text is in English; version and full title remain embedded in the ZIP root and output names.

**No claim of asymptotic/material recovery is made from a finite horizon.**
