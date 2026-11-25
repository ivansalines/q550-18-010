# Q550 — v18.53v
## Velocity-Occupation Metric Continuation, Conditioned Hessian & Oriented-Valley Branch-Persistence Audit

**Classification before server run: PREPARED / NOT YET A PHYSICAL RESULT.**

This is a reproducible **post-processing audit** of the authenticated v18.53u/v18.53t physical traces. It runs **no new 11D or 21D IVPs**. It varies lambda continuously in two metric families previously introduced by v18.53u. Every scientific assertion is contingent on an authentic server sendback.

## Commands on the Q550 server

```bash
cd /root/elab1/hysteresis
B=q550_velocity_occupation_metric_continuation_conditioned_hessian_oriented_valley_branch_persistence_v18_53v_bundle
unzip -q "$B.zip"
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash "$B/scripts/run_on_server.sh" --require-parent \
  2>&1 | tee v18_53v_server.log
```

To preflight only (no geometrical scans):

```bash
PYTHON_BIN=python3.13 bash "$B/scripts/run_on_server.sh" --preflight-only
```

To smoke-test one authentic physical case with 3 lambda points, without sending or overwriting a scientific report:

```bash
PYTHON_BIN=python3.13 bash "$B/scripts/run_on_server.sh" --smoke
```

## Experiment

Six physical cases (two seeds × three kappa rules), two velocity-occupation calibration epochs, 17 lambda values from 0 to 4 in steps of 0.25: **204 continuation nodes**, plus 141-point oriented-valley comparisons per node (**28,764 curve rows**). Metric-conditioned weights depend on the actual kinetic velocity covariance sampled along the authenticated model; they are **not** unique or fundamental Lagrangian metrics.

At each continuation node:
- Reoptimize the complete nonzero 10D residual minimum in origin and width.
- Evaluate the exact gradient and Hessian of the *archived Cubic Hermite interpolant* and the scaled Hessian eigenvalues and condition number.
- Calculate the implicit derivative of the minimizer with respect to lambda and compare it to centered secants along the continuation family.
- Independently check the numerical minimum with multiple starting points at alternating nodes.
- Evaluate the oriented 45D polygon in the principal SPD square-root identification, with its metric-dependent angle from lambda=0.
- Compare lambda 0, 1 and 3 to all 30 authenticated v18.53u cases.

Original immutable historical source: `inputs/v18_53u_authenticated_server_sendback.zip`. The helper `inputs/u_metric_geometry_authenticated.py` is byte-identical to the parent v18.53u `scripts/metric_geometry.py`, checked against its SHA-256 before use. The v18.53t sendback remains embedded inside v18.53u. All source and result hashes are verified.

## Output

`results/metric_lambda_continuation.csv`, `metric_lambda_valley_curves.csv`, `metric_lambda_family_summary.csv`, `metric_lambda_parent_u_replay.csv`, `metric_lambda_oriented_45d.npz`, two figures, `report.json`, `SCIENTIFIC_REPORT_SERVER.md`, `OUTPUT_MANIFEST.json`, and the full ZIP sendback next to the extracted directory.

## Limits

The canonical kinetic 5D metric remains Euclidean in mass-whitened Cartesian coordinates. The lambda-family is **noncanonical** and its continuity does not establish arbitrary-metric invariance. The Hessian is exact for the source-interpolated finite numerical data, not the continuum physical ODE. Far endpoints are in formal clock continuation beyond xi=25.6; no material closure, physical holonomy, 11D recurrence or memory certificate is asserted.
