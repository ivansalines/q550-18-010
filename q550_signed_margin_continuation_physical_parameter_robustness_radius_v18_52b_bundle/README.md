# Q550 — v18.52b

**Signed-Margin Continuation & Physical Parameter Robustness Radius Audit**

This is the next executable audit after v18.52a. Its job is to search for a *boundary*, not manufacture a robustness certificate.

## Source conventions

- `inputs/q550_dynamically_admissible_cross_channel_turning_robustness_v18_52a_bundle_server_sendback.zip` is included and hash-authenticated.
- The lighter edition requires the complete, original `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` in `inputs/` or the parent `/root/elab1/hysteresis` workspace. This is the exact historic input used by v18.52a.
- The full edition also includes that large v18.51s parent within `inputs/`. The two editions run the **same scripts**.
- No original model source is replaced; the audited parent 52a script and the deeper v18.51t → v18.51s → v18.51r → v18.51q ancestry are staged in a temporary directory.

## Execution (research server)

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_signed_margin_continuation_physical_parameter_robustness_radius_v18_52b_bundle.zip
cd q550_signed_margin_continuation_physical_parameter_robustness_radius_v18_52b_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_52b_server.log
```

To extend the search by explicitly choosing the normalized radii, use the audit command directly after verifying the source manifest:

```bash
OPENBLAS_NUM_THREADS=1 python3.13 scripts/run_audit.py \
  --levels 2 4 8 16 32 64 --bisect-steps 5
```

Larger parameter ranges are not guaranteed to be dynamically admissible. The default predeclared scan is deliberately modest.

## Outputs

- `results/physical_parameter_continuation_scan.csv`: all anchor/new branches with status and resolution evidence.
- `results/physical_parameter_continuation_curves.npz`: full section-wise responses and turning budgets.
- `results/physical_kappa_signed_margin_continuation.png`: positive-margin scan and failure-region visualization.
- `results/report.json`: separate numeric audit status, scientific verdict, per-side sampled lower radii, failure classification, and continuation brackets.
- `results/SCIENTIFIC_REPORT_SERVER.md`: interpretation and limitations after the physical server calculation.
- `SHA256SUMS_OUTPUT.json`: output integrity.

## Interpretation

**`SAMPLED_STABILITY_UP_TO_TESTED_LIMIT__NO_RADIUS_BOUNDARY_FOUND`** means no failure at sampled values. It is not a certified parameter interval. **`SAMPLED_SIGN_FAILURE_BRACKET_FOUND`** means a candidate crossing has been bracketed by truly integrated paths, but neither interpolation nor bisection makes it a rigorous proof of a physical bifurcation. **`INCOMPLETE_DUE_TO_SOLVER_LIMIT`** must never be rebranded as physical loss of dominance.

## Local limits

Only source verification, historic-band preflight and synthetic unit tests are run locally. Nonlinear and tangent continuation is exclusively entrusted to the research server.
