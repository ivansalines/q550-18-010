# Q550 — v18.52a: Dynamically Admissible Cross-Channel Turning Robustness Audit

**Purpose.** Test whether the positive v2–omega2 mixed turning term, observed in v18.51z, persists when the *historical nonlinear model parameter* kappa0 is varied across the original/autonomous range. This is not another passive observer or algebraic source-recombination experiment.

## Scientific protocol

- Three already integrated, independent nonlinear kappa0 trajectories from v18.51r are used for **each** of the two physical seed families.
- On each trajectory, the structural 10D field Jacobian and partial derivative of the field with respect to kappa0 are recomputed from the historical nonlinear engine, at fixed state and physical drive phase.
- The source channels are strictly the model-defined velocity slots `v2 = state[5]`, `omega2 = state[7]`, plus an *unaltered* `other` residual.
- Three independent, nonhomogeneous forced tangent response channels are integrated and recomposed with the original total sensitivity.
- Signed angular turning terms are monitored on all 487 qualified sections (`xi >= 0.01`) through `xi=0.20`, on both seed families.
- The midpoint response is checked against the independent v18.51z server sendback; a coarser Jacobian sampling grid and optional separate Green propagator check give numerical controls.
- A finding that the physical dominance fails is a **valid scientific result**, not a numerical error; the report distinguishes these outcomes.

## Download variants

**Light archive:** contains this source tree, the authenticated Z sendback and the SHA-256 of the historical S sendback. The user must have `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` at `/root/elab1/hysteresis/` or place it into `inputs/` before running.

**Full archive:** includes both authenticated Z and S sendbacks, along with all scripts/tests/reference implementation; it has no dependency on an external copied archive.

## On the user's research server

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_dynamically_admissible_cross_channel_turning_robustness_v18_52a_bundle.zip
cd q550_dynamically_admissible_cross_channel_turning_robustness_v18_52a_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_52a_server.log
```

For the light archive, substitute the light ZIP filename; make sure the S parent is in the parent workspace first.

Optional faster branch survey: `python3.13 scripts/run_audit.py --skip-green` skips the independent fundamental-matrix check at the midpoint only; its absence is documented by the output checks.

## Outputs

- `results/physical_parameter_robustness_sections.csv`: signed turning budget and geometric invariants, all seeds/branches/qualified sections.
- `results/physical_branch_turning_arrays.npz`: numerical channel response/turning arrays.
- `results/physical_kappa_branch_cross_channel_turning_robustness.png`: visual comparison.
- `results/report.json`: numerical validity and **separate scientific robustness decision**.
- `results/SCIENTIFIC_REPORT_SERVER.md`: short interpretation and limitations.
- `SHA256SUMS_OUTPUT.json`: output SHA-256 integrity report.
- `q550_dynamically_admissible_cross_channel_turning_robustness_v18_52a_bundle_server_sendback.zip` in the parent workspace.

## Provenance and limitations

No synthetic physics, tunable source-channel mixing or off-orbit forcing is used in the real Q550 audit. The only perturbed dynamical parameter is the historical kappa0. Reusing six existing nonlinear IVP branches avoids re-running their more expensive integrations; **fresh structural jets and forced tangent IVPs still run on the user's server**.

Success within this narrow parameter band does **not** imply robustness to arbitrary admissible changes in the generator, material metric, initial data or history, and does not certify physical return or a material memory constitutive law. The algebraic turning-budget identity is not independent physics; its signed dominance and stability across original model branches are the empirical questions.
