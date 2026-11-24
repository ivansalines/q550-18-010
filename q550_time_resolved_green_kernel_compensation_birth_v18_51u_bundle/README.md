# Q550 — v18.51u: Time-Resolved Green Kernel & Compensation-Birth Audit

This is the next reproducible numerical **caliper**, built on the fully authenticated `v18.51t` server sendback. It asks when a *source-time* perturbation contributes to a *later observation*, when radial opposition becomes detectable, and which past windows generate the measured signed cancellation.

**This bundle does not rerun the historical nonlinear 10D Q550 solver.** The expensive nonlinear trajectory, original analytic state Jacobian, and parameter-source samples are already archived and audited in the nested v18.51t genealogy. The new work integrates only the stored tangent propagator and several independent small impulse-restart linear IVPs, then constructs the time-resolved two-time kernels and their projected integrals. All heavy audit calculations are intended to run on the research server.

## Full archive: server commands

```bash
cd /root/elab1/hysteresis
unzip q550_time_resolved_green_kernel_compensation_birth_v18_51u_bundle.zip
cd q550_time_resolved_green_kernel_compensation_birth_v18_51u_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_51u_server.log
```

The script creates an adjacent `q550_time_resolved_green_kernel_compensation_birth_v18_51u_bundle_server_sendback.zip` with the scripts, verified 19-ancestor input chain through t, results, report, and SHA-256 manifests.

## Lightweight archive (if the main download misbehaves)

The `..._source_only_bundle.zip` retains the identical scientific root folder and all scripts/tests/docs, but **omits the 86 MB already-delivered v18.51t server sendback**. Put the following existing file in the extracted root's `inputs/` directory before running:

`q550_source_to_response_green_transport_radial_azimuthal_causality_v18_51t_bundle_server_sendback.zip`

For example, from the extracted root directory:

```bash
cp ../q550_source_to_response_green_transport_radial_azimuthal_causality_v18_51t_bundle_server_sendback.zip inputs/
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_51u_server.log
```

The source manifest strictly rejects a missing or modified ancestor. Nothing is silently fetched or substituted.

## Separate execution modes

- `python3.13 tests/test_green_kernel_synthetic.py`: tiny, exactly soluble causal-kernel benchmark; no physical Q550 integration.
- `python3.13 scripts/run_audit.py --mode preflight`: verify nested genealogy, source arrays, and scientific parent checks. **This does not calculate a new kernel.**
- `bash run_on_server.sh`: perform the full v18.51u kernel audit, validate all tests, write generated report/output hashes and pack sendback.

## Expected outputs **after the full server run**

- `results/green_causal_two_time_kernels.npz`: two-seed, 513×513 triangular causal kernels for radial/azimuthal projections and `v2`, `omega2`, `other`, total source channels.
- `results/compensation_birth_section_ledger.csv`: 1,026 observation-section records of signed cumulative contributions and radial cancellation.
- `results/causal_source_window_contributions.csv`: non-overlapping source-time windows (0.025 event-time width) whose signed Simpson sums recompose the observations.
- `results/green_kernel_causal_map_original.png` and `..._autonomous.png`: four-panel causal time maps with source time horizontal and observation time vertical.
- `results/compensation_birth_and_channel_history.png`: accumulated radial response and cancellation thresholds.
- `results/report.json`: quantitative error budgets, channel reconstruction, independently repropagated impulse probes, window/refinement diagnostics, source-age measures, operational birth markers and all checks.
- `results/SCIENTIFIC_REPORT_SERVER.md`: English scientific report **written after** measurements, not a preliminary report mislabeled as completed.

## Main scientific distinction

The two-time kernel `K(xi,u)` is a model-conditioned, moving-frame, first-variation response. A negative kernel sector, a cancellation threshold or a response-weighted source age is **not yet** an intrinsic physical memory mechanism, unique causality, material torque, or completed 10D return. The `other` channel preserves first-oscillator and any additional original parameter-force components. The original model and mathematical context are unchanged.

Requirements: Python 3.13, NumPy 1.26+, SciPy 1.11+, Matplotlib 3.8+. Source and report content, as well as all code comments, are in English.
