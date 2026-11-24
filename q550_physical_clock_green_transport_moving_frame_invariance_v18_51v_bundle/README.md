# Q550 — v18.51v: Physical-Clock Green Transport & Moving-Frame Invariance Audit

## Purpose

This audit continues the authenticated v18.51u two-time Green-kernel genealogy. It distinguishes:

1. A **frozen-history** change of independent variable `xi -> tau`, which changes a Green source density by `1/c(tau_source)` while preserving the same integrated response.
2. The **physical-clock sensitivity at fixed tau** when the clock itself depends on kappa and the physical trajectory. This requires an explicit clock-shift correction and direct matching of the two distinct historical kappa trajectories.
3. Passive changes of the moving second-plane observer frame. Full Cartesian vectors and their Euclidean norms are invariant; individual radial/azimuthal components, source cancellation percentages, and plain derivatives generally are not. Covariant moving-frame derivatives include the connection term.

It does **not** certify a new material constitutive memory, thermodynamic dissipation, full return closure, globally invariant source-channel attribution, or an independent physical actuator.

## Inputs and provenance

`inputs/q550_time_resolved_green_kernel_compensation_birth_v18_51u_bundle_server_sendback.zip` contains the full genealogical chain U -> T -> S -> R. The small `authenticated_*_v18_51t.npz` and `authenticated_*_v18_51r.npz` caches are compared **byte-for-byte** against the corresponding sources extracted from this nested ancestry during the server preflight. The source manifest is checked before any numerical computation.

The archived v18.51r nonlinear trajectories already contain `xi`, `midpoint_Y`, `midpoint_Q`, physical `clock` and `tau`, `Qprime`, and all original/autonomous kappa branches for the original and autonomous seed families. The v18.51t archive contains independent tangent sensitivity and channel-resolved responses. The v18.51u archive contains the two-time Green kernels.

The audit does **not** reload or reintegrate the original expensive nonlinear 10D solver. It processes the saved dynamics and source histories; the two-time arrays and nonuniform-time reconstructions are best run on the research server.

## Run on the server

```bash
cd /root/elab1/hysteresis
unzip /root/Downloads/q550_physical_clock_green_transport_moving_frame_invariance_v18_51v_bundle.zip
cd q550_physical_clock_green_transport_moving_frame_invariance_v18_51v_bundle

set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_51v_server.log
```

`run_on_server.sh` verifies the source manifest, runs a fast analytic synthetic suite, verifies the nested provenance, builds the full physical-clock and observer-frame data, runs independent output checks, writes SHA256 output manifests, and creates:

`../q550_physical_clock_green_transport_moving_frame_invariance_v18_51v_bundle_server_sendback.zip`

## Light package

The source-only ZIP contains all scripts, protocols, tests, authenticated cached R/T arrays, and preliminary metadata, but **not** the ~103 MB v18.51u sendback. Place a copy of the already-available U server sendback under `inputs/` with its exact historical filename before running. The expected SHA256 is recorded in `SOURCE_MANIFEST.json`; a mismatch stops the audit before analysis.

## Output files

- `results/report.json`: actual numerical result and scientific checks.
- `results/SCIENTIFIC_REPORT_SERVER.md`: generated scientific summary and limits.
- `results/physical_clock_frame_section_ledger.csv`: two-seed section-by-section record including physical time, channel responses and clock correction.
- `results/physical_clock_frame_green_tensors.npz`: clock-density kernels, reconstructed physical vector responses, moving-frame fields and fixed-tau sensitivity comparisons.
- `results/physical_clock_vs_fixed_tau_sensitivity.png`: same-physical-time correction.
- `results/moving_frame_components_vs_physical_vectors.png`: frame-dependent components versus invariant plane-vector norm.
- `results/physical_clock_causal_green_density.png`: the causal two-time kernel expressed in units of physical source time.
- `SHA256SUMS_OUTPUT.json`: hashes of actual server output files.

The pre-distributed `SCIENTIFIC_REPORT.md` is only a protocol. The real scientific report is generated during the server run; do not mistake one for the other.

## Dependencies

Python 3.11+, NumPy, SciPy and Matplotlib. Python 3.13 is recommended on the established server.
