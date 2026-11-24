# Q550 — v18.52x: Fixed-Charge Relative Poincaré Return Map & Transverse Identifiability Audit

This is a **reproducible, no-new-IVP** audit of the authentic v18.52l native long return trajectories in the authentic v18.51s phase clock. All script comments and reports are in English.

## Core question

Does the physical 10D state return after one external-phase winding, **modulo counterrotation** and on its conserved angular-charge leaf? Can the available perturbation genealogy identify the transverse derivative of a return map on that leaf? A limit cycle, an attracting cycle, or a Floquet spectrum must never be asserted from a phase projection alone.

## Findings anticipated from verified input structure, not assumed outputs

- The 20 long native series comprise 2 seeds × (reference, four signed kick amplitudes) × (base and refined). **Refined traces are numerical controls, not extra independent perturbation directions.**
- Historical kicks change the exact angular charge `J = exp(2 logr1)*omega1 - exp(2 logr2)*omega2`. Each trajectory conserves its **own** J. Mixing kick trajectories therefore samples different invariant leaves.
- The horizon xi=[12.8,25.6] covers approximately 1.7 external-phase turns. It supports **one section return**, not two iterates at the same section.
- Four signed kick data points per seed cannot span an 8D fixed-J transverse tangent space; the fixed-J Floquet spectrum is structurally **not identifiable** from these native kicks.

## What is computed

1. **Authenticated ancestry:** SHA-256-pin the original v18.51s sendback (including nested v18.51r) and v18.52l native 10D NPZ; reconstruct authentic source clock `ctx['theta']` through source loader, not from the 11th native column (which is material time, not phase).
2. **Late signed phase:** certify that `cubic_clock_sign == -1` on the sampled late branch and that the reconstructed phase is strictly decreasing. This makes the section at external phase modulo `2*pi` well-defined on the sampled branch.
3. **Physical quotient chart:** reduce 10D physical state by counterrotation. At fixed `J` the eight chart coordinates are `(logr1,v1,omega1,logr2,v2,z,vz,alpha)` with `alpha = (theta1+theta2) mod 2*pi`, and `omega2` reconstructed via the exact charge. This is a chart, not a canonical Euclidean metric.
4. **Section crossings:** for each of 17 phase anchors per seed, find start and end phase-separated by `-2*pi`, interpolate all genuine native traces and compare all eight reduced coordinates. Report channel-wise discrepancies and a **descriptive** norm. All 20 × 17 = 340 observed pairs should be available.
5. **Charge stratification:** for every trace record its initial/final/max-drift charge. For each kick compare its charge with the same-seed reference on the same section. **Never project a kick to another charge leaf and claim its observed endpoint is physically propagated there.**
6. **Input-rank and identifiability:** compute singular values of four native kick directions and mixed-charge central secant gains at ±1e-4, ±3e-4. These are **one-parameter family directional responses**, explicitly *not* fixed-charge Floquet multipliers.
7. **Quality controls:** verify 10 refined/base pairs, compare local quadratic and linear interpolation, test gauge invariance, hidden-posture false returns, and a charge-mixed false stability certificate.

## Execution on the user's server

```bash
cd /root/elab1/hysteresis
unzip q550_fixed_charge_relative_poincare_return_map_and_transverse_identifiability_v18_52x_bundle.zip
cd q550_fixed_charge_relative_poincare_return_map_and_transverse_identifiability_v18_52x_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh --require-parent --require-native \
  2>&1 | tee ../v18_52x_server.log
```

Default authentic inputs expected in `/root/elab1/hysteresis`:

- `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` (source-sign and original clock ancestry; SHA pinned in `inputs_config.json`)
- `q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz` (native 20×65,537×11, last column **tau**; SHA pinned)

If stored elsewhere, pass `--parent-v1851s /full/path/...zip --native-v1852l /full/path/...npz`. No source paths or result metrics are fit to obtain a desired answer. `run_on_server.sh` always packages `..._server_sendback.zip` after a successfully completed set of preflight tests, even on authentic audit failure.

## Outputs

- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`
- `results/one_step_fixed_phase_return_pairs.csv`
- `results/per_trajectory_charge_conservation.csv`
- `results/section_input_rank_and_identifiability.json`
- `results/mixed_charge_kick_family_directional_response.csv`
- `results/refined_vs_base_return_precision.csv`, `results/section_interpolation_precision.csv`
- `results/original_relative_return_map.png`, `results/autonomous_relative_return_map.png` if matplotlib is available
- `OUTPUT_MANIFEST.json`, and a self-contained server sendback ZIP

## Interpretation

Even if the first-return residual is tiny, **one return is not a stable limit cycle**. Even if a kick perturbation appears to contract, it does **not** imply transverse attraction *on a fixed-J leaf*, because the historical perturbations change J and span fewer than eight independent tangent directions. A future complete Floquet diagnosis requires independent, model-consistent same-J tangent vectors and genuine repeated return/monodromy measurements, not a forced periodicity inferred from one projection.
