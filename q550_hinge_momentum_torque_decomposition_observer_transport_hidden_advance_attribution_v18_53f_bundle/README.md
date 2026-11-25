# Q550 v18.53f — Hinge-Momentum Torque Decomposition, Observer Transport & Hidden-Advance Attribution Audit

## Purpose

Distinguish genuine kinetic-plane hinge torque, model-source force attribution, and apparent torque due solely to a rotating observer, using the authenticated v18.52l absolute 10D states and relative 5D traces. Preserve the v18.53e fixed endpoint selection (24 pairs) and observer genealogy.

**Scientific boundary:** the flat-kinetic `J_hinge^(T)` is NOT an established Noether invariant of the full driven Q550. `J1-J2` is the independently certified charge. Do not promote a minimized projected endpoint gap to an absolute orbit closure.

## Execution (server)

```bash
cd /root/elab1/hysteresis
unzip q550_hinge_momentum_torque_decomposition_observer_transport_hidden_advance_attribution_v18_53f_bundle.zip
cd q550_hinge_momentum_torque_decomposition_observer_transport_hidden_advance_attribution_v18_53f_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  --require-parent --require-native --require-trace \
  2>&1 | tee ../v18_53f_server.log
```

The script locates the SHA-pinned v18.51s sendback and v18.52l NPZ files in the normal hysteresis workspace. You can override them with `--parent-v1851s PATH --native-npz PATH --trace-npz PATH` and `--base PATH`.

**No new Q550 orbit integration.** The historic `stage_engine` source reconstruction may perform a short inherited context-generation step; the native 10D states and relative 5D trajectories are not recomputed.

## Inputs

- Authenticated `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` (expected SHA-256 in config). Mandatory for full source torque attribution.
- Native `extended_horizon_10d_states.npz` (20 series, full original length 65,537 samples each, subsampled to stride 64). The 11th column is material time `tau`, **never** the external phase.
- Relative `extended_horizon_5d_vectors_and_traces.npz` (8 series), subsampled with stride 16, with 24 fixed-endpoint pairs inherited and SHA-pinned from v18.53e.
- In-archive actual historical source genealogy: v18.51t loader, historical support physical RHS, v18.52v-r1 source audit, v18.53e source audit, v18.49h kinetic embedding and v18.52n fixed observer basis.

All original source references are retained. Large NPZ data are **not duplicated** in the ZIP: they are external SHA-pinned dependencies.

## Scientific outputs

- `results/authenticated_original_source_metadata.json`: verified parent identity, source parameters, fixed native kappa and phase genealogy.
- `results/authenticated_rhs_verification.csv`: agreement between source-exact 11D RHS (material-time normalized) and the independent six-force reconstruction.
- `results/native_absolute_hinge_torque_summary.json`: per-series force-channel time integrals, static and rotating-plane torque identity residuals, charge drift.
- `results/native_hinge_torque_and_observer_transport.csv`: phase-space source torque and observer torque at 201 shown times per native series.
- `results/fixed_endpoint_hidden_advance_force_vs_observer.csv`: 24 physical endpoint pairs × 5 **preselected** observer angles; projected and hidden 5D advance, velocity-integral check.
- `results/relative_moving_observer_hidden_ledger.csv`: the signed hidden coordinate change split into physical velocity integral and frame-transport integral, using only relative separation vectors.
- `results/original_hinge_torque_observer_decomposition.png`, `results/fixed_endpoint_hidden_advance.png` (if matplotlib installed).
- `results/SCIENTIFIC_REPORT_SERVER.md`, `results/report.json`, `OUTPUT_MANIFEST.json`, and sendback ZIP.

## Independent checks

1. Source integrity and strict parent authentication; the v18.53e output and v18.52n frame are independently checked by SHA.
2. Native 5D embedding cross-checked against historical `kinetic_metric_v18_49h.py` and flat kinetic metric `P=V`.
3. Source RHS checked at native samples against independent coordinate equations and 5D force decomposition.
4. `J1-J2` charge conserved while fixed-plane kinetic hinge moment can vary.
5. Within-plane observer basis rotation contributes exactly **zero** to change of oriented moment. Tilting the *plane itself* produces nonzero observer transport.
6. Both static and moving-plane moment changes compared with integrated torques; changing the observer is not treated as a physical interaction.
7. Preselected fixed endpoint pairs compared in all observer angles, and hidden 5D advance is checked by integrating the stored relative velocities.
8. Synthetic tests falsify the false assertion that a closed 2D helical shadow is an absolute orbit closure.

If any required source is missing, the strict run exits nonzero **while packaging an honest partial diagnostic sendback**. It never invents an original phase. The optional local run without `--require-parent` is explicitly classified as partial.

## Install

Python 3.13 preferred; `numpy`, `scipy` and optional `matplotlib`. Standard library only otherwise. No required `sympy` for this audit.
