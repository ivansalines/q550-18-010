# Q550 — v18.52w: Drive-Zero Crossing, Posture Persistence & Relative-Cycle Discrimination

**Purpose:** Test whether the physical Q550 10D state continues to move as the historical external phase drive reverses, and determine whether long-horizon returns after one external phase winding are genuine full-state recurrences modulo counterrotation. This is **not** an attracting limit-cycle proof or a Hamiltonian conservation theorem.

## Authenticated input genealogy

1. `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` — historical source tree, original `make_context`, `theta`, `cubic_clock_sign`, and **nested v18.51r midpoint 10D early return states** for `xi in [0,0.2]`. SHA-pinned by `inputs_config.json`. No fake forward branch is generated.
2. `q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz` — 20 full-length 10D return series; `common_xi in [12.8,25.6]`, 65,537 samples per series. SHA-pinned by `inputs_config.json`. **This archive does NOT include the signed drive's zero crossing.** Column 10 in this archive is accumulated material time `tau`, **not** the phase `Theta`.
3. `reference/v18_52v_r1_run_audit.py`, `reference/parent_v18_51t_run_audit.py`, `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py`: authentic functional source lineage used for re-staging and verification.
4. Reference summaries v18.52v-r1 and v18.52j included with their retained SHA-256 manifest.

## Run on server

```bash
cd /root/elab1/hysteresis
unzip q550_drive_zero_crossing_posture_persistence_and_relative_cycle_discrimination_v18_52w_bundle.zip
cd q550_drive_zero_crossing_posture_persistence_and_relative_cycle_discrimination_v18_52w_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh --require-parent --require-native 2>&1 | tee ../v18_52w_server.log
```

Optional explicit paths:

```bash
python3.13 scripts/run_audit.py \
  --parent-v1851s ../q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip \
  --native-v1852l ../q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz \
  --require-parent --require-native
```

`run_on_server.sh` verifies the source manifest, executes synthetic tests, runs the audit, and **packages a sendback even after a failure** (with `results/FAILURE_DIAGNOSTIC.txt`). It does not integrate new native Q550 orbits; it re-stages the historical short ancestor environment to inspect original sources.

## Outputs

- `results/actual_drive_zero_posture_events.csv`: physical 10D and Cartesian 5D speed, tangent direction before/after zero, counterrotation charge, drive work, energy balance across windows.
- `results/drive_zero_local_curves.csv`: phase command, material speed, independent angular/gauge rates, work channels around source zero.
- `results/one_turn_full_quotient_recurrence_scout.csv`: one external-phase turn, source/target support time, all 8 non-angle state coordinates plus angle-sum mismatch (mod 2pi), charge consistency; no projected-only surrogate.
- `results/authenticated_source_genealogy.json`: recovered executable source for clock, phase, original drive and parent hash.
- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, up to four plots and `OUTPUT_MANIFEST.json`.

**Interpretation guardrails:** A nonzero 5D tangent at drive zero establishes only local persistence. A one-phase-turn approximate match modulo counterrotation is not proof of a periodic or attracting relative orbit. Even an exact periodic orbit can have a conserved angular Noether charge in a driven/dissipative dynamical system without globally conserved energy. Full limit-cycle certification requires repeated Poincare return, independent stability diagnostics, and a clear autonomous or periodically driven setting.
