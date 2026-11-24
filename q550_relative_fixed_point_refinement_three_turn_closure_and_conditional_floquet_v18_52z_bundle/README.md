# Q550 — v18.52z: Relative Fixed-Point Refinement, Three-Turn Closure & Conditional Floquet Audit

## Purpose

Continue the exact `u` candidates preserved in the authentic v18.52y sendback, which finished with an ~1.5e-5 **uncertified** relative fixed-point residual after only three Newton updates. Do not repeat the original search from scratch and do not substitute the mixed-charge historical kicks.

This audit asks (i) whether source-authenticated fixed-J Newton continuation actually closes the 8D quotient state; (ii) whether independently integrated 1T/2T/3T trajectories confirm closure; (iii) whether the fixed-J quotient return map has a numerically converged spectrum that may be called conditional Floquet multipliers; and (iv) if so, how drive work and dissipation balance over the closed relative cycle.

## Directory contents

* `scripts/run_audit.py`: new authentic server-only IVP/root/spectral/energy audit, source SHA checks, conditional decisions.
* `scripts/return_map.py`: exact fixed-J lift, 8D quotient, wrapped differential, finite-difference 8-column Jacobian; source retained verbatim from v18.52y.
* `scripts/energy_math.py`, `scripts/phase_energy.py`: original/source-consistent energy and signed-work formulas.
* `reference/v18_52y_authenticated_report.json`: exact source for both candidate vectors, not manually transcribed.
* `reference/*.py`: historical executable Python source, including the v18.52y server flow, native engine parent staging and original physical RHS ancestor.
* `tests/test_synthetic.py`: quick falsification tests, **not** authentic Q550 outcomes.
* `inputs_config.json`: source and native SHA-256, exact numerical acceptance tolerances.
* `MATHEMATICAL_PROTOCOL.md`, `SCIENTIFIC_DECISION_RULES.md`, `SOURCE_PROVENANCE.md`, `requirements.txt`: English scientific and executable genealogy.
* `SOURCE_MANIFEST.json`, `OUTPUT_MANIFEST.json`: SHA-256 verification of inputs/sources and resulting sendback outputs.

## Server command

```bash
cd /root/elab1/hysteresis
unzip q550_relative_fixed_point_refinement_three_turn_closure_and_conditional_floquet_v18_52z_bundle.zip
cd q550_relative_fixed_point_refinement_three_turn_closure_and_conditional_floquet_v18_52z_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh --require-parent --require-native 2>&1 | tee ../v18_52z_server.log
```

For partial/debug use `--seeds original`, or `--max-newton 7`, but retain all source SHA controls. No GUI or user editing is required.

**External inputs expected in parent workspace**:

- `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` SHA `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`.
- `q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz` SHA `92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12`.

The authentic phase comes from the parent `ctx['theta']`, not the 11th native NPZ column (which is tau). The fixed native kappa is `0.1666331385421536` in BOTH families.

## Outputs

Server writes `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`, optionally two `*_three_turn_fixed_charge_closure.png` figures, plus `OUTPUT_MANIFEST.json` and the ZIP `q550_relative_fixed_point_refinement_three_turn_closure_and_conditional_floquet_v18_52z_bundle_server_sendback.zip` adjacent to the extracted root. If a required ancestor is missing, the strict run exits nonzero and still packages an explicit partial/failure report.

## Scientific caveats

No user has yet produced a v18.52z authentic return closure; local tests validate **only** software mechanics. A periodic fixed point in the late externally driven return section modulo counterrotation is not proof of a globally conservative limit cycle or a full historical forward-return orbit. The carica J is exactly conserved while the total energy can still be driven and dissipated. See `SCIENTIFIC_DECISION_RULES.md`.
