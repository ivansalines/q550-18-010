# Q550 — v18.52y: Fixed-Charge Tangent Basis, Relative-Closure Search & Return-Map Stability Audit

This server-run bundle continues authentic Q550 v18.52x. Unlike v18.52x, it **launches new authentic 10D IVPs on the user's server** and forms an eight-column Jacobian using eight *independent tangent vectors on the same exact angular-charge leaf*. It explicitly rejects historic mixed-charge kick secants as a full Floquet spectrum.

## The question

On the **late prescribed negative-drive return branch**, is there a relative fixed point of the phase-stroboscopic map after one `-2*pi` external phase turn? If a root is genuinely found, is the fixed-J reduced return map locally attracting? In other words, can the observed near recurrence be promoted to a relative periodic orbit and a conditional fixed-J limit-cycle claim?

## Exact chart and constraint

`Y=(logr1,v1,theta1,omega1,logr2,v2,theta2,omega2,z,vz)`; quotient by counterrotation `(theta1,theta2)->(theta1+eps,theta2-eps)` and fix `J = exp(2*logr1)*omega1 - exp(2*logr2)*omega2`. The eight coordinates are `u=(logr1,v1,omega1,logr2,v2,z,vz,alpha)`, with `alpha=theta1+theta2 (mod 2*pi)`. Lifting always enforces

`omega2 = (exp(2*logr1)*omega1 - J_fixed)/exp(2*logr2)`.

The 10x8 analytic lift Jacobian is independently compared with centered finite differences, and its columns must satisfy `grad J · column = 0`. Counterrotation gauge covariance is directly falsified on a re-integrated orbit. Chart Euclidean norms and singular values are **descriptive**, not kinetic-metric invariants.

## Source and phase genealogy

- Native source `v18.52l/results/extended_horizon_10d_states.npz` (20×65,537 state arrays); **column 11 is material time tau, not phase**. Use native original and autonomous reference starting at `xi=12.8`, not the mixed-charge perturbations.
- SHA-256 authenticated `v18.51s` sendback containing the authentic stage loader, original physical RHS, original external theta clock, and nested historical ancestors.
- Continue the actual 10D RHS via `engine.physical.rhs_from_original([Y, theta(event+xi)], kappa_fixed, engine.params, engine.slope, cubic_clock_sign(event+xi), ...)`. On this branch sign is constant `-1` and the phase decreases by `2*pi` in `2*pi/engine.omega` of support time. **Never substitute the autonomous +x 11th clock or an invented phase.**
- The same `kappa_fixed=0.1666331385421536` is used for both families, not the historical seed-construction kappas.

## Falsification gates and outputs

1. Verify both SHA-256 input archives, authentic historical engine parameters, clock sign, native time range and first native reference state.
2. Build and verify the rank-eight fixed-charge tangent lift; test counterrotation by changing the gauge before propagation.
3. **Replay one authentic physical turn and compare with the independently archived native endpoint**; abort if incompatible.
4. Propagate one and two full phase turns from the same original native state, and check `P(P(u))` against the direct two-period integration. These are actual extensions beyond the historical one-step observation.
5. Compute each of the eight return-map derivatives via **two independent centered finite-difference scales**, with true same-J re-integrations. Record their convergence, singular values and directional columns. At a non-fixed point, Jacobian eigenvalues are descriptive one-step map diagnostics, **not Floquet multipliers**.
6. Try a limited, damped Newton search for `P(u)-u=0` on the same J leaf. Never force closure. Recheck the candidate's one-turn and two-turn closure independently. **Only** if both close, recompute its Jacobian and conditionally report fixed-charge Floquet multipliers if the two-step discretizations agree.
7. Emit compact scientific report, per-column table, per-seed charts, hashes, run failure diagnostics if necessary, and a full sendback archive containing actual scripts.

## Command on the server

```bash
cd /root/elab1/hysteresis
unzip q550_fixed_charge_tangent_basis_relative_closure_search_and_return_map_stability_v18_52y_bundle.zip
cd q550_fixed_charge_tangent_basis_relative_closure_search_and_return_map_stability_v18_52y_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
  bash run_on_server.sh --require-parent --require-native \
  2>&1 | tee ../v18_52y_server.log
```

Expected authentic inputs, both automatically discovered in `/root/elab1/hysteresis`:

- `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`
- `q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz`

If paths differ, use `--parent-v1851s /absolute/path/parent.zip --native-v1852l /absolute/path/states.npz`. Dependencies: Python 3.13, NumPy, SciPy; Matplotlib is optional for two plots.

**This audit now requires numerical 10D integrations.** They are intentionally deferred to the user's server. Default effort is two seeds, 32+ finite-difference integrations plus a bounded closure search on each; configure using `--seeds original`, `--skip-closure`, `--max-newton 1` for a smaller exploratory run, or tighten `--rtol`, `--atol`, `--maxstep` for independent solver convergence. Because matrix conditioning and dynamical stiffness can make this demanding, failures are retained and no automatic success is assumed.

## Return package

`q550_fixed_charge_tangent_basis_relative_closure_search_and_return_map_stability_v18_52y_bundle_server_sendback.zip` (next to the extracted root). Key files:

- `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`
- `results/eight_direction_tangent_jacobian.csv`
- `results/relative_return_summary.csv`
- `results/original_fixed_charge_return_map.png`, `results/autonomous_fixed_charge_return_map.png` (if plotting available)
- `OUTPUT_MANIFEST.json`; all primary/reference scripts, config, SHA lineage, failure trace when relevant

## What is not implied

A period of the external forcing is not proof of full physical periodicity. A single near return is not an attracting cycle. One eigenvalue calculation at a non-fixed point is not a Floquet certificate. Attraction, even if eventually established, is restricted to the specified **fixed-J, gauge-reduced, externally driven return model**; it does not prove total-energy conservation or a global conservative limit cycle, nor a complete forward–return cycle.
