# Q550 v18.53a — Phase-Locked Relative Equilibrium, Cross-Seed Equivalence & Instantaneous RHS Audit

## Scientific question

The source-authenticated v18.52z three-turn return map has a near-machine-precision fixed point modulo counter-rotation, with a numerical 8D quotient spectrum inside the unit disk. Is that fixed point an **instantaneously stationary relative equilibrium** of the Q550 equations on the constant-negative external return drive, rather than a nontrivial periodic orbit? Do the original/autonomous seeds represent the same locked state after correcting their external phase origins? Does the **instantaneous reduced generator** reproduce the measured conditional Floquet spectrum?

## Independent tests

1. SHA-256 verify original v18.52z server report and authenticate native v18.52l states and v18.51s executable ancestry. No guessed historical parameters and no substitution of the 11th native `.npz` column (it is material time tau, not Theta).
2. From the v18.52z candidate and conserved J, solve three independent stationary equations at `Delta=pi, x=1/2`: two radial force balances plus `omega1+omega2=-x`; omega2 is eliminated exactly via J. The local algebraic root does not integrate a period.
3. Reconstruct authentic external phase from `ctx['theta']` at the native xi=12.8 anchor; evaluate 11D `engine.physical.rhs_from_original` with the negative late-return clock. Verify all eight instantaneous reduced components, including the clock-compensated angle sum and auxiliary x dynamics. Check counter-rotation gauge equivariance and J.
4. Compare candidate and rooted instantaneous RHS independently. Check alternative `Delta=0` phase branch is not a stationary radial solution (synthetic falsification).
5. Compare original and autonomous **after removing their separate phase origins** and the counter-rotation gauge. Raw absolute theta sums may differ; this is not evidence of two physical solutions.
6. Compute 8x8 reduced Jacobian by three finite-difference scales and compare `exp(lambda_i*T_xi)` against the archived v18.52z one-turn Floquet eigenvalues using optimal complex matching. Do NOT interpret this as the complete 11D spectrum or a full historical forward-return cycle.

## Server execution

```bash
cd /root/elab1/hysteresis
unzip q550_phase_locked_relative_equilibrium_cross_seed_equivalence_instantaneous_rhs_v18_53a_bundle.zip
cd q550_phase_locked_relative_equilibrium_cross_seed_equivalence_instantaneous_rhs_v18_53a_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh --require-parent --require-native 2>&1 | tee ../v18_53a_server.log
```

No new long-horizon 10D IVPs; authenticated stage may execute a brief historical make_context IVP. This bundle deliberately remains small and does not embed the 97MB native data or the multi-ancestor parent. Both are required on the server, with their original SHA-256 fingerprints.

Expected external inputs:
- `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` — SHA-256 `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`.
- `q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz` — SHA-256 `92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12`.

The sendback ZIP lands beside the extracted bundle and contains executable scripts, input provenance, JSON and scientific interpretation. If a required source is missing, the strict run fails explicitly and still packages a partial report. Run `PYTHON_BIN=python3 bash run_on_server.sh` for a nonstrict offline algebraic smoke check (not a Q550 server certificate).
