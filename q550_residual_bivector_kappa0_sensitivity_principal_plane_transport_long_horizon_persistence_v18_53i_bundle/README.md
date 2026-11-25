# Q550 — v18.53i: Residual-Bivector κ₀-Sensitivity, Principal-Plane Transport & Long-Horizon Persistence Audit

Reproducible, source-authenticated late-return Q550 numerical audit. Designed to run genuine parameter-distinct ODE integrations **on the user's server**; locally tested without expensive production Q550 integrations.

## Scientific starting point

The authenticated v18.53h fix1 sendback certified 10/10 real IVPs on ξ∈[12.8,25.6], exact generator-source comparison, shared-κ historical replay, and local agreement of principal residual-torque planes between **original κ₀=1/6** and **autonomous κ₀≈0.16659961041764054**. The 5D residual is the material-time integral of the radial-potential plus phase-coupling torque bivectors. This audit asks whether the alignment persists over a wider κ range and three consecutive late-time windows. It does NOT presume it persists.

## Parameter sweep and integration design

- Two authenticated seed families (`original`, `autonomous`), both initialized from original native v18.52l at **ξ=12.8**; within each family all κ variants have identical initial states and the same source historical external phase prescription.
- 11 genuinely distinct κ₀ choices: original, autonomous, shared κ₀, and shared ±{0.25, 1, 3, 10}×Δκ₀, where Δκ₀=κ₀,original−κ₀,autonomous. Thus 22 primary IVPs.
- 4 additional higher-accuracy IVPs, refining both historical κ values for both seed families: **26 IVPs** total.
- Default endpoint ξ=51.2, with stages [12.8,25.6], [25.6,38.4], [38.4,51.2]. Also calculate cumulative intervals from 12.8 to 25.6, 38.4, and 51.2.
- The native historical archive ends at **ξ=25.6**. Subsequent forcing is a *formal continuation of the authenticated original-source clock function*, **not** additional authenticated historical observations nor the historical pre-12.8 ramp. No full-cycle stability or asymptotic-in-time theorem is inferred.
- Historical ξ-grid and material-clock convention preserved on [12.8,25.6]; extended samples at step 0.025.
- Source agreement, conserved counterrotation charge, two-level integrator refinement, historical shared-κ replay, reference v18.53h residual reproduction, and bivector quadrature controls are acceptance checks.

## Observables

For each accepted run and temporal interval:

`R = ∫ [Q∧F_radial + Q∧F_phase] dτ`, `R` a 5×5 skew tensor. A real Schur decomposition produces singular strengths σ₁,σ₂ and oriented rank-two projectors P₁,P₂. The algorithm logs **singular-gap** warnings; individual principal planes are not stable observables at degeneracy. Gauge-invariant distances, projectors and two-dimensional principal angles quantify geometry.

The parameter response is measured by the complete-bivector finite difference `dR/dκ`, singular-value derivatives, and `dP_i/dκ` in Frobenius norm for four predetermined κ half-widths. A true perturbative derivative requires convergence of the shrinking-step estimates, not merely small finite differences. Consecutive interval projectors measure **temporal transport**, kept distinct from the cumulative-window results.

**Cancellation scores are integrated signed bivector diagnostics, not pointwise torque-balance identities.** The numerical κ range (±10×Δκ about shared) is an exploratory finite range, not a proof of a parameter-invariant law.

## Server execution

```bash
cd /root/elab1/hysteresis
unzip q550_residual_bivector_kappa0_sensitivity_principal_plane_transport_long_horizon_persistence_v18_53i_bundle.zip
cd q550_residual_bivector_kappa0_sensitivity_principal_plane_transport_long_horizon_persistence_v18_53i_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 \
 bash run_on_server.sh --require-parent --require-native --require-dynamic \
 2>&1 | tee ../v18_53i_server.log
```

For explicit source locations:

```bash
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
 --parent-v1851s /root/elab1/hysteresis/q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip \
 --native-npz /root/elab1/hysteresis/q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz \
 --require-parent --require-native --require-dynamic
```

**Do not overwrite v18.53h fix1.** ZIP root, output folder, README, source manifest and all code use this complete unique v18.53i name. The large SHA-pinned ancestors are external and are automatically discovered on the server using the authenticated v18.53h loader. They are intentionally not duplicated in the portable ZIP.

To change the continuation horizon, pass `--xi-max 38.4` or `--xi-max 76.8` (the default is 51.2). Larger horizons increase computational load and extrapolate the forcing farther beyond the archived window. For the standard 51.2 run, expect exactly 26 successful IVPs. The sendback records the **actual count**, even on a later failure.

## Output files

- `results/report.json`: classification, SHA/source statuses, true IVP count, material and numerical checks; `SCIENTIFIC_REPORT_SERVER.md` is readable interpretation and scope.
- `integrated_kappa_runs.csv`: all IVPs, conservation and RHS checks, historical replay diagnostics.
- `residual_bivector_windows.csv`: full 5D skew residual, principal singular values, cancellation, conditioning per interval.
- `kappa_principal_plane_separation.csv`: principal angles, cosine, relative residual change versus shared κ₀ at *matching endpoints*.
- `kappa_local_sensitivity.csv`: κ₀ finite differences of bivector, principal projectors, strengths, cancellation.
- `sensitivity_resolution_comparison.csv`: finite-step stability diagnostic across ±0.25Δκ and ±Δκ.
- `principal_plane_time_transport.csv`: consecutive-segment plane motion for original/autonomous/shared κ₀, calculated on each self-consistent trajectory.
- `baseline_53h_replay.json`: independent reproduction of authenticated v18.53h on historical interval and refined orbit checks.
- `sampled_native_radius_charge.csv`: small sample of radius and charge for plotting/review; not full 10D saved orbit.
- `principal_plane_sensitivity_original.png`, `principal_plane_sensitivity_autonomous.png`, `principal_plane_transport.png`: figures if matplotlib is available.
- `OUTPUT_MANIFEST.json` plus parent-level `_server_sendback.zip`: source + output hashes.

In partial mode with no large ancestors, the script emits an explicitly incomplete report and returns nonzero under `--require-*` flags. **Partial sendback is not a scientific certificate**. Validate version and actual IVP count before interpreting results.

## Dependencies

Python 3.13 recommended, NumPy, SciPy, Matplotlib. Core script and reference scripts are bundled. Read `SOURCE_PROVENANCE.md`, `MATHEMATICAL_PROTOCOL.md`, and `SCIENTIFIC_DECISION_RULES.md` for audit boundaries.
