# Q550 — v18.53b
## Relative-Equilibrium Branch Continuation, Fixed-Charge Stability Boundary & Bifurcation Audit

This bundle continues the **single authenticated phase-locked Δ=π relative-equilibrium branch** from v18.53a, not the unclosed historical trajectory. The original and autonomous seeds are phase/gauge equivalent at the common anchor; they are not counted as separate branches.

The original 11D physical RHS is staged from the historical v18.51s SHA-pinned sendback on the server. Its instantaneous field and 8D Jacobian are compared at **every accepted branch root** against the independently reconstructed equations. No extended orbital integration is needed.

### Server run

From `/root/elab1/hysteresis`:

```bash
unzip q550_relative_equilibrium_branch_continuation_fixed_charge_stability_boundary_and_bifurcation_v18_53b_bundle.zip
cd q550_relative_equilibrium_branch_continuation_fixed_charge_stability_boundary_and_bifurcation_v18_53b_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh --require-parent 2>&1 | tee ../v18_53b_server.log
```

You may pass `--parent-v1851s /absolute/path/to/q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip` if the ancestor is elsewhere. A local-only smoke test may use `bash run_on_server.sh --offline-only`; its outputs are **analytic screening only**, not a source-authenticated Q550 result.

### Independent scans

- `charge`: J from −2.2 to −0.7, holding kappa and constant signed drive s fixed. **Each point belongs to a different J leaf; the 8D stability is computed tangent to that point's leaf.**
- `kappa`: kappa from 0.085 to 0.25 at fixed J and s.
- `signed_drive`: constant s from −1 to +1 at fixed J and kappa. This is a **counterfactual constant external-control family**, not a replay of the historical cubic forward/return ramp. At s=0 phase does not turn: the phase-turn multiplier is deliberately undefined.

The ordinary scans follow the Δ=π sheet by predicted Newton correctors. **A separate pseudo-arclength continuation** then resumes from the ordinary signed-drive branch tip, solving a four-dimensional augmented root problem and following the same sheet through a possible parameter turn. The turning point and a real eigenvalue crossing are candidate saddle-node evidence, not a nondegenerate bifurcation theorem. The script does not claim global branch exhaustiveness.

### Outputs

`results/report.json`, `results/branch_continuation.csv`, `results/pseudo_arclength_drive.csv`, `results/pseudo_arclength_fold_candidates.json`, `results/refined_fold_nondegeneracy.json`, `results/branch_detailed_spectrum.json`, `results/stability_boundary_candidates.json`, `results/corrector_stops_not_proven_folds.json`, `results/branch_radius_and_spectral_abscissa.png`, `results/SCIENTIFIC_REPORT_SERVER.md`, plus the full source bundle, source/output manifests and upstream provenance. The wrapper packages a complete sendback even on a strict failure.

### Scientific exclusions

No global conservative/Hamiltonian claim; no global radial noncollapse theorem; no physical full forward–return cycle; no quantized J; no bifurcation certificate from a missing root or sign change. The genuine v18.53a anchor is pinned by SHA-256.
