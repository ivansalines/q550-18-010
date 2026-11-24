# Q550 v18.53b — Relative-Equilibrium Branch Continuation, Fixed-Charge Stability Boundary & Bifurcation Audit

**AUTHENTIC_FIXED_CHARGE_RELATIVE_EQUILIBRIUM_BRANCHES_WITH_LOCAL_SPECTRAL_BOUNDARY_DIAGNOSTICS__NO_GLOBAL_BIFURCATION_THEOREM**

The physical source has been checked at each root only if the authenticated v18.51s parent was available.

## Branch scans

### charge
- Continuation points: 62
- Stable reduced 8D: 62; unstable: 0
- Candidate stability-boundary intervals: 0
- Corrector stops: 0 (not proved folds)

### kappa
- Continuation points: 62
- Stable reduced 8D: 62; unstable: 0
- Candidate stability-boundary intervals: 0
- Corrector stops: 0 (not proved folds)

### signed_drive
- Continuation points: 78
- Stable reduced 8D: 78; unstable: 0
- Candidate stability-boundary intervals: 0
- Corrector stops: 1 (not proved folds)

## Interpretation

Continuation in charge compares different invariant charge leaves, with stability always computed tangent to each leaf.
Continuation in kappa or constant signed-drive holds charge fixed. A signed drive value is a counterfactual constant protocol parameter, not an actual reconstructed ramp trajectory.
A local spectral sign crossing does not alone prove a Hopf, pitchfork, saddle-node or global bifurcation. A failed corrector is not evidence of a fold.
No Poincare orbit was reintegrated here. The inherited v18.53a one-turn Floquet match applies only to its authenticated anchor.
At s=0 the external phase does not turn: no one-phase-turn Floquet multiplier is defined.
The two seed offsets are gauge/phase equivalent at the anchor and not separately continued as different physical branches.
The signed-drive pseudo-arclength branch is followed through the apparent turn; its numerical parameter reversal and real eigenvalue crossing are candidates for a saddle-node, NOT a nondegenerate bifurcation certificate.
Pseudo-arclength candidate: {"status": "NUMERICAL_PSEUDOARCLENGTH_FOLD_AND_REAL_STABILITY_CROSSING_CANDIDATE", "n_steps": 60, "n_failed_steps": 0, "max_signed_drive_including_tip": 0.9268409019848072, "max_signed_drive_position_in_sequence": 3, "parameter_turn_observed": true, "min_equilibrium_3x3_singular": 0.0008212623455905393, "min_spectral_abscissa": -0.020869318826671472, "max_spectral_abscissa": 0.4428854274006488, "native_all_checked": true}
Refined local fold diagnostic (not interval certificate): {"status": "NUMERICALLY_NONDEGENERATE_LOCAL_SADDLE_NODE_CANDIDATE__NOT_INTERVAL_CERTIFICATE", "J": -1.441016445824766, "kappa": 0.1666331385421536, "s_fold": 0.9268428336365817, "core_a": [0.5323827558021118, -1.1851911989057509, -0.46689946765986373], "r1": 1.702985276132668, "r2": 0.3056877267519828, "root_augmented_residual_max": 7.448163635837926e-11, "jacobian_equilibrium_singular_values": [89.8538967252151, 1.1369690330736044, 2.4274967852109883e-09], "parametric_transversality_left_null_dot_Fs": -0.4052573899870457, "quadratic_fold_curvature_left_null_dot_Faa_vv_by_stencil": [-0.43152706436967314, -0.4315243370232702, -0.43152410195735946], "smallest_8d_eigenvalue": {"real": 3.132888730289945e-08, "imag": 0.0}, "other_7_max_real": -0.04036235847711806, "conditions": {"root_augmented_max_abs_lt_1e_7": true, "3x3_null_singular_value_lt_1e_6": true, "parametric_transversality_nonzero": true, "quadratic_fold_curvature_nonzero": true, "quadratic_stencil_consistent": true, "8d_one_near_zero_real_mode": true, "other_7_modes_left_half_plane": true, "authentic_source_match": true}, "authentic_original_rhs_at_fold": {"available": true, "native_equilibrium_rhs_max": 1.9932976885231687e-15, "source_reconstructed_field_relative": 7.717299992935316e-15, "source_reconstructed_jacobian_relative": 4.54770043378145e-13, "native_jacobian_step_relative": 6.650331550023252e-09, "native_local_eigen_max_abs_difference": 1.364946419779578e-10, "source_match": true, "native_spectral_abscissa": 3.132888730289945e-08}, "classification": "NUMERICAL_LOCAL_SADDLE_NODE_NONDGENERACY_DIAGNOSTIC_NOT_GLOBAL_PROOF"}
