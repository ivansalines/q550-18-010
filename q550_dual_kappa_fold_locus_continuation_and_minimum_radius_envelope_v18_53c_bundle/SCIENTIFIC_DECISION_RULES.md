# Scientific decision rules — Q550 v18.53c

`AUTHENTIC_DUAL_KAPPA_FOLD_LOCUS_AND_SAMPLED_CRITICAL_RADIUS_ENVELOPE__NOT_GLOBAL_MINIMUM_OR_INTERVAL_PROOF` requires all of:

- SHA-256 matching v18.53a authentic report and v18.53b authenticated report and refined fold, with v18.53b `parent_v1851s_verified`, `source_checked_refined_fold` and native fold match set.
- The original v18.51s generator archive is staged and its SHA-256 matches; original model parameters match the established Q550 configuration.
- Each of the high-resolution 102 narrow fold samples satisfies the augmented stationary/fold equations and the nondegeneracy and 8D spectral tests.
- Every accepted narrow and extended root, plus the two original/autonomous endpoint roots and same-drive control, agrees with the authentic original physical RHS and 8D Jacobian within configured thresholds.
- The two κ endpoints and intermediate original shared κ anchor are explicitly identified as **parameter prescriptions**, not independent physical seed branches.

`OFFLINE_FOLD_LOCUS_SCREENING__ORIGINAL_ENGINE_NOT_VERIFIED` is the only permitted successful local-only classification.
`AUTHENTIC_PARENT_PRESENT_BUT_FOLD_LOCUS_SOURCE_CHECK_INCOMPLETE` reports any missing or failed authentic checks and must fail strict mode.

For any corrector loss outside the narrow interval record `CONTINUATION_STOP_NOT_FOLD_DISAPPEARANCE`: do not infer a physical branch ending.

**The following claims are explicitly excluded even on full PASS:**

- `r2_fold` is a guaranteed global lower bound of any actual trajectory.
- The two critical fold radii occur at one and the same fixed signed-drive parameter.
- A sampled monotonic derivative proves a continuous global envelope between endpoints, or uniqueness of the entire fold locus.
- The v18.53b saddle-node on the constant-drive family occurred in the historical forward–return trajectory.
- The system has globally conserved energy or a globally conservative Hamiltonian formulation.
