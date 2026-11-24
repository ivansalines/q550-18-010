# Decision rules — Q550 v18.52v-r1

- `AUTHENTICATED_CLOCK_SIGN_AND_FIXED_KAPPA_SIGNED_ENERGY_FULL_RESOLUTION`: v18.51s source SHA and v18.52l state SHA verified; all 20 10D series included; clock source differential sign certified (`theta_prime-omega*sign` max <1e-7); correct common fixed-kappa energy closure max <1e-8 at stride 1. The status means *numerical energy ledger and source genealogy*, not physical forward-return path equivalence or all-time theorem.
- `PARTIAL_AUTHENTICATION_OR_NUMERICAL_AUDIT_FAILED`: at least one numeric criterion not met, even if some source identities hold. Report full diagnostics.
- `NATIVE_AUTHENTICATED_PARENT_MISSING_NO_FORCED_PHASE_IMPUTATION`: no ancestor; neither external phase nor energy certificate invented. Strict mode fails.
- `INCOMPLETE_EXCEPTION`: runtime failure; preserve partial report and package sendback.

No forward physical trajectory is generated; frozen-state sign inversion never counts as a measured forward work integral. A true global noncollapse claim requires separate global-in-time estimates.
