# Q550 v18.53n-fix1 — Minimal runtime repair

## Authentic failure
- Original server sendback: `provenance/q550_full_10d_residual_two_parameter_near_return_minimum_boundary_source_robustness_v18_53n_bundle_server_sendback.zip`.
- Server report: `FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE`.
- Physical IVPs completed: 1 of 10.
- Python traceback at `scripts/full_net_geometry.py`, inside `endpoint_source_diff`: `net_j` was undefined; available function: `net_jac`.

## Changes
1. Replace exactly `_,J=net_j(p,at,slope)` with `_,J=net_jac(p,at,slope)`.
2. Add a synthetic test calling `endpoint_source_diff` with 10D primitive and endpoint source; verify returned Jacobian agreement.
3. Mark new run version v18.53n-fix1 and update README.
4. Include original failing ZIP byte-for-byte, parent m-fix1 ZIP byte-for-byte, and all scripts.

No changes to physical equations, κ0, domain box, objective, minimization tolerances, ancestor inputs, numerical integrations, or scientific decision rules.

## Interpretation
Original v18.53n had **no scientific certificate**. Fix1 must run the full Q550 experiment on the user's server; synthetic/preflight tests do not replace physical results.
