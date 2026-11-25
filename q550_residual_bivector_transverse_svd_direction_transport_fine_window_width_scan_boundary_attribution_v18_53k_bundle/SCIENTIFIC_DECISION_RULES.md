# Scientific decisions — v18.53k

- `PREFLIGHT_PASS__NO_NEW_IVP_NO_NEW_WIDTH_EVIDENCE`: source/manifest checks and synthetic tests only. Do not interpret as physical scan.
- `PARTIAL_AUTHENTIC_J_REFERENCE__NEW_PHYSICAL_IVPS_NOT_EXECUTED`: v18.53j ancestor passes SHA/authentication, but required historical parent or native NPZ is missing. No new physical evidence.
- `AUDIT_FAIL_NOT_CERTIFIED`: parent mismatch, failed source/replay, failed numerical or integrator check. Inspect report/error, do not promote results.
- `AUTHENTIC_FINE_WIDTH_SCAN__NO_RESOLVED_SIGN_CROSSING`: all ten required IVPs complete and 756 triplets computed, but no signed third-direction coefficient changes sign under the continuity/conditioning gates.
- `AUTHENTIC_FINE_WIDTH_SCAN__SIGNED_TRANSVERSE_CROSSING_CANDIDATES_ONLY`: all integrations and reconcilations complete; at least one **sampled sign change** survives orientation continuity gates. It remains a conditional numerical candidate, NOT a certified root or dynamical flip.

An apparent minimum located exactly at a scanned grid value is only a *sampled* minimum. No continuous extremum, recurrence period, or physical holonomy follows. Report rank gaps, quadrature sensitivity and interval endpoint derivative checks; do not treat exact interval-additivity as evidence of robustness.
