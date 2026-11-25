# Q550 v18.53n-fix1 — Local validation (no new Q550 physical integrations)

## Observed authentic failure (original sendback)

- ZIP CRC: PASS (19 members).
- Source manifest: 14 source files authenticated by SHA-256.
- Output manifest: 3 result files authenticated by SHA-256.
- Server decision: `FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE`.
- Physical integrations: 1/10; refined 0/4.
- Traceback cause: `endpoint_source_diff` called undefined `net_j`; the implemented derivative routine is `net_jac`.

## Fix scope

Only the function reference and direct regression coverage were changed. No changes to physical Q550 equations, kappa rules, time domain, optimizer parameters, dynamic sources, or scientific acceptance criteria. The failing original archive is preserved byte-identically in `provenance/`.

## Verification performed locally

1. Full parent m-fix1 archive CRC and manifest authentication: PASS (20 source + 34 output files).
2. Source manifest validation for this fix1 distribution: PASS (see generated manifest for count).
3. `python3 tests/test_synthetic.py`: PASS, including direct call to `endpoint_source_diff` and boundary-source agreement to <1e-10 for a consistent synthetic source and primitive.
4. `python3 tests/test_analyze_case_smoke.py`: PASS: entire `analyze_case` path on a declared synthetic 10D history, 2279 grid points, 22 multistart minimizations, parent-baseline comparison, PCHIP sensitivity and endpoint attribution.
5. `python3 scripts/run_audit.py --preflight`: PASS, no Q550 integration performed.
6. Python compilation, ZIP CRC, source/output manifest and extraction-root check: PASS in final packaging validation.

## Scientific status

**No v18.53n-fix1 physical conclusion yet.** The original sendback cannot certify a full-10D minimizer or material return, because it aborted after one physical integration. Running fix1 on the user's Q550 server remains necessary.
