# Q550 v18.53m-fix1 — Numerical trace-domain repair

## Observed physical sendback failure (2026-11-09)

- Original authentic Q550 v18.53m server sendback: SHA-256 `a9348cf4eb3b9de75d75396dd3635e5ad68e2094272a5ad8b0e2e9c3b7b7c5dd`.
- ZIP CRC PASS; original SOURCE_MANIFEST.json 16/16 and OUTPUT_MANIFEST.json 19/19 PASS.
- Parent v18.53l authentic source and outputs verified; original external physical inputs authenticated.
- Actual Q550 IVPs=1/10; error `ENDPOINT_OUTSIDE_TRACE` inside scipy.optimize.least_squares.
- Source: the data start at xi=12.8; optimizer was permitted to attempt origins down to 12.6. Its trial called primitive_at below xi[0] and correctly refused extrapolation.
- The failed sendback cannot establish the existence or nonexistence of any projected root.

## Changes in fix1

1. Correct rectangular optimizer origin lower bound from 12.6 to 12.8. The 25x29 grid, 3 anchor planes, physical equations, root tolerances and all pinned ancestors are unchanged.
2. `validate_trace_coverage` now rejects any optimizer bounds whose start or final endpoint fall outside the actual xi trace. It runs **before** all Q550 IVPs and at every root optimization; no clipping of the physical integral occurs.
3. The optional finite-difference Jacobian *diagnostic* uses one-sided differences if a centered step would leave the trace; it does not replace the exact piecewise-linear endpoint Jacobian used by the optimizer.
4. New `tests/test_trace_domain_regression.py` reproduces the old invalid bounds on a 12.8–56.0 synthetic trace, asserts rejection, then confirms safe optimization with the new bounds.
5. The original failed server ZIP is kept unmodified in `provenance/v18_53m_failed_server_sendback.zip` for audit traceability. All historical v18.53l sources and inputs remain byte-identical.

## Scientific boundaries

- The fix has **not** generated new physical Q550 integrations locally; all numerical scientific results remain pending server execution.
- The new search box excludes 12.6<=o<12.8, a region outside the actually available trace; it does not exclude previously sampled grid points.
- A 2D projected zero remains observer-dependent and is not evidence of full 10D closure, material memory, or temporal holonomy.
- All new comments, scientific docs and test output are in English.
