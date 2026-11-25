# Local preflight — Q550 v18.53v

**Status:** PASS for source authentication and synthetic/postprocessing tests. The complete 204-node scientific calculation has deliberately **not** been executed here; it belongs on the user's Q550 server.

- Historical authentic v18.53u server sendback: ZIP CRC PASS; source SHA manifests 14/14; output SHA manifests 8/8. Original parent bytes retained as an immutable nested ZIP.
- Historical v18.53t embedded server sendback: SHA-256 pinned in `inputs_config.json`; source manifest 13/13; output manifest 7/7.
- Historical v18.53u metric helper: byte-for-byte authenticated SHA-256 matched.
- Mathematical synthetic tests: induced metric differentiation; determinant normalization; analytic 2-parameter gradient/Hessian; mixed parameter derivative; implicit response of a separate synthetic quadratic least-squares model. PASS.
- Authentic physical-data smoke replay: original/shared seed-rule, preformal metric, λ = 0, 0.5, 1; minimum continuation, 2x2 Hessian, independent competitor, oriented curve. PASS (3 nodes). No new physical IVPs.
- Synthetic packaging regression: ZIP entries built from `Path` objects, CRC and SHA manifests reverified. PASS.
- Historic missing `results/` directory regression avoided by explicit `OUT.mkdir(parents=True, exist_ok=True)`.
- Preflight cannot certify the actual branch-persistence outcome or a material closure. The genuine 204-point server run determines the output classification.
