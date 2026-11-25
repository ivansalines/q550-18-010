# Local preflight validation — Q550 v18.53y

Executed with Python 3.13 and single-threaded BLAS:

- Source manifest: PASS, 15 authenticated files in the tested build before distribution rehash.
- Minimal-turn rotation, exterior-square generator and analytic SO(5) gauge tests: PASS.
- Path-safe package/sendback mock: PASS with CRC + manifest checks.
- Authentic x source/output manifests: 14/14 + 9/9, plus authenticated w (13+7), v (14+10), u (14+8), t (13+7): PASS.
- Authentic Q550 smoke (original seed, shared kappa, preformal calibration, lambda in {0,1,4}): 3 parent paths, 9 start/midpoint/end lifts, PASS.
- Maximum parent x fixed physical-frame angle reconstruction in smoke: 2.159e-12 degrees.
- Maximum time-dependent-gauge pullback error in smoke: 5.857e-17.

No new physical 11D/21D integration, and no full 204-path audit performed locally. Full scientific output must come from the server.
