# Q550 — v18.51s

**Parameter-Sensitivity Transport & Transverse-Posture Attribution Audit**

This is a two-level executable and genealogically authenticated study of the true Q550 nonlinear model. It follows the v18.51r server sendback, with its full nested q→p→... historical ancestry and the original seed/parameter identities.

**Level A (completed before packaging).** Authenticate the 87 MB v18.51r parent and reproduce the 513-section-centered secant parameter sensitivity at fixed original/autonomous seed, with independent physical Cartesian immersion differentiation and decomposition into radial/azimuthal/posture components. Read `results/LOCAL_CACHED_BASELINE.json` for the verified local output. The cache-only Level A passes 10/10 tests; it is NOT an independent parameter tangent ODE.

**Level B (server task).** Reactivate the authentic nonlinear RHS/analytic jet, compute `∂κf` on fixed midpoint states, and integrate `S'=J S+∂κf`, `S(0)=0`; independently check its state jet, parameter stencil, grid refinement and its agreement with the two already-integrated κ branches. The Level B server execution has **not** been completed in the packaging environment. It must NOT be reported as PASS until run and checked.

Run on the user's Python 3.13 server from the unpacked root folder:

```bash
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_51s_server.log
```

A completed run generates `../q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`, containing the manifest-authenticated parent, the exact source code, new results, logs produced externally when redirected, tests and reports.

For a quick independent baseline-only check, without staging the legacy engine:

```bash
OPENBLAS_NUM_THREADS=1 python3.13 scripts/run_audit.py --mode cached
python3.13 tests/test_results.py
```

This intentionally does NOT certify the proposed full parameter tangent test. Re-run `run_on_server.sh` to restore and complete the full scientific report. The entire source history is bundled in `inputs/...v18_51r_bundle_server_sendback.zip`, and all source comments/reports are in English. Dependencies: numpy, scipy, matplotlib; pandas may be required by inherited historical staging modules. No arbitrary replacement potential, viscosity, metric, observer gauge or external κ threshold has been added.
