# Q550 v18.51u — restart-axis repair patch

Issue: `ValueError: operands could not be broadcast together with shapes (2,3) (3,2)` in `restart_impulse_checks`.

Root cause: `actual = L[-1]@z` has (projection, source-channel) orientation `(2,3)`, as does `K[-1,j,:,:3]`. The original code incorrectly transposed the latter.

Repair: remove that transpose, assert shape equality, and add a regression to the independent synthetic impulse-restart benchmark. Add the optional `restart_indices` argument so the test can use its short grid. Default production restart indices `(0,128,256,384)` remain unchanged. No model parameters, sources, input data or solver tolerances are altered.

The enclosed `SOURCE_MANIFEST.json` includes the recomputed hashes of modified files and also authenticates this note. All other hashes match the original bundle. This is a source-only overlay applied to the **existing extracted full v18.51u root directory**, preserving the 86 MB ancestor input.

From `/root/elab1/hysteresis`:

```bash
unzip -o /root/Downloads/q550_time_resolved_green_kernel_compensation_birth_v18_51u_restart_axes_patch.zip
cd q550_time_resolved_green_kernel_compensation_birth_v18_51u_bundle
OPENBLAS_NUM_THREADS=1 python3.13 scripts/verify_source_manifest.py
OPENBLAS_NUM_THREADS=1 python3.13 tests/test_green_kernel_synthetic.py
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh 2>&1 | tee ../v18_51u_server_repaired.log
```

Do not run the patched `run_audit.py` against the original unpatched `SOURCE_MANIFEST.json`, or the source verification will correctly refuse it.
