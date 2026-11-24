# Executable provenance

This audit is based on byte-preserved archives supplied by the Q550 research workflow:

- v18.51r nonlinear state tensors, recovered from the exact v18.51r server sendback and cross-verified against the same nested bytes in the authenticated v18.51s parent.
- v18.52c signed-margin physical midpoint sendback (24 completed genuine midpoint probes), retained as an immutable ZIP.
- `parent_v18_51t_run_audit.py` from v18.52a, providing the original nested-model staging function without reimplementing a historical solver.
- `kinetic_metric_v18_49h.py` from the authenticated v18.51p lineage, defining the ten-dimensional Cartesian kinetic embedding and exact Jacobian.

For the full stage, the user provides the original authenticated v18.51s server sendback. `source_manifest.py` checks all included source files, and the server script checks the exact parent ZIP SHA256 and compares the embedded v18.51r NPZ bytes with the separately supplied cached copy. All output files receive their own SHA256 manifest after the run.

An identical computational purpose and mathematical map, not an ungrounded reconstruction from plot pixels, define the executable genealogy.
