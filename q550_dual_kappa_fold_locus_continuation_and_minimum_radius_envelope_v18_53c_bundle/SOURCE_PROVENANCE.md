# Source provenance and executable genealogy — v18.53c

- **`reference/v18_53b_authenticated_report.json` and `reference/v18_53b_authenticated_fold.json`** are extracted byte-for-byte from the user-supplied v18.53b server sendback. Their exact SHA-256 digests are checked against `inputs_config.json`.
- **`reference/v18_53a_authenticated_report.json`** is the SHA-pinned v18.53a original-equilibrium certificate inherited from v18.53b, including fixed charge and underlying equilibrium state.
- **`reference/v18_53b_run_audit.py`** preserves the executable stationary/reduced field and v18.53b original-engine verification used in this step. The v18.53c script uses these functions by importing that exact file, rather than reimplementing the RHS.
- **`reference/v18_52v_r1_run_audit.py`** supplies `find_parent` and `stage_parent` to stage the *genuine* v18.51s source on the user's server.
- **`reference/parent_v18_51t_run_audit.py`**, **`reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py`**, **`reference/v18_53a_run_audit.py`** preserve older executable lineage.
- **`scripts/energy_math.py`, `scripts/return_map.py`, `scripts/phase_energy.py`** preserve mathematical and quotient structure.
- **`scripts/run_audit.py`** implements a new analytical fold Jacobian, four-variable fold correction, implicit κ derivatives, dual/extended continuation, sampled radius envelope, original-engine crosschecks, and a distinct same-drive comparison.
- **`tests/test_synthetic.py`** tests false sign changes, derivative agreement, fixed-J fold calculation, and the analytic 3×3 Jacobian against independent finite differences. **`tests/test_archived_rhs_crosscheck.py`** probes the separate historical support RHS at 64 arbitrary states.

The large v18.51s source ZIP is NOT embedded. Its expected SHA-256 is pinned to `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`, and the server MUST supply it to receive a genuine `AUTHENTIC...` conclusion.

The generated offline screening preview (if present) is a model-derivation diagnostic only, NOT a substitute for native original-engine confirmation.
