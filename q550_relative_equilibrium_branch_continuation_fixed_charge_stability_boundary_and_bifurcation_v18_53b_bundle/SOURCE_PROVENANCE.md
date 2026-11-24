# Source provenance and executable genealogy

- `reference/v18_53a_authenticated_report.json`: actual **server sendback** scientific report from `q550_phase_locked_relative_equilibrium_cross_seed_equivalence_instantaneous_rhs_v18_53a_bundle_server_sendback.zip`. It contains the generator-certified equilibrium and spectral comparison. SHA fixed by `inputs_config.json`.
- `reference/v18_53a_run_audit.py`, `reference/v18_53a_inputs_config.json`: preceding executable audit and its model config.
- `reference/v18_52v_r1_run_audit.py`: historically validated staging interface `find_parent` / `stage_parent` for the authenticated v18.51s archive.
- `reference/parent_v18_51t_run_audit.py`: original engine stage loader preserved as historically used.
- `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py`: independent archived RHS used by source-crosscheck tests.
- `scripts/energy_math.py`, `scripts/return_map.py`, `scripts/phase_energy.py`: inherited mathematical support and fixed-charge quotient implementations.
- `scripts/run_audit.py`: actual new v18.53b continuation/generator comparison.
- `tests/*`: independent archive RHS comparison and intentionally failing false boundary/false drift controls.

The large historical v18.51s source ZIP is **not embedded**. The server must provide it at `/root/elab1/hysteresis` or via `--parent-v1851s`. The exact file SHA and v18.53a report hash are machine enforced. Historical native 10D arrays are not required for this stationary branch audit, and no invented orbit is used.
