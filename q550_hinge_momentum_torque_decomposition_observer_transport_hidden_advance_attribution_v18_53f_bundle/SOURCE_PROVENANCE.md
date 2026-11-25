# Provenance and executable genealogy

- Native original 10D v18.52l: `extended_horizon_10d_states.npz`; SHA-256 pinned in `inputs_config.json` (never embedded in ZIP).
- Relative physical 5D v18.52l: `extended_horizon_5d_vectors_and_traces.npz`; SHA-256 pinned in config (never embedded).
- Original source v18.51s nested archive: `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`; SHA-256 pinned in config; used to rehydrate true 11D physical RHS and historical `ctx['theta']`.
- `reference/parent_v18_51t_run_audit.py`: historical engine loader for nested parent.
- `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py`: historical independent RHS source.
- `reference/v18_52v_r1_run_audit.py`: original signed clock and shared-kappa reconciliation implementation.
- `scripts/energy_math.py`, `scripts/phase_energy.py`: original model-energy and clock work math references, included unchanged.
- `reference/kinetic_metric_v18_49h.py`: SHA-pinned exact 10D->5D kinetic embedding.
- `reference/v18_52n_authenticated_report.json`: SHA-pinned observer frame and hidden PCA direction.
- `reference/v18_53e_authenticated_report.json` and `reference/v18_53e_fixed_endpoint_pair_summary.json`: SHA-pinned original v18.53e output ancestry. They are also verified against the carried `reference/v18_53e_OUTPUT_MANIFEST.json`.
- `reference/v18_53e_run_audit.py`: actual source that created endpoint selections and absolute hinge geometry.

The v18.53f Python code is the complete **new** execution source; historical scripts are retained to audit functional genealogy. The runtime preserves both reported evidence and input SHA-256, instead of relying on historical filenames alone.
