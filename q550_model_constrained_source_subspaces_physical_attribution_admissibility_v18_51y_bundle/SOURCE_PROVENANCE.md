# Q550 v18.51y — Genealogy

- Immediate ancestor: actual Q550 v18.51x server sendback (`q550_physical_source_channel_recombination_geometric_attribution_robustness_v18_51x_bundle_server_sendback.zip`). The full v18.51y archive includes this exact parent, preserving all nested ancestral ZIPs. It is optional in the lightweight bundle only because T/W/X snapshots are SHA-256 locked.
- `inputs/t_green_channel_transport_states_v18_51t.npz`: original model-derived 10D Jacobian and parameter source, forced and Green responses.
- `inputs/w_channel_vector_geometric_invariants_v18_51w.npz`: physical second-plane source-vector geometry.
- `inputs/x_channel_recombination_robustness_arrays_v18_51x.npz`: prior X mixing curves and critical eta.
- `reference/parent_v18_51t_run_audit.py`, `reference/parent_v18_51w_run_audit.py`, `reference/parent_v18_51x_run_audit.py` and matching reports are preserved for functional genealogy.
- `scripts/run_audit.py`: new v18.51y generating script with strict-source-support and chart-covariance checks.
- `tests/`: a synthetic unit-level model with negative controls, and scientific checks on real archive data.
- This stage does not repeat any expensive physical IVP. It examines the well-formed recorded tangent transport and its coordinate meaning.
