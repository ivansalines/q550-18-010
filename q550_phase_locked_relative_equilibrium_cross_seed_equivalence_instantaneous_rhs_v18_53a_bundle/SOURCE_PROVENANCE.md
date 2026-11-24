# Source provenance

- **Directly authenticated historical report:** `reference/v18_52z_authenticated_report.json` copied verbatim from user server sendback, SHA-256 pinned in `inputs_config.json`. It contains 8D final Newton candidates, J, period, 3-turn certificate, conditional Floquet eigenspectra, model kappa and scientific boundaries.
- **Original 11D physical source:** historical `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py` and `reference/v18_52v_r1_run_audit.py` + `reference/parent_v18_51t_run_audit.py` source loaders, with complete staged v18.51s server sendback SHA-verified when authentic runs occur.
- **Native state:** v18.52l `extended_horizon_10d_states.npz`, SHA-256 pinned and loaded for the gauge reference and time anchor; the 11th NPZ column is accumulated material time, not Theta.
- **Retained executable mathematics:** source-faithful fixed-J chart `scripts/return_map.py` and energy scripts from v18.52z, plus newly written independent algebraic stationary root and instantaneous native RHS verifier.
- **Local preflight:** computes independent algebraic fixed-J roots and the predicted linear spectrum using reconstructed equations. It does **not** stage the inaccessible historical 19-ancestor engine, and cannot certify the original instantaneous RHS, phase offset or full cross-seed equivalence.
