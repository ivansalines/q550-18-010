# Genealogy

v18.52v extends the authenticated results of v18.52u (phase-clock genealogy), v18.52t (signed energy ledger), and v18.52s-r2 (model-parameter reconciliation); it loads original `v18.51s` from the workspace, with the exact source archive hash specified in `inputs_config.json`.

Historical source code is **not replaced** by an approximation. The included `reference/parent_v18_51t_run_audit.py` stages original nested `v18.51s` -> `v18.51r` -> deep original engine after authenticating the outer SHA-256. The per-state source functions are inspected and compared with independently bundled energy formulas. The native v18.52l states are verified against the recorded exact SHA-256. All dynamically generated outputs receive SHA-256 in the server sendback.
