# Source genealogy — Q550 v18.52z

1. The immediately preceding **v18.52y server sendback** supplies *both* source-authenticated eight-dimensional Newton candidate vectors, full search history, fixed charge, original phase, and native replay checks. Its `results/report.json` is retained verbatim as `reference/v18_52y_authenticated_report.json`, SHA-256 pinned independently in `inputs_config.json`.
2. The identical historical `reference/v18_52y_run_audit.py`, `reference/v18_52v_r1_run_audit.py`, `reference/parent_v18_51t_run_audit.py`, and original support source are retained for function-level provenance. No other version's model is silently substituted.
3. The authentic full historical **v18.51s sendback** remains an external SHA-pinned input; it should be in `/root/elab1/hysteresis`. It is not reproduced in this small bundle because of its large nested ancestry.
4. The native **v18.52l 10D state file** likewise remains external, SHA pinned. It supplies the first original native reference state and independent one-turn endpoint, not a synthetic cycle.
5. The new `scripts/run_audit.py` implements source-authenticated Newton continuation from the preserved candidate, independent 1/2/3-turn closure and conditional Floquet and energy checks. `scripts/return_map.py`, `energy_math.py` and `phase_energy.py` are retained functional Python sources, not merely results.

The new audit performs genuine native 10D IVPs on the user's server only. File names and SHA-256 in the source/output manifests provide executable ancestry. Full paths may differ; functional source identity is enforced by hashes.
