# Source genealogy — v18.52x

- `reference/v18_52w_run_audit.py`: actual predecessor v18.52w engine/ancestry loader and authentic drive-zero/one-turn scout logic (included **as executable source**, not just a report).
- `reference/v18_52v_r1_run_audit.py`, `reference/parent_v18_51t_run_audit.py`, `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py`: original engine staging and the historical Q550 physical RHS required by predecessor.
- `scripts/energy_math.py`, `scripts/phase_energy.py`: needed source-level dependency modules of the historical engine loader.
- `reference/v18_52w_authenticated_report.json`: exact previous server report; its SHA256 is recorded in `inputs_config.json`.
- `reference/v18_52w_inputs_config.json`: historical source parameters, SHA pins and protocol.
- `inputs_config.json`: new scientific protocol, pinned upstream hashes and output policies.
- `scripts/run_audit.py`: this audit's executable implementation.
- `tests/test_synthetic.py`: independent gauge/charge/rank/false-recurrence falsifiers.
- `SOURCE_MANIFEST.json`: SHA256 for every bundled reproducibility source; `verify_manifest.py` checks it before the real audit.
- `OUTPUT_MANIFEST.json`: generated at runtime with SHA256 for every output, packaged together with the scripts in the server sendback.

External large inputs are intentionally **not copied** into the ~100 KB runner bundle to avoid retransferring an authentic 97 MB native dataset. They are recovered from original filenames and verified against SHA pins: the v18.51s sendback and the v18.52l native NPZ. This is reproducibility via exact-byte pinned external inputs, not omission of provenance.

No result claimed to arise from fresh numerical integrations. All work uses the established authentic trajectories.
