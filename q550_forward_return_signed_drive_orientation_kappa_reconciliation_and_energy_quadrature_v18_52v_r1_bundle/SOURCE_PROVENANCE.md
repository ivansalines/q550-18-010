# Provenance — v18.52v-r1

Source references preserved without modification:
- Archived v18.52v `scripts/run_audit.py` copied to `reference/v18_52v_original_run_audit.py`; its original server `report.json` and native ledger copied under `reference` for explicit before/after.
- v18.52t authenticated 20-series baseline CSV copied to `reference/v18_52t_authenticated_20_series_baseline.csv`.
- Historical v18.51s engine staging code retained in `reference/parent_v18_51t_run_audit.py` (same as previous bundle).
- Historical scalar RHS retained in `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py` (same as previous bundle).
- Original v18.51s ZIP and v18.52l NPZ are external SHA-256 pinned inputs; not duplicated because they are large and live on the user's server.

Novel r1 scripts are functional reconstructions of the v18.52v audit with two explicit scientific corrections: the source sign and the native shared kappa. All calculation scripts and data formats are in this ZIP. Manifest hashes pin every distributed source.
