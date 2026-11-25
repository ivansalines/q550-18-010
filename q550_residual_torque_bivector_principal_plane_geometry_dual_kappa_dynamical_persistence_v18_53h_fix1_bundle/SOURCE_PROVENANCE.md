# Genealogy / provenance

- Historical original v18.51s nested generator sendback: exact SHA-256 `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b` (external server input, required for new dynamics).
- Native v18.52l 10D trajectories: exact SHA-256 `92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12` (external server input, required).
- v18.53f: pinned generator-stage loader, kinetic embedding, source force decomposition and authenticated report/manifest, shipped in `reference/` and `scripts/`.
- v18.53g: authenticated full server `report.json`, integrated bivector ledger, parameter-plane projections and output manifest shipped in `reference/`, individually pinned and cross-verified.
- Independent generated code: `scripts/bivector_geometry.py` real Schur 5D principal-plane geometry, `scripts/run_audit.py` full source-driven historical-clock IVPs, `tests/test_synthetic.py`, manifest validation and packaging scripts.
- This is a *functional extension* of the 53g scientific operation, keeping all keys for SHA-verifiable ancestors; scripts are included rather than results-only archives.
- The generator uses the historical drive and original stage chain. Large external antecedent bytes are not duplicated in every ZIP, but hashes and explicit paths are always provided.

## v18.53h fix1 provenance

- Recovered the actual server failure from the submitted v18.53h sendback:
  `NameError: name 'clock_error' is not defined` at `scripts/run_audit.py:125`.
- Original report status `AUDIT_FAIL_NOT_CERTIFIED`; `n_new_Q550_ivps=0`.
- This patch computes the material clock derivative at full native resolution and
  adds synthetic positive/negative tests; no historical reference sources changed.
- All old results from the failed sendback are excluded from this reissued input ZIP.
- `n_new_Q550_ivps` now counts successfully completed IVPs incrementally, including refinements, so interrupted audits never overclaim or under-report the work done.
