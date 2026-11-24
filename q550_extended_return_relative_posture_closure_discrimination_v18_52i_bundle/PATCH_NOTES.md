# v18.52i Patch 1 — Ledger h-seam angle read alias

Original failure: `KeyError: 'angle_from_h_seam_deg'` at the CSV writer.
The numerical traces store `angle_from_seam_deg`, while the public CSV names
its column `angle_from_h_seam_deg`. Patched the read through `ledger_value`;
the CSV schema and the in-memory plotting key remain unchanged.

- No changes to nonlinear 10D dynamics, material clock, kappa or 5D geometry.
- New smoke test: `tests/test_ledger_hotfix.py`.
- Original parent h and nested historical archives remain byte-for-byte intact.
- `SOURCE_MANIFEST.json` refreshed for changed sources.
- Existing two NPZ files from failed attempt are partial results, not a full sendback.
