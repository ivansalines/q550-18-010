# Provenance — Q550 v18.52n

- Original model generator, ODE integration and historical trajectory construction: inherited from v18.52l; **not recreated, altered or claimed as part of v18.52n**.
- Physical 5D NPZ source: v18.52l authentic server output, SHA-256 `d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195`.
- Paired observer original experiment: v18.52m authentic server sendback report and output manifest pinned in `reference/PARENT_v18_52m_*`; prior Q1-axis observer used m=0.
- This v18.52n source: `scripts/run_audit.py` (full geometry, diagnostics and figures), `tests/test_synthetic.py`, `scripts/verify_lineage.py`, `scripts/write_manifest.py`, `scripts/pack_sendback.py`.
- User-specified: `m≈-0.308`, nominal visible direction `atan(m)≈-17.119°`; observer tilt pair `±17.04°`. These are **not** inferred from the v18.52l trace without the actual four observational pairs.
- No historical script is silently replaced, no source data modified, no new ODE integration.
- `inputs/four_point_template.csv` is deliberately **header-only**: we do not invent the coordinates of four user-observed points.
