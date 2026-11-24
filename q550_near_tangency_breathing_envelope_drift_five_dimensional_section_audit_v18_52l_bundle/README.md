# Q550 — v18.52l: Near-Tangency, Breathing-Envelope Drift & Five-Dimensional Section Audit

Authentic nonlinear 10D continuation from v18.52k (global xi=12.8) to xi=25.6, with no parameter or kappa changes. Physical 5D horizontal section Q2=0 is sampled with nonlinear dense evaluation, not inferred from projected pixels. Full 5D crossing pairs quantify deviations from exact reversal or mirror-like recurrence; they are not closure certificates.

Run on the research server:

```bash
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh --xi-end 25.6 --samples 65537 --refine
```

Bundled outputs: all historical 10D restart states, 5D scalar and vector diagnostics, horizontal_section_crossings_5d.csv, horizontal_section_pairs_5d.csv, horizontal_section_report.json, horizontal_section_projected_5d.png, figures, report and hashes. The k sendback and historical generator are included byte-for-byte. Requires the historical v18.51s external source path as in the prior audits.
