#!/usr/bin/env python3
"""Pinned physical ancestry and prior audited interpretations, no user guesswork."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
l=json.loads((ROOT/'reference/PARENT_v18_52l_OUTPUT_MANIFEST.json').read_text())
m=json.loads((ROOT/'reference/PARENT_v18_52m_report.json').read_text())
n=json.loads((ROOT/'reference/PARENT_v18_52n_report.json').read_text())
manifestn=json.loads((ROOT/'reference/PARENT_v18_52n_OUTPUT_MANIFEST.json').read_text())
assert l['version']=='v18.52l' and m['version']=='v18.52m' and n['version']=='v18.52n'
assert manifestn['version']=='v18.52n'
sha=l['file_sha256']['results/extended_horizon_5d_vectors_and_traces.npz']
assert sha=='d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195'
assert m['input_sha256']==sha
assert abs(n['observed_m_reference']+.308)<1e-14
assert abs(n['theta_deg']-17.04)<1e-14
print('LINEAGE PASS: v18.52l physical data; v18.52m/v18.52n observer audits')
