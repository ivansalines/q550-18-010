#!/usr/bin/env python3
"""Check locally retained v18.52l/v18.52m output histories, without remote calls."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
l=json.loads((ROOT/'reference/PARENT_v18_52l_OUTPUT_MANIFEST.json').read_text())
m=json.loads((ROOT/'reference/PARENT_v18_52m_OUTPUT_MANIFEST.json').read_text())
r=json.loads((ROOT/'reference/PARENT_v18_52m_report.json').read_text())
assert l['version']=='v18.52l' and r['version']=='v18.52m'
sha=l['file_sha256']['results/extended_horizon_5d_vectors_and_traces.npz']
assert sha=='d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195'
assert r['input_sha256']==sha
assert len(r['PCA_hidden_axis_Q3_Q4_Q5'])==3
assert len(m['file_sha256'])>=8
print('PARENT LINEAGE PASS: v18.52l physical vectors and v18.52m server report pinned')
