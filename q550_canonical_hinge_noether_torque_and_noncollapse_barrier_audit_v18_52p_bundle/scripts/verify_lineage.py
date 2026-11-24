#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
a=json.loads((ROOT/'reference/PARENT_v18_52l_OUTPUT_MANIFEST.json').read_text())
b=json.loads((ROOT/'reference/PARENT_v18_52o_report.json').read_text())
c=json.loads((ROOT/'reference/PARENT_v18_52o_OUTPUT_MANIFEST.json').read_text())
assert a['version']=='v18.52l' and b['version']=='v18.52o' and c['version']=='v18.52o'
assert a['file_sha256']['results/extended_horizon_5d_vectors_and_traces.npz']=='d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195'
assert a['file_sha256']['results/extended_horizon_10d_states.npz']=='92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12'
assert b['source_sha256']==a['file_sha256']['results/extended_horizon_5d_vectors_and_traces.npz']
assert b['stage']=='AUTHENTIC_FULL_RESOLUTION_POSTPROCESS'
assert abs(b['nodal_axis_slope']+.308)<1e-14
assert hashlib.sha256((ROOT/'reference/kinetic_metric_v18_49h.py').read_bytes()).hexdigest()=='29bef2c91899a9da46bc8c9159b89c70362f505e0adead3f814c487a9a17bb1a'
# Check every file that was sealed before the scientific server execution.
source=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
assert source['version']=='v18.52p'
for rel,digest in source['sha256'].items():
 f=ROOT/rel
 assert f.is_file(),f'Missing sealed source: {rel}'
 assert hashlib.sha256(f.read_bytes()).hexdigest()==digest,f'Source drift: {rel}'
print('SEALED SOURCE PASS',len(source['sha256']))
print('LINEAGE PASS: authentic v18.52l native 5D/10D ancestor; v18.52o observer; v18.49h kinetic embedding')
