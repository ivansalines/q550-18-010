#!/usr/bin/env python3
"""Check immutable source files, locked model input and external-parent lineage."""
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
m=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
assert m['version']=='v18.52a' and m['root']==ROOT.name
for name,h in m['files'].items():
    p=ROOT/name
    assert p.is_file(),f'Not found: {name}'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==h, f'Source drift: {name}'
from run_audit import find_input,sha256_file,SNAME,ZNAME
for path,tag in ((find_input(SNAME),'v18_51s_sha256'),(find_input(ZNAME),'v18_51z_sha256')):
    assert sha256_file(path)==m['parents'][tag],f'Parent mismatch {path}'
    with zipfile.ZipFile(path) as z:assert z.testzip() is None
print('Q550 v18.52a authenticated source files:',len(m['files']),'locked parents: 2')
