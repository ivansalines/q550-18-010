#!/usr/bin/env python3
"""Verify immutable source and embedded t-parent SHA256."""
import hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
assert m['version']=='v18.51u' and m['root']==r.name
for path,h in m['files'].items():
    p=r/path
    assert p.is_file(),f'Missing source dependency {path}'
    with p.open('rb') as f:
        calc=hashlib.file_digest(f,'sha256').hexdigest()
    assert calc==h,f'Source SHA256 mismatch: {path}'
print('Q550 v18.51u source verified',len(m['files']))
