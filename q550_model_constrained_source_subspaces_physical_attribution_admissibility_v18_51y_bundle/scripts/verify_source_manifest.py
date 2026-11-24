#!/usr/bin/env python3
"""Check v18.51y SHA-256 snapshot integrity including optional full X parent."""
from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1]
m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
assert m['version']=='v18.51y' and m['root']==r.name

def sha(f):
    h=hashlib.sha256()
    with f.open('rb') as b:
        for chunk in iter(lambda:b.read(4*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
for name,expected in m['files'].items():
    f=r/name
    assert f.is_file() and sha(f)==expected,name
optional=m['optional_full_x_parent']
f=r/optional['path']
if f.exists():
    assert sha(f)==optional['sha256'],'Optional full parent X SHA256 mismatch'
print('Q550 v18.51y SOURCE VERIFIED',len(m['files']), 'optional full parent',f.is_file())
