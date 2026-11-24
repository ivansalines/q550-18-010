#!/usr/bin/env python3
"""Strict source file SHA256 verification."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
assert m['root']==root.name and m['version']=='v18.51t'
for name,h in m['files'].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,name
print('Q550 v18.51t source SHA256 verified',len(m['files']))
