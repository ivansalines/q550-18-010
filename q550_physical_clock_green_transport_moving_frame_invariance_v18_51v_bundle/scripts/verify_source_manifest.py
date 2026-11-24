#!/usr/bin/env python3
"""Verify source SHA256 entries, including the complete U ancestor sendback."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
assert m['version']=='v18.51v' and m['root']==root.name
for rel,sha in m['files'].items():
    p=root/rel
    assert p.exists(),f'Missing source: {rel}'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==sha,f'Source hash mismatch: {rel}'
print('Q550 v18.51v verified',len(m['files']),'source files')
