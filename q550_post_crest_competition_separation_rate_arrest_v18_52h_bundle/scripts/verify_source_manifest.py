#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
assert m['version']=='v18.52h'
for name,expected in m['sha256'].items():
    digest=hashlib.sha256((root/name).read_bytes()).hexdigest()
    if digest!=expected:
        raise RuntimeError(f'Source file SHA-256 mismatch: {name}')
print('SOURCE MANIFEST PASS',len(m['sha256']),'files')
