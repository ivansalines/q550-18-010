#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
assert m['version']=='v18.52l'
for name,digest in m['sha256'].items():
    got=hashlib.sha256((root/name).read_bytes()).hexdigest()
    if got!=digest:raise RuntimeError(f'Source SHA mismatch: {name}')
print('SOURCE MANIFEST PASS',len(m['sha256']),'files')
