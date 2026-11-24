#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text());n=0
for name,want in m['files'].items():
 p=root/name
 if not p.is_file():raise RuntimeError('MISSING '+name)
 got=hashlib.sha256(p.read_bytes()).hexdigest()
 if got!=want:raise RuntimeError('SHA256_MISMATCH '+name)
 n+=1
print('SOURCE_MANIFEST_PASS',n)
