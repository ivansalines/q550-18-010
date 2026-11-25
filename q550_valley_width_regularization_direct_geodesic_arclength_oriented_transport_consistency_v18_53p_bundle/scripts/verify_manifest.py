#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
for rel,dig in m['files'].items():
 p=root/rel
 if not p.is_file() or sha(p)!=dig:raise SystemExit('SOURCE_HASH_MISMATCH '+rel)
print('V1853P_SOURCE_MANIFEST_PASS',len(m['files']))
