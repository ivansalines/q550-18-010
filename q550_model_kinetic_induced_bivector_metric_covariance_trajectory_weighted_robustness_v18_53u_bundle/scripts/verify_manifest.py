#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
a=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
for path,expect in a.items():
 p=root/path
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expect:
  raise RuntimeError('SOURCE_HASH_MISMATCH '+path)
print('SOURCE_MANIFEST_PASS',len(a),'files')
