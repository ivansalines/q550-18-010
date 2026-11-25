#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
files=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
for rel,dig in files.items():
 p=root/rel
 if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=dig:raise SystemExit('V1853Q_HASH_FAIL '+rel)
print('V1853Q_SOURCE_MANIFEST_PASS',len(files))
