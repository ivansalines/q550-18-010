#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for rel,h in m['files'].items():
 p=root/rel
 if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise SystemExit('SOURCE_SHA_MISMATCH: '+rel)
print('SOURCE_MANIFEST_PASS',len(m['files']))
