#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
r=Path(__file__).resolve().parents[1];m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
for fn,h in m['files'].items():
 p=r/fn
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise RuntimeError('SOURCE_HASH_MISMATCH '+fn)
print('SOURCE_MANIFEST_PASS',len(m['files']))
