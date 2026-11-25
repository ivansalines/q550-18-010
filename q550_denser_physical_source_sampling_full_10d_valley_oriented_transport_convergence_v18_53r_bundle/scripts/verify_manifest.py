#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
r=Path(__file__).resolve().parents[1]
m=json.loads((r/'SOURCE_MANIFEST.json').read_text())['files']
for name,h in m.items():
 p=r/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=h:
  raise SystemExit('SOURCE_MANIFEST_FAILED '+name)
print('SOURCE_MANIFEST_PASS',len(m),'files',flush=True)
