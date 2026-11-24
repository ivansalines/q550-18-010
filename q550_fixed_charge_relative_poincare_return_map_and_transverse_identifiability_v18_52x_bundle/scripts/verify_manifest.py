#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for name,expected in m['files'].items():
 p=root/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
  print('SOURCE_MANIFEST_FAIL',name,flush=True);sys.exit(2)
print('SOURCE_MANIFEST_PASS',len(m['files']))
