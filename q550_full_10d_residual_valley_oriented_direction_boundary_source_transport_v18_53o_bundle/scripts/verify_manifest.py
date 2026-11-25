#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
for rel,expected in m.items():
    p=root/rel
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
        raise SystemExit('SOURCE_MANIFEST_FAIL '+rel)
print('SOURCE_MANIFEST_PASS',len(m),'files')
