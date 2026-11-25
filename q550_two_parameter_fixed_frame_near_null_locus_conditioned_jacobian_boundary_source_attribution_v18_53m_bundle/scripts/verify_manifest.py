#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for rel,digest in m['files'].items():
    p=root/rel
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
        raise SystemExit('V1853M_SOURCE_MANIFEST_FAIL '+rel)
print('V1853M_SOURCE_MANIFEST_PASS',len(m['files']))
