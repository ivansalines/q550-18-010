#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for name,expected in m['sha256'].items():
    p=root/name
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
        raise SystemExit('Source file missing or drifted: '+name)
print('PASS: source provenance hashes',len(m['sha256']))
