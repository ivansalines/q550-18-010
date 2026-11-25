#!/usr/bin/env python3
"""Verify v18.53l source immutability before every server run."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for f,digest in m['files'].items():
    p=root/f
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
        raise SystemExit('V1853L_SOURCE_MANIFEST_FAIL '+f)
print('V1853L_SOURCE_MANIFEST_PASS',len(m['files']))
