#!/usr/bin/env python3
"""Check exact source-file SHA-256, including the authenticated ancestor ZIP."""
import hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
entry=json.loads((r/'SOURCE_MANIFEST.json').read_text())['files']
for rel,correct in entry.items():
    p=r/rel
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=correct:
        raise RuntimeError('SOURCE_MANIFEST_FAILURE '+rel)
print('SOURCE_MANIFEST_PASS',len(entry),'files')
