#!/usr/bin/env python3
"""Verify every pinned v18.53v source, including authenticated parent bytes."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
for rel,expected in manifest.items():
    p=root/rel
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
        raise RuntimeError('SOURCE_HASH_MISMATCH '+rel)
print('SOURCE_MANIFEST_PASS',len(manifest),'files')
