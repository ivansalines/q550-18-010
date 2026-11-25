#!/usr/bin/env python3
"""Verify immutable source manifest entries before any execution."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
entries=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
for rel,digest in entries.items():
    path=root/rel
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
        raise RuntimeError('SOURCE_MANIFEST_MISMATCH '+rel)
print('V1853S_SOURCE_MANIFEST_PASS',len(entries),'files')
