#!/usr/bin/env python3
"""Authenticate the distributable files of v18.53y; excludes computed results."""
from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
for file,digest in manifest.items():
    p=root/file
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
        raise RuntimeError('SOURCE_MANIFEST_HASH_MISMATCH '+file)
print('SOURCE_MANIFEST_PASS',len(manifest),'files')
