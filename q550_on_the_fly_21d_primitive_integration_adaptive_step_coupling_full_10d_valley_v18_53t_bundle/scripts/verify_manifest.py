#!/usr/bin/env python3
"""Verify distribution source hashes; no generated results are trusted."""
import json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
files=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
for name,sha in files.items():
    p=root/name
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=sha:
        raise RuntimeError('SOURCE_HASH_FAILURE '+name)
print('V1853T_SOURCE_MANIFEST_PASS',len(files),'files')
