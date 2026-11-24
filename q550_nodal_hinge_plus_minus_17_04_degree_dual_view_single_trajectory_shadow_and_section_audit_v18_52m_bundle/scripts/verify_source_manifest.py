#!/usr/bin/env python3
"""Verify the exact Q550 v18.52m analysis source bundle before running."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text())
assert manifest['version']=='v18.52m'
for rel,digest in manifest['sha256'].items():
    p=root/rel
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
        raise RuntimeError(f'SOURCE_MANIFEST mismatch: {rel}')
print('SOURCE MANIFEST PASS',len(manifest['sha256']),'files')
