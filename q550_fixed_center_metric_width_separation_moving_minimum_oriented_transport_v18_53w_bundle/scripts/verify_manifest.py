#!/usr/bin/env python3
"""Verify the immutable distribution source manifest."""
import hashlib,json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
listing=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
for name,digest in listing.items():
    p=root/name
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
        raise SystemExit('SOURCE_MANIFEST_MISMATCH '+name)
print('SOURCE_MANIFEST_PASS',len(listing),'files')
