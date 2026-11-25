#!/usr/bin/env python3
"""Validate the exact bundled source files before doing any physics."""
import json,hashlib
from pathlib import Path
base=Path(__file__).resolve().parents[1]
files=json.loads((base/'SOURCE_MANIFEST.json').read_text())['files']
for name,digest in files.items():
    p=base/name
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
        raise SystemExit('SOURCE_MANIFEST_FAILED: '+name)
print('SOURCE_MANIFEST_PASS',len(files),'files')
