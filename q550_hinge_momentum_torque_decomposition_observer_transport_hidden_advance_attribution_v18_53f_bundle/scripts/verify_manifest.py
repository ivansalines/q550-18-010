#!/usr/bin/env python3
"""Fail-closed SHA-256 verification of all portable execution inputs."""
import json,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
files=manifest['files']
for rel,digest in files.items():
 p=ROOT/rel
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
  raise SystemExit('SOURCE_MANIFEST_FAIL: '+rel)
print('SOURCE_MANIFEST_PASS',len(files))
