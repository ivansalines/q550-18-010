#!/usr/bin/env python3
"""Verify all source archive hashes (including the authenticated V sendback)."""
import hashlib
import json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
for k,h in m['files'].items():
 p=r/k
 assert p.is_file(),f'Missing source: {k}'
 assert hashlib.sha256(p.read_bytes()).hexdigest()==h,f'Invalid source SHA-256: {k}'
print('Q550 v18.51w source hashes PASS',len(m['files']))
