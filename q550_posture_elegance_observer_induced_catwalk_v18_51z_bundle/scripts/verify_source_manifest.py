#!/usr/bin/env python3
"""Check exact input and source files against the fixed manifest."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for rel,expected in manifest['sha256'].items():
    f=root/rel
    assert f.is_file(),f'Missing: {rel}'
    actual=hashlib.sha256(f.read_bytes()).hexdigest()
    assert actual==expected, f'SHA256 mismatch: {rel}'
print('Q550 v18.51z verified source files',len(manifest['sha256']))
