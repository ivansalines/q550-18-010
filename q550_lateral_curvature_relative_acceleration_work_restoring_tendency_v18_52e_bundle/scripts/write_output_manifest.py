#!/usr/bin/env python3
"""Hash every produced result; do not reuse incomplete preliminary manifests."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
items={}
for p in sorted((ROOT/'results').rglob('*')):
    if p.is_file() and p.name!='OUTPUT_MANIFEST.json':
        with open(p,'rb') as f:items[p.relative_to(ROOT).as_posix()]=hashlib.file_digest(f,'sha256').hexdigest()
(ROOT/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.52e','file_sha256':items},indent=2)+'\n')
print('Q550 v18.52e results hashed:',len(items))
