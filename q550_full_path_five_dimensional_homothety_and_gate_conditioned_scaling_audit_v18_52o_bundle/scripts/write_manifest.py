#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
a={}
for f in sorted(ROOT.rglob('*')):
 if f.is_file() and f.name!='OUTPUT_MANIFEST.json' and '__pycache__' not in f.parts:
  with f.open('rb') as h:a[str(f.relative_to(ROOT))]=hashlib.file_digest(h,'sha256').hexdigest()
(ROOT/'OUTPUT_MANIFEST.json').write_text(json.dumps(dict(version='v18.52o',files_sha256=a),indent=2)+'\n')
print('MANIFEST PASS',len(a))
