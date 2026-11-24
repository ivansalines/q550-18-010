#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
files={}
for f in sorted(ROOT.rglob('*')):
 if f.is_file() and f.name!='OUTPUT_MANIFEST.json' and '__pycache__' not in f.parts:
  with f.open('rb') as stream:files[str(f.relative_to(ROOT))]=hashlib.file_digest(stream,'sha256').hexdigest()
(ROOT/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.52p','files_sha256':files},indent=2,sort_keys=True)+'\n')
print('OUTPUT MANIFEST PASS',len(files))
