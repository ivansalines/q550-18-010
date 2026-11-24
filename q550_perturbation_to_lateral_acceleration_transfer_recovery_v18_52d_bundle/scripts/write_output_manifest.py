#!/usr/bin/env python3
"""Hash all generated results and server report, never only preflight outputs."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
res={}
for f in sorted((ROOT/'results').rglob('*')):
    if f.is_file():
        rel=f.relative_to(ROOT).as_posix()
        with f.open('rb') as obj:res[rel]=hashlib.file_digest(obj,'sha256').hexdigest()
manifest={'version':'v18.52d','root':ROOT.name,'file_sha256':res}
(ROOT/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Q550 v18.52d output hashes:',len(res))
