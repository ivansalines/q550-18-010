#!/usr/bin/env python3
"""Write v18.51y output hashes after completed server replay."""
from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1]
files=sorted(f for f in (r/'results').iterdir() if f.is_file())
hashes={}
for f in files:
    h=hashlib.sha256()
    with f.open('rb') as b:
        for chunk in iter(lambda:b.read(4*1024*1024),b''):h.update(chunk)
    hashes[str(f.relative_to(r))]=h.hexdigest()
(r/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'root':r.name,'version':'v18.51y','file_sha256':hashes},indent=2)+'\n')
print('Q550 v18.51y OUTPUT SHA256',len(hashes))
