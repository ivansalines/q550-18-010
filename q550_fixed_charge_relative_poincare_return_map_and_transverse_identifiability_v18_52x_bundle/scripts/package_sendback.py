#!/usr/bin/env python3
import hashlib,json,zipfile
from pathlib import Path
r=Path(__file__).resolve().parents[1];out=r/'results';out.mkdir(exist_ok=True)
files={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}
(r/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.52x','root':r.name,'files':files},indent=2)+'\n')
p=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(r.rglob('*')):
  if f.is_file() and '__pycache__' not in f.parts:z.write(f,arcname=r.name+'/'+str(f.relative_to(r)))
print('SENDBACK',p,'results',len(files))
