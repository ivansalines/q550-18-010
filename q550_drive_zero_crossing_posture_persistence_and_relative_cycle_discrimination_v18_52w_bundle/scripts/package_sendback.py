#!/usr/bin/env python3
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1];out=root/'results';out.mkdir(exist_ok=True)
output_files={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.52w','root':root.name,'files':output_files},indent=2)+'\n')
p=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(root.rglob('*')):
  if f.is_file():z.write(f,arcname=root.name+'/'+str(f.relative_to(root)))
print('SENDBACK',p,'outputs',len(output_files))
