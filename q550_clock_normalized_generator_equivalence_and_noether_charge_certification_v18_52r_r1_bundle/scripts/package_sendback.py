#!/usr/bin/env python3
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root.parent/(root.name+'_server_sendback.zip')
sha=lambda p: hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
files=[p for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts]
output_checks={str(p.relative_to(root)):sha(p) for p in files if 'results'==p.relative_to(root).parts[0]}
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'version':'v18.52r','root':root.name,'file_sha256':output_checks},indent=2)+'\n')
files=[p for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts]
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,arcname=root.name+'/'+str(p.relative_to(root)))
with zipfile.ZipFile(out) as z:assert z.testzip() is None
print('SENDBACK:',out,'bytes:',out.stat().st_size)
