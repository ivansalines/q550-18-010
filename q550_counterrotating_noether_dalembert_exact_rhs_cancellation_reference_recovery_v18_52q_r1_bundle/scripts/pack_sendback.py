#!/usr/bin/env python3
import zipfile,hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1];zipdest=r.parent/(r.name+'_server_sendback.zip')
paths=sorted(p for p in r.rglob('*') if p.is_file() and not p.name.endswith(('.pyc',)) and '__pycache__' not in p.parts)
result={str(p.relative_to(r)):hashlib.file_digest(open(p,'rb'),'sha256').hexdigest() for p in paths if 'results' in p.parts}
(r/'OUTPUT_MANIFEST.json').write_text(json.dumps(dict(version='v18.52q',file_sha256=result),indent=2)+'\n')
with zipfile.ZipFile(zipdest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(p for p in r.rglob('*') if p.is_file() and '__pycache__' not in p.parts):z.write(p,arcname=r.name+'/'+str(p.relative_to(r)))
with zipfile.ZipFile(zipdest) as z:
 corrupt=z.testzip()
 if corrupt:raise RuntimeError('Archive CRC failed: '+corrupt)
print('SENDBACK READY:',zipdest,zipdest.stat().st_size,'bytes')
