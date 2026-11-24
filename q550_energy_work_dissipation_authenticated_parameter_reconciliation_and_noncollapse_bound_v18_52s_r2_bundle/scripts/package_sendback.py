#!/usr/bin/env python3
"""Package results without source NPZ; include full executable scripts."""
import zipfile,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
archive=root.parent/(root.name+'_server_sendback.zip')
paths=[p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc']
sha={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.52s-r2','sha256':sha},indent=2)+'\n')
paths.append(root/'OUTPUT_MANIFEST.json')
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:
 for p in sorted(paths):z.write(p,str(p.relative_to(root.parent)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
print('SENDBACK',archive,archive.stat().st_size)
