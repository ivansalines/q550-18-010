#!/usr/bin/env python3
"""Write a reproducible sendback with source and result SHA-256 hashes."""
import json,hashlib,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
files=sorted(p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='OUTPUT_MANIFEST.json')
mani={'version':'v18.53h','root_folder':root.name,'files':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(mani,indent=2)+'\n')
files.append(root/'OUTPUT_MANIFEST.json')
out=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
 for p in files:z.write(p,arcname=root.name+'/'+str(p.relative_to(root)))
print('SENDBACK',out,'FILES',len(files))
