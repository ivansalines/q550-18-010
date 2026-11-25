#!/usr/bin/env python3
"""Save complete executed source+outputs+ancestors in one full-name ZIP."""
from pathlib import Path
import hashlib,json,zipfile
root=Path(__file__).resolve().parents[1]
files=sorted((p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts
 and p.suffix!='.pyc' and p.name!='OUTPUT_MANIFEST.json'),
 key=lambda p:str(p.relative_to(root)))
m={'version':'v18.53m','root_folder':root.name,
  'files':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(m,indent=2)+'\n')
files.append(root/'OUTPUT_MANIFEST.json')
out=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:
    for f in files:z.write(f,arcname=root.name+'/'+str(f.relative_to(root)))
print('V1853M_SENDBACK',out,'MEMBERS',len(files))
