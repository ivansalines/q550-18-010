#!/usr/bin/env python3
"""Make a self-contained sendback with sources, provenance and results."""
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
files=sorted([p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='OUTPUT_MANIFEST.json'])
manifest={'version':'v18.53f','root_folder':ROOT.name,'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(ROOT/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
files.append(ROOT/'OUTPUT_MANIFEST.json')
out=ROOT.parent/(ROOT.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
 for p in sorted(files):z.write(p,arcname=ROOT.name+'/'+str(p.relative_to(ROOT)))
print('SENDBACK',out,'FILES',len(files))
