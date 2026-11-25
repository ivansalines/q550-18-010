#!/usr/bin/env python3
"""Generate a non-forgeable-by-filename source/output ledger and full sendback."""
from pathlib import Path
import zipfile,hashlib,json
root=Path(__file__).resolve().parents[1]
files=sorted(p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='OUTPUT_MANIFEST.json')
ledger={'version':'v18.53j','root_folder':root.name,'files':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(ledger,indent=2)+'\n')
files.append(root/'OUTPUT_MANIFEST.json')
target=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
 for p in files:z.write(p,arcname=root.name+'/'+str(p.relative_to(root)))
print('SENDBACK',target,'FILES',len(files))
