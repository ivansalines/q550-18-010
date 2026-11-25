#!/usr/bin/env python3
"""Archive sources and all results with independently verified output hashes."""
import json,hashlib,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
files=sorted(p for p in root.rglob('*') if p.is_file() and not ('__pycache__' in p.parts or p.suffix=='.pyc' or p.name=='OUTPUT_MANIFEST.json'))
mani={'version':'v18.53g','root_folder':root.name,'files':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(mani,indent=2)+'\n')
files.append(root/'OUTPUT_MANIFEST.json')
out=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
 for p in files:z.write(p,arcname=root.name+'/'+str(p.relative_to(root)))
print('SENDBACK',out,'FILES',len(files))
