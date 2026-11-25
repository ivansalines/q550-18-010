#!/usr/bin/env python3
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'files':{str(p.relative_to(root)):sha(p) for p in sorted((root/'results').rglob('*')) if p.is_file() and p.name!='.keep'}}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
out=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and p!=out:z.write(p,str(root.name+'/'+p.relative_to(root)))
print('SEND_BACK_CREATED',out,'OUTPUT_HASHES',len(manifest['files']))
