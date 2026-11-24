#!/usr/bin/env python3
import zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
output=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and not p.name.endswith('.pyc'):
            z.write(p,arcname=root.name+'/'+str(p.relative_to(root)))
with zipfile.ZipFile(output) as z:
    assert z.testzip() is None
print('PACKED',output)
