#!/usr/bin/env python3
"""Pack fully reproducible v18.51u sendback including authenticated t ancestor."""
from pathlib import Path
import zipfile
r=Path(__file__).resolve().parents[1]
p=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1,allowZip64=True) as z:
    for f in sorted(r.rglob('*')):
        if not f.is_file() or '__pycache__' in f.parts:continue
        compress=zipfile.ZIP_STORED if f.suffix=='.zip' else zipfile.ZIP_DEFLATED
        z.write(f,f.relative_to(r.parent),compress_type=compress)
with zipfile.ZipFile(p) as z:assert z.testzip() is None
print('Q550 v18.51u sendback',p,'bytes',p.stat().st_size)
