#!/usr/bin/env python3
"""Make genealogically named Q550 v18.51y server sendback."""
from pathlib import Path
import zipfile
r=Path(__file__).resolve().parents[1]
output=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(output,'w',allowZip64=True) as z:
    for p in sorted(r.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            z.write(p,arcname=str(Path(r.name)/p.relative_to(r)),
                    compress_type=zipfile.ZIP_STORED if p.suffix in ('.zip','.npz','.png') else zipfile.ZIP_DEFLATED)
print('Q550 v18.51y SENDBACK',output)
