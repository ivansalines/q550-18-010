#!/usr/bin/env python3
"""Create a full, executable, self-contained Q550 server sendback archive."""
from pathlib import Path
import zipfile
r=Path(__file__).resolve().parents[1]
target=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
 for p in sorted(r.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:
   z.write(p,arcname=r.name+'/'+str(p.relative_to(r)))
with zipfile.ZipFile(target) as z: assert z.testzip() is None
print('Q550 v18.51w sendback:',target,'bytes:',target.stat().st_size)
