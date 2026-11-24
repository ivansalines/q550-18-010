#!/usr/bin/env python3
import zipfile
from pathlib import Path
r=Path(__file__).resolve().parents[1];o=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(o,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
 for p in sorted(r.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:
   z.write(p,arcname=r.name+'/'+p.relative_to(r).as_posix())
with zipfile.ZipFile(o) as z:assert z.testzip() is None
print('SENDBACK PACKED',o,'bytes',o.stat().st_size)
