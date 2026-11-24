#!/usr/bin/env python3
import json,zipfile
from pathlib import Path
r=Path(__file__).resolve().parents[1]
out=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
 for f in sorted(r.rglob('*')):
  if f.is_file() and '__pycache__' not in f.parts:z.write(f,arcname=r.name+'/'+f.relative_to(r).as_posix())
with zipfile.ZipFile(out) as z:
 assert z.testzip() is None
print('SENDBACK PACKED',out,'bytes',out.stat().st_size)
