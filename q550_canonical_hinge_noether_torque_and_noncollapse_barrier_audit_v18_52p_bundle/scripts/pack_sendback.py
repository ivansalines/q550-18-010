#!/usr/bin/env python3
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
ROOT=Path(__file__).resolve().parents[1]
output=ROOT.parent/(ROOT.name+'_server_sendback.zip')
with ZipFile(output,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(ROOT.rglob('*')):
  if f.is_file() and '__pycache__' not in f.parts:z.write(f,arcname=ROOT.name+'/'+f.relative_to(ROOT).as_posix())
with ZipFile(output) as z:assert z.testzip() is None
print('SENDBACK PACK PASS',output,output.stat().st_size)
