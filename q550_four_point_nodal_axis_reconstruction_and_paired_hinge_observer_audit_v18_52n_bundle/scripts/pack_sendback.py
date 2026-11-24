#!/usr/bin/env python3
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
ROOT=Path(__file__).resolve().parents[1]
dst=ROOT.parent/(ROOT.name+'_server_sendback.zip')
with ZipFile(dst,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(ROOT.rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts:continue
  z.write(p,arcname=ROOT.name+'/'+p.relative_to(ROOT).as_posix())
with ZipFile(dst) as z: assert z.testzip() is None
print('SENDBACK ZIP PASS',dst, dst.stat().st_size,'bytes')
