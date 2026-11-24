#!/usr/bin/env python3
"""Pack small audited source/output return, excluding large already-known parents."""
import zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(root.rglob('*')):
  if not f.is_file() or '__pycache__' in f.parts:continue
  if 'inputs' in f.parts:continue
  if 'reference' in f.parts:continue
  z.write(f,arcname=root.name+'/'+f.relative_to(root).as_posix())
print('Q550 v18.52b sendback:',out,flush=True)
