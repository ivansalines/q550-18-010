#!/usr/bin/env python3
"""Build a compact, reproducible server sendback excluding trusted large parents."""
from pathlib import Path
import zipfile
root=Path(__file__).resolve().parents[1]
out=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:
    for f in sorted(root.rglob('*')):
        if not f.is_file() or '__pycache__' in f.parts or 'inputs' in f.parts:continue
        z.write(f,arcname=root.name+'/'+f.relative_to(root).as_posix())
print('Q550 v18.52c sendback:',out,flush=True)
