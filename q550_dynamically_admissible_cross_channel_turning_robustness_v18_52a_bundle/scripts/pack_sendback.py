#!/usr/bin/env python3
"""Pack result, scripts and provenance from actual executed input tree."""
from pathlib import Path
import zipfile
R=Path(__file__).resolve().parents[1]
P=R.parent/(R.name+'_server_sendback.zip')
with zipfile.ZipFile(P,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4) as z:
    for x in sorted(R.rglob('*')):
        if not x.is_file() or '__pycache__' in x.parts:continue
        z.write(x,R.name+'/'+x.relative_to(R).as_posix())
with zipfile.ZipFile(P) as z:assert z.testzip() is None
print('Created sendback',P)
