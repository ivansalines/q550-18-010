#!/usr/bin/env python3
"""Package all scripts, archived provenance, and generated results in sendback."""
from pathlib import Path
import zipfile
r=Path(__file__).resolve().parents[1]
out=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1,allowZip64=True) as z:
    for p in sorted(r.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            z.write(p,p.relative_to(r.parent))
with zipfile.ZipFile(out) as z:
    assert z.testzip() is None
print('Q550 v18.51t sendback',out,'bytes',out.stat().st_size)
