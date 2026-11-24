#!/usr/bin/env python3
"""Pack executable analysis genealogy and small derived outputs, NOT parent NPZ."""
import zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
archive=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
    for p in sorted(root.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.suffix=='.pyc': continue
        z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
with zipfile.ZipFile(archive) as z:
    bad=z.testzip()
    if bad:raise RuntimeError('Sendback zip CRC failed '+bad)
print('PACKED SERVER SENDBACK',archive)
