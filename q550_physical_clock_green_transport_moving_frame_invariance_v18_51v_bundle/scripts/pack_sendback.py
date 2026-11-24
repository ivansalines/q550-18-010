#!/usr/bin/env python3
"""Preserve the exact full scientific bundle root and genealogy in sendback."""
from pathlib import Path
import zipfile
root=Path(__file__).resolve().parents[1]
target=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(target,'w',allowZip64=True) as z:
    for p in sorted(root.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        # A nested ZIP was already compressed; avoid expensively recompressing it.
        compression=zipfile.ZIP_STORED if p.suffix=='.zip' else zipfile.ZIP_DEFLATED
        z.write(p,arcname=str(p.relative_to(root.parent)),compress_type=compression,compresslevel=6)
with zipfile.ZipFile(target) as z:
    bad=z.testzip()
    if bad:raise RuntimeError('Sendback CRC failure: '+bad)
print('Q550 v18.51v sendback:',target)
