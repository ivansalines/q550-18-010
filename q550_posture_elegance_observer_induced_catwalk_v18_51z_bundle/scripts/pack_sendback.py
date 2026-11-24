#!/usr/bin/env python3
"""Create sendback preserving the exact descriptive root and source genealogy."""
import zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
name=root.name+'_server_sendback.zip'
out=root.parent/name
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as archive:
    for f in sorted(root.rglob('*')):
        if f.is_file() and '__pycache__' not in f.parts:
            archive.write(f,arcname=str(Path(root.name)/f.relative_to(root)))
with zipfile.ZipFile(out) as archive:
    assert archive.testzip() is None
print('Q550 v18.51z sendback:',out)
