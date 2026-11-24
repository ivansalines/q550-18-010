import zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root.parent/f'{root.name}_server_sendback.zip'
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
 for p in sorted(root.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:
   z.write(p,arcname=(root.name+'/'+p.relative_to(root).as_posix()))
with zipfile.ZipFile(out) as z:
 assert z.testzip() is None
print('Q550 v18.52f sendback',out)
