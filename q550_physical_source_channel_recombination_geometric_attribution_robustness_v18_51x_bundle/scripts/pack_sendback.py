from pathlib import Path
import zipfile
r=Path(__file__).resolve().parents[1]
output=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(output,'w',allowZip64=True) as z:
 for f in sorted(r.rglob('*')):
  if f.is_file() and '__pycache__' not in f.parts:
   method=zipfile.ZIP_STORED if f.suffix in ('.zip','.npz','.png') else zipfile.ZIP_DEFLATED
   z.write(f,arcname=str(Path(r.name)/f.relative_to(r)),compress_type=method)
print('Q550 v18.51x sendback:',output)
