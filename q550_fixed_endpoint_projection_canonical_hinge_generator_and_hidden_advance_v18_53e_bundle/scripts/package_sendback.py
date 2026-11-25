"""Archive complete runnable genealogy and all scientific outputs with hashes."""
from pathlib import Path
import json,hashlib,zipfile
r=Path(__file__).resolve().parents[1];o=r/'results'
f={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(o.rglob('*')) if p.is_file()}
(r/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.53e','root':r.name,'files':f},indent=2)+'\n')
a=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(a,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(r.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:z.write(p,r.name+'/'+str(p.relative_to(r)))
print('SENDBACK',a,len(f),'output SHA256')
