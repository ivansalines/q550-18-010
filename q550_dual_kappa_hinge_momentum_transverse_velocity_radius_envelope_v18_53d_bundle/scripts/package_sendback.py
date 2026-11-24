import zipfile,json,hashlib
from pathlib import Path
r=Path(__file__).resolve().parents[1]
f={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (r/'results').rglob('*') if p.is_file()}
(r/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.53d','root':r.name,'files':f},indent=2)+'\n')
p=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for q in r.rglob('*'):
  if q.is_file() and '__pycache__' not in q.parts:z.write(q,r.name+'/'+str(q.relative_to(r)))
print('SENDBACK',p,len(f),'output hashes')
