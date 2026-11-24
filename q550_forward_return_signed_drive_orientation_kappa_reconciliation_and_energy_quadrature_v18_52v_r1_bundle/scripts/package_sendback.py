from pathlib import Path
import hashlib,json,zipfile
r=Path(__file__).resolve().parents[1]
files=sorted(p for p in r.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='OUTPUT_MANIFEST.json')
m={'version':'v18.52v-r1','root':r.name,'file_sha256':{str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(r/'OUTPUT_MANIFEST.json').write_text(json.dumps(m,indent=2)+'\n');files.append(r/'OUTPUT_MANIFEST.json')
p=r.parent/(r.name+'_server_sendback.zip')
with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:
 for f in files:z.write(f,f.relative_to(r.parent))
with zipfile.ZipFile(p) as z:assert z.testzip() is None
print('SENDBACK',p,p.stat().st_size)
