from pathlib import Path
import hashlib,json,zipfile
root=Path(__file__).resolve().parents[1]
files=sorted(p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='OUTPUT_MANIFEST.json')
manifest={'version':'v18.52u','root':root.name,'file_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
files.append(root/'OUTPUT_MANIFEST.json')
path=root.parent/(root.name+'_server_sendback.zip')
with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:
 for p in files:z.write(p,p.relative_to(root.parent))
with zipfile.ZipFile(path) as z:
 assert z.testzip() is None
print('SENDBACK',path,path.stat().st_size)
