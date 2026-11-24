from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1]
files=sorted(x for x in (r/'results').iterdir() if x.is_file() and x.name!='LOCAL_PREFLIGHT.json')
m={'root':r.name,'version':'v18.51x','file_sha256':{str(x.relative_to(r)):hashlib.sha256(x.read_bytes()).hexdigest() for x in files}}
(r/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(m,indent=2)+'\n')
print('Q550 v18.51x output hashes written:',len(files))
