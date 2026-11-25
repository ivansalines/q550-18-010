#!/usr/bin/env python3
"""Regression check for the former str/Path ZIP arcname failure."""
import hashlib,json,subprocess,sys,tempfile
from pathlib import Path
from zipfile import ZipFile
script=Path(__file__).resolve().parents[1]/'scripts/package_sendback.py'
with tempfile.TemporaryDirectory(prefix='v1853v_pack_') as tmp:
    base=Path(tmp)/'synthetic_bundle';(base/'scripts').mkdir(parents=True);(base/'results').mkdir()
    (base/'scripts/package_sendback.py').write_bytes(script.read_bytes())
    (base/'README.md').write_text('synthetic packaging test\n')
    (base/'SOURCE_MANIFEST.json').write_text(json.dumps({'files':{'README.md':hashlib.sha256((base/'README.md').read_bytes()).hexdigest()}}))
    (base/'results/report.json').write_text(json.dumps({'n_cases':204,'scientific_decision':'SYNTHETIC_TEST_ONLY'}))
    for i in range(7):(base/'results'/f'dummy_{i}.csv').write_text('a,b\n1,2\n')
    subprocess.run([sys.executable,str(base/'scripts/package_sendback.py')],check=True,capture_output=True,text=True)
    with ZipFile(base.parent/(base.name+'_server_sendback.zip')) as z:
        assert z.testzip() is None
        p=base.name+'/'
        for m in ['SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json']:
            for rel,digest in json.loads(z.read(p+m))['files'].items():
                assert hashlib.sha256(z.read(p+rel)).hexdigest()==digest
print('SYNTHETIC_SENDBACK_PACKAGING_PASS')
