#!/usr/bin/env python3
"""Small synthetic algebraic tests only. No physical Q550 integration."""
from pathlib import Path
import sys, tempfile, shutil, json, zipfile, subprocess
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'scripts'))
from run_audit import refine
from source_valley_math import Primitive,optimized_width,curve,wedge
x=np.r_[np.linspace(12.8,25.6,31),np.linspace(25.6,56,55)[1:]]
for factor in (2,4):
    y=refine(x,factor)
    assert len(y)==(len(x)-1)*factor+1
    assert np.max(abs(y[::factor]-x))<1e-13
# A primitive whose source is analytically known; no Q550 physical data.
c=np.zeros((len(x),10));c[:,0]=np.sin(x/4);c[:,1]=np.cos(x/4);c[:,2]=x/400
s=np.zeros_like(c);s[:,0]=np.cos(x/4)/4;s[:,1]=-np.sin(x/4)/4;s[:,2]=1/400
p=Primitive(x,c,s,'hermite')
q=p.eval(17.3);dq=p.eval(17.3,1)
assert q.shape==(10,) and dq.shape==(10,)
minimum=optimized_width(p,14.3,12.75,12.85,25)
assert np.isfinite(minimum['norm'])
origin=np.linspace(14.25,14.35,9);w=np.full(len(origin),12.8)
result=curve(p,origin,w)
assert np.isfinite(result['direct']) and result['direct']>0
assert np.isfinite(result['chain']) and result['chain']>0
assert result['oriented_polygon'].shape==(45,)
assert np.linalg.norm(wedge(np.eye(10)[0],np.eye(10)[1]))==1
# A disposable packaging regression for the historical Path concatenation bug.
with tempfile.TemporaryDirectory() as td:
    td=Path(td);fake=td/'test_root';(fake/'scripts').mkdir(parents=True);(fake/'results').mkdir()
    shutil.copy2(root/'scripts/package_sendback.py',fake/'scripts/package_sendback.py')
    (fake/'SOURCE_MANIFEST.json').write_text(json.dumps({'files':{}}))
    (fake/'results/report.json').write_text(json.dumps({'scientific_decision':'SYNTHETIC_TEST_ONLY'}))
    (fake/'results/check.csv').write_text('a\n1\n')
    subprocess.check_call([sys.executable,str(fake/'scripts/package_sendback.py')],stdout=subprocess.DEVNULL)
    with zipfile.ZipFile(td/'test_root_server_sendback.zip') as z:
        assert z.testzip() is None
        assert 'test_root/results/check.csv' in z.namelist()
print('V1853R_SYNTHETIC_GEOMETRY_AND_PACKAGING_PASS')
