#!/usr/bin/env python3
"""Mock global-xi drive and material-time seam, no Q550 nonlinear calculation."""
from pathlib import Path
import importlib.util,sys
import numpy as np
from types import SimpleNamespace
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
s=importlib.util.spec_from_file_location('v1852j_mock',p);m=importlib.util.module_from_spec(s)
sys.modules[s.name]=m;s.loader.exec_module(m)
def field(x,z):
    f=np.zeros(11);f[0]=x;f[10]=2.;return f
fake=SimpleNamespace(model_rhs=lambda *args:field)
z=np.zeros(11);z[0]=10.;z[10]=5.4
sol=m.resume(fake,{},None,None,None,None,z,3.6,1e-11,1e-12,.04)
assert abs(sol.y[0,-1]-(10.+.5*(3.6**2-3.2**2)))<1e-9
assert abs(sol.y[10,-1]-(5.4+2.*(3.6-3.2)))<1e-9
assert np.allclose(sol.sol(3.2),z,atol=1e-11)
print('MOCK GLOBAL-XI + MATERIAL CLOCK SEAM PASS')
