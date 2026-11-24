#!/usr/bin/env python3
"""Catch accidental global drive-phase reset and lost clock continuity."""
from pathlib import Path
import importlib.util,sys
from types import SimpleNamespace
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
s=importlib.util.spec_from_file_location('v1852i_restart',p);m=importlib.util.module_from_spec(s)
sys.modules[s.name]=m;s.loader.exec_module(m)
def field(x,z):
    v=np.zeros(11);v[0]=x;v[10]=2.;return v
fake=SimpleNamespace(model_rhs=lambda *args:field)
z=np.zeros(11);z[0]=10.;z[10]=2.7
sol=m.resume(fake,{},None,None,None,None,z,3.2,1e-11,1e-12,.04)
assert abs(sol.y[0,-1]-(10+.5*(3.2**2-1.6**2)))<1e-9
assert abs(sol.y[10,-1]-(2.7+2*(3.2-1.6)))<1e-9
assert np.allclose(sol.sol(1.6),z,atol=1e-11)
print('MOCK GLOBAL-XI + 11D CLOCK SEAM PASS')
