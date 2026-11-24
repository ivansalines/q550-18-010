#!/usr/bin/env python3
"""Check global time is not reset and the integrated material clock is continuous."""
import importlib.util,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
spec=importlib.util.spec_from_file_location('h52restart',p)
h=importlib.util.module_from_spec(spec);sys.modules[spec.name]=h;spec.loader.exec_module(h)
def f(x,z):
    # This test detects accidental restart at xi=0, not the physical Q550 field.
    a=np.zeros(11);a[0]=x;a[10]=1.;return a
mock=SimpleNamespace(model_rhs=lambda *args:f)
z=np.zeros(11);z[0]=10.;z[10]=1.3
sol=h.resume_ivp(mock,{},None,None,None,None,z,1.6,1e-11,1e-12,.05)
assert abs(sol.y[0,-1]-(10.+.5*(1.6**2-.8**2)))<1e-10
assert abs(sol.y[10,-1]-2.1)<1e-11
assert np.allclose(sol.sol(.8),z,atol=1e-11)
print('MOCK RESTART PASS: xi global, 11D material clock continuous')
