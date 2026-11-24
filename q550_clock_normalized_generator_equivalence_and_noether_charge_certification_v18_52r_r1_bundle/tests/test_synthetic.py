#!/usr/bin/env python3
"""Tests of clock normalization and charge rejection, using synthetic states."""
import importlib.util,sys
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts'/'clock_audit.py'
spec=importlib.util.spec_from_file_location('q550_clock_audit_synthetic',p)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
rng=np.random.default_rng(1852)
for _ in range(40):
 y=m.sample_state(rng)
 a=rng.normal(size=10);scale=0.134383735972165
 c,res=m.align(scale*a,a)
 assert abs(c-scale)<1e-14 and res<1e-14
 broken=scale*a.copy();broken[3]+=.2
 _,res2=m.align(broken,a)
 assert res2>1e-4
 d=m.dcharge(y,broken)
 assert np.isfinite(d[2])
 assert abs(m.charge(m.anti_rotate(y,.8))-m.charge(y))<1e-12
assert m.archived_self_test(80)['pass']
print('SYNTHETIC PASS: constant clock, mismatch rejected, charge and archived ODE')
