#!/usr/bin/env python3
"""Synthetic section and full-5D comparison regression; no Q550 integration."""
import importlib.util,sys
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
spec=importlib.util.spec_from_file_location('v1852l_sections_test',p)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
def sample(t):
    q=np.array([1.0, np.sin(t), .4*np.cos(t), 0., 0.])
    v=np.array([0., np.cos(t), -.4*np.sin(t), 0., 0.])
    return dict(dQ=q,dV=v,G=float(q@v),D=.5*float(q@q),xi_base=t,xi_kick=t)
tt=np.linspace(.15,4*np.pi+.15,401)
result=m.section_audit(tt,[sample(float(t)) for t in tt],sample)
assert result['summary']['n_horizontal_section_crossings']==4,result['summary']
assert result['summary']['n_opposite_orientation_pairs']>=1,result['summary']
assert all(abs(x['q1'])<1e-9 for x in result['crossings'])
assert all(x['full_5d_distance']>=0 for x in result['pairs'])
assert all(x['tangency_index']>0 for x in result['crossings'])
print('SYNTHETIC HORIZONTAL 5D SECTION TEST PASS')
