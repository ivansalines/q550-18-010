#!/usr/bin/env python3
"""Exact event classification, 5D orientation and no-false-closure controls."""
from pathlib import Path
import importlib.util,sys
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
s=importlib.util.spec_from_file_location('v1852l_synthetic',p)
m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
assert m.SEAM==12.8
assert m.event_discrimination(dict(rate_zero=[dict(crossing='negative_to_positive')],convexity_zero=[])).startswith('NUMERICAL_LOCAL_D_MINIMUM')
assert not m.event_discrimination(dict(rate_zero=[dict(crossing='positive_to_negative')],convexity_zero=[])).startswith('NUMERICAL_LOCAL_D_MINIMUM')
a=np.array([1.,0,0,0,0]);b=np.array([0.,1,0,0,0]);c=2*a
x=m.compare_posture_vectors(a,b)
assert abs(x['angle_between_vectors_deg']-90.)<1e-10
assert abs(x['relative_norm_ratio']-1.)<1e-10
assert abs(x['distance_in_5d']-np.sqrt(2))<1e-10
x=m.compare_posture_vectors(a,c)
assert abs(x['angle_between_vectors_deg'])<1e-10 and abs(x['relative_norm_ratio']-2)<1e-10
assert abs(x['distance_in_5d']-1)<1e-10
for theta in (0.,np.pi/4,np.pi/2,np.pi):
    q=np.array([np.cos(theta),np.sin(theta),0,0,0]);v=np.array([-np.sin(theta),np.cos(theta),0,0,0])
    z=m.vector_geometry(dict(dQ=q,dV=v),a,1e-12)
    assert abs(z['radial_rate'])<1e-12 and abs(z['angular_speed']-1.)<1e-12
    assert abs(z['angle_from_h_seam_deg']-np.degrees(theta))<1e-8
assert m.vector_geometry(dict(dQ=np.zeros(5),dV=np.ones(5)),a,1e-12)['orientation_defined'] is False
print('SYNTHETIC EVENT + PHYSICAL 5D POSTURE TEST PASS')
