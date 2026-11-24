#!/usr/bin/env python3
"""Lightweight exact event and physical geometry control, no Q550 integration."""
from pathlib import Path
import importlib.util,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def import_m(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m
m=import_m(ROOT/'scripts/run_audit.py','v1852i_synthetic')
h=import_m(ROOT/'reference/parent_v18_52h_run_audit.py','v1852h_events')
# Genuine negative-to-positive zero of G at tau=2.0: D minimum, not orbit closure.
def evaluate(t):
    g=t-2.0
    return dict(tau=t,G=g,C=1.,D=1.+.5*g*g,P=2.,I=1.,R=.5,
                xi_base=t,xi_kick=t+1e-4,delta_xi=1e-4)
times=np.linspace(1.,3.,43);vals=[evaluate(t) for t in times]
roots=h.crossing_events(times,vals,evaluate,'G','rate_zero')
assert len(roots)==1 and roots[0]['crossing']=='negative_to_positive'
assert abs(roots[0]['tau']-2)<1e-10
assert m.event_discrimination(dict(rate_zero=roots,convexity_zero=[])).startswith('NUMERICAL_LOCAL_D_MINIMUM')
assert not m.event_discrimination(dict(rate_zero=[],convexity_zero=[])).startswith('NUMERICAL_LOCAL_D_MINIMUM')
# Exact rotating q and v in a plane embedded in R^5: G=0 but transverse speed nonzero.
a0=np.array([1,0,0,0,0],float)
for theta in (0.,np.pi/4,np.pi/2,np.pi):
    q=np.array([np.cos(theta),np.sin(theta),0,0,0])
    v=np.array([-np.sin(theta),np.cos(theta),0,0,0])
    p=dict(dQ=q,dV=v,G=float(q@v))
    geom=m.vector_geometry(p,a0,1e-12)
    assert abs(geom['radial_rate'])<1e-12
    assert abs(geom['tangential_speed']-1)<1e-12
    assert abs(geom['angular_speed']-1)<1e-12
    assert abs(geom['angle_from_h_seam_deg']-np.degrees(theta))<1e-8
assert not m.vector_geometry(dict(dQ=np.zeros(5),dV=np.ones(5)),a0,1e-12)['orientation_defined']
print('SYNTHETIC EVENT AND 5D POSTURE TEST PASS')
