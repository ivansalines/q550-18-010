#!/usr/bin/env python3
"""Synthetic regression tests for gauge covariance and fixed-frame geometry."""
import importlib.util, pathlib, numpy as np
root=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('audit_x',root/'scripts/gauge_frame_audit.py')
x=importlib.util.module_from_spec(spec);spec.loader.exec_module(x)
rng=np.random.default_rng(53)
for lam in [0.,1.,4.]:
    for family in x.CFG['gauge_families']:
        R5=x.frame5(family,lam)
        assert np.linalg.norm(R5.T@R5-np.eye(5))<1e-12
        assert abs(np.linalg.det(R5)-1)<1e-12
        R10=x.mg.wedge_map(R5);R45=x.exterior_square(R10)
        assert np.linalg.norm(R10.T@R10-np.eye(10))<1e-12
        assert np.linalg.norm(R45.T@R45-np.eye(45))<1e-12
        source=np.cumsum(rng.normal(size=(31,10))*.02,axis=0)+np.linspace(0,1,31)[:,None]*rng.normal(size=(1,10))
        source/=np.linalg.norm(source,axis=1)[:,None]
        arc,P,_=x.mg.oriented_polygon(source,np.eye(10))
        arc_rot,Protat,_=x.mg.oriented_polygon(source@R10.T,np.eye(10))
        assert abs(arc-arc_rot)<1e-12
        assert np.linalg.norm(Protat-R45@P)<5e-12
        assert np.linalg.norm(R45.T@Protat-P)<5e-12
    if lam==0:
        for family in x.CFG['gauge_families']:
            assert np.linalg.norm(x.frame5(family,lam)-np.eye(5))<1e-14
N=np.column_stack([np.cos(np.linspace(0,1,141)),np.sin(np.linspace(0,1,141)),np.ones(141)*.1,np.zeros((141,7))])
d,r=x.moving_frame_test(N,np.eye(10),0.16)
assert d['covariant_arc_error']<1e-12
assert d['covariant_polygon_error']<1e-12
assert d['max_connection_pair_error']<1e-12
print('SYNTHETIC_5D_INDUCED_GAUGE_AND_COVARIANT_MOVING_FRAME_PASS')
