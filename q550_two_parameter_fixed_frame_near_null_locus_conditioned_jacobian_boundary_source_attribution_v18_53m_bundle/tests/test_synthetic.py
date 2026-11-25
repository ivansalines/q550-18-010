#!/usr/bin/env python3
"""Pure synthetic algebraic sanity checks. Never substitute for a Q550 IVP."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from near_null_geometry import (projected,residual,boundary_jacobian,determinant_from_edges,
    optimize_projected_root,primitive_slope,primitive_at,fixed_grid,spaced_seeds)

xi=np.arange(0.,61.00001,.025)
anchor_o,anchor_w=14.4,12.8
f1=np.interp(xi,[0,12,18,48,50,56,60],[0,-2.4,3.6,0,0,0,0])
f2=np.interp(xi,[0,12,18,48,50,56,60],[0,0,0,0,-2.8,3.2,0])
C=np.zeros((len(xi),10));C[:,0]=f1;C[:,1]=f2
C[:,4]=np.interp(xi,[0,14.4,52.8,60],[0,3,4,5])  # real complement not obliged to vanish
basis=np.eye(10)[:,:2]
p=np.array([anchor_o,anchor_w]);f=projected(xi,C,basis,p)
assert np.linalg.norm(f)<1e-12
assert np.linalg.norm(residual(xi,C,*p))>0.9
J=boundary_jacobian(xi,C,basis,p)
assert np.linalg.norm(J-np.array([[-1.,0.],[1.,3.]]))<1e-9,J
source_start=basis.T@primitive_slope(xi,C,p[0]);source_end=basis.T@primitive_slope(xi,C,p[0]+3*p[1])
assert abs(np.linalg.det(J)-determinant_from_edges(source_start,source_end))<1e-10
for h in [1e-4,1e-5]:
    # Here the root hits interpolation knots: one-sided slopes can change.
    for col in [0,1]:
        q=p.copy();q[col]+=h
        pp=p.copy();pp[col]-=h
        fd=(projected(xi,C,basis,q)-projected(xi,C,basis,pp))/(2*h)
        assert np.linalg.norm(J[:,col]-fd)<1e-8,(col,fd,J)
sol=optimize_projected_root(xi,C,basis,(14.1,12.81),(12.6,15.4),(12.70,12.91))
assert sol['projected_norm']<1e-8,sol
assert abs(sol['origin_xi']-anchor_o)<1e-7
assert abs(sol['width_xi']-anchor_w)<1e-7
assert sol['full_bivector_norm']>0.9
assert sol['sigma_min']>0.1
assert sol['boundary_det_identity_error']<1e-9
origin=np.linspace(12.8,15.2,25);width=np.linspace(12.74,12.88,29)
grid=fixed_grid(xi,C,basis,origin,width)
assert len(grid)==725
assert min(x[2] for x in grid)<1e-9
assert spaced_seeds(grid,12,(12.6,15.4),(12.70,12.91))
rot=np.array([[np.cos(.637),-np.sin(.637)],[np.sin(.637),np.cos(.637)]])
for testpoint in [(14.1,12.79),(14.4,12.80),(15.1,12.84)]:
    a=projected(xi,C,basis,testpoint)
    b=projected(xi,C,basis@rot,testpoint)
    assert abs(np.linalg.norm(a)-np.linalg.norm(b))<1e-12
    q=boundary_jacobian(xi,C,basis,testpoint)
    r=boundary_jacobian(xi,C,basis@rot,testpoint)
    assert np.linalg.norm(np.sort(np.linalg.svd(q,compute_uv=False))-np.sort(np.linalg.svd(r,compute_uv=False)))<1e-12
print('SYNTHETIC_2D_ROOT_BOUNDARY_JACOBIAN_SO2_GAUGE_AND_FULL_NET_NONCLOSURE_PASS')
