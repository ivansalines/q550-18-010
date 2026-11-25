#!/usr/bin/env python3
"""Small deterministic 5D skew/SVD, phase and FTC falsification tests."""
import sys
from pathlib import Path
import numpy as np
r=Path(__file__).resolve().parents[1];sys.path.insert(0,str(r/'scripts'))
from fine_svd_geometry import to10,from10,metrics,transported_sign_frame,derivative_diagnostic,boundary_derivative
from three_window_geometry import accumulators,interval

def close(x,y,tol=5e-11):
    if not np.allclose(x,y,rtol=tol,atol=tol):raise AssertionError(f'mismatch {x} versus {y}')

# Isometric coordinate embedding and antisymmetry under 5D rotations.
rng=np.random.default_rng(551853)
B=rng.normal(size=(5,5));B=B-B.T
close(np.linalg.norm(to10(B)),np.linalg.norm(B));close(from10(to10(B)),B)
Q,_=np.linalg.qr(rng.normal(size=(5,5)));close(np.linalg.norm(to10(Q@B@Q.T)),np.linalg.norm(to10(B)))

# An equilateral planar triad and small third-coordinate drift can produce
# ~100% transverse residual share without a special dynamical mechanism.
E=np.zeros((3,10))
E[0,0]=1.;E[1,0]=-.5;E[1,1]=np.sqrt(3)/2;E[2,0]=-.5;E[2,1]=-np.sqrt(3)/2
for i in range(3):E[i,2]=.04
m=metrics(E)
assert m['direction_defined'] and m['transverse_share_squared']>.999999
close(m['cancellation_score'],1-np.linalg.norm(E.sum(0))/sum(np.linalg.norm(E,axis=1)))
# Arbitrary O(10) observer-coordinate change does not alter intrinsic SVD metrics.
W,_=np.linalg.qr(rng.normal(size=(10,10)));rot=metrics(E@W)
for k in ('net_norm','cancellation_score','transverse_norm','planar_norm','transverse_share_squared','rank23_gap'):
    close(m[k],rot[k])
# Relabeling three windows does not change the aggregate SVD quantities.
perm=metrics(E[[2,0,1],:]);close(perm['transverse_norm'],m['transverse_norm'])

# Gauge-continuation test: sign of an SVD axis is not itself an observable.
rows=[]
for i,w in enumerate([12.70,12.75,12.80,12.85]):
    D=E.copy();D[:,2]=w-12.775
    rows.append({'width':w,'origin':12.8,'seed':'synthetic','kn':'shared','m':metrics(D)})
transported_sign_frame(rows,12.8)
assert all(x['transport_ok'] for x in rows)
p=[x['signed_transverse_coefficient'] for x in rows]
assert p[1]*p[2]<0,('Expected a gauge-continuous synthetic crossing',p)
assert all(x['v3_projector_cosine_to_anchor']>.999 for x in rows)

# Constant skew torque density and affine material clock: FTC endpoint derivative
# can be checked without approximations of a nonlinear model.
xi=np.linspace(0,60,2401);tau=1.7*xi
C=from10(np.array([1.,0,0,0,0,0,0,0,0,0]));Bu=np.repeat(C[None,:,:],len(xi),axis=0);Bc=np.zeros_like(Bu)
acc=accumulators(xi,tau,Bu,Bc)
for o,w in [(12.8,12.8),(14.4,12.775),(15.2,13.6)]:
    v=sum(interval(acc,o,o+3*w))
    close(v,1.7*3*w*C,1e-9)
    pred=boundary_derivative(Bu,Bc,xi,np.full(len(xi),1.7),o,w)
    close(pred,3*1.7*to10(C),1e-9)
print('V1853K_SYNTHETIC_PASS: skew-isometry, O(10) invariance, SVD gauge, artificial bias, sign crossing, FTC')
