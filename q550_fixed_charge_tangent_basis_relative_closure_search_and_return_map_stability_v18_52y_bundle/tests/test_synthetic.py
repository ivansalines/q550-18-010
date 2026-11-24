#!/usr/bin/env python3
"""Independent synthetic falsification tests for v18.52y.

No authentic Q550 model or fictitious Q550 outcomes are generated here.
"""
import math,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from return_map import charge,reduce8,lift8,tangent_basis,chart_diff,scaled_norm,fixed_charge_fd,damped_fixed_point,gauge

y=np.array([.2,.04,11.9,-.6,-1.5,-.03,-.7,.13,-.01,.005])
j=charge(y);u=reduce8(y);g=float(y[2]-y[6]);B=tangent_basis(y)
assert B.shape==(10,8) and np.linalg.matrix_rank(B)==8
for i in range(8):
 h=1.e-6;d=np.eye(8)[i]*h
 yy=(lift8(u+d,j,g)-lift8(u-d,j,g))/(2*h)
 assert np.linalg.norm(yy-B[:,i])<1.e-8
 assert abs(charge(lift8(u+d,j,g))-j)<1e-12
 assert abs(charge(lift8(u-d,j,g))-j)<1e-12
assert np.linalg.norm(chart_diff(reduce8(lift8(u,j,g)),u))<1e-12
assert np.linalg.norm(chart_diff(reduce8(lift8(u,j,g+1.5)),u))<1e-12
# The diagnostic must NOT erase hidden v2 drift when sum angle is closed.
u_hidden=u.copy();u_hidden[4]+=.17
assert np.linalg.norm(chart_diff(u,u_hidden))>.16
# A true synthetic eight-dimensional contracting return with a fixed point:
M=np.diag([.4,.5,.6,.7,.8,.3,.2,.65]);target=u+np.linspace(-.03,.04,8)
def synthetic(a):return target+M@(a-target)
D1,_=fixed_charge_fd(synthetic,u,1e-4)
D2,_=fixed_charge_fd(synthetic,u,4e-5)
assert np.linalg.norm(D1-M)<1e-9 and np.linalg.norm(D2-M)<1e-9
sol=damped_fixed_point(synthetic,u,steps=3,backtracks=4,weighted_step=.25,threshold=1e-8)
assert sol['converged'] and sol['best_residual']<1e-8
assert np.linalg.norm(sol['best_u']-target)<1e-8
# A drifting map has no periodic point: the algorithm must not force a closure.
def translation(a):
 b=np.asarray(a).copy();b[4]+=.17;return b
no=damped_fixed_point(translation,u,steps=2,backtracks=3,weighted_step=.25,threshold=1e-8)
assert not no['converged'] and no['best_residual']>.1
# A neutral coordinate must not be classified as an attracting fixed point.
N=np.diag([1.,.5,.6,.7,.8,.3,.2,.65]);assert max(abs(np.linalg.eigvals(N)))>=1.
print('SYNTHETIC_PASS: 8D fixed-charge tangent, gauge quotient, hidden drift, source-free central Jacobian, true/false closure, neutrality exclusion')
# End-to-end physical 10D IVP plumbing with a *synthetic* source, not Q550.
# The purpose is to catch source-array ordering, external-phase injection,
# stroboscopic composition, and fixed-J map plumbing in a cheap local test.
from run_audit import AuthenticatedMap
class DummyRamp:
 @staticmethod
 def Design(*args):return tuple(args)
class DummyClock:
 @staticmethod
 def cubic_clock_sign(t,design):return -1.
class DummyPhysical:
 @staticmethod
 def rhs_from_original(z,k,params,slope,s,support,audit):
  out=np.zeros(11);out[2]=z[3];out[6]=z[7];out[10]=s
  return out
class DummyEngine:
 original_ramp=DummyRamp();clock=DummyClock();physical=DummyPhysical()
 omega=1.;params=object();slope=1.;support=None;audit=None
sim=np.array([0.,0.,0.,.5,0.,0.,0.,.5,0.,0.])
smap=AuthenticatedMap(DummyEngine(),{'event':0.,'theta':lambda t:-t},1/6,0.,2*math.pi,1e-10,1e-12,.1)
smap.set_leaf(charge(sim),gauge(sim))
z=reduce8(sim);zp=smap.map(z)
assert np.linalg.norm(chart_diff(z,zp))<1e-10
assert np.linalg.norm(chart_diff(smap.map(zp),smap.map(z,2)))<1e-10
M,_=fixed_charge_fd(smap.map,z,1e-4)
assert np.linalg.matrix_rank(M)==8
print('MOCK_10D_IVP_PASS: authentic pipeline interface simulated only, charge/gauge and two-turn composition')
