#!/usr/bin/env python3
"""Check local reduced reconstruction against historical support_rhs_float on off-root states.

This test does not authenticate the full ancestor archive; strict server
verification later stages rhs_from_original and the external phase protocol.
"""
import math,sys
from pathlib import Path
import numpy as np
from run_audit import xi_rhs_local,K,P
from return_map import lift8
from energy_math import sigmoid,dU
from test_support_q550_oriented_coupled_bulge_cone_audit_v1 import support_rhs_float
class Audit:
 @staticmethod
 def sigmoid(z):return sigmoid(z)
 @staticmethod
 def dU(r,e):return dU(r,e)
 @staticmethod
 def dV_barrier(x,B):return B*(2*x-1)/(x*(1-x))

omega=0.844357915384206
rng=np.random.default_rng(1853)
maxerr=0.
for _ in range(64):
 u=np.array([rng.uniform(.3,.8),rng.normal(0,.02),rng.uniform(-.6,-.4),
             rng.uniform(-2.35,-1.9),rng.normal(0,.03),rng.normal(0,.08),
             rng.normal(0,.1),rng.uniform(-math.pi,math.pi)])
 J=-1.441016445824766
 y=lift8(u,J,.73)
 phase=(y[2]+y[6])-u[-1]
 source=support_rhs_float(np.r_[y,phase],kappa0=K,params=P,audit=Audit())
 # Signed 11D historical wrapper sets the final derivative to -omega.
 F=source*(omega/(2*math.pi)); F[10]=-omega
 f=np.array([F[0],F[1],F[3],F[4],F[5],F[8],F[9],F[2]+F[6]-F[10]])
 error=np.max(abs(f-xi_rhs_local(u,J)))
 maxerr=max(maxerr,error)
assert maxerr<2e-11, maxerr
print('ARCHIVED_RHS_INDEPENDENT_RANDOM_STATES_PASS',64,'max_error',maxerr)
