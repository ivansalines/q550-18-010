#!/usr/bin/env python3
"""Independent historical RHS check for the negative signed-drive source.
Only compares to reference support RHS; full authentication is on server.
"""
import math
import numpy as np
from run_audit import B,C
reduced_field=B.reduced_field
from return_map import lift8
from energy_math import sigmoid,dU,Params
from test_support_q550_oriented_coupled_bulge_cone_audit_v1 import support_rhs_float
class Audit:
 @staticmethod
 def sigmoid(z):return sigmoid(z)
 @staticmethod
 def dU(r,e):return dU(r,e)
 @staticmethod
 def dV_barrier(x,B):return B*(2*x-1)/(x*(1-x))
p=Params(**C['model_params']);omega=.844357915384206
rng=np.random.default_rng(1853);worst=0.
for _ in range(64):
 u=np.array([rng.uniform(.3,.8),rng.normal(0,.02),rng.uniform(-.6,-.4),
   rng.uniform(-2.35,-1.9),rng.normal(0,.03),rng.normal(0,.08),
   rng.normal(0,.1),rng.uniform(-math.pi,math.pi)])
 J=-1.441016445824766;y=lift8(u,J,.73);phase=(y[2]+y[6])-u[-1]
 source=support_rhs_float(np.r_[y,phase],kappa0=C['v1852l_fixed_model_kappa'],params=p,audit=Audit())
 F=source*(omega/(2*math.pi));F[10]=-omega
 fd=np.array([F[0],F[1],F[3],F[4],F[5],F[8],F[9],F[2]+F[6]-F[10]])
 worst=max(worst,float(np.max(abs(fd-reduced_field(u,J,C['v1852l_fixed_model_kappa'],-1,omega)))))
assert worst<2e-11,worst
print('HISTORICAL_SUPPORT_RHS_INDEPENDENT_64_STATES_PASS',worst)
