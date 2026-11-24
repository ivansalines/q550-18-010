#!/usr/bin/env python3
"""Synthetic, independent tests of the variational equation and Cartesian chart."""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import chart_jacobian,project_components
# Exact toy system q'=v, v'=kappa: S_q=t^2/2, S_v=t.
def rhs(t,z):return np.array([z[1],1.])
sol=solve_ivp(rhs,(0.,.2),np.zeros(2),method='DOP853',rtol=1e-12,atol=1e-14,dense_output=True)
for t in (0.,.02,.1,.2):
 assert np.linalg.norm(sol.sol(t)-np.array([t*t/2,t]))<1e-11
# Analytic Cartesian immersion derivative checks independent of Q550 solver.
y=np.array([.2,-.4,1.2,.3,-.3,.1,-.8,.2,.07,.1],float)
mx=2.4
def Q(y):
 x=1/(1+np.exp(-y[8]));r1=np.exp(y[0]);r2=np.exp(y[4]);
 return np.array([r1*np.cos(y[2]),r1*np.sin(y[2]),r2*np.cos(y[6]),r2*np.sin(y[6]),np.sqrt(mx)*x])
q=Q(y);J=chart_jacobian(y,q)
for i in range(10):
 d=np.zeros(10);d[i]=1e-6
 numerical=(Q(y+d)-Q(y-d))/(2e-6)
 assert np.linalg.norm(numerical-J[:,i])<3e-10
# Planar polar decomposition preserves Cartesian squared norm.
dq=np.array([.3,-.2,.15,-.28,.03]);r1,t1,r2,t2=project_components(dq,y)
assert abs(r1*r1+t1*t1-float(dq[:2]@dq[:2]))<1e-15
assert abs(r2*r2+t2*t2-float(dq[2:4]@dq[2:4]))<1e-15
# A nonzero angular difference does not imply a zero radial difference.
assert abs(r2)>1e-3 and abs(t2)>1e-3
print('SYNTHETIC MODEL-INDEPENDENT VARIATION / EMBEDDING TESTS PASS')
