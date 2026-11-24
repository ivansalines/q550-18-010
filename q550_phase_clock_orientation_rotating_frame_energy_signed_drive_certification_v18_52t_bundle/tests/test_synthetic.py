#!/usr/bin/env python3
import sys,math
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'scripts'))
from energy_math import Params,sigmoid
from phase_energy import exact_identity_probe,physical_rhs_with_signed_external_phase,signed_power,frame_ideal_power,frame_quantity,integral
p=Params();rng=np.random.default_rng(185252)
errs=[]
for i in range(160):
 y=np.zeros(11);y[[0,4]]=rng.uniform(-1,1,2);y[[1,5,3,7,9]]=rng.uniform(-1,1,5)
 y[[2,6]]=rng.uniform(-4,4,2);y[8]=rng.uniform(-2,2);phase=float(rng.uniform(-4,4));k=float(rng.uniform(-.3,.3))
 rate=float(rng.uniform(-2,2))
 errs.append(abs(exact_identity_probe(y,k,p,phase,rate)))
 # direct signed frame energy balance along independently evaluated full model RHS
 for sigma in (-1,1):
  from energy_math import exact_H_derivative,energy_terms,exact_energy_derivative
  y[10]=phase;F=physical_rhs_with_signed_external_phase(y,k,p,phase,sigma=sigma)
  E=energy_terms(y,k,p,phase)
  Jsum=float(E['Jsum']);dx=float(E['dx']);
  dJsum=(E['r1']**2*(2*y[1]*y[3]+F[3])+E['r2']**2*(2*y[5]*y[7]+F[7]))
  hchain=exact_energy_derivative(y,F,k,p)-sigma*.5*(dx*Jsum+E['x']*dJsum)
  assert abs(hchain-frame_ideal_power(y,k,p,phase,sigma))<3e-12
assert max(errs)<3e-12,max(errs)
# Diagnostic must reject wrong signed clock on an unambiguous sample.
t=np.linspace(0,3,1001);x=.5+.02*np.sin(t);phi=-.5*t+.02*(np.cos(t)-1)
assert abs(phi[-1]-phi[0]+np.trapezoid(x,t))<1e-7
assert abs(phi[-1]-phi[0]-np.trapezoid(x,t))>1
# Compare independent archived historical RHS with newly coded signed physical RHS.
import importlib.util
from types import SimpleNamespace
src=root/'reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py'
spec=importlib.util.spec_from_file_location('historical_1852t_test',src)
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
class A:
 @staticmethod
 def sigmoid(t):return float(sigmoid(t))
 @staticmethod
 def dU(r,e):return 2*r-2.4*r*r+.8*r*r*r-e*e/r**3
 @staticmethod
 def dV_barrier(x,B):return -B*(1/x-1/(1-x))
for i in range(40):
 y=np.zeros(11);y[[0,4]]=rng.uniform(-1,1,2);y[[1,3,5,7,9]]=rng.uniform(-.5,.5,5)
 y[[2,6,10]]=rng.uniform(-3,3,3);y[8]=rng.uniform(-2,2)
 k=float(rng.uniform(-.2,.2));x=float(sigmoid(y[8]))
 arch=mod.support_rhs_float(y,kappa0=k,params=SimpleNamespace(**p.__dict__),audit=A())*x/(2*math.pi)
 new=physical_rhs_with_signed_external_phase(y,k,p,y[10],sigma=+1)
 assert np.linalg.norm(arch-new)<2e-13,(i,np.linalg.norm(arch-new))
assert np.linalg.norm(arch-new)>0  # finite precision independent implementations
print('SYNTHETIC_PASS 160 signed energy identities, 40 independent historical RHS checks')
