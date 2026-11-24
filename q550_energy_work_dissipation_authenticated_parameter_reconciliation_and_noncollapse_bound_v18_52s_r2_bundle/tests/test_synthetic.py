#!/usr/bin/env python3
"""Compact synthetic tests: fail on reversed phase-drive sign or missing work."""
import sys,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from energy_math import Params,energy_terms,a_priori_bound,exact_energy_derivative
import numpy as np
p=Params();k=1/6
assert p.phi0==0 and p.e_reserve==0.01 and p.mu_x==0.05 and p.m_x==1.0
assert p.gamma_x==0.3 and p.barrier_B==0.3 and p.gamma_rho==0.3
z=np.array([.5,.2,.6,-.1,-.1,.4,1.4,.25,0.,.3,.7])
t=energy_terms(z,k,p)
assert abs(float(t['E_dot'])-float(-t['D_rho']-t['D_x']+t['W_x']+t['W_Theta']))<1e-13
assert abs(float(t['H'])-float(t['E']-.5*t['x']*t['Jsum']))<1e-13
assert abs(float(t['H_dot'])-float(-t['D_rho']-t['D_x']+t['W_x']-.5*t['dx']*t['Jsum']))<1e-13
assert a_priori_bound(k,p)['K']>0
assert a_priori_bound(3,p)['status']=='NOT_APPLICABLE'
assert a_priori_bound(k,Params(e_reserve=0))['status']=='NOT_APPLICABLE'
print('SYNTHETIC_PASS')
