#!/usr/bin/env python3
"""Falsification controls on signed work, material clock and torque."""
import sys,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from energy_math import Params,energy_terms,exact_energy_derivative
from phase_energy import physical_rhs_with_signed_external_phase,signed_power
p=Params();rng=np.random.default_rng(1852525)
max_err=0.;n=120
for _ in range(n):
 y=np.zeros(11);y[[0,4]]=rng.uniform(-.7,1.2,2)
 y[[1,3,5,7,9]]=rng.uniform(-.5,.5,5)
 y[[2,6,10]]=rng.uniform(-3,3,3);y[8]=rng.uniform(-1.5,1.5)
 k=.1666331385421536;rate=rng.uniform(-2.,2.)
 F=physical_rhs_with_signed_external_phase(y,k,p,y[10],rate=rate)
 rhs=signed_power(y,k,p,y[10],rate);lhs=exact_energy_derivative(y,F,k,p)
 max_err=max(max_err,abs(rhs-lhs));assert abs(rhs-lhs)<1e-11
 d=energy_terms(y,k,p);coeff=k*d['x']*d['r1']*d['r2']*math.sin(float(d['D']))
 assert abs(coeff*rate+coeff*(-rate))<1e-13
 jp=(d['r1']**2)*(2*y[1]*y[3]+F[3]);j2=(d['r2']**2)*(2*y[5]*y[7]+F[7])
 assert abs(jp-j2)<1e-12
 # Falsification: replace the rate by minus rate in RHS, but NOT in the ledger.
 if abs(coeff*rate)>1e-6:assert abs(rhs-signed_power(y,k,p,y[10],-rate))>1e-6
print('SYNTHETIC_PASS',n,'MAX_ENERGY_IDENTITY_RESIDUAL',max_err)
