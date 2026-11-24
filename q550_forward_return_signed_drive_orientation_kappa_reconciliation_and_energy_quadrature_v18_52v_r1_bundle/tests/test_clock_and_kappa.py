#!/usr/bin/env python3
"""Independent synthetic falsification tests; no parent engine required."""
import sys, math, numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from energy_math import energy_terms,Params
from phase_energy import exact_identity_probe,integral
C=__import__('json').loads((ROOT/'inputs_config.json').read_text())
ko=C['kappa_seeds']['original'];ka=C['kappa_seeds']['autonomous']
kf=C['v1852l_fixed_model_kappa']
assert abs(kf-(ko+ka)/2)<5e-11
assert abs(kf-ko)>1e-5 and abs(kf-ka)>1e-5
omega=0.844357915384206
# Synthetic smooth ramp, independent of historical beta implementation:
# sign = 1-2*(3u²-2u³), excursion = -t+2*R*(u³-0.5u⁴)
R=.24
for t in np.linspace(.005,R-.005,64):
 u=t/R
 sign=1-2*(3*u*u-2*u**3)
 def exc(tt):
  z=tt/R
  return -tt+2*R*(z**3-0.5*z**4)
 num=(exc(t+1e-6)-exc(t-1e-6))/(2e-6)
 assert abs(num+sign)<1e-8
 theta_rate=-omega*num
 assert abs(theta_rate-omega*sign)<1e-8
 # Deliberately wrong diagnostic should reject the sign for internal ramp points.
 if abs(sign)>.2:assert abs(theta_rate+omega*sign)>.3*omega
# Fixed-kappa exact identity independently checked for arbitrary states and phase rates.
rng=np.random.default_rng(18525261);p=Params(**C['model_params'])
for i in range(40):
 y=np.zeros(11)
 y[[0,4]]=rng.uniform(-1,1,2)
 y[[1,3,5,7,9]]=rng.uniform(-.8,.8,5)
 y[[2,6,10]]=rng.uniform(-3,3,3)
 y[8]=rng.uniform(-1.3,1.3)
 rate=rng.uniform(-.75,.75)
 assert abs(exact_identity_probe(y,kf,p,y[10],rate))<1e-12
print('CLOCK_AND_FIXED_KAPPA_SYNTHETIC_PASS 40 + 64')
