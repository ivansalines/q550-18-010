#!/usr/bin/env python3
"""Fast adversarial event/kinematic tests; no Q550 nonlinear solve."""
from pathlib import Path
import importlib.util
import sys
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
spec=importlib.util.spec_from_file_location('v1852g_unit',p)
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)

def parabola(t):
    # Physical trajectory q=t-t², v=1-2t, a=-2 in one Cartesian direction.
    # D=.5*q², G=q*v, C=v²+q*a.
    q=t-t*t;v=1-2*t;a=-2.
    D=.5*q*q;G=q*v;P=v*v;I=-q*a;C=P-I
    return {'D':D,'G':G,'C':C,'P':P,'I':I,'R':I/P if P!=0 else None,
      'xi_base':t,'xi_kick':t,'delta_xi':0.}

t=np.linspace(0,1.2,1002)
r=[parabola(x) for x in t]
ev=mod.find_events(t,r,parabola)
assert len(ev['convexity_zero'])==2,ev
assert abs(ev['convexity_zero'][0]['tau']-(3-np.sqrt(3))/6)<1e-9
assert abs(ev['convexity_zero'][1]['tau']-(3+np.sqrt(3))/6)<1e-9
assert any(abs(x['tau']-.5)<1e-9 for x in ev['rate_zero'])
assert any(abs(x['tau']-1.)<1e-9 for x in ev['rate_zero'])
assert len(ev['position_separation_maxima'])==1,ev
assert abs(ev['position_separation_maxima'][0]['tau']-.5)<1e-9
assert abs(parabola(.5)['G'])<1e-14
assert parabola(.21)['G']>0 and parabola(.75)['G']<0
# Distinct thresholds: at tau=.3 C<0 but D'=G>0.
assert parabola(.3)['C']<0 and parabola(.3)['G']>0
# An inward acceleration may coexist with strictly positive convexity.
assert parabola(.1)['I']>0 and parabola(.1)['C']>0
# Pure translation has G>0, C>0, and I=0 after coincidence.
def translation(t):
    return {'D':.5*t*t, 'G':t, 'C':1., 'P':1., 'I':0., 'R':0.,
            'xi_base':t,'xi_kick':t,'delta_xi':0.}
no=mod.find_events(t,[translation(x) for x in t],translation)
assert not no['convexity_zero'] and not no['rate_zero'] and not no['position_separation_maxima']
print('PASS: synthetic exact convexity/rate/maxima order and absence-of-return tests')
