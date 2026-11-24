#!/usr/bin/env python3
"""Independent synthetic tests: transverse displacement, phase, and returns."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import find_roots_on_grid,lag_fit
# Along-curve shift plus perpendicular change are distinguished exactly.
w=np.array([2.,-.4,1.,0.,0.]);perp=np.array([.4,2.,0.,0.,0.])
ell=-.0013;delta=ell*w+.003*perp/np.linalg.norm(perp)
meas=float(delta@w/(w@w)); rem=delta-meas*w
assert abs(meas-ell)<1e-14 and abs(rem@w)<1e-14
# A bona fide synthetic x turnaround and an equal-x return are detected.
t=np.linspace(0,.2,513);x=0.5-(t-.1)**2
turn=find_roots_on_grid(t,np.gradient(x,t,edge_order=2))
assert any(abs(a-.1)<.0005 for a in turn)
x2=.5+(t-.05)*(t-.15)
returns=find_roots_on_grid(t,x2,float(x2[0]))
assert abs(returns[0])<1e-12 and abs(returns[-1]-.2)<.0005
# A monotone x has no turn and no spurious initial recurrence.
assert not find_roots_on_grid(t,np.ones_like(t))
# One global shift is recovered for a translated straight synthetic trajectory.
P=np.column_stack([t,3*t,t*0,t*0,t*0]); ds=.0003
H=P+ds*np.array([1.,3.,0.,0.,0.])[None,:]
fit=lag_fit(t,P,H)
assert abs(fit['optimal_constant_delay']-ds)<1e-8
assert fit['rms_after_shift']<1e-10
print('Q550 v18.51r independent geometry, projected returns and phase-lag tests PASS')
