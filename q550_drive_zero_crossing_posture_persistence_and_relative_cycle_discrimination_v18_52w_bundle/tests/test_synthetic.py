#!/usr/bin/env python3
"""Independent tests for posture metrics and gauge-quotient recurrence invariants."""
import sys,math
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import project_Q_and_tangent,unit_angle,unwrap_angle,full_one_turn_relative_map
from energy_math import sigmoid

y=np.array([math.log(.7),.2,.3,.4,math.log(.8),-.15,-.35,-.2,.1,.03]);f=np.array([.2,.01,.4,.05,-.15,.2,-.2,.02,.03,.01,0.])
q,u=project_Q_and_tangent(y,f)
h=1e-7
p,_=project_Q_and_tangent(y+h*f[:10],f);m,_=project_Q_and_tangent(y-h*f[:10],f)
assert np.linalg.norm((p-m)/(2*h)-u)<1e-8
assert unit_angle(np.array([1.,0.]),np.array([0.,1.]))==90.
assert abs(unwrap_angle(2*math.pi))<1e-13
# Synthetic relative periodic orbit: rotation angles have opposite slopes and invariant sum.
xi=np.linspace(0,4*math.pi,257);phase=-xi;Y=np.tile(np.r_[y,0.],(len(xi),1));Y[:,2]=.25*xi;Y[:,6]=-.25*xi
Y[:,10]=xi
rr=full_one_turn_relative_map(xi,{'original__reference':Y},{'original':phase})
assert len(rr)==17
assert max(abs(r['quotient_error_absolute_mixed_units']) for r in rr)<1e-12
assert max(abs(r['phase_match_error_rad']) for r in rr)<1e-12
# Counterexample: a projected angle can return while a hidden velocity drifts.
Y[:,3]+=.03*xi
bad=full_one_turn_relative_map(xi,{'original__reference':Y},{'original':phase})
assert min(r['quotient_error_absolute_mixed_units'] for r in bad)>.1
print('SYNTHETIC_PASS projection, gauge invariance, hidden-drift falsifier, one-turn alignment')

# End-to-end synthetic zero-drive event: the clock stops, but material tangent does not.
from run_audit import zero_crossing
class DummyClock:
 @staticmethod
 def cubic_clock_sign(t,design):return 1.-2.*t/.24
class DummyRamp:
 @staticmethod
 def Design(*args):return args
class DummyPhysical:
 def rhs_from_original(self,z,k,params,slope,s,source,audit):return np.r_[f[:10],.85*s]
class DummyEngine:
 clock=DummyClock();original_ramp=DummyRamp();physical=DummyPhysical()
 omega=.85;slope=.85/(2*math.pi);support=None;audit=None;params=None
xis=np.linspace(0,.2,513)
base=y[None,:]+xis[:,None]*f[None,:10]
ctx={'original':{'event':.091,'theta':lambda t:.85*(t-t*t/.24)},
     'autonomous':{'event':.09102,'theta':lambda t:.85*(t-t*t/.24)}}
events,evrep,curves=zero_crossing(xis,{k:base.copy() for k in ctx},1/6,DummyEngine(),ctx)
assert len(events)==6 and len(curves)==402
assert all(abs(q['theta_xi_rate_at_zero'])<1e-12 and q['physical_Q5_speed_xi_at_zero']>0 for q in evrep.values())
assert all(q['early_Y_derivative_vs_source_rhs_relative']<1e-8 for q in evrep.values())
assert all(r['sign_before']>0 and r['sign_after']<0 for r in events)
print('SYNTHETIC_DRIVE_ZERO_PASS nonzero material tangent, signed command, pre/post sampling')
