#!/usr/bin/env python3
"""Falsification tests: gauge-only return, hidden drift, charge mixing, rank identifiability."""
import math,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import J,reduce8,lift8,diff8,row_one,rank_audit,phase_gates,interpolate_state,modangle

base=np.array([.5,.04,-1.5,-.6,-2.1,-.03,.7,.13,-.01,.005])
j=J(base);u=reduce8(base)
y=lift8(u,j)
assert abs(J(y)-j)<1e-12
assert np.linalg.norm(diff8(u,reduce8(y)))<1e-13
# Gauge transformation cannot move physical quotient at fixed charge.
g=base.copy();g[2]+=.77;g[6]-=.77
assert np.linalg.norm(diff8(reduce8(base),reduce8(g)))<1e-13
assert abs(J(g)-j)<1e-13
same=row_one('original__reference','original',0,0.,2*math.pi,0.,base,g)
assert same['relative_8D_chart_norm']<1e-12
# Hidden omega1/v2 drift survives the quotient even with charge fixed.
bad=base.copy();bad[3]+=.04
bad[7]=(math.exp(2*bad[0])*bad[3]-j)/math.exp(2*bad[4])
assert abs(J(bad)-j)<1e-12
assert row_one('original__reference','original',0,0.,2*math.pi,0.,base,bad)['relative_8D_chart_norm']>.03
# Synthetic phase inversion defines one-turn gates in a 1.7-turn archive; NOT a twice-iterated map.
xi=np.linspace(12.8,25.6,513)
ph=-xi*11./12.8 # 11 rad total
pp=phase_gates(xi,ph,17)
assert len(pp)==17 and all(xb>xa for xa,xb,_ in pp)
assert ph[0]-ph[-1]<4*math.pi
Y=np.tile(np.r_[base,0.],(len(xi),1))
assert np.linalg.norm(interpolate_state(xi,Y,15.)-base)<1e-14
assert np.linalg.norm(interpolate_state(xi,Y,15.,'quadratic')-base)<1e-13
# Four perturbations on DIFFERENT J leaves; no transverse fixed-J Jacobian certificate.
records={'original__reference':{'start':base,'end':base}}
for amp in (1e-4,3e-4):
 for sign in (-1,1):
  n=f'original_{"positive" if sign>0 else "negative"}_{amp:.1e}__kick'
  yy=base.copy();yy[0]+=sign*amp*.1;yy[3]+=sign*amp*.7;yy[7]+=sign*amp*.2
  zz=yy.copy();zz[1]+=.01*sign*amp
  records[n]={'start':yy,'end':zz}
secants,rank=rank_audit('original',0,records)
assert len(secants)==2
assert rank['fixed_charge_8D_Jacobian_identifiable'] is False
assert rank['max_abs_charge_offset']>0
assert not rank['all_kicks_on_reference_charge']
assert rank['input_span_max_possible_rank']==4<rank['chart_dimension_after_fixed_J']
assert not any(s['is_fixed_charge_floquet_multiplier'] for s in secants)
# Honest fixed-J perturbations still give <=4 source vectors; rank 8 never inferred.
for rec in records.values():rec['start'][7]=(math.exp(2*rec['start'][0])*rec['start'][3]-j)/math.exp(2*rec['start'][4])
_,r_fixed=rank_audit('original',0,records)
assert r_fixed['all_kicks_on_reference_charge']
assert not r_fixed['full_floquet_spectrum_identifiable']
print('SYNTHETIC_PASS gauge quotient, J leaf, hidden posture, phase gates, one-turn, mixed-charge exclusion, rank ceiling')
