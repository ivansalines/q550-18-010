#!/usr/bin/env python3
"""Synthetic event tests: known C-zero, later G-zero, and independent ratio crest."""
import importlib.util,sys
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
spec=importlib.util.spec_from_file_location('h52test',p)
h=importlib.util.module_from_spec(spec);sys.modules[spec.name]=h;spec.loader.exec_module(h)
T=np.linspace(.1,3.,600)
def physical(t):
    t=float(t)
    # D'=G=2t-t², D''=C=2-2t. D maximum at t=2; G crest at t=1.
    G=2*t-t*t;C=2-2*t;P=5.;I=P-C
    return dict(tau=t,D=t*t-t**3/3,G=G,C=C,P=P,I=I,
                R=I/P,xi_base=t,xi_kick=t,delta_xi=0)
values=[physical(t) for t in T]
c=h.crossing_events(T,values,physical,'C','convexity_zero')
g=h.crossing_events(T,values,physical,'G','rate_zero')
assert len(c)==1 and len(g)==1,(c,g)
assert abs(c[0]['tau']-1)<1e-9
assert abs(g[0]['tau']-2)<1e-9
assert c[0]['crossing']=='positive_to_negative'
assert g[0]['crossing']=='positive_to_negative'
assert c[0]['tau']<g[0]['tau']
assert len([e for e in g if e['crossing']=='positive_to_negative'])==1
def ratio_physical(t):
    p=physical(t);p['P']=4.;p['I']=4.*(1.7-.15*(t-1.9)**2)
    p['R']=p['I']/p['P'];p['C']=p['P']-p['I']
    return p
R=[ratio_physical(t) for t in T]
crests=h.ratio_crest_events(T,R,ratio_physical,1e-12)
assert len(crests)==1 and abs(crests[0]['tau']-1.9)<2e-4,crests
fine=h.event_convergence(c+g,(float(T[0]),float(T[-1])),physical,1e-12)
assert len(fine)==2 and all(x['passed'] for x in fine)
print('SYNTHETIC PASS: convexity first, rate arrest second, independent ratio crest')
