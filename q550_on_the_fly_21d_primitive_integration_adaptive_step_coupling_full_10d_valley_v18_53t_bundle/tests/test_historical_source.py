#!/usr/bin/env python3
"""Test authentic force construction on synthetic states; no physical IVP."""
from pathlib import Path
from contextlib import ExitStack
import sys
import numpy as np
from scipy.special import expit
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_audit import load_ancestors
from physics_21d import source_10d,source_array
with ExitStack() as stack:
    s,r,h,m,k,o,j,centers,sroot,rr,meta=load_ancestors(stack)
    t=np.array([0.,0.21,.43]);a=np.zeros((len(t),11))
    a[:,0]=np.log(1.2);a[:,4]=np.log(1.1)
    a[:,2]=np.array([0.0,.02,.04]);a[:,6]=.35
    a[:,3]=.6;a[:,7]=-.45;a[:,8]=0.03;a[:,9]=-.1
    a[:,10]=np.array([0.,.35,.72])
    ctx={'event':0.,'theta':lambda x:.07*x}
    kp=.1666666666666; omega=.85
    historical=h.primary_metrics(t,a,ctx,kp)
    historical_source=np.array([m.to10(matrix) for matrix in historical['B_U']+historical['B_C']])
    historical_source *= (omega/expit(a[:,8]))[:,None]
    derived=source_array(h,m,ctx,kp,t,a,omega)
    err=np.max(np.abs(derived-historical_source))
    assert err<2e-13,err
    for idx in range(len(t)):
        item=source_10d(h,m,ctx,kp,float(t[idx]),a[idx],omega/expit(a[idx,8]))
        assert np.linalg.norm(item-derived[idx])<2e-13
print('AUTHENTIC_HISTORICAL_SOURCE_PASS original radial_potential+phase_coupling and xi-clock conversion')
