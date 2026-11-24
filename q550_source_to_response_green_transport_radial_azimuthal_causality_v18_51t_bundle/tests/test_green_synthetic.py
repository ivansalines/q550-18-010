#!/usr/bin/env python3
"""Independently soluble triangular system tests source->Green transport."""
import importlib.util,sys
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts'/'run_audit.py'
s=importlib.util.spec_from_file_location('q550_v1851t_synthetic_green',p)
m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
xi=np.linspace(0.,.2,513)
J=np.zeros((10,10));J[4,5]=1.3;J[6,7]=.9;J[6,5]=.4;J[4,7]=-.5;J[4,9]=.12
b=np.zeros(10);b[5]=10;b[7]=2;b[9]=.03
assert np.allclose(J@J,np.zeros((10,10)))
idx=np.arange(0,513,2)
a=m.transport(xi,idx,np.repeat(J[None],len(idx),axis=0),np.repeat(b[None],len(idx),axis=0))
sc=a['sc'];g=a['green']
exact=xi[:,None]*b[None]+.5*xi[:,None]**2*(J@b)[None]
assert np.max(np.abs(sc[:,:,3]-exact))<1e-9
assert np.max(np.abs(g[2::2,:,3]-exact[2::2]))<1e-9
assert np.max(np.abs(sc[:,:,:3].sum(axis=2)-sc[:,:,3]))<1e-10
assert np.max(np.abs(sc[:,:,2]))>0
assert abs(sc[-1,6,0])>0 and abs(sc[-1,4,1])>0, 'Cross channels must survive'
assert np.max(np.abs(sc[2::2]-g[2::2]))<1e-9
print('Q550 v18.51t synthetic Green causal-channel algebra 7/7 PASS')
