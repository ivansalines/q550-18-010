#!/usr/bin/env python3
"""Fix1 regression: original optimizer box would leave native xi=[12.8,56]."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from near_null_geometry import (validate_trace_coverage, optimize_projected_root,
    projected)
xi=np.linspace(12.8,56.0,1729)
try:
    validate_trace_coverage(xi,(12.6,15.4),(12.7,12.91))
except ValueError as ex:
    assert 'START_BELOW_TRACE' in str(ex)
else:
    raise AssertionError('original unsafe box was accepted')
valid=validate_trace_coverage(xi,(12.8,15.4),(12.7,12.91))
assert valid['box_start_min']==12.8 and abs(valid['box_end_max']-54.13)<1e-10
C=np.zeros((len(xi),10))
C[:,0]=np.sin((xi-12.8)/7.)
C[:,1]=np.cos((xi-12.8)/8.)
B=np.eye(10)[:,:2]
assert np.all(np.isfinite(projected(xi,C,B,(12.8,12.8))))
for seed in [(12.80000001,12.8),(15.399999,12.90999),(14.4,12.8)]:
    result=optimize_projected_root(xi,C,B,seed,(12.8,15.4),(12.7,12.91))
    assert np.isfinite(result['projected_norm'])
    assert result['origin_xi']>=xi[0]-1e-10
    assert result['end_xi']<=xi[-1]+1e-10
print('FIX1_TRACE_DOMAIN_GUARD_AND_BOUNDARY_SAFE_FD_PASS')
