#!/usr/bin/env python3
"""Synthetic controls for fixed-hinge dual observer (no source required)."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import run_audit as audit
import numpy as np
x=audit.null_self_test()
assert all(v=='PASS' for v in x.values()),x
u=np.array([1.,0.,0.])
bplus=audit.observer_basis(+17.04,u)[1]
bminus=audit.observer_basis(-17.04,u)[1]
assert np.linalg.norm(bplus-bminus)>0
assert abs(bplus[0])<1e-14 and abs(bminus[0])<1e-14
assert abs(bplus[1]-bminus[1])<1e-14
assert abs(bplus[2]+bminus[2])<1e-14
print('SYNTHETIC HINGE / TWO SHADOWS TEST PASS')
