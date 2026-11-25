#!/usr/bin/env python3
"""Fast synthetic algebra tests: no physical integrations or Q550 claims."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from quadrature_math import (cumulative_trapezoid,baseline_quadrature,
                             source_pchip_integral,decompose_nested_2x4x,max_norm)

x1=np.linspace(1.2,3.2,21)
x2=np.linspace(1.2,3.2,41)
x4=np.linspace(1.2,3.2,81)
def fn(x):return np.stack([x*x,3*x+2,np.sin(x)],axis=-1)
y1=fn(x1);y2=fn(x2);y4=fn(x4)
# Simpson exact for polynomial <= degree 3, with a constant g-component and a sine diagnostic.
q2=baseline_quadrature(x2,y2,2,'simpson')
q4=baseline_quadrature(x4,y4,4,'simpson')
exact=np.stack([(x1**3-x1[0]**3)/3,
                1.5*(x1*x1-x1[0]*x1[0])+2*(x1-x1[0])],axis=-1)
assert np.max(abs(q4[:,:2]-exact))<1e-13
assert np.max(abs(q2[:,:2]-exact))<1e-13
T4=baseline_quadrature(x4,y4,4,'trapezoid')
assert max_norm(T4-q4)>0
c=cumulative_trapezoid(x4,y4)
assert np.max(abs(T4-c[::4]))<5e-14
split=decompose_nested_2x4x(x2,y2+.001,x4,y4)
assert split['identity_max_error']<5e-14
assert np.max(abs(split['total']-split['quadrature_on_common_g4']-split['physical_sample_change']))<5e-14
# For exact linear g, source-integrated PCHIP matches analytical integral.
linear=np.stack([2*x4+1,3*x4-2],axis=1)
pint=source_pchip_integral(x4,linear)
truth=np.stack([(x4*x4+x4)-(x4[0]**2+x4[0]),
                1.5*(x4*x4-x4[0]**2)-2*(x4-x4[0])],axis=1)
assert np.max(abs(pint-truth))<1e-12
# Perturbed source, identical knots; decomposition MUST attribute exact T2 term.
split2=decompose_nested_2x4x(x2,y2+.003,x4,y4)
assert max_norm(split2['physical_sample_change'])>0
print('V1853S_SYNTHETIC_QUADRATURE_AND_TELESCOPING_TESTS_PASS')
