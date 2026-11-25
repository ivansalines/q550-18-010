#!/usr/bin/env python3
"""Synthetic algebra, gauge, and boundary tests only. Not Q550 physical evidence."""
import numpy as np,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from oriented_valley_geometry import net_geometry,wedge,cos,derive_curves,fd_derivative_check
from scipy.interpolate import CubicSpline
x=np.linspace(0,20,4001)
np.random.seed(123)
raw=np.column_stack([np.sin(x*(.15+i*.014)+i*.3)+.002*i*x for i in range(10)])
c=CubicSpline(x,raw,axis=0,extrapolate=False)
at=lambda t:np.asarray(c(t))
slope=lambda t:np.asarray(c(t,1))
o,w=2.301,2.0003
q=net_geometry(o,w,at,slope)
h=1e-6
np.testing.assert_allclose(q['jo'],((at(o+h+3*w)-at(o+h))-(at(o-h+3*w)-at(o-h)))/(2*h),atol=3e-9)
np.testing.assert_allclose(q['jw'],((at(o+3*(w+h))-at(o))-(at(o+3*(w-h))-at(o)))/(2*h),atol=3e-9)
np.testing.assert_allclose(q['du'], (np.eye(10)-np.outer(q['u'],q['u']))@q['jo']/q['r'],rtol=1e-12,atol=1e-12)
assert abs(np.linalg.norm(q['form'])-q['ang_speed'])<1e-10
assert abs(q['radial']-q['start_radial']-q['end_radial'])<1e-10
Q,_=np.linalg.qr(np.random.randn(10,10))
assert abs(np.linalg.norm(Q@q['u'])-1)<1e-10
np.testing.assert_allclose(np.linalg.norm(wedge(Q@q['u'],Q@q['du'])),np.linalg.norm(q['form']),atol=1e-10)
np.testing.assert_allclose(float((Q@q['u'])@(Q@q['jo'])),q['radial'],atol=1e-10)
# Stationary 10D chord can be nonzero and perpendicular to both endpoint sources.
n=np.zeros(10);n[0]=.07868;g0=np.zeros(10);g1=np.zeros(10);g0[1]=2;g1[2]=-3
assert np.dot(n,g0)==0 and np.dot(n,g1)==0 and np.linalg.norm(n)>0
# Exercise the entire curve builder on smooth synthetic data in a wider domain.
curve=derive_curves(x,at,slope,6.5,2.,halfwidth=.15,npoints=31,width_range=(1.8,2.2))
assert len(curve['fixed_rows'])==31 and len(curve['valley_rows'])==31
assert np.all(np.isfinite(curve['fixed_unit']))
fd=fd_derivative_check(o,w,at,slope,h,x)
assert fd['evaluated'] and fd['relative_defect']<1e-7
print('SYNTHETIC_ORIENTATION_AND_BOUNDARY_TESTS_PASS')
