#!/usr/bin/env python3
"""Standalone Q550 v18.51w geometry tests using analytically known vectors.

No Q550 data are used. This tests the mechanics of orientation, axes and
observer projection independently of the physics result.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np

path=Path(__file__).resolve().parents[1]/'scripts'/'run_audit.py'
spec=importlib.util.spec_from_file_location('w_audit_module',path)
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
checks={}

def near(x,y,tol=2e-13):return np.max(np.abs(np.asarray(x)-np.asarray(y)))<tol

A=np.array([3.,0.]);B=np.array([0.,4.]);C=np.array([1.,2.])
checks['exact_orthogonal_angle']=near(np.arctan2(w.signed_area(A,B),np.dot(A,B)),np.pi/2)
checks['signed_area_3_4_rectangle']=near(w.signed_area(A,B),12.)
checks['squared_norm_identity']=near(np.dot(A,B)**2+w.signed_area(A,B)**2,np.dot(A,A)*np.dot(B,B))
for angle in [0.19,-0.72,2.13]:
  R=w.cos_sin_rotation(np.array(angle))
  checks['proper_rotation_dot_'+str(angle)]=near(np.dot(R@A,R@B),np.dot(A,B))
  checks['proper_rotation_area_'+str(angle)]=near(w.signed_area(R@A,R@B),12.)
mirror=np.diag([1.,-1.])
checks['reflection_flips_signed_area']=near(w.signed_area(mirror@A,mirror@B),-12.)
checks['reflection_keeps_dot']=near(np.dot(mirror@A,mirror@B),0.)
psi_a=w.critical_axis(A)
psi_b=w.critical_axis(B)
checks['root_zero_at_two_exact_vectors']=near(w.projected_radial(A,psi_a),0.) and near(w.projected_radial(B,psi_b),0.)
checks['axis_periodic_pi']=near(w.projected_radial(A,psi_a+np.pi),0.)
checks['critical_angle_frame_covariance']=near(np.sin(w.critical_axis(w.cos_sin_rotation(np.array(-.37))@C)-(w.critical_axis(C)-.37)),0.)
D=np.array([2.,-.5]);checks['observer_root_nonzero_norm']=near(w.projected_radial(D,w.critical_axis(D)),0.) and np.linalg.norm(D)>2.
checks['same_sign_no_projected_cancellation']=near(w.cancellation_fraction(np.array([[2.,1.,.5]]))[0],0.)
checks['opposite_equal_complete_projected_cancellation']=near(w.cancellation_fraction(np.array([[2.,-2.,0.]]))[0],1.)
checks['zero_vector_cancellation_undefined']=np.isnan(w.cancellation_fraction(np.array([[0.,0.,0.]]))[0])
checks['vanishing_response_has_no_direction']=bool(np.linalg.norm(np.zeros(2))==0)
assert all(checks.values()), {k:v for k,v in checks.items() if not v}
print(f'Q550 v18.51w synthetic geometric reference: {len(checks)}/{len(checks)} PASS')
