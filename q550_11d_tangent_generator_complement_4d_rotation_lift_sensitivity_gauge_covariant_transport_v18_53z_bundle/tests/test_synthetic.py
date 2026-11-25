#!/usr/bin/env python3
"""Independent analytic fixtures for recovered forcing and projected tangent rotation."""
from pathlib import Path
import importlib.util,sys,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from tangent_source import biv,parametric_forces,reconstruct_phase,config_jacobian,physical_velocity_jacobian,tangent_lift,jacobian_11d,rhs_material_xi
from connection_math import wedge5,exterior_generator
from complement_frame import skew,cayley
spec=importlib.util.spec_from_file_location('torque_fixture',ROOT/'reference/torque_math_v18_53f_authenticated.py')
torq=importlib.util.module_from_spec(spec);spec.loader.exec_module(torq)
p={'phi0':0.,'e_reserve':.01,'gamma_rho':.3,'barrier_B':.3,'mu_x':.05,'gamma_x':.3,'m_x':1.}
k=.1666331385421536;omega=.84435791663
xi=np.linspace(12.8,13.8,41)
z=np.zeros((len(xi),11));z[:,0]=.05;z[:,1]=.01;z[:,2]=.13;z[:,3]=.05
z[:,4]=-.05;z[:,5]=-.02;z[:,6]=-.09;z[:,7]=-.04;z[:,8]=.01;z[:,9]=.06
z[:,10]=np.linspace(0,1,len(xi))
external=.03+.06*(xi-xi[0]);D=z[:,2]+z[:,6]-external
q,pot,fc,fs=parametric_forces(z,k,p,torq)
g=omega/(1/(1+np.exp(-z[:,8])))[:,None]*(biv(q,pot)+np.cos(D)[:,None]*biv(q,fc)+np.sin(D)[:,None]*biv(q,fs))
phase,diag=reconstruct_phase(xi,z,g,omega,k,p,torq)
assert max(np.abs(phase(xi)-external))<1e-10,diag
assert diag['max_source_10d_fit_error']<1e-12,diag
for state in z[[0,20,40]]:
    L1=tangent_lift(state,1.,'fixed_generalized_velocity')
    L2=tangent_lift(state,1.,'fixed_cartesian_velocity')
    assert np.linalg.norm(config_jacobian(state)@L1-np.eye(5))<1e-12
    assert np.linalg.norm(config_jacobian(state)@L2-np.eye(5))<1e-12
    assert np.linalg.norm(physical_velocity_jacobian(state)@L2)<1e-12
    A1=jacobian_11d(state,0.03,omega,k,p,torq,2e-5)
    A2=jacobian_11d(state,0.03,omega,k,p,torq,1e-5)
    assert np.linalg.norm(A1-A2)<1e-6
    assert np.linalg.norm(A1[:,10])==0
A=np.zeros((5,5));A[1,2]=.3;A[2,1]=-.3
assert np.linalg.norm(exterior_generator(A)+exterior_generator(A).T)<1e-12
assert np.linalg.norm(cayley(.01*A).T@cayley(.01*A)-np.eye(5))<1e-12
print('SYNTHETIC_11D_SOURCE_PHASE_LIFT_AND_EXTERIOR_TRANSPORT_PASS')
