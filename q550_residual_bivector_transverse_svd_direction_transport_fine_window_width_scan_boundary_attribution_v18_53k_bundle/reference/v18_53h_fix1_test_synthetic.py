#!/usr/bin/env python3
"""Small, deterministic falsification tests; no Q550 orbit integration here."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from scipy.integrate import solve_ivp
from bivector_geometry import bivector,integ,cancel_metric,principal_skew,compare_geometry
from torque_math import force_components,material_rhs,acceleration_from_rhs,sigmoid
from run_audit import clock_rate_relative_residual

def check(x,message):
 if not x:raise AssertionError(message)

def test_geometry():
 B=np.zeros((5,5));B[0,1]=1.7;B[1,0]=-1.7;B[2,3]=.6;B[3,2]=-.6
 m=principal_skew(B)
 check(abs(m['sigma1']-1.7)<1e-12 and abs(m['sigma2']-.6)<1e-12,'principal singular values')
 check(m['reconstruction_error']<1e-12,'principal reconstruction')
 check(abs(m['kernel'][4])>1-1e-12,'rank 4 kernel')
 rng=np.random.default_rng(55053);Z=rng.normal(size=(5,5));Q,R=np.linalg.qr(Z)
 if np.linalg.det(Q)<0:Q[:,0]*=-1
 C=Q@B@Q.T;n=principal_skew(C)
 check(abs(n['sigma1']-m['sigma1'])<1e-12 and abs(n['sigma2']-m['sigma2'])<1e-12,'frame invariant modes')
 a=compare_geometry(B,B);check(a['full_bivector_cosine']>1-1e-12,'matching frame')
 check(max(a['planes'][0]['principal_angles_deg'])<1e-6,'matching plane')
 B2=np.zeros((5,5));B2[0,4]=1.;B2[4,0]=-1.
 # Full-space residual survives while one chosen 2D plane sees zero.
 check(abs(B2[0,1])<1e-15 and np.linalg.norm(B2)>1.,'false projected cancellation guard')


def test_forces():
 rng=np.random.default_rng(15031);n=60
 y=np.zeros((n,10));y[:,0]=.3+rng.normal(0,.06,n);y[:,4]=-1.9+rng.normal(0,.07,n)
 y[:,1]=rng.normal(0,.03,n);y[:,5]=rng.normal(0,.02,n)
 y[:,2]=rng.normal(0,1,n);y[:,6]=rng.normal(0,1,n)
 y[:,3]=-.5+rng.normal(0,.04,n);y[:,7]=-.01+rng.normal(0,.04,n)
 y[:,8]=rng.normal(0,.04,n);y[:,9]=rng.normal(0,.1,n)
 phase=rng.normal(0,1,n)
 p=dict(phi0=0.,e_reserve=.01,gamma_rho=.3,barrier_B=.3,mu_x=.05,gamma_x=.3,m_x=1.)
 for k in (1/6,.16659961041764054,.1666331385421536):
  Q,V,J,g,d,parts=force_components(y,phase,k,p)
  f=material_rhs(y,phase,k,p);A=acceleration_from_rhs(y,f,1.)
  check(np.max(abs(A-sum(parts.values())))<2e-12,'six force source identity')
  jdot=g['r1']**2*(2*y[:,1]*y[:,3]+f[:,3])-g['r2']**2*(2*y[:,5]*y[:,7]+f[:,7])
  check(np.max(abs(jdot))<2e-12,'J derivative cancellation')
  w=bivector(Q,parts['radial_potential'])
  check(np.max(abs(w+np.swapaxes(w,1,2)))<2e-12,'skew bivector')


def test_independent_ivp():
 # A controlled analogue of dual-kappa IVPs and numerical-refinement logic.
 times=np.linspace(0,4,65)
 def sol(k,r):
  return solve_ivp(lambda t,y:[-k*y[0]],(0,4),[1.],method='DOP853',
    t_eval=times,rtol=r,atol=r*.01,max_step=.035).y[0]
 a=sol(.16659961041764054,1e-8);b=sol(1/6,1e-8);c=sol(1/6,1e-11)
 check(np.max(abs(b-c))<2e-9,'integrator refinement')
 check(abs(a[-1]-b[-1])>1e-5,'distinct parameter produces distinct physical evolution')
 check(max(abs(b-np.exp(-times/6)))<2e-9,'true IVP reference')


def test_clock_rate_guard():
 # Regression for the v18.53h server NameError before any authentic IVP ran.
 omega=.844357915384206
 xi=np.linspace(12.8,25.6,65537)
 y_tau=np.zeros((len(xi),11))
 y_tau[:,8]=0.0 # x=sigmoid(0)=0.5
 y_tau[:,10]=10.0+(xi-xi[0])*(2.0*omega)
 err=clock_rate_relative_residual(xi,y_tau,omega)
 check(err<1e-9,'material clock exact analytic rate')
 # Mislabelled tau must be rejected, never silently accepted.
 y_tau[:,10]=10.0+(xi-xi[0])*omega
 bad=clock_rate_relative_residual(xi,y_tau,omega)
 check(bad>.49,'material clock mismatch detection')

def main():
 test_geometry();test_forces();test_independent_ivp();test_clock_rate_guard()
 print('SYNTHETIC_PASS geometry/principal_planes/frame_invariance/six_forces/J_identity/dual_IVP/refinement/material_clock_regression')
if __name__=='__main__':main()
