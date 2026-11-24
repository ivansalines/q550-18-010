#!/usr/bin/env python3
"""Eight-dimensional chart and genuinely fixed-angular-charge perturbations.

Physical state order in the authenticated Q550 10D RHS:
  [logr1, v1, theta1, omega1, logr2, v2, theta2, omega2, z, vz].
Only the simultaneous angular difference theta1-theta2 is quotiented.  J is
held exactly fixed by eliminating omega2, not adjusted after propagating.
Neither chart Euclidean norm nor singular values are kinetic metric invariants.
"""
import math
import numpy as np

TWOPI=2*math.pi
LABELS=('logr1','v1','omega1','logr2','v2','z','vz','alpha')

def wrap(a): return math.atan2(math.sin(float(a)),math.cos(float(a)))

def charge(y):
 y=np.asarray(y,float)
 return float(np.exp(2*y[0])*y[3]-np.exp(2*y[4])*y[7])

def reduce8(y):
 y=np.asarray(y,float)
 return np.array([y[0],y[1],y[3],y[4],y[5],y[8],y[9],wrap(y[2]+y[6])],float)

def gauge(y):
 y=np.asarray(y,float);return float(y[2]-y[6])

def lift8(u,J,g):
 u=np.asarray(u,float);y=np.zeros(10,float)
 y[0]=u[0];y[1]=u[1];y[3]=u[2];y[4]=u[3];y[5]=u[4];y[8]=u[5];y[9]=u[6]
 y[2]=(u[7]+g)/2;y[6]=(u[7]-g)/2
 y[7]=(math.exp(2*y[0])*y[3]-J)/math.exp(2*y[4])
 if not np.all(np.isfinite(y)):raise FloatingPointError('INVALID_FIXED_CHARGE_LIFT')
 return y

def tangent_basis(y):
 """10x8 analytic D(lift8) at fixed charge and fixed angular gauge."""
 y=np.asarray(y,float);B=np.zeros((10,8),float)
 for i,j in enumerate([0,1,3,4,5,8,9]):B[j,i]=1.
 B[2,7]=.5;B[6,7]=.5
 r=math.exp(2*(float(y[0])-float(y[4])))
 B[7,0]=2*r*y[3]
 B[7,2]=r
 B[7,3]=-2*y[7]
 return B

def chart_diff(a,b):
 d=np.asarray(b,float)-np.asarray(a,float);d=d.copy();d[7]=wrap(d[7]);return d

def closure_delta(u,p):return chart_diff(u,p)

def scaled_norm(v,scales):return float(np.linalg.norm(np.asarray(v,float)/np.asarray(scales,float)))

def fixed_charge_fd(fun,u,eps,scales=None):
 """Central finite differences of an EIGHT-dimensional map on a fixed J leaf.

 'fun' must already lift each perturbed chart point on the same exact leaf;
 it must return the eight coordinates of the evolved state.
 """
 u=np.asarray(u,float);scales=np.maximum(1.,np.abs(u)) if scales is None else np.asarray(scales,float)
 D=np.empty((8,8),float);columns=[]
 for col in range(8):
  h=float(eps)*float(scales[col]);du=np.zeros(8);du[col]=h
  up=fun(u+du);um=fun(u-du)
  D[:,col]=chart_diff(um,up)/(2*h)
  columns.append({'index':col,'name':LABELS[col],'step':h})
 return D,columns

def damped_fixed_point(fun,start,steps=3,backtracks=4,weighted_step=.25,threshold=1e-7,fd_eps=1e-4):
 """Bounded damped Newton search for P(u)-u=0 on a fixed-charge chart.

 A closure proof is NOT returned unless an actual independently recomputed
 P(u) residual is small.  On failure return the best evaluated candidate.
 """
 x=np.asarray(start,float).copy();history=[];best=None
 scales=np.maximum(1.,np.abs(x))
 for it in range(int(steps)+1):
  p=np.asarray(fun(x),float);delta=closure_delta(x,p)
  residual=scaled_norm(delta,scales)
  rec={'iteration':it,'residual_normalized':residual,'u':x.tolist(),'P_u':p.tolist()}
  history.append(rec)
  if best is None or residual<best[0]:best=(residual,x.copy())
  if residual<threshold:return {'converged':True,'best_residual':best[0],'best_u':best[1], 'history':history}
  if it==steps:break
  D,_=fixed_charge_fd(fun,x,fd_eps)
  A=D-np.eye(8)
  s=np.linalg.svd(A,compute_uv=False)
  rec['newton_operator_singular_values']=s.tolist()
  try: step=np.linalg.lstsq(A,-delta,rcond=1e-10)[0]
  except (np.linalg.LinAlgError,ValueError):break
  length=scaled_norm(step,scales)
  if length>weighted_step:step*=weighted_step/length
  accepted=False
  for back in range(backtracks+1):
   factor=2.**(-back);proposal=x+factor*step
   try:
    dp=closure_delta(proposal,fun(proposal));rr=scaled_norm(dp,scales)
   except (ValueError,RuntimeError,FloatingPointError,OverflowError):continue
   if rr<residual*(1-1e-4*factor):
    x=proposal;rec['accepted_backtrack']=back;rec['accepted_residual']=rr;accepted=True;break
  if not accepted:rec['stagnated']=True;break
 return {'converged':False,'best_residual':best[0], 'best_u':best[1], 'history':history}
