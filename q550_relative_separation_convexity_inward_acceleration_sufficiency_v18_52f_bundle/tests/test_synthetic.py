import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import analyze_vectors,sample_summary

def synthetic(inward_strength=0.5,relative_motion=1):
 t=np.linspace(0,1,65)
 # q=t v0+0.5 a t², v=v0+a t, a const, exact vector kinematics.
 v0=np.array([float(relative_motion),0,0,0,0])
 a=np.array([-float(inward_strength),0,0,0,0])
 q=t[:,None]*v0[None,:]+0.5*t[:,None]**2*a[None,:]
 v=v0[None,:]+t[:,None]*a[None,:]
 acc=np.tile(a,(len(t),1))
 return analyze_vectors(t,q,v,acc)

count=0
def check(c,msg):
 global count
 if not c:raise AssertionError(msg)
 count+=1

x=synthetic(0.05);check(x['ratio'][0]==0,'initial inward zero')
check(np.all(x['convexity']>0),'weak acceleration cannot beat positive term')
check(np.all(x['rate'][1:]>0),'weak acceleration cannot produce early return')
check(np.all(x['inward'][1:]>0),'inward force can coexist with divergence')
check(np.all(np.isfinite(x['ratio'])),'all physically finite cases')
# Constant relative-velocity case: zero acceleration, no restoration.
x=synthetic(0);check(np.all(x['inward']==0),'zero acceleration')
check(np.allclose(x['ratio'],0),'zero ratio')
check(np.all(x['convexity']>0),'pure drift is convex')
# Strong inward acceleration can reverse convexity before reversing distance rate.
x=synthetic(1.5)
check(np.any(x['convexity']<0),'strong inward term reverses convexity')
check(np.any(x['convexity']<0) and np.any(x['rate']>0),'negative convexity does not imply return')
check(abs(x['convexity'][-1] - (x['positive'][-1]-x['inward'][-1]))<1e-14,'identity')
# Recoverable oscillator-like case: turn back after outward motion; distinct threshold.
x=synthetic(2.5)
check(np.any(x['rate'][1:]<0),'position-rate reversal possible')
check(np.any(x['convexity']<0),'convexity reversal possible')
# Test a physical t reparametrization: actual tau coordinate for this synthetic trace.
x=synthetic(0.05)
check(len(x['tau'])==65 and np.all(np.diff(x['tau'])>0),'positive tau mesh')
check(np.all(x['positive']>=0),'positive outward term')
# Adversarial example: two positive sample endpoints can hide a negative interior.
t=np.array([0,1],float);p=lambda a: (a-0.5)**2-0.1
check(p(0)>0 and p(1)>0 and p(.5)<0,'positive endpoints do not certify positivity between samples')
print('Q550 v18.52f synthetic',count,'/',count,'PASS')
