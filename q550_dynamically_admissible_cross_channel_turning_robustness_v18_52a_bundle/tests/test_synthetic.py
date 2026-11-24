#!/usr/bin/env python3
"""Light analytic algebra, 2D orientation, and numerical turning-budget checks."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import angular_budget,wedge

xi=np.linspace(0.,.2,513)
A=np.column_stack([2.+xi,-.1-xi+2*xi**2]);B=np.column_stack([.2+xi,1.+.7*xi]);C=np.column_stack([np.zeros_like(xi)+.01,xi**2])
v=np.stack([A,B,C,A+B+C],axis=2)
bud=angular_budget(xi,v)
checks={}
checks['total_vector_closure']=np.allclose(v[:,:,:3].sum(axis=2),v[:,:,3],atol=1e-14)
checks['wedge_antisymmetric']=np.allclose(wedge(A,B),-wedge(B,A))
checks['wedge_formula']=np.allclose(wedge(A,B),A[:,0]*B[:,1]-A[:,1]*B[:,0])
checks['turning_budget']=bud['wedge_reconstruction_max_abs']<1e-12
checks['same_total_under_explicit_recombination']=np.allclose(.9*A+.1*B+.1*A+.9*B,A+B)
checks['mixing_changes_cross_wedge']=not np.allclose(wedge(.9*A+.1*B,.1*A+.9*B),wedge(A,B))
checks['angle_wedge_positive']=np.all(bud['oriented_area']>0)
checks['finite_rates']=np.isfinite(bud['rate_total']).all()
phi=.1*np.sin(2*np.pi*xi/.2)
c=np.cos(phi);s=np.sin(phi)
rot=lambda z:np.column_stack([c*z[:,0]+s*z[:,1],-s*z[:,0]+c*z[:,1]])
checks['passive_angle_invariance']=np.allclose(wedge(rot(A),rot(B)),wedge(A,B),atol=1e-13)
checks['passive_dot_invariance']=np.allclose(np.sum(rot(A)*rot(B),axis=1),np.sum(A*B,axis=1),atol=1e-13)
checks['passive_norm_invariance']=np.allclose(np.linalg.norm(rot(A),axis=1),np.linalg.norm(A,axis=1))
checks['sum_of_rates']=np.allclose(bud['rate_total'],bud['rate_self_ab']+bud['rate_mixed']+bud['rate_remaining'],atol=1e-12)
checks['qualified_nonzero']=np.count_nonzero(xi>=.01)==487
print('Synthetic checks:',sum(checks.values()),'/',len(checks))
for key,ok in checks.items():print(('PASS' if ok else 'FAIL'),key)
if not all(checks.values()):raise SystemExit(1)
