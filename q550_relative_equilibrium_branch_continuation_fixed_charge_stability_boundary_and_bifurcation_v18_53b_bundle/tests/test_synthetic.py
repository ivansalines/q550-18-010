#!/usr/bin/env python3
"""Fast local invariance tests and deliberate falsification. No physics IVP."""
import math,json,sys
from pathlib import Path
import numpy as np
from run_audit import (C,ROOT,radial_residual,solve,reduced_field,correct_and_audit,
  native_field,find_candidates,sha,pseudo_arclength_drive,refine_algebraic_fold,make_row)
from return_map import lift8,reduce8,tangent_basis,charge
ref=ROOT/'reference/v18_53a_authenticated_report.json'
assert sha(ref)==C['v1853a_report_sha256']
R=json.loads(ref.read_text()); r=R['results']['original']; J=r['J']; k=r['kappa']; s=-1.; w=0.844357915384206
u=np.array(r['equilibrium_reduced_chart']);a=np.array(r['algebraic_branch_root']['core_root'])
assert max(abs(radial_residual(a,J,k,s)))<1e-9
assert max(abs(reduced_field(u,J,k,s,w)))<1e-9
qa=correct_and_audit(a,J,k,s,w)
assert qa['spectral_abscissa'] < 0,qa['spectral_abscissa']
assert abs(qa['spectral_abscissa'] - r['floquet_vs_exp_instantaneous_jacobian']['generator_spectral_abscissa'])<2e-7
# Perturbations preserve J exactly by the lift, and provide 8 independent chart columns.
y=lift8(u,J,0.);B=tangent_basis(y)
assert np.linalg.matrix_rank(B)==8
rng=np.random.default_rng(53)
for i in range(40):
 du=rng.normal(0,1e-4,8);test=lift8(u+du,J,0.)
 assert abs(charge(test)-J)<1e-13
# Test continuation in three directions (without source claims).
for dJ,dk,ds in [(-.05,0.,0.),(0.,.01,0.),(0.,0.,.04)]:
 ok,aa,err=solve(a,J+dJ,k+dk,s+ds)
 assert ok and err<1e-9,(dJ,dk,ds,err)
 qa=correct_and_audit(aa,J+dJ,k+dk,s+ds,w)
 assert qa['root_8d_max_abs']<1e-8
# A zero phase clock must have no one-phase-turn multiplier, regardless of stability.
ok,aa,err=solve(a,J,k,0.)
assert ok
assert correct_and_audit(aa,J,k,0.,w)['phase_turn_xi'] is None
# Falsification: synthetic eigenvalue crossing is a candidate, never a proof of Hopf.
dummy=lambda x:{'parameter':x,'spectral_abscissa':x,'smallest_equilibrium_singular':1.}
can=find_candidates('charge',[dummy(-.1),dummy(.1)])
assert len(can)==1 and can[0]['spectral_sign_change'] and 'NOT_BIFURCATION_THEOREM' in can[0]['classification']
# Independent local saddle-node test from two predicted signed-drive roots.
rows=[]
for ss in (.9,.925):
 ok,aa,err=solve(a,J,k,ss)
 assert ok
 qa=correct_and_audit(aa,J,k,ss,w)
 rows.append(make_row('signed_drive',ss,J,k,ss,aa,err,qa))
pa,pd,pf,pst=pseudo_arclength_drive(rows,(J,k,-1.),w,None,steps=12,h=.02)
assert pst['parameter_turn_observed'] and pst['max_signed_drive_including_tip']<1.
fold=refine_algebraic_fold(pa,(J,k,-1.),w,None)
assert fold['status']=='NUMERICALLY_NONDEGENERATE_LOCAL_SADDLE_NODE_CANDIDATE__NOT_INTERVAL_CERTIFICATE'
assert all(val for key,val in fold['conditions'].items() if key!='authentic_source_match')
assert not fold['conditions']['authentic_source_match']
print('SYNTHETIC_FIXED_CHARGE_AND_BOUNDARY_FALSIFIERS_PASS')
