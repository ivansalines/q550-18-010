#!/usr/bin/env python3
"""v18.53a lightweight algebraic checks, NOT an authenticated engine test."""
import json,math,sys
from pathlib import Path
import numpy as np
from run_audit import solvecandidate,xi_rhs_local,jac_fd,spectrum_match,C,K,SEEDS
from return_map import lift8,charge,wrap
R=Path(__file__).resolve().parents[1]
z=json.loads((R/'reference/v18_52z_authenticated_report.json').read_text())
for seed in SEEDS:
 src=z['results'][seed];u=np.asarray(src['newton_final_candidate'],float);J=float(src['J_reference'])
 c=solvecandidate(u,J)
 assert c['success'] and c['root_residual_max']<1e-11
 assert c['candidate_core_shift']<1e-7
 y=lift8(np.array([*c['chart_u_delta'][:7],0.]),J,0.)
 assert abs(charge(y)-J)<1e-12
 f=xi_rhs_local(c['chart_u_delta'],J)
 assert max(abs(f))<1e-11
 bad=np.array(c['chart_u_delta'],float);bad[-1]=0.
 assert max(abs(xi_rhs_local(bad,J)))>.1, 'False zero-phase branch must be rejected'
 ms=[jac_fd(lambda q:xi_rhs_local(q,J),c['chart_u_delta'],eps) for eps in C['jacobian_steps']]
 rel=np.linalg.norm(ms[1]-ms[2])/max(1.,np.linalg.norm(ms[2]))
 assert rel<1e-4, rel
 sp=spectrum_match(ms[-1],float(src['late_source_period_xi']),src['conditional_floquet'])
 assert sp['rel_eigen_discrepancy_max']<.002, sp
 # Displacing the hidden oscillator velocity is NOT a relative equilibrium.
 shifted=np.asarray(c['chart_u_delta']);shifted[4]=.07
 assert max(abs(xi_rhs_local(shifted,J)))>.001
 print(seed,'root',c['root_residual_max'],'local_spectrum_vs_parent',sp['rel_eigen_discrepancy_max'])
print('SYNTHETIC_AND_ARCHIVED_ALGEBRA_PASS')
