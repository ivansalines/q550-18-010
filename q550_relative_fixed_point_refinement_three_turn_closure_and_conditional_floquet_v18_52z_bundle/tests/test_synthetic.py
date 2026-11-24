#!/usr/bin/env python3
"""v18.52z local synthetic falsification (not authentic Q550 evidence)."""
from pathlib import Path
import json, math,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from return_map import charge, reduce8, lift8, gauge, tangent_basis, chart_diff
from run_audit import newton_refine, three_turn_gate, spectrum
from types import SimpleNamespace

root=Path(__file__).resolve().parents[1]
cfg=json.loads((root/'inputs_config.json').read_text())
r=json.loads((root/'reference/v18_52y_authenticated_report.json').read_text())
assert r['scientific_decision'].endswith('NO_RELATIVE_CYCLE_CERTIFICATE')
for seed in ('original','autonomous'):
 x=np.asarray(r['results'][seed]['closure_search']['candidate_u'],float)
 assert x.size==8 and np.isfinite(x).all()
 assert r['results'][seed]['closure_search']['accepted_as_relative_fixed_point'] is False

Y=np.array([.2,.01,3.8,-.44,-1.3,.05,-1.1,.13,.04,-.02]);J=charge(Y);u=reduce8(Y);g=gauge(Y)
B=tangent_basis(Y); assert np.linalg.matrix_rank(B)==8
for q in range(8):
 v=np.eye(8)[q]*1e-5
 assert abs(charge(lift8(u+v,J,g))-J)<1e-12
 assert abs(charge(lift8(u-v,J,g))-J)<1e-12

A=np.diag([.35,.42,.5,.61,.69,.72,.83,.91]);target=u+np.array([.01,.02,-.03,.015,-.022,.013,.017,-.02])
args=SimpleNamespace(max_newton=7,closure_tol=1e-9,newton_fd=5e-5,max_weighted_step=.12,max_backtracks=10,
 two_turn_tol=3e-9,three_turn_tol=5e-9,composition_tol=1e-8,charge_tol=1e-8)
class SyntheticMap:
 def __init__(self):self.j=J;self.max_j_drift=0.
 def map(self,x,turns=1):return target+np.linalg.matrix_power(A,turns)@(np.asarray(x)-target)
sm=SyntheticMap();q,err,h=newton_refine(sm,u,args)
assert err<1e-9 and np.linalg.norm(q-target)<1e-8
v=three_turn_gate(sm,q,1,args);assert v['pass']

# A false phase recurrence with a hidden velocity translation must NEVER pass.
class DriftMap(SyntheticMap):
 def map(self,x,turns=1):
  y=np.asarray(x,float).copy();y[4]+=.17*turns;return y
bad=DriftMap();z,b,history=newton_refine(bad,u,args)
assert b>.1 and not three_turn_gate(bad,z,1,args)['pass']

# Stability is a numerical eigenvalue question only AFTER exact fixed-point acceptance.
N=A.copy();N[0,0]=1.0
class Neutral(SyntheticMap):
 def map(self,x,turns=1):return target+np.linalg.matrix_power(N,turns)@(np.asarray(x)-target)
ind=Neutral();sp=spectrum(ind,ind,target,SimpleNamespace(jac_steps=[1e-4,3e-5,1e-5],jac_tol=.005,attraction_margin=.001))
assert sp['jacobian_numerically_converged']
assert sp['conditional_leaf_attractor'] is False

from phase_energy import exact_identity_probe
from energy_math import Params
par=Params();q=np.array([.3,.02,.1,-.45,-1.9,.01,.2,.12,-.04,.015,0.])
for rate in (-.3,.1,.7):
 assert abs(exact_identity_probe(q,.1666331385421536,par,.15,rate))<2e-12
print('SYNTHETIC_PASS: historical candidate lineage, fixed-J 8D basis, three-turn composition, Newton root, hidden-drift falsification, neutral Floquet exclusion')
