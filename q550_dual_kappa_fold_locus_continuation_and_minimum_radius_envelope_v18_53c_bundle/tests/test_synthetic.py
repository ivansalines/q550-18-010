#!/usr/bin/env python3
"""Q550 v18.53c synthetic falsifiers, source-pinned algebra and derivative checks.
Do not interpret these local checks as original-engine authentication.
"""
import math,json
import numpy as np
from run_audit import (C,ROOT,sha,B,core_jacobian,augmented,fold_correct,fold_metrics,
   kappa_sensitivity,scan_values,sampled_envelope)
J=-1.441016445824766
k=float(C['v1852l_fixed_model_kappa']);om=.844357915384206
p=ROOT/'reference/v18_53b_authenticated_fold.json'
assert sha(p)==C['v1853b_fold_sha256']
fold=json.loads(p.read_text());z0=np.r_[fold['core_a'],fold['s_fold']]
assert np.max(abs(augmented(z0,J,k)))<2e-7
# Analytic core Jacobian must agree with an independent finite-difference derivative.
rng=np.random.default_rng(1853)
worst=0.
for i in range(64):
 a=z0[:3]+rng.normal(size=3)*np.array([.08,.06,.04]);kk=k+rng.uniform(-.008,.008);ss=z0[3]+rng.uniform(-.03,.03)
 M=core_jacobian(a,J,kk,ss);h=1e-5;N=np.zeros((3,3))
 for j in range(3):
  ee=np.zeros(3);ee[j]=h
  N[:,j]=(B.radial_residual(a+ee,J,kk,ss)-B.radial_residual(a-ee,J,kk,ss))/(2*h)
 worst=max(worst,float(np.max(abs(M-N))))
assert worst<1e-6,worst
z,d=fold_correct(z0,J,k,om)
assert d['radial_core_nondegenerate'] and d['one_real_mode_near_zero'] and d['other_seven_stable']
# Both distinct prescriptions must have independently solved fold roots.
kap=C['kappa_seeds'];oa,od=fold_correct(z,J,float(kap['original']),om);aa,ad=fold_correct(z,J,float(kap['autonomous']),om)
assert abs(oa[3]-aa[3])>1e-5
assert abs(math.exp(oa[1])-math.exp(aa[1]))>1e-6
# Tangent differential is checked against two physically re-solved folds.
sens=kappa_sensitivity(z,J,k)
for h in [1e-5,3e-6]:
 zp,_=fold_correct(z,J,k+h,om);zm,_=fold_correct(z,J,k-h,om)
 empirical=(zp-zm)/(2*h)
 assert abs(empirical[3]-sens['ds_dk'])<3e-4,(h,empirical[3],sens['ds_dk'])
 assert abs((math.exp(zp[1])-math.exp(zm[1]))/(2*h)-sens['dr2_dk'])<3e-4
values=np.r_[np.linspace(float(kap['autonomous']),float(kap['original']),11),k]
rows,states,stops=scan_values(values,k,z,J,om,None,'synthetic')
assert not stops and len(rows)==12,len(rows)
env=sampled_envelope(rows,float(kap['original']),float(kap['autonomous']))
assert env['monotonic_sampled'] and env['radius_component_switches']==0
assert env['sampled_min_radius']<env['sampled_max_of_min_radius']
# Deliberate falsification: a stationary phase does not ensure no hidden-radius derivative.
fake=[dict(item) for item in rows]
for i,r in enumerate(fake):r['dr2_dk']=math.sin(2*math.pi*i/(len(fake)-1))
env2=sampled_envelope(fake,float(kap['original']),float(kap['autonomous']))
assert not env2['monotonic_sampled'] and len(env2['detected_derivative_zero_brackets'])>=1
print('SYNTHETIC_DUAL_KAPPA_FOLD_AND_MINMAX_FALSIFIERS_PASS',len(rows),'analytic_jac_err',worst)
