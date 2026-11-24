#!/usr/bin/env python3
"""Q550 v18.53a: authenticated instantaneous rotating equilibrium test.

The local algebraic root is independent of the preceding three-turn closure.
Authentic 11D RHS and clock phases are loaded only from the SHA-verified v18.51s
ancestor. This is not a global Hamiltonian or complete forward-return theorem.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,json,math,sys,traceback
from pathlib import Path
import numpy as np
from scipy.optimize import root,linear_sum_assignment
from return_map import charge,lift8,reduce8,wrap
from energy_math import sigmoid,dU,Params,energy_terms
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;OUT=ROOT/'results'
C=json.loads((ROOT/'inputs_config.json').read_text())
SEEDS=('original','autonomous');P=Params(**C['model_params']);K=float(C['v1852l_fixed_model_kappa'])

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 return h.hexdigest()

def save_report(r):
 OUT.mkdir(exist_ok=True,parents=True)
 (OUT/'report.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')

def importer(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod)
 return mod

def solvecandidate(u,J,k=K):
 """Exact stationary-branch equations for delta=pi and x=1/2, fixed J.
 The charge is imposed algebraically and never patched after integration.
 """
 def residual(a):
  l1,l2,w1=np.asarray(a,float)
  if not np.all(np.isfinite(a)) or abs(l1)>20 or abs(l2)>20:return np.array([1e6,1e6,1e6])
  r1,r2=math.exp(l1),math.exp(l2)
  w2=(r1*r1*w1-J)/(r2*r2); x=.5
  f1=w1*w1-dU(r1,P.e_reserve)/r1+k*x*r2/r1
  f2=w2*w2-dU(r2,P.e_reserve)/r2+k*x*r1/r2
  return np.array([f1,f2,w1+w2+x])
 result=root(residual,[u[0],u[3],u[2]],method='hybr',tol=1e-12)
 a=np.asarray(result.x,float); rr=residual(a)
 l1,l2,w1=a;rad1,rad2=math.exp(l1),math.exp(l2)
 w2=(rad1**2*w1-J)/rad2**2
 ueq=np.asarray([l1,0.,w1,l2,0.,0.,0.,math.pi],float)
 return {'success':bool(result.success or np.max(np.abs(rr))<C['root_residual_tol']),
         'root_residual':rr.tolist(),'root_residual_max':float(np.max(np.abs(rr))),
         'equilibrium_core':{'logr1':float(l1),'logr2':float(l2),'r1':rad1,'r2':rad2,
                             'x':.5,'omega1':float(w1),'omega2':float(w2),'J':float(J),
                             'delta':'pi_mod_2pi'},
         'core_root':a.tolist(),'chart_u_delta':ueq.tolist(),
         'candidate_core_shift':float(np.linalg.norm(a-np.array([u[0],u[3],u[2]])))}

def xi_rhs_local(y,J,k=K,delta=math.pi,omega=.844357915384206):
 """Local algebraic reduced field for synthetic & off-server self-check only.
 xi source f=(omega/x) d(physical state)/d tau, while Theta'_xi=-omega.
 The 10D source-authenticated server test is independent of this expression.
 """
 l1,v1,w1,l2,v2,z,vz,D=np.asarray(y,float)
 r1,r2=math.exp(l1),math.exp(l2);x=float(sigmoid(z));g=x*(1-x)
 w2=(r1*r1*w1-J)/r2**2
 sin,cos=math.sin(D),math.cos(D)
 dr1=r1*v1;dr2=r2*v2; dx=g*vz
 a1=(r1*w1*w1-dU(r1,P.e_reserve)-k*x*r2*cos-P.gamma_rho*dr1)/r1-v1*v1
 a2=(r2*w2*w2-dU(r2,P.e_reserve)-k*x*r1*cos-P.gamma_rho*dr2)/r2-v2*v2
 dw1=k*x*r2*sin/r1-2*v1*w1
 xdd=(-P.barrier_B*(2*x-1)/(x*(1-x))+P.mu_x*r1*r2*sin-P.gamma_x*dx)/P.m_x
 # dV_barrier/dx = B*(2x-1)/(x*(1-x)); so -dV is -B*(2x-1)/g.
 dvz=xdd/g-(1-2*x)*vz*vz
 fac=omega/x
 return fac*np.asarray([v1,a1,dw1,v2,a2,vz,dvz,w1+w2+x])

def direct_native_rhs(eng,ctx,xi,u_delta,J,g):
 """Return RHS at one instantaneous state in the physical return protocol.
 Chart uses delta=theta1+theta2-phi0-Theta_ext, not absolute theta sum.
 """
 time=float(ctx['event'])+float(xi)
 phase=float(ctx['theta'](time))
 a=np.asarray(u_delta,float).copy()
 a[7]=phase+float(eng.params.phi0)+a[7]
 y=lift8(a,J,g)
 des=eng.original_ramp.Design('v18_50p','beta:2:2',.24)
 sign=float(eng.clock.cubic_clock_sign(time,des))
 if abs(sign+1.)>1e-12:raise RuntimeError('PHASE_IS_NOT_IN_LATE_NEGATIVE_DRIVE')
 F=np.asarray(eng.physical.rhs_from_original(np.r_[y,phase],K,eng.params,eng.slope,sign,eng.support,eng.audit),float)
 if F.shape!=(11,):raise RuntimeError('UNEXPECTED_NATIVE_RHS_DIMENSION')
 f=np.asarray([F[0],F[1],F[3],F[4],F[5],F[8],F[9],F[2]+F[6]-F[10]])
 return f,F,y,phase

def jac_fd(fn,u,eps):
 x=np.asarray(u,float);D=np.empty((len(x),len(x)))
 for i in range(len(x)):
  d=np.zeros(len(x));d[i]=eps*max(1.,abs(x[i]))
  D[:,i]=(fn(x+d)-fn(x-d))/(2*d[i])
 return D

def spectrum_match(instant_D,T,archived):
 lam=np.linalg.eigvals(instant_D); multipliers=np.exp(lam*T)
 old=np.array([complex(e['real'],e['imag']) for e in archived['eigenvalues']])
 row,col=linear_sum_assignment(abs(multipliers[:,None]-old[None,:]))
 errors=np.abs(multipliers[row]-old[col]);rel=errors/np.maximum(1.,abs(old[col]))
 return {'instantaneous_generator_eigenvalues':[{'re':float(v.real),'im':float(v.imag)} for v in lam],
         'predicted_monodromy':[{'re':float(v.real),'im':float(v.imag),'mod':float(abs(v))} for v in multipliers],
         'archived_monodromy':[{'re':float(v.real),'im':float(v.imag),'mod':float(abs(v))} for v in old],
         'best_permutation':[int(i) for i in col],
         'abs_eigen_discrepancy_max':float(np.max(errors)),
         'rel_eigen_discrepancy_max':float(np.max(rel)),
         'generator_spectral_abscissa':float(np.max(lam.real)),
         'prediction_spectral_radius':float(np.max(abs(multipliers))),
         'archived_spectral_radius':float(np.max(abs(old)))}

def audit_seed(seed,eng,ctx,xi0,g,old):
 u=np.asarray(old['newton_final_candidate'],float);J=float(old['J_reference']);period=float(old['late_source_period_xi'])
 c=solvecandidate(u,J)
 ue=np.array(c['chart_u_delta'],float)
 # Absolute phase: residual of the v18.52z root on delta=pi branch.
 p0=float(ctx['theta'](ctx['event']+xi0))
 delta=wrap(u[7]-float(eng.params.phi0)-p0)
 phase_branch=wrap(delta-math.pi)
 native_u=u.copy();native_u[7]=delta
 candidate_f,F_c,y_c,phase=direct_native_rhs(eng,ctx,xi0,native_u,J,g)
 root_f,F_e,y_e,_=direct_native_rhs(eng,ctx,xi0,ue,J,g)
 e=energy_terms(np.r_[y_e,phase],K,P)
 ener=float(e['E'])
 ch=charge(y_e)
 gp=[-.91,.51,1.62]
 gauge_errors=[]
 for eps in gp:
  f,_,_,_=direct_native_rhs(eng,ctx,xi0,ue,J,g+2*eps)
  gauge_errors.append(float(np.max(np.abs(root_f-f))))
 h=C['jacobian_steps'];M=[jac_fd(lambda v:direct_native_rhs(eng,ctx,xi0,v,J,g)[0],ue,st) for st in h]
 jm=float(np.linalg.norm(M[0]-M[1])/max(1.,np.linalg.norm(M[1])))
 js=float(np.linalg.norm(M[1]-M[2])/max(1.,np.linalg.norm(M[2])))
 sp=spectrum_match(M[-1],period,old['conditional_floquet'])
 return {'seed':seed,'J':J,'kappa':K,'t0_xi':xi0,'period_xi':period,
         'phase_external_at_anchor':p0, 'candidate_alpha':float(u[7]),
         'candidate_delta':float(delta),'candidate_phase_branch_distance_pi':float(abs(phase_branch)),
         'candidate_instantaneous_rhs':candidate_f.tolist(),
         'candidate_instantaneous_rhs_max':float(np.max(np.abs(candidate_f))),
         'algebraic_branch_root':c, 'native_root_reduced_rhs':root_f.tolist(),
         'native_root_rhs_max':float(np.max(np.abs(root_f))),
         'native_root_charge':ch,'native_root_charge_error':float(abs(ch-J)),
         'native_root_F11':F_e.tolist(),
         'native_root_angvel_sum_vs_clock':float(F_e[2]+F_e[6]-F_e[10]),
         'native_root_energy':ener,'gauge_errors':gauge_errors,
         'signed_clock_11d_residual':float(abs(F_e[10]+float(eng.omega))),
         'source_phase_finite_difference_rate':float((ctx['theta'](ctx['event']+xi0+1e-5)-ctx['theta'](ctx['event']+xi0-1e-5))/(2e-5)),
         'jacobian_fd_rel_disagreement_large_mid':jm,
         'jacobian_fd_rel_disagreement_mid_small':js,
         'floquet_vs_exp_instantaneous_jacobian':sp,
         'equilibrium_reduced_chart':ue.tolist(),
         'gauge_at_parent_native_start':g,
         'instantaneous_root_pass':bool(c['success'] and c['root_residual_max']<C['root_residual_tol'] and
           c['candidate_core_shift']<C['root_change_tol'] and abs(phase_branch)<C['alpha_phase_alignment_tol'] and
           np.max(np.abs(candidate_f))<C['instantaneous_candidate_rhs_tol'] and
           abs(F_e[10]+float(eng.omega))<1e-10 and
           np.max(np.abs(root_f))<C['instantaneous_algebraic_root_rhs_tol'] and abs(ch-J)<C['charge_tol'] and max(gauge_errors)<C['gauge_tol']),
         'local_spectrum_pass':bool(jm<C['jacobian_rel_convergence_tol'] and js<C['jacobian_rel_convergence_tol'] and sp['rel_eigen_discrepancy_max']<C['monodromy_match_tol'])}

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--native-v1852l',type=Path,default=BASE/'q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results/extended_horizon_10d_states.npz')
 ap.add_argument('--require-parent',action='store_true');ap.add_argument('--require-native',action='store_true')
 args=ap.parse_args(); OUT.mkdir(exist_ok=True,parents=True)
 zfile=ROOT/'reference/v18_52z_authenticated_report.json'
 if sha(zfile)!=C['v1852z_report_sha256']:raise RuntimeError('AUTHENTICATED_Z_REPORT_HASH_MISMATCH')
 z=json.loads(zfile.read_text())
 if z.get('scientific_decision')!=C['parent_archived_v1852z_status']:raise RuntimeError('UNEXPECTED_AUTHENTICATED_Z_STATUS')
 report={'version':'v18.53a','title':C['title'],'parent_z_sha256':sha(zfile),'source_requirements_met':False,
         'status':'AUTHENTIC_INPUTS_UNAVAILABLE__OFFLINE_ANALYTIC_DIAGNOSTICS_ONLY',
         'global_conservative_dynamics_certified':False,'full_historical_forward_return_cycle_certified':False,
         'results':{}}
 for seed in SEEDS:
  old=z['results'][seed]
  cc=solvecandidate(np.asarray(old['newton_final_candidate'],float),float(old['J_reference']))
  report['results'][seed]={'offline_equilibrium_core':cc,'parent_three_turn_closed':bool(old['relative_fixed_point_accepted'])}
 save_report(report)
 stage=importer(ROOT/'reference/v18_52v_r1_run_audit.py','q550_stage_53a')
 parent=stage.find_parent(args.parent_v1851s)
 if parent is None or not args.native_v1852l.is_file():
  if (args.require_parent and parent is None) or (args.require_native and not args.native_v1852l.is_file()):raise RuntimeError('AUTHENTIC_Q550_INPUTS_REQUIRED')
  return
 if sha(parent)!=C['v1851s_sha256'] or sha(args.native_v1852l)!=C['v1852l_10d_sha256']:
  raise RuntimeError('AUTHENTIC_SOURCE_OR_NATIVE_SHA_MISMATCH')
 with np.load(args.native_v1852l,allow_pickle=False) as p:
  xi=np.asarray(p['common_xi'],float)
  Y={s:np.asarray(p[s+'__reference__Y_tau'][0,:10],float) for s in SEEDS}
  if len(xi)!=65537 or abs(xi[0]-12.8)>1e-10:raise RuntimeError('NATIVE_GRID_CHANGED')
 temp=None
 try:
  temp,eng,_,ks,contexts,source=stage.stage_parent(parent,xi[:2]);t0=float(xi[0])
  if any(abs(float(ks[s])-K)>1e-14 for s in SEEDS):raise RuntimeError('MODEL_KAPPA_CHANGED')
  params={key:float(getattr(eng.params,key)) for key in C['model_params']}
  if any(abs(params[k]-C['model_params'][k])>1e-12 for k in params):raise RuntimeError('AUTHENTIC_MODEL_PARAMETERS_MISMATCH')
  report.update({'source_requirements_met':True,'v1851s_sha256':sha(parent),'v1852l_sha256':sha(args.native_v1852l),
                 'verified_model_params':params,'source_genealogy':{k:source[k] for k in ('authenticated_parent','rhs_from_original','cubic_clock_sign','cubic_phase_excursion','make_context')}})
  for seed in SEEDS:
   g=float(Y[seed][2]-Y[seed][6]);r=audit_seed(seed,eng,contexts[seed],t0,g,z['results'][seed]);
   report['results'][seed]=r;save_report(report)
   print(seed,'native root residual',r['native_root_rhs_max'],'local eig error',r['floquet_vs_exp_instantaneous_jacobian']['rel_eigen_discrepancy_max'],flush=True)
  a,b=(report['results'][s] for s in SEEDS)
  # Both families share J,kappa; compare intrinsic delta-locked state rather than absolute phase/gauge.
  u1=np.asarray(a['equilibrium_reduced_chart']);u2=np.asarray(b['equilibrium_reduced_chart'])
  diff=u1-u2;diff[-1]=wrap(diff[-1]);equiv=float(np.linalg.norm(diff))
  source_delta=wrap((a['candidate_alpha']-b['candidate_alpha'])-(a['phase_external_at_anchor']-b['phase_external_at_anchor']))
  Jdiff=abs(a['J']-b['J'])
  report['cross_seed']={'J_difference':Jdiff,'intrinsic_locked_chart_difference':equiv,
    'candidate_phase_offset_invariance_defect':float(source_delta),
    'phase_external_difference':a['phase_external_at_anchor']-b['phase_external_at_anchor'],
    'candidate_alpha_difference':a['candidate_alpha']-b['candidate_alpha'],
    'equivalent_on_fixed_charge_and_phase_frame':bool(equiv<C['cross_seed_equivalence_tol'] and Jdiff<C['charge_tol'] and abs(source_delta)<C['cross_seed_equivalence_tol'])}
  both=all(report['results'][s]['instantaneous_root_pass'] for s in SEEDS)
  spect=all(report['results'][s]['local_spectrum_pass'] for s in SEEDS)
  report['status']=('AUTHENTIC_INSTANTANEOUS_PHASE_LOCKED_RELATIVE_EQUILIBRIUM_AND_CROSS_SEED_MATCH_WITH_GENERATOR_FLOQUET_AGREEMENT' if both and spect and report['cross_seed']['equivalent_on_fixed_charge_and_phase_frame'] else
   'AUTHENTIC_SOURCE_AUDITED__PARTIAL_OR_FAILED_PHASE_LOCK_AND_LOCAL_SPECTRAL_MATCH')
  report['instantaneous_equilibrium_pass']=both;report['generator_floquet_agreement_pass']=spect
  lines=['# Q550 v18.53a — phase-locked relative equilibrium audit','', '**'+report['status']+'**','',
   'Instantaneous source RHS is authenticated against the v18.51s ancestor, and the candidate is SHA-pinned to v18.52z.','',
   'The local algebraic root uses Delta=pi, x=1/2, fixed J, and a late negative external clock. The parent map and this root are independent computations.','']
  for s in SEEDS:
   rr=report['results'][s];p=rr['algebraic_branch_root']['equilibrium_core']
   lines+=['## '+s,'',f"- Radial equilibrium r1={p['r1']:.14g}, r2={p['r2']:.14g}; omega1={p['omega1']:.14g}, omega2={p['omega2']:.14g}",
    f"- Candidate instantaneous RHS residual {rr['candidate_instantaneous_rhs_max']:.12g}",
    f"- Independently rooted RHS residual {rr['native_root_rhs_max']:.12g}",
    f"- Candidate delta-minus-pi residual {rr['candidate_phase_branch_distance_pi']:.12g}",
    f"- Relative Floquet-vs-exp(local generator) max error {rr['floquet_vs_exp_instantaneous_jacobian']['rel_eigen_discrepancy_max']:.12g}",
    f"- Instantaneous algebraic equilibrium certified: {rr['instantaneous_root_pass']}",f"- Local spectral comparison passed: {rr['local_spectrum_pass']}",'']
  lines+=['## Cross-seed equivalence',json.dumps(report['cross_seed'],indent=2),'',
   '## Strict interpretation','The result is conditional on the *constant negative externally prescribed return drive*, exact fixed J and quotient of counter-rotation. It is neither a full historical forward-return cycle nor a global conservative Hamiltonian statement.']
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
  save_report(report)
 finally:
  if temp is not None:temp.cleanup()

if __name__=='__main__':
 try:main()
 except Exception as exc:
  OUT.mkdir(parents=True,exist_ok=True);(OUT/'audit_failure.txt').write_text(traceback.format_exc());print('V18_53A_FAILURE',repr(exc),flush=True);sys.exit(1)
