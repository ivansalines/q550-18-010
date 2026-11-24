#!/usr/bin/env python3
"""Q550 v18.52y: authenticated fixed-J 8D return map, real IVPs on server.

CRITICAL distinctions:
  • phase is the EXTERNAL authenticated late-return forcing, not native Y[:,10];
  • the eight columns are actual fixed-J tangent directions, not historical kicks;
  • one-period DP at a non-fixed point is NOT a Floquet matrix;
  • Newton failure and a single apparent phase recurrence are NOT a cycle;
  • attractivity is judged only after physical quotient closure and converged DP.
No long integrations are performed by ChatGPT: script to be run on user's server.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,sys,traceback
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from return_map import TWOPI,LABELS,wrap,charge,reduce8,gauge,lift8,tangent_basis,chart_diff,scaled_norm,fixed_charge_fd,damped_fixed_point

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;OUT=ROOT/'results'
C=json.loads((ROOT/'inputs_config.json').read_text())
SEEDS=('original','autonomous')
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
L52='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'

def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()

def save(name,obj):
 (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def rows(name,data):
 if not data:return
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,list(data[0]));w.writeheader();w.writerows(data)

def module_at(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

class AuthenticatedMap:
 def __init__(self,eng,context,k,t0,period,rtol,atol,maxstep,cache=True):
  self.eng=eng;self.ctx=context;self.k=float(k);self.t0=float(t0);self.period=float(period)
  self.rtol=float(rtol);self.atol=float(atol);self.maxstep=float(maxstep)
  self.design=eng.original_ramp.Design('v18_50p','beta:2:2',.24)
  self.cache={} if cache else None;self.eval_count=0;self.nfev=0;self.max_j_drift=0.
  self.g=0.;self.j=0.
 def set_leaf(self,J,g):
  self.j=float(J);self.g=float(g)
  if self.cache is not None:self.cache.clear()
 def source(self,xi,y):
  t=float(self.ctx['event'])+float(xi)
  theta=float(self.ctx['theta'](t))
  s=float(self.eng.clock.cubic_clock_sign(t,self.design))
  if abs(s+1.)>1e-9:raise RuntimeError('LEFT_CONSTANT_NEGATIVE_LATE_DRIVE')
  y=np.asarray(y,float)
  if not np.all(np.isfinite(y)):raise FloatingPointError('NONFINITE_10D_STATE')
  with np.errstate(over='raise',invalid='raise',divide='raise'):
   F=np.asarray(self.eng.physical.rhs_from_original(np.r_[y,theta],self.k,self.eng.params,self.eng.slope,s,self.eng.support,self.eng.audit),float)
  if F.shape!=(11,) or not np.all(np.isfinite(F)):raise RuntimeError('NATIVE_PHYSICAL_RHS_INVALID')
  return F[:10]
 def integrate(self,y0,turns=1):
  t1=self.t0+float(turns)*self.period
  sol=solve_ivp(self.source,(self.t0,t1),np.asarray(y0,float),method='DOP853',rtol=self.rtol,atol=self.atol,max_step=self.maxstep,t_eval=[t1])
  self.eval_count+=1;self.nfev+=sol.nfev
  if not sol.success or sol.y.shape!=(10,1):raise RuntimeError('SOLVER_FAILURE '+str(sol.message))
  y1=sol.y[:,-1]
  drift=abs(charge(y1)-charge(y0));self.max_j_drift=max(drift,self.max_j_drift)
  if drift>2e-7:raise RuntimeError('CHARGE_LOSS_IN_IVP '+str(drift))
  return y1
 def map(self,u,turns=1):
  u=np.asarray(u,float);key=None
  if self.cache is not None:
   key=(int(turns),u.tobytes())
   if key in self.cache:return self.cache[key].copy()
  y0=lift8(u,self.j,self.g)
  q=reduce8(self.integrate(y0,turns))
  if self.cache is not None:self.cache[key]=q.copy()
  return q

def phase_period(eng,ctx,t0):
 des=eng.original_ramp.Design('v18_50p','beta:2:2',.24)
 omega=float(eng.omega)
 if omega<=0:raise RuntimeError('SOURCE_OMEGA_NOT_POSITIVE')
 period=TWOPI/omega
 ts=np.linspace(t0,t0+2*period,9)
 signs=[float(eng.clock.cubic_clock_sign(float(ctx['event'])+float(t),des)) for t in ts]
 if any(abs(s+1.)>1e-10 for s in signs):raise RuntimeError('NOT_LATE_CONSTANT_MINUS_ONE_CLOCK')
 ph=lambda t:float(ctx['theta'](float(ctx['event'])+float(t)))
 first=ph(t0+period)-ph(t0)
 second=ph(t0+2*period)-ph(t0)
 if abs(first+TWOPI)>1e-7 or abs(second+2*TWOPI)>1e-7:
  raise RuntimeError('SOURCE_EXTERNAL_PHASE_PERIOD_INCONSISTENT '+repr((first,second)))
 return float(period),{'source_omega':omega,'source_sign_samples':signs,'phase_one_turn':first,'phase_two_turn':second,'period_xi':period}

def native_y_at(xi,Y,t):
 idx=np.searchsorted(xi,t)
 if idx==0:return Y[0,:10].copy()
 if idx>=len(xi):return Y[-1,:10].copy()
 a=(t-xi[idx-1])/(xi[idx]-xi[idx-1])
 return (1-a)*Y[idx-1,:10]+a*Y[idx,:10]

def basis_checks(y):
 j=charge(y);u=reduce8(y);g=gauge(y)
 B=tangent_basis(y)
 tangent_test=[]
 for q in range(8):
  eps=1e-6
  v=np.zeros(8);v[q]=eps
  finite=(lift8(u+v,j,g)-lift8(u-v,j,g))/(2*eps)
  tangent_test.append(float(np.linalg.norm(finite-B[:,q])))
 grad=np.zeros(10)
 grad[0]=2*math.exp(2*y[0])*y[3];grad[3]=math.exp(2*y[0])
 grad[4]=-2*math.exp(2*y[4])*y[7];grad[7]=-math.exp(2*y[4])
 charge_tangent=grad@B
 sv=np.linalg.svd(B,compute_uv=False)
 return {'basis_shape':[10,8],'basis_rank':int(np.linalg.matrix_rank(B)),
  'basis_singular_values':sv.tolist(), 'max_analytic_vs_fd_column_error':max(tangent_test),
  'max_abs_charge_derivative_on_basis':float(np.max(np.abs(charge_tangent))),
  'J_reference':j,'gauge_reference':g}

def complex_json(z):return {'real':float(z.real),'imag':float(z.imag),'abs':float(abs(z))}

def analyze_seed(seed,engine,ctx,xi,Y,cfg,args,summary_rows,deriv_rows):
 t0=float(xi[0]);period,periodcheck=phase_period(engine,ctx,t0)
 u0=reduce8(Y[0,:10]);y0=Y[0,:10];j=charge(y0)
 b=basis_checks(y0)
 if b['basis_rank']!=8 or b['max_abs_charge_derivative_on_basis']>1e-9 or b['max_analytic_vs_fd_column_error']>2e-8:
  raise RuntimeError('FIXED_CHARGE_TANGENT_BASIS_NOT_VERIFIED '+seed)
 obj=AuthenticatedMap(engine,ctx,cfg['v1852l_fixed_model_kappa'],t0,period,args.rtol,args.atol,args.maxstep)
 obj.set_leaf(j,gauge(y0))
 # First check: genuinely reproduce historical one-turn native trajectory.
 p0=obj.map(u0); ynat=native_y_at(xi,Y,t0+period);pnative=reduce8(ynat)
 native_delta=chart_diff(pnative,p0)
 native_error=float(np.linalg.norm(native_delta)/max(1.,np.linalg.norm(pnative)))
 if native_error>cfg['native_one_turn_validation_rel_tol']:
  raise RuntimeError('INTEGRATOR_DOES_NOT_REPRODUCE_AUTHENTIC_ONE_TURN '+seed+' err='+str(native_error))
 # Verify quotient gauge invariance using original exact J at same physical leaf.
 control=AuthenticatedMap(engine,ctx,cfg['v1852l_fixed_model_kappa'],t0,period,args.rtol,args.atol,args.maxstep,cache=False)
 control.set_leaf(j,gauge(y0)+.341)
 gauged=control.map(u0)
 gauge_error=float(np.linalg.norm(chart_diff(p0,gauged)))
 if gauge_error>cfg['gauge_invariance_tol']:
  raise RuntimeError('COUNTERROTATION_NOT_A_SYMMETRY_FOR_NATIVE_RETURN '+seed+' err='+str(gauge_error))
 # First vs second authentic extended return, not inferred by naive phase interpolation.
 p2=obj.map(u0,turns=2)
 p1_iter=obj.map(p0)
 semigroup_err=float(np.linalg.norm(chart_diff(p1_iter,p2)))
 # Note: map P(u) is a stroboscopic map with phase periodicity. P^2=P o P.
 if semigroup_err>5e-6:raise RuntimeError('ONE_TURN_MAP_NOT_AUTONOMOUS_ON_LATE_BRANCH '+seed+' err='+str(semigroup_err))
 res0=chart_diff(u0,p0);res1=chart_diff(p0,p2)
 # Full independent 8-column fixed-leaf derivative at historical NONCLOSED point.
 D1,c1=fixed_charge_fd(obj.map,u0,args.fd1)
 D2,c2=fixed_charge_fd(obj.map,u0,args.fd2)
 dnorm=float(np.linalg.norm(D2))
 convergence=float(np.linalg.norm(D1-D2)/max(1.,dnorm))
 singular=np.linalg.svd(D2,compute_uv=False)
 eigs=np.linalg.eigvals(D2)
 base={'seed':seed,'phase_period':periodcheck,'J_fixed':j,'native_one_turn_relative_error':native_error,
  'counterrotation_gauge_equivariance_error':gauge_error,'P_squared_vs_extended_integration_error':semigroup_err,
  'first_return_residual_chart_norm':float(np.linalg.norm(res0)),
  'second_return_residual_chart_norm':float(np.linalg.norm(chart_diff(u0,p2))),
  'consecutive_return_displacement_1':float(np.linalg.norm(res0)),
  'consecutive_return_displacement_2':float(np.linalg.norm(res1)),
  'basis':b,'two_actual_phase_turns_integrated':True,'DP_at_reference_exists':True,
  'historical_reference_is_closed':bool(np.linalg.norm(res0)<args.closure_tol),
  'reference_D_P_singular_values_descriptive':singular.tolist(),
  'reference_DP_eigenvalues_not_floquet':[complex_json(q) for q in eigs],
  'finite_difference_steps':[args.fd1,args.fd2],
  'relative_two_step_jacobian_difference':convergence,
  'fixed_charge_perturbation_directions':8,
  'initial_chart':u0.tolist(),'P_initial_chart':p0.tolist(),'P2_initial_chart':p2.tolist(),
  'chart_displacement_after_one_turn':{s:float(v) for s,v in zip(LABELS,res0)},
  'chart_displacement_after_two_turns':{s:float(v) for s,v in zip(LABELS,chart_diff(u0,p2))}}
 for col in range(8):
  h=args.fd2*max(1.,abs(u0[col]));v=np.zeros(8);v[col]=h
  chargep=charge(lift8(u0+v,j,obj.g));chargem=charge(lift8(u0-v,j,obj.g))
  deriv_rows.append({'seed':seed,'column':col,'label':LABELS[col],'fd_h':h,
   'fixed_leaf_charge_plus_error':chargep-j,'fixed_leaf_charge_minus_error':chargem-j,
   'column_DP_norm':float(np.linalg.norm(D2[:,col])),
   'column_DP_difference_two_steps':float(np.linalg.norm(D1[:,col]-D2[:,col])),
   'not_floquet_until_fixed_point':True})
 # Bounded search, NEVER retroactively treat near-miss as a cycle.
 if args.skip_closure:
  search={'converged':False,'best_residual':float(np.linalg.norm(res0)), 'best_u':u0,'history':[], 'skipped':True}
 else:
  search=damped_fixed_point(obj.map,u0,steps=args.max_newton,backtracks=cfg['max_line_search_backtracks'],
   weighted_step=cfg['newton_max_weighted_step'],threshold=args.closure_tol,fd_eps=args.fd1)
 uc=np.asarray(search['best_u'],float)
 pp=obj.map(uc);cc=chart_diff(uc,pp)
 closure_res=float(np.linalg.norm(cc))
 closure_relative=float(np.linalg.norm(cc)/max(1.,np.linalg.norm(uc)))
 repeat=obj.map(uc,turns=2)
 repeat_err=float(np.linalg.norm(chart_diff(uc,repeat))/max(1.,np.linalg.norm(uc)))
 # Independent second step of the same candidate (no reused native trajectory).
 two_step_consistency=float(np.linalg.norm(chart_diff(obj.map(pp),repeat)))
 closure_pass=bool(closure_res<args.closure_tol and repeat_err<2*args.closure_tol and two_step_consistency<5e-6 and obj.max_j_drift<1e-8)
 base['closure_search']={'search_converged':bool(search['converged']),
  'closure_8D_absolute_norm':closure_res,'closure_8D_normalized_norm':closure_relative,
  'two_turn_reclosure_normalized_norm':repeat_err,
  'P_composed_vs_integrated_two_turn':two_step_consistency,
  'search_history':search['history'],'candidate_u':uc.tolist(),
  'accepted_as_relative_fixed_point':closure_pass,
  'physical_10D_full_orbit_not_asserted':True}
 if closure_pass:
  Dc1,_=fixed_charge_fd(obj.map,uc,args.fd1)
  Dc2,_=fixed_charge_fd(obj.map,uc,args.fd2)
  jacerr=float(np.linalg.norm(Dc1-Dc2)/max(1.,np.linalg.norm(Dc2)))
  eig=np.linalg.eigvals(Dc2)
  numerically_converged=jacerr<cfg['floquet_jacobian_convergence_rel_threshold']
  base['closure_search']['fixed_point_jacobian_two_step_rel_error']=jacerr
  base['closure_search']['conditional_floquet_multipliers']=[complex_json(q) for q in eig] if numerically_converged else None
  base['closure_search']['conditional_leaf_attracting']=bool(max(abs(eig))<1-1e-3) if numerically_converged else None
  base['closure_search']['spectral_radius']=float(max(abs(eig))) if numerically_converged else None
  base['closure_search']['floquet_certificate_conditional']=numerically_converged
 else:
  base['closure_search']['conditional_floquet_multipliers']=None
  base['closure_search']['floquet_certificate_conditional']=False
 base['statistics']={'physical_10D_ivps':obj.eval_count+control.eval_count,'total_solver_nfev':obj.nfev+control.nfev,
  'max_native_integration_J_drift':max(obj.max_j_drift,control.max_j_drift),'fixed_charge_jacobian_probe_evaluations':32}
 summary_rows.append({'seed':seed,'one_turn_norm':float(np.linalg.norm(res0)),'two_turn_norm':float(np.linalg.norm(chart_diff(u0,p2))),
  'native_replay_error':native_error,'gauge_error':gauge_error,
  'max_abs_J_drift':base['statistics']['max_native_integration_J_drift'],
  'eight_directions_certified':True,'jacobian_convergence':convergence,
  'closure_search_residual':closure_res,'relative_fixed_point_certified':closure_pass,
  'floquet_spectrum_reported':base['closure_search']['conditional_floquet_multipliers'] is not None})
 return base

def render(results):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return []
 plots=[]
 for seed,r in results.items():
  his=r['closure_search']['search_history']
  fig,ax=plt.subplots(2,1,figsize=(8.5,6.5),constrained_layout=True)
  vals=[r['first_return_residual_chart_norm'],r['second_return_residual_chart_norm'],r['closure_search']['closure_8D_absolute_norm']]
  ax[0].bar(['Native P(u)-u','Native P²(u)-u','Best closure'],vals)
  ax[0].set_ylabel('Descriptive 8D chart norm (mixed coordinates)');ax[0].grid(axis='y',alpha=.3)
  ax[1].plot(range(8),r['reference_D_P_singular_values_descriptive'],'o-',label='One-step DP singular values (not Floquet)')
  ax[1].set_yscale('log');ax[1].set_xlabel('Singular value index');ax[1].set_ylabel('Singular value');ax[1].grid(alpha=.3);ax[1].legend(fontsize=8)
  fig.suptitle('Q550 v18.52y '+seed+' — real two-turn fixed-J return map')
  name=seed+'_fixed_charge_return_map.png';fig.savefig(OUT/name,dpi=130);plt.close(fig);plots.append(name)
 return plots

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--native-v1852l',type=Path,default=BASE/L52/'results/extended_horizon_10d_states.npz')
 ap.add_argument('--require-parent',action='store_true');ap.add_argument('--require-native',action='store_true')
 ap.add_argument('--seeds',choices=['both','original','autonomous'],default='both')
 ap.add_argument('--skip-closure',action='store_true')
 ap.add_argument('--max-newton',type=int,default=C['max_newton_iterations'])
 ap.add_argument('--fd1',type=float,default=C['finite_difference_epsilon_1'])
 ap.add_argument('--fd2',type=float,default=C['finite_difference_epsilon_2'])
 ap.add_argument('--rtol',type=float,default=C['solver_rtol'])
 ap.add_argument('--atol',type=float,default=C['solver_atol'])
 ap.add_argument('--maxstep',type=float,default=C['solver_max_step_xi'])
 ap.add_argument('--closure-tol',type=float,default=C['closure_residual_threshold'])
 args=ap.parse_args();OUT.mkdir(exist_ok=True,parents=True)
 if args.fd1<=0 or args.fd2<=0 or args.rtol<=0 or args.maxstep<=0:raise ValueError('INVALID_SOLVER_CONFIGURATION')
 candidates=([args.parent_v1851s] if args.parent_v1851s else [])+[BASE/(S51+'_server_sendback.zip'),Path('/root/elab1/hysteresis')/(S51+'_server_sendback.zip'),Path('/root/Downloads')/(S51+'_server_sendback.zip')]
 parent=next((p for p in candidates if p and p.is_file()),None)
 report={'version':'v18.52y','title':C['title'],'parent_available':bool(parent),
  'native_available':args.native_v1852l.is_file(),'solver_configuration':{'rtol':args.rtol,'atol':args.atol,'max_step':args.maxstep,'fd_eps':[args.fd1,args.fd2], 'max_newton':args.max_newton},
  'full_global_conservative_dynamics_claimed':False,'global_limit_cycle_certificate':False,'original_historical_mixed_charge_kicks_not_used_for_DP':True}
 if not parent or not args.native_v1852l.is_file():
  report['scientific_decision']='AUTHENTIC_ANCESTOR_MISSING__SYNTHETIC_ONLY'
  save('report.json',report)
  if args.require_parent or args.require_native:raise RuntimeError('REQUIRED_AUTHENTIC_ANCESTOR_MISSING')
  return
 if digest(parent)!=C['v1851s_sha256']:raise RuntimeError('AUTHENTIC_PARENT_SHA256_MISMATCH')
 if digest(args.native_v1852l)!=C['v1852l_10d_sha256']:raise RuntimeError('AUTHENTIC_NATIVE_SHA256_MISMATCH')
 v=module_at(ROOT/'reference/v18_52v_r1_run_audit.py','q550_v1852y_auth_parent')
 temp=None
 try:
  native=np.load(args.native_v1852l,allow_pickle=False)
  with native as z:
   xi=np.asarray(z['common_xi'],float)
   if len(xi)!=65537 or abs(xi[0]-12.8)>1e-8 or abs(xi[-1]-25.6)>1e-8:
    raise RuntimeError('UNEXPECTED_NATIVE_XI_GRID')
   sel=SEEDS if args.seeds=='both' else (args.seeds,)
   series={seed:np.asarray(z[seed+'__reference__Y_tau'],float) for seed in sel}
  temp,engine,_,kappas,contexts,source=v.stage_parent(parent,xi[:2])
  report['parent_sha256']=digest(parent);report['native_sha256']=digest(args.native_v1852l)
  report['original_source_provenance']={k:source[k] for k in ('make_context','rhs_from_original','cubic_clock_sign','cubic_phase_excursion')}
  report['fixed_native_kappa']=C['v1852l_fixed_model_kappa']
  results={};summary_rows=[];deriv_rows=[]
  for seed in sel:
   print('AUTHENTIC FIXED-CHARGE RETURN MAP',seed,flush=True)
   results[seed]=analyze_seed(seed,engine,contexts[seed],xi,series[seed],C,args,summary_rows,deriv_rows)
   save('report.json',{**report,'partial_results':results,'scientific_decision':'RUNNING__NO_FINAL_CERTIFICATE'})
  rows('eight_direction_tangent_jacobian.csv',deriv_rows)
  rows('relative_return_summary.csv',summary_rows)
  report['results']=results;report['plots']=render(results)
  report['full_eight_direction_J_tangent_basis_verified']=True
  report['all_selected_native_one_turn_replays_validated']=True
  report['all_two_turn_propagations_completed']=True
  report['relative_fixed_point_found']=any(r['closure_search']['accepted_as_relative_fixed_point'] for r in results.values())
  report['full_8D_fixed_J_floquet_certified']=any(r['closure_search']['floquet_certificate_conditional'] for r in results.values())
  report['scientific_decision']=('RELATIVE_FIXED_POINT_AND_LOCAL_FIXED_J_RETURN_SPECTRUM_CONDITIONALLY_CERTIFIED' if report['full_8D_fixed_J_floquet_certified'] else
                                'RELATIVE_FIXED_POINT_CANDIDATE__FLOQUET_UNRESOLVED' if report['relative_fixed_point_found'] else
                                'AUTHENTIC_FIXED_CHARGE_8D_TANGENT_MAP_AND_TWO_TURN_PROPAGATION__NO_RELATIVE_CYCLE_CERTIFICATE')
  save('report.json',report)
  lines=['# Q550 — v18.52y authenticated server report','','**'+report['scientific_decision']+'**','',
   'Every finite-difference probe is lifted with exactly the same Noether charge; historical kicks are not used as tangent vectors.','',
   'This audits the late *externally forced return model*, not a demonstrated autonomous 11D closed cycle.','',
   '## Numerical diagnostics','']
  for seed,r in results.items():
   lines.extend(['### '+seed,
    '- First-return chart norm: '+str(r['first_return_residual_chart_norm']),
    '- Two-turn chart norm: '+str(r['second_return_residual_chart_norm']),
    '- Native one-turn replay relative error: '+str(r['native_one_turn_relative_error']),
    '- Fixed-J tangent-basis rank: '+str(r['basis']['basis_rank']),
    '- Two-epsilon Jacobian difference: '+str(r['relative_two_step_jacobian_difference']),
    '- Closure search residual: '+str(r['closure_search']['closure_8D_absolute_norm']),
    '- Relative fixed point certified: '+str(r['closure_search']['accepted_as_relative_fixed_point']),
    '- Floquet spectrum conditionally available: '+str(r['closure_search']['floquet_certificate_conditional']),''])
  lines.extend(['## Interpretation','A nonzero return residual excludes a certified relative periodic orbit at that candidate. At non-fixed points DP singular values are local map responses, never Floquet multipliers.','Conserving the counterrotation charge does not establish global energy conservation; transverse attraction, if found, is restricted to the J leaf and the specified driven return protocol.'])
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
 finally:
  if temp is not None:temp.cleanup()

if __name__=='__main__':
 try:main()
 except Exception as exc:
  OUT.mkdir(exist_ok=True,parents=True)
  (OUT/'audit_failure.txt').write_text(traceback.format_exc())
  print('AUTHENTIC_AUDIT_FAILURE',repr(exc),flush=True)
  sys.exit(1)
