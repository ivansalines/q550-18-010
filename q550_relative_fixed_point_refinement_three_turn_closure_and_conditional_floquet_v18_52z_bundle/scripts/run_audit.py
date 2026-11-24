#!/usr/bin/env python3
"""Q550 v18.52z: source-authenticated fixed-J Newton refinement and conditional Floquet.

All authentic 10D integrations belong on the user's server, NOT in local preflight.
No root, period or stability conclusion is asserted without independent
recomputation and integrity checking of the historical source and native states.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, math, sys, traceback
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, simpson
from return_map import charge, reduce8, lift8, gauge, chart_diff, fixed_charge_fd, LABELS, TWOPI
from energy_math import energy_terms, Params

ROOT=Path(__file__).resolve().parents[1]; BASE=ROOT.parent; OUT=ROOT/'results'
CFG=json.loads((ROOT/'inputs_config.json').read_text())
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
L52='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
SEEDS=('original','autonomous')

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 return h.hexdigest()

def report_save(report):
 OUT.mkdir(parents=True,exist_ok=True)
 (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')

def module(path, name):
 spec=importlib.util.spec_from_file_location(name,path)
 m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def cplx(a):return {'real':float(np.real(a)),'imag':float(np.imag(a)),'modulus':float(abs(a))}

def desc(u):return {k:float(v) for k,v in zip(LABELS,u)}

def normdelta(a,b):return float(np.linalg.norm(chart_diff(a,b)))

def strict_native_replay(engine,ctx,xi,Y,J,g,args,old):
 Y0=Y[0,:10];u0=reduce8(Y0)
 ya=module(ROOT/'reference/v18_52y_run_audit.py','q550_v52y_authenticated_map')
 t0=float(xi[0]);period,phase=ya.phase_period(engine,ctx,t0)
 # Existing historical native gate: candidate must be based on an authentic one-turn map.
 obj=ya.AuthenticatedMap(engine,ctx,CFG['v1852l_fixed_model_kappa'],t0,period,args.rtol,args.atol,args.maxstep)
 obj.set_leaf(J,g)
 native_u=ya.reduce8(ya.native_y_at(xi,Y,t0+period))
 replay=obj.map(u0)
 mismatch=normdelta(replay,native_u)/max(1.,float(np.linalg.norm(native_u)))
 if mismatch>2e-6:raise RuntimeError('NATIVE_ONE_TURN_REPLAY_MISMATCH '+repr(mismatch))
 check=float(old['native_one_turn_relative_error'])
 if abs(mismatch-check)>2e-6:raise RuntimeError('V1852Y_REPLAY_REGRESSION_MISMATCH')
 return ya,period,phase,{'native_replay_error':mismatch,'previous_v1852y_replay_error':check}

def newton_refine(obj,start,args):
 x=np.asarray(start,float).copy(); scales=np.maximum(1.,np.abs(x)); hist=[]; best_x=x.copy();best_rr=float('inf')
 for it in range(args.max_newton+1):
  p=obj.map(x);d=chart_diff(x,p); rr=float(np.linalg.norm(d)); rnorm=float(np.linalg.norm(d/scales))
  rec={'iteration':it,'residual_absolute':rr,'residual_scaled':rnorm,'chart':x.tolist()}
  hist.append(rec)
  if rr<best_rr:best_rr=rr;best_x=x.copy()
  print(f'  Newton {it:02d} |P(u)-u|={rr:.9e} scaled={rnorm:.9e}',flush=True)
  if rr<args.closure_tol:
   rec['stop_reason']='CLOSURE_THRESHOLD_REACHED';break
  if it>=args.max_newton:
   rec['stop_reason']='NEWTON_ITERATION_LIMIT';break
  D,_=fixed_charge_fd(obj.map,x,args.newton_fd)
  A=D-np.eye(8);s=np.linalg.svd(A,compute_uv=False)
  rec['min_singular_DP_minus_I']=float(s[-1]); rec['max_singular_DP_minus_I']=float(s[0]);rec['jacobian_cond']=float(s[0]/max(s[-1],1e-300))
  try: step=np.linalg.lstsq(A,-d,rcond=1e-11)[0]
  except np.linalg.LinAlgError:
   rec['stop_reason']='SINGULAR_NEWTON_OPERATOR';break
  ls=float(np.linalg.norm(step/scales))
  if ls>args.max_weighted_step:step*=args.max_weighted_step/ls
  accepted=False
  for b in range(args.max_backtracks+1):
   factor=2.**(-b); trial=x+factor*step
   try:
    dd=chart_diff(trial,obj.map(trial)); rt=float(np.linalg.norm(dd))
   except (ValueError,RuntimeError,FloatingPointError,OverflowError):continue
   if math.isfinite(rt) and rt<rr*(1-1e-4*factor):
    x=trial;rec['accepted_backtrack']=b;rec['accepted_residual_absolute']=rt;accepted=True;break
  if not accepted:
   rec['stop_reason']='BACKTRACK_NO_DECREASE';break
 return best_x,best_rr,hist

def three_turn_gate(obj,u,omega,args):
 one=obj.map(u,1);two=obj.map(u,2);three=obj.map(u,3)
 R=[normdelta(u,p) for p in (one,two,three)]
 sem2=normdelta(obj.map(one),two)
 sem3=normdelta(obj.map(two),three)
 extra3=normdelta(obj.map(obj.map(one)),three)
 C={'one_turn_residual_absolute':R[0], 'two_turn_residual_absolute':R[1],
    'three_turn_residual_absolute':R[2], 'P_squared_vs_direct_two_turn':sem2,
    'P_cubed_vs_direct_three_turn':extra3, 'P_of_direct_two_turn_vs_direct_three_turn':sem3,
    'one_turn_residual_relative':R[0]/max(1.,np.linalg.norm(u)),
    'two_turn_residual_relative':R[1]/max(1.,np.linalg.norm(u)),
    'three_turn_residual_relative':R[2]/max(1.,np.linalg.norm(u)),
    'charge_J':obj.j, 'max_integration_charge_drift':obj.max_j_drift,
    'phase_three_turns':-3*TWOPI}
 C['pass']=bool(R[0]<args.closure_tol and R[1]<args.two_turn_tol and R[2]<args.three_turn_tol
            and sem2<args.composition_tol and sem3<args.composition_tol and extra3<args.composition_tol
            and obj.max_j_drift<args.charge_tol)
 return C

def independent_gate(engine,ctx,u,J,g,t0,period,args,coarse,ya):
 obj=ya.AuthenticatedMap(engine,ctx,CFG['v1852l_fixed_model_kappa'],t0,period,
                         args.verify_rtol,args.verify_atol,args.verify_maxstep,cache=False)
 obj.set_leaf(J,g)
 fine=three_turn_gate(obj,u,engine.omega,args)
 fine['independent_solve_one_turn_difference']=normdelta(coarse.map(u,1),obj.map(u,1))
 fine['independent_solve_three_turn_difference']=normdelta(coarse.map(u,3),obj.map(u,3))
 shifted=ya.AuthenticatedMap(engine,ctx,CFG['v1852l_fixed_model_kappa'],t0,period,args.verify_rtol,args.verify_atol,args.verify_maxstep,cache=False)
 shifted.set_leaf(J,g+0.341)
 fine['candidate_counterrotation_gauge_equivariance_error']=normdelta(obj.map(u,1),shifted.map(u,1))
 fine['pass']=bool(fine['pass'] and fine['independent_solve_one_turn_difference']<args.solver_agree_tol
                  and fine['independent_solve_three_turn_difference']<args.solver_agree_tol
                  and fine['candidate_counterrotation_gauge_equivariance_error']<args.gauge_tol)
 return obj,fine

def spectrum(obj,ind,u,args):
 D={}
 for eps in args.jac_steps:
  mat,_=fixed_charge_fd(obj.map,u,eps)
  D[str(eps)]=mat
 def comp(a,b):return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(b)))
 e1,e2,e3=args.jac_steps
 convergence12=comp(D[str(e1)],D[str(e2)])
 convergence23=comp(D[str(e2)],D[str(e3)])
 Din,_=fixed_charge_fd(ind.map,u,e2)
 diff_independent=comp(D[str(e2)],Din)
 eig=np.linalg.eigvals(D[str(e2)])
 radius=float(np.max(np.abs(eig)))
 svals=np.linalg.svd(D[str(e2)],compute_uv=False)
 fdpass=convergence12<args.jac_tol and convergence23<args.jac_tol and diff_independent<args.jac_tol
 return {'finite_difference_steps':args.jac_steps,'jacobian_disagreement_large_mid':convergence12,
         'jacobian_disagreement_mid_small':convergence23,
         'independent_solver_jacobian_disagreement':diff_independent,
         'eigenvalues':[cplx(x) for x in eig], 'spectral_radius':radius,
         'singular_values_descriptive_coordinate_dependent':svals.tolist(),
         'jacobian_numerically_converged':fdpass,
         'conditional_leaf_attractor':(radius<1-args.attraction_margin) if fdpass else None,
         'attraction_decision_margin':args.attraction_margin,
         'no_global_conservative_claim':True}

def energy_cycle_ledger(obj,engine,u,args):
 from energy_math import energy_terms
 y0=lift8(u,obj.j,obj.g)
 p=Params(**CFG['model_params']);k=obj.k
 sol=solve_ivp(obj.source,(obj.t0,obj.t0+obj.period),y0,method='DOP853',rtol=args.verify_rtol,atol=args.verify_atol,max_step=args.verify_maxstep,dense_output=True)
 if not sol.success:raise RuntimeError('CYCLE_ENERGY_IVP_FAILURE')
 recs=[]
 for n in CFG['energy_ledger_quadrature_n']:
  xi=np.linspace(obj.t0,obj.t0+obj.period,int(n))
  y=np.asarray(sol.sol(xi),float).T
  phases=np.array([float(obj.ctx['theta'](float(obj.ctx['event'])+float(x))) for x in xi])
  e=energy_terms(y,k,p,phase=phases)
  xs=e['x'];dtau_dxi=float(engine.omega)/xs
  # On the authenticated late-return branch dTheta/dtau = -x, while dTheta/dxi = -omega.
  rate=-xs
  Wtheta=k*xs*e['r1']*e['r2']*np.sin(e['D'])*rate
  wphase=float(simpson(Wtheta*dtau_dxi,x=xi));wx=float(simpson(e['W_x']*dtau_dxi,x=xi))
  d_rad=float(simpson(e['D_rho']*dtau_dxi,x=xi));d_x=float(simpson(e['D_x']*dtau_dxi,x=xi))
  Dt=float(simpson(dtau_dxi,x=xi))
  deltaE=float(e['E'][-1]-e['E'][0]);balance=wphase+wx-d_rad-d_x
  recs.append({'grid_points':int(n),'delta_energy':deltaE, 'work_drive':wphase,'work_x':wx,
               'dissipation_rho':d_rad,'dissipation_x':d_x, 'net_signed_work_minus_dissipation':balance,
               'energy_balance_residual':deltaE-balance, 'material_time_per_turn':Dt,
               'energy_initial':float(e['E'][0]),'energy_final':float(e['E'][-1])})
 return {'quadrature_replicates':recs,
         'ledger_two_grid_difference':abs(recs[0]['energy_balance_residual']-recs[1]['energy_balance_residual']),
         'not_global_energy_conservation':True,
         'no_historical_forward_branch_claim':True}

def analyze_seed(seed,engine,ctx,xi,Y,old,args):
 Y0=Y[0,:10];J=charge(Y0);g=gauge(Y0);u0=reduce8(Y0)
 assert abs(J-float(old['J_fixed']))<5e-12,'V1852Y_CHARGE_ANCHOR_MISMATCH'
 cand=np.asarray(old['closure_search']['candidate_u'],float)
 assert cand.shape==(8,) and np.all(np.isfinite(cand)),'INVALID_PARENT_CANDIDATE'
 ya,T,phase,replay=strict_native_replay(engine,ctx,xi,Y,J,g,args,old)
 t0=float(xi[0]);reference_dist=normdelta(u0,cand)
 # 3T timing is certified against the source phase, not assumed from a clock change.
 phase0=float(ctx['theta'](float(ctx['event'])+t0))
 drift=[float(ctx['theta'](float(ctx['event'])+t0+i*T))-phase0 for i in (1,2,3)]
 if any(abs(drift[i]+(i+1)*TWOPI)>1e-7 for i in range(3)):raise RuntimeError('THREE_TURN_PHASE_PERIOD_FAILURE')
 obj=ya.AuthenticatedMap(engine,ctx,CFG['v1852l_fixed_model_kappa'],t0,T,args.rtol,args.atol,args.maxstep)
 obj.set_leaf(J,g)
 start_res=normdelta(cand,obj.map(cand))
 print(f'{seed}: authentic y candidate residual {start_res:.9e}',flush=True)
 if abs(start_res-float(old['closure_search']['closure_8D_absolute_norm']))>1e-6:
  raise RuntimeError('PARENT_CANDIDATE_REPLAY_REGRESSION_FAILURE '+seed)
 xf,best,history=newton_refine(obj,cand,args)
 closure=three_turn_gate(obj,xf,engine.omega,args)
 independent,ver=independent_gate(engine,ctx,xf,J,g,t0,T,args,obj,ya)
 accept=bool(closure['pass'] and ver['pass'])
 print(f'{seed}: 3-turn gate accepted={accept} abs1={closure["one_turn_residual_absolute"]:.3e} abs3={closure["three_turn_residual_absolute"]:.3e}',flush=True)
 result={'seed':seed,'J_reference':J,'fixed_kappa':obj.k,'late_source_period_xi':T,
         'authenticated_source_phase_three_turns':drift,'native_replay':replay,
         'v1852y_candidate_report_sha256':CFG['v1852y_report_sha256'],
         'initial_candidate':cand.tolist(),'initial_candidate_residual_absolute':start_res,
         'initial_historical_chart_distance':reference_dist,
         'newton_history':history,'newton_final_candidate':xf.tolist(),
         'newton_smallest_residual':best, 'one_two_three_turns':closure,
         'independent_tight_solver_validation':ver,
         'relative_fixed_point_accepted':accept,
         'conditional_floquet':None, 'energy_ledger':None,
         'solver_effort':{'coarse_ivps':obj.eval_count,'tight_ivps':independent.eval_count,
                          'coarse_nfev':obj.nfev,'tight_nfev':independent.nfev}}
 if accept:
  sp=spectrum(obj,independent,xf,args)
  result['conditional_floquet']=sp
  result['energy_ledger']=energy_cycle_ledger(independent,engine,xf,args)
  # Full certificate requires 8D root AND FD/solver-converged spectral data.
  result['conditional_floquet_certified']=bool(sp['jacobian_numerically_converged'])
 else:
  result['conditional_floquet_certified']=False
 return result

def plot(results):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return []
 fns=[]
 for seed,r in results.items():
  history=r['newton_history'];xx=[h['iteration'] for h in history]; yy=[h['residual_absolute'] for h in history]
  fig,axes=plt.subplots(2,1,figsize=(9,7),constrained_layout=True)
  axes[0].semilogy(xx,np.maximum(yy,1e-16),marker='o',label='Fixed-J relative return residual')
  axes[0].axhline(CFG['closure_chart_abs_tol'],ls='--',label='Acceptance threshold')
  axes[0].set_xlabel('Newton iteration');axes[0].set_ylabel('Chart 8D residual');axes[0].legend();axes[0].grid(alpha=.2)
  d=r['one_two_three_turns']
  axes[1].bar(['1 turn','2 turns','3 turns'],[d['one_turn_residual_absolute'],d['two_turn_residual_absolute'],d['three_turn_residual_absolute']])
  axes[1].set_yscale('log');axes[1].set_ylabel('Direct reclosure residual');axes[1].grid(alpha=.2)
  fig.suptitle('Q550 v18.52z — '+seed+' | fixed-J closure and three-turn audit')
  fn=seed+'_three_turn_fixed_charge_closure.png';fig.savefig(OUT/fn,dpi=145);plt.close(fig);fns.append(fn)
 return fns

def main():
 ap=argparse.ArgumentParser(description='Authenticated late-return fixed-J Newton and three-turn audit')
 ap.add_argument('--parent-v1851s',type=Path);ap.add_argument('--native-v1852l',type=Path,default=BASE/L52/'results/extended_horizon_10d_states.npz')
 ap.add_argument('--require-parent',action='store_true');ap.add_argument('--require-native',action='store_true');ap.add_argument('--seeds',choices=['both','original','autonomous'],default='both')
 ap.add_argument('--max-newton',type=int,default=CFG['max_newton_steps']);ap.add_argument('--rtol',type=float,default=CFG['solver_rtol']);ap.add_argument('--atol',type=float,default=CFG['solver_atol']);ap.add_argument('--maxstep',type=float,default=CFG['solver_max_step_xi'])
 ap.add_argument('--verify-rtol',type=float,default=CFG['independent_solver_rtol']);ap.add_argument('--verify-atol',type=float,default=CFG['independent_solver_atol']);ap.add_argument('--verify-maxstep',type=float,default=CFG['independent_max_step_xi'])
 ap.add_argument('--closure-tol',type=float,default=CFG['closure_chart_abs_tol']);ap.add_argument('--two-turn-tol',type=float,default=CFG['two_turn_chart_abs_tol']);ap.add_argument('--three-turn-tol',type=float,default=CFG['three_turn_chart_abs_tol']);ap.add_argument('--composition-tol',type=float,default=CFG['composed_vs_direct_three_turn_tol']);ap.add_argument('--solver-agree-tol',type=float,default=CFG['independent_solver_one_turn_tol']);ap.add_argument('--charge-tol',type=float,default=CFG['maximum_charge_drift'])
 ap.add_argument('--gauge-tol',type=float,default=CFG['gauge_equivariance_tol']);ap.add_argument('--newton-fd',type=float,default=CFG['newton_fd_epsilon']);ap.add_argument('--jac-steps',type=float,nargs=3,default=CFG['jacobian_fd_epsilon']);ap.add_argument('--jac-tol',type=float,default=CFG['jacobian_rel_step_convergence_tol']);ap.add_argument('--attraction-margin',type=float,default=CFG['leaf_attraction_radius_margin']);ap.add_argument('--max-backtracks',type=int,default=CFG['max_backtracks']);ap.add_argument('--max-weighted-step',type=float,default=CFG['max_weighted_step'])
 args=ap.parse_args();OUT.mkdir(exist_ok=True,parents=True)
 if not (args.rtol>0 and args.verify_rtol>0 and args.closure_tol>0 and args.max_newton>=0 and all(x>0 for x in args.jac_steps)):
  raise ValueError('INVALID_AUDIT_TOLERANCES')
 path=ROOT/'reference/v18_52y_authenticated_report.json'
 if sha(path)!=CFG['v1852y_report_sha256']:raise RuntimeError('V1852Y_CANDIDATE_SOURCE_SHA256_MISMATCH')
 yreport=json.loads(path.read_text())
 if yreport['scientific_decision']!='AUTHENTIC_FIXED_CHARGE_8D_TANGENT_MAP_AND_TWO_TURN_PROPAGATION__NO_RELATIVE_CYCLE_CERTIFICATE':
  raise RuntimeError('UNEXPECTED_V1852Y_PARENT_STATUS')
 candidates=([args.parent_v1851s] if args.parent_v1851s else [])+[BASE/(S51+'_server_sendback.zip'),Path('/root/elab1/hysteresis')/(S51+'_server_sendback.zip'),Path('/root/Downloads')/(S51+'_server_sendback.zip')]
 parent=next((p for p in candidates if p and p.is_file()),None)
 report={'version':'v18.52z','title':CFG['title'],'v1852y_candidate_sha256':sha(path),
         'v1851s_available':bool(parent),'v1852l_available':args.native_v1852l.is_file(),
         'source_authorization_required':True,'global_conservative_system_certified':False,
         'full_historical_forward_return_cycle_certified':False,'no_unclosed_floquet_claim':True,
         'seeds':args.seeds,'configuration':vars(args).copy()}
 report['configuration']={k:str(v) if isinstance(v,Path) else v for k,v in report['configuration'].items()}
 if not parent or not args.native_v1852l.is_file():
  report['scientific_decision']='AUTHENTIC_INPUTS_UNAVAILABLE__NO_Q550_INTEGRATION'
  report_save(report)
  if (args.require_parent and not parent) or (args.require_native and not args.native_v1852l.is_file()):
   raise RuntimeError('REQUIRED_AUTHENTIC_ANCESTOR_NOT_FOUND')
  return
 if sha(parent)!=CFG['v1851s_sha256']:raise RuntimeError('V1851S_SOURCE_SHA256_MISMATCH')
 if sha(args.native_v1852l)!=CFG['v1852l_10d_sha256']:raise RuntimeError('V1852L_NATIVE_SHA256_MISMATCH')
 stage=module(ROOT/'reference/v18_52v_r1_run_audit.py','q550_v52z_v_parent')
 temp=None
 try:
  with np.load(args.native_v1852l,allow_pickle=False) as z:
   xi=np.asarray(z['common_xi'],float)
   if len(xi)!=65537 or abs(xi[0]-12.8)>1e-8 or abs(xi[-1]-25.6)>1e-8:raise RuntimeError('UNEXPECTED_NATIVE_SECTION_GRID')
   selected=SEEDS if args.seeds=='both' else (args.seeds,)
   Y={seed:np.asarray(z[seed+'__reference__Y_tau'],float) for seed in selected}
  temp,engine,_,_,contexts,source=stage.stage_parent(parent,xi[:2])
  if abs(float(engine.omega)-float(yreport['results']['original']['phase_period']['source_omega']))>1e-10:raise RuntimeError('AUTHENTIC_SOURCE_OMEGA_CHANGED')
  report['v1851s_sha256']=sha(parent);report['v1852l_sha256']=sha(args.native_v1852l)
  report['authenticated_source_provenance']={k:source[k] for k in ('make_context','rhs_from_original','cubic_clock_sign','cubic_phase_excursion')}
  report['fixed_kappa']=CFG['v1852l_fixed_model_kappa'];report['results']={}
  for seed in selected:
   print('SOURCE-AUTHENTICATED THREE TURN FIXED-J AUDIT',seed,flush=True)
   report['results'][seed]=analyze_seed(seed,engine,contexts[seed],xi,Y[seed],yreport['results'][seed],args)
   report['scientific_decision']='PARTIAL_AUTHENTIC_RESULTS__NOT_FINAL_CERTIFICATE'
   report_save(report)
  report['plots']=plot(report['results'])
  report['relative_fixed_points_certified']={s:r['relative_fixed_point_accepted'] for s,r in report['results'].items()}
  report['conditional_floquet_certified']={s:r['conditional_floquet_certified'] for s,r in report['results'].items()}
  passed=list(report['relative_fixed_points_certified'].values())
  floquet=list(report['conditional_floquet_certified'].values())
  if all(passed) and all(floquet):status='TWO_AUTHENTIC_RELATIVE_FIXED_POINTS_THREE_TURN_AND_CONDITIONAL_FLOQUET_CERTIFIED'
  elif any(floquet):status='PARTIAL_RELATIVE_FIXED_POINT_AND_CONDITIONAL_FLOQUET_CERTIFICATE'
  elif any(passed):status='RELATIVE_FIXED_POINT_REACHED__FULL_FLOQUET_UNCERTIFIED'
  else:status='SOURCE_AUTHENTIC_REFINEMENT__NO_RELATIVE_CYCLE_CERTIFICATE'
  report['scientific_decision']=status
  report_save(report)
  lines=['# Q550 v18.52z — authenticated three-turn fixed-J return audit','', '**'+status+'**','',
   'The source late-return protocol is externally driven. This certificate cannot be extended to the full historical forward/return 11D dynamics without further verification.','',
   'All root searches are restricted to a J-fixed eight-coordinate quotient; energy is not a globally conserved quantity.','']
  for seed,r in report['results'].items():
   lines+=['## '+seed, '',f'- Reference charge J: {r["J_reference"]:.16g}',
           f'- Parent v18.52y candidate |P(u)-u|: {r["initial_candidate_residual_absolute"]:.12g}',
           f'- Refined |P(u)-u|: {r["one_two_three_turns"]["one_turn_residual_absolute"]:.12g}',
           f'- Direct 2-turn reclosure: {r["one_two_three_turns"]["two_turn_residual_absolute"]:.12g}',
           f'- Direct 3-turn reclosure: {r["one_two_three_turns"]["three_turn_residual_absolute"]:.12g}',
           f'- Independent tight-solver validation: {r["independent_tight_solver_validation"]["pass"]}',
           f'- Relative fixed point accepted: {r["relative_fixed_point_accepted"]}',
           f'- Floquet conditionally certified: {r["conditional_floquet_certified"]}','']
   if r['conditional_floquet']:
    lines+=['- Conditional quotient spectral radius: '+str(r['conditional_floquet']['spectral_radius']),
            '- Within-fixed-J attracting: '+str(r['conditional_floquet']['conditional_leaf_attractor']),'']
   if r['energy_ledger']:
    last=r['energy_ledger']['quadrature_replicates'][-1]
    lines+=['- Energy change on the candidate relative cycle: '+str(last['delta_energy']),
            '- Drive work: '+str(last['work_drive']), '- x work: '+str(last['work_x']),
            '- Radial dissipation: '+str(last['dissipation_rho']), '- x dissipation: '+str(last['dissipation_x']),
            '- Energy ledger residual: '+str(last['energy_balance_residual']),'']
  lines+=['## Strict scientific boundaries','Even a converged root is a relative cycle of the source-authenticated late, constant-negative external-drive protocol. Floquet eigenvalues are conditional on the exact charge leaf and quotient. A neutral direction or absent spectral gap must not be called attracting. A charge conservation law alone is not total-energy conservation.']
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
 finally:
  if temp is not None:temp.cleanup()

if __name__=='__main__':
 try:main()
 except Exception as exc:
  OUT.mkdir(parents=True,exist_ok=True)
  (OUT/'audit_failure.txt').write_text(traceback.format_exc())
  print('AUTHENTIC_V1852Z_AUDIT_FAILURE',repr(exc),flush=True)
  sys.exit(1)
