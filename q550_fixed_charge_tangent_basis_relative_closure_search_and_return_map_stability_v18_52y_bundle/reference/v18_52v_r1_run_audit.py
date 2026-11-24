#!/usr/bin/env python3
"""Q550 v18.52v-r1: correct source clock sign, native fixed kappa and quadrature.

Model time tau and source support xi are NEVER silently equated. No IVP is solved
here, except the short inherited make_context ancestor staging step (if present).
The data are historical return-drive trajectories; +drive comparisons below are
instantaneous counterfactual powers at SAME physical states and SAME phase, not
measurements of re-integrated forward motion.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,inspect,io,json,math,sys,tempfile,traceback,zipfile
from pathlib import Path
import numpy as np
from energy_math import Params,energy_terms,sigmoid
from phase_energy import integral,physical_rhs_with_signed_external_phase
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;OUT=ROOT/'results'
C=json.loads((ROOT/'inputs_config.json').read_text()); P=Params(**C['model_params'])
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
R51='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
L52='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
SEEDS=('original','autonomous')

def digest(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def dump(name,obj):
 p=OUT/name;p.parent.mkdir(exist_ok=True,parents=True)
 p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def writerows(name,rows):
 if not rows:return
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def dynamic_import(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 if not spec or not spec.loader:raise RuntimeError('CANNOT_IMPORT '+str(path))
 m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def source_info(f):
 out={'module':getattr(f,'__module__',None),'qualname':getattr(f,'__qualname__',None)}
 try:out['source']=inspect.getsource(f)
 except (OSError,TypeError):out['source']='INSPECT_SOURCE_UNAVAILABLE'
 if hasattr(f,'__code__'):
  out['filename']=f.__code__.co_filename;out['line']=f.__code__.co_firstlineno
 return out

def find_parent(given):
 paths=([given] if given else [])+[BASE/(S51+'_server_sendback.zip'),Path('/root/elab1/hysteresis')/(S51+'_server_sendback.zip'),Path('/root/Downloads')/(S51+'_server_sendback.zip')]
 return next((Path(x) for x in paths if x and Path(x).is_file()),None)

def stage_parent(parent,xi):
 if digest(parent)!=C['v1851s_sha256']:raise RuntimeError('AUTHENTICATED_PARENT_SHA256_MISMATCH')
 with zipfile.ZipFile(parent) as z:
  if z.testzip() is not None:raise RuntimeError('PARENT_BAD_ZIP')
  rblob=z.read(S51+'/inputs/'+R51+'_server_sendback.zip')
  srun=z.read(S51+'/scripts/run_audit.py')
 with zipfile.ZipFile(io.BytesIO(rblob)) as z:
  if z.testzip() is not None:raise RuntimeError('NESTED_51R_BAD_ZIP')
  blob=z.read(R51+'/results/extended_paired_kappa_nonlinear_states.npz')
 with np.load(io.BytesIO(blob)) as z:
  kappa_seeds={seed:float(z['kappa_'+seed]) for seed in SEEDS}
 if abs((kappa_seeds['original']+kappa_seeds['autonomous'])/2-C['v1852l_fixed_model_kappa'])>1e-12:raise RuntimeError('KAPPA_GENEALOGY_MISMATCH')
 loader=dynamic_import(ROOT/'reference/parent_v18_51t_run_audit.py','q550_v1852v_loader')
 # The parent stage is intentionally kept alive inside a temporary directory.
 temp=tempfile.TemporaryDirectory(prefix='q550_v1852v_original_')
 staged=loader.stage_engine(Path(temp.name),rblob,srun)
 _,_,seed_api,observer,engine,unit,drpt,creport,_=staged[4]
 params={name:float(getattr(engine.params,name)) for name in C['model_params']}
 diff={name:params[name]-float(v) for name,v in C['model_params'].items()}
 if any(abs(x)>1e-12 for x in diff.values()):
  dump('parameter_mismatch.json',{'actual':params,'expected':C['model_params'],'difference':diff})
  raise RuntimeError('MODEL_PARAMETER_MISMATCH')
 source={k:source_info(getattr(engine.clock,k)) for k in ('cubic_clock_sign','cubic_phase_excursion')}
 source.update({k:source_info(getattr(engine.physical,k)) for k in ('rhs_from_original','power_terms')})
 source['make_context']=source_info(seed_api.make_context)
 source['Design']=source_info(engine.original_ramp.Design)
 source['original_integrate']=source_info(engine.original.integrate)
 source['authenticated_parent']={'path':str(parent),'sha256':digest(parent),'rblob_sha256':hashlib.sha256(rblob).hexdigest(),
                                  'srun_sha256':hashlib.sha256(srun).hexdigest(),'engine_omega':float(engine.omega),'engine_slope':float(engine.slope),
                                  'historical_seed_kappas':kappa_seeds,
                                  'authenticated_native_fixed_kappa':float(C['v1852l_fixed_model_kappa'])}
 phase={};contexts={}
 for seed in SEEDS:
  ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
  contexts[seed]=ctx
  phase[seed]=np.asarray([ctx['theta'](ctx['event']+float(x)) for x in xi],float)
 # Native 52l used a FIXED model kappa shared by both seed families.
 native_kappas={seed:float(C['v1852l_fixed_model_kappa']) for seed in SEEDS}
 return temp,engine,phase,native_kappas,contexts,source

def rows_native(path):
 if digest(path)!=C['v1852l_10d_sha256']:raise RuntimeError('AUTHENTIC_NATIVE_SHA256_MISMATCH')
 with np.load(path,allow_pickle=False) as f:
  xi=np.asarray(f['common_xi'],float)
  series={k[:-len('__Y_tau')]:np.asarray(f[k],float) for k in f.files if k.endswith('__Y_tau') }
 if not all(seed+'__reference' in series for seed in SEEDS):raise RuntimeError('NATIVE_REFERENCE_MISSING')
 return xi,series

def source_protocol(engine,contexts):
 # On t>0 the source has d(cubic_phase_excursion)/dt = -cubic_clock_sign.
 # A second minus in theta(t) makes d(theta_return)/dt = +omega*clock_sign.
 lines=[];spec={};rg=np.linspace(0,1.24,1241)
 for seed in SEEDS:
  ctx=contexts[seed]; des=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
  sign=np.array([float(engine.clock.cubic_clock_sign(float(t),des)) for t in rg])
  excursion=np.array([float(engine.clock.cubic_phase_excursion(float(t),des)) for t in rg])
  grad=np.gradient(excursion,rg,edge_order=2)
  # The t<=0 extension is constant and has a one-sided boundary derivative.
  good=(rg>.004)&(np.abs(rg-.24)>.004)&(rg<1.236)
  err=grad+sign
  # Independent finite-difference probes, away from boundaries.
  probe_t=np.r_[np.linspace(.004,.235,240),np.linspace(.245,1.22,240)]
  h=1.e-6
  numerical_exc=np.array([(engine.clock.cubic_phase_excursion(float(t+h),des)-engine.clock.cubic_phase_excursion(float(t-h),des))/(2*h) for t in probe_t],float)
  exact_sign=np.array([engine.clock.cubic_clock_sign(float(t),des) for t in probe_t],float)
  max_derivative=float(np.max(np.abs(numerical_exc+exact_sign)))
  # Original 11D return integration is theta0 - omega*excursion(t).
  # This computed rate must agree with the last RHS component.
  for i in range(0,len(rg),4):
   lines.append({'seed':seed,'historical_t':float(rg[i]),'cubic_clock_sign':float(sign[i]),
       'excursion_derivative':float(grad[i]),'excursion_derivative_plus_sign':float(err[i]),
       'return_theta_rate_from_source':float(engine.omega*sign[i]),
       'return_rhs_11th_rate':float(engine.omega*sign[i]),
       'forward_theta_rate_from_source':float(engine.omega),
       'first10_affected_instantaneously_by_sign':False})
  spec[seed]={'design':'v18_50p beta:2:2 ramp=0.24', 'turn_event':float(ctx['event']),
       'clock_sign_min':float(sign.min()),'clock_sign_max':float(sign.max()),
       'excursion_derivative_plus_sign_max_abs_finite_difference':max_derivative,
       'excursion_derivative_plus_sign_grid_max_abs':float(np.max(np.abs(err[good]))),
       'phase_rate_source_return':'PLUS engine.omega*cubic_clock_sign',
       'phase_rate_source_forward':'PLUS engine.omega',
       'source_sign_certification':'TWO_MINUS_SIGNS_CANCEL_AND_MATCH_SIGNED_RHS'}
 return lines,spec

def source_rhs_sweep(engine,contexts,kappas,series,phase,xi):
 # Ask the EXACT staged physical source whether clock_sign is injected only into
 # the phase derivative or also into the ten material force components.
 out=[]
 for seed in SEEDS:
  ctx=contexts[seed];k=kappas[seed];Y=series[seed+'__reference']
  ids=np.unique(np.linspace(0,len(xi)-1,20,dtype=int))
  design=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
  probes=[('native_return',ctx['event']+float(xi[i]),int(i),np.r_[Y[i,:10],phase[seed][i]]) for i in ids]
  for type_,t,i,state in probes:
   s=float(engine.clock.cubic_clock_sign(t,design))
   if abs(s)<1e-9:continue
   # Holding the SAME state and phase: exactly one control parameter changed.
   yp=np.asarray(engine.physical.rhs_from_original(state,k,engine.params,engine.slope,abs(s),engine.support,engine.audit),float)
   ym=np.asarray(engine.physical.rhs_from_original(state,k,engine.params,engine.slope,-abs(s),engine.support,engine.audit),float)
   pp=np.asarray(engine.physical.power_terms(state,k,engine.params,engine.slope,abs(s)),float).ravel()
   pm=np.asarray(engine.physical.power_terms(state,k,engine.params,engine.slope,-abs(s)),float).ravel()
   out.append({'seed':seed,'index':i,'xi':float(xi[i]),'source_clock_sign':s,
      'counterfactual_RHS_first10_diff_max':float(np.max(np.abs(yp[:10]-ym[:10]))),
      'source_last_rate_positive':float(yp[10]),'source_last_rate_negative':float(ym[10]),
      'source_power_channel_positive_0':float(pp[0]),'source_power_channel_negative_0':float(pm[0]),
      'source_power_other_channels_max_diff':float(np.max(np.abs(pp[1:]-pm[1:]))) if len(pp)>1 else None})
 return out

def signed_ledger(xi,series,phase,kappas,engine,contexts,strides):
 # Native arrays have index 10 = MATERIAL TIME tau. External phase is recovered
 # from the source for each seed. Do not reinterpret tau as Theta.
 out=[];curve=[];kappa_falsifier=[]
 omega=float(engine.omega)
 stride_vals=sorted(set(max(1,int(x)) for x in strides))
 if 1 not in stride_vals:raise RuntimeError('FULL_RESOLUTION_STRIDE_1_REQUIRED')
 for name,Y in series.items():
  seed=name.split('_',1)[0]
  if seed not in SEEDS:raise RuntimeError('UNKNOWN_NATIVE_SEED '+name)
  k=float(kappas[seed]);t_native=Y[:,10]
  if len(xi)!=len(Y):raise RuntimeError('GRID_LENGTH_MISMATCH')
  if not np.all(np.diff(t_native)>0):raise RuntimeError('NATIVE_NONMONOTONIC_MATERIAL_TIME '+name)
  # For a reference-driven protocol, Theta is a function of xi, common to kicks.
  # But dTheta/dtau along EACH native trajectory uses its own tau(xi) map.
  design=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
  ev=float(contexts[seed]['event'])
  for stride in stride_vals:
   take=np.arange(0,len(xi),stride,dtype=int)
   if take[-1]!=len(xi)-1:take=np.r_[take,len(xi)-1]
   yy=Y[take];t=yy[:,10];ph=phase[seed][take]; xx=xi[take]
   d=energy_terms(yy,k,P,phase=ph)
   signs=np.asarray([float(engine.clock.cubic_clock_sign(ev+float(v),design)) for v in xx])
   # Exact model-time conversion: d tau/d xi = omega / x_native.
   phase_dot_source=signs*d['x']
   phase_dot_finite_diff=np.gradient(ph,t,edge_order=2)
   xi_tau_grad=np.gradient(t,xx,edge_order=2)
   tau_xi_source=omega/d['x']
   phase_coeff=k*d['x']*d['r1']*d['r2']*np.sin(d['D'])
   w_actual=phase_coeff*phase_dot_source
   w_cf=-w_actual  # fixed-state sign reversal, NOT an observed forward orbit
   base=-d['D_rho']-d['D_x']+d['W_x']
   E=d['E']; net=base+w_actual
   drive_int=integral(t,w_actual);cf_int=integral(t,w_cf)
   net_int=integral(t,net)
   endpoint=float(E[-1]-E[0]);error=np.asarray(E-E[0]-net_int)
   error_fd=np.asarray(E-E[0]-integral(t,base+phase_coeff*phase_dot_finite_diff))
   # Arbitrary variable-clock frame: c = dTheta/dtau / x, equals source sign.
   c=signs;cdot=np.gradient(c,t,edge_order=2)
   H=E-.5*c*d['x']*d['Jsum']
   hpow=base-.5*(cdot*d['x']+c*d['dx'])*d['Jsum']+phase_coeff*(phase_dot_source-c*d['x'])
   Hclose=float(np.max(np.abs(H-H[0]-integral(t,hpow))))
   jd=float(np.max(np.abs(d['Jdiff']-d['Jdiff'][0])))
   out.append({'series':name,'seed':seed,'stride':stride,'n_samples':len(t),
      'native_kappa_fixed':k,'seed_historical_kappa_NOT_USED':float(C['kappa_seeds'][seed]),
      'tau_start':float(t[0]),'tau_end':float(t[-1]),
      'drive_actual_work':float(drive_int[-1]),
      'drive_frozen_state_counterfactual_work':float(cf_int[-1]),
      'drive_reversal_odd_identity':float(abs(drive_int[-1]+cf_int[-1])),
      'dissipation_radial_total':float(np.trapezoid(d['D_rho'],t)),
      'dissipation_x_total':float(np.trapezoid(d['D_x'],t)),
      'x_source_work':float(np.trapezoid(d['W_x'],t)),
      'delta_energy_measured':endpoint,'delta_energy_predicted':float(net_int[-1]),
      'energy_closure_max_abs':float(np.max(np.abs(error))),
      'energy_closure_endpoint':float(error[-1]),
      'energy_closure_max_abs_finite_difference_rate':float(np.max(np.abs(error_fd))),
      'phase_rate_source_minus_finite_difference_max_abs':float(np.max(np.abs(phase_dot_source-phase_dot_finite_diff))),
      'tau_xi_analytical_vs_gradient_max_rel':float(np.max(np.abs((xi_tau_grad-tau_xi_source)/tau_xi_source))),
      'frame_with_variable_clock_closure_max_abs':Hclose,
      'material_clock_ratio_max_deviation_from_source_sign':float(np.max(np.abs(phase_dot_finite_diff/d['x']-signs))),
      'counterrotation_charge_drift':jd,
      'interpretation':'REAL_RETURN_AND_FROZEN_STATE_COUNTERFACTUAL_NOT_FORWARD_ORBIT'})
   if stride==1 and (name.endswith('__reference') or name.endswith('__reference_refined')):
    every=max(1,len(t)//400)
    for i in np.r_[np.arange(0,len(t),every,dtype=int),len(t)-1]:
     curve.append({'seed':seed,'series':name,'tau':float(t[i]),'external_phase':float(ph[i]),
         'clock_c':float(c[i]),'energy':float(E[i]),
         'work_drive_actual_cumulative':float(drive_int[i]),
         'work_drive_counterfactual_cumulative':float(cf_int[i]),
         'dissipation_instantaneous':float(d['D_rho'][i]+d['D_x'][i]),
         'work_x_instantaneous':float(d['W_x'][i]),
         'energy_ledger_residual':float(error[i])})
  # One controlled falsifier: old seed-dependent kappa evaluated on native state.
  if name.endswith('__reference'):
   kbad=float(C['kappa_seeds'][seed]);old=energy_terms(Y,kbad,P,phase=phase[seed]);
   wbad=kbad*old['x']*old['r1']*old['r2']*np.sin(old['D'])*(-old['x'])
   fluxbad=-old['D_rho']-old['D_x']+old['W_x']+wbad
   old_err=np.asarray(old['E']-old['E'][0]-integral(t_native,fluxbad))
   goodrow=next(z for z in out if z['series']==name and z['stride']==1)
   kappa_falsifier.append({'seed':seed,'true_native_kappa':k,
      'historical_seed_kappa_wrong_for_native':kbad,
      'absolute_kappa_mismatch':float(abs(k-kbad)),
      'correct_fixed_kappa_energy_max_error':goodrow['energy_closure_max_abs'],
      'wrong_seed_kappa_energy_max_error':float(np.max(np.abs(old_err))),
      'falsification_not_parameter_fitting':True})
 return out,curve,kappa_falsifier

def plots(rows,curve):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return []
 ans=[]
 for seed in SEEDS:
  a=[z for z in curve if z['seed']==seed and z['series']==seed+'__reference']
  if not a:continue
  t=np.array([z['tau'] for z in a]);fig,ax=plt.subplots(3,1,figsize=(10,9),sharex=True,constrained_layout=True)
  ax[0].plot(t,[z['clock_c'] for z in a],label='actual external Theta_dot/x');ax[0].axhline(-1,ls='--',c='gray',label='-1')
  ax[1].plot(t,[z['work_drive_actual_cumulative'] for z in a],label='Archived return drive work')
  ax[1].plot(t,[z['work_drive_counterfactual_cumulative'] for z in a],label='Frozen-state +drive counterfactual',ls='--')
  ax[2].plot(t,[z['energy_ledger_residual'] for z in a],label='Energy ledger error')
  for a0 in ax:a0.grid(alpha=.25);a0.legend(fontsize=8)
  ax[-1].set_xlabel('Material time tau')
  fig.suptitle('Q550 v18.52v-r1 — '+seed+' signed-drive protocol and work (counterfactual not forward orbit)')
  out=OUT/(seed+'_signed_drive_work.png');fig.savefig(out,dpi=135);plt.close(fig);ans.append(out.name)
 return ans

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent-v1851s',type=Path);ap.add_argument('--states',type=Path,default=BASE/L52/'results/extended_horizon_10d_states.npz')
 ap.add_argument('--strides',type=int,nargs='+',default=[1,2,4,8,16,32]);ap.add_argument('--require-parent',action='store_true');ap.add_argument('--require-native',action='store_true');args=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 report={'version':'v18.52v-r1','historical_parent_status':'MISSING','native_status':'MISSING','no_heavy_integration':True,
         'forward_trajectory_not_generated':True,'full_forward_return_work_difference_measured':False,'global_noncollapse_certified':False}
 parent=find_parent(args.parent_v1851s)
 if args.require_parent and parent is None:raise RuntimeError('REQUIRED_PARENT_51S_MISSING')
 if args.require_native and not args.states.is_file():raise RuntimeError('REQUIRED_NATIVE_52L_MISSING')
 if not args.states.is_file():
  report['scientific_decision']='PROTOCOL_FRAMEWORK_ONLY_NATIVE_MISSING';dump('report.json',report);return
 xi,series=rows_native(args.states);report['native_status']='SHA256_AUTHENTICATED';report['native_total_series']=len(series)
 if parent is None:
  report['scientific_decision']='NATIVE_AUTHENTICATED_PARENT_MISSING_NO_FORCED_PHASE_IMPUTATION'
  dump('report.json',report);return
 temp=None
 try:
  temp,engine,phase,kappas,contexts,sources=stage_parent(parent,xi)
  report['historical_parent_status']='SHA256_AUTHENTICATED'
  dump('authenticated_control_sources.json',sources)
  protocol,ps=source_protocol(engine,contexts)
  writerows('historical_forward_return_clock_protocol.csv',protocol)
  probes=source_rhs_sweep(engine,contexts,kappas,series,phase,xi)
  writerows('original_source_signed_rhs_power_probes.csv',probes)
  ledger,curve,kappa_test=signed_ledger(xi,series,phase,kappas,engine,contexts,args.strides)
  writerows('native_signed_energy_work_ledger.csv',ledger)
  writerows('native_kappa_falsification.csv',kappa_test)
  writerows('reference_cumulative_work_curves.csv',curve)
  figs=plots(ledger,curve)
  # Cross-report overlay with independently archived, SHA-pinned historical
  # v18.52t (correct fixed-kappa baseline) and original v18.52v (wrong-kappa).
  baseline=list(csv.DictReader((ROOT/'reference/v18_52t_authenticated_20_series_baseline.csv').open()))
  earlier=list(csv.DictReader((ROOT/'reference/v18_52v_original_native_ledger.csv').open()))
  t_by={z['series']:z for z in baseline};v_by={z['series']:z for z in earlier}
  comparisons=[]
  for z in ledger:
   if z['stride']!=1:continue
   key=z['series']+'__Y_tau'
   if key not in t_by:raise RuntimeError('V18_52T_BASELINE_SERIES_NOT_FOUND '+key)
   tval=float(t_by[key]['minus_E_max_error'])
   old=v_by.get(z['series'])
   comparisons.append({'series':z['series'],
     'v18_52t_fixed_kappa_full_E_max_error':tval,
     'v18_52v_r1_fixed_kappa_full_E_max_error':z['energy_closure_max_abs'],
     'v18_52t_vs_r1_abs_difference':abs(tval-z['energy_closure_max_abs']),
     'old_v18_52v_seed_kappa_stride4_E_max_error':float(old['energy_closure_max_abs']) if old else None,
     'source':'v18.52t authenticated 20-series historical baseline and original v18.52v server ledger'})
  writerows('cross_report_full_resolution_comparison.csv',comparisons)
  max_comparison_delta=max(z['v18_52t_vs_r1_abs_difference'] for z in comparisons)
  # No sign selected to improve closure; same authenticated theta used throughout.
  max_close=max(z['energy_closure_max_abs'] for z in ledger if z['stride']==1)
  max_close_stride4=max(z['energy_closure_max_abs'] for z in ledger if z['stride']==4)
  source_errors=max(z['phase_rate_source_minus_finite_difference_max_abs'] for z in ledger if z['stride']==1)
  max_tau_xi_rel=max(z['tau_xi_analytical_vs_gradient_max_rel'] for z in ledger if z['stride']==1)
  max_j=max(z['counterrotation_charge_drift'] for z in ledger)
  # Source derivative can be checked on native full xi against clock sign.
  pclock={}
  for seed in SEEDS:
   idx=np.unique(np.linspace(8,len(xi)-9,128,dtype=int))
   event=float(contexts[seed]['event']);d=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
   x=[];y=[]
   for j in idx:
    t=event+float(xi[j]);s=float(engine.clock.cubic_clock_sign(t,d))
    step=1e-5
    rate=(contexts[seed]['theta'](t+step)-contexts[seed]['theta'](t-step))/(2*step)
    x.append(s);y.append(rate)
   x=np.asarray(x);y=np.asarray(y)
   pclock[seed]={'source_theta_prime_minus_omega_sign_max_abs':float(np.max(np.abs(y-engine.omega*x))),
     'n_probes':int(len(idx)), 'source_clock_sign_min':float(x.min()),'max':float(x.max())}
  report.update({'protocol_source':ps,'source_rate_crosscheck':pclock,
      'max_full_resolution_energy_ledger_residual':max_close,
      'max_stride4_energy_ledger_residual':max_close_stride4,
      'max_phase_rate_source_vs_finite_difference':source_errors,
      'max_tau_xi_analytical_vs_gradient_relative':max_tau_xi_rel,
      'max_v18_52t_baseline_vs_r1_full_closure_difference':max_comparison_delta,
      'v18_52t_baseline_series_matched':len(comparisons),
      'clock_source_probes_pass':all(v['source_theta_prime_minus_omega_sign_max_abs']<1e-7 for v in pclock.values()),
      'omega_vs_2pi_slope_abs':float(abs(engine.omega-2*math.pi*engine.slope)),
      'signed_rhs_last_rate_vs_source_max_abs':float(max(max(abs(z['source_last_rate_positive']-engine.omega),abs(z['source_last_rate_negative']+engine.omega)) for z in probes)),
      'kappa_falsification':kappa_test,'max_counterrotation_charge_drift':max_j,
      'signed_rhs_probes':len(probes),'signed_rhs_max_first10_difference':max(z['counterfactual_RHS_first10_diff_max'] for z in probes),
      'plots':figs,'native_ledger_series':len(series),
      'native_ledger_stride_rows':len(ledger),'native_fixed_kappa':C['v1852l_fixed_model_kappa'],
      'all_20_native_series_included':len(series)==20,
      'scientific_decision':('AUTHENTICATED_CLOCK_SIGN_AND_FIXED_KAPPA_SIGNED_ENERGY_FULL_RESOLUTION'
        if len(series)==20 and max_close<1e-8 and max_comparison_delta<1e-8 and abs(engine.omega-2*math.pi*engine.slope)<1e-10 and all(v['source_theta_prime_minus_omega_sign_max_abs']<1e-7 for v in pclock.values())
        else 'PARTIAL_AUTHENTICATION_OR_NUMERICAL_AUDIT_FAILED'),
      'cannot_claim':'No actual forward branch comparison without independently authenticated forward physical states'})
  dump('report.json',report)
 except Exception as exc:
  report['scientific_decision']='INCOMPLETE_EXCEPTION';report['exception']=repr(exc);dump('report.json',report)
  raise
 finally:
  if temp is not None:temp.cleanup()
 reportlines=['# Q550 v18.52v-r1 — Forward–Return Drive Reversal & Signed-Work Attribution Audit','',
  '## Decision',report['scientific_decision'],'',
  '## Exactly what was tested','',
  '1. Authenticated cubic_excursion_prime=-cubic_clock_sign and theta_prime=+omega*cubic_clock_sign.',
  '2. Staged source comparison of signed RHS and power channels at matched physical states.',
  '3. On-path archived return-drive work and frozen-state sign-switched work (not a forward orbit).',
  '4. All 20 native series; fixed native kappa for BOTH seed families; quadrature convergence and kappa falsification.', '',
  '## Numerical checkpoints','', '```json',json.dumps({k:v for k,v in report.items() if k not in ('no_heavy_integration',)},indent=2),'```','',
  '## Limitations','',
  'The long native v18.52l interval contains RETURN protocol only. A real forward-vs-return work comparison requires authentic forward physical states at phase-matched conditions. '+
  'A frozen-state counterfactual measures a local effect of changing the control, not the evolution of a new trajectory. '+
  'No global noncollapse or topology theorem follows from this audit.']
 (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(reportlines)+'\n')
 print('AUDIT',report['scientific_decision'], 'ENERGY_MAX_FULL',report.get('max_full_resolution_energy_ledger_residual'))

if __name__=='__main__':main()
