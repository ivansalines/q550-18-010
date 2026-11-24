#!/usr/bin/env python3
"""Q550 v18.52v: actual source forward-return drive protocol vs frozen-state work.

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
  kappas={seed:float(z['kappa_'+seed]) for seed in SEEDS}
 if abs((kappas['original']+kappas['autonomous'])/2-C['v1852l_fixed_model_kappa'])>1e-12:raise RuntimeError('KAPPA_GENEALOGY_MISMATCH')
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
                                  'srun_sha256':hashlib.sha256(srun).hexdigest(),'engine_omega':float(engine.omega),'engine_slope':float(engine.slope)}
 phase={};contexts={}
 for seed in SEEDS:
  ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
  contexts[seed]=ctx
  phase[seed]=np.asarray([ctx['theta'](ctx['event']+float(x)) for x in xi],float)
 return temp,engine,phase,kappas,contexts,source

def rows_native(path):
 if digest(path)!=C['v1852l_10d_sha256']:raise RuntimeError('AUTHENTIC_NATIVE_SHA256_MISMATCH')
 with np.load(path,allow_pickle=False) as f:
  xi=np.asarray(f['common_xi'],float)
  series={k[:-len('__Y_tau')]:np.asarray(f[k],float) for k in f.files if k.endswith('__Y_tau') and '__refined__' not in k and '_refined__' not in k}
 if not all(seed+'__reference' in series for seed in SEEDS):raise RuntimeError('NATIVE_REFERENCE_MISSING')
 return xi,series

def source_protocol(engine,contexts):
 # Source law versus derivative: independent, sample forward and return control,
 # without asserting any numerically integrated forward material orbit exists.
 lines=[];spec={};rg=np.linspace(0,1.24,1241)
 for seed in SEEDS:
  ctx=contexts[seed]
  des=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
  f=engine.clock.cubic_clock_sign
  u=engine.clock.cubic_phase_excursion
  sign=np.asarray([float(f(float(t),des)) for t in rg],float)
  ex=np.asarray([float(u(float(t),des)) for t in rg],float)
  dx=np.gradient(ex,rg,edge_order=2)
  err=dx-sign
  # Original source return phase has the hard-coded minus sign.
  for i in range(len(rg)):
   if i%4==0 or (abs(sign[i])<.2 and i%2==0):
    lines.append({'seed':seed,'historical_t':float(rg[i]),'cubic_clock_sign':float(sign[i]),
       'cubic_excursion_derivative':float(dx[i]),'excursion_rate_error':float(err[i]),
       'return_phase_rate_from_source':float(-engine.omega*sign[i]),
       'forward_source_phase_rate':float(engine.omega)})
  spec[seed]={'design':'v18_50p beta:2:2 ramp=0.24','turn_event':float(ctx['event']),
       'clock_sign_min':float(sign.min()),'clock_sign_max':float(sign.max()),
       'excursion_derivative_vs_sign_max_abs':float(np.max(np.abs(err[4:-4]))),
       'phase_rate_source_return':'-engine.omega*cubic_clock_sign',
       'phase_rate_source_forward':'engine.omega'}
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

def signed_ledger(xi,series,phase,kappas,stride):
 out=[];curve=[]
 for name,Y in series.items():
  seed=name.split('_',1)[0];k=kappas[seed]
  take=np.arange(0,len(xi),max(1,stride),dtype=int)
  if take[-1]!=len(xi)-1:take=np.r_[take,len(xi)-1]
  yy=Y[take]; t=yy[:,10];ph=phase[seed][take]
  if not np.all(np.diff(t)>0):raise RuntimeError('NATIVE_NONMONOTONIC_MATERIAL_TIME '+name)
  phase_dot=np.gradient(ph,t,edge_order=2)
  d=energy_terms(yy,k,P,phase=ph)
  phase_coeff=k*d['x']*d['r1']*d['r2']*np.sin(d['D'])
  w_actual=phase_coeff*phase_dot
  w_cf=-w_actual # At the SAME q and Theta, no new trajectory.
  base=-d['D_rho']-d['D_x']+d['W_x']
  net=base+w_actual
  net_cf=base+w_cf
  actual_int=integral(t,w_actual);cf_int=integral(t,w_cf)
  net_int=integral(t,net)
  E=d['E']; ecl=float(np.max(np.abs((E-E[0])-net_int)))
  # Check the frame corrected by a LOCAL instantaneous c(t) while keeping
  # c_dot term; critical at ramp where c varies.
  c=phase_dot/d['x'];dc=np.gradient(c,t,edge_order=2)
  H=E-.5*c*d['x']*d['Jsum']
  hpow=base-.5*(dc*d['x']+c*d['dx'])*d['Jsum']+phase_coeff*(phase_dot-c*d['x'])
  Hclose=float(np.max(np.abs((H-H[0])-integral(t,hpow))))
  charges=d['Jdiff'];jdrift=float(np.max(np.abs(charges-charges[0])))
  tmin=float(t[0]);tmax=float(t[-1]);n=len(t)
  out.append({'series':name,'seed':seed,'n':n,'tau_start':tmin,'tau_end':tmax,
   'drive_actual_work':float(actual_int[-1]),'drive_frozen_state_counterfactual_work':float(cf_int[-1]),
   'drive_reversal_odd_identity':float(abs(actual_int[-1]+cf_int[-1])),
   'dissipation_radial_total':float(np.trapezoid(d['D_rho'],t)),
   'dissipation_x_total':float(np.trapezoid(d['D_x'],t)),
   'x_source_work':float(np.trapezoid(d['W_x'],t)),
   'delta_energy_measured':float(E[-1]-E[0]),
   'delta_energy_predicted':float(net_int[-1]),
   'energy_closure_max_abs':ecl,
   'frame_with_variable_clock_closure_max_abs':Hclose,
   'material_clock_ratio_max_deviation_from_minus_one':float(np.max(np.abs(c+1))),
   'counterrotation_charge_drift':jdrift,
   'interpretation':'actual_return_vs_frozen_state_counterfactual_NOT_forward_trajectory'})
  if name.endswith('__reference'):
   every=max(1,n//400)
   for i in np.r_[np.arange(0,n,every,dtype=int),n-1]:
    curve.append({'seed':seed,'tau':float(t[i]),'external_phase':float(ph[i]),
        'clock_c':float(c[i]),'energy':float(E[i]),'work_drive_actual_cumulative':float(actual_int[i]),
        'work_drive_counterfactual_cumulative':float(cf_int[i]),
        'dissipation_instantaneous':float(d['D_rho'][i]+d['D_x'][i]),
        'work_x_instantaneous':float(d['W_x'][i]),
        'energy_ledger_residual':float((E[i]-E[0])-net_int[i])})
 return out,curve

def plots(rows,curve):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return []
 ans=[]
 for seed in SEEDS:
  a=[z for z in curve if z['seed']==seed]
  if not a:continue
  t=np.array([z['tau'] for z in a]);fig,ax=plt.subplots(3,1,figsize=(10,9),sharex=True,constrained_layout=True)
  ax[0].plot(t,[z['clock_c'] for z in a],label='actual external Theta_dot/x');ax[0].axhline(-1,ls='--',c='gray',label='-1')
  ax[1].plot(t,[z['work_drive_actual_cumulative'] for z in a],label='Archived return drive work')
  ax[1].plot(t,[z['work_drive_counterfactual_cumulative'] for z in a],label='Frozen-state +drive counterfactual',ls='--')
  ax[2].plot(t,[z['energy_ledger_residual'] for z in a],label='Energy ledger error')
  for a0 in ax:a0.grid(alpha=.25);a0.legend(fontsize=8)
  ax[-1].set_xlabel('Material time tau')
  fig.suptitle('Q550 v18.52v — '+seed+' signed-drive protocol and work (counterfactual not forward orbit)')
  out=OUT/(seed+'_signed_drive_work.png');fig.savefig(out,dpi=135);plt.close(fig);ans.append(out.name)
 return ans

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent-v1851s',type=Path);ap.add_argument('--states',type=Path,default=BASE/L52/'results/extended_horizon_10d_states.npz')
 ap.add_argument('--stride',type=int,default=4);ap.add_argument('--require-parent',action='store_true');ap.add_argument('--require-native',action='store_true');args=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 report={'version':'v18.52v','historical_parent_status':'MISSING','native_status':'MISSING','no_heavy_integration':True,
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
  ledger,curve=signed_ledger(xi,series,phase,kappas,args.stride)
  writerows('native_signed_energy_work_ledger.csv',ledger)
  writerows('reference_cumulative_work_curves.csv',curve)
  figs=plots(ledger,curve)
  # No sign selected to improve closure; same authenticated theta used throughout.
  max_close=max(z['energy_closure_max_abs'] for z in ledger)
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
   pclock[seed]={'source_theta_prime_plus_omega_sign_max_abs':float(np.max(np.abs(y+engine.omega*x))),
     'n_probes':int(len(idx)), 'source_clock_sign_min':float(x.min()),'max':float(x.max())}
  report.update({'protocol_source':ps,'source_rate_crosscheck':pclock,
      'max_return_energy_ledger_residual':max_close,'max_counterrotation_charge_drift':max_j,
      'signed_rhs_probes':len(probes),'signed_rhs_max_first10_difference':max(z['counterfactual_RHS_first10_diff_max'] for z in probes),
      'plots':figs,'native_ledger_series':len(ledger),
      'scientific_decision':'AUTHENTIC_SOURCE_PROTOCOL_AND_ARCHIVED_SIGNED_WORK_WITH_FROZEN_STATE_COUNTERFACTUAL_NOT_FORWARD_ORBIT',
      'cannot_claim':'No actual forward branch comparison without independently authenticated forward physical states'})
  dump('report.json',report)
 except Exception as exc:
  report['scientific_decision']='INCOMPLETE_EXCEPTION';report['exception']=repr(exc);dump('report.json',report)
  raise
 finally:
  if temp is not None:temp.cleanup()
 reportlines=['# Q550 v18.52v — Forward–Return Drive Reversal & Signed-Work Attribution Audit','',
  '## Decision',report['scientific_decision'],'',
  '## Exactly what was tested','',
  '1. Historical forward and return source expressions; signed cubic source rate and derivative.',
  '2. Staged source comparison of signed RHS and power channels at matched physical states.',
  '3. On-path archived return-drive work and frozen-state sign-switched work (not a forward orbit).',
  '4. Mechanical energy ledger, variable-clock frame and counter-rotation charge.', '',
  '## Numerical checkpoints','', '```json',json.dumps({k:v for k,v in report.items() if k not in ('no_heavy_integration',)},indent=2),'```','',
  '## Limitations','',
  'The long native v18.52l interval contains RETURN protocol only. A real forward-vs-return work comparison requires authentic forward physical states at phase-matched conditions. '+
  'A frozen-state counterfactual measures a local effect of changing the control, not the evolution of a new trajectory. '+
  'No global noncollapse or topology theorem follows from this audit.']
 (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(reportlines)+'\n')
 print('AUDIT',report['scientific_decision'], 'ENERGY_MAX',report.get('max_return_energy_ledger_residual'))

if __name__=='__main__':main()
