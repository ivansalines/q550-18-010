#!/usr/bin/env python3
"""Q550 v18.53h — genuine paired-kappa IVPs and principal-plane residual torque.

Source-driven IVPs begin from authenticated v18.52l native states at xi=12.8,
with identical initial conditions *within each seed family* and identical
historical external phase forcing. They do NOT redo the pre-xi=12.8 ramp.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,sys,traceback
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1];REF=ROOT/'reference';OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'scripts'))
from torque_math import embed,force_components,material_rhs,sigmoid,axes
from bivector_geometry import bivector,integ,cancel_metric,principal_skew,compare_geometry
C=json.loads((ROOT/'inputs_config.json').read_text());P=C['params'];K=C['kappas']

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def write(name,data):
 (OUT/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def rows_csv(name,items):
 if not items:return
 names=list(dict.fromkeys(k for row in items for k in row))
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=names);w.writeheader();w.writerows(items)
def provenance():
 for f,h in C['reference_sha'].items():
  if sha(REF/f)!=h:raise RuntimeError('G_REFERENCE_SHA_MISMATCH '+f)
 previous=json.loads((REF/'v18_53g_report.json').read_text())
 if previous.get('scientific_decision')!='AUTHENTIC_FIXED_STATE_DUAL_KAPPA_PERSISTENCE_MAPPED_WITH_OBSERVER_OR_INVARIANT_FAILURES__NOT_NEW_ORBITS':
  raise RuntimeError('G_AUTHENTIC_STATUS_INVALID')
 full=json.loads((REF/'v18_53g_OUTPUT_MANIFEST.json').read_text())
 for old,new in [('results/report.json','v18_53g_report.json'),('results/persistence_metrics.json','v18_53g_persistence_metrics.json'),
   ('results/coordinate_free_bivector_persistence.csv','v18_53g_coordinate_free_bivector_persistence.csv'),
   ('results/observer_kappa_torque_persistence.csv','v18_53g_observer_kappa_torque_persistence.csv')]:
  if full['files'][old]!=sha(REF/new):raise RuntimeError('G_OUTPUT_MANIFEST_DRIFT '+new)
 if previous.get('native_status')!='SHA256_AUTHENTICATED' or previous.get('parent_status')!='SHA256_AUTHENTICATED':
  raise RuntimeError('G_PARENT_NATIVE_AUTHENTICATION_MISSING')
 return previous

def loader():
 script=REF/'v18_53f_original_run_audit.py'
 if sha(script)!=C['source_sha256']['v18_53f_source']:raise RuntimeError('OLD_ENGINE_LOADER_SHA_MISMATCH')
 spec=importlib.util.spec_from_file_location('q550_v1853h_authenticated_loader',str(script))
 mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
 return mod

def rhs(engine,context,kappa,xi,z):
 th=float(context['theta'](context['event']+float(xi)))
 design=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
 clock=float(engine.clock.cubic_clock_sign(context['event']+float(xi),design))
 val=np.asarray(engine.physical.rhs_from_original(np.r_[z[:10],th],kappa,engine.params,
   engine.slope,clock,engine.support,engine.audit),float)
 if val.shape!=(11,) or not np.all(np.isfinite(val)):raise RuntimeError('BAD_ORIGINAL_11D_RHS')
 if abs(val[10]-engine.omega*clock)>2e-10:raise RuntimeError('EXTERNAL_CLOCK_RHS_DRIFT')
 # Native Y_tau's 11th column is material time; not the external clock.
 tau_rate=float(engine.omega)/float(sigmoid(z[8]))
 return np.r_[val[:10],tau_rate]

def solve(engine,context,kappa,xi,z0,tight=False):
 rt=C['refined_rtol'] if tight else C['integrator_rtol'];at=C['refined_atol'] if tight else C['integrator_atol']
 sol=solve_ivp(lambda t,z:rhs(engine,context,kappa,t,z),(float(xi[0]),float(xi[-1])),z0,
  method='DOP853',t_eval=xi,rtol=rt,atol=at,max_step=C['integrator_max_step_xi'])
 if not sol.success or not np.all(np.isfinite(sol.y)):raise RuntimeError('DUAL_KAPPA_IVP_FAILED: '+str(sol.message))
 return sol.y.T,dict(nfev=int(sol.nfev),steps_evaluated=int(len(xi)),rtol=rt,atol=at,max_step=C['integrator_max_step_xi'])

def validated_rhs(engine,ctx,kappa,xi,z):
 idx=np.linspace(0,len(xi)-1,C['source_probe_count_per_trajectory'],dtype=int)
 errors=[]
 for i in idx:
  phi=ctx['theta'](ctx['event']+float(xi[i]));Fx=rhs(engine,ctx,kappa,xi[i],z[i])[:10]
  fmat=material_rhs(z[i,:10],phi,kappa,P)[0]
  src=Fx*float(sigmoid(z[i,8]))/engine.omega
  err=np.max(abs(src-fmat))/max(1.,np.max(abs(src)))
  errors.append(float(err))
 return max(errors)

def primary_metrics(xi,z,ctx,kappa):
 t=z[:,10];Y=z[:,:10];Q,V,J,ge=embed(Y,P['m_x'])
 phase=np.asarray([ctx['theta'](ctx['event']+float(a)) for a in xi])
 _,_,_,_,_,force=force_components(Y,phase,kappa,P)
 BU=bivector(Q,force['radial_potential']);BC=bivector(Q,force['phase_coupling'])
 U=integ(t,BU);Coup=integ(t,BC);R=U+Coup
 others=sum((bivector(Q,F) for k,F in force.items() if k not in ('radial_potential','phase_coupling')),np.zeros_like(BU))
 All=integ(t,BU+BC+others)
 met=cancel_metric(U,Coup);rec=principal_skew(R);allpr=principal_skew(All)
 met.update(residual_principal_sigma1=rec['sigma1'],residual_principal_sigma2=rec['sigma2'],
   residual_gap=rec['relative_gap'],residual_kernel_component_5=float(abs(rec['kernel'][4])),
   total_six_principal_sigma1=allpr['sigma1'],total_six_principal_sigma2=allpr['sigma2'],
   charge_initial=float(J[0]),charge_range=float(np.ptp(J)),r1_range=[float(np.min(ge['r1'])),float(np.max(ge['r1']))],
   r2_range=[float(np.min(ge['r2'])),float(np.max(ge['r2']))],
   material_duration=float(t[-1]-t[0]))
 return dict(Q=Q,V=V,J=J,U=U,C=Coup,R=R,All=All,B_U=BU,B_C=BC,t=t,metric=met)

def window_metrics(seed,kn,kappa,xi,d):
 rows=[];n=len(xi);intervals=[(0,n-1,'full')]
 for i in range(C['window_count']):
  a=round((n-1)*i/C['window_count']);b=round((n-1)*(i+1)/C['window_count'])
  intervals.append((a,b,'window_'+str(i+1)))
 for i,j,title in intervals:
  a=integ(d['t'][i:j+1],d['B_U'][i:j+1]);b=integ(d['t'][i:j+1],d['B_C'][i:j+1]);R=a+b
  geom=principal_skew(R);m=cancel_metric(a,b)
  rows.append(dict(seed=seed,kappa_rule=kn,kappa=kappa,window=title,xi_start=float(xi[i]),xi_end=float(xi[j]),
    **m,sigma1=geom['sigma1'],sigma2=geom['sigma2'],singular_relative_gap=geom['relative_gap'],
    kernel_abs_x_component=float(abs(geom['kernel'][4])),
    principal_plane_degenerate=bool(geom['relative_gap']<C['principal_gap_warning_relative'])))
 return rows

def native_53g_reconciliation(npz,xi,seed,ctx,omega):
 label=seed+'__reference__Y_tau';arr=np.asarray(npz[label][::C['native_stride']],float)
 y=arr[:,:10];tau=arr[:,10]
 phase=np.asarray([ctx['theta'](ctx['event']+float(v)) for v in xi])
 Q,V,J,g=embed(y,P['m_x']);_,_,_,_,_,F=force_components(y,phase,K['shared'],P)
 U=integ(tau,bivector(Q,F['radial_potential']));Cp=integ(tau,bivector(Q,F['phase_coupling']))
 import csv
 with (REF/'v18_53g_coordinate_free_bivector_persistence.csv').open() as f:
  row=next(r for r in csv.DictReader(f) if r['series']==label and r['kappa_rule']=='shared')
 err=max(abs(np.linalg.norm(U)-float(row['bivector_norm_U'])),abs(np.linalg.norm(Cp)-float(row['bivector_norm_C'])))
 if err>4e-10:raise RuntimeError('G_BIVECTOR_BASELINE_RECONCILIATION_FAIL '+seed+' '+str(err))
 return float(err),clock_error,arr

def analyse(native,engine,contexts):
 from scipy import __version__ as scipy_version
 import numpy as np
 all_windows=[];runs=[];pairs=[];refcheck=[];diverged=[];provenance=[];plotdata={}
 with np.load(native,allow_pickle=False) as f:
  allxi=np.asarray(f['common_xi'],float);xi=allxi[::C['native_stride']]
  if abs(xi[0]-12.8)>1e-10 or abs(xi[-1]-25.6)>1e-10:raise RuntimeError('NATIVE_XI_BOUNDS_CHANGED')
  if len(xi)<128:raise RuntimeError('BAD_SAMPLE_COUNT')
  for seed in ('original','autonomous'):
   ctx=contexts[seed];baseline,clock_error,archived=native_53g_reconciliation(f,xi,seed,ctx,engine.omega)
   refcheck.append(dict(seed=seed,baseline_53g_bivector_error=baseline,material_clock_rate_relative_residual=clock_error))
   initial=np.asarray(archived[0],float)
   data={}
   for kn in ('shared','original','autonomous'):
    kap=K[kn]
    z,info=solve(engine,ctx,kap,xi,initial)
    source_err=validated_rhs(engine,ctx,kap,xi,z)
    if source_err>C['source_rhs_rel_tolerance']:raise RuntimeError('AUTHENTIC_SOURCE_RHS_DRIFT '+seed+' '+kn)
    d=primary_metrics(xi,z,ctx,kap);data[kn]=d
    err_native=float(np.max(np.abs(z-archived)/np.maximum(1,np.abs(archived)))) if kn=='shared' else None
    if kn=='shared' and err_native>C['shared_replay_max_scaled_error']:
     raise RuntimeError('SHARED_KAPPA_NATIVE_REPLAY_FAILURE '+seed+' '+str(err_native))
    refinement_err=None
    if kn in ('original','autonomous'):
     fine,info2=solve(engine,ctx,kap,xi,initial,tight=True)
     refinement_err=float(np.max(abs(fine-z)/np.maximum(1,abs(fine))))
     if refinement_err>C['integrator_refinement_max_scaled_error']:
      raise RuntimeError('DUAL_KAPPA_INTEGRATOR_CONVERGENCE_FAILURE '+seed+' '+kn+' '+str(refinement_err))
    if d['metric']['charge_range']>C['initial_charge_max_drift']:
     raise RuntimeError('DYNAMIC_COUNTERROTATION_CHARGE_DRIFT '+seed+' '+kn)
    runs.append(dict(seed=seed,kappa_rule=kn,kappa=kap,integrator_nfev=info['nfev'],source_rhs_max_relative=source_err,
       native_midpoint_replay_max_scaled_error=err_native,integrator_refinement_error=refinement_err,
       initial_state_is_historical_shared_kappa_seed=True,dynamic_kappa_integrated=True,
       historical_pre_xi_12_8_ramp_reintegrated=False,**d['metric']))
    all_windows.extend(window_metrics(seed,kn,kap,xi,d))
    plotdata[seed,kn]=[xi.copy(),z[:,:10].copy(),d['metric']]
   a,b=data['original'],data['autonomous']
   comparison=compare_geometry(a['R'],b['R'],C['principal_gap_warning_relative'])
   # Pointwise temporal state separation at same xi: true dynamical response.
   ya,yb=plotdata[seed,'original'][1],plotdata[seed,'autonomous'][1]
   Qa,Qb=a['Q'],b['Q']
   qdiff=np.linalg.norm(Qa-Qb,axis=1)
   pairs.append(dict(seed=seed,kappa_original=K['original'],kappa_autonomous=K['autonomous'],
     kappa_gap=K['original']-K['autonomous'],max_5d_position_separation=float(np.max(qdiff)),
     final_5d_position_separation=float(qdiff[-1]),max_native_10d_scaled_separation=float(np.max(abs(ya-yb)/np.maximum(1,abs(ya)))),
     endpoint_5d_separation=float(np.linalg.norm(Qa[-1]-Qb[-1])),
     original_cancellation=a['metric']['invariant_cancellation'],autonomous_cancellation=b['metric']['invariant_cancellation'],
     **comparison))
   for window in ('full','window_1','window_2','window_3','window_4'):
    aw=next(v for v in all_windows if v['seed']==seed and v['kappa_rule']=='original' and v['window']==window)
    bw=next(v for v in all_windows if v['seed']==seed and v['kappa_rule']=='autonomous' and v['window']==window)
    diverged.append(dict(seed=seed,window=window,original_cancellation=aw['invariant_cancellation'],
      autonomous_cancellation=bw['invariant_cancellation'],delta=aw['invariant_cancellation']-bw['invariant_cancellation'],
      same_cancel_class_at_90_percent=(aw['invariant_cancellation']>=.9)==(bw['invariant_cancellation']>=.9)))
 rows_csv('dynamically_integrated_dual_kappa_orbits.csv',runs)
 rows_csv('bivector_principal_window_geometry.csv',all_windows)
 rows_csv('paired_kappa_dynamical_geometry.csv',[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in pairs])
 rows_csv('window_cancellation_persistence.csv',diverged)
 write('paired_kappa_principal_planes.json',pairs)
 write('baseline_53g_reconciliation.json',refcheck)
 return dict(runs=runs,pairs=pairs,window_rows=all_windows,window_persistence=diverged,baseline_check=refcheck),plotdata

def plot(data):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return
 fig,ax=plt.subplots(2,1,figsize=(9.6,7.6),layout='constrained')
 for seed in ('original','autonomous'):
  for kn in ('original','autonomous'):
   xi,Y,met=data[seed,kn]
   r2=np.exp(Y[:,4]);ax[0].plot(xi,r2,label=f'{seed} initial / κ {kn}')
   Q,V,J,g=embed(Y,P['m_x']);ax[1].plot(xi,J-J[0],label=f'{seed} / κ {kn}')
 ax[0].set_ylabel('r2 (self-consistent late-window evolution)');ax[1].set_ylabel('J - J(start)')
 ax[1].set_xlabel('xi');
 for a in ax:a.grid(alpha=.2);a.legend(fontsize=7)
 fig.suptitle('Q550 v18.53h — independently integrated dual-κ late-window orbits')
 fig.savefig(OUT/'dual_kappa_dynamic_radius_and_charge.png',dpi=150);plt.close(fig)

def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--base',type=Path,default=ROOT.parent)
 parser.add_argument('--native-npz',type=Path)
 parser.add_argument('--parent-v1851s',type=Path)
 parser.add_argument('--require-parent',action='store_true')
 parser.add_argument('--require-native',action='store_true')
 parser.add_argument('--require-dynamic',action='store_true')
 args=parser.parse_args();status=dict(version='v18.53h',scientific_decision='UNCLASSIFIED',
  parent_status='MISSING',native_status='MISSING',n_new_Q550_ivps=0,
  no_full_historical_pre_12_8_ramp_claim=True,no_global_bifurcation_or_hinge_noether_certificate=True,
  kappa_rules=K)
 temp=None;code=0
 try:
  provenance();status['v18_53g_reference']='SHA256_AND_PREVIOUS_OUTPUT_MANIFEST_AUTHENTICATED'
  L=loader();parent=L.find_parent(args.parent_v1851s);native=args.native_npz or L.find_npz('10d',args.base)
  if parent and Path(parent).is_file():
   if sha(parent)!=C['parent_v1851s_sha256']:raise RuntimeError('HISTORICAL_51S_PARENT_SHA_MISMATCH')
   status['parent_status']='SHA256_AUTHENTICATED'
  if native and Path(native).is_file():
   if sha(native)!=C['native_10d_sha256']:raise RuntimeError('HISTORICAL_52L_NATIVE_SHA_MISMATCH')
   status['native_status']='SHA256_AUTHENTICATED'
  if status['parent_status']=='SHA256_AUTHENTICATED' and status['native_status']=='SHA256_AUTHENTICATED':
   temp,eng,contexts,meta=L.stage_original(parent)
   status['engine_provenance']=meta
   result,plotdata=analyse(native,eng,contexts)
   status['n_new_Q550_ivps']=10
   status['metrics']=dict(n_traces=len(result['runs']),n_windows=len(result['window_rows']),
    max_charge_drift=max(v['charge_range'] for v in result['runs']),
    max_source_rel=max(v['source_rhs_max_relative'] for v in result['runs']),
    max_replay_error=max(v['native_midpoint_replay_max_scaled_error'] for v in result['runs'] if v['kappa_rule']=='shared'),
    max_refinement_error=max(v['integrator_refinement_error'] for v in result['runs'] if v['integrator_refinement_error'] is not None),
    seed_pair_comparison=result['pairs'],sampled_all_geometric_comparisons=True)
   status['scientific_decision']='AUTHENTIC_SOURCE_SELF_CONSISTENT_DUAL_KAPPA_LATE_WINDOW_DYNAMICAL_GEOMETRY__NOT_FULL_RAMP_OR_NOETHER_HINGE'
   plot(plotdata)
  else:status['scientific_decision']='PARTIAL_AUTHENTIC_G_FROZEN_STATE_REFERENCE__DYNAMIC_IVPS_NOT_RUN'
 except Exception as exc:
  status['scientific_decision']='AUDIT_FAIL_NOT_CERTIFIED';status['failure']=str(exc)
  status['traceback']=traceback.format_exc();code=2
 finally:
  if temp is not None:temp.cleanup()
  write('report.json',status)
  summary=['# Q550 v18.53h — Residual-Torque Principal Planes & Dual-κ Dynamic Persistence','',
   '## Scientific decision', '`'+status['scientific_decision']+'`','',
   '## Actual source availability',f"Parent: {status['parent_status']}; Native: {status['native_status']}",'',
   '## Scientific scope',
   'The frozen-state v18.53g provenance is authenticated, but does not certify distinct trajectories.',
   'On full execution, two independent native 10D IVPs per initial seed are run at κ-original and κ-autonomous, plus a shared-κ baseline replay; a stricter replay checks convergence.',
   'The external historical phase is preserved from the authenticated context; the initial state at xi=12.8 is common between κ variants. This is NOT a complete re-integration of pre-xi=12.8 historical ramp.',
   'Integrated skew bivectors are decomposed into two oriented principal planes. Principal-plane tracking near spectral degeneracy must be interpreted only at the combined-subspace level.',
   'Bivector compensation score uses signed integrals and is not time-local cancellation. No globally conserved hinge charge or new bifurcation certificate is asserted.','',
   '## Metrics',json.dumps(status.get('metrics',{'diagnostic':'not computed'}),indent=2)]
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(summary)+'\n')
  print('V18_53H_SCIENTIFIC_DECISION',status['scientific_decision'])
  if args.require_parent and status['parent_status']!='SHA256_AUTHENTICATED':code=max(code,3)
  if args.require_native and status['native_status']!='SHA256_AUTHENTICATED':code=max(code,3)
  if args.require_dynamic and status['scientific_decision']!='AUTHENTIC_SOURCE_SELF_CONSISTENT_DUAL_KAPPA_LATE_WINDOW_DYNAMICAL_GEOMETRY__NOT_FULL_RAMP_OR_NOETHER_HINGE':code=max(code,3)
 return code
if __name__=='__main__':sys.exit(main())
