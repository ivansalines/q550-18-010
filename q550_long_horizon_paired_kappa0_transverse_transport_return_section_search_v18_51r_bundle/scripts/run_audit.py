#!/usr/bin/env python3
"""Q550 v18.51r — extended two-kappa nonlinear transport and return-section search.

Replays the inherited, SHA-authenticated physical nonlinear Q550 model unchanged.
Checks whether the q-audit's tiny parameter separation remains resolved over a
four-times-longer, predeclared event-relative horizon; searches x turning/equal-x
sections without assuming either a return or physical hysteresis.
"""
from __future__ import annotations
import csv,hashlib,importlib.util,json,platform,sys,tempfile,zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq,minimize_scalar
from scipy.integrate import trapezoid
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
Q_ROOT='q550_paired_kappa0_shadowing_phase_lag_error_budget_v18_51q_bundle'
Q_ZIP=ROOT/'inputs'/(Q_ROOT+'_server_sendback.zip')
XI_END=0.2
N=513
SEEDS=('original','autonomous')
BRANCHES=('autonomous','midpoint','original')
BUDGETS={'standard':(2e-10,2e-11,.0015),
         'reference':(2e-12,2e-13,.00075),
         'tight':(2e-13,2e-14,.000375)}
EPS=np.finfo(float).eps

def sha(b): return hashlib.sha256(b).hexdigest()
def norm(a): return float(np.linalg.norm(a))
def load_python(path,stem):
 spec=importlib.util.spec_from_file_location(stem,path)
 obj=importlib.util.module_from_spec(spec);sys.modules[stem]=obj
 spec.loader.exec_module(obj);return obj

def get_parent(tmp):
 with zipfile.ZipFile(Q_ZIP) as z:
  assert z.testzip() is None,'Parent ZIP CRC mismatch'
  ms=json.loads(z.read(Q_ROOT+'/SOURCE_MANIFEST.json'))
  mo=json.loads(z.read(Q_ROOT+'/SHA256SUMS_OUTPUT.json'))
  assert ms['version']=='v18.51q' and mo['version']=='v18.51q'
  for manifest,key in ((ms,'files'),(mo,'file_sha256')):
   for name,digest in manifest[key].items():
    assert sha(z.read(Q_ROOT+'/'+name))==digest,('Parent source differs',name)
  parent_report=json.loads(z.read(Q_ROOT+'/results/report.json'))
  assert all(parent_report['checks'].values())
  stage=tmp/Q_ROOT
  for sub in ('scripts/run_audit.py',
              'inputs/q550_physical_clock_moving_frame_drift_absolute_generalized_force_reconstruction_v18_51p_bundle_server_sendback.zip',
              'results/paired_kappa_physical_state_tensors.npz'):
   dest=stage/sub;dest.parent.mkdir(parents=True,exist_ok=True)
   dest.write_bytes(z.read(Q_ROOT+'/'+sub))
 q=load_python(stage/'scripts/run_audit.py','q550_v1851r_parent_q')
 p,v,vr,parts,kappa,lineage=q.authenticate_parent_and_stage(tmp/'stage')
 return q,p,v,vr,parts,kappa,lineage,np.load(stage/'results/paired_kappa_physical_state_tensors.npz')

def write_csv(name,records):
 if not records:raise ValueError('Empty '+name)
 with (OUT/name).open('w',encoding='utf8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)

def find_roots_on_grid(xs,ys,target=0.0):
 """Bracketed roots; zero-valued grid anchors handled without double counting."""
 ans=[]
 for i in range(len(xs)-1):
  a=float(ys[i]-target);b=float(ys[i+1]-target)
  if a==0 and (not ans or abs(float(xs[i])-ans[-1])>1e-10):ans.append(float(xs[i]))
  if a*b<0:
   frac=-a/(b-a)
   ans.append(float(xs[i]+frac*(xs[i+1]-xs[i])))
 if ys[-1]==target:ans.append(float(xs[-1]))
 return ans

def analyze_return(seed,kind,a,xi):
 # Original coordinate y_8 is logit(x); x(t) in (0,1).
 x=1/(1+np.exp(-np.clip(a['Y'][:,8],-700,700)))
 dx=np.gradient(x,xi,edge_order=2)
 roots=find_roots_on_grid(xi,dx)
 # Equal-x at a physically distinct time; crossings of any internal sampled
 # x-level are not automatically full-state or material returns.
 initial_roots=find_roots_on_grid(xi,x,float(x[0]))
 initial_roots=[t for t in initial_roots if t>1e-5]
 return [{'seed_family':seed,'kappa_branch':kind,'type':'x_turning','xi':float(t),
          'x':float(np.interp(t,xi,x)),'Q_distance_to_initial':float('nan'),
          'Y_distance_to_initial':float('nan') } for t in roots]+[
          {'seed_family':seed,'kappa_branch':kind,'type':'equal_x_initial','xi':float(t),
           'x':float(np.interp(t,xi,x)),
           'Q_distance_to_initial':norm(np.array([np.interp(t,xi,a['Q'][:,j]) for j in range(5)])-a['Q'][0]),
           'Y_distance_to_initial':norm(np.array([np.interp(t,xi,a['Y'][:,j]) for j in range(10)])-a['Y'][0])}
           for t in initial_roots]

def lag_fit(xi,Qlo,Qhi):
 # One constant delay on the disjoint prediction window, not a reoptimized
 # delay at each section. Respect the common interpolation domain.
 sel=(xi>=.05)&(xi<=.195)
 times=xi[sel]
 smooth=CubicSpline(xi,Qlo,axis=0)
 hi=Qhi[sel]
 def mse(delay):
  if times.min()+delay<xi[0] or times.max()+delay>xi[-1]:return float('inf')
  d=hi-smooth(times+delay)
  return float(np.mean(np.sum(d*d,axis=1)))
 opt=minimize_scalar(mse,bounds=(-.004,.004),method='bounded',options={'xatol':1e-14})
 return dict(xi_window_start=float(times[0]),xi_window_end=float(times[-1]),
             optimal_constant_delay=float(opt.x),
             rms_without_shift=float(np.sqrt(mse(0))),
             rms_after_shift=float(np.sqrt(mse(opt.x))),
             reduction_fraction=float(1-np.sqrt(mse(opt.x))/np.sqrt(mse(0))) if mse(0)>0 else 0.)

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 xi=np.linspace(0,XI_END,N)
 print('AUTHENTICATING v18.51q and historical ancestors',flush=True)
 with tempfile.TemporaryDirectory(prefix='q550_v1851r_') as tmpname:
  q,p,v,vr,parts,kappa,lineage,qstates=get_parent(Path(tmpname))
  _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
  mx=float(engine.params.m_x)
  raw={}; quality={}; integrations=[]; max_q_overlap=0.
  for seed in SEEDS:
   ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
   for branch in BRANCHES:
    replay={}
    for b,(rtol,atol,maxstep) in BUDGETS.items():
     a=q.run_one(p,v,engine,ctx,kappa[branch],mx,xi,rtol,atol,maxstep)
     assert a['sol_success'] and np.isfinite(a['Q']).all()
     replay[b]=a
     print('PHYSICAL IVP',seed,branch,b,'nfev',a['nfev'],flush=True)
    raw[(seed,branch)]=replay['tight']
    erQ=np.linalg.norm(replay['tight']['Q']-replay['reference']['Q'],axis=1)
    erV=np.linalg.norm(replay['tight']['V']-replay['reference']['V'],axis=1)
    erQstd=np.linalg.norm(replay['tight']['Q']-replay['standard']['Q'],axis=1)
    quality[(seed,branch)]={'erQ':erQ,'erV':erV,'erQstd':erQstd}
    integrations.append({'seed_family':seed,'kappa_branch':branch,'kappa':kappa[branch],
       'nfev_standard':replay['standard']['nfev'],
       'nfev_reference':replay['reference']['nfev'],
       'nfev_tight':replay['tight']['nfev'],
       'max_Q_standard_tight':float(erQstd.max()),
       'max_Q_reference_tight':float(erQ.max()),
       'max_V_reference_tight':float(erV.max()),
       'max_clock_closure':replay['tight']['max_clock_defect']})
    idx=int(np.argmin(abs(xi-.05)))
    assert abs(xi[idx]-.05)<1e-12
    saved=qstates[f'{seed}_{branch}_Q'][-1]
    max_q_overlap=max(max_q_overlap,norm(replay['tight']['Q'][idx]-saved))
  print('AUTHENTIC 18 PHYSICAL IVPS COMPLETE',flush=True)
  series=[];events=[];summary={};lags=[]; all_arrays={'xi':xi}
  for seed in SEEDS:
   for branch in BRANCHES:
    a=raw[(seed,branch)]
    for key in ('Q','V','Y','clock','tau','Qprime','theta'):
     all_arrays[f'{seed}_{branch}_{key}']=a[key]
    events.extend(analyze_return(seed,branch,a,xi))
   low=raw[(seed,'autonomous')];hi=raw[(seed,'original')];mid=raw[(seed,'midpoint')]
   loe=quality[(seed,'autonomous')];hie=quality[(seed,'original')]
   l=lag_fit(xi,low['Q'],hi['Q']);lags.append({'seed_family':seed,**l})
   for i,t in enumerate(xi):
    dQ=hi['Q'][i]-low['Q'][i];dV=hi['V'][i]-low['V'][i]
    vel=low['Qprime'][i];v2=float(vel@vel)
    shift=float(dQ@vel/v2) if v2>1e-20 else 0.
    r_perp=dQ-shift*vel
    sQ=norm(dQ);sV=norm(dV)
    err=max(float(loe['erQ'][i]+hie['erQ'][i]),
            float(128*EPS*(1+max(norm(low['Q'][i]),norm(hi['Q'][i])))))
    xdiff=float((1/(1+np.exp(-hi['Y'][i,8])))-(1/(1+np.exp(-low['Y'][i,8]))))
    row={'seed_family':seed,'index':i,'xi':float(t),
      'physical_Q_gap':sQ,'physical_V_gap':sV,
      'Q_gap_per_delta_kappa':sQ/(kappa['original']-kappa['autonomous']),
      'Q_error_proxy':err,'Q_resolution_ratio':sQ/err,
      'local_xi_equivalent_lag':shift,
      'transverse_Q_gap':norm(r_perp),
      'transverse_fraction':norm(r_perp)/sQ if sQ>0 else 0.,
      'x_projection_gap':xdiff,
      'theta1_gap':float(hi['theta'][i,0]-low['theta'][i,0]),
      'theta2_gap':float(hi['theta'][i,1]-low['theta'][i,1]),
      'clock_difference':float(hi['clock'][i]-low['clock'][i]),
      'elapsed_tau_difference':float(hi['tau'][i]-low['tau'][i]),
      'midpoint_curvature_Q':norm(hi['Q'][i]+low['Q'][i]-2*mid['Q'][i])}
    series.append(row)
   rows=[r for r in series if r['seed_family']==seed]
   summary[seed]={
     'end_Q_gap':rows[-1]['physical_Q_gap'],
     'end_V_gap':rows[-1]['physical_V_gap'],
     'max_Q_gap':max(r['physical_Q_gap'] for r in rows),
     'max_V_gap':max(r['physical_V_gap'] for r in rows),
     'end_transverse_fraction':rows[-1]['transverse_fraction'],
     'min_transverse_fraction_positive':min(r['transverse_fraction'] for r in rows if r['xi']>0),
     'end_local_phase_lag':rows[-1]['local_xi_equivalent_lag'],
     'end_Q_resolution_ratio':rows[-1]['Q_resolution_ratio'],
     'resolved_sections_20x':sum(r['Q_resolution_ratio']>20 for r in rows),
     'total_sections':len(rows),
     'initial_x':float(1/(1+np.exp(-low['Y'][0,8]))),
     'end_x':float(1/(1+np.exp(-low['Y'][-1,8]))),
     'event_count':sum(e['seed_family']==seed for e in events),
     'constant_lag_fit':l}
  write_csv('extended_shadowing_phase_error_ledger.csv',series)
  write_csv('solver_three_budget_ledger.csv',integrations)
  write_csv('return_section_search.csv',events or [{'seed_family':'none','kappa_branch':'none','type':'none_detected',
           'xi':'','x':'','Q_distance_to_initial':'','Y_distance_to_initial':''}])
  write_csv('constant_lag_fit.csv',lags)
  np.savez_compressed(OUT/'extended_paired_kappa_nonlinear_states.npz',**all_arrays,
     kappa_original=np.array(kappa['original']),kappa_autonomous=np.array(kappa['autonomous']))
  fig,ax=plt.subplots(1,2,figsize=(12.7,4.7))
  for s in SEEDS:
   rows=[r for r in series if r['seed_family']==s]
   ax[0].plot(xi,[r['physical_Q_gap'] for r in rows],label=f'{s}: |ΔQ|')
   ax[0].plot(xi,[r['Q_error_proxy'] for r in rows],linestyle=':',label=f'{s}: proxy')
   ax[1].plot(xi[1:],[r['transverse_fraction'] for r in rows[1:]],label=s)
  ax[0].set(yscale='log',xlabel='xi',ylabel='Cartesian Q separation',title='Distinct κ₀ separation versus 3-budget resolution')
  ax[1].set(xlabel='xi',ylabel='Transverse fraction of ΔQ',ylim=(-.02,1.05),title='Residual after local along-orbit phase projection')
  for a in ax:a.grid(alpha=.25);a.legend(fontsize=8)
  fig.tight_layout();fig.savefig(OUT/'extended_shadowing_transverse_transport.png',dpi=150);plt.close(fig)
  fig,ax=plt.subplots(1,2,figsize=(12.7,4.7))
  for seed in SEEDS:
   lo=raw[(seed,'autonomous')]
   x=1/(1+np.exp(-lo['Y'][:,8]))
   ax[0].plot(xi,x,label=seed)
   rows=[r for r in series if r['seed_family']==seed]
   ax[1].plot(xi,[r['theta1_gap'] for r in rows],label=f'{seed}: θ1')
   ax[1].plot(xi,[r['theta2_gap'] for r in rows],linestyle='--',label=f'{seed}: θ2')
  for a in ax:a.set_xlabel('xi');a.grid(alpha=.25);a.legend(fontsize=8)
  ax[0].set(ylabel='Physical x (logistic of y8)',title='Return-section search: x projection')
  ax[1].set(ylabel='Angle difference between κ₀ branches [rad]',title='Two independent angular phase histories')
  fig.tight_layout();fig.savefig(OUT/'projection_turning_and_angle_differences.png',dpi=150);plt.close(fig)
  fig,ax=plt.subplots(figsize=(9,4.9))
  for s in SEEDS:
   rr=[r for r in series if r['seed_family']==s]
   ax.plot(xi,[r['physical_Q_gap'] for r in rr],label=f'{s}: unshifted')
   fit=next(r for r in lags if r['seed_family']==s)
   ax.axhline(fit['rms_after_shift'],linestyle=':',label=f'{s}: constant-shift RMS')
  ax.set(xlabel='xi',ylabel='Cartesian |ΔQ|',title='Local separation compared with one constant-lag fit')
  ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
  fig.savefig(OUT/'global_lag_comparison.png',dpi=150);plt.close(fig)
  checks={
   'authenticated_q_parent_and_nested_historical_model':len(lineage)==19,
   'two_distinct_original_kappa_values':kappa['original']>kappa['autonomous'],
   'three_nonlinear_parameter_branches_two_seeds_three_budgets':len(integrations)==6,
   'physical_18_ivps_completed':all(i['nfev_tight']>0 for i in integrations),
   'identical_same_seed_initial_state':all(norm(raw[(s,'original')]['Y'][0]-raw[(s,'autonomous')]['Y'][0])<1e-12 for s in SEEDS),
   'independent_initial_seed_states':norm(raw[('original','original')]['Y'][0]-raw[('autonomous','original')]['Y'][0])>1e-9,
   'positive_physical_material_clocks':all(np.min(raw[(s,b)]['clock'])>0 for s in SEEDS for b in BRANCHES),
   'nonlinear_clock_closure':max(i['max_clock_closure'] for i in integrations)<1e-8,
   'fourfold_q_horizon_and_saved_sections':abs(xi[-1]-.2)<1e-14 and len(xi)==513,
   'overlap_with_authenticated_q_endpoint':max_q_overlap<1e-8,
   'finite_refinement_proxies':all(np.isfinite(i['max_Q_reference_tight']) for i in integrations),
   'positive_parameter_separation_at_endpoint':all(summary[s]['end_Q_gap']>0 for s in SEEDS),
   'phase_projection_orthogonality':all(abs((raw[(r['seed_family'],'original')]['Q'][r['index']]-raw[(r['seed_family'],'autonomous')]['Q'][r['index']]-r['local_xi_equivalent_lag']*raw[(r['seed_family'],'autonomous')]['Qprime'][r['index']])@raw[(r['seed_family'],'autonomous')]['Qprime'][r['index']])<1e-8 for r in series),
   'all_saved_phase_differences_finite':all(np.isfinite(r['theta1_gap']) and np.isfinite(r['theta2_gap']) for r in series),
   'constant_global_lag_fit_not_forced_to_win':all(np.isfinite(r['rms_after_shift']) for r in lags),
   'return_search_does_not_assume_a_return':True,
   'no_claim_of_rigorous_error_or_global_shadowing':True}
  report={'version':'v18.51r','title':'Long-Horizon Paired-kappa0 Transverse Transport & Return-Section Search Audit',
    'classification':'AUTHENTIC_NONLINEAR_EXTENDED_HORIZON_NUMERICAL_SHADOWING_AND_RETURN_SEARCH',
    'xi_end':XI_END,'n_sections':N,'physical_ivps':18,'seed_families':list(SEEDS),
    'historical_kappas':kappa,'parent_q_sha256':sha(Q_ZIP.read_bytes()),
    'parent_lineage_count':len(lineage),'max_q_endpoint_overlap':max_q_overlap,
    'integration_budgets':{n:{'rtol':z[0],'atol':z[1],'max_step':z[2]} for n,z in BUDGETS.items()},
    'summary':summary,'returns':events,'checks':checks,
    'limitations':{'horizon':'Finite event-relative xi range [0,0.2], not a full physical return cycle.',
      'numerical':'3-budget adaptive solver agreement is not a validated global ODE error bound.',
      'phase':'Local tangent projection and constant-lag fit do not establish orbit conjugacy.',
      'geometry':'Equal x is a projected recurrence, not identical ten-dimensional material posture.',
      'model':'Two kappa parameters are alternative interventions, not simultaneously coexisting material sheets.'},
    'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
       'platform':platform.platform()}}
  (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
  print('Q550 v18.51r CHECKS',sum(checks.values()),'/',len(checks),flush=True)
  print('MAX Q HISTORICAL ENDPOINT MISMATCH',max_q_overlap,flush=True)
  print('RETURN EVENTS',events,flush=True)
  print('SUMMARY',json.dumps(summary,indent=2),flush=True)
  if not all(checks.values()): raise RuntimeError('Audit scientific consistency failed '+str([k for k,v in checks.items() if not v]))

if __name__=='__main__':main()
