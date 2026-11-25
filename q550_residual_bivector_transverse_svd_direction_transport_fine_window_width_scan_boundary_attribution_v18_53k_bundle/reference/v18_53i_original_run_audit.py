#!/usr/bin/env python3
"""Q550 v18.53i — parameter sensitivity, principal-plane transport, late-horizon IVPs.

Authentic v18.51s historical generator, v18.52l native conditions, and the
v18.53h fix1 code/outputs are pinned by SHA256. Beyond xi=25.6 the forcing
is a *model continuation* of the historical external clock, not measured history.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';REF=ROOT/'reference';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'scripts'))
from bivector_geometry import principal_skew,compare_geometry
from torque_math import embed
CONFIG=json.loads((ROOT/'inputs_config.json').read_text());C=CONFIG['long_horizon'];K=CONFIG['kappas']

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def write_json(p,x):
 (OUT/p).write_text(json.dumps(x,indent=2,allow_nan=False,ensure_ascii=False)+'\n')

def csv_rows(p,rows):
 if not rows:return
 names=list(dict.fromkeys(k for r in rows for k in r))
 with (OUT/p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=names);w.writeheader();w.writerows(rows)

def import_history():
 p=REF/'v18_53h_fix1_run_audit.py'
 sp=importlib.util.spec_from_file_location('q550_53h_authenticated_ancestor',p)
 m=importlib.util.module_from_spec(sp);sys.modules[sp.name]=m;sp.loader.exec_module(m)
 return m

def assert_lineage(h):
 h.provenance()
 for n,expected in CONFIG['v1853h_references_sha256'].items():
  if sha(REF/n)!=expected:raise RuntimeError('H_BASELINE_REFERENCE_SHA256_FAIL '+n)
 manifest=json.loads((REF/'v18_53h_OUTPUT_MANIFEST.json').read_text())
 pairs={
  'v18_53h_authenticated_report.json':'results/report.json',
  'v18_53h_dynamically_integrated_dual_kappa_orbits.csv':'results/dynamically_integrated_dual_kappa_orbits.csv',
  'v18_53h_paired_kappa_principal_planes.json':'results/paired_kappa_principal_planes.json',
  'v18_53h_bivector_principal_window_geometry.csv':'results/bivector_principal_window_geometry.csv'}
 for loc,old in pairs.items():
  if manifest['files'][old]!=sha(REF/loc):raise RuntimeError('H_OUTPUT_MANIFEST_FAIL '+loc)
 old=json.loads((REF/'v18_53h_authenticated_report.json').read_text())
 if old.get('scientific_decision')!='AUTHENTIC_SOURCE_SELF_CONSISTENT_DUAL_KAPPA_LATE_WINDOW_DYNAMICAL_GEOMETRY__NOT_FULL_RAMP_OR_NOETHER_HINGE' or old.get('n_new_Q550_ivps')!=10:
  raise RuntimeError('H_DYNAMIC_ANCESTOR_NOT_CERTIFIED')
 if old.get('parent_status')!='SHA256_AUTHENTICATED' or old.get('native_status')!='SHA256_AUTHENTICATED':
  raise RuntimeError('H_PARENT_OR_NATIVE_NOT_CERTIFIED')
 if K!=old.get('kappa_rules'):raise RuntimeError('H_KAPPA_RULE_GENEALOGY_DRIFT')
 return old

def grid_kappa():
 delta=K['original']-K['autonomous'];v=[('autonomous',K['autonomous']),('original',K['original']),('shared',K['shared'])]
 for mag in (0.25,1,3,10):
  v.extend([(f'center_minus_{mag:g}_gap',K['shared']-mag*delta),(f'center_plus_{mag:g}_gap',K['shared']+mag*delta)])
 if len(set(k for _,k in v))!=len(v):raise RuntimeError('DUPLICATE_KAPPA_GRID')
 return v

def xi_grid(native_xi,xi_max):
 old=np.asarray(native_xi,float)[::CONFIG['native_stride']]
 if len(old)<100 or abs(old[0]-C['start_xi'])>1e-9 or abs(old[-1]-C['historical_endpoint_xi'])>1e-9:
  raise RuntimeError('ARCHIVED_NATIVE_WINDOW_CHANGED')
 dur=float(xi_max-C['historical_endpoint_xi'])
 if dur<=0:return old,len(old)
 n=int(round(dur/C['extension_grid_step_xi']))
 extra=np.linspace(C['historical_endpoint_xi'],xi_max,n+1)[1:]
 xx=np.r_[old,extra]
 if np.any(np.diff(xx)<=0):raise RuntimeError('NONINCREASING_XI_GRID')
 return xx,len(old)

def checkpoints(xi_max):
 return [v for v in C['checkpoints_xi'] if v<=xi_max+1e-9]

def slices_for(xi,xi_max):
 pts=[float(C['start_xi'])]+checkpoints(xi_max)
 if abs(pts[-1]-xi_max)>1e-7:pts.append(xi_max)
 records=[]
 for i,(a,b) in enumerate(zip(pts[:-1],pts[1:])):
  records.append((f'segment_{i+1}',float(a),float(b)))
 for i,b in enumerate(pts[1:]):records.append((f'cumulative_{i+1}',pts[0],float(b)))
 return records

def extract_segment(xi,z,lo,hi):
 i=int(np.argmin(abs(xi-lo)));j=int(np.argmin(abs(xi-hi)))
 if abs(xi[i]-lo)>1e-7 or abs(xi[j]-hi)>1e-7 or j-i<50:raise RuntimeError('MISSING_CHECKPOINT_GRID')
 return xi[i:j+1],z[i:j+1]

def independent_53h_reconcile(se,kn,metric):
 with (REF/'v18_53h_dynamically_integrated_dual_kappa_orbits.csv').open() as f:
  hist=next((r for r in csv.DictReader(f) if r['seed']==se and r['kappa_rule']==kn),None)
 if hist is None:raise RuntimeError('MISSING_53H_BASELINE_ROW')
 err_s=abs(metric['residual_principal_sigma1']-float(hist['residual_principal_sigma1']))
 err_c=abs(metric['invariant_cancellation']-float(hist['invariant_cancellation']))
 if err_s>C['baseline_53h_sigmas_abs_tolerance'] or err_c>C['baseline_53h_cancellation_abs_tolerance']:
  raise RuntimeError(f'REPLAY_53H_RESIDUAL_MISMATCH {se} {kn} sigma={err_s} cancellation={err_c}')
 return dict(baseline_sigma1_abs_error=float(err_s),baseline_cancellation_abs_error=float(err_c))

def tensor_geometry(d):
 r=d['R'];p=principal_skew(r)
 return dict(residual_norm=float(np.linalg.norm(r)),radial_norm=float(np.linalg.norm(d['U'])),
  phase_norm=float(np.linalg.norm(d['C'])),cancellation=float(d['metric']['invariant_cancellation']),
  residual_s1=p['sigma1'],residual_s2=p['sigma2'],gap=p['relative_gap'])

def numerical_quadrature_defect(d):
 from bivector_geometry import integ
 rr=d['B_U']+d['B_C'];fine=d['R'];half=integ(d['t'][::2],rr[::2])
 # All our segments have even sample intervals; guard a missing endpoint explicitly.
 if len(d['t'])%2==0:return None
 return float(np.linalg.norm(fine-half)/max(1.,np.linalg.norm(fine)))

def plane_projector(d,i):
 p=principal_skew(d['R']);block=p['blocks'][i]
 return block['projector'] if block['projector'] is not None else None

def pair_metrics(a,b):
 g=compare_geometry(a['R'],b['R'],C['finite_difference_projector_gap_warn'])
 angles=[max(v['principal_angles_deg']) if v['principal_angles_deg'] is not None else None for v in g['planes']]
 return dict(bivector_cosine=float(g['full_bivector_cosine']),relative_change=float(g['relative_residual_change']),
  plane1_angle_deg=angles[0],plane2_angle_deg=angles[1],
  plane_pairing_warning=bool(g['principal_plane_pairing_ill_conditioned']))

def derivative_metrics(lo,mid,hi,step):
 A,B,D=lo['R'],mid['R'],hi['R']
 v=(D-A)/(2*step)
 curvature=(D-2*B+A)/(step*step)
 m=dict(dR_dk_fro=float(np.linalg.norm(v)),dR_dk_relative=float(np.linalg.norm(v)/max(np.linalg.norm(B),1e-15)),
  curvature_R_fro=float(np.linalg.norm(curvature)),nonlinearity_ratio=float(np.linalg.norm(D-2*B+A)/max(np.linalg.norm(D-A),1e-15)),
  d_sigma1_dk=float((hi['metric']['residual_principal_sigma1']-lo['metric']['residual_principal_sigma1'])/(2*step)),
  d_sigma2_dk=float((hi['metric']['residual_principal_sigma2']-lo['metric']['residual_principal_sigma2'])/(2*step)),
  d_cancellation_dk=float((hi['metric']['invariant_cancellation']-lo['metric']['invariant_cancellation'])/(2*step)))
 for i in (0,1):
  Pm=plane_projector(mid,i);Pl=plane_projector(lo,i);Ph=plane_projector(hi,i)
  if all(p is not None for p in (Pm,Pl,Ph)) and min(principal_skew(d['R'])['relative_gap'] for d in (lo,mid,hi))>=C['finite_difference_projector_gap_warn']:
   m[f'dP{i+1}_dk_fro']=float(np.linalg.norm((Ph-Pl)/(2*step)))
   m[f'plane{i+1}_sensitivity_angle_full_deg']=float(max(pair_metrics(lo,hi)[f'plane{i+1}_angle_deg'],0))
  else:m[f'dP{i+1}_dk_fro']=None;m[f'plane{i+1}_sensitivity_angle_full_deg']=None
 return m

def plot_results(separations,sensitivities,transport):
 try:
  import matplotlib;matplotlib.use('Agg')
  import matplotlib.pyplot as plt
 except ImportError:return
 # Separate plots; every point is an independently integrated IVP, not an interpolation.
 for seed in ('original','autonomous'):
  rows=[r for r in separations if r['seed']==seed and r['window'].startswith('cumulative_')]
  if not rows:continue
  fig,ax=plt.subplots(2,1,figsize=(9,8),layout='constrained')
  for label in sorted(set(r['window'] for r in rows)):
   v=sorted([r for r in rows if r['window']==label],key=lambda z:z['kappa'])
   x=np.array([r['kappa']-K['shared'] for r in v]);ax[0].plot(x,[r['plane1_angle_deg'] for r in v],marker='o',label=label)
   ax[1].plot(x,[r['plane2_angle_deg'] for r in v],marker='o',label=label)
  for a in ax:a.legend(fontsize=8);a.grid(alpha=.3);a.set_xlabel('kappa0 minus shared')
  ax[0].set_ylabel('1st plane max angle vs shared (deg)');ax[1].set_ylabel('2nd plane max angle vs shared (deg)')
  fig.suptitle(f'Q550 v18.53i — true IVP principal-plane sensitivity ({seed} historical seed)')
  fig.savefig(OUT/f'principal_plane_sensitivity_{seed}.png',dpi=150);plt.close(fig)
 if transport:
  fig,ax=plt.subplots(figsize=(9,5.6),layout='constrained')
  for seed in ('original','autonomous'):
   for k in ('original','autonomous'):
    v=[r for r in transport if r['seed']==seed and r['kappa_rule']==k]
    v.sort(key=lambda z:z['later_segment'])
    ax.plot([r['later_segment'] for r in v],[r['plane2_angle_deg'] for r in v],marker='o',label=f'{seed} seed / {k} kappa')
  ax.grid(alpha=.25);ax.set_ylabel('2nd principal-plane transport angle (deg)')
  ax.set_xlabel('later segment');ax.legend(fontsize=8);fig.savefig(OUT/'principal_plane_transport.png',dpi=150);plt.close(fig)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,default=ROOT.parent)
 ap.add_argument('--native-npz',type=Path);ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--xi-max',type=float,default=51.2)
 ap.add_argument('--require-parent',action='store_true');ap.add_argument('--require-native',action='store_true')
 ap.add_argument('--require-dynamic',action='store_true')
 args=ap.parse_args()
 status=dict(version='v18.53i',scientific_decision='UNCLASSIFIED',parent_status='MISSING',native_status='MISSING',
  actual_completed_Q550_ivps=0,historical_window=[12.8,25.6],model_continuation_endpoint_xi=args.xi_max,
  continued_history_beyond_xi25_6_is_not_archived_history=True,no_full_ramp_or_noether_claim=True,
  kappa_parameter_grid=grid_kappa())
 temp=None;rc=0
 try:
  if not 25.6<=args.xi_max<=76.8:raise RuntimeError('XI_MAX_OUTSIDE_VALIDATED_AUDIT_REQUEST_RANGE')
  h=import_history();old=assert_lineage(h);status['authenticated_v1853h_fix1_reference']='OUTPUT_MANIFEST_AND_REPORT_SHA256'
  L=h.loader();parent=L.find_parent(args.parent_v1851s)
  native=args.native_npz or L.find_npz('10d',args.base)
  if parent and Path(parent).is_file():
   if sha(parent)!=CONFIG['parent_v1851s_sha256']:raise RuntimeError('V1851S_PARENT_SHA_DRIFT')
   status['parent_status']='SHA256_AUTHENTICATED'
  if native and Path(native).is_file():
   if sha(native)!=CONFIG['native_10d_sha256']:raise RuntimeError('V1852L_NATIVE_SHA_DRIFT')
   status['native_status']='SHA256_AUTHENTICATED'
  if status['parent_status']=='SHA256_AUTHENTICATED' and status['native_status']=='SHA256_AUTHENTICATED':
   temp,eng,ctxs,metadata=L.stage_original(parent);status['engine_provenance']=metadata
   print('PARENT_AND_NATIVE_AUTHENTICATED',flush=True)
   checked=[];dyn=[];sep=[];windows=[];sensitivity=[];transport=[];replay=[];refinement=[];trace=[]
   from bivector_geometry import bivector,integ
   with np.load(native,allow_pickle=False) as n:
    native_xi=np.asarray(n['common_xi'],float);xi,hcount=xi_grid(native_xi,args.xi_max)
    # Historical prefix is exactly the same grid as v18.53h (stride 64).
    assert abs(xi[hcount-1]-25.6)<1e-8
    for seed in ('original','autonomous'):
     ctx=ctxs[seed]
     old_reconcile,clockerr,archived=h.native_53g_reconciliation(n,xi[:hcount],seed,ctx,eng.omega)
     if clockerr>C['material_clock_rate_rel_tolerance']:raise RuntimeError('MATERIAL_CLOCK_RATE_FAILURE')
     replay.append(dict(seed=seed,source_g_reference_bivector_error=old_reconcile,material_clock_rel_err=clockerr))
     y0=np.asarray(archived[0],float)
     # Explicitly inspect extrapolated phase/drive: no silent flat extension.
     design=eng.original_ramp.Design('v18_50p','beta:2:2',.24)
     probe_xi=np.linspace(12.8,args.xi_max,31)
     theta=[float(ctx['theta'](ctx['event']+float(u))) for u in probe_xi]
     sign=[float(eng.clock.cubic_clock_sign(ctx['event']+float(u),design)) for u in probe_xi]
     if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(sign)) or max(abs(v) for v in sign)>1.00001:
      raise RuntimeError('HISTORICAL_CLOCK_CONTINUATION_UNDEFINED_OR_INVALID')
     # A source-generated forcing extension is not a measurement of an archived future.
     status.setdefault('clock_continuation',{})[seed]=dict(theta_start=theta[0],theta_end=theta[-1],
      clock_sign_start=sign[0],clock_sign_end=sign[-1],
      beyond_historical_endpoint_is_model_extrapolation=args.xi_max>25.6)
     dbase={}
     for kn,kappa in grid_kappa():
      z,info=h.solve(eng,ctx,kappa,xi,y0)
      status['actual_completed_Q550_ivps']+=1
      if not np.all(np.isfinite(z)) or np.max(np.abs(z[:,[0,4]]))>40 or np.max(np.abs(z[:,8]))>700:
       raise RuntimeError('DYNAMIC_STATE_NONPHYSICAL_OR_LOGIT_OVERFLOW '+seed+' '+kn)
      sourceerr=h.validated_rhs(eng,ctx,kappa,xi,z)
      if sourceerr>C['source_field_rel_tolerance']:raise RuntimeError('AUTHENTIC_SOURCE_MISMATCH '+seed+' '+kn)
      A=h.primary_metrics(xi,z,ctx,kappa)
      charge_range=A['metric']['charge_range']
      if charge_range>C['charge_abs_drift_tolerance']:raise RuntimeError('CHARGE_CONSERVATION_DRIFT '+seed+' '+kn)
      archive_err=None
      if kn=='shared':
       archive_err=float(np.max(abs(z[:hcount]-archived)/np.maximum(1,abs(archived))))
       if archive_err>C['native_replay_scaled_max']:raise RuntimeError('SHARED_NATIVE_REPLAY_FAIL '+seed+' '+str(archive_err))
      run=dict(seed=seed,kappa_rule=kn,kappa=float(kappa),xi_max=args.xi_max,
       authentic_replay_to_xi25_6=archive_err,charge_range=charge_range,
       source_rel_error=sourceerr,nfev=info['nfev'],true_integrated_IVP=True,
       historic_external_clock_extended_beyond_native=args.xi_max>25.6)
      if kn in ('original','autonomous'):
       zfine,ii=h.solve(eng,ctx,kappa,xi,y0,tight=True)
       status['actual_completed_Q550_ivps']+=1
       ferr=float(np.max(abs(zfine-z)/np.maximum(1,abs(zfine))))
       refinement.append(dict(seed=seed,kappa_rule=kn,max_scaled_error=ferr,
        nfev_normal=info['nfev'],nfev_refined=ii['nfev']))
       if ferr>C['integrator_refinement_scaled_max']:raise RuntimeError('LONG_IVP_REFINEMENT_FAILED '+seed+' '+kn+' '+str(ferr))
      if kn in ('shared','original','autonomous'):
       rgeo=embed(z[:,:10],CONFIG['params']['m_x'])[3]
       i=np.linspace(0,len(xi)-1,min(125,len(xi)),dtype=int)
       for j in i:
        trace.append(dict(seed=seed,kappa_rule=kn,xi=float(xi[j]),r1=float(rgeo['r1'][j]),
         r2=float(rgeo['r2'][j]),x=float(rgeo['x'][j]),J=float(A['J'][j])))
      dyn.append(run)
      for key,lo,hi in slices_for(xi,args.xi_max):
       xs,zz=extract_segment(xi,z,lo,hi);D=h.primary_metrics(xs,zz,ctx,kappa)
       qerr=numerical_quadrature_defect(D)
       if qerr is not None and qerr>C['residual_quadrature_relative_tolerance']:
        raise RuntimeError('BIVECTOR_QUADRATURE_MISMATCH '+seed+' '+kn+' '+key+' '+str(qerr))
       key2=(kn,key)
       dbase[key2]=D
       met=tensor_geometry(D)
       windows.append(dict(seed=seed,kappa_rule=kn,kappa=float(kappa),window=key,
        xi_start=lo,xi_end=hi,quadrature_halfgrid_rel=qerr,
        principal_degenerate=bool(met['gap']<C['finite_difference_projector_gap_warn']),**met))
       if key=='cumulative_1' and kn in ('shared','original','autonomous'):
        checks=independent_53h_reconcile(seed,kn,D['metric']);checked.append(dict(seed=seed,kappa_rule=kn,**checks))
      print('IVP_DONE',seed,kn,'k',f'{kappa:.12f}','source',f'{sourceerr:.3g}',
       'ivp_count',status['actual_completed_Q550_ivps'],flush=True)
     # Compare every parameter grid point against the same true shared-k trajectory.
     for key,lo,hi in slices_for(xi,args.xi_max):
      center=dbase['shared',key]
      for kn,kappa in grid_kappa():
       if kn=='shared':continue
       d=dbase[kn,key];pm=pair_metrics(center,d)
       sep.append(dict(seed=seed,kappa_rule=kn,kappa=float(kappa),kappa_minus_shared=float(kappa-K['shared']),
        window=key,xi_start=lo,xi_end=hi,
        sigma1=d['metric']['residual_principal_sigma1'],sigma2=d['metric']['residual_principal_sigma2'],
        cancellation=d['metric']['invariant_cancellation'],**pm))
      for m in (0.25,1,3,10):
       left=dbase[f'center_minus_{m:g}_gap',key];right=dbase[f'center_plus_{m:g}_gap',key]
       val=derivative_metrics(left,center,right,m*(K['original']-K['autonomous']))
       sensitivity.append(dict(seed=seed,window=key,step_gap_multiples=m,xi_start=lo,xi_end=hi,**val))
     # Principal-plane transport inside each same physical evolution.
     segments=[r for r in slices_for(xi,args.xi_max) if r[0].startswith('segment_')]
     for kn in ('shared','original','autonomous'):
      for (k0,a0,b0),(k1,a1,b1) in zip(segments[:-1],segments[1:]):
       left=dbase[kn,k0];right=dbase[kn,k1]
       transport.append(dict(seed=seed,kappa_rule=kn,earlier_segment=k0,later_segment=k1,
        xi_join=b0,**pair_metrics(left,right)))
   if args.xi_max==51.2 and status['actual_completed_Q550_ivps']!=C['q550_integrations_expected']:
    raise RuntimeError('UNEXPECTED_IVP_COUNT '+str(status['actual_completed_Q550_ivps']))
   # Identify derivative stability: central slope at adjacent step scales is an observable,
   # not a theorem of uniform parameter robustness.
   stiffness=[]
   for seed in ('original','autonomous'):
    for window in sorted(set(r['window'] for r in sensitivity)):
     v={r['step_gap_multiples']:r for r in sensitivity if r['seed']==seed and r['window']==window}
     a,b=v[.25],v[1]
     stiffness.append(dict(seed=seed,window=window,
      small_vs_medium_dR_slope_rel_defect=abs(a['dR_dk_fro']-b['dR_dk_fro'])/max(1,a['dR_dk_fro']),
      small_vs_medium_mode1_dP_slope_rel_defect=abs(a['dP1_dk_fro']-b['dP1_dk_fro'])/max(1,a['dP1_dk_fro']) if a['dP1_dk_fro'] is not None and b['dP1_dk_fro'] is not None else None,
      small_vs_medium_mode2_dP_slope_rel_defect=abs(a['dP2_dk_fro']-b['dP2_dk_fro'])/max(1,a['dP2_dk_fro']) if a['dP2_dk_fro'] is not None and b['dP2_dk_fro'] is not None else None))
   # Strictly descriptive geometry: a parameter-extreme can refute persistence without failing authentic integration.
   far=[r for r in sep if abs(r['kappa_minus_shared'])>5*(K['original']-K['autonomous'])]
   weak=[r for r in far if r['bivector_cosine']<C['strong_geometry_cosine_threshold'] or
      (r['plane1_angle_deg'] is not None and r['plane1_angle_deg']>C['strong_plane_max_angle_deg']) or
      (r['plane2_angle_deg'] is not None and r['plane2_angle_deg']>C['strong_plane_max_angle_deg'])]
   status['scientific_decision']=('AUTHENTIC_DUAL_KAPPA_PARAMETER_SENSITIVITY_AND_EXTENDED_MODEL_CLOCK_IVPS__STRONG_GEOMETRY_NOT_UNIFORM'
     if weak else 'AUTHENTIC_DUAL_KAPPA_PARAMETER_SENSITIVITY_AND_EXTENDED_MODEL_CLOCK_IVPS__SAMPLED_GEOMETRY_PERSISTENT')
   csv_rows('integrated_kappa_runs.csv',dyn)
   csv_rows('residual_bivector_windows.csv',windows)
   csv_rows('kappa_principal_plane_separation.csv',sep)
   csv_rows('kappa_local_sensitivity.csv',sensitivity)
   csv_rows('principal_plane_time_transport.csv',transport)
   csv_rows('sensitivity_resolution_comparison.csv',stiffness)
   csv_rows('sampled_native_radius_charge.csv',trace)
   write_json('baseline_53h_replay.json',dict(source_reconciliation=replay,residual_reconciliation=checked,refinements=refinement))
   status['summary']=dict(number_parameter_values=len(grid_kappa()),number_seeds=2,
    number_integrations=status['actual_completed_Q550_ivps'],
    number_parameter_plane_comparisons=len(sep),number_horizon_windows=len(windows),
    number_projector_derivatives=len(sensitivity),number_time_transport_comparisons=len(transport),
    weak_geometry_wider_kappa_cases=len(weak),
    max_charge_drift=max(r['charge_range'] for r in dyn),max_authentic_rhs_error=max(r['source_rel_error'] for r in dyn),
    max_numeric_refinement_error=max(r['max_scaled_error'] for r in refinement),
    max_historical_shared_replay_error=max(r['authentic_replay_to_xi25_6'] for r in dyn if r['kappa_rule']=='shared'),
    max_historical_53h_residual_replay_error=max(max(r['baseline_cancellation_abs_error'],r['baseline_sigma1_abs_error']) for r in checked),
    max_quadrature_defect=max(r['quadrature_halfgrid_rel'] or 0 for r in windows),
    min_principal_gap=min(r['gap'] for r in windows),
    source_history_expires_xi=25.6,continued_clock_end_xi=args.xi_max)
   plot_results(sep,sensitivity,transport)
  else:
   status['scientific_decision']='PARTIAL_VERIFIED_53H_REFERENCE__NO_NEW_Q550_IVPS'
 except Exception as exc:
  status['scientific_decision']='AUDIT_FAIL_NOT_CERTIFIED';status['failure']=str(exc)
  status['traceback']=traceback.format_exc();rc=2
 finally:
  if temp is not None:temp.cleanup()
  write_json('report.json',status)
  lines=['# Q550 — v18.53i: κ₀ sensitivity, principal-plane transport & long-horizon persistence','',
   '## Verdict',f"`{status['scientific_decision']}`",'',
   '## Completed integrations',str(status['actual_completed_Q550_ivps']),'',
   '## Provenance',f"v18.51s: {status['parent_status']}; native v18.52l: {status['native_status']}",'',
   '## Scientific interpretation',
   'The 12.8–25.6 interval is independently checked against the authenticated v18.53h fix1 and original native trajectory.',
   f"The 25.6–{args.xi_max} extension, if computed, is a continued ORIGINAL-SOURCE forcing protocol, NOT independently authenticated historical observations.",
   'A source-matched IVP may diverge parametrically without invalidating the original model. Report angles, derivative conditioning and time transport independently.',
   'Principal-plane projectors are meaningful away from singular-value degeneracy; Schur bases themselves have arbitrary gauge.',
   'Time-integrated bivector compensation does not certify instantaneous cancellation, global stability, or hinge Noether charge.',
   '', '## Summary',json.dumps(status.get('summary',{'not_computed':status.get('failure','inputs absent')}),indent=2)]
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
  print('V18_53I_SCIENTIFIC_DECISION',status['scientific_decision'],'IVPs',status['actual_completed_Q550_ivps'],flush=True)
  if args.require_parent and status['parent_status']!='SHA256_AUTHENTICATED':rc=max(3,rc)
  if args.require_native and status['native_status']!='SHA256_AUTHENTICATED':rc=max(3,rc)
  if args.require_dynamic and not status['scientific_decision'].startswith('AUTHENTIC_DUAL_KAPPA_PARAMETER_SENSITIVITY_AND_EXTENDED_MODEL_CLOCK_IVPS'):rc=max(3,rc)
 return rc

if __name__=='__main__':sys.exit(main())
