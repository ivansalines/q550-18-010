#!/usr/bin/env python3
"""v18.53g — observer-plane vs physical-parameter persistence of torque cancellation.

No trajectories are re-integrated. A kappa change is an instantaneous,
fixed-state source counterfactual, never a self-consistent historical orbit.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,sys,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'reference';OUT=ROOT/'results'
sys.path.insert(0,str(ROOT/'scripts'))
from torque_math import axes,embed,force_components,material_rhs,wedge,integral
C=json.loads((ROOT/'inputs_config.json').read_text())
P=C['params'];K=C['kappas']

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for part in iter(lambda:f.read(1<<20),b''):h.update(part)
 return h.hexdigest()
def put(name,data):
 OUT.mkdir(parents=True,exist_ok=True)
 (OUT/name).write_text(json.dumps(data,indent=2,allow_nan=False,ensure_ascii=False)+'\n')
def csvput(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def old_loader():
 p=REF/'v18_53f_original_run_audit.py'
 spec=importlib.util.spec_from_file_location('q550_v1853g_parent_loader',str(p))
 m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m
 spec.loader.exec_module(m)
 return m

def provenance():
 for key,filename in (
  ('v18_53f_report','v18_53f_authenticated_report.json'),
  ('v18_53f_summary','v18_53f_authenticated_native_summary.json'),
  ('v18_53f_output_manifest','v18_53f_OUTPUT_MANIFEST.json'),
  ('v18_52n_report','v18_52n_authenticated_report.json'),
  ('v18_53f_source','v18_53f_original_run_audit.py'),
  ('historical_kinetic','kinetic_metric_v18_49h.py')):
  if sha(REF/filename)!=C['source_sha256'][key]:raise RuntimeError('REFERENCE_SHA_MISMATCH: '+filename)
 report=json.loads((REF/'v18_53f_authenticated_report.json').read_text())
 prev=json.loads((REF/'v18_53f_OUTPUT_MANIFEST.json').read_text())
 for rel,fn in (
  ('results/report.json','v18_53f_authenticated_report.json'),
  ('results/native_absolute_hinge_torque_summary.json','v18_53f_authenticated_native_summary.json')):
  if prev['files'][rel]!=sha(REF/fn):raise RuntimeError('53F_REPORT_OUTPUT_MANIFEST_MISMATCH')
 if report.get('scientific_decision')!='AUTHENTIC_ORIGINAL_HINGE_TORQUE_AND_OBSERVER_TRANSPORT_LEDGER__HINGE_NOT_NOETHER_CHARGE':
  raise RuntimeError('53F_PREVIOUS_SCIENTIFIC_CERTIFICATE_NOT_ACCEPTED')
 if not report.get('absolute_hinge',{}).get('original_source_authenticated'):
  raise RuntimeError('53F_NOT_PARENT_AUTHENTICATED')
 meta=report['original_engine_provenance']
 if abs(meta['native_fixed_kappa']-K['shared'])>1e-13:raise RuntimeError('53F_KAPPA_DRIFT')
 for seed in ('autonomous','original'):
  if abs(meta['historical_seed_kappas'][seed]-K[seed])>C['kappa_relative_difference_tolerance']:
   raise RuntimeError('53F_HISTORICAL_KAPPA_DRIFT')
 summ=json.loads((REF/'v18_53f_authenticated_native_summary.json').read_text())
 if len(summ)!=C['sampled_series_expected']:raise RuntimeError('53F_SUMMARY_COUNT_DRIFT')
 return report,{z['series']:z for z in summ}

def wedge_bivector(Q,F,t):
 b=Q[:,:,None]*F[:,None,:]-F[:,:,None]*Q[:,None,:]
 return np.trapezoid(b,t,axis=0)
def cscore(a,b):
 den=abs(a)+abs(b)
 if den<1e-15:return None
 return float(np.clip(1-abs(a+b)/den,0,1))
def binvariants(A,B):
 na=float(np.linalg.norm(A));nb=float(np.linalg.norm(B));den=na+nb
 return dict(
  bivector_norm_U=na,bivector_norm_C=nb,
  bivector_cancel_score=(float(np.clip(1-np.linalg.norm(A+B)/den,0,1)) if den>1e-15 else None),
  bivector_cosine=(float(np.sum(A*B)/(na*nb)) if na*nb>1e-18 else None),
  bivector_residual_norm=float(np.linalg.norm(A+B)))

def analyse(native,loader,engine,contexts,previous):
 if sha(native)!=C['native_10d_sha256']:raise RuntimeError('NATIVE_10D_SHA_MISMATCH')
 observer=json.loads((REF/'v18_52n_authenticated_report.json').read_text())
 e,n,h=axes(C['observer_slope'],observer['pca_hidden'])
 from torque_math import acceleration_from_rhs
 table=[];intrinsic=[];source=[];correlations=[]
 max_baseline=0.;max_source=0.;max_charge=0.
 with np.load(native,allow_pickle=False) as npz:
  xi=np.asarray(npz['common_xi'],float)
  labels=sorted(k for k in npz.files if k.endswith('__Y_tau'))
  if len(labels)!=C['sampled_series_expected']:raise RuntimeError('NATIVE_SERIES_COUNT_INVALID')
  for label in labels:
   a=np.asarray(npz[label][::C['native_stride']],float)
   t=a[:,10];y=a[:,:10];xx=xi[::C['native_stride']]
   if len(t)!=len(xx) or np.min(np.diff(t))<=0:raise RuntimeError('INVALID_TIME_OR_SIZE: '+label)
   seed=label.split('_',1)[0]
   if seed not in contexts:raise RuntimeError('UNRECOGNIZED_SEED: '+label)
   context=contexts[seed]
   phase=np.asarray([context['theta'](context['event']+float(z)) for z in xx],float)
   Q,V,J,g=embed(y,P['m_x']);charge_range=float(np.ptp(J));max_charge=max(max_charge,charge_range)
   if charge_range>C['charge_drift_abs_tolerance']:raise RuntimeError('AUTHENTIC_ORBIT_CHARGE_DRIFT: '+label)
   for kn,kappa in K.items():
    _,_,_,_,D,force=force_components(y,phase,kappa,P)
    F=material_rhs(y,phase,kappa,P)
    A=acceleration_from_rhs(y,F,P['m_x'])
    if np.max(abs(A-sum(force.values())))>5e-12:raise RuntimeError('SIX_FORCE_ACCELERATION_INCONSISTENT')
    # Counterrotation charge is an algebraic identity for all kappa variants;
    # do not mistake it for trajectory-level conservation for counterfactual kappa.
    r1=g['r1'];r2=g['r2']
    jdot=(r1*r1*(2*y[:,1]*y[:,3]+F[:,3])-r2*r2*(2*y[:,5]*y[:,7]+F[:,7]))
    if np.max(abs(jdot))>1e-11:raise RuntimeError('COUNTERROTATION_DERIVATIVE_IDENTITY_FAIL')
    # Parent RHS is independent of rewritten materials: verify each variant.
    picked=np.unique(np.linspace(0,len(t)-1,C['source_probe_count_per_series'],dtype=int))
    source_error=0.
    design=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
    for j in picked:
     clock=float(engine.clock.cubic_clock_sign(context['event']+float(xx[j]),design))
     original=np.asarray(engine.physical.rhs_from_original(
      np.r_[y[j],phase[j]],kappa,engine.params,engine.slope,clock,engine.support,engine.audit),float)
     converted=original[:10]*g['x'][j]/engine.omega
     error=float(np.max(abs(converted-F[j]))/max(1.,float(np.max(abs(converted)))))
     source_error=max(source_error,error)
     if abs(float(original[10])-engine.omega*clock)>1e-10:raise RuntimeError('EXTERNAL_CLOCK_SOURCE_DRIFT')
    max_source=max(max_source,source_error)
    if source_error>C['source_rel_tolerance']:raise RuntimeError('AUTHENTIC_PARENT_RHS_DRIFT: '+label+' '+kn)
    source.append(dict(series=label,kappa_rule=kn,kappa=kappa,
        max_parent_rhs_rel_residual=source_error,n_spotchecks=len(picked),charge_derivative_max_abs=float(np.max(abs(jdot)))))
    MU=wedge_bivector(Q,force['radial_potential'],t)
    MC=wedge_bivector(Q,force['phase_coupling'],t)
    bi=binvariants(MU,MC)
    bi.update(series=label,kappa_rule=kn,kappa=kappa,charge_range_on_archived_shared_kappa_trajectory=charge_range)
    bi['invariant_cancel_strong']=(bi['bivector_cancel_score'] is not None and bi['bivector_cancel_score']>=C['cancel_score_threshold'] and bi['bivector_cosine']<=C['cancel_cosine_threshold'])
    intrinsic.append(bi)
    for angle in C['observer_angles_degrees']:
     rad=math.radians(angle)
     na=math.cos(rad)*n+math.sin(rad)*e
     # Observer plane is held physically fixed while time evolves.
     components={key:float(np.trapezoid(wedge(Q,fv,na,h),t)) for key,fv in force.items()}
     pu=components['radial_potential'];pc=components['phase_coupling'];sc=cscore(pu,pc)
     exact_u=float(na@MU@h);exact_c=float(na@MC@h)
     projection_error=max(abs(exact_u-pu),abs(exact_c-pc))
     if projection_error>1e-10:raise RuntimeError('BIVECTOR_PROJECTION_IDENTITY_FAIL')
     u=wedge(Q,force['radial_potential'],na,h)
     v=wedge(Q,force['phase_coupling'],na,h)
     local_opposition=float(np.mean((u*v)<0))
     local_integrated_abs_ratio=(float(np.trapezoid(abs(u+v),t))/max(1e-14,float(np.trapezoid(abs(u)+abs(v),t))))
     pscore=1-local_integrated_abs_ratio
     same_sign=bool(pu*pc<0)
     strong=bool(same_sign and sc is not None and sc>=C['cancel_score_threshold'])
     row=dict(series=label,seed=seed,kappa_rule=kn,kappa=float(kappa),angle_deg=angle,
       source_type=('SHARED_KAPPA_TRUE_NATIVE_RHS' if kn=='shared' else 'FROZEN_STATE_KAPPA_COUNTERFACTUAL_NOT_NATIVE_TRAJECTORY'),
       integral_radial_potential=pu,integral_phase_coupling=pc,
       integral_other_4_channels=sum(v for key,v in components.items() if key not in ('radial_potential','phase_coupling')),
       integral_all_6_channels=sum(components.values()),
       signed_pair_cancel_score=sc,
       signed_pair_opposite_sign=same_sign,
       strong_projected_cancellation=strong,
       fraction_samples_opposite_pointwise=local_opposition,
       integrated_absolute_local_cancel_score=pscore,
       bivector_projection_max_error=projection_error)
     for k2,v2 in components.items():row['integral_'+k2]=v2
     table.append(row)
     if kn=='shared' and angle==0:
      old=previous[label]['force_torque_components_integrated']
      err=max(abs(pu-old['radial_potential']),abs(pc-old['phase_coupling']))
      max_baseline=max(max_baseline,err)
      if err>C['baseline_integral_abs_tolerance']:
       raise RuntimeError('V18_53F_BASELINE_RECONCILIATION_FAILED: '+label+' '+str(err))
    # At frozen state and identical external phase, force C is linear in kappa.
    correlations.append(dict(series=label,kappa_rule=kn,kappa=kappa,
      old_shared_kappa=K['shared'],source_phase_authenticated=True,
      observer_plane_and_state_held_fixed=True))
 if len(table)!=len(labels)*len(K)*len(C['observer_angles_degrees']):raise RuntimeError('ROWS_COUNT_WRONG')
 csvput('observer_kappa_torque_persistence.csv',table)
 csvput('coordinate_free_bivector_persistence.csv',intrinsic)
 csvput('authenticated_parent_rhs_spotchecks.csv',source)
 orig=[r for r in table if r['kappa_rule'] in ('original','autonomous')]
 fails=[dict(series=r['series'],kappa_rule=r['kappa_rule'],angle_deg=r['angle_deg'],score=r['signed_pair_cancel_score'],
   pair_opposite_sign=r['signed_pair_opposite_sign']) for r in orig if not r['strong_projected_cancellation']]
 fails_invariant=[dict(series=r['series'],kappa_rule=r['kappa_rule'],score=r['bivector_cancel_score'],cosine=r['bivector_cosine'])
   for r in intrinsic if r['kappa_rule'] in ('original','autonomous') and not r['invariant_cancel_strong']]
 allproj=bool(len(fails)==0)
 allinv=bool(len(fails_invariant)==0)
 metrics=dict(n_native_series=len(labels),n_projection_rows=len(table),n_invariant_rows=len(intrinsic),
  original_autonomous_projected_strong_count=len(orig)-len(fails),
  original_autonomous_projected_total=len(orig),
  original_autonomous_bivector_strong_count=len([r for r in intrinsic if r['kappa_rule'] in ('original','autonomous')])-len(fails_invariant),
  max_source_field_relative_error=max_source,max_previous_baseline_integral_abs_error=max_baseline,max_native_charge_drift=max_charge,
  persist_all_predeclared_planes_both_kappas=allproj,
  persist_coordinate_free_bivector_both_kappas=allinv,
  persistence_fails_sample=fails[:50],bivector_fails_sample=fails_invariant[:50],
  counterfactual_kappa_not_new_trajectory=True,
  no_global_hinge_Noether_claim=True)
 put('persistence_metrics.json',metrics)
 return table,intrinsic,metrics

def plot(table,intrinsic):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return
 selected=[r for r in table if r['series'] in ('original__reference__Y_tau','autonomous__reference__Y_tau')]
 if not selected:return
 fig,ax=plt.subplots(2,1,figsize=(10.4,8.0),layout='constrained')
 for seed in ('original','autonomous'):
  for kn in ('autonomous','original'):
   z=sorted((r for r in selected if r['seed']==seed and r['kappa_rule']==kn),key=lambda r:r['angle_deg'])
   if not z:continue
   ax[0].plot([r['angle_deg'] for r in z],[r['signed_pair_cancel_score'] for r in z],'.-',label=f'{seed} state / κ {kn}')
 ax[0].axhline(C['cancel_score_threshold'],color='gray',ls='--',label='predeclared strong-cancellation threshold')
 ax[0].set_ylabel('Integrated signed pair cancellation score');ax[0].set_ylim(-.03,1.03)
 for seed in ('original','autonomous'):
  for kn in ('autonomous','original'):
   z=next((r for r in intrinsic if r['series']==seed+'__reference__Y_tau' and r['kappa_rule']==kn),None)
   if z:ax[1].scatter([K[kn]],[z['bivector_cancel_score']],s=80,label=f'{seed} state / κ {kn}')
 ax[1].set_xlabel('κ₀');ax[1].set_ylabel('Coordinate-free 5D bivector cancellation')
 ax[1].axhline(C['cancel_score_threshold'],color='gray',ls='--')
 for a in ax:a.legend(fontsize=8);a.grid(alpha=.25)
 fig.suptitle('Q550 v18.53g — plane-dependent vs coordinate-free cancellation')
 fig.savefig(OUT/'torque_cancellation_observer_vs_invariant.png',dpi=155);plt.close(fig)

def preview_from_previous(prev):
 # Linear in kappa on identical states; this is not a recalculated orbit or parent verification.
 rows=[]
 for series,old in prev.items():
  base=old['force_torque_components_integrated'];u=base['radial_potential'];v=base['phase_coupling']
  for kn,k in K.items():
   vc=v*k/K['shared'];rows.append(dict(series=series,kappa_rule=kn,kappa=k,
     integrated_radial_from_v1853f=u,scaled_integrated_coupling=vc,
     predicted_projected_cancel_score=cscore(u,vc),
     status='ALGEBRAIC_KAPPA_RESCALE_ON_ARCHIVED_SHARED_KAPPA_STATES_NOT_NEW_ORBIT'))
 csvput('preview_53f_fixed_plane_analytic_kappa_rescale.csv',rows)
 return rows

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--base',type=Path,default=ROOT.parent)
 ap.add_argument('--native-npz',type=Path)
 ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--require-native',action='store_true')
 ap.add_argument('--require-parent',action='store_true')
 args=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 try:
  old,previous=provenance()
  preview=preview_from_previous(previous)
  l=old_loader()
  native=args.native_npz or l.find_npz('10d',args.base)
  parent=l.find_parent(args.parent_v1851s)
  status=dict(version=C['version'],scientific_decision='NOT_CLASSIFIED',
    native_status='MISSING',parent_status='MISSING',
    source_v18_53f_status='AUTHENTICATED_SHA_AND_MANIFEST',
    frozen_kappa_counterfactual_not_new_orbits=True,
    no_noether_hinge_charge_certificate=True,
    kappa_rules=K,observer_angles_degrees=C['observer_angles_degrees'],
    n_prev_analytic_preview_rows=len(preview),new_Q550_orbit_integrations=0)
  temp=None
  try:
   if native and native.is_file():
    if sha(native)!=C['native_10d_sha256']:raise RuntimeError('NATIVE_SHA_MISMATCH')
    status['native_status']='SHA256_AUTHENTICATED'
   if parent and Path(parent).is_file():
    if sha(parent)!=C['parent_v1851s_sha256']:raise RuntimeError('PARENT_SHA_MISMATCH')
    status['parent_status']='SHA256_AUTHENTICATED'
   if status['native_status']=='SHA256_AUTHENTICATED' and status['parent_status']=='SHA256_AUTHENTICATED':
    temp,eng,contexts,metadata=l.stage_original(parent)
    rows,intrinsic,metrics=analyse(native,l,eng,contexts,previous)
    put('original_generator_provenance.json',metadata)
    status['persistence']=metrics
    status['scientific_decision']=('AUTHENTIC_FIXED_STATE_DUAL_KAPPA_ALL_PREDECLARED_PROJECTED_AND_5D_BIVECTOR_CANCELLATION_PERSISTENCE__NOT_NEW_ORBITS' if
       metrics['persist_all_predeclared_planes_both_kappas'] and metrics['persist_coordinate_free_bivector_both_kappas']
       else 'AUTHENTIC_FIXED_STATE_DUAL_KAPPA_PERSISTENCE_MAPPED_WITH_OBSERVER_OR_INVARIANT_FAILURES__NOT_NEW_ORBITS')
    plot(rows,intrinsic)
   else:
    status['scientific_decision']='PARTIAL_MISSING_AUTHENTIC_SOURCE_OR_NATIVE_STATES'
   put('report.json',status)
   report=['# Q550 v18.53g — Torque-Cancellation Persistence','',
    '## Classification',f"`{status['scientific_decision']}`",'',
    '## Source status',f"Native 10D: {status['native_status']}; original v18.51s parent: {status['parent_status']}",
    '','## Scope and scientific boundary',
    'The historical original and autonomous series are native to SHARED κ₀. All evaluations at historical original/autonomous κ₀ are fixed-state fixed-external-phase counterfactuals, not self-consistent trajectories.',
    'Projected cancellation can depend on the chosen observer plane. An additional coordinate-free integrated 5D bivector test is used and is not implied by a high single-plane score.',
    'Temporal pointwise cancellation and signed integrated cancellation are separately reported.',
    'No global Noether hinge charge or true orbital closure is certified.','',
    '## Results',json.dumps(status.get('persistence',{'status':'not evaluated because parent or native states missing'}),indent=2)]
   (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(report)+'\n')
   print('V18_53G_STATUS',status['scientific_decision'])
   if args.require_parent and status['parent_status']!='SHA256_AUTHENTICATED':raise RuntimeError('REQUIRED_PARENT_MISSING')
   if args.require_native and status['native_status']!='SHA256_AUTHENTICATED':raise RuntimeError('REQUIRED_NATIVE_MISSING')
  finally:
   if temp is not None:temp.cleanup()
 except Exception:
  (OUT/'audit_failure.txt').write_text(traceback.format_exc())
  put('failure_summary.json',{'version':C['version'],'classification':'FAILED_INPUT_INTEGRITY_OR_REQUIRED_SOURCE',
    'traceback_file':'audit_failure.txt'})
  raise
if __name__=='__main__':main()
