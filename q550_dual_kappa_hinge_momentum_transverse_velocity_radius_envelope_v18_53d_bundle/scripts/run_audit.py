#!/usr/bin/env python3
"""v18.53d. Authenticated fold kinematics plus relative 5D observer-depth audit.

Crucial separation: fold rotor J1-J2 is a real model first integral. A wedge
computed on matched 5D relative dQ,dV is a kinematic observer diagnostic,
NOT canonical J_hinge, not Noether, and not absolute orbital momentum.
"""
import argparse,csv,json,math,hashlib,traceback
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parents[1];REF=R/'reference';OUT=R/'results'
C=json.loads((R/'inputs_config.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def emit(name,data):
 OUT.mkdir(exist_ok=True,parents=True)
 (OUT/name).write_text((data if isinstance(data,str) else json.dumps(data,indent=2,allow_nan=False,ensure_ascii=False))+'\n')
def csvout(name,rows):
 if not rows:return
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def unit(q):
 q=np.asarray(q,float);n=float(np.linalg.norm(q))
 if n<1e-15 or not np.isfinite(n):raise ValueError('Degenerate observer frame')
 return q/n
def frame(slope,hidden,angle):
 e=unit([1,slope,0,0,0]);n=unit([-slope,1,0,0,0]);h=unit([0,0,*hidden])
 a=math.radians(angle);w=math.cos(a)*n+math.sin(a)*h;depth=-math.sin(a)*n+math.cos(a)*h
 assert np.linalg.norm(np.stack([e,w,depth])@np.stack([e,w,depth]).T-np.eye(3))<1e-13
 return e,n,h,w,depth
def wedge(q,v,a,b):return (q@a)*(v@b)-(q@b)*(v@a)
def rotor(r1,r2,w1,w2):
 j1=r1*r1*w1;j2=r2*r2*w2
 return dict(vtheta1=r1*w1,vtheta2=r2*w2,J1=j1,J2=j2,J=j1-j2)
def lineage():
 m=json.loads((REF/'v18_53c_output_manifest.json').read_text())
 for f,p in [('v18_53c_authenticated_report.json','results/report.json'),('v18_53c_authenticated_fold_locus.csv','results/fold_locus_continuation.csv'),('v18_53c_authenticated_fold_states.json','results/fold_locus_selected_states.json')]:
  if sha(REF/f)!=m['files'][p]:raise RuntimeError('AUTHENTIC_FOLD_HASH_MISMATCH '+f)
 m=json.loads((REF/'v18_52n_output_manifest.json').read_text())
 if sha(REF/'v18_52n_authenticated_report.json')!=m['files_sha256']['results/report.json']:raise RuntimeError('OBSERVER_ANCESTOR_HASH_MISMATCH')
 p=json.loads((REF/'v18_53c_authenticated_report.json').read_text())
 a=json.loads((REF/'v18_52n_authenticated_report.json').read_text())
 if not p['scientific_decision'].startswith('AUTHENTIC_DUAL_KAPPA_FOLD_LOCUS') or a['stage']!='AUTHENTIC_FULL_POSTPROCESS':raise RuntimeError('ANCESTOR_NOT_AUTHENTICATED')
 if abs(p['authenticated_v18_53b_anchor']['J']-C['charge_reference'])>1e-13:raise RuntimeError('CHARGE_LINEAGE_MISMATCH')
 return p,a
def fold_analysis(parent):
 with (REF/'v18_53c_authenticated_fold_locus.csv').open() as f:source=list(csv.DictReader(f))
 rows=[]
 for x in source:
  r1,r2,w1,w2=map(float,(x['r1_fold'],x['r2_fold'],x['omega1_fold'],x['omega2_fold']))
  mech=rotor(r1,r2,w1,w2)
  if abs(mech['J']-float(x['J']))>2e-12 or x['native_checked']!='True' or x['native_source_match']!='True':raise RuntimeError('UNCERTIFIED_FOLD_ROW')
  rows.append(dict(region=x['region'],kappa=float(x['kappa']),signed_drive=float(x['signed_drive_fold']),r1=r1,r2=r2,omega1=w1,omega2=w2,r_min=min(r1,r2),dr2_dk=float(x['dr2_dk']),**mech))
 csvout('fold_rotor_kinematics.csv',rows)
 dual=sorted([x for x in rows if x['region']=='dual'],key=lambda y:y['kappa'])
 if len(dual)!=parent['continuation_counts']['dual']:raise RuntimeError('INCOMPLETE_DUAL_FOLD_LOCUS')
 a,b=dual[0],dual[-1]
 for key,v in [('autonomous',a),('original',b)]:
  if abs(parent['endpoints'][key]['r2_fold']-v['r2'])>1e-12:raise RuntimeError('FOLD_ENDPOINT_MISMATCH')
 delta={k:b[k]-a[k] for k in ('r1','r2','omega1','omega2','vtheta1','vtheta2','J1','J2','J')}
 # At common s=-1, the exact phase-locking relation is omega1+omega2=-x=-1/2.
 # Since the authenticated v18.53c comparison stores radii, recover the
 # two angular rates algebraically from the SAME charge (no new model fit).
 same_drive={}
 for label in ('autonomous','original'):
  entry=parent['fixed_common_drive_comparison'][label]
  r1=float(entry['r1']);r2=float(entry['r2']);J=C['charge_reference']
  w1=(J-.5*r2*r2)/(r1*r1+r2*r2);w2=-.5-w1
  mec=rotor(r1,r2,w1,w2)
  if abs(mec['J']-J)>2e-14:raise RuntimeError('COMMON_DRIVE_CHARGE_MISMATCH')
  same_drive[label]=dict(r1=r1,r2=r2,omega1=w1,omega2=w2,**mec,source='v18.53c authenticated radii + analytic fixed-drive phase-locking identity')
 spans={k:[min(v[k] for v in dual),max(v[k] for v in dual)] for k in ('r2','omega2','vtheta2','J1','J2','J')}
 rep=dict(status='AUTHENTICATED_V18_53C_FOLD_POSTPROCESS',n_all=len(rows),n_dual=len(dual),
  kappa_autonomous=a['kappa'],kappa_original=b['kappa'],fold_autonomous=a,fold_original=b,
  difference_original_minus_autonomous=delta,sampled_spans=spans,
  same_fixed_drive_s_minus1=parent['fixed_common_drive_comparison'],
  same_drive_rotor_mechanics_reconstructed_from_charge= same_drive,
  note='Each fold occurs at its own signed drive. These radii are not dynamical minimum radii or global noncollapse bounds.',
  J_hinge_canonical_status='NOT_IDENTIFIABLE_FROM_FOLD_R1_R2_OMEGA1_OMEGA2')
 emit('fold_rotor_report.json',rep);return rep,dual

def closest_separated_pair(t,q,v,e,w,view_depth,scale):
 N=len(t);n=max(10,int(N*C['early_late_fraction']))
 early=np.arange(n);late=np.arange(N-n,N)
 scr=np.column_stack([q@e,q@w]);tree=cKDTree(scr[late]/scale)
 d,ix=tree.query(scr[early]/scale)
 for pos in np.argsort(d):
  i=int(early[pos]);j=int(late[ix[pos]])
  if t[j]-t[i] < C['min_separation_material_time']:continue
  dq=q[j]-q[i];full=float(np.linalg.norm(dq));shadow=float(np.linalg.norm([dq@e,dq@w]));unseen=math.sqrt(max(0,full*full-shadow*shadow))
  assert abs(full**2-shadow**2-unseen**2)<1e-14
  return dict(i=i,j=j,tau_a=float(t[i]),tau_b=float(t[j]),delta_tau=float(t[j]-t[i]),
   screen_gap=shadow,full_gap=full,unseen_gap=unseen,depth_signed=float(dq@view_depth),
   screen_gap_rel=shadow/scale,full_gap_rel=full/scale,unseen_gap_rel=unseen/scale,
   depth_signed_rel=float(dq@view_depth)/scale,full_velocity_mismatch=float(np.linalg.norm(v[j]-v[i])))
 raise RuntimeError('NO_TIME_SEPARATED_PAIR')

def trace_analysis(source,observer):
 if sha(source)!=C['v18_52l_5d_expected_sha256']:raise RuntimeError('5D_TRACE_SHA_MISMATCH')
 hidden=observer['pca_hidden'];slope=C['observer_slope_hypothesis']
 if abs(observer['used_nodal_slope']-slope)>1e-13:raise RuntimeError('OBSERVER_SLOPE_MISMATCH')
 rows=[];hist=[]
 with np.load(source,allow_pickle=False) as src:
  labels=sorted({k.split('__')[0] for k in src.files if k.endswith('__dQ')})
  if len(labels)!=C['v18_52l_expected_kicks']:raise RuntimeError('UNEXPECTED_TRACE_LABELS')
  for label in labels:
   sl=slice(None,None,int(C['observer_stride']))
   tau=np.asarray(src[label+'__tau'][sl],float);q=np.asarray(src[label+'__dQ'][sl],float);v=np.asarray(src[label+'__dV'][sl],float)
   if q.shape!=(len(tau),5) or q.shape!=v.shape or not np.all(np.diff(tau)>0):raise RuntimeError('BAD_TRACE_SHAPE')
   scale=float(np.sqrt(np.mean(np.sum(q*q,axis=1))))
   e,n,h,w,dd=frame(slope,hidden,0)
   L=wedge(q,v,n,h);rad=np.hypot(q@n,q@h);vt=np.divide(L,rad,out=np.full_like(L,np.nan),where=rad>1e-14)
   hist.append(dict(label=label,source_is_relative_dQ_dV=True,n_samples=len(tau),relative_q_rms=scale,
    kinematic_hinge_wedge_span=float(np.ptp(L)),kinematic_hinge_wedge_min=float(L.min()),kinematic_hinge_wedge_max=float(L.max()),
    hinge_plane_radius_range=[float(rad.min()),float(rad.max())],hinge_plane_transverse_speed_finite_fraction=float(np.isfinite(vt).mean())))
   for angle in C['observer_tilts_degrees']:
    e,n,h,w,depth=frame(slope,hidden,angle)
    c=closest_separated_pair(tau,q,v,e,w,depth,scale)
    i,j=c.pop('i'),c.pop('j')
    screenL=wedge(q,v,e,w)
    flag=bool(c['screen_gap_rel']<C['approximate_shadow_closure_threshold'] and c['unseen_gap_rel']>C['min_hidden_residual_threshold'])
    rows.append(dict(label=label,angle_deg=angle,q_rms=scale,**c,
     kinematic_plane_wedge_at_a=float(L[i]),kinematic_plane_wedge_at_b=float(L[j]),
     screen_kinematic_wedge_at_a=float(screenL[i]),screen_kinematic_wedge_at_b=float(screenL[j]),
     hinge_plane_radius_at_a=float(rad[i]),hinge_plane_radius_at_b=float(rad[j]),
     hinge_transverse_speed_at_a=float(vt[i]),hinge_transverse_speed_at_b=float(vt[j]),
     prespecified_view_false_closure_flag=flag))
 csvout('physical_5d_relative_projection_shadow_hidden_depth.csv',rows)
 emit('relative_5d_hinge_wedge_summary.json',hist)
 rep=dict(status='AUTHENTICATED_V18_52L_RELATIVE_5D_POSTPROCESS',trace_hash=sha(source),n_labels=len(labels),n_observer_comparisons=len(rows),
  n_near_shadow_but_not_full=sum(int(x['prespecified_view_false_closure_flag']) for x in rows),
  source_is_relative_5d_perturbation_minus_reference=True,
  not_a_test_of_absolute_orbit_closure=True,
  observer_geometry_from_v18_52n=True,
  h_pca_is_an_observer_direction_not_dynamically_proven_hinge=True,
  kinetic_wedge_is_not_canonical_J_hinge=True,
  projection_ambiguous_even_when_plane_shadow_closes=True)
 return rep,rows

def plots(dual,shadows):
 try:
  import matplotlib;matplotlib.use('Agg')
  import matplotlib.pyplot as plt
 except ImportError:return
 k=[d['kappa'] for d in dual]
 fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained')
 for a,key,title in zip(ax.flat,('r2','omega2','vtheta2','J2'),('Critical smaller radius','Rotor-2 angular rate','Rotor-2 transverse speed','Rotor-2 angular momentum')):
  a.plot(k,[r[key] for r in dual]);a.set(title=title,xlabel='kappa0, each at own fold',ylabel=key);a.grid(alpha=.2)
 fig.suptitle('Authentic v18.53c fold kinematics — distinct kappa prescriptions');fig.savefig(OUT/'dual_kappa_rotor_radius_speed_momentum.png',dpi=140);plt.close(fig)
 if not shadows:return
 sel=sorted((r for r in shadows if r['label']=='original_negative_1.0e-04'),key=lambda r:r['angle_deg'])
 fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained');angles=[r['angle_deg'] for r in sel]
 for name,label in [('screen_gap_rel','Shadow (2D)'),('full_gap_rel','Full (5D)')]:ax[0].plot(angles,[r[name] for r in sel],'.-',label=label)
 for name,label in [('unseen_gap_rel','Unseen 5D'),('depth_signed_rel','Signed viewing depth')]:ax[1].plot(angles,[r[name] for r in sel],'.-',label=label)
 for a in ax:a.set(xlabel='Prespecified inclination (degrees)',ylabel='Gap / RMS ||dQ||');a.grid(alpha=.2);a.legend()
 fig.suptitle('One authentic Q550 relative 5D trajectory, several observers (not absolute orbit)');fig.savefig(OUT/'tilted_observer_shadow_vs_hidden_5d_gap.png',dpi=140);plt.close(fig)

def locate_trace(base):
 p=Path(base);f='extended_horizon_5d_vectors_and_traces.npz'
 return next((x for x in [p/'q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle/results'/f,p/'q550_v18_52l_review'/f,p/f] if x.is_file()),None)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--base',default=str(R.parent));ap.add_argument('--trace-npz',default='');ap.add_argument('--require-trace',action='store_true');args=ap.parse_args()
 try:
  report,observer=lineage();fold,dual=fold_analysis(report)
  source=Path(args.trace_npz) if args.trace_npz else locate_trace(args.base)
  obs,shadows=(trace_analysis(source,observer) if source is not None and source.is_file() else (None,[]))
  if args.require_trace and obs is None:raise RuntimeError('AUTHENTIC_v18_52l_5D_TRACE_UNAVAILABLE')
  plots(dual,shadows)
  decision='AUTHENTIC_FOLD_KINEMATICS_AND_RELATIVE_5D_OBSERVER_AUDIT__CANONICAL_HINGE_MOMENTUM_NOT_CERTIFIED' if obs else 'AUTHENTIC_FOLD_KINEMATICS__RELATIVE_5D_TRACE_MISSING'
  summary=dict(version='v18.53d',scientific_decision=decision,fold=fold,relative_5d_observer=obs,
    canonical_J_hinge='NOT_IDENTIFIABLE_WITH_AVAILABLE_MOMENTUM_AND_METRIC',
    Noether_J1_minus_J2_is_distinct_from_5d_observer_wedge=True,
    global_radius_lower_bound=False,full_physical_orbit_closure=False,
    dual_fold_kappas_are_not_the_historical_shared_kappa_of_5d_traces=True)
  emit('report.json',summary)
  emit('SCIENTIFIC_REPORT_SERVER.md',f'''# Q550 v18.53d — Scientific report\n\nDecision: `{decision}`.\n\nThe {fold['n_dual']} authenticated sampled dual-kappa fold points provide absolute rotor radii, angular rates, tangential velocities and angular momenta; J1-J2 stays equal to the fixed charge. Fold drive is different at each kappa. No global dynamical radius minimum is claimed.\n\n5D relative observations: {('verified physical v18.52l traces' if obs else 'trace missing')}. This dataset contains perturbation-minus-reference dQ/dV, NOT an absolute orbit. Prescribed angles can change projected closure even though the full 5D state is untouched. Projection cannot prove full return.\n\nThe kinematic wedge (dQ.n)(dV.h)-(dQ.h)(dV.n) is an observer-conditioned area rate. It is NOT canonical hinge momentum without the genuine (q,p) and kinetic metric; conservation is not claimed. The 5D inherited runs used one shared kappa0, not two independent dynamical kappa prescriptions.\n''')
  print('SCIENTIFIC_DECISION',decision);return 0
 except Exception as ex:
  emit('report.json',dict(version='v18.53d',scientific_decision='AUDIT_INCOMPLETE',error=str(ex),traceback=traceback.format_exc()))
  print('AUDIT_INCOMPLETE',ex);return 2
if __name__=='__main__':raise SystemExit(main())
