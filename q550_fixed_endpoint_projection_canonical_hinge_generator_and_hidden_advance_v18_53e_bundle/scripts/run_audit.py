#!/usr/bin/env python3
"""Q550 v18.53e: fixed-time-pair observer test and absolute kinetic hinge moment.

Never confuse a relative dQ/dV wedge with absolute moment, or assume a
kinetic candidate is a conserved Noether charge of the driven Q550 model.
No Q550 ODE is integrated here.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,traceback
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'reference'; OUT=ROOT/'results'
CONFIG=json.loads((ROOT/'inputs_config.json').read_text())

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(1<<20),b''):h.update(x)
 return h.hexdigest()
def writejson(n,obj):
 OUT.mkdir(exist_ok=True,parents=True)
 (OUT/n).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def writecsv(n,rows):
 if not rows:return
 with (OUT/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def unit(v):
 v=np.asarray(v,float);l=np.linalg.norm(v)
 if l<1e-14 or not np.isfinite(l):raise ValueError('DEGENERATE_FRAME')
 return v/l

def axes(slope,hh,alpha):
 e=unit([1,slope,0,0,0]);n=unit([-slope,1,0,0,0]);h=unit([0,0,*hh]);a=np.deg2rad(alpha)
 w=np.cos(a)*n+np.sin(a)*h;d=-np.sin(a)*n+np.cos(a)*h
 assert np.max(abs(np.stack([e,w,d])@np.stack([e,w,d]).T-np.eye(3)))<1e-13
 return e,n,h,w,d

def signed_wedge(q,p,n,h):
 return (q@n)*(p@h)-(q@h)*(p@n)

def lineage():
 # Two-level proof that the frozen previous candidate rows, fold results and
 # observer frame are unchanged from their immediately previous sendbacks.
 mani=json.loads((REF/'v18_53d_OUTPUT_MANIFEST.json').read_text())['files']
 pairs={'v18_53d_report.json':'results/report.json',
        'v18_53d_projection_candidates.csv':'results/physical_5d_relative_projection_shadow_hidden_depth.csv',
        'v18_53d_fold_rotor_kinematics.csv':'results/fold_rotor_kinematics.csv',
        'v18_53d_fold_rotor_report.json':'results/fold_rotor_report.json'}
 for local,original in pairs.items():
  if sha(REF/local)!=mani[original]:raise RuntimeError('V18_53D_OUTPUT_LINEAGE_MISMATCH: '+local)
  if sha(REF/local)!=CONFIG['v18_53d_reference_sha256'][local]:raise RuntimeError('V18_53D_PIN_MISMATCH: '+local)
 if sha(REF/'kinetic_metric_v18_49h.py')!=CONFIG['kinetic_metric_v18_49h_expected_sha256']:
  raise RuntimeError('KINETIC_METRIC_SOURCE_SHA_MISMATCH')
 old=json.loads((REF/'v18_53d_report.json').read_text())
 if not old['scientific_decision'].startswith('AUTHENTIC_FOLD_KINEMATICS_AND_RELATIVE_5D_OBSERVER_AUDIT'):
  raise RuntimeError('PARENT_53D_NOT_AUTHENTIC')
 m52n=json.loads((REF/'v18_52n_output_manifest.json').read_text())['files_sha256']
 if sha(REF/'v18_52n_authenticated_report.json')!=m52n['results/report.json']:
  raise RuntimeError('V18_52N_OBSERVER_GENEALOGY_MISMATCH')
 ob=json.loads((REF/'v18_52n_authenticated_report.json').read_text())
 if abs(ob['used_nodal_slope']-CONFIG['observer_slope_hypothesis'])>1e-13:
  raise RuntimeError('OBSERVER_SLOPE_DIFFERENT')
 return old,ob

def load_candidate_pairs():
 ans={}
 with (REF/'v18_53d_projection_candidates.csv').open(newline='') as f:
  for row in csv.DictReader(f):
   ang=float(row['angle_deg'])
   if ang in CONFIG['anchor_angles_degrees']:ans[(row['label'],ang)]=row
 return ans

def lookup(t,v):
 i=int(np.argmin(abs(t-v)))
 if abs(float(t[i])-float(v))>max(1e-9,2e-7*float(np.median(np.diff(t)))):
  raise RuntimeError(f'PREVIOUS_ENDPOINT_NOT_ON_NATIVE_SUBSAMPLE: {v}, got {t[i]}')
 return i

def intrinsic_pair(t,q,scale):
 N=len(t);n=max(10,int(N*CONFIG['early_late_fraction']))
 early=np.arange(n);late=np.arange(N-n,N)
 dist,ix=cKDTree(q[late]/scale).query(q[early]/scale)
 for k in np.argsort(dist):
  i=int(early[k]);j=int(late[ix[k]])
  if t[j]-t[i]>=CONFIG['min_separation_material_time']:return i,j
 raise RuntimeError('NO_INTRINSIC_PAIR')

def fixed_endpoint(trace,ob):
 if sha(trace)!=CONFIG['v18_52l_5d_expected_sha256']:raise RuntimeError('5D_TRACE_SHA_MISMATCH')
 pairs=load_candidate_pairs(); rows=[]; comparisons=[]
 slope=CONFIG['observer_slope_hypothesis'];h=ob['pca_hidden']
 with np.load(trace,allow_pickle=False) as z:
  labels=sorted({k.split('__')[0] for k in z.files if k.endswith('__dQ')})
  if len(labels)!=8:raise RuntimeError('UNEXPECTED_RELATIVE_TRACE_COUNT')
  for label in labels:
   sl=slice(None,None,int(CONFIG['observer_stride']))
   t=np.asarray(z[f'{label}__tau'][sl],float)
   q=np.asarray(z[f'{label}__dQ'][sl],float);v=np.asarray(z[f'{label}__dV'][sl],float)
   if len(t)<100 or q.shape!=(len(t),5) or v.shape!=q.shape or np.min(np.diff(t))<=0:
    raise RuntimeError('BAD_RELATIVE_TRACE_ARRAY')
   scale=float(np.sqrt(np.mean(np.einsum('ij,ij->i',q,q))))
   intrinsic_i,intrinsic_j=intrinsic_pair(t,q,scale)
   anchor=[('intrinsic_5d',intrinsic_i,intrinsic_j)]
   for a in CONFIG['anchor_angles_degrees']:
    row=pairs.get((label,a))
    if row is None:raise RuntimeError('MISSING_V18_53D_ANCHOR: '+str((label,a)))
    anchor.append((f'v18_53d_{a:+g}deg',lookup(t,float(row['tau_a'])),lookup(t,float(row['tau_b']))))
   for kind,i,j in anchor:
    if t[j]-t[i]<CONFIG['min_separation_material_time']:raise RuntimeError('BAD_TIME_SEPARATION')
    dq=q[j]-q[i];dv=v[j]-v[i]
    D=float(np.linalg.norm(dq));Vm=float(np.linalg.norm(dv))
    e,n,hh,w0,d0=axes(slope,h,0)
    inplane_e=float(dq@e);inplane_n=float(dq@n);inplane_h=float(dq@hh)
    other=math.sqrt(max(0,D*D-inplane_e**2-inplane_n**2-inplane_h**2))
    # Analytically optimized observer is explicitly data-selected and cannot
    # be a pre-registered physical hinge angle.
    a_opt=float(math.degrees(math.atan2(-inplane_n,inplane_h))) if np.hypot(inplane_n,inplane_h)>1e-15 else 0.
    for a in [*CONFIG['observer_tilts_degrees'],a_opt]:
     e,n,hh,w,d=axes(slope,h,a)
     dx=float(dq@e);dy=float(dq@w);depth=float(dq@d)
     screen=math.hypot(dx,dy)
     hidden=math.sqrt(max(0,D*D-screen*screen))
     depth_beyond=math.sqrt(max(0,hidden*hidden-depth*depth))
     if abs((screen*screen+depth*depth+depth_beyond*depth_beyond)-D*D)>5e-13*max(1,D*D):
      raise RuntimeError('PROJECTION_PYTHAGOREAN_IDENTITY_FAIL')
     la=signed_wedge(q[i],v[i],n,hh);lb=signed_wedge(q[j],v[j],n,hh)
     rows.append(dict(label=label,anchor=kind,angle_deg=float(a),angle_selection='OPTIMIZED_AFTER_ENDPOINT_SELECTION' if a==a_opt and a not in CONFIG['observer_tilts_degrees'] else ('OPTIMIZED_AFTER_ENDPOINT_SELECTION' if abs(a-a_opt)<1e-12 and a_opt in CONFIG['observer_tilts_degrees'] else 'PRE_REGISTERED'),
       tau_a=float(t[i]),tau_b=float(t[j]),delta_tau=float(t[j]-t[i]),
       q_rms=scale,screen_gap=screen,full_5d_gap=D,hidden_gap=hidden,signed_depth=depth,other_hidden_gap=depth_beyond,
       screen_gap_rel=screen/scale,full_gap_rel=D/scale,hidden_gap_rel=hidden/scale,velocity_5d_gap=Vm,
       signed_hinge_kinematic_wedge_at_a=float(la),signed_hinge_kinematic_wedge_at_b=float(lb),
       gap_identity_error=float(screen*screen+hidden*hidden-D*D),
       apparent_closure_flag=bool(screen/scale<CONFIG['near_shadow_threshold_rms'] and hidden/scale>CONFIG['far_hidden_threshold_rms']),
       relative_separation_not_absolute_orbit=True))
   # For every anchor, the full state discrepancy must be invariant across the observer sweep.
   for kind,i,j in anchor:
    vals=[r for r in rows if r['label']==label and r['anchor']==kind]
    if np.ptp([r['full_5d_gap'] for r in vals])>1e-13:raise RuntimeError('FIXED_ENDPOINT_FULL_GAP_CHANGED')
    rr=[r for r in vals if r['angle_selection']=='PRE_REGISTERED']
    comparisons.append(dict(label=label,anchor=kind,n_views=len(rr),tau_a=float(t[i]),tau_b=float(t[j]),
      invariant_full_5d_gap=float(np.linalg.norm(q[j]-q[i])),
      screen_gap_range=[min(x['screen_gap_rel'] for x in rr),max(x['screen_gap_rel'] for x in rr)],
      hidden_gap_range=[min(x['hidden_gap_rel'] for x in rr),max(x['hidden_gap_rel'] for x in rr)],
      any_pre_registered_false_closure=any(x['apparent_closure_flag'] for x in rr),
      optimized_screen_gap_rel=min(x['screen_gap_rel'] for x in vals if x['angle_selection']=='OPTIMIZED_AFTER_ENDPOINT_SELECTION') if any(x['angle_selection']=='OPTIMIZED_AFTER_ENDPOINT_SELECTION' for x in vals) else None))
 writecsv('fixed_endpoint_observer_projection.csv',rows)
 writejson('fixed_endpoint_pair_summary.json',comparisons)
 return dict(source_sha256=sha(trace),n_relative_series=len(labels),n_fixed_endpoint_anchors=len(comparisons),n_observer_measurements=len(rows),
   n_pre_registered_false_closures=sum(int(x['apparent_closure_flag']) for x in rows if x['angle_selection']=='PRE_REGISTERED'),
   n_data_optimized_false_closures=sum(int(x['apparent_closure_flag']) for x in rows if x['angle_selection']=='OPTIMIZED_AFTER_ENDPOINT_SELECTION'),
   absolute_orbit_closure_tested=False,endpoint_pair_is_fixed_across_inclinations=True,
   observer_slope_is_a_predeclared_hypothesis_not_a_fit=True),rows,comparisons

def embed10(y,mx):
 """Historical v18.49h kinetic 10D -> (Q,V) embedding, vectorized."""
 y=np.asarray(y,float);l1,v1,t1,w1,l2,v2,t2,w2,z,vz=[y[:,k] for k in range(10)]
 r1=np.exp(l1);r2=np.exp(l2)
 x=np.where(z>=0,1/(1+np.exp(-np.clip(z,0,700))),np.exp(np.clip(z,-700,0))/(1+np.exp(np.clip(z,-700,0))))
 sx=math.sqrt(mx)
 q=np.column_stack([r1*np.cos(t1),r1*np.sin(t1),r2*np.cos(t2),r2*np.sin(t2),sx*x])
 v=np.column_stack([r1*(v1*np.cos(t1)-w1*np.sin(t1)),r1*(v1*np.sin(t1)+w1*np.cos(t1)),
                    r2*(v2*np.cos(t2)-w2*np.sin(t2)),r2*(v2*np.sin(t2)+w2*np.cos(t2)),sx*x*(1-x)*vz])
 jr=r1*r1*w1-r2*r2*w2
 return q,v,jr

def load_historical_kinetic():
 spec=importlib.util.spec_from_file_location('kinetic_metric_v18_49h',str(REF/'kinetic_metric_v18_49h.py'))
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 return module

def absolute_kinetic_hinge(native,ob):
 if sha(native)!=CONFIG['v18_52l_10d_expected_sha256']:raise RuntimeError('10D_NATIVE_SHA_MISMATCH')
 hist=load_historical_kinetic()
 e,n,h,_,_=axes(CONFIG['observer_slope_hypothesis'],ob['pca_hidden'],0)
 rows=[];series=[];model_checks=[]
 with np.load(native,allow_pickle=False) as z:
  keys=sorted(k for k in z.files if k.endswith('__Y_tau'))
  if len(keys)!=20:raise RuntimeError('UNEXPECTED_NATIVE_10D_SERIES')
  for k in keys:
   y=np.asarray(z[k][::int(CONFIG['native_stride'])],float)
   if y.shape[1]!=11 or len(y)<100 or np.min(np.diff(y[:,10]))<=0:raise RuntimeError('BAD_NATIVE_STATES')
   q,v,jrot=embed10(y,CONFIG['mx_kinetic']);tau=y[:,10]
   hmom=signed_wedge(q,v,n,h)
   kdelta=np.ptp(jrot)
   for idx in np.linspace(0,len(y)-1,7,dtype=int):
    H,J=hist.embedding_and_jac(y[idx,:10],CONFIG['mx_kinetic'])
    model_checks.append(max(float(np.max(abs(q[idx]-H[:5]))),float(np.max(abs(v[idx]-H[5:])))))
   vel=np.gradient(hmom,tau,edge_order=2)
   series.append(dict(series=k,n_samples=len(tau),
      J_rot_charge_mean=float(np.mean(jrot)),J_rot_charge_range=float(kdelta),
      J_hinge_kinetic_min=float(hmom.min()),J_hinge_kinetic_max=float(hmom.max()),
      J_hinge_kinetic_range=float(np.ptp(hmom)),
      J_hinge_kinetic_initial=float(hmom[0]),J_hinge_kinetic_final=float(hmom[-1]),
      J_hinge_kinetic_endpoint_change=float(hmom[-1]-hmom[0]),
      finite_difference_hinge_torque_abs_max=float(np.max(abs(vel))),
      kinetic_momentum_assumption='P=V from declared flat 5D T; full driven Lagrangian not authenticated'))
   for idx in np.linspace(0,len(tau)-1,145,dtype=int):
    rows.append(dict(series=k,tau=float(tau[idx]),J_rot_charge=float(jrot[idx]),
        J_hinge_kinetic=float(hmom[idx]),Q_dot_n=float(q[idx]@n),Q_dot_h=float(q[idx]@h),
        V_dot_n=float(v[idx]@n),V_dot_h=float(v[idx]@h),
        derivative_of_hinge_kinetic_fd=float(vel[idx])))
 if max(model_checks)>3e-13:raise RuntimeError('HISTORICAL_KINETIC_EMBEDDING_MISMATCH')
 writecsv('native_absolute_kinetic_hinge_moment.csv',rows)
 writejson('native_absolute_kinetic_hinge_summary.json',series)
 return dict(source_sha256=sha(native),n_native_series=len(series),
       max_kinetic_embedding_error=max(model_checks),
       max_J_rot_charge_range=max(x['J_rot_charge_range'] for x in series),
       max_J_hinge_kinetic_range=max(x['J_hinge_kinetic_range'] for x in series),
       reference_original=next(x for x in series if x['series']=='original__reference__Y_tau'),
       kinetic_candidate_evaluated=True,
       full_Q550_canonical_generator_authenticated=False,
       Noether_hinge_conservation_certificate=False,
       full_10d_integrations_performed=0),rows,series

def plots(shadows,native):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return
 if shadows:
  label='original_negative_1.0e-04'; f=[r for r in shadows if r['label']==label and r['anchor']=='v18_53d_+0deg' and r['angle_selection']=='PRE_REGISTERED']
  f.sort(key=lambda a:a['angle_deg'])
  if f:
   fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
   ang=[r['angle_deg'] for r in f]
   ax[0].plot(ang,[r['screen_gap_rel'] for r in f],'.-',label='2D shadow');ax[0].plot(ang,[r['hidden_gap_rel'] for r in f],'.-',label='Unseen 5D')
   ax[0].plot(ang,[r['full_gap_rel'] for r in f],':',label='Full 5D (same endpoints)')
   ax[0].set(xlabel='Observer inclination (degrees)',ylabel='Gap / RMS relative dQ');ax[0].legend();ax[0].grid(alpha=.3)
   ax[1].plot(ang,[r['signed_depth'] for r in f],'.-',label='Signed depth');ax[1].plot(ang,[r['other_hidden_gap'] for r in f],'.-',label='Other hidden 5D')
   ax[1].set(xlabel='Same observer inclination (degrees)',ylabel='Displacement (dQ units)');ax[1].legend();ax[1].grid(alpha=.3)
   fig.suptitle('Q550 authentic relative separation: SAME time endpoints at every observer angle')
   fig.savefig(OUT/'fixed_endpoint_shadow_hidden_observer.png',dpi=140);plt.close(fig)
 if native:
  f=[r for r in native if r['series']=='original__reference__Y_tau'];f.sort(key=lambda a:a['tau'])
  if f:
   fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
   ax[0].plot([r['tau'] for r in f],[r['J_hinge_kinetic'] for r in f],label='Kinetic hinge candidate J_nh')
   ax[0].axhline(0,ls=':',color='gray');ax[0].set(xlabel='Material tau',ylabel='Kinetic hinge moment');ax[0].grid(alpha=.3);ax[0].legend()
   ax[1].plot([r['tau'] for r in f],[r['J_rot_charge'] for r in f],label='Conserved J1-J2')
   ax[1].set(xlabel='Material tau',ylabel='Rotor counterrotation charge');ax[1].grid(alpha=.3);ax[1].legend()
   fig.suptitle('Q550 authentic ABSOLUTE 5D: hinge candidate varies, counterrotation is conserved')
   fig.savefig(OUT/'absolute_kinetic_hinge_vs_conserved_rotor_charge.png',dpi=140);plt.close(fig)

def locate(base,which):
 name=('extended_horizon_5d_vectors_and_traces.npz' if which=='5d' else 'extended_horizon_10d_states.npz')
 b=Path(base)
 return next((x for x in [b/'q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'/'results'/name,
  b/'q550_v18_52l_review'/name,b/'q550_v18_52p_validation'/name,b/name] if x.is_file()),None)

def main():
 pa=argparse.ArgumentParser();pa.add_argument('--base',default=str(ROOT.parent));pa.add_argument('--trace-npz',default='');pa.add_argument('--native-npz',default='');pa.add_argument('--require-trace',action='store_true');pa.add_argument('--require-native',action='store_true')
 args=pa.parse_args()
 try:
  parent,ob=lineage()
  p5=Path(args.trace_npz) if args.trace_npz else locate(args.base,'5d')
  p10=Path(args.native_npz) if args.native_npz else locate(args.base,'10d')
  relative,obs_rows,pairs=(fixed_endpoint(p5,ob) if p5 and p5.is_file() else (None,[],[]))
  native,native_rows,series=(absolute_kinetic_hinge(p10,ob) if p10 and p10.is_file() else (None,[],[]))
  plots(obs_rows,native_rows)
  missing=[]
  if relative is None:missing.append('5D relative trace')
  if native is None:missing.append('10D native states')
  if not missing:
   decision='AUTHENTIC_FIXED_ENDPOINT_OBSERVER_GEOMETRY_AND_ABSOLUTE_KINETIC_HINGE__NO_FULL_NOETHER_HINGE_CHARGE'
  else:
   decision='PARTIAL_AUTHENTIC_GENEALOGY_ONLY__MISSING_'+','.join(missing)
  summary=dict(version='v18.53e',scientific_decision=decision,
   v18_53d_parent_status=parent['scientific_decision'],relative_5d=relative,absolute_10d_kinetic=native,
   fixed_endpoints_across_observers=True,hinge_geometry_observer_hypothesis_not_physical_symmetry=True,
   canonical_hinge_generator_full_Q550='NOT_CERTIFIED',
   kinetic_flat_Lagrangian_hinge_candidate='EVALUATED_FROM_AUTHENTIC_NATIVE_Q_AND_V' if native else 'NOT_EVALUATED',
   rotor_counterrotation_charge='J1-J2 IS CONSERVED',
   separate_dual_kappa_fold_and_shared_kappa_traces=True,
   global_orbit_closure_certificate=False,relative_orbit_closure_certificate=False,
   missing_inputs=missing)
  writejson('report.json',summary)
  writejson('SCIENTIFIC_REPORT_SERVER.json',dict(version='v18.53e',decision=decision,relative_5d=relative,absolute_10d_kinetic=native))
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text(f'''# Q550 v18.53e — Scientific report

**Decision:** `{decision}`.

## Fixed endpoint observers

{json.dumps(relative,indent=2) if relative else '5D trace unavailable'}

For each of eight authenticated 5D perturbation-minus-reference traces, three independently anchored time pairs are held FIXED while observing at five predeclared inclinations. The sixth viewpoint, where present, is optimized analytically only AFTER the endpoints were selected and cannot serve as independent physical evidence. Distance in complete 5D is unchanged. Hidden distance comprises signed viewing depth and orthogonal residual components. Velocity mismatches are retained. This is **not** absolute orbit closure.

## Absolute native kinetic hinge generator

{json.dumps(native,indent=2) if native else 'Native 10D states unavailable'}

The historical v18.49h exact flat kinetic embedding supplies absolute (Q,V) and kinetic Legendre P=V. The observer-plane SO(2) moment candidate is J_nh=(Q.n)(V.h)-(Q.h)(V.n). Unlike J1-J2, this candidate varies. It is a valid canonical moment of the declared flat kinetic sector, but a conserved Noether hinge charge of the **full driven Q550 action** is NOT demonstrated. The observer PCA plane is not a proven dynamical symmetry axis.

## Scope

No new Q550 integrations, no claim of true helical absolute orbit, no interval proof of a radius bound, no original-vs-autonomous dual-kappa 5D trajectory comparison. The native long 5D traces were generated with a shared kappa.
''')
  print('SCIENTIFIC_DECISION',decision)
  if (args.require_trace and relative is None) or (args.require_native and native is None):return 2
  return 0
 except Exception as exc:
  writejson('report.json',dict(version='v18.53e',scientific_decision='AUDIT_INCOMPLETE',error=str(exc),traceback=traceback.format_exc()))
  print('AUDIT_INCOMPLETE',repr(exc));return 3
if __name__=='__main__':raise SystemExit(main())
