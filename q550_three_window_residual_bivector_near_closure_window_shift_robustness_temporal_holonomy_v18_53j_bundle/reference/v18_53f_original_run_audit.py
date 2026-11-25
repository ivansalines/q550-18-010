#!/usr/bin/env python3
"""Q550 v18.53f — source-exact hinge torque and observer/hidden advance.

Scientific boundary: 10D native states require the authenticated v18.51s
external phase to evaluate the physical RHS. No guesswork from tau (NPZ col 10).
The full model has a conserved counterrotation charge, NOT a conserved hinge
moment. No local Q550 IVPs are run; only parent-source context staging may run.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,math,sys,tempfile,traceback,zipfile
from pathlib import Path
import numpy as np
from torque_math import (axes,embed,force_components,material_rhs,acceleration_from_rhs,
 wedge,inplane_basis,tilted_plane,observer_transport,integral,unit)
ROOT=Path(__file__).resolve().parents[1];REF=ROOT/'reference';OUT=ROOT/'results';BASE=ROOT.parent
C=json.loads((ROOT/'inputs_config.json').read_text());P=C['params']
S51=C['legacy_parent_v1851s_name'];R51=C['legacy_nested_v1851r_name']
L52='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
SEEDS=('original','autonomous')

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for data in iter(lambda:f.read(1<<20),b''):h.update(data)
 return h.hexdigest()
def put(name,obj):
 OUT.mkdir(parents=True,exist_ok=True);(OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False,ensure_ascii=False)+'\n')
def csvput(name,rows):
 if not rows:return
 OUT.mkdir(parents=True,exist_ok=True)
 with (OUT/name).open('w',newline='') as f:
  fields=list(dict.fromkeys(k for r in rows for k in r));w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def mod(path,name):
 spec=importlib.util.spec_from_file_location(name,str(path));m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def preflight():
 parent=REF/'v18_53e_authenticated_report.json'; f=REF/'v18_53e_fixed_endpoint_pair_summary.json'
 mani=REF/'v18_53e_OUTPUT_MANIFEST.json'
 for path,key in ((parent,'v18_53e_authenticated_sha256'),(f,'v18_53e_fixed_endpoints_sha256'),(mani,'v18_53e_output_manifest_sha256'),
  (REF/'v18_52n_authenticated_report.json','v18_52n_observer_sha256'),(REF/'kinetic_metric_v18_49h.py','historical_kinetic_sha256')):
  if sha(path)!=C[key]:raise RuntimeError('REFERENCE_SHA256_MISMATCH: '+str(path))
 m=json.loads(mani.read_text()); old=json.loads(parent.read_text())
 if not old['scientific_decision'].startswith('AUTHENTIC_FIXED_ENDPOINT_OBSERVER_GEOMETRY_AND_ABSOLUTE_KINETIC_HINGE'):
  raise RuntimeError('PREVIOUS_53E_NOT_AUTHENTICATED')
 for local,source in (('v18_53e_authenticated_report.json','results/report.json'),
   ('v18_53e_fixed_endpoint_pair_summary.json','results/fixed_endpoint_pair_summary.json')):
  if sha(REF/local)!=m['files'][source]:raise RuntimeError('53E_OUTPUT_NOT_LINEAGE_VERIFIED: '+local)
 return old,json.loads((REF/'v18_52n_authenticated_report.json').read_text())

def find_npz(which,base):
 name=('extended_horizon_10d_states.npz' if which=='10d' else 'extended_horizon_5d_vectors_and_traces.npz')
 return next((p for p in [base/L52/'results'/name,base/'q550_v18_52l_review'/name,
  base/'q550_v18_52p_validation'/name,base/name] if p.exists()),None)

def find_parent(given):
 return next((p for p in ([given] if given else [])+[BASE/(S51+'_server_sendback.zip'),
   Path('/root/elab1/hysteresis')/(S51+'_server_sendback.zip'),Path('/root/Downloads')/(S51+'_server_sendback.zip')]
  if p and Path(p).is_file()),None)

def stage_original(parent):
 if sha(parent)!=C['parent_v1851s_sha256']:raise RuntimeError('PARENT_V1851S_SHA256_MISMATCH')
 with zipfile.ZipFile(parent) as z:
  if z.testzip() is not None:raise RuntimeError('PARENT_BAD_ZIP')
  rblob=z.read(S51+'/inputs/'+R51+'_server_sendback.zip')
  script=z.read(S51+'/scripts/run_audit.py')
 with zipfile.ZipFile(io.BytesIO(rblob)) as z:
  if z.testzip() is not None:raise RuntimeError('NESTED_PARENT_BAD_ZIP')
  states=z.read(R51+'/results/extended_paired_kappa_nonlinear_states.npz')
 with np.load(io.BytesIO(states)) as z:
  original=float(z['kappa_original']);autonomous=float(z['kappa_autonomous'])
 if abs(.5*(original+autonomous)-C['fixed_kappa'])>1e-12:raise RuntimeError('NATIVE_KAPPA_SHARED_GENEALOGY_MISMATCH')
 loader=mod(REF/'parent_v18_51t_run_audit.py','q550_v1853f_original_stage_loader')
 td=tempfile.TemporaryDirectory(prefix='q550_v18_53f_original_')
 staged=loader.stage_engine(Path(td.name),rblob,script)
 _,_,seed_api,observer,eng,unit,drpt,creport,_=staged[4]
 mismatch={k:float(getattr(eng.params,k))-v for k,v in P.items()}
 if max(abs(z) for z in mismatch.values())>1e-12:
  put('parameter_disagreement.json',mismatch);raise RuntimeError('AUTHENTICATED_MODEL_PARAMETER_MISMATCH')
 ctx={seed:seed_api.make_context(observer,seed,eng,unit,drpt,creport) for seed in SEEDS}
 return td,eng,ctx,dict(parent_sha256=sha(parent),source_omega=float(eng.omega),model_parameters=P,
  historical_seed_kappas={'original':original,'autonomous':autonomous},native_fixed_kappa=C['fixed_kappa'],
  external_phase_source='authenticated v18.51s make_context/theta + cubic_clock_sign')

def moments_tau(t,Q,V,n,h):
 J=wedge(Q,V,n,h);return J,np.gradient(J,t,edge_order=2)

def native_audit(path,engine=None,contexts=None):
 if sha(path)!=C['native_10d_sha256']:raise RuntimeError('AUTHENTIC_NATIVE_10D_SHA256_MISMATCH')
 observer=json.loads((REF/'v18_52n_authenticated_report.json').read_text());e,n,h=axes(C['observer_slope'],observer['pca_hidden'])
 hist=mod(REF/'kinetic_metric_v18_49h.py','q550_v1853f_historical_metric')
 all_rows=[];summaries=[];full_rows=[];source_probes=[]
 with np.load(path,allow_pickle=False) as f:
  xi=np.asarray(f['common_xi'],float);labels=sorted(k for k in f.files if k.endswith('__Y_tau'))
  if len(labels)!=C['sampled_series_expected']:raise RuntimeError('UNEXPECTED_NATIVE_SERIES_COUNT')
  for label in labels:
   a=np.asarray(f[label][::C['native_stride']],float);tt=a[:,10];ee=xi[::C['native_stride']]
   if a.ndim!=2 or a.shape[1]!=11 or np.min(np.diff(tt))<=0:raise RuntimeError('BAD_NATIVE_CHART')
   Q,V,Jchg,g=embed(a[:,:10],P['m_x']);jfix,jfd=moments_tau(tt,Q,V,n,h)
   seed=label.split('_',1)[0]
   if seed not in SEEDS:raise RuntimeError('BAD_SEED')
   # Independent verify of historical exact embedding; full 5D absolute state.
   errors=[]
   for j in np.linspace(0,len(a)-1,8,dtype=int):
    H,jac=hist.embedding_and_jac(a[j,:10],P['m_x'])
    errors.append(float(np.max(np.abs(np.r_[Q[j],V[j]]-H))))
   if max(errors)>2e-13:raise RuntimeError('HISTORICAL_5D_EMBEDDING_MISMATCH')
   T=max(float(tt[-1]-tt[0]),1e-15);u=(tt-tt[0])/T
   alpha=np.deg2rad(C['observer_in_plane_amplitude_degrees'])*np.sin(2*np.pi*C['observer_transport_cycles']*u)
   alphadot=np.deg2rad(C['observer_in_plane_amplitude_degrees'])*2*np.pi*C['observer_transport_cycles']/T*np.cos(2*np.pi*C['observer_transport_cycles']*u)
   beta=np.deg2rad(C['observer_out_of_plane_amplitude_degrees'])*np.sin(2*np.pi*C['observer_transport_cycles']*u)
   betadot=np.deg2rad(C['observer_out_of_plane_amplitude_degrees'])*2*np.pi*C['observer_transport_cycles']/T*np.cos(2*np.pi*C['observer_transport_cycles']*u)
   na=np.cos(alpha)[:,None]*n+np.sin(alpha)[:,None]*h;ha=-np.sin(alpha)[:,None]*n+np.cos(alpha)[:,None]*h
   nad=alphadot[:,None]*ha;had=-alphadot[:,None]*na
   nm=np.cos(beta)[:,None]*n+np.sin(beta)[:,None]*e;hm=np.broadcast_to(h,nm.shape)
   nmd=betadot[:,None]*(-np.sin(beta)[:,None]*n+np.cos(beta)[:,None]*e);hmd=np.zeros_like(hm)
   j_in=wedge(Q,V,na,ha);j_m=wedge(Q,V,nm,hm)
   obs_in=observer_transport(Q,V,na,ha,nad,had);obs_m=observer_transport(Q,V,nm,hm,nmd,hmd)
   if np.max(np.abs(j_in-jfix))>5e-13 or np.max(abs(obs_in))>5e-13:raise RuntimeError('IN_PLANE_INVARIANCE_FAILURE')
   row=dict(series=label,n_samples=len(a),charge_range=float(np.ptp(Jchg)),
    static_hinge_moment_range=float(np.ptp(jfix)),in_plane_rotation_J_difference_max=float(np.max(abs(j_in-jfix))),
    in_plane_observer_transport_max=float(np.max(abs(obs_in))),moving_plane_J_range=float(np.ptp(j_m)),
    moving_plane_observer_transport_max=float(np.max(abs(obs_m))),
    finite_difference_static_torque_max=float(np.max(abs(jfd))),historical_embedding_max_abs=max(errors),
    physical_RHS_source_used=bool(engine is not None),source_field_max_relative=None,source_force_max_abs=None,
    static_torque_derivative_max_error=None,static_torque_integral_max_error=None,
    moving_plane_torque_integral_max_error=None,
    full_Q550_Noether_hinge_charge_certified=False)
   if np.ptp(Jchg)>C['counterrotation_drift_tolerance']:raise RuntimeError('COUNTERROTATION_CHARGE_DRIFT')
   if engine is not None:
    ctx=contexts[seed];phase=np.asarray([ctx['theta'](ctx['event']+float(xx)) for xx in ee],float)
    design=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
    sign=np.asarray([float(engine.clock.cubic_clock_sign(ctx['event']+float(xx),design)) for xx in ee])
    Q1,V1,J1,g1,Delta,parts=force_components(a[:,:10],phase,C['fixed_kappa'],P)
    FF_model=material_rhs(a[:,:10],phase,C['fixed_kappa'],P)
    Afr=acceleration_from_rhs(a[:,:10],FF_model,P['m_x'])
    A=sum(parts.values());ferr=float(np.max(np.abs(Afr-A)))
    source_rel=0.;source_force=0.;rate_err=0.
    for j in range(len(a)):
     Fxi=np.asarray(engine.physical.rhs_from_original(np.r_[a[j,:10],phase[j]],
       C['fixed_kappa'],engine.params,engine.slope,sign[j],engine.support,engine.audit),float)
     if Fxi.shape!=(11,):raise RuntimeError('PARENT_RHS_BAD_SHAPE')
     Fsrc=Fxi[:10]*g1['x'][j]/engine.omega
     source_rel=max(source_rel,float(np.max(abs(Fsrc-FF_model[j])))/max(1.,float(np.max(abs(Fsrc)))))
     source_force=max(source_force,float(np.max(abs(acceleration_from_rhs(a[j,:10],Fsrc,P['m_x'])[0]-A[j]))))
     rate_err=max(rate_err,abs(float(Fxi[10])-float(engine.omega*sign[j])))
    if source_rel>C['native_source_field_rel_tolerance'] or source_force>C['source_force_abs_tolerance'] or rate_err>1e-10:
     raise RuntimeError(f'PARENT_RHS_NOT_MATCHED {label}: {source_rel}, {source_force}, {rate_err}')
    ft={name:wedge(Q,force,n,h) for name,force in parts.items()}
    tau0=sum(ft.values());taum=wedge(Q,A,nm,hm);totalm=taum+obs_m
    errfd=float(np.max(abs(jfd-tau0)))
    integ0=integral(tt,tau0);integm=integral(tt,totalm)
    errint=float(np.max(abs(jfix-jfix[0]-integ0)))
    errintm=float(np.max(abs(j_m-j_m[0]-integm)))
    # Analytical in-plane torque under rotation must reproduce physical torque.
    torque_in=wedge(Q,A,na,ha)
    if np.max(abs(torque_in-tau0))>5e-13:raise RuntimeError('IN_PLANE_TORQUE_INVARIANCE_FAILURE')
    row.update(source_field_max_relative=source_rel,source_force_max_abs=max(source_force,ferr),
      signed_external_drive_source_max_abs=rate_err,static_torque_derivative_max_error=errfd,
      static_torque_integral_max_error=errint,moving_plane_torque_integral_max_error=errintm,
      static_total_torque_change=float(integ0[-1]),moving_total_torque_change=float(integm[-1]),
      moving_physical_torque_integral=float(integral(tt,taum)[-1]),
      moving_observer_transport_integral=float(integral(tt,obs_m)[-1]),
      static_endpoint_moment_change=float(jfix[-1]-jfix[0]),moving_endpoint_moment_change=float(j_m[-1]-j_m[0]),
      force_torque_components_integrated={name:float(integral(tt,z)[-1]) for name,z in ft.items()})
    source_probes.append(dict(series=label,phase_source='authenticated_v1851s',n_samples=len(a),
      original_rhs_field_max_relative=source_rel,physical_acceleration_max_error=max(source_force,ferr),
      original_signed_clock_max_error=rate_err))
    if max(errint,errintm)>C['kinematic_derivative_abs_tolerance']:
     raise RuntimeError(f'INTEGRATED_TORQUE_BALANCE_FAILURE {label}: {errint}, {errintm}')
   for j in np.linspace(0,len(tt)-1,201,dtype=int):
    item=dict(series=label,tau=float(tt[j]),J_counterrotation=float(Jchg[j]),
      J_hinge_static=float(jfix[j]),J_hinge_inplane=float(j_in[j]),J_hinge_tilted_plane=float(j_m[j]),
      apparent_observer_transport_torque=float(obs_m[j]),finite_difference_J_hinge_dot=float(jfd[j]),
      beta_deg=float(np.rad2deg(beta[j])))
    if engine is not None:
     item.update(torque_radial_potential=float(ft['radial_potential'][j]),
       torque_phase_coupling=float(ft['phase_coupling'][j]),torque_radial_dissipation=float(ft['radial_dissipation'][j]),
       torque_x_barrier=float(ft['x_barrier'][j]),torque_x_activity=float(ft['x_activity'][j]),
       torque_x_dissipation=float(ft['x_dissipation'][j]),torque_total_model=float(tau0[j]),
       torque_moving_physical=float(taum[j]),torque_moving_total=float(totalm[j]))
    all_rows.append(item)
   summaries.append(row)
 if len(summaries)!=C['sampled_series_expected']:raise RuntimeError('NATIVE_SERIES_MISSING')
 csvput('native_hinge_torque_and_observer_transport.csv',all_rows)
 csvput('authenticated_rhs_verification.csv',source_probes)
 put('native_absolute_hinge_torque_summary.json',summaries)
 return dict(n_native_series=len(summaries),max_charge_drift=max(z['charge_range'] for z in summaries),
  max_static_hinge_excursion=max(z['static_hinge_moment_range'] for z in summaries),
  max_in_plane_rotation_spurious_torque=max(z['in_plane_observer_transport_max'] for z in summaries),
  max_moving_plane_transport_torque=max(z['moving_plane_observer_transport_max'] for z in summaries),
  original_source_authenticated=bool(engine is not None),
  max_original_RHS_relative_residual=(max(z['source_field_max_relative'] for z in summaries) if engine is not None else None),
  max_torque_ledger_integral_residual=(max(z['static_torque_integral_max_error'] for z in summaries) if engine is not None else None),
  max_moving_plane_ledger_integral_residual=(max(z['moving_plane_torque_integral_max_error'] for z in summaries) if engine is not None else None),
  kinetic_plane_momentum_not_Noether_charge=True),all_rows,summaries

def sample_index(t,t0):
 j=int(np.argmin(abs(t-t0)))
 if abs(float(t[j])-float(t0))>max(1e-9,2e-7*float(np.median(np.diff(t)))):
  raise RuntimeError(f'FIXED_ENDPOINT_NOT_FOUND {t0}')
 return j

def relative_audit(path):
 if sha(path)!=C['relative_5d_sha256']:raise RuntimeError('AUTHENTIC_RELATIVE_5D_SHA256_MISMATCH')
 observer=json.loads((REF/'v18_52n_authenticated_report.json').read_text());e,n,h=axes(C['observer_slope'],observer['pca_hidden'])
 pairs=json.loads((REF/'v18_53e_fixed_endpoint_pair_summary.json').read_text());rows=[];summaries=[]
 with np.load(path,allow_pickle=False) as z:
  labels=sorted({k[:-len('__dQ')] for k in z.files if k.endswith('__dQ')})
  if len(labels)!=C['relative_series_expected']:raise RuntimeError('BAD_RELATIVE_SERIES_COUNT')
  cache={}
  for label in labels:
   tau=np.asarray(z[label+'__tau'][::C['relative_stride']],float)
   dq=np.asarray(z[label+'__dQ'][::C['relative_stride']],float);dv=np.asarray(z[label+'__dV'][::C['relative_stride']],float)
   if dq.shape!=(len(tau),5) or dv.shape!=dq.shape:raise RuntimeError('BAD_RELATIVE_ARRAYS')
   cache[label]=(tau,dq,dv)
  for pair in pairs:
   label=pair['label'];t,q,v=cache[label];i=sample_index(t,float(pair['tau_a']));j=sample_index(t,float(pair['tau_b']))
   if j<=i:raise RuntimeError('NON_FORWARD_PAIR')
   Dq=q[j]-q[i];Dq_int=np.trapezoid(v[i:j+1],t[i:j+1],axis=0)
   eps=float(np.linalg.norm(Dq-Dq_int));scale=max(1e-13,float(np.linalg.norm(Dq)))
   re=eps/scale
   if re>C['relative_endpoint_integration_rel_tolerance']:raise RuntimeError(f'RELATIVE_Q_AND_V_INCONSISTENT: {label} {re}')
   T=float(t[-1]-t[0]);omega=2*np.pi*C['observer_transport_cycles']/T
   a=np.deg2rad(C['observer_in_plane_amplitude_degrees'])*np.sin(omega*(t[i:j+1]-t[0]))
   ad=np.deg2rad(C['observer_in_plane_amplitude_degrees'])*omega*np.cos(omega*(t[i:j+1]-t[0]))
   depth=(-np.sin(a)[:,None]*n+np.cos(a)[:,None]*h)
   w=(np.cos(a)[:,None]*n+np.sin(a)[:,None]*h)
   velocity_depth=np.einsum('ij,ij->i',v[i:j+1],depth)
   moving_depth=-ad*np.einsum('ij,ij->i',q[i:j+1],w)
   actual=float(q[j]@depth[-1]-q[i]@depth[0]);phys=float(np.trapezoid(velocity_depth,t[i:j+1]));obs=float(np.trapezoid(moving_depth,t[i:j+1]))
   mob_error=abs(actual-(phys+obs))
   one=dict(label=label,anchor=pair['anchor'],tau_a=float(t[i]),tau_b=float(t[j]),
    displacement_full_5d=float(np.linalg.norm(Dq)),velocity_integral_5d=float(np.linalg.norm(Dq_int)),
    velocity_5d_displacement_absolute_mismatch=eps,relative_integral_mismatch=re,
    moving_observer_hidden_coordinate_advance=actual,moving_observer_physical_hidden_integral=phys,
    moving_observer_frame_transport_integral=obs,moving_observer_transport_balance_abs_error=mob_error,
    data_type='RELATIVE_KICK_MINUS_REFERENCE_NOT_ABSOLUTE_ORBIT')
   for angle in C['observer_tilts_degrees']:
    rad=math.radians(angle);wfix=math.cos(rad)*n+math.sin(rad)*h;d=-math.sin(rad)*n+math.cos(rad)*h
    screen=np.hypot(float(Dq@e),float(Dq@wfix));hidden=math.sqrt(max(0,float(Dq@Dq)-screen**2))
    advance=float(Dq@d);advance_integrated=float(Dq_int@d)
    rows.append(dict(label=label,anchor=pair['anchor'],angle_deg=angle,tau_a=float(t[i]),tau_b=float(t[j]),
      full_5d_advance=float(np.linalg.norm(Dq)),visible_2d_advance=screen,hidden_5d_advance=hidden,
      signed_hidden_plane_advance=advance,integrated_signed_hidden_velocity=advance_integrated,
      signed_hidden_ledger_error=float(advance-advance_integrated),
      fixed_endpoints=True,relative_separation_only=True))
   summaries.append(one)
 if len(summaries)!=C['fixed_endpoint_pairs_expected']:raise RuntimeError('NOT_24_FIXED_ENDPOINTS')
 csvput('fixed_endpoint_hidden_advance_force_vs_observer.csv',rows)
 csvput('relative_moving_observer_hidden_ledger.csv',summaries)
 return dict(n_fixed_endpoint_pairs=len(summaries),n_fixed_angle_measurements=len(rows),
  max_relative_velocity_integral_mismatch=max(z['relative_integral_mismatch'] for z in summaries),
  max_moving_frame_hidden_ledger_error=max(z['moving_observer_transport_balance_abs_error'] for z in summaries),
  relative_motion_not_absolute_orbit=True,angles_predeclared=True,optimized_post_selection_excluded=True),rows,summaries

def plots(rows,relative,summary):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return
 target=[z for z in rows if z['series']=='original__reference__Y_tau']
 if target:
  t=np.array([a['tau'] for a in target]); fig,ax=plt.subplots(3,1,figsize=(11,10),sharex=True,layout='constrained')
  ax[0].plot(t,[a['J_hinge_static'] for a in target],label='Kinetic hinge moment: fixed plane')
  ax[0].plot(t,[a['J_hinge_tilted_plane'] for a in target],label='Moving plane',alpha=.8)
  ax[0].plot(t,[a['J_hinge_inplane'] for a in target],':',label='Rotated axes in same plane')
  ax[1].plot(t,[a['finite_difference_J_hinge_dot'] for a in target],label='dJ/dtau (numerical)')
  if 'torque_total_model' in target[0]:
   ax[1].plot(t,[a['torque_total_model'] for a in target],':',label='Sum of genuine force torques')
   for key in ['torque_phase_coupling','torque_radial_potential','torque_radial_dissipation']:
    ax[2].plot(t,[a[key] for a in target],label=key)
  else:
   ax[2].plot(t,[a['apparent_observer_transport_torque'] for a in target],label='Plane transport (kinematic)')
  for a in ax:a.legend(fontsize=8);a.grid(alpha=.2)
  ax[-1].set_xlabel('Material time tau');fig.suptitle('Q550 v18.53f: fixed-plane physical torque versus moving observer')
  fig.savefig(OUT/'original_hinge_torque_observer_decomposition.png',dpi=140);plt.close(fig)
 if relative:
  r=[z for z in relative if z['label']=='original_negative_1.0e-04' and z['anchor']=='v18_53d_+0deg']
  if r:
   r.sort(key=lambda z:z['angle_deg']);angle=[z['angle_deg'] for z in r]
   fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
   ax.plot(angle,[z['visible_2d_advance'] for z in r],'.-',label='2D visible endpoint gap')
   ax.plot(angle,[z['hidden_5d_advance'] for z in r],'.-',label='Hidden 5D endpoint gap')
   ax.plot(angle,[z['full_5d_advance'] for z in r],':',label='Full 5D (same endpoints)')
   ax.set(xlabel='Predeclared observer angle [deg]',ylabel='Relative displacement',title='Q550 authentic relative separation: observer decomposition')
   ax.legend();ax.grid(alpha=.2);fig.savefig(OUT/'fixed_endpoint_hidden_advance.png',dpi=140);plt.close(fig)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,default=BASE)
 ap.add_argument('--native-npz',type=Path);ap.add_argument('--trace-npz',type=Path);ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--require-native',action='store_true');ap.add_argument('--require-trace',action='store_true');ap.add_argument('--require-parent',action='store_true')
 args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
 try:
  parent_report,observer=preflight()
  npath=args.native_npz or find_npz('10d',args.base)
  tpath=args.trace_npz or find_npz('5d',args.base)
  ppath=find_parent(args.parent_v1851s)
  report=dict(version=C['version'],scientific_decision='NOT_YET_CLASSIFIED',
   source_manifest='SOURCE_MANIFEST.json',previous_53e_authenticated_sha256=C['v18_53e_authenticated_sha256'],
   native_status='MISSING',relative_status='MISSING',parent_status='MISSING',
   new_Q550_orbit_integrations=0,no_global_Noether_hinge_charge_certificate=True,
   torque_source_decomposition_status='NOT_AUTHENTICATED',
   external_drive_torque_at_frozen_state='Clock sign enters phase rate; first 10 RHS components are not directly signed',
   moving_observer_is_observer_choice_not_physical_hinge_transport=True)
  temp=None;eng=None;ctx=None
  try:
   if ppath:
    temp,eng,ctx,prov=stage_original(ppath)
    report['parent_status']='SHA256_AUTHENTICATED';report['original_engine_provenance']=prov
    put('authenticated_original_source_metadata.json',prov)
   if npath and npath.is_file():
    n,nrows,ns=native_audit(npath,eng,ctx);report['native_status']='SHA256_AUTHENTICATED';report['absolute_hinge']=n
    report['torque_source_decomposition_status']='AUTHENTICATED' if eng is not None else 'PARTIAL_NATIVE_FINITE_DIFFERENCE_ONLY'
   else:nrows=[]
   if tpath and tpath.is_file():
    rel,rr,prs=relative_audit(tpath);report['relative_status']='SHA256_AUTHENTICATED';report['relative_hidden_advance']=rel
   else:rr=[]
   if not npath or not tpath or not ppath:
    report['scientific_decision']='PARTIAL_INPUTS_NO_FULL_FORCE_TORQUE_CERTIFICATE'
   else:
    report['scientific_decision']='AUTHENTIC_ORIGINAL_HINGE_TORQUE_AND_OBSERVER_TRANSPORT_LEDGER__HINGE_NOT_NOETHER_CHARGE'
   put('report.json',report)
   plots(nrows,rr,report)
   msg=['# Q550 v18.53f — Hinge-Momentum Torque Decomposition, Observer Transport & Hidden-Advance Attribution','',
    '## Scope','Fixed plane kinetic hinge moment, moving-plane observer correction, authenticated archived 10D source force decomposition, and pre-fixed 5D relative separation endpoint transport.',
    '',f"**Classification:** `{report['scientific_decision']}`",'','## Input status',
    f"- Original source: {report['parent_status']}",f"- Native 10D: {report['native_status']}",f"- Relative 5D: {report['relative_status']}",
    '','## Scientific boundaries','The six torque channels are literal terms of the documented physical material-time RHS, not uniquely separable forms of energy transfer. A clock-sign switch at fixed physical state and phase does not directly change the first 10 components. Changing clock history changes phase and therefore later coupling torque.',
    'The 5D hinge momentum is canonical for the DECLARED flat kinetic sector P=V, not a proven full-system Noether charge.',
    'A rotating basis WITHIN the same plane must contribute zero geometric torque. Tilting the physical observer plane in the surrounding 5D space adds a frame-transport term; this is a choice of observer, not physical transport of a uniquely defined hinge.',
    'The 5D endpoint experiments use kicked minus reference separation, not absolute orbit closure.',
    'Do not identify one torque channel with a fundamental physical interaction without further derivation.',
    '','## Results',json.dumps({k:report.get(k) for k in ('absolute_hinge','relative_hidden_advance')},indent=2)]
   (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(msg)+'\n')
   if args.require_native and report['native_status']=='MISSING':raise RuntimeError('REQUIRE_NATIVE_MISSING')
   if args.require_trace and report['relative_status']=='MISSING':raise RuntimeError('REQUIRE_TRACE_MISSING')
   if args.require_parent and report['parent_status']=='MISSING':raise RuntimeError('REQUIRE_PARENT_MISSING')
   print('V18_53F_STATUS',report['scientific_decision'])
  finally:
   if temp is not None:temp.cleanup()
 except Exception:
  (OUT/'audit_failure.txt').write_text(traceback.format_exc())
  put('failure_summary.json',{'version':C['version'],'scientific_decision':'FAILED_OR_MISSING_REQUIRED_SOURCE','traceback_file':'audit_failure.txt'})
  raise

if __name__=='__main__':main()
