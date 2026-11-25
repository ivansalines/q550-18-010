#!/usr/bin/env python3
"""Q550 v18.52p: physical canonical-hinge *candidate* and collapse-barrier audit.

Numerical postprocessing only. Do not infer a Noether charge from coordinate rotation.
No original Q550 Hamiltonian or conjugate-momentum law is invented. The native
kinetic chart gives a justified P=V LEGENDRE KINETIC PROXY, not a certified
canonical action for the driven 10D engine. All source records authenticated.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
PARENT='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
BASE=ROOT.parent/PARENT/'results'
VECTORS='extended_horizon_5d_vectors_and_traces.npz'
STATES='extended_horizon_10d_states.npz'
SHA5='d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195'
SHA10='92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12'
KIN_SHA='29bef2c91899a9da46bc8c9159b89c70362f505e0adead3f814c487a9a17bb1a'
SEEDS=('original','autonomous')
SIGNS=('negative','positive')
AMPS=('1.0e-04','3.0e-04')
LABELS=[f'{s}_{p}_{a}' for s in SEEDS for p in SIGNS for a in AMPS]
M=-.308
ALPHA=17.04


def digest(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def saved_json(path,item):
 Path(path).write_text(json.dumps(item,indent=2,sort_keys=True,allow_nan=False)+'\n')

def saved_csv(path,rows,cols):
 with Path(path).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)

def norms(v):return np.linalg.norm(v,axis=-1)

def stats(z):
 z=np.asarray(z,float);z=z[np.isfinite(z)]
 if len(z)==0:return {'n':0}
 return dict(n=int(len(z)),minimum=float(z.min()),maximum=float(z.max()),median=float(np.median(z)),mean=float(z.mean()),p01=float(np.quantile(z,.01)),p99=float(np.quantile(z,.99)))

def cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]

def canonical_rotation(angle_rad):
 c,s=math.cos(angle_rad),math.sin(angle_rad)
 return np.array([[c,s],[-s,c]],dtype=float)

def passive_nodal_basis(m=M):
 """The FIRST ROW is (1,m)/sqrt(1+m^2); theta=atan(m), *negative*."""
 return canonical_rotation(math.atan(m))

def canonical_block(R):
 S=np.zeros((4,4));S[:2,:2]=R;S[2:,2:]=R;return S

def canonical_identity(q,p,R):
 # Genuine symplectic coordinate change, independent of dynamics.
 qr=q@R.T;pr=p@R.T
 return dict(J_before=cross(q,p),J_after=cross(qr,pr),norm_before=norms(p),norm_after=norms(pr))

def count_sign_flips(v,zero_floor=0.):
 v=np.asarray(v,float);sg=np.sign(v);sg[abs(v)<=zero_floor]=0;sg=sg[sg!=0]
 return int(np.count_nonzero(sg[1:]!=sg[:-1])) if len(sg)>1 else 0

def bivector(q,v):
 return np.stack([q[:,i]*v[:,j]-q[:,j]*v[:,i] for i in range(5) for j in range(i+1,5)],axis=1)

def derivative(y,t):return np.gradient(y,t,edge_order=2,axis=0)

def safe_residual(a,b,scale):
 return float(np.max(abs(a-b))/max(1e-30,float(scale)))

def pca_hidden_from_o():
 rep=json.loads((ROOT/'reference/PARENT_v18_52o_report.json').read_text())
 h=np.asarray(rep['pca_hidden'],float)
 if abs(float(np.linalg.norm(h))-1)>1e-10:raise RuntimeError('v18.52o historical PCA axis not unit')
 return h

def hinge_gate(q,h,m=M,tilt=ALPHA):
 R=passive_nodal_basis(m)
 t=q[:,:2]@R[1]
 z=q[:,2:]@h
 a=math.radians(tilt)
 return math.cos(a)*t+math.sin(a)*z,math.cos(a)*t-math.sin(a)*z

def roots(t,f):
 indices=np.flatnonzero((f[:-1]<0)&(f[1:]>=0) | (f[:-1]>0)&(f[1:]<=0))
 out=[]
 for i in indices:
  f0,f1=f[i],f[i+1]
  u=float(-f0/(f1-f0));r=float(t[i]+u*(t[i+1]-t[i]))
  if not out or r-out[-1][0]>1e-9:out.append((r,1 if f1>f0 else -1))
 return out

def source_preflight(f5,f10,require_states):
 parent=json.loads((ROOT/'reference/PARENT_v18_52l_OUTPUT_MANIFEST.json').read_text())
 prior=json.loads((ROOT/'reference/PARENT_v18_52o_report.json').read_text())
 if parent['file_sha256']['results/'+VECTORS]!=SHA5 or parent['file_sha256']['results/'+STATES]!=SHA10:
  raise ValueError('Pinned v18.52l ancestry mismatch')
 if prior['version']!='v18.52o' or prior['stage']!='AUTHENTIC_FULL_RESOLUTION_POSTPROCESS':
  raise ValueError('v18.52o predecessor is not genuine full-resolution audit')
 if digest(ROOT/'reference/kinetic_metric_v18_49h.py')!=KIN_SHA:
  raise ValueError('Historical kinetic embedding file changed')
 if not f5.exists() or digest(f5)!=SHA5:raise FileNotFoundError('Authentic v18.52l physical 5D data missing or SHA failure: '+str(f5))
 state_ok=f10.exists() and digest(f10)==SHA10
 if require_states and not state_ok:raise FileNotFoundError('Authentic v18.52l native 10D states missing or SHA failure: '+str(f10))
 return {'v18_52l_vectors_sha256':SHA5,'v18_52l_10d_states_sha256':SHA10 if state_ok else None,
         'historical_kinetic_embedding_sha256':KIN_SHA,'ancestry_v18_52o_version':prior['version'],
         'native_10d_states_available':bool(state_ok)}

def geometry(q,v,a,t):
 r5=norms(q);p=q[:,:2];vp=v[:,:2];ap=a[:,:2]
 r2=norms(p);j=cross(p,vp);tau=cross(p,ap)
 floor=max(1e-18,float(np.median(r2))*1e-11)
 if np.min(r2)<=floor:raise RuntimeError('Planar projected radius near zero: continuous phase / barrier not identifiable')
 vnorm=norms(vp)
 # Cauchy identity in transverse plane: ||v||² = rdot² + (J/r)².
 rdot=np.sum(p*vp,axis=1)/r2
 T=.5*np.sum(vp*vp,axis=1)
 Trot=.5*(j/r2)**2
 Trr=.5*rdot**2
 radforce=np.sum(p*ap,axis=1)/r2
 centrifugal=j*j/(r2**3)
 rddot=centrifugal+radforce
 # NOTE: Jdot is torque, for P=V kinetic proxy, no H or Noether claim.
 dJ=derivative(j,t)
 abs_floor=max(1e-15,1e-7*float(np.max(abs(tau))))
 torque_consistency=safe_residual(dJ[2:-2],tau[2:-2],max(float(np.max(abs(tau))),abs_floor))
 qv=np.sum(q*v,axis=1);biv=bivector(q,v);biv_a=bivector(q,a)
 bnorm=norms(biv)
 bividentity=safe_residual(bnorm**2,r5**2*np.sum(v*v,axis=1)-qv*qv,max(float(np.max(bnorm**2)),1e-28))
 # Continuous unwinding exists on open segment; not topological integer.
 phase=np.unwrap(np.arctan2(p[:,1],p[:,0]));turns=float((phase[-1]-phase[0])/(2*math.pi))
 eps=1e-35
 rel_J_span=float((j.max()-j.min())/max(float(np.median(abs(j))),eps))
 jmin=float(np.min(abs(j)));vmax=float(np.max(vnorm))
 # A posteriori bound conditioned on observed j_min and finite observed speed upper bound.
 bound=jmin/max(vmax,eps)
 barrier_ok=bool(jmin>1e-14 and count_sign_flips(j)==0)
 info={'n_samples':len(t),'tau_start':float(t[0]),'tau_end':float(t[-1]),
   'r2_stats':stats(r2),'r5_stats':stats(r5),'relative_kinetic_J_stats':stats(j),
   'relative_J_relative_range_to_median':rel_J_span,'relative_J_sign_flips':count_sign_flips(j),
   'torque_cross_q_a_stats':stats(tau),'dJ_dt_vs_torque_relative_max':torque_consistency,
   'two_dimensional_kinetic_energy_decomposition_max_abs':float(np.max(abs(T-Trot-Trr))),
   'five_dimensional_bivector_squared_identity_relative_max':bividentity,
   'bivector_norm_stats':stats(bnorm),'bivector_derivative_norm_stats':stats(norms(biv_a)),
   'turning_angle_open_path_turns':turns,'turning_angle_integer_topological_certificate':False,
   'centrifugal_kinematic_term_stats':stats(centrifugal),
   'radial_acceleration_from_q_v_a_stats':stats(rddot),
   'inward_radial_component_stats':stats(radforce),
   'centrifugal_exceeds_inward_force_fraction':float(np.mean(centrifugal>np.maximum(-radforce,0))),
   'kinematic_barrier_energy_stats':stats(Trot),'plane_kinetic_energy_stats':stats(T),
   'j_min_observed':jmin,'v2_max_observed':vmax,
   'conditional_sampled_interval_lower_bound_jmin_over_vmax':bound,
   'measured_r2_min_over_conditional_bound':float(r2.min()/bound) if bound>0 else None,
   'no_J_zero_seen_on_sampled_grid':barrier_ok,
   'certified_noether_charge':False,'certified_future_nocollapse':False}
 signals=dict(t=t,r2=r2,r5=r5,j=j,torque=tau,dJ=dJ,bnorm=bnorm,radforce=radforce,
  centrifugal=centrifugal,rddot=rddot,Trot=Trot,T=T,turns=(phase-phase[0])/(2*math.pi))
 return info,signals

def extract_qv_polar(y):
 """Native archived y=(l1,v1,theta1,w1,l2,v2,theta2,w2,z,vz)."""
 if y.ndim!=2 or y.shape[1]<10:raise ValueError('Not native 10D polar chart')
 qq=np.empty((len(y),4));vv=np.empty_like(qq)
 for off,k in ((0,0),(4,2)):
  l,vr,th,omega=(y[:,off+j] for j in range(4))
  r=np.exp(l);c=np.cos(th);s=np.sin(th)
  qq[:,k]=r*c;qq[:,k+1]=r*s
  vv[:,k]=r*(vr*c-omega*s);vv[:,k+1]=r*(vr*s+omega*c)
 return qq,vv

def abs_axis(y,t,key):
 q,v=extract_qv_polar(y)
 j1=cross(q[:,:2],v[:,:2]);j2=cross(q[:,2:4],v[:,2:4])
 polar1=np.exp(2*y[:,0])*y[:,3];polar2=np.exp(2*y[:,4])*y[:,7]
 if max(abs(j1-polar1))>1e-10 or max(abs(j2-polar2))>1e-10:raise RuntimeError('Kinetic polar/angular momentum identity broke')
 jac={'J_first':j1,'J_second':j2}
 r={}
 for nm,j in jac.items():
  dJ=derivative(j,t)
  a=derivative(v[:,:2] if nm=='J_first' else v[:,2:4],t)
  qp=q[:,:2] if nm=='J_first' else q[:,2:4]
  torque=cross(qp,a)
  scale=max(float(np.max(abs(dJ))),float(np.max(abs(torque))),1e-25)
  r[nm]=dict(stats=stats(j),relative_range_to_median=float((j.max()-j.min())/max(float(np.median(abs(j))),1e-25)),
     sign_flips=count_sign_flips(j),derivative_stats=stats(dJ),torque_stats=stats(torque),
     derivative_torque_relative_max=float(np.max(abs(dJ[2:-2]-torque[2:-2]))/scale))
 return r,dict(t=t,J_first=j1,J_second=j2)

def pair_absolute_delta(ref,kick,trel):
 """Different solution clocks! PCHIP onto the actual common MATERIAL tau only."""
 vr=(trel>=max(ref['t'][0],kick['t'][0]))&(trel<=min(ref['t'][-1],kick['t'][-1]))
 if np.count_nonzero(vr)<100:raise RuntimeError('Native 10D reference and kick have no time overlap')
 t=trel[vr]
 out={}
 for plane in ('J_first','J_second'):
  Jr=PchipInterpolator(ref['t'],ref[plane])(t)
  Jk=PchipInterpolator(kick['t'],kick[plane])(t)
  out[plane]=dict(t=t,delta=Jk-Jr,ref=Jr,kick=Jk)
 return vr,out

def process(args):
 require_states=not args.allow_missing_10d
 ancestry=source_preflight(args.source5,args.source10,require_states)
 stage='FULL_RESOLUTION_NATIVE_10D_AND_RELATIVE_5D_KINEMATIC_AUDIT' if args.sample_stride==1 and ancestry['native_10d_states_available'] else 'DEVELOPMENT_DIAGNOSTIC_ONLY'
 if args.mode=='preflight':
  saved_json(ROOT/'PREFLIGHT.json',dict(version='v18.52p',stage=stage,source_authentication=ancestry,labels=LABELS));print('PREFLIGHT PASS',stage);return
 out=ROOT/'results';out.mkdir(exist_ok=True)
 h=pca_hidden_from_o();R=passive_nodal_basis(M);theta=math.degrees(math.atan(M));
 # Strict canonical angle and rotation checks do NOT certify any dynamical Noether symmetry.
 S=canonical_block(R);I=np.eye(2);O=np.block([[np.zeros((2,2)),I],[-I,np.zeros((2,2))]])
 symp_err=float(np.max(abs(S.T@O@S-O)))
 report={'version':'v18.52p','title':'Canonical Hinge, Noether-Torque Candidate & Conditional Noncollapse Barrier Audit',
  'stage':stage,'source_authentication':ancestry,'sample_stride':args.sample_stride,
  'passive_nodal_axis_angle_deg':theta,'observed_visual_tilt_deg':ALPHA,
  'm_slope':M,'R_nodal':R.tolist(),
  'positive_angle_user_matrix_aligns_slope_minus_0p308':False,
  'canonical_symplectic_block_error':symp_err,
  'canonical_statement':'S=diag(R,R) is a symplectic coordinate change; J invariant under R algebraically, not necessarily under dynamics.',
  'historical_kinetic_chart':'Q=(r1 cos theta1,r1 sin theta1,r2 cos theta2,r2 sin theta2,sqrt(mx)*x); 2D kinetic conjugate proxy P1:2=V1:2.',
  'independent_full_hamiltonian_and_legendre_action_verified':False,
  'noether_conservation_certified':False,'topological_charge_certified':False,'global_no_collapse_certified':False,
  'physical_absolute_trajectories':{},'physical_relative_trajectories':{},
  'warning':'Q_delta and V_delta are DIFFERENCES of two flows. q cross delta-P is not J_kick-J_ref; a coordinate transformation cannot establish a dynamical Noether symmetry.'}
 if symp_err>1e-13:raise AssertionError('Canonical matrix not symplectic')
 # Verify rotation invariance on exact same relative q,v for all trajectory times.
 emax=0.;time_rows=[];gates=[];per_branch={};all_signal={}
 with np.load(args.source5,allow_pickle=False) as f:
  for label in LABELS:
   pref=label+'__';t=f[pref+'tau'][::args.sample_stride].astype(float)
   q=f[pref+'dQ'][::args.sample_stride].astype(float)
   v=f[pref+'dV'][::args.sample_stride].astype(float)
   a=f[pref+'da'][::args.sample_stride].astype(float)
   if len(t)<100 or q.shape!=(len(t),5) or a.shape!=q.shape or not np.all(np.diff(t)>0):
    raise ValueError('Invalid authenticated 5D physical trajectory: '+label)
   rr=canonical_identity(q[:,:2],v[:,:2],R)
   emax=max(emax,float(np.max(abs(rr['J_after']-rr['J_before']))))
   geom,signals=geometry(q,v,a,t)
   wp,wm=hinge_gate(q,h)
   ge=[]
   for view,w in (('+',wp),('-',wm)):
    for idx,(event_tau,direction) in enumerate(roots(t,w)):
     obj=dict(label=label,view=view,event_index=idx,tau=event_tau,orientation=direction,
      J_kinetic=float(np.interp(event_tau,t,signals['j'])),
      r2=float(np.interp(event_tau,t,signals['r2'])),
      r5=float(np.interp(event_tau,t,signals['r5'])),
      radial_second_derivative=float(np.interp(event_tau,t,signals['rddot'])),
      torque=float(np.interp(event_tau,t,signals['torque'])))
     gates.append(obj);ge.append(obj)
   geom['n_observer_gate_events']=len(ge)
   report['physical_relative_trajectories'][label]=geom
   per_branch[label]=(t,geom,signals)
   if label=='original_negative_1.0e-04':all_signal=signals
   for i in range(0,len(t),args.output_stride):
    time_rows.append(dict(label=label,tau=float(t[i]),r2=float(signals['r2'][i]),r5=float(signals['r5'][i]),
        J_kinetic=float(signals['j'][i]),torque=float(signals['torque'][i]),
        J_time_derivative=float(signals['dJ'][i]),bivector_norm_5d=float(signals['bnorm'][i]),
        radial_second_derivative=float(signals['rddot'][i]),centrifugal_term=float(signals['centrifugal'][i]),
        radial_force_term=float(signals['radforce'][i]),
        kinetic_rotation_barrier_energy=float(signals['Trot'][i])))
 report['max_J_invariance_under_nodal_coordinate_rotation_abs']=emax
 # Native 10D full-state kinetic proxy, optionally disabled only for local DEVELOPMENT.
 delta_rows=[];absplot={}
 if ancestry['native_10d_states_available']:
  with np.load(args.source10,allow_pickle=False) as states:
   refs={}
   for seed in SEEDS:
    key=seed+'__reference__Y_tau'
    Y=states[key][::args.sample_stride].astype(float)
    t=Y[:,10]
    if len(t)<100 or np.any(np.diff(t)<=0):raise ValueError('Invalid native reference material clock '+key)
    geom,js=abs_axis(Y[:,:10],t,key)
    report['physical_absolute_trajectories'][key]=geom
    refs[seed]=js
   for label in LABELS:
    key=label+'__kick__Y_tau'
    Y=states[key][::args.sample_stride].astype(float)
    t=Y[:,10]
    if len(t)<100 or np.any(np.diff(t)<=0):raise ValueError('Invalid native kick material clock '+key)
    geom,js=abs_axis(Y[:,:10],t,key)
    report['physical_absolute_trajectories'][key]=geom
    seed=label.split('_')[0]
    tr=per_branch[label][0]
    mask,delta=pair_absolute_delta(refs[seed],js,tr)
    jrel=per_branch[label][2]['j'][mask]
    d=delta['J_first']['delta']
    deriv=dict(valid_count=int(mask.sum()),norm_difference_between_absolute_delta_J_and_relative_J=stats(d-jrel),
      delta_absolute_J_stats=stats(d),relative_J_at_matched_times_stats=stats(jrel),
      comparison_not_identity=True)
    report['physical_relative_trajectories'][label]['absolute_Jkick_minus_Jreference_comparison']=deriv
    for i in range(0,len(d),max(1,args.output_stride)):
     delta_rows.append(dict(label=label,tau=float(delta['J_first']['t'][i]),
        Jkick_minus_Jreference=float(d[i]),Jrelative_cross_deltaQ_deltaV=float(jrel[i]),
        difference=float(d[i]-jrel[i])))
    if label=='original_negative_1.0e-04':absplot={'ref':refs[seed],'kick':js,'label':label,'delta':delta['J_first']}
 report['global_observation']={
  'all_relative_J_sign_flips':int(sum(x['relative_J_sign_flips'] for x in report['physical_relative_trajectories'].values())),
  'all_relative_J_nonconstant':bool(all(x['relative_J_relative_range_to_median']>1e-3 for x in report['physical_relative_trajectories'].values())),
  'n_relative_branches':len(report['physical_relative_trajectories']),
  'n_absolute_kinetic_trajectories':len(report['physical_absolute_trajectories']),
  'n_observer_gate_events':len(gates),
  'observed_finite_horizon_positive_radius':bool(all(x['r5_stats']['minimum']>0 for x in report['physical_relative_trajectories'].values())),
  'noether_hinge_classification':'NOT_CERTIFIED_FROM_DRIVEN_ENGINE_WITHOUT_ACTION_SYMMETRY_TEST',
  'topological_hinge_classification':'OPEN_PATH_WINDING_NOT_INTEGER_TOPOLOGICAL_CERTIFICATE',
  'collapse_barrier_classification':'FINITE_HORIZON_NONZERO_J_AND_SPEED_OBSERVED; GLOBAL NO-COLLAPSE NOT CERTIFIED'}
 saved_json(out/'report.json',report)
 saved_csv(out/'relative_5d_hinge_kinematics_decimated.csv',time_rows,list(time_rows[0]))
 saved_csv(out/'observer_gate_kinetic_torque_events.csv',gates,list(gates[0]) if gates else ['label','view','event_index','tau','orientation','J_kinetic','r2','r5','radial_second_derivative','torque'])
 saved_csv(out/'native_absolute_vs_relative_J_matched_time.csv',delta_rows,list(delta_rows[0]) if delta_rows else ['label','tau','Jkick_minus_Jreference','Jrelative_cross_deltaQ_deltaV','difference'])
 plot(all_signal,absplot,out)
 lines=[
 '# Q550 v18.52p — Scientific server report','',f'**Stage:** `{stage}`','',
 '## Exact canonical coordinate identity',
 f'Fixed observation-based nodal slope `m={M}`; passive nodal angle `theta={theta:.9f} deg` (negative). Visual hinge tilt is a DIFFERENT angle, `±{ALPHA} deg`.',
 f'Symplectic block matrix error `{symp_err:.3e}`; max coordinate-transform J error `{emax:.3e}`. These are algebraic identities, NOT tests of dynamical conservation.','',
 '## Finite-horizon physical relative kinetic observations',
 '| Physical branch | J min | J max | sign flips | relative drift range/median | min r2 | min r5 | open path turns |',
 '|---|---:|---:|---:|---:|---:|---:|---:|']
 for label,z in report['physical_relative_trajectories'].items():
  j=z['relative_kinetic_J_stats'];lines.append(f"| {label} | {j['minimum']:.6g} | {j['maximum']:.6g} | {z['relative_J_sign_flips']} | {z['relative_J_relative_range_to_median']:.4g} | {z['r2_stats']['minimum']:.6g} | {z['r5_stats']['minimum']:.6g} | {z['turning_angle_open_path_turns']:.6g} |")
 lines.extend(['','## Absolute canonical kinetic proxy versus relative two-flow quantity',
  f"Native original 10D states available: `{ancestry['native_10d_states_available']}`.",
  'From exact historical kinetic embedding, the first polar pair gives `J1_kinetic=r1^2*w1`, and second pair `J2_kinetic=r2^2*w2`. Both are reconstructed and numerically differentiated along EACH physical material clock. The original action / true canonical momentum in presence of forcing is not supplied, so the proxy is NOT a certified Noether charge.',
  '| Native physical trajectory | J1 min | J1 max | J1 relative range | J2 sign flips |',
  '|---|---:|---:|---:|---:|'])
 for key,quantity in report['physical_absolute_trajectories'].items():
  j1=quantity['J_first'];j2=quantity['J_second']
  lines.append(f"| {key} | {j1['stats']['minimum']:.6g} | {j1['stats']['maximum']:.6g} | {j1['relative_range_to_median']:.5g} | {j2['sign_flips']} |")
 lines.extend(['',
  'The quantity `J_kick-J_reference` is not the same as `(Q_kick-Q_reference) x (V_kick-V_reference)`; an explicit matched-time comparison is written to CSV.','',
  '## Noncollapse: necessary conditions rather than a false theorem',
  'In the fixed physical 2D plane, `T2 = 1/2 (dr/dtau)^2 + Jkin^2/(2r^2)` and `rddot = Jkin^2/r^3 + (q_plane·a_plane)/r` are exact kinematic identities. Observed `J_min>0` and finite sampled speed yield a posteriori `r >= min|J|/max|v2|` on the sampled domain only; no global speed bound, exact Noether conservation or continuous-time zero-exclusion proof is present.',
  'Open-path unwrapped angular turns are not integral topological charges.','',
  '## Scientific disposition',
  '**Noether charge:** not certified. **Topological charge:** not certified. **Global impossibility of collapse:** not certified. A physically meaningful finite-horizon geometric obstruction may still be present and is measured here.',''])
 (out/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines))
 print('Q550 v18.52p',stage,'relative branches',len(per_branch),'absolute physical states',len(report['physical_absolute_trajectories']),
       'gates',len(gates),'Noether NOT certified; collapse NOT certified',flush=True)

def plot(s,absplot,out):
 if not s:return
 t=s['t'];fig,ax=plt.subplots(3,1,figsize=(12.5,10),sharex=True)
 ax[0].plot(t,s['j'],label='relative J_kinetic = q1 v2 - q2 v1')
 ax[0].plot(t,s['dJ'],alpha=.72,label='numerical dJ/dtau')
 ax[0].plot(t,s['torque'],ls='--',alpha=.7,label='q1 a2 - q2 a1')
 ax[0].axhline(0,color='.5',lw=.7);ax[0].legend(fontsize=8,ncol=3);ax[0].set_ylabel('Angular quantities (mixed units)')
 ax[1].plot(t,s['r2'],label='r2 of physical relative Q1,Q2');ax[1].plot(t,s['r5'],label='full physical r5',alpha=.7)
 ax[1].set_ylabel('Physical separation radius');ax[1].legend(fontsize=8)
 ax[2].plot(t,s['centrifugal'],label='positive J²/r³ kinematic term')
 ax[2].plot(t,s['radforce'],label='(q·a)/r radial acceleration')
 ax[2].plot(t,s['rddot'],label='sum: radial acceleration',alpha=.75)
 ax[2].axhline(0,color='.5',lw=.7);ax[2].legend(fontsize=8);ax[2].set_ylabel('Radial acceleration');ax[2].set_xlabel('Material time τ')
 fig.suptitle('Q550 v18.52p — relative 5D physical hinge: torque and radial balance (NOT Noether certificate)')
 fig.tight_layout();fig.savefig(out/'relative_hinge_torque_and_collapse_balance.png',dpi=150);plt.close(fig)
 if absplot:
  ref=absplot['ref'];kick=absplot['kick'];dd=absplot['delta']
  fig,axs=plt.subplots(2,1,figsize=(12,7),sharex=False)
  axs[0].plot(ref['t'],ref['J_first'],label='reference absolute kinetic J1')
  axs[0].plot(kick['t'],kick['J_first'],label='kick absolute kinetic J1',alpha=.7)
  axs[0].set_ylabel('Absolute planar kinetic J');axs[0].legend(fontsize=8)
  axs[1].plot(dd['t'],dd['delta'],label='Jkick - Jreference')
  axs[1].set(xlabel='Physical common τ',ylabel='Difference of absolute J');axs[1].legend(fontsize=8)
  fig.suptitle('Q550 v18.52p — absolute native 10D physical kinetic momentum proxies')
  fig.tight_layout();fig.savefig(out/'native_absolute_kinetic_moment_comparison.png',dpi=150);plt.close(fig)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--mode',choices=('preflight','audit'),default='audit')
 p.add_argument('--source5',type=Path,default=BASE/VECTORS)
 p.add_argument('--source10',type=Path,default=BASE/STATES)
 p.add_argument('--sample-stride',type=int,default=1)
 p.add_argument('--output-stride',type=int,default=32)
 p.add_argument('--allow-missing-10d',action='store_true',help='DEVELOPMENT ONLY; server run must supply genuine 10D states')
 args=p.parse_args()
 if args.sample_stride<1 or args.output_stride<1:raise ValueError('Invalid strides')
 if args.allow_missing_10d and args.sample_stride==1:raise ValueError('Missing 10D is allowed only for development stride >1')
 process(args)

if __name__=='__main__':main()
