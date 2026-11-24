#!/usr/bin/env python3
"""Q550 v18.53c: dual-kappa fold locus and sampled minimum-radius envelope.

The historical RHS is the arbiter only when the SHA-pinned ancestor is staged.
The analytic fold curve is a numerical family, NOT a global/interval theorem.
The 'minimum radius' is min(r1,r2) of a critical relative equilibrium,
not the minimum radius of arbitrary dynamical trajectories.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,sys,traceback
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.optimize import linear_sum_assignment

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'
C=json.loads((ROOT/'inputs_config.json').read_text());F=C['dual_fold']


def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for data in iter(lambda:f.read(1<<20),b''):h.update(data)
 return h.hexdigest()


def load_module(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec)
 sys.modules[name]=mod;spec.loader.exec_module(mod)
 return mod

# The actual v18.53b code is included, not silently recreated.
B=load_module(ROOT/'reference/v18_53b_run_audit.py','authentic_ancestor_v1853b_reconstructed_equations')


def output(name,obj):
 OUT.mkdir(exist_ok=True,parents=True);p=OUT/name
 p.write_text(obj if isinstance(obj,str) else json.dumps(obj,indent=2,allow_nan=False,ensure_ascii=False)+'\n')


def write_csv(name,rows):
 if not rows:return
 fields=list(dict.fromkeys(k for row in rows for k in row))
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def core_jacobian(a,J,k,s):
 """Analytic d(F1,F2,F3)/d(log r1, log r2, omega1) for delta=pi."""
 l1,l2,w1=map(float,a)
 r1,r2=math.exp(l1),math.exp(l2)
 w2=(r1*r1*w1-J)/(r2*r2)
 e=B.P.e_reserve; x=.5
 ratio1=r2/r1;ratio2=r1/r2;ratio_sq=r1*r1/(r2*r2)
 dw2_dl1=2*ratio_sq*w1
 dw2_dl2=-2*w2
 dw2_dw1=ratio_sq
 return np.array([
  [2.4*r1-1.6*r1*r1-4*e*e/r1**4-k*x*ratio1, k*x*ratio1,2*w1],
  [2*w2*dw2_dl1+k*x*ratio2,2*w2*dw2_dl2+2.4*r2-1.6*r2*r2-4*e*e/r2**4-k*x*ratio2,2*w2*dw2_dw1],
  [dw2_dl1,dw2_dl2,1+dw2_dw1]],dtype=float)


def augmented(z,J,k):
 a=np.asarray(z[:3],float);s=float(z[3]);f=B.radial_residual(a,J,k,s)
 if not np.all(np.isfinite(f)) or max(np.abs(f))>1e5:return np.array([1e6]*4,float)
 A=core_jacobian(a,J,k,s)
 sv=np.linalg.svd(A,compute_uv=False)
 den=max(1e-18,float(sv[0]*sv[1]))
 return np.r_[f,np.linalg.det(A)/den]


def fold_metrics(z,J,k,omega,eng=None):
 a=np.asarray(z[:3],float);s=float(z[3]);A=core_jacobian(a,J,k,s)
 U,sv,Vh=np.linalg.svd(A);w=U[:,-1];v=Vh[-1,:]
 Fs=np.array([0.,0.,-.5]);Fk=np.array([.5*math.exp(a[1]-a[0]),.5*math.exp(a[0]-a[1]),0.])
 trans=float(w@Fs);kappa_trans=float(w@Fk)
 hlist=[1e-3,3e-4,1e-4]
 quad=[float(w@((B.radial_residual(a+h*v,J,k,s)+B.radial_residual(a-h*v,J,k,s)-2*B.radial_residual(a,J,k,s))/h**2)) for h in hlist]
 qa=B.correct_and_audit(a,J,k,s,omega,eng)
 eigs=np.array([complex(e['re'],e['im']) for e in qa['eigenvalues']]);ix=np.argsort(abs(eigs))
 e0=eigs[ix[0]];others=eigs[ix[1:]]
 native=qa['native']
 return {'aug_max':float(np.max(np.abs(augmented(z,J,k)))),
   'radial_max':float(np.max(abs(B.radial_residual(a,J,k,s)))),
   'smallest_singular_3x3':float(sv[-1]),'singular_values_3x3':[float(q) for q in sv],
   'fold_drive_transversality':trans,'fold_kappa_projection':kappa_trans,
   'fold_quadratic_curvature_stencils':quad,
   'smallest_eigenvalue_complex':{'re':float(e0.real),'im':float(e0.imag)},
   'other_seven_max_real':float(max(others.real)),
   'spectral_abscissa':qa['spectral_abscissa'],
   'native_available':native.get('available',False),'native_source_match':native.get('source_match',False),
   'native_rhs_max':native.get('native_equilibrium_rhs_max'),
   'native_field_relative':native.get('source_reconstructed_field_relative'),
   'native_jacobian_relative':native.get('source_reconstructed_jacobian_relative'),
   'native_jacobian_step_relative':native.get('native_jacobian_step_relative'),
   'radial_core_nondegenerate':bool(abs(trans)>1e-5 and min(abs(q) for q in quad)>1e-5 and sv[-1]<5e-7),
   'one_real_mode_near_zero':bool(abs(e0)<float(F['local_fold_mode_tol']) and abs(e0.imag)<float(F['local_fold_mode_tol'])),
   'other_seven_stable':bool(max(others.real)<-1e-4)}


def fold_correct(zseed,J,k,omega,eng=None):
 sol=root(lambda z:augmented(z,J,k),np.asarray(zseed,float),method='hybr',tol=2e-12)
 z=np.asarray(sol.x,float)
 resid=float(np.max(abs(augmented(z,J,k))))
 if not np.all(np.isfinite(z)) or not np.isfinite(resid):raise RuntimeError('NONFINITE_FOLD_CORRECTOR')
 if resid>float(F['fold_augmented_root_tol']):raise RuntimeError('FOLD_CORRECTOR_RESIDUAL '+str(resid))
 if min(math.exp(z[0]),math.exp(z[1]))<float(F['min_fold_radius']):raise RuntimeError('FOLD_RADIUS_CHART_FLOOR')
 diag=fold_metrics(z,J,k,omega,eng)
 if not all(diag[key] for key in ('radial_core_nondegenerate','one_real_mode_near_zero','other_seven_stable')):
  raise RuntimeError('FOLD_LOCAL_CONDITIONS_NOT_MET '+str({k:diag[k] for k in ('smallest_singular_3x3','one_real_mode_near_zero','other_seven_stable')}))
 if eng is not None and not (diag['native_source_match'] and diag['native_rhs_max']<float(F['native_max_rhs_tol'])):
  raise RuntimeError('AUTHENTIC_RHS_OR_JACOBIAN_NOT_MATCHED')
 return z,diag


def aug_jacobian(z,J,k,scale=2e-5):
 A=np.zeros((4,4));z=np.asarray(z,float)
 for j in range(4):
  h=scale*max(1.,abs(z[j]));e=np.zeros(4);e[j]=h
  A[:,j]=(augmented(z+e,J,k)-augmented(z-e,J,k))/(2*h)
 return A


def kappa_sensitivity(z,J,k,h=2e-6):
 # The 4D fold locus is defined by G(z,k)=0; derivative from implicit differentiation.
 Jz=aug_jacobian(z,J,k)
 Gk=(augmented(z,J,k+h)-augmented(z,J,k-h))/(2*h)
 dz=np.linalg.solve(Jz,-Gk)
 return {'ds_dk':float(dz[3]),'dr1_dk':float(math.exp(z[0])*dz[0]),
  'dr2_dk':float(math.exp(z[1])*dz[1]),'domega1_dk':float(dz[2]),
  'four_dimensional_fold_condition':float(np.linalg.cond(Jz)),
  'augmented_linear_residual':float(np.max(abs(Jz@dz+Gk)))}


def row_from_root(k,z,diag,J,source,omega):
 r1=math.exp(z[0]);r2=math.exp(z[1]);s=float(z[3]);w2=(r1*r1*z[2]-J)/(r2*r2)
 deriv=kappa_sensitivity(z,J,k)
 return dict(region=source,kappa=float(k),signed_drive_fold=s,r1_fold=r1,r2_fold=r2,
   r_min_fold=min(r1,r2),r_max_fold=max(r1,r2),min_radius_component=('r1' if r1<r2 else 'r2'),
   omega1_fold=float(z[2]),omega2_fold=float(w2),J=J,fold_equilibrium_augmented_residual=diag['aug_max'],
   radial_residual=diag['radial_max'],smallest_singular_equilibrium=diag['smallest_singular_3x3'],
   zero_mode_re=diag['smallest_eigenvalue_complex']['re'],zero_mode_im=diag['smallest_eigenvalue_complex']['im'],
   other_seven_max_real=diag['other_seven_max_real'],transversality_drive=diag['fold_drive_transversality'],
   quadratic_fold_curvature=diag['fold_quadratic_curvature_stencils'][-1],
   within_signed_drive_minus1_plus1=bool(-1<=s<=1),
   native_checked=diag['native_available'],native_source_match=diag['native_source_match'],
   native_rhs_max=diag['native_rhs_max'],native_jacobian_rel=diag['native_jacobian_relative'],
   **deriv)


def scan_values(values,seed_k,seed_z,J,omega,eng,label):
 rows=[];failures=[];states=[]
 central_value=min(values,key=lambda v:abs(v-seed_k))
 for direction in (-1,1):
  subseq=sorted((float(v) for v in values if (float(v)-seed_k)*direction>1e-13),reverse=(direction<0))
  prev=seed_z.copy();prev_k=seed_k
  for k in subseq:
   try:
    step=k-prev_k;seed=prev.copy()
    # implicit tangent predictor includes all 4 state components; solve for it directly
    jac=aug_jacobian(prev,J,prev_k);hh=2e-6
    gk=(augmented(prev,J,prev_k+hh)-augmented(prev,J,prev_k-hh))/(2*hh)
    seed=prev-step*np.linalg.solve(jac,gk)
    z,diag=fold_correct(seed,J,k,omega,eng)
    if np.linalg.norm(z-prev)>max(0.2,35*abs(step)):
     raise RuntimeError('FOLD_BRANCH_JUMP_GUARD')
    rows.append(row_from_root(k,z,diag,J,label,omega));states.append({'kappa':k,'core_z':z.tolist(),'diag':diag})
    prev,prev_k=z,k
   except Exception as e:
    failures.append({'region':label,'kappa':k,'last_kappa':prev_k,'status':'CONTINUATION_STOP_NOT_FOLD_DISAPPEARANCE','error':str(e)})
    break
 z,diag=fold_correct(seed_z,J,seed_k,omega,eng)
 rows.append(row_from_root(seed_k,z,diag,J,label,omega));states.append({'kappa':seed_k,'core_z':z.tolist(),'diag':diag})
 rows.sort(key=lambda v:v['kappa']);states.sort(key=lambda v:v['kappa'])
 return rows,states,failures


def sampled_envelope(rows,original_k,autonomous_k):
 kmin=min(original_k,autonomous_k);kmax=max(original_k,autonomous_k)
 rr=[r for r in rows if kmin-1e-12<=r['kappa']<=kmax+1e-12]
 if len(rr)<3:return {'status':'INSUFFICIENT_GRID_POINTS'}
 vals=np.array([r['r_min_fold'] for r in rr]);dv=np.array([r['dr2_dk'] if r['min_radius_component']=='r2' else r['dr1_dk'] for r in rr]);kk=np.array([r['kappa'] for r in rr])
 jmin=int(np.argmin(vals));jmax=int(np.argmax(vals))
 signs=np.sign(dv)
 zeros=[{'from':float(kk[i]),'to':float(kk[i+1])} for i in range(len(rr)-1) if signs[i]*signs[i+1]<0]
 return {'status':'SAMPLED_ENVELOPE_NOT_INTERVAL_CERTIFICATE',
  'n_grid':len(rr),'interval':[float(kk[0]),float(kk[-1])],
  'sampled_min_radius':float(vals[jmin]),'sampled_min_location_kappa':float(kk[jmin]),
  'sampled_max_of_min_radius':float(vals[jmax]),'sampled_max_location_kappa':float(kk[jmax]),
  'monotonic_sampled':bool(np.all(dv>0) or np.all(dv<0)),
  'derivative_range':[float(np.min(dv)),float(np.max(dv))],
  'radius_component_switches':int(sum(rr[i]['min_radius_component']!=rr[i+1]['min_radius_component'] for i in range(len(rr)-1))),
  'detected_derivative_zero_brackets':zeros,
  'continuous_extremum_excluded':False,'global_dynamic_radius_bound_certified':False}


def plot(rows,env,orig,auto):
 try:
  import matplotlib;matplotlib.use('Agg')
  import matplotlib.pyplot as plt
 except ImportError:return
 fig,ax=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
 for name,group,color in [('dual',[r for r in rows if r['region']=='dual'],'tab:blue'),('extended',[r for r in rows if r['region']=='extended'],'tab:gray')]:
  if not group:continue
  kk=np.array([r['kappa'] for r in group]);alpha=.85 if name=='dual' else .5
  ax[0,0].plot(kk,[r['r1_fold'] for r in group],label=f'{name} r1',alpha=alpha,color=color,lw=1.2)
  ax[0,0].plot(kk,[r['r2_fold'] for r in group],label=f'{name} r2',alpha=alpha,color=color,lw=1.6,ls='--')
  ax[0,1].plot(kk,[r['signed_drive_fold'] for r in group],alpha=alpha,color=color,lw=1.6,label=name)
  ax[1,0].plot(kk,[r['r_min_fold'] for r in group],alpha=alpha,color=color,lw=1.7,label=name)
  ax[1,1].plot(kk,[r['dr2_dk'] for r in group],alpha=alpha,color=color,lw=1.6,label=name)
 for a in ax.flat:
  for k,c in [(orig,'tab:orange'),(auto,'tab:green')]:a.axvline(k,color=c,lw=.8,ls=':')
  a.grid(alpha=.15)
 ax[0,0].set(title='Radii on the fold locus',xlabel='kappa0',ylabel='r1,r2')
 ax[0,1].set(title='Fold signed-drive location',xlabel='kappa0',ylabel='s_fold')
 ax[1,0].set(title='Critical equilibrium minimum radius (not orbital bound)',xlabel='kappa0',ylabel='min(r1,r2)')
 ax[1,1].set(title='Implicit slope along the fold',xlabel='kappa0',ylabel='dr2_fold / dkappa0')
 for a in ax.flat:a.legend(fontsize=8)
 fig.suptitle('Q550 v18.53c — Dual-kappa fold-radius envelope (constant-drive equilibria)')
 fig.savefig(OUT/'dual_kappa_fold_radius_envelope.png',dpi=150);plt.close(fig)
 # zoom directly compares the two prescribed kappa values
 dual=[q for q in rows if q['region']=='dual']
 fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
 ax.plot([r['kappa'] for r in dual],[r['r_min_fold'] for r in dual],lw=2)
 for k,c,n in [(orig,'tab:orange','original'),(auto,'tab:green','autonomous')]:
  nearest=min(dual,key=lambda r:abs(r['kappa']-k));ax.scatter([nearest['kappa']],[nearest['r_min_fold']],c=c,zorder=3,label=n)
 ax.grid(alpha=.2);ax.legend();ax.set(xlabel='kappa0',ylabel='min(r1_fold,r2_fold)',title='Q550 v18.53c — sampled min/max of the critical minimum radius')
 fig.savefig(OUT/'dual_kappa_minimum_radius_zoom.png',dpi=150);plt.close(fig)


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent-v1851s',type=Path);ap.add_argument('--require-parent',action='store_true');ap.add_argument('--offline-only',action='store_true');args=ap.parse_args()
 OUT.mkdir(exist_ok=True,parents=True)
 if args.offline_only and args.require_parent:raise RuntimeError('CONTRADICTORY_REQUIRE_PARENT_OFFLINE_ONLY')
 for nm,entry in [('v18_53b_authenticated_report.json','v1853b_report_sha256'),('v18_53b_authenticated_fold.json','v1853b_fold_sha256'),('v18_53a_authenticated_report.json','v1853a_report_sha256')]:
  if sha(ROOT/'reference'/nm)!=C[entry]:raise RuntimeError('HISTORICAL_REFERENCE_SHA_MISMATCH '+nm)
 repb=json.loads((ROOT/'reference/v18_53b_authenticated_report.json').read_text())
 fold=json.loads((ROOT/'reference/v18_53b_authenticated_fold.json').read_text())
 if not (repb['source_checked_refined_fold'] and repb['parent_v1851s_verified'] and fold['conditions']['authentic_source_match']):
  raise RuntimeError('V1853B_ANCHOR_NOT_AUTHENTICATED')
 J=float(fold['J']);k0=float(fold['kappa']);z0=np.r_[np.array(fold['core_a'],float),fold['s_fold']]
 if float(np.max(abs(augmented(z0,J,k0))))>1e-7:raise RuntimeError('V1853B_ANCHOR_FOLD_EQUATIONS_MISMATCH')
 ko=float(C['kappa_seeds']['original']);ka=float(C['kappa_seeds']['autonomous'])
 omega=.844357915384206;eng=None;temp=None;source={}
 report={'version':'v18.53c','title':C['title'],'scientific_decision':'OFFLINE_FOLD_LOCUS_SCREENING__ORIGINAL_ENGINE_NOT_VERIFIED',
  'authenticated_v18_53b_anchor':{'J':J,'kappa':k0,'signed_drive_fold':float(z0[3]),'report_sha256':sha(ROOT/'reference/v18_53b_authenticated_report.json'),'fold_sha256':sha(ROOT/'reference/v18_53b_authenticated_fold.json')},
  'kappa_values':{'original':ko,'autonomous':ka,'shared_v18_52l':k0},'parent_v1851s_staged':False,'warnings':[],
  'radius_kind':'RADIUS_OF_CRITICAL_RELATIVE_EQUILIBRIUM__NOT_MINIMUM_OVER_TRAJECTORIES',
  'full_ramp_studied':False,'global_radial_bound_certified':False,'fold_locus_interval_theorem':False}
 try:
  if not args.offline_only:
   stage=load_module(ROOT/'reference/v18_52v_r1_run_audit.py','historical_v1852v_stage_v1853c')
   parent=stage.find_parent(args.parent_v1851s)
   if parent is not None:
    if sha(parent)!=C['v1851s_sha256']:raise RuntimeError('GENEALOGICAL_PARENT_SHA_MISMATCH')
    temp,eng,_,_,_,source=stage.stage_parent(parent,np.array([12.8,12.8001],float))
    for n,value in C['model_params'].items():
     if abs(float(getattr(eng.params,n))-value)>1e-12:raise RuntimeError('AUTHENTIC_PARAMETER_MISMATCH '+n)
    omega=float(eng.omega);report['parent_v1851s_staged']=True
    report['source_provenance']={'parent_sha256':sha(parent),'omega':omega,'original_rhs':source['rhs_from_original'],'original_clock':source['cubic_clock_sign']}
   elif args.require_parent:raise RuntimeError('V1851S_PARENT_REQUIRED_FOR_STRICT_AUDIT')
  # Recompute the authenticated central fold from the analytical equation, then its physical RHS.
  z0,d0=fold_correct(z0,J,k0,omega,eng)
  report['shared_anchor_recheck']=row_from_root(k0,z0,d0,J,'shared_anchor',omega)
  vals=np.unique(np.r_[np.linspace(min(ka,ko),max(ka,ko),int(F['n_grid'])),[ka,ko,k0]])
  dual,dual_states,dfail=scan_values(vals,k0,z0,J,omega,eng,'dual')
  extended=[];extended_states=[];efail=[]
  if F.get('extended_enabled',False):
   ex=np.unique(np.r_[np.linspace(float(F['extended_lo']),float(F['extended_hi']),int(F['extended_n'])),[ka,ko,k0]])
   extended,extended_states,efail=scan_values(ex,k0,z0,J,omega,eng,'extended')
  rows=dual+extended
  write_csv('fold_locus_continuation.csv',rows)
  output('fold_locus_selected_states.json',{'dual':dual_states,'extended':extended_states})
  output('continuation_stops_not_proven_fold_termination.json',dfail+efail)
  env=sampled_envelope(dual,ko,ka)
  # Validate explicit actual kappa roots, not nearest-grid interpolation.
  endpoint={}
  for name,k in [('original',ko),('autonomous',ka),('shared',k0)]:
   z,d=fold_correct(z0,J,k,omega,eng)
   endpoint[name]=row_from_root(k,z,d,J,'anchor_'+name,omega)
  report['endpoints']=endpoint
  # Separate a same-drive comparison from the different-drive fold endpoints.
  # The two fold radii occur at DIFFERENT s_fold; they are not a simultaneous pair.
  anchored=json.loads((ROOT/'reference/v18_53a_authenticated_report.json').read_text())
  fixed_a=np.asarray(anchored['results']['original']['algebraic_branch_root']['core_root'],float)
  same_drive={}
  for name,k in [('original',ko),('autonomous',ka)]:
   good,a,eps=B.solve(fixed_a,J,k,-1.)
   if not good:raise RuntimeError('SAME_DRIVE_ANCHOR_CORRECTOR_FAILED '+name)
   qa=B.correct_and_audit(a,J,k,-1.,omega,eng)
   if eng is not None and (not qa['native']['source_match'] or qa['native']['native_equilibrium_rhs_max']>F['native_max_rhs_tol']):
    raise RuntimeError('SAME_DRIVE_ANCHOR_NATIVE_MISMATCH '+name)
   same_drive[name]={'kappa':k,'signed_drive':-1.,'r1':float(math.exp(a[0])),'r2':float(math.exp(a[1])),
     'min_radius':float(min(math.exp(a[0]),math.exp(a[1]))),
     'root_max_abs':float(np.max(abs(B.radial_residual(a,J,k,-1.)))),
     'native_checked':qa['native']['available'],'native_match':qa['native'].get('source_match',False)}
  report['fixed_common_drive_comparison']=same_drive
  report['sampled_envelope']=env
  report['continuation_counts']={'dual':len(dual),'extended':len(extended),'dual_failures':len(dfail),'extended_failures':len(efail)}
  report['native_checks']={'checked_all':bool(eng is not None and all(r['native_checked'] and r['native_source_match'] and r['native_rhs_max']<F['native_max_rhs_tol'] for r in rows)),
    'checked_endpoints':bool(eng is not None and all(r['native_checked'] and r['native_source_match'] for r in endpoint.values()))}
  report['sampled_envelope']['full_grid_completed']=bool(len(dfail)==0 and len(dual)==len(vals))
  report['sampled_envelope']['within_two_kappa_only']=True
  report['sampled_envelope']['observed_min_at_endpoint']=bool(abs(env['sampled_min_location_kappa']-min(ka,ko))<1e-12 or abs(env['sampled_min_location_kappa']-max(ka,ko))<1e-12)
  report['sampled_envelope']['observed_max_at_endpoint']=bool(abs(env['sampled_max_location_kappa']-min(ka,ko))<1e-12 or abs(env['sampled_max_location_kappa']-max(ka,ko))<1e-12)
  report['cross_seed_r2_fold_difference_original_minus_autonomous']=float(endpoint['original']['r2_fold']-endpoint['autonomous']['r2_fold'])
  report['cross_seed_signed_drive_fold_difference_original_minus_autonomous']=float(endpoint['original']['signed_drive_fold']-endpoint['autonomous']['signed_drive_fold'])
  # Avoid upgrading a sampled sign survey into an interval proof.
  if report['native_checks']['checked_all'] and report['native_checks']['checked_endpoints'] and report['sampled_envelope']['full_grid_completed']:
   report['scientific_decision']='AUTHENTIC_DUAL_KAPPA_FOLD_LOCUS_AND_SAMPLED_CRITICAL_RADIUS_ENVELOPE__NOT_GLOBAL_MINIMUM_OR_INTERVAL_PROOF'
  elif eng is not None:report['scientific_decision']='AUTHENTIC_PARENT_PRESENT_BUT_FOLD_LOCUS_SOURCE_CHECK_INCOMPLETE'
  plot(rows,env,ko,ka)
  summary=['# Q550 v18.53c — '+C['title'],'', '**'+report['scientific_decision']+'**','',
   '## Dual-kappa radius at the *critical relative equilibrium*','',
   f"Original kappa = {ko:.16g}, r2_fold = {endpoint['original']['r2_fold']:.15g}, s_fold = {endpoint['original']['signed_drive_fold']:.15g}",
   f"Autonomous kappa = {ka:.16g}, r2_fold = {endpoint['autonomous']['r2_fold']:.15g}, s_fold = {endpoint['autonomous']['signed_drive_fold']:.15g}",
   '', '**Crucial:** the two fold radii are evaluated at distinct signed drives. They do NOT form a pair at a single simultaneous drive setting.',
   f"At common drive s=-1: original r2 = {same_drive['original']['r2']:.15g}, autonomous r2 = {same_drive['autonomous']['r2']:.15g}.",
   '',f"Sampled min of min(r1,r2): {env['sampled_min_radius']:.15g} at k={env['sampled_min_location_kappa']:.15g}",
   f"Sampled max of min(r1,r2): {env['sampled_max_of_min_radius']:.15g} at k={env['sampled_max_location_kappa']:.15g}",
   f"Monotonic on sampled grid: {env['monotonic_sampled']}, derivative range: {env['derivative_range']}",
   '', '## Scientific limitations','',
   '- This is a narrow envelope over prescribed kappa rules at a **fixed charge** and **counterfactual constant signed drive fold**. It is not a global lower bound on dynamical radii.',
   '- Endpoints bound *sampled* extrema only if the scan is monotonic on the sampled grid. A rigorous continuous enclosure requires interval arithmetic or a validated implicit derivative sign.',
   '- The two original/autonomous kappa values are two parameter prescriptions, not two independent seed branches at the shared anchor.',
   '- A local saddle-node normal-form diagnostic is not a global bifurcation theorem, and the historical forward-return ramp is not simulated.']
  output('SCIENTIFIC_REPORT_SERVER.md','\n'.join(summary)+'\n')
  output('report.json',report)
  print('Q550_V1853C',report['scientific_decision'],'dual_points',len(dual),'extended_points',len(extended),'native',eng is not None,flush=True)
  if args.require_parent and report['scientific_decision']!='AUTHENTIC_DUAL_KAPPA_FOLD_LOCUS_AND_SAMPLED_CRITICAL_RADIUS_ENVELOPE__NOT_GLOBAL_MINIMUM_OR_INTERVAL_PROOF':
   raise RuntimeError('STRICT_AUTHENTIC_DUAL_KAPPA_FOLD_CHECK_INCOMPLETE')
 finally:
  if temp is not None:temp.cleanup()

if __name__=='__main__':
 try:main()
 except Exception as exc:
  OUT.mkdir(exist_ok=True,parents=True);output('audit_failure.txt',traceback.format_exc())
  print('V18_53C_AUDIT_FAILURE',repr(exc),flush=True);sys.exit(1)
