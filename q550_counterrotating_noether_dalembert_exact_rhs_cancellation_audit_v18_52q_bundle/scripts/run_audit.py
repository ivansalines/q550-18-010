#!/usr/bin/env python3
"""Q550 v18.52q: exact counter-rotation charge; archived and current engine tests.

Only evaluates existing RHS and precomputed trajectories. No nonlinear IVP.
Numerical structural checks supplement (but do not replace) the analytic proof.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,math,sys,tempfile,zipfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
R51='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
S51NAME=S51+'_server_sendback.zip'
V52L='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
DEFAULT_STATES=BASE/V52L/'results'/'extended_horizon_10d_states.npz'
PARAMS=SimpleNamespace(phi0=.28,e_reserve=.26,gamma_rho=.3,barrier_B=.48,mu_x=.5,gamma_x=.02,m_x=1.2)

def sha(path):
 with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def module_from(path,tag):
 spec=importlib.util.spec_from_file_location(tag,str(path));mod=importlib.util.module_from_spec(spec)
 sys.modules[tag]=mod;spec.loader.exec_module(mod);return mod

def reference_source():
 return module_from(ROOT/'reference'/'test_support_q550_oriented_coupled_bulge_cone_audit_v1.py','q550_historical_rhs')

class IndependentPotential:
 @staticmethod
 def sigmoid(z):
  return float(1/(1+math.exp(-z))) if z>=0 else float(math.exp(z)/(1+math.exp(z)))
 @staticmethod
 def dU(r,e):return 2*r - 2.4*r*r + .8*r**3 - e*e/r**3
 @staticmethod
 def dV_barrier(x,B):return -B*(1/x-1/(1-x))


def charge(y):
 y=np.asarray(y,float);return np.exp(2*y[...,0])*y[...,3]-np.exp(2*y[...,4])*y[...,7]

def transfers(y,f):
 """Exact derivative of native C along arbitrary vector field F, any clock."""
 y=np.asarray(y,float);f=np.asarray(f,float)
 t1=np.exp(2*y[...,0])*(2*y[...,3]*f[...,0]+f[...,3])
 t2=np.exp(2*y[...,4])*(2*y[...,7]*f[...,4]+f[...,7])
 return t1,t2,t1-t2

def random_state(rng):
 y=np.zeros(11);y[0]=rng.uniform(-.15,1.25);y[4]=rng.uniform(-.2,1.25)
 y[[1,5]]=rng.uniform(-.7,.7,2);y[[3,7]]=rng.uniform(-1.1,1.1,2)
 y[[2,6,10]]=rng.uniform(-5,5,3);y[8]=rng.uniform(-1.2,1.2);y[9]=rng.uniform(-.8,.8)
 return y

def antirotate(z,angle):
 a=np.array(z,float,copy=True);a[2]+=angle;a[6]-=angle;return a

def archived_physical(historical,z,k,params,audit):
 support=np.asarray(historical.support_rhs_float(z,kappa0=k,params=params,audit=audit),float)
 x=audit.sigmoid(float(z[8]))
 return support*(x/(2*math.pi))

def archived_probes(n=240):
 rng=np.random.default_rng(1852);historical=reference_source();a=IndependentPotential();rows=[]
 for k in (0.,1/6, .1665996104,.37):
  for _ in range(n//4):
   z=random_state(rng);alpha=float(rng.uniform(-2.4,2.4))
   f=archived_physical(historical,z,k,PARAMS,a)
   f2=archived_physical(historical,antirotate(z,alpha),k,PARAMS,a)
   support=historical.support_rhs_float(z,kappa0=k,params=PARAMS,audit=a)
   t1,t2,d=transfers(z,f);t1s,t2s,ds=transfers(z,support)
   expected=k*a.sigmoid(z[8])*np.exp(z[0]+z[4])*np.sin(z[2]+z[6]-PARAMS.phi0-z[10])
   rows.append(dict(k=float(k),angle=alpha,torque1=t1,torque2=t2,
    cancellation_error=abs(d),source_torque_error=max(abs(t1-expected),abs(t2-expected)),
    support_cancellation_error=abs(ds),equivariance_error=float(np.max(abs(f2-f)))))
 err=max(x['cancellation_error'] for x in rows)
 equiv=max(x['equivariance_error'] for x in rows)
 serr=max(x['support_cancellation_error'] for x in rows)
 torque=max(x['source_torque_error'] for x in rows)
 # Counterfactual: unbalanced forcing in first angular acceleration
 z=random_state(rng);f=archived_physical(historical,z,1/6,PARAMS,a)
 broken=np.array(f);broken[3]+=.003
 broken_derivative=float(transfers(z,broken)[2]);wanted=float(math.exp(2*z[0])*.003)
 proof=abs(broken_derivative-wanted)<1e-12 and abs(broken_derivative)>1e-5
 return dict(n_probes=len(rows),max_dC_physical=err,max_dC_support=serr,
  max_anti_rotation_equivariance=equiv,max_common_torque_error=torque,
  counterfactual_break_dC=broken_derivative,counterfactual_predicted=wanted,
  counterfactual_pass=bool(proof),passes=bool(err<2e-12 and equiv<2e-12 and serr<2e-11 and torque<2e-12 and proof))

def native_traces(path,expected_hash,output):
 if not path.is_file():raise FileNotFoundError('v18.52l 10D original states missing: '+str(path))
 actual=sha(path)
 if actual!=expected_hash:raise RuntimeError('v18.52l 10D SHA mismatch: '+actual)
 out=[]
 with np.load(path,allow_pickle=False) as arr:
  for k in arr.files:
   if k=='common_xi':continue
   y=np.asarray(arr[k],float)
   if y.shape[1]<10:raise ValueError('Bad 10D array '+k)
   c=charge(y);j1=np.exp(2*y[:,0])*y[:,3];j2=np.exp(2*y[:,4])*y[:,7]
   vals=dict(series=k,n=int(len(y)),C0=float(c[0]),Cfinal=float(c[-1]),C_range=float(np.ptp(c)),
    J1_range=float(np.ptp(j1)),J2_range=float(np.ptp(j2)))
   out.append(vals)
 write_csv(output,out)
 return dict(n_series=len(out),trace_n_samples=min(t['n'] for t in out),
  max_range=max(t['C_range'] for t in out),max_individual_J_variation=max(max(t['J1_range'],t['J2_range']) for t in out),
  source_sha256=actual,passed=bool(all(t['C_range']<2e-10 for t in out))),out

def make_reference_plot(path,output):
 """No recomputation of dynamics; show native reference charge vs individual transfers."""
 try:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
 except ImportError:
  return False
 with np.load(path,allow_pickle=False) as a:
  y=np.asarray(a['original__reference__Y_tau'],float)
 t=y[:,10];j1=np.exp(2*y[:,0])*y[:,3];j2=np.exp(2*y[:,4])*y[:,7];c=j1-j2
 fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True,constrained_layout=True)
 axes[0].plot(t,j1-j1[0],label='J1(t) - J1(0)',lw=1.3)
 axes[0].plot(t,j2-j2[0],ls='--',label='J2(t) - J2(0)',lw=1.15)
 axes[0].set_ylabel('Individual angular-momentum transfers');axes[0].legend();axes[0].grid(alpha=.18)
 axes[1].plot(t,c-c[0],color='black',lw=.7,label='C(t) - C(0) = J1-J2 change')
 axes[1].set_ylabel('Charge drift (roundoff)');axes[1].set_xlabel('Original matched material time tau')
 axes[1].legend();axes[1].grid(alpha=.18)
 fig.suptitle('Q550 v18.52q — equal dynamic angular torques, exact counter-rotating charge')
 fig.savefig(output,dpi=155);plt.close(fig)
 return True

def current_engine_probes(source_path,expected_hash,n_per_seed=48):
 """Resolve historical nonlinear generator. No integration, random off-orbit tests."""
 if not source_path.is_file():raise FileNotFoundError('Original v18.51s sendback missing: '+str(source_path))
 real=sha(source_path)
 if real!=expected_hash:raise RuntimeError('v18.51s SHA changed: '+real)
 with zipfile.ZipFile(source_path) as z:
  if z.testzip():raise RuntimeError('v18.51s ZIP CRC failure')
  nested=z.read(S51+'/inputs/'+R51+'_server_sendback.zip')
  srun=z.read(S51+'/scripts/run_audit.py')
 parent=module_from(ROOT/'reference'/'parent_v18_51t_run_audit.py','q550_v1851t_original_loader')
 orig=reference_source()
 rng=np.random.default_rng(18521108)
 results=[]
 with tempfile.TemporaryDirectory(prefix='q550_v1852q_original_generator_') as tmp:
  parts=parent.stage_engine(Path(tmp),nested,srun)[4]
  _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
  for seed in ('original','autonomous'):
   ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
   for i in range(n_per_seed):
    z=random_state(rng)
    # Use the genuine original drive phase at a random event-clock point.
    xi=float(rng.uniform(0.,.8));z[10]=float(ctx['theta'](ctx['event']+xi))
    k=1/6 if seed=='original' else .1665996104
    f=np.asarray(engine.jet.structural_jet(z,k,engine.params,engine.slope)[0][:10],float)
    phi=float(rng.uniform(-2.0,2.0));rot=antirotate(z,phi)
    fr=np.asarray(engine.jet.structural_jet(rot,k,engine.params,engine.slope)[0][:10],float)
    t1,t2,dc=transfers(z,f)
    # A direct archived RHS comparison is attempted if its potential functions
    # use the same callable signatures as the original generated engine.
    comp=None;comp_full=None;comp_note=''
    try:
     ff=orig.support_rhs_float(z,kappa0=k,params=engine.params,audit=IndependentPotential())[:10]
     angular_indices=np.array([0,2,3,4,6,7],int)
     af=f[angular_indices]; bf=np.asarray(ff)[angular_indices]
     comp=float(np.max(abs(af-bf))/max(1.,float(np.max(abs(af))),float(np.max(abs(bf)))))
     comp_full=float(np.max(abs(f-ff))/max(1.,float(np.max(abs(f))),float(np.max(abs(ff)))))
     comp_note='historical_support_rhs_direct_angular_and_full'
    except Exception as ex:
     comp_note='historical_support_rhs_unavailable: '+str(ex)[:180]
    results.append(dict(seed=seed,k=k,xi=xi,torque1=float(t1),torque2=float(t2),
       charge_derivative=float(dc),equivariance_error=float(np.max(abs(fr-f))),
       relative_archived_angular_rhs_disagreement=comp,relative_archived_full_rhs_disagreement=comp_full,source_comparison_note=comp_note))
 maxdc=max(abs(x['charge_derivative']) for x in results)
 maxeq=max(x['equivariance_error'] for x in results)
 maxrelative=max((x['relative_archived_angular_rhs_disagreement'] for x in results if x['relative_archived_angular_rhs_disagreement'] is not None),default=None)
 compare_n=sum(x['relative_archived_angular_rhs_disagreement'] is not None for x in results)
 # scaled relative tolerance accounts for large time-parameterized fields
 scale=max(1.,max(abs(x['torque1'])+abs(x['torque2']) for x in results))
 return dict(input_v1851s_sha256=real,n_probes=len(results),max_absolute_dC=maxdc,
  max_absolute_anti_rotation_equivariance_error=maxeq,
  max_relative_archived_angular_rhs_disagreement=maxrelative,n_original_rhs_comparisons=compare_n,
  relative_dC_scale=maxdc/scale,passes=bool(maxdc/scale<1e-11 and maxeq<2e-10 and compare_n==len(results) and maxrelative is not None and maxrelative<2e-10)),results

def write_csv(path,rows):
 path.parent.mkdir(parents=True,exist_ok=True)
 with open(path,'w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['empty']);w.writeheader();w.writerows(rows)

def save_report(root,summary):
 (root/'report.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
 a=summary['archived_rhs'];n=summary.get('native_traces',{});e=summary.get('actual_original_engine',{})
 lines=['# Q550 v18.52q — Counter-Rotating Exact RHS Conservation Audit','',
 '**Analytical identity for archived angular RHS:** d/dt(r1² ω1 - r2² ω2) ≡ 0. This is a symbolic result from the explicit angular equations, NOT inferred from small numerical drift.','',
 '## Archived nonlinear RHS probes','',
 f"- Random off-orbit states: {a['n_probes']}; charge-derivative max |residual|: {a['max_dC_physical']:.3e} (physical clock), {a['max_dC_support']:.3e} (support clock).",
 f"- Anti-rotation equivariance max |ΔF|: {a['max_anti_rotation_equivariance']:.3e}.",
 f"- Deliberately unbalanced angular forcing produces nonzero dC/dt: {a['counterfactual_break_dC']:.9g} (predicted {a['counterfactual_predicted']:.9g}).",
 '', '## Authenticated native Q550 traces',
 f"- Series: {n.get('n_series','not run')}; largest C range: {n.get('max_range','not run')}.",
 '', '## Authenticated ORIGINAL structural_jet probe (separate from 2026-Aug archived RHS)',
 f"- Status: {e.get('status','RUN')}; tested off-trajectory states: {e.get('n_probes','not run')}; max charge derivative: {e.get('max_absolute_dC','not run')}.",
 f"- Direct original-vs-archived RHS comparisons: {e.get('n_original_rhs_comparisons','not run')}; max relative mismatch: {e.get('max_relative_archived_angular_rhs_disagreement','not run')}.",
 '', '## Interpretation','',
 '- Counter-rotation leaves Delta=theta1+theta2-phi0-Theta exactly unchanged. Both angular momenta receive the same torque. Their difference C is an exact first integral in the documented RHS.',
 '- The generalized Noether–d’Alembert interpretation is justified for the angular sector, including radial/auxiliary forces with zero angular work along X=∂theta1−∂theta2; a purely conservative full 11D Hamiltonian is NOT asserted.',
 '- The charge is real-valued, not a topological winding integer.',
 '- A positive finite-horizon collective radius is not a theorem of global no-collapse: this requires an all-future kinetic-energy upper bound and C≠0.',
 '- The invariant is in native absolute polar coordinates, not in the observer 17.04° tilted plane nor in the difference of two trajectory separations.',
 '', '## Decision',f"**{summary['status']}**",'']
 (root/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines))

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--mode',choices=('archived','server'),default='server')
 ap.add_argument('--source-states',type=Path,default=DEFAULT_STATES)
 ap.add_argument('--source-s-path',type=Path,default=BASE/S51NAME)
 ap.add_argument('--engine-probes',type=int,default=48)
 ap.add_argument('--archive-probes',type=int,default=240)
 args=ap.parse_args()
 config=json.loads((ROOT/'inputs_config.json').read_text());out=ROOT/'results';out.mkdir(exist_ok=True)
 a=archived_probes(args.archive_probes)
 n,rows=native_traces(args.source_states,config['original_10d_sha256'],out/'all_native_absolute_charge_series.csv')
 make_reference_plot(args.source_states,out/'counterrotating_angular_transfer_and_charge.png')
 report=dict(version='v18.52q',mode=args.mode,analytical_identity='EXACT_FOR_ARCHIVED_EXPLICIT_Q550_ANGULAR_RHS',
  archived_rhs=a,native_traces=n,actual_original_engine={'status':'NOT_RUN'},
  generated_from_source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'reference').glob('*.py'))})
 if args.mode=='server':
  try:
   engine,eng_rows=current_engine_probes(args.source_s_path,config['v1851s_sha256'],args.engine_probes)
   write_csv(out/'current_engine_off_orbit_probes.csv',eng_rows)
   report['actual_original_engine']=dict(status='RUN',**engine)
  except Exception as exc:
   report['actual_original_engine']=dict(status='FAILED_OR_UNAVAILABLE',error_type=type(exc).__name__,error=str(exc))
 report['checks']={'archived_rhs_random_identity':bool(a['passes']),
  'native_all_series_constant':bool(n['passed']),
  'actual_original_engine_off_orbit_verified':bool(report['actual_original_engine'].get('passes',False)),
  'deliberate_unbalanced_torque_is_detected':bool(a['counterfactual_pass'])}
 report['status']=('ORIGINAL_ENGINE_NUMERICALLY_VERIFIED_WITH_EXACT_ARCHIVED_RHS_DERIVATION' if
  all(report['checks'].values()) else 'ARCHIVED_EXACT_DERIVATION_ONLY' if args.mode=='archived' and a['passes'] and n['passed']
  else 'UNRESOLVED_OR_FAILED')
 save_report(out,report)
 print(json.dumps({k:report[k] for k in ('version','mode','status','checks')},indent=2),flush=True)
 if args.mode=='server' and not all(report['checks'].values()):raise SystemExit('Current engine audit has unresolved checks')
 if not a['passes'] or not n['passed']:raise SystemExit('Archived/reference checks failed')

if __name__=='__main__':main()
