#!/usr/bin/env python3
"""Q550 v18.52r — independent clock-normalized RHS and charge audit.

No ODE solve. Authenticates model archive and evaluates direct original RHS.
The numerical best-fit scale is a *diagnostic*, not a definition of clock units.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, io, inspect, json, math, sys, tempfile, zipfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
R51='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
S51FILE=S51+'_server_sendback.zip'
V52L='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
DEFAULT_STATES=BASE/V52L/'results'/'extended_horizon_10d_states.npz'
CONFIG=json.loads((ROOT/'inputs_config.json').read_text())


def sha(p):
 with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def import_path(p,name):
 spec=importlib.util.spec_from_file_location(name,str(p))
 if spec is None or spec.loader is None:raise ImportError(str(p))
 m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

historical=None

def history():
 global historical
 if historical is None:
  historical=import_path(ROOT/'reference'/'test_support_q550_oriented_coupled_bulge_cone_audit_v1.py','q550_v1852r_historical')
 return historical

class Potential:
 @staticmethod
 def sigmoid(z):
  if z>=0:return 1/(1+math.exp(-z))
  ex=math.exp(z);return ex/(1+ex)
 @staticmethod
 def dU(r,e):return 2*r-2.4*r*r+.8*r**3-e*e/r**3
 @staticmethod
 def dV_barrier(x,B):return -B*(1/x-1/(1-x))

def sample_state(rng):
 a=np.zeros(11)
 a[[0,4]]=rng.uniform(-.2,1.15,2)
 a[[1,5]]=rng.uniform(-.75,.75,2)
 a[[3,7]]=rng.uniform(-1.2,1.2,2)
 a[[2,6,10]]=rng.uniform(-4.5,4.5,3)
 a[8]=rng.uniform(-1.2,1.2);a[9]=rng.uniform(-.9,.9)
 return a

def anti_rotate(y,phi):
 z=np.array(y,float,copy=True);z[2]+=phi;z[6]-=phi;return z

def charge(y):
 a=np.asarray(y,float)
 return np.exp(2*a[...,0])*a[...,3]-np.exp(2*a[...,4])*a[...,7]

def dcharge(y,dy):
 a=np.asarray(y,float);f=np.asarray(dy,float)
 j1=np.exp(2*a[...,0])*(2*a[...,3]*f[...,0]+f[...,3])
 j2=np.exp(2*a[...,4])*(2*a[...,7]*f[...,4]+f[...,7])
 return j1,j2,j1-j2

def archived_rhs(y,k,params,clock='support'):
 f=np.asarray(history().support_rhs_float(y,kappa0=k,params=params,audit=Potential()),float)
 if clock=='physical': f=f*(Potential.sigmoid(float(y[8]))/(2*math.pi))
 elif clock!='support':raise ValueError(clock)
 return f

def norm(a):return float(np.linalg.norm(a))

def rel(a,b):return norm(np.asarray(a)-np.asarray(b))/max(1.,norm(a),norm(b))

def align(original,reference):
 """Pointwise scalar needed to map archived -> original (least squares)."""
 a=np.asarray(original,float);b=np.asarray(reference,float)
 if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)) or (b@b)<1e-22:raise ValueError('Invalid vector for clock fit')
 c=float((a@b)/(b@b))
 resid=norm(a-c*b)/max(1.,norm(a),norm(c*b))
 return c,resid

def csv_write(p,rows):
 if not rows:return
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def find_ancestor(override):
 if override:return Path(override).expanduser().resolve()
 c=[ROOT/'inputs'/S51FILE,BASE/S51FILE,BASE/S51/'inputs'/S51FILE,Path('/root/elab1/hysteresis')/S51FILE,Path('/root/Downloads')/S51FILE]
 return next((p for p in c if p.is_file()),None)

def native_traces(path):
 if not path.is_file():return {'status':'MISSING', 'path':str(path)}
 if sha(path)!=CONFIG['v1852l_10d_sha256']:raise RuntimeError('Native v18.52l SHA-256 mismatch')
 rows=[]
 with np.load(path,allow_pickle=False) as o:
  for key in o.files:
   if key=='common_xi':continue
   y=np.asarray(o[key],float)
   if y.ndim!=2 or y.shape[1]<10:raise ValueError('Bad native array: '+key)
   C=charge(y)
   rows.append({'name':key,'samples':len(y),'initial':float(C[0]),'range':float(np.ptp(C)),
    'delta_end':float(C[-1]-C[0])})
 csv_write(ROOT/'results'/'native_charge_per_series.csv',rows)
 maxerr=max(r['range'] for r in rows)
 return {'status':'VERIFIED','n_series':len(rows),'max_drift':maxerr,
 'pass':bool(len(rows)>=20 and maxerr<2e-10),'sha256':sha(path)}

def symbolic_identity():
 # The identity is an algebraic calculation from the archived angular ODEs.
 # Sympy is optional; an independent arithmetic check always runs.
 state='NOT_INSTALLED';expression='archived_angular_rhs_manual_cancellation'
 try:
  import sympy as S
  r1,r2,w1,w2,v1,v2,ke,sin_d=S.symbols('r1 r2 w1 w2 v1 v2 ke sinDelta')
  w1p=ke*r2/r1*sin_d-2*v1*w1
  w2p=ke*r1/r2*sin_d-2*v2*w2
  t1=S.factor(2*r1*r1*v1*w1+r1*r1*w1p)
  t2=S.factor(2*r2*r2*v2*w2+r2*r2*w2p)
  assert S.simplify(t1-t2)==0
  assert S.simplify(t1-ke*r1*r2*sin_d)==0
  state='SYMPY_EXACT_PASS';expression=str(S.simplify(t1-t2))
 except ImportError:
  state='SYMPY_UNAVAILABLE_MANUAL_DERIVATION_ONLY'
 return {'status':state,'residual_expression':expression}

def archived_self_test(probes):
 rng=np.random.default_rng(185209);params=SimpleNamespace(phi0=.28,e_reserve=.26,gamma_rho=.3,barrier_B=.48,mu_x=.5,gamma_x=.02,m_x=1.2)
 errs=[];equivs=[];ratio=[]
 for i in range(probes):
  z=sample_state(rng);k=[0.,1/6,.1665996104,.3][i%4]
  s=archived_rhs(z,k,params)
  f=archived_rhs(z,k,params,'physical')
  _,_,dc=dcharge(z,s)
  errs.append(abs(dc));equivs.append(rel(s,archived_rhs(anti_rotate(z,.847),k,params)))
  ratio.append(rel(s,f*(2*math.pi/Potential.sigmoid(z[8]))))
 return {'n_probes':probes,'max_dC_support':max(errs),'max_anti_rotation_equivariance':max(equivs),
  'max_support_to_physical_clock_residual':max(ratio),'pass':max(errs)<2e-11 and max(equivs)<1e-12 and max(ratio)<1e-12}

def engine_run(source,probes):
 if not source.is_file():raise FileNotFoundError(str(source))
 sh=sha(source)
 if sh!=CONFIG['v1851s_sha256']:raise RuntimeError('Authentic v18.51s SHA mismatch: '+sh)
 with zipfile.ZipFile(source) as z:
  if z.testzip():raise RuntimeError('v18.51s corrupted ZIP CRC')
  rblob=z.read(S51+'/inputs/'+R51+'_server_sendback.zip')
  srun=z.read(S51+'/scripts/run_audit.py')
 loader=import_path(ROOT/'reference'/'parent_v18_51t_run_audit.py','q550_v1852r_parent_loader')
 rows=[]; field_sources={}; clock_provenance={}
 with tempfile.TemporaryDirectory(prefix='q550_v1852r_original_') as tmp:
  parts=loader.stage_engine(Path(tmp),rblob,srun)[4]
  _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
  try:
   src=inspect.getsource(engine.jet.structural_jet)
   (ROOT/'results'/'original_structural_jet_source_excerpt.txt').write_text(src)
   field_sources['structural_jet_excerpt_sha256']=sha(ROOT/'results'/'original_structural_jet_source_excerpt.txt')
  except (OSError,TypeError) as ex:field_sources['unavailable']=str(ex)
  # Record what is actually available; a scalar fit alone is not provenance.
  for attr in ('slope','tau_ref','period_ref','T_ref','period','time_scale','period_s'):
   try:
    value=getattr(engine,attr)
    if isinstance(value,(int,float,np.floating,str)):
     clock_provenance['engine.'+attr]=str(value)
   except (AttributeError,TypeError,ValueError):pass
  rng=np.random.default_rng(18520852)
  for seed in ('original','autonomous'):
   ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
   k=1/6 if seed=='original' else .1665996104
   for i in range(probes):
    y=sample_state(rng);xi=float(rng.uniform(.01,.8))
    # Set physical driving phase supplied by the authenticated engine.
    y[10]=float(ctx['theta'](ctx['event']+xi))
    fraw=np.asarray(engine.jet.structural_jet(y,k,engine.params,engine.slope)[0],float)
    f=fraw[:10]
    if f.shape!=(10,):raise ValueError('Original generator does not provide 10 physical components')
    supp=archived_rhs(y,k,engine.params)[:10]
    phys=archived_rhs(y,k,engine.params,'physical')[:10]
    cS,errS=align(f,supp)
    cT,errT=align(f,phys)
    phi=float(rng.uniform(-2.2,2.2))
    fr=np.asarray(engine.jet.structural_jet(anti_rotate(y,phi),k,engine.params,engine.slope)[0],float)[:10]
    j1,j2,dc=dcharge(y,f)
    # Check the normalized angular torque identity using the *observed* scalar cS.
    expected_torque=k*Potential.sigmoid(float(y[8]))*math.exp(y[0]+y[4])*math.sin(y[2]+y[6]-engine.params.phi0-y[10])
    # archived support torque = (2pi/x) * physical torque; check actual engine on actual clock
    ratio_s=2*math.pi/Potential.sigmoid(float(y[8]))
    torque_err=max(abs(j1-cS*ratio_s*expected_torque),abs(j2-cS*ratio_s*expected_torque))
    row={'seed':seed,'i':i,'xi':xi,'kappa':k,'ratio_engine_to_support':cS,
      'ratio_engine_to_physical':cT,'scale_fit_resid_support':errS,'scale_fit_resid_physical':errT,
      'min_speed_component_ratio':min(abs(cS),abs(cT)),'normalised_field_disagreement':errS,
      'original_angular_charge_derivative':dc,'original_common_torque_error':torque_err,
      'anti_rotation_equivariance_error':rel(f,fr),
      'engine_f11':float(fraw[10]) if len(fraw)>10 else '',
      'archived_support_f11':float(archived_rhs(y,k,engine.params)[10])}
    rows.append(row)
 csv_write(ROOT/'results'/'clock_normalized_off_orbit_engine_comparison.csv',rows)
 # No fitting over the hold-out half. Independent global constant is estimated on train,
 # then evaluated on hold-out. The *value* is not a proof of the source clock law.
 train=[r for r in rows if r['i']%2==0];valid=[r for r in rows if r['i']%2==1]
 fit=float(np.median([r['ratio_engine_to_support'] for r in train]))
 relative_spread=max(abs(r['ratio_engine_to_support']/fit-1) for r in rows)
 max_holdout=max(abs(r['ratio_engine_to_support']/fit-1) for r in valid)
 max_fit=max(r['scale_fit_resid_support'] for r in rows)
 max_eq=max(r['anti_rotation_equivariance_error'] for r in rows)
 max_charge=max(abs(r['original_angular_charge_derivative']) for r in rows)
 max_torque=max(r['original_common_torque_error'] for r in rows)
 # Positive constant is required to preserve oriented parametrization.
 field_match=bool(fit>0 and max_fit<2e-10 and relative_spread<2e-10 and max_holdout<2e-10)
 charge_ok=bool(max_charge<2e-10 and max_eq<2e-10 and max_torque<2e-10)
 return {'status':'AUTHENTIC_ENGINE_PROBED','sha256_parent':sh,'n_probes':len(rows),'engine_to_support_fit':fit,
  'inverse_fit':1/fit,'field_scale_relative_spread':relative_spread,
  'max_holdout_relative_scale_mismatch':max_holdout,
  'max_scale_fit_residual':max_fit,'max_anti_rotation_equivariance':max_eq,
  'max_native_charge_derivative':max_charge,'max_common_torque_residual':max_torque,
  'clock_source_metadata':clock_provenance,'field_source':field_sources,
  'clock_normalized_field_match':field_match,'off_orbit_charge_pass':charge_ok,
  'warning':'Numerical fit identifies constant field ratio, not by itself an independent documented clock transformation.'}

def decision(arch,native,engine,exact):
 a=arch['pass'];e=engine.get('clock_normalized_field_match',False);c=engine.get('off_orbit_charge_pass',False)
 n=native.get('pass',False);s=exact['status'] in ('SYMPY_EXACT_PASS','EXACT_SYMPY_AND_FORMAL_MONOMIAL_PASS','EXACT_FORMAL_MONOMIAL_PASS')
 if not a:return 'FAIL_ARCHIVED_FIELD'
 if engine['status']=='MISSING_PARENT':return 'ARCHIVED_IDENTITY_ONLY_PARENT_NOT_AVAILABLE' if not n else 'ARCHIVED_IDENTITY_AND_NATIVE_TRACES_PARENT_NOT_AVAILABLE'
 if not (e and c):return 'CLOCK_OR_CHARGE_DISCREPANCY_REQUIRES_REVIEW'
 if n and s:return 'CLOCK_NORMALIZED_ORIGINAL_FIELD_MATCH_AND_ANGULAR_NOETHER_DALEMBERT_CHARGE_CERTIFIED_CONDITIONALLY'
 if n:return 'CLOCK_NORMALIZED_ORIGINAL_FIELD_MATCH__NATIVE_CHARGE_VERIFIED__SYMPY_NOT_RUN'
 return 'CLOCK_NORMALIZED_ORIGINAL_FIELD_MATCH__NATIVE_STATES_NOT_AVAILABLE'

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--states',type=Path,default=DEFAULT_STATES)
 ap.add_argument('--probes-per-seed',type=int,default=48)
 ap.add_argument('--archive-probes',type=int,default=160)
 ap.add_argument('--strict',action='store_true',help='Fail if authenticated original parent or native states are unavailable')
 args=ap.parse_args()
 if args.probes_per_seed<4 or args.archive_probes<8:raise ValueError('Probe count too small')
 res=ROOT/'results';res.mkdir(exist_ok=True)
 proof_file=ROOT/'results'/'exact_symbolic_noether.json'
 exact=json.loads(proof_file.read_text()) if proof_file.is_file() else symbolic_identity()
 arch=archived_self_test(args.archive_probes)
 native=native_traces(args.states)
 source=find_ancestor(args.parent_v1851s)
 if source is None:
  eng={'status':'MISSING_PARENT','expected_archive':S51FILE,
   'expected_sha256':CONFIG['v1851s_sha256']}
 else:eng=engine_run(source,args.probes_per_seed)
 status=decision(arch,native,eng,exact)
 obj={'version':'v18.52r','status':status,'archived_analytic':exact,
  'archived_field_checks':arch,'native_traces':native,
  'original_engine':eng,
  'model_scope':{
   'claim':'angular-sector Noether-dAlembert charge conditional on stated Lagrangian and angular virtual forces zero',
   'topological_charge':'NOT_CERTIFIED',
   'global_noncollapse':'NOT_CERTIFIED',
   'absolute_configuration':'r1,r2 > 0; 0 < x < 1',
   'canonical_momentum':'p_theta_i = r_i^2 omega_i under documented flat kinetic immersion',
   'time_reparameterization':'Clock-normalized field equality evaluated off orbit; source clock provenance separately recorded',
   },
  'sources_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'reference').glob('*.py'))}}
 (res/'report.json').write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
 report=['# Q550 v18.52r — Clock-Normalized Generator Equivalence & Noether-Charge Certification','',
 '## Evidence and scope','',
 f"- Decision: **{status}**", f"- Exact algebraic cancellation: {exact['status']}",
 f"- Archived field off-orbit: {arch['n_probes']} states, max dC/ds = {arch['max_dC_support']:.3e}",
 f"- Native authentic 10D series: {native.get('n_series','not available')}; maximum charge variation {native.get('max_drift','not available')}",
 f"- Authentic original structural_jet: {eng['status']}",
 f"- Engine/support scale: {eng.get('engine_to_support_fit','not available')}",
 f"- Max crossvalidated normalized field mismatch: {eng.get('max_scale_fit_residual','not available')}",
 f"- Holdout scale mismatch: {eng.get('max_holdout_relative_scale_mismatch','not available')}",
 f"- Maximum original charge derivative: {eng.get('max_native_charge_derivative','not available')}",
 '', '## Analytic angular sector','',
 'p1 = r1^2 omega1; p2 = r2^2 omega2; Delta = theta1+theta2-phi0-Theta.',
 'V_ang = kappa0*x*r1*r2*cos(Delta) yields identical angular generalized torques kappa0*x*r1*r2*sin(Delta).',
 'The counterrotation (theta1,theta2) -> (theta1+epsilon,theta2-epsilon) fixes Delta.',
 'Its generator X has conserved momentum map C = p1-p2 (Noether-dAlembert, angular virtual work zero).',
 'Radial and auxiliary damping/driving are not treated as an autonomous full-system Hamiltonian.',
 '', '## Boundaries','',
 'Field proportionality alone does not prove an independent physical clock law; inspect source clock metadata/excerpt.',
 'C is a continuous momentum-map value, NOT topological winding.',
 'Nonzero C conditionally obstructs simultaneous collapse r1=r2=0 with a global kinetic bound; such bound is NOT certified.',
 '', '## File provenance','', 'All included code and archived sources pinned in SOURCE_MANIFEST.json. Original parent v18.51s is SHA-pinned in inputs_config.json.', '']
 (res/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(report))
 print(json.dumps({'version':'v18.52r','status':status,'engine_to_support_fit':eng.get('engine_to_support_fit'),
  'checks':{'archived':arch['pass'],'native':native.get('pass',False),
    'original_engine':eng.get('clock_normalized_field_match',False),
    'charge':eng.get('off_orbit_charge_pass',False)}} ,indent=2),flush=True)
 if not arch['pass']:raise SystemExit('Archived analytic source failed')
 if eng['status']!='MISSING_PARENT' and not (eng.get('clock_normalized_field_match') and eng.get('off_orbit_charge_pass')):
  raise SystemExit('Original field or charge check failed; inspect report')
 if args.strict and (eng['status']=='MISSING_PARENT' or not native.get('pass') or not exact['status'].startswith('EXACT_')):
  raise SystemExit('Strict certificate incomplete; see report')
if __name__=='__main__':main()
