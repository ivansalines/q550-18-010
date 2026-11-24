#!/usr/bin/env python3
"""Q550 v18.52s-r2: exact energy ledger, native traces, finite-time vs uniform bounds.
No ODE solver. Archives remain untouched. No use of other generated v18.52s outputs.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,sys,tempfile,zipfile,io
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from energy_math import Params,sigmoid,geometry,energy_terms,exact_energy_derivative,exact_H_derivative,a_priori_bound

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent
CONFIG=json.loads((ROOT/'inputs_config.json').read_text())
P=Params(**CONFIG['model_params'])
L='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
DEFAULT=BASE/L/'results/extended_horizon_10d_states.npz'
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
R51='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
S51FILE=S51+'_server_sendback.zip'
DEFAULT_S51=BASE/S51FILE


def sha(path):
 with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def csv_out(path,rows):
 if not rows:return
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)

def historical_source():
 p=ROOT/'reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py'
 spec=importlib.util.spec_from_file_location('q550_v1852s_historical',p)
 m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
 return m

class Derivative:
 @staticmethod
 def sigmoid(z):return float(sigmoid(z))
 @staticmethod
 def dU(r,e):return float(2*r-2.4*r*r+.8*r*r*r-e*e/r**3)
 @staticmethod
 def dV_barrier(x,B):return float(-B*(1/x-1/(1-x)))

def rhs_physical(historical,y,k):
 F=np.asarray(historical.support_rhs_float(np.asarray(y,float),kappa0=k,params=SimpleNamespace(**CONFIG['model_params']),audit=Derivative()),float)
 x=float(sigmoid(y[8]));return F*(x/(2*math.pi))

def independent_off_orbit(n):
 hist=historical_source();rng=np.random.default_rng(18525201)
 vals=[]; hvals=[]; finite_checks=[]; charge=[]
 for i in range(n):
  y=np.zeros(11);y[[0,4]]=rng.uniform(-1,1.3,2);y[[1,5]]=rng.uniform(-.6,.6,2)
  y[[2,6,10]]=rng.uniform(-4,4,3);y[[3,7]]=rng.uniform(-1,1,2)
  y[8]=rng.uniform(-2,2);y[9]=rng.uniform(-.7,.7)
  k=[0.,1/6,.1665996104,.37][i%4]
  F=rhs_physical(hist,y,k);t=energy_terms(y,k,P)
  e=exact_energy_derivative(y,F,k,P);h=exact_H_derivative(y,F,k,P)
  vals.append(float(e-t['E_dot']));hvals.append(float(h-t['H_dot']))
  _,_,x,_,_,_,_,_=geometry(y,P);r1,r2,_,_,_,_,_,_=geometry(y,P)
  if abs(k)<2:
   G=float(t['E'])+86.4
   lower=.1*(r1**4+r2**4)+(1-abs(k)/2)*(r1*r1+r2*r2)+.5*P.e_reserve**2*(r1**-2+r2**-2)+P.barrier_B*math.log(4)
   finite_checks.append(float(lower-G))
  # independent derivative of counterrotation charge
  charge.append(float(r1*r1*(2*y[3]*F[0]+F[3])-r2*r2*(2*y[7]*F[4]+F[7])))
 a=max(map(abs,vals));b=max(map(abs,hvals));c=max(finite_checks);d=max(map(abs,charge))
 return dict(n=n,max_abs_E_balance=a,max_abs_H_balance=b,max_lower_bound_violation=max(0.,c),max_dJdiff=d,pass_all=bool(a<1e-10 and b<1e-10 and c<1e-12 and d<1e-11))

def integrate_on_tau(tau,flux):
 return np.r_[0.,np.cumsum(.5*(flux[1:]+flux[:-1])*np.diff(tau))]

def import_path(p,name):
 spec=importlib.util.spec_from_file_location(name,str(p))
 m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def locate_engine_source(override):
 candidates=([override] if override else [])+[DEFAULT_S51,ROOT/'inputs'/S51FILE,
  Path('/root/elab1/hysteresis')/S51FILE,Path('/root/Downloads')/S51FILE]
 return next((p for p in candidates if p and p.is_file()),None)

def recover_authenticated_global_phase(source,xi):
 """The NPZ 11th column is elapsed material tau, NEVER global Theta!
 The phase Theta(xi) is externally supplied by the authenticated original generator.
 """
 if sha(source)!=CONFIG['v1851s_sha256']:
  raise RuntimeError('Authentic v18.51s SHA-256 mismatch: '+sha(source))
 with zipfile.ZipFile(source) as z:
  if z.testzip():raise RuntimeError('Original v18.51s archive CRC error')
  rblob=z.read(S51+'/inputs/'+R51+'_server_sendback.zip')
  run_source=z.read(S51+'/scripts/run_audit.py')
 with zipfile.ZipFile(io.BytesIO(rblob)) as rzip:
  if rzip.testzip():raise RuntimeError('Nested v18.51r ZIP CRC error')
  native_ancestor=rzip.read(R51+'/results/extended_paired_kappa_nonlinear_states.npz')
 with np.load(io.BytesIO(native_ancestor),allow_pickle=False) as n:
  k_orig=float(n['kappa_original']);k_auto=float(n['kappa_autonomous'])
 k_common=(k_orig+k_auto)/2
 if abs(k_common-CONFIG['v1852l_fixed_model_kappa'])>1e-12:
  raise RuntimeError('Fixed model kappa disagrees with authenticated v18.52l report: '+str(k_common))
 parent=import_path(ROOT/'reference/parent_v18_51t_run_audit.py','q550_v1852s_parent_loader')
 with tempfile.TemporaryDirectory(prefix='q550_v1852s_original_clock_') as temp:
  parts=parent.stage_engine(Path(temp),rblob,run_source)[4]
  _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
  phase={}
  for seed in ('original','autonomous'):
   ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
   event=float(ctx['event'])
   phase[seed]=np.array([float(ctx['theta'](event+float(t))) for t in xi])
  # The original model parameters cannot be silently replaced by the fixed configuration.
  params=CONFIG['model_params']
  authenticated={key:float(getattr(engine.params,key)) for key in params}
  differences={key:authenticated[key]-float(value) for key,value in params.items()}
  parameter_audit={'status':'MATCH' if all(abs(v)<=1e-12 for v in differences.values()) else 'MISMATCH',
   'parent_v1851s_sha256':sha(source),'engine_authenticated_parameters':authenticated,
   'audit_declared_parameters':params,'original_v1852s_unverified_parameters':CONFIG.get('supersedes_wrong_parameter_config',{}),
   'differences':differences,'note':'Never substitute parameters merely to force an energy ledger match'}
  (ROOT/'results').mkdir(parents=True,exist_ok=True)
  (ROOT/'results'/'parameter_audit.json').write_text(json.dumps(parameter_audit,indent=2)+'\n')
  if parameter_audit['status']!='MATCH':
   raise RuntimeError('AUTHENTICATED_MODEL_PARAMETER_MISMATCH (see results/parameter_audit.json): '+str(differences))
  clock_meta={'parent_sha256':sha(source),'original_engine_slope':float(engine.slope),
   'kappa_original':k_orig,'kappa_autonomous':k_auto,'kappa_shared':k_common,
   'xi_interval':[float(xi[0]),float(xi[-1])],
   'phase_interval':{key:[float(v[0]),float(v[-1])] for key,v in phase.items()},
   'theta_is_external':True,'Y_tau_column_11_is_material_time':True,
   'parameter_authentication':'SHA256_SOURCE_AND_EXACT_SEVEN_PARAM_MATCH',
   'authenticated_model_parameters':authenticated}
 return phase,k_common,clock_meta

def evaluate_trace(Y,k,stride,seed,name,allrows,plot_source,phase):
 indices=np.arange(0,len(Y),stride)
 if indices[-1]!=len(Y)-1:indices=np.r_[indices,len(Y)-1]
 y=np.asarray(Y[indices],float)
 phi=np.asarray(phase[indices],float)
 t=energy_terms(y,k,P,phase=phi)
 tau=y[:,10];x=t['x'];E=t['E'];H=t['H'];G=E+86.4
 assert np.all(np.isfinite(tau)) and np.all(np.diff(tau)>0), 'Native 11th column tau must increase'
 assert np.all((x>0)&(x<1))
 # A global phase rate is reconstructed from the lineage function, not from Y[:,10].
 phase_tau=np.gradient(phi,tau,edge_order=2)
 theta_external_work=k*x*t['r1']*t['r2']*np.sin(t['D'])*phase_tau
 extra=theta_external_work-t['W_Theta']
 energy_rate=t['E_dot']+extra
 h_rate=t['H_dot']+extra
 physical_elapsed=tau-tau[0]
 work_x=integrate_on_tau(tau,t['W_x']);work_th=integrate_on_tau(tau,theta_external_work)
 loss_rad=integrate_on_tau(tau,t['D_rho']);loss_x=integrate_on_tau(tau,t['D_x'])
 ledger=work_x+work_th-loss_rad-loss_x
 hwork=integrate_on_tau(tau,t['W_x']-.5*t['dx']*t['Jsum']+extra)
 hledger=hwork-loss_rad-loss_x
 b=a_priori_bound(k,P)
 bound=G[0]*np.exp(b['K']*physical_elapsed)
 e_bound_violation=float(np.max(np.maximum(G-bound,0.)/np.maximum(bound,1e-30)))
 rmin=np.minimum(t['r1'],t['r2']);r_total=t['rhoR'];J=t['Jdiff']
 check_r=abs(P.e_reserve)/np.sqrt(2*bound)
 check_charge=np.abs(J[0])/np.sqrt(2*bound)
 phase_error=float(np.max(abs(phase_tau-x)))
 rows=dict(series=name,seed=seed,steps=int(len(y)),kappa0=k,model_parameter_source='authenticated v18.51s engine.params',
  xi_phase_source='authenticated v18.51s ctx[theta](event+xi)',
  tau_interval=float(tau[-1]-tau[0]),phase_interval=float(phi[-1]-phi[0]),
  phase_rate_minus_x_max_abs=phase_error,
  E0=float(E[0]),E_end=float(E[-1]),E_range=float(np.ptp(E)),
  H0=float(H[0]),H_end=float(H[-1]),H_range=float(np.ptp(H)),
  energy_delta=float(E[-1]-E[0]),budget_delta=float(ledger[-1]),
  energy_closure_absolute=float((E[-1]-E[0])-ledger[-1]),
  energy_closure_relative_to_budget_scale=float(abs((E[-1]-E[0])-ledger[-1])/max(1.,np.ptp(E),abs(ledger[-1]))),
  H_closure_absolute=float((H[-1]-H[0])-hledger[-1]),
  net_drive=float(work_x[-1]+work_th[-1]),total_dissipation=float(loss_rad[-1]+loss_x[-1]),
  net_Wx=float(work_x[-1]),net_WTheta=float(work_th[-1]),
  observed_min_r1=float(np.min(t['r1'])),observed_min_r2=float(np.min(t['r2'])),
  observed_max_Tang=float(np.max(t['T_ang'])),observed_max_G=float(np.max(G)),
  bound_K=b['K'],bound_final_G=float(bound[-1]),energy_bound_max_relative_violation=e_bound_violation,
  min_radial_inequality_margin=float(np.min(rmin-check_r)),
  min_combined_charge_inequality_margin=float(np.min(r_total-check_charge)),
  charge_range=float(np.ptp(J)),
  observed_energy_upper_bound_is_global=False,global_uniform_noncollapse_certified=False)
 allrows.append(rows)
 if name in ('original__reference__Y_tau','autonomous__reference__Y_tau'):
  plot_source[seed]={'theta':tau,'elapsed':physical_elapsed,'E':E,'H':H,'energy_budget':ledger,
   'H_budget':hledger,'work_x':work_x,'work_Theta':work_th,'diss':loss_rad+loss_x,
   'r1':t['r1'],'r2':t['r2'],'T_ang':t['T_ang'],'Wdot':t['W_x']+theta_external_work,
   'Ddot':t['D_rho']+t['D_x']}
 return rows

def make_plots(d):
 try:
  import matplotlib
  matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return False
 for seed,v in d.items():
  th=v['theta']; E=v['E'];H=v['H']
  fig,axes=plt.subplots(3,1,figsize=(11,9),sharex=True,constrained_layout=True)
  axes[0].plot(th,E-E[0],lw=1.2,label='Measured E - E0')
  axes[0].plot(th,v['energy_budget'],lw=1,ls='--',label='Integrated work - dissipation')
  axes[0].set_ylabel('Energy change');axes[0].legend();axes[0].grid(alpha=.2)
  axes[1].plot(th,H-H[0],lw=1.2,label='Rotating-frame H - H0')
  axes[1].plot(th,v['H_budget'],lw=1,ls='--',label='Integrated H budget')
  axes[1].set_ylabel('Rotating-frame energy');axes[1].legend();axes[1].grid(alpha=.2)
  axes[2].plot(th,v['Wdot'],label='Work injection rate',lw=1)
  axes[2].plot(th,v['Ddot'],label='Dissipated power (positive)',lw=1)
  axes[2].set_ylabel('Material-time power');axes[2].set_xlabel('Material time tau');axes[2].legend();axes[2].grid(alpha=.2)
  fig.suptitle('Q550 v18.52s-r2 — authentic '+seed+' reference; exact energy ledger')
  fig.savefig(ROOT/'results'/f'{seed}_energy_work_dissipation_ledger.png',dpi=145);plt.close(fig)
 return True

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--states',type=Path,default=DEFAULT)
 ap.add_argument('--stride',type=int,default=1);ap.add_argument('--off-orbit',type=int,default=160)
 ap.add_argument('--require-traces',action='store_true');ap.add_argument('--parent-v1851s',type=Path);args=ap.parse_args()
 if args.stride<1 or args.off_orbit<16:raise ValueError('stride>=1 off-orbit>=16')
 (ROOT/'results').mkdir(exist_ok=True)
 probes=independent_off_orbit(args.off_orbit)
 traces={'status':'MISSING','path':str(args.states)};rows=[];plot_source={}
 parent_source=locate_engine_source(args.parent_v1851s)
 phase_context={'status':'MISSING_OR_UNAVAILABLE','required':S51FILE,'expected_sha256':CONFIG['v1851s_sha256']}
 if args.states.is_file():
  actual=sha(args.states)
  expected=CONFIG['v1852l_10d_sha256']
  if actual!=expected:raise RuntimeError('v18.52l original 10D SHA mismatch: '+actual)
  if parent_source is None:
   traces=dict(status='NATIVE_STATES_AUTHENTIC_BUT_PHASE_CONTEXT_MISSING',file_sha256=actual,
    warning='Native column 11 stores material tau, NOT global Theta; cannot evaluate the energy ledger without authenticated phase source')
   if args.require_traces:raise FileNotFoundError('v18.51s original phase source missing; need '+S51FILE)
  else:
   with np.load(args.states,allow_pickle=False) as z:
    xi=np.asarray(z['common_xi'],float)
    phase_by_seed,k_shared,clock_meta=recover_authenticated_global_phase(parent_source,xi)
    phase_context=dict(status='AUTHENTIC_ORIGINAL_PHASE_RECONSTRUCTED',**clock_meta)
    for name in z.files:
     if name=='common_xi':continue
     seed='original' if name.startswith('original') else 'autonomous'
     evaluate_trace(z[name],k_shared,args.stride,seed,name,rows,plot_source,phase_by_seed[seed])
   csv_out(ROOT/'results/native_20_series_energy_ledger.csv',rows)
   make_plots(plot_source)
   maxerr=max(r['energy_closure_relative_to_budget_scale'] for r in rows)
   phase_error=max(r['phase_rate_minus_x_max_abs'] for r in rows)
   traces=dict(status='AUTHENTIC_NATIVE_AND_PHASE_VERIFIED',file_sha256=actual,
     phase_source_sha256=sha(parent_source),n_series=len(rows),min_count=min(r['steps'] for r in rows),
     max_energy_ledger_relative_error=maxerr,max_phase_rate_minus_x=phase_error,
     max_charge_range=max(r['charge_range'] for r in rows),
     max_energy_upper_bound_violation=max(r['energy_bound_max_relative_violation'] for r in rows),
     ledger_pass=bool(maxerr<5e-3),phase_clock_pass=bool(phase_error<5e-3),
     bound_pass=bool(max(r['energy_bound_max_relative_violation'] for r in rows)<1e-11),
     plots_generated=list(plot_source))
 elif args.require_traces:raise FileNotFoundError('v18.52l 10D input missing: '+str(args.states))
 theorem={k:a_priori_bound(v,P) for k,v in CONFIG['kappa_seeds'].items()}
 report={'version':'v18.52s-r2','model':'Q550 archived physical-time RHS',
  'proof_scope':'Exact parameter-dependent model energy identity; verified parameter lineage required for authentic native ledger; finite-time bound only',
  'parameter_audit_file':'results/parameter_audit.json' if (ROOT/'results/parameter_audit.json').is_file() else None,
  'model_parameters':CONFIG['model_params'],
  'off_orbit':probes,'native':traces,'phase_context':phase_context,'finite_time_bounds':theorem,
  'scientific_decision':('ENERGY_IDENTITY_AND_FINITE_TIME_NONCOLLAPSE_BOUND_WITH_AUTHENTIC_LEDGER' if probes['pass_all'] and traces['status']=='AUTHENTIC_NATIVE_AND_PHASE_VERIFIED' and traces['ledger_pass'] and traces['phase_clock_pass'] and traces['bound_pass'] else 'ANALYTIC_PROOF_WITH_PARTIAL_OR_FAILED_TRACE_SUPPORT'),
  'no_global_uniform_bound_certified':True,
  'not_proved':['global uniform energy upper bound','uniform-in-time radius lower bound','no collapse of 5D perturbation separation','topological charge interpretation']}
 (ROOT/'results/report.json').write_text(json.dumps(report,indent=2)+'\n')
 md=['# Q550 v18.52s — Energy–Work–Dissipation & Noncollapse Scientific Report',
 '','## Classification',report['scientific_decision'],'',
 '## Exact physical-time identity',
 'E= 1/2 Σ(rdot_i²+r_i²ω_i²)+1/2*m_x*xdot²+ΣU(r_i)+V_B(x)+κ0*x*r1*r2*cos(Δ).',
 'dE/dt = -γρ Σrdot_i²-γx xdot² + r1*r2*(μx*sinΔ+κ0*cosΔ)*xdot + κ0*x²*r1*r2*sinΔ.',
 'Δ=θ1+θ2−φ0−Θ; Θdot=x. The last term is a phase-drive work rate, NOT dissipation.',
 '','## Rotating frame',
 'H=E - (x/2)(J1+J2), Hdot= -γρ Σrdot_i²-γx xdot² +xdot*[r1*r2*(μx*sinΔ+κ0*cosΔ)-(J1+J2)/2].',
 '','## Finite-time theorem vs uniform global theorem',
 'For |κ0|<2, e_reserve>0, B>0, mx>0 and nonnegative damping, let G=E+86.4.',
 'An explicit lower bound gives G ≥ (1/10)Σr_i^4+(1-|κ0|/2)Σr_i²+(e²/2)Σr_i^(-2)+(mx/2)xdot²+B ln4+other nonnegative kinetic terms.',
 'Thus Gdot≤K G, K=sqrt(μx²+κ0²)*(2.5+1/mx)+|κ0|/[2(1-|κ0|/2)].',
 'Therefore on any finite regular physical-time interval G(t)≤G(0)exp(Kt), and r_i(t)≥|e|/sqrt(2G(0)exp(Kt)).',
 'For C=J1−J2 constant, r1²+r2²≥C²/[2G(0)exp(Kt)].',
 '**This excludes finite-time collapse of either absolute radius, conditional on continuation and the stated model domain; it does NOT yield a time-independent lower bound as t→∞.**',
 'No theorem concerning collapse of relative 5D separations follows from this radial statement.',
 '','## Independently evaluated RHS probes',json.dumps(probes,indent=2),'','## Native traces',json.dumps(traces,indent=2),'','## External phase provenance',json.dumps(phase_context,indent=2),
 '','## Provenance', 'Historical RHS: reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py; exact SHA 10D native v18.52l in inputs_config.json.',
 'No nonlinear ODE integration was performed by this audit.']
 (ROOT/'results/SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(md)+'\n')
 print(json.dumps({'decision':report['scientific_decision'],'off_orbit':probes,'native':traces,'K':{s:v['K'] for s,v in theorem.items()}},indent=2))
 if not probes['pass_all'] or (traces['status']=='AUTHENTIC_NATIVE_AND_PHASE_VERIFIED' and (not traces['ledger_pass'] or not traces['bound_pass'] or not traces['phase_clock_pass'])):
  raise SystemExit(2)

if __name__=='__main__':main()
