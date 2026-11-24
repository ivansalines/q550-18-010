#!/usr/bin/env python3
"""Q550 v18.52t: observed external clock, signed-drive and moving-frame audit.

No IVP integration. The 11th column of the v18.52l NPZ is physical material
elapsed time, NOT the Theta variable of the eleven-component original RHS.
Model source and native vector hashes must match authenticated ancestors.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,math,sys,tempfile,zipfile
from pathlib import Path
import numpy as np
from energy_math import Params,energy_terms,sigmoid
from phase_energy import (integral,signed_power,frame_quantity,frame_ideal_power,
                          phase_work_per_unit,summary_distribution,physical_rhs_with_signed_external_phase)
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent
CONFIG=json.loads((ROOT/'inputs_config.json').read_text());P=Params(**CONFIG['model_params'])
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
R51='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
L52='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
DEFAULT_STATES=BASE/L52/'results/extended_horizon_10d_states.npz'
DEFAULT_PARENT=BASE/(S51+'_server_sendback.zip')


def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def dynamic_import(path,name):
 s=importlib.util.spec_from_file_location(name,path)
 m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def csvwrite(path,rows):
 if not rows:return
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def authenticate_parent(path,xi):
 if sha(path)!=CONFIG['v1851s_sha256']:raise RuntimeError('v18.51s SHA256 mismatch: '+sha(path))
 with zipfile.ZipFile(path) as z:
  if z.testzip() is not None:raise RuntimeError('v18.51s zip CRC failure')
  rblob=z.read(S51+'/inputs/'+R51+'_server_sendback.zip')
  run_source=z.read(S51+'/scripts/run_audit.py')
 with zipfile.ZipFile(io.BytesIO(rblob)) as z:
  if z.testzip() is not None:raise RuntimeError('nested v18.51r CRC failure')
  a=z.read(R51+'/results/extended_paired_kappa_nonlinear_states.npz')
 with np.load(io.BytesIO(a),allow_pickle=False) as nz:
  ko,ka=float(nz['kappa_original']),float(nz['kappa_autonomous'])
 k=(ko+ka)/2
 if abs(k-CONFIG['v1852l_fixed_model_kappa'])>1e-12:
  raise RuntimeError('Lineage kappa mismatch '+str(k))
 parent=dynamic_import(ROOT/'reference/parent_v18_51t_run_audit.py','q550_v1852t_staging_parent')
 with tempfile.TemporaryDirectory(prefix='q550_v1852t_staging_') as d:
  parts=parent.stage_engine(Path(d),rblob,run_source)[4]
  _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
  params={name:float(getattr(engine.params,name)) for name in CONFIG['model_params']}
  mismatch={n:params[n]-float(p) for n,p in CONFIG['model_params'].items()}
  (ROOT/'results/parameter_audit.json').write_text(json.dumps({'status':'MATCH' if all(abs(v)<1e-12 for v in mismatch.values()) else 'MISMATCH',
      'authenticated_parameters':params,'expected_parameters':CONFIG['model_params'],'differences':mismatch},indent=2)+'\n')
  if any(abs(v)>=1e-12 for v in mismatch.values()):raise RuntimeError('AUTHENTICATED_MODEL_PARAMETER_MISMATCH '+str(mismatch))
  phase={};ctxs={}
  for seed in ('original','autonomous'):
   ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
   phase[seed]=np.array([float(ctx['theta'](float(ctx['event'])+float(t))) for t in xi])
   ctxs[seed]=dict(event=float(ctx['event']))
  # Independent original source: report what its 11th component actually asserts,
  # never overwrite it with the externally reconstructed phase rate.
  klist=[ko,ka]
  first=[]
  for seed,ks in zip(('original','autonomous'),klist):
   ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
   for ii in np.linspace(1,len(xi)-2,16,dtype=int):
    th=float(phase[seed][ii]);z=np.zeros(11)
    z[[0,4]]=[.3,-.8];z[[1,5]]=[.01,-.03]
    z[[2,6]]=[.2,-.7];z[[3,7]]=[-.4,.2];z[8]=0.;z[9]=.02;z[10]=th
    f=np.asarray(engine.jet.structural_jet(z,ks,engine.params,engine.slope)[0],float)
    x=float(sigmoid(z[8]))
    normalized=f*x/(2*math.pi*float(engine.slope))
    first.append({'seed':seed,'x':x,'normalized_engine_theta_rate':float(normalized[10]),
       'engine_theta_rate_over_x':float(normalized[10]/x)})
  # Cross-check ALL 11 components against the historical physical RHS after
  # removal of engine.slope and support-clock factor; independent of native data.
  engine_errors=[]; engine_last=[]
  rng=np.random.default_rng(18525217)
  for seed,ks in zip(('original','autonomous'),klist):
   for i in range(48):
    y=np.zeros(11);y[[0,4]]=rng.uniform(-1,1,2)
    y[[1,3,5,7,9]]=rng.uniform(-.8,.8,5)
    y[[2,6,10]]=rng.uniform(-4,4,3)
    y[8]=rng.uniform(-2,2)
    f=np.asarray(engine.jet.structural_jet(y,ks,engine.params,engine.slope)[0],float)
    x=float(sigmoid(y[8]))
    normed=f*x/(2*math.pi*float(engine.slope))
    archived=physical_rhs_with_signed_external_phase(y,ks,P,y[10],sigma=+1)
    engine_errors.append(float(np.linalg.norm(normed-archived)/max(1.,np.linalg.norm(archived))))
    engine_last.append(float(normed[10]/x))
  # Evaluate original off-orbit generator at genuine native coordinates outside
  # temporary dir, while its loaded code may rely on its staged source path.
  def engine_probe(y,th,ks):
   z=np.r_[np.asarray(y[:10],float),float(th)]
   f=np.asarray(engine.jet.structural_jet(z,ks,engine.params,engine.slope)[0],float)
   x=float(sigmoid(z[8]))
   return f*x/(2*math.pi*float(engine.slope))
  # Need all probes while temporary genealogy is staged.
  origin_meta={'v1851s_sha256':sha(path),'kappa_original':ko,'kappa_autonomous':ka,'kappa_fixed':k,
    'engine_slope':float(engine.slope),'clock_11th_component_probes':first,
    'engine_versus_historical_normalized_max_relative':max(engine_errors),
    'engine_last_rate_over_x_range':[min(engine_last),max(engine_last)],
    'source':'v18.51s staged authenticated engine'}
  return phase,k,origin_meta,engine_probe

# Staged generator function must not be called after its TemporaryDirectory is removed.
# Instead, a controlled off-orbit source sample is recorded in authenticate_parent.

def trace_metrics(Y,phi,k,stride,name,seed):
 tau=np.asarray(Y[:,10],float)
 if np.any(np.diff(tau)<=0) or not np.all(np.isfinite(tau)):raise RuntimeError('Native tau invalid '+name)
 ids=np.arange(0,len(tau),stride)
 if ids[-1]!=len(tau)-1:ids=np.r_[ids,len(tau)-1]
 y=np.asarray(Y[ids],float);T=tau[ids];PH=np.asarray(phi[ids],float)
 # theta function is shared by all kicks in each seed; kick x may differ from reference x.
 x=np.asarray(sigmoid(y[:,8]),float);thdot=np.gradient(PH,T,edge_order=2)
 d=energy_terms(y,k,P,phase=PH)
 E=np.asarray(d['E'],float);c=np.asarray(phase_work_per_unit(y,k,P,PH),float)
 diss=np.asarray(d['D_rho']+d['D_x']);wx=np.asarray(d['W_x']);J=np.asarray(d['Jsum'])
 exact_flux=-diss+wx+c*thdot
 E_delta=E-E[0]
 # Independent angular-equation replay from the SAVED trajectory derivatives.
 # Use saved sample time (do not match array rows across distinct kicks).
 dw1=np.gradient(np.asarray(Y[:,3],float),tau,edge_order=2)[ids]
 dw2=np.gradient(np.asarray(Y[:,7],float),tau,edge_order=2)[ids]
 rr1,rr2=np.asarray(d['r1']),np.asarray(d['r2'])
 sinD=np.sin(d['D'])
 f1=k*x*(rr2/rr1)*sinD-2*y[:,1]*y[:,3]
 f2=k*x*(rr1/rr2)*sinD-2*y[:,5]*y[:,7]
 angular_rms_rel_1=float(np.linalg.norm(dw1-f1)/max(1e-10,np.linalg.norm(dw1)))
 angular_rms_rel_2=float(np.linalg.norm(dw2-f2)/max(1e-10,np.linalg.norm(dw2)))
 abs_closure=E_delta-integral(T,exact_flux)
 terms={}
 for sigma in (-1,1):
  es=(-diss+wx+c*sigma*x)
  Eledger=integral(T,es)
  H=np.asarray(frame_quantity(y,k,P,PH,sigma),float)
  Hledger=integral(T,np.asarray(frame_ideal_power(y,k,P,PH,sigma),float))
  ideal_mismatch=c*(thdot-sigma*x)
  terms['minus' if sigma==-1 else 'plus']={
    'energy_endpoint_error':float(E_delta[-1]-Eledger[-1]),
    'energy_max_abs_cumulative_error':float(np.max(np.abs(E_delta-Eledger))),
    'energy_integrated_absolute_mismatch':float(np.trapezoid(np.abs(ideal_mismatch),T)),
    'H_endpoint_error':float(H[-1]-H[0]-Hledger[-1]),
    'H_max_abs_cumulative_error':float(np.max(np.abs(H-H[0]-Hledger))),
    'H_total_change':float(H[-1]-H[0]),
    'clock_rms_rate_error':float(np.sqrt(np.mean((thdot-sigma*x)**2))),
    'clock_integrated_signed_mismatch':float(PH[-1]-PH[0]-sigma*np.trapezoid(x,T))}
 # For reference traces, higher fidelity clock test is integral of x at full grid.
 row={'series':name,'seed':seed,'n_samples':int(len(y)),'tau_duration':float(T[-1]-T[0]),
    'theta_delta':float(PH[-1]-PH[0]),'x_integral':float(np.trapezoid(x,T)),
    'reference_clock_preferred_sign':int(-1 if terms['minus']['clock_rms_rate_error']<terms['plus']['clock_rms_rate_error'] else 1),
    'native_angular1_rhs_rms_rel':angular_rms_rel_1,
    'native_angular2_rhs_rms_rel':angular_rms_rel_2,
    'arbitrary_rate_E_max_abs_closure':float(np.max(np.abs(abs_closure))),
    'arbitrary_rate_E_endpoint_closure':float(abs_closure[-1]),
    'clock_external_rate_minus_x_rms':terms['minus']['clock_rms_rate_error'],
    'clock_external_rate_plus_x_rms':terms['plus']['clock_rms_rate_error'],
    'plus_E_max_error':terms['plus']['energy_max_abs_cumulative_error'],
    'minus_E_max_error':terms['minus']['energy_max_abs_cumulative_error'],
    'plus_H_max_error':terms['plus']['H_max_abs_cumulative_error'],
    'minus_H_max_error':terms['minus']['H_max_abs_cumulative_error'],
    'signed_plus_E_power_integrated_abs_error':terms['plus']['energy_integrated_absolute_mismatch'],
    'signed_minus_E_power_integrated_abs_error':terms['minus']['energy_integrated_absolute_mismatch'],
    'native_Theta_is_external':True,'global_time_uniform_bound_certified':False}
 return row,dict(T=T,E=E,DE=E_delta,theta=PH,x=x,rate=thdot,
    Eledger_actual=integral(T,exact_flux),d=d,c=c,terms=terms,J=J,wd=wx-diss)


def plots(ref):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return []
 out=[]
 for seed,rr in ref.items():
  T=rr['T'];d=rr['d'];X=rr['x'];D=rr['DE'];rate=rr['rate'];J=rr['J'];c=rr['c']
  fig,axes=plt.subplots(3,1,figsize=(11,10),sharex=True,constrained_layout=True)
  axes[0].plot(T,rate,label='Authenticated external dTheta/dtau',lw=1.1)
  axes[0].plot(T,-X,'--',label='-x',lw=1)
  axes[0].plot(T,+X,':',label='+x',lw=1)
  axes[0].set_ylabel('Clock rate');axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
  axes[1].plot(T,D,label='Observed delta E',lw=1.2)
  for sig in (-1,1):
   es=-d['D_rho']-d['D_x']+d['W_x']+c*sig*X
   axes[1].plot(T,integral(T,es),'--' if sig==-1 else ':',label=('Signed -x' if sig==-1 else 'Signed +x')+' budget')
  axes[1].set_ylabel('Energy change');axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
  for sig in (-1,1):
   H=rr['E']-sig*.5*X*J
   # H ledger uses already computed model energy terms, no artificial reconstruction.
   hb=integral(T,-d['D_rho']-d['D_x']+d['W_x']-sig*.5*d['dx']*J)
   axes[2].plot(T,H-H[0]-hb,label=('H- ' if sig==-1 else 'H+ ')+'frame closure')
  axes[2].set_ylabel('Rotating-frame closure');axes[2].set_xlabel('Material time tau')
  axes[2].legend(fontsize=8);axes[2].grid(alpha=.2)
  fig.suptitle('Q550 v18.52t — '+seed+' authentic reference clock, signed drive, frame energy')
  p=ROOT/'results'/f'{seed}_signed_drive_rotating_frame.png';fig.savefig(p,dpi=140);plt.close(fig);out.append(str(p.relative_to(ROOT)))
 return out


def main():
 a=argparse.ArgumentParser()
 a.add_argument('--states',type=Path,default=DEFAULT_STATES)
 a.add_argument('--parent-v1851s',type=Path)
 a.add_argument('--stride',type=int,default=1)
 a.add_argument('--require-original-engine',action='store_true')
 a.add_argument('--require-native',action='store_true')
 args=a.parse_args()
 if args.stride<1:raise ValueError('--stride must be >=1')
 OUT=ROOT/'results';OUT.mkdir(parents=True,exist_ok=True)
 status={'native':'UNAVAILABLE','original_engine':'UNAVAILABLE','scientific_decision':'NOT_EVALUATED'}
 filename=args.states;parent=args.parent_v1851s or next((p for p in [DEFAULT_PARENT,ROOT/'inputs'/(S51+'_server_sendback.zip'),Path('/root/elab1/hysteresis')/(S51+'_server_sendback.zip'),Path('/root/Downloads')/(S51+'_server_sendback.zip')] if p.is_file()),None)
 if args.require_native and not filename.is_file():raise RuntimeError('AUTHENTIC_NATIVE_STATES_MISSING '+str(filename))
 if args.require_original_engine and (parent is None or not parent.is_file()):raise RuntimeError('AUTHENTIC_ENGINE_PARENT_MISSING '+str(parent))
 hist_source=ROOT/'reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py'
 # Strong analytical probe independent of native genealogy.
 from phase_energy import exact_identity_probe
 rng=np.random.default_rng(18525203);analytic=[]
 for i in range(160):
  q=np.zeros(11);q[[0,4]]=rng.uniform(-1.5,1.5,2);q[[1,3,5,7,9]]=rng.uniform(-.8,.8,5)
  q[[2,6]]=rng.uniform(-3,3,2);q[8]=rng.uniform(-2,2);phase=float(rng.uniform(-5,5))
  analytic.append(abs(exact_identity_probe(q,float(rng.uniform(-.3,.3)),P,phase,float(rng.uniform(-2,2)))))
 status['arbitrary_external_rate_exact_identity_max_numerical_residual']=max(analytic)
 rows=[];ref={};details={}
 if filename.is_file():
  if sha(filename)!=CONFIG['v1852l_10d_sha256']:raise RuntimeError('NATIVE_STATES_SHA_MISMATCH '+sha(filename))
  status['native']='AUTHENTIC_SHA256_VERIFIED'
  if parent is not None and parent.is_file():
   with np.load(filename,allow_pickle=False) as z:
    xi=np.asarray(z['common_xi'],float)
    phase,k,orig,probe=authenticate_parent(parent,xi)
    status['original_engine']='AUTHENTIC_AND_PARAMETER_MATCH'
    status['engine_metadata']=orig
    # Parent provides a single external phase function per seed across all kicks.
    for name in z.files:
     if name=='common_xi':continue
     seed='original' if name.startswith('original') else 'autonomous'
     Y=np.asarray(z[name],float)
     if Y.shape!=(len(xi),11):raise RuntimeError('Invalid native shape '+name)
     rr,diag=trace_metrics(Y,phase[seed],k,args.stride,name,seed)
     rows.append(rr)
     if name.endswith('reference__Y_tau'):ref[seed]=diag
   status['native']='AUTHENTIC_20D_SERIES_AND_EXTERNAL_PHASE' if len(rows)==20 else 'AUTHENTIC_WRONG_SERIES_COUNT'
   csvwrite(OUT/'signed_phase_energy_per_series.csv',rows)
   status['figures']=plots(ref)
   status['max_native_angular_rhs_rms_rel']=max(max(r['native_angular1_rhs_rms_rel'],r['native_angular2_rhs_rms_rel']) for r in rows)
   status['n_series']=len(rows)
   status['max_exact_E_energy_closure']=max(r['arbitrary_rate_E_max_abs_closure'] for r in rows)
   status['clock_reference_preferred_signs']={seed:next(r['reference_clock_preferred_sign'] for r in rows if r['series']==seed+'__reference__Y_tau') for seed in ('original','autonomous')}
   status['clock_reference_phase_integral_residual']={seed:next(r['theta_delta']+r['x_integral'] for r in rows if r['series']==seed+'__reference__Y_tau') for seed in ('original','autonomous')}
   status['signed_frame_reference_max_closures']={seed:{'minus':next(r['minus_H_max_error'] for r in rows if r['series']==seed+'__reference__Y_tau'),
       'plus':next(r['plus_H_max_error'] for r in rows if r['series']==seed+'__reference__Y_tau')} for seed in ('original','autonomous')}
   status['non_reference_kicks_have_external_phase_clock_inherited_from_reference']=True
   # Do not certify global uniform noncollapse from a finite trace.
   status['global_uniform_noncollapse_bound_certified']=False
  else:status['native']='AUTHENTIC_STATES_BUT_EXTERNAL_PHASE_LINEAGE_MISSING'
 else:status['native']='MISSING_SOURCE'
 certified=(status['native']=='AUTHENTIC_20D_SERIES_AND_EXTERNAL_PHASE'
     and all(x==-1 for x in status.get('clock_reference_preferred_signs',{}).values())
     and max(abs(x) for x in status.get('clock_reference_phase_integral_residual',{}).values())<1e-6
     and status['max_exact_E_energy_closure']<1e-6
    and status.get('max_native_angular_rhs_rms_rel',1)>-1 and status.get('max_native_angular_rhs_rms_rel',1)<1e-3
    and status['engine_metadata']['engine_versus_historical_normalized_max_relative']<1e-10)
 original_plus=(certified and abs(status['engine_metadata']['engine_last_rate_over_x_range'][0]-1)<1e-10
    and abs(status['engine_metadata']['engine_last_rate_over_x_range'][1]-1)<1e-10)
 status['original_engine_vs_external_clock_discordance']=bool(original_plus)
 status['scientific_decision']=('AUTHENTIC_EXTERNAL_NEGATIVE_CLOCK_WITH_INTERNAL_POSITIVE_CLOCK_DISCORDANCE_AND_SIGNED_ENERGY_LEDGER' if original_plus
    else ('AUTHENTIC_EXTERNAL_PHASE_NEGATIVE_ORIENTATION_AND_SIGNED_ENERGY_LEDGER_VERIFIED' if certified
    else 'ANALYTIC_SIGNED_IDENTITY_WITH_INCOMPLETE_OR_FAILED_AUTHENTIC_CLOCK_CHECK'))
 status['original_11D_clock_identity_claim']='NOT_CERTIFIED_FROM_EXTERNAL_CLOCK_ALONE'
 status['notes']=['Original engine field may have an independent 11th component; compare source law and external source clock before equating them.',
    'Native external phase was constructed per seed using the reference clock; it need not satisfy dTheta/dtau=+/-x_kick for each perturbation.',
    'Moving-frame H_sigma cancellation requires signed clock law; phase rotation alone is not a physical energy conservation theorem.',
    'No uniform-in-time energy or 5D separation noncollapse theorem derived.']
 (OUT/'report.json').write_text(json.dumps(status,indent=2)+'\n')
 md=['# Q550 v18.52t — Phase-Clock Orientation and Signed-Drive Certification','',
     '## Scientific decision',status['scientific_decision'],'',
     '## Exact identity for arbitrary external Theta',
     'E_dot = -gamma_rho*(r1dot^2+r2dot^2)-gamma_x*xdot^2 + r1*r2*(mu_x*sin(Delta)+kappa0*cos(Delta))*xdot + kappa0*x*r1*r2*sin(Delta)*Theta_dot.',
     'Jsum_dot=2*kappa0*x*r1*r2*sin(Delta). For sigma=+/-1 and Theta_dot=sigma*x:',
     'H_sigma=E-sigma*x*(J1+J2)/2; Hsigma_dot = -D + W_x - sigma*xdot*(J1+J2)/2.',
     '', '## Authenticated clock and native measurements','```json',json.dumps(status,indent=2),'```','',
     '## Caveats','The v18.52l NPZ 11th column is tau, NOT Theta. External Theta comes only from SHA-authenticated v18.51s contexts.',
     'An external reference clock is shared across kicks within a seed; a signed clock law is not exactly enforced on each kick x(t).',
     'Full original 11D source-clock consistency is not proven merely by native energy-closure.',
     'Finite-time bounds must not be described as uniform global bounds.',
     'No new IVP integrations are performed.']
 (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(md)+'\n')
 print(json.dumps({'decision':status['scientific_decision'],'native':status['native'],'clock':status.get('clock_reference_preferred_signs'),
     'phase_residual':status.get('clock_reference_phase_integral_residual'),
     'max_energy_closure':status.get('max_exact_E_energy_closure'),
     'engine_last_component':status.get('engine_metadata',{}).get('clock_11th_component_probes',[])[:1]},indent=2))
 if args.require_native and status['native']!='AUTHENTIC_20D_SERIES_AND_EXTERNAL_PHASE':raise SystemExit(2)
 if args.require_original_engine and not certified:raise SystemExit(3)

if __name__=='__main__':main()
