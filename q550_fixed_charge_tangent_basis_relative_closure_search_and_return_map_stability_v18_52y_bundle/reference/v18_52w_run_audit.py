#!/usr/bin/env python3
"""Q550 v18.52w. Authentic zero-drive posture and one-turn relative recurrence.

No new Q550 orbit integrations. Early 10D states come from authentic nested v18.51r
midpoint return protocol, late 20-series states from authentic v18.52l.
Never identify an externally imposed phase halt with a halt in material motion.
No limit cycle or Hamiltonian/global-energy claim from recurrence scouting.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,math,sys,traceback,zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from energy_math import Params,energy_terms,sigmoid
ROOT=Path(__file__).resolve().parents[1]; BASE=ROOT.parent; OUT=ROOT/'results'
C=json.loads((ROOT/'inputs_config.json').read_text()); P=Params(**C['model_params'])
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
R51='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
L52='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
SEEDS=('original','autonomous')

def sha(b):return hashlib.sha256(b).hexdigest()
def writejson(name,obj):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,allow_nan=False,ensure_ascii=False)+'\n')
def csvwrite(name,rows):
 if not rows:return
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
def loadmodule(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def unwrap_angle(d):return math.atan2(math.sin(float(d)),math.cos(float(d)))
def unit_angle(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float);na=float(np.linalg.norm(a));nb=float(np.linalg.norm(b))
 return math.degrees(math.acos(float(np.clip(a@b/(na*nb),-1.,1.)))) if min(na,nb)>1e-13 else None

def project_Q_and_tangent(y,F):
 """5D physical embedding Q=(r1cos theta1,r1sin theta1,r2cos theta2,r2sin theta2,x).
 F is dY/dxi from source RHS; projected tangent is dQ/dxi.
 """
 r1=math.exp(float(y[0]));r2=math.exp(float(y[4]));x=float(sigmoid(y[8]));a,b=float(y[2]),float(y[6])
 q=np.array([r1*math.cos(a),r1*math.sin(a),r2*math.cos(b),r2*math.sin(b),x])
 u=np.array([r1*(F[0]*math.cos(a)-F[2]*math.sin(a)),r1*(F[0]*math.sin(a)+F[2]*math.cos(a)),
             r2*(F[4]*math.cos(b)-F[6]*math.sin(b)),r2*(F[4]*math.sin(b)+F[6]*math.cos(b)),
             x*(1-x)*F[8]])
 return q,u

def extract_early(parent):
 assert sha(parent.read_bytes())==C['v1851s_sha256'],'v18.51s SHA256 mismatch'
 with zipfile.ZipFile(parent) as z:
  assert z.testzip() is None,'v18.51s zip CRC failure'
  nested=z.read(S51+'/inputs/'+R51+'_server_sendback.zip')
 assert sha(nested)=='542fdb636e8e5722cbe5ffa42ca78dad4c6e9aa4cb8c0bd19742b29ca4bc10d2','v18.51r nested SHA256 mismatch'
 with zipfile.ZipFile(io.BytesIO(nested)) as z:
  assert z.testzip() is None,'v18.51r zip CRC failure'
  rep=json.loads(z.read(R51+'/results/report.json'))
  assert rep['version']=='v18.51r' and all(rep['checks'].values()),'v18.51r authentication failed'
  blob=z.read(R51+'/results/extended_paired_kappa_nonlinear_states.npz')
 with np.load(io.BytesIO(blob),allow_pickle=False) as f:
  xi=f['xi'].copy();ys={seed:f[f'{seed}_midpoint_Y'].copy() for seed in SEEDS}
  kp={seed:float(f['kappa_'+seed]) for seed in SEEDS}
 assert len(xi)==513 and abs(xi[0])<1e-13 and abs(xi[-1]-.2)<1e-12
 assert all(y.shape==(513,10) for y in ys.values())
 km=.5*(kp['original']+kp['autonomous'])
 assert abs(km-C['v1852l_fixed_model_kappa'])<1e-12
 return xi,ys,kp,km,{'parent_sha256':sha(parent.read_bytes()),'nested_sha256':sha(nested),'nested_states_sha256':sha(blob),'branch':'v18.51r MIDPOINT RETURN (NOT FORWARD)','early_xi':[float(xi[0]),float(xi[-1])]}

def stage_engine(parent,xi):
 old=loadmodule(ROOT/'reference/v18_52v_r1_run_audit.py','q550_v1852w_parent_loader')
 # Original authentication code, with native kappa checked against the nested ancestor.
 return old.stage_parent(parent,xi)

def sample_source(engine,ctx,y,k,xi):
 t=float(ctx['event']+xi);des=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
 s=float(engine.clock.cubic_clock_sign(t,des));th=float(ctx['theta'](t))
 z=np.r_[np.asarray(y,float),th]
 F=np.asarray(engine.physical.rhs_from_original(z,k,engine.params,engine.slope,s,engine.support,engine.audit),float)
 if F.shape!=(11,):raise RuntimeError('SOURCE_RHS_UNEXPECTED_SHAPE '+str(F.shape))
 q,u=project_Q_and_tangent(y,F)
 dx=energy_terms(z,k,P,phase=th)
 work= k*float(dx['x'])*float(dx['r1'])*float(dx['r2'])*math.sin(float(dx['D']))*(s*float(dx['x']))
 return {'xi':float(xi),'t':t,'sign':s,'theta':th,'F':F,'Q5':q,'V5':u,
  'speed10':float(np.linalg.norm(F[:10])),'speed5':float(np.linalg.norm(u)),
  'phase_rate_xi':float(F[10]),'E':float(dx['E']),'D_rho':float(dx['D_rho']),
  'D_x':float(dx['D_x']),'W_x':float(dx['W_x']),'W_theta':work,
  'Jdiff':float(dx['Jdiff']),'Jsum':float(dx['Jsum']),
  'alpha_rate_xi':float(F[2]+F[6]),'gauge_rate_xi':float((F[2]-F[6])/2),
  'r1':float(dx['r1']),'r2':float(dx['r2']),'x':float(dx['x'])}

def zero_crossing(xi,early_ys,k,engine,contexts):
 rows=[];events={};deriv_rows=[]
 des=engine.original_ramp.Design('v18_50p','beta:2:2',.24)
 signfun=lambda t:float(engine.clock.cubic_clock_sign(float(t),des))
 zt=brentq(signfun,.04,.20,xtol=1e-14)
 if abs(zt-.12)>2e-9:raise RuntimeError('EXPECTED_DESIGN_ZERO_DRIFT')
 for seed in SEEDS:
  ctx=contexts[seed];root_xi=zt-float(ctx['event'])
  if not (xi[8]<root_xi<xi[-9]):raise RuntimeError('DRIVE_ZERO_NOT_BRACKETED_BY_EARLY_ARCHIVE '+seed)
  interp=PchipInterpolator(xi,early_ys[seed],axis=0)
  at=lambda d:sample_source(engine,ctx,interp(root_xi+d),k,root_xi+d)
  middle=at(0.)
  # Original 51r orbit velocities are checked against source RHS, without refitting.
  dY=interp.derivative()(root_xi)
  rel=float(np.linalg.norm(dY-middle['F'][:10])/max(1e-12,np.linalg.norm(dY),np.linalg.norm(middle['F'][:10])))
  # Whether the zero of the drive is also a stoppage of physical motion is a data question.
  halfwidths=C['drive_zero_window_halfwidths_xi']
  for h in halfwidths:
   a=at(-h);b=at(h)
   # Real energy transfer across the zero on the actual archived return path.
   grid=np.linspace(root_xi-h,root_xi+h,129)
   sample=[at(float(v-root_xi)) for v in grid]
   flux=np.array([(-q['D_rho']-q['D_x']+q['W_x']+q['W_theta'])*float(engine.omega)/q['x'] for q in sample])
   phase_flux=np.array([q['W_theta']*float(engine.omega)/q['x'] for q in sample])
   diss_flux=np.array([(q['D_rho']+q['D_x'])*float(engine.omega)/q['x'] for q in sample])
   deltaE=float(b['E']-a['E'])
   integrated=float(np.trapezoid(flux,grid))
   rows.append({'seed':seed,'halfwidth_xi':h,'zero_xi':root_xi,'zero_source_t':zt,
    'actual_delta_E_across_zero':deltaE,'signed_work_and_dissipation_predicted_deltaE':integrated,
    'energy_window_closure_error':deltaE-integrated,
    'integrated_signed_drive_work_across_zero':float(np.trapezoid(phase_flux,grid)),
    'integrated_dissipation_across_zero':float(np.trapezoid(diss_flux,grid)),
    'sign_before':a['sign'],'sign_at_zero':middle['sign'],'sign_after':b['sign'],
    'phase_rate_at_zero_xi':middle['phase_rate_xi'],'RHS10_speed_at_zero':middle['speed10'],
    'physical_Q5_speed_at_zero':middle['speed5'],'RHS10_speed_before':a['speed10'],
    'RHS10_speed_after':b['speed10'],'Q5_speed_before':a['speed5'],'Q5_speed_after':b['speed5'],
    'Q5_posture_angle_before_after_deg':unit_angle(a['V5'],b['V5']),
    'Q5_posture_angle_before_zero_deg':unit_angle(a['V5'],middle['V5']),
    'Q5_posture_angle_zero_after_deg':unit_angle(middle['V5'],b['V5']),
    'Jdiff_change_before_after':b['Jdiff']-a['Jdiff'],
    'alpha_rate_zero':middle['alpha_rate_xi'],'counterrotation_gauge_rate_zero':middle['gauge_rate_xi'],
    'work_drive_power_at_zero':middle['W_theta'],'x_source_work_at_zero':middle['W_x'],
    'dissipation_total_at_zero':middle['D_rho']+middle['D_x'],
    'source_vs_interpolated_early_Y_tangent_relative':rel,
    'interpretation':'ZERO_OF_EXTERNALLY_PRESCRIBED_DRIVE; NOT AUTOMATIC MATERIAL_VELOCITY_ZERO'})
  events[seed]={'source_design_zero_t':float(zt),'context_event':float(ctx['event']),
    'early_xi_crossing':float(root_xi),'native_midpoint_kappa':float(k),
    'source_rhs_zero_sign':middle['sign'],'theta_xi_rate_at_zero':middle['phase_rate_xi'],
    'physical_Q5_speed_xi_at_zero':middle['speed5'],'physical_RHS10_norm_xi_at_zero':middle['speed10'],
    'early_Y_derivative_vs_source_rhs_relative':rel,'Jdiff_at_zero':middle['Jdiff'],
    'Wtheta_at_zero':middle['W_theta'],'D_total_at_zero':middle['D_rho']+middle['D_x'],
    'W_x_at_zero':middle['W_x'],'alpha_rate_at_zero':middle['alpha_rate_xi'],
    'counterrotation_gauge_rate_at_zero':middle['gauge_rate_xi'],
    'source_phase_derivative_minus_omega_sign_abs':abs((ctx['theta'](zt+1e-6)-ctx['theta'](zt-1e-6))/(2e-6)-float(engine.omega)*middle['sign'])}
  grid=np.linspace(root_xi-.025,root_xi+.025,201)
  for v in grid:
   o=sample_source(engine,ctx,interp(v),k,v)
   deriv_rows.append({'seed':seed,'xi_minus_zero':float(v-root_xi),'xi':float(v),
    'source_clock_sign':o['sign'],'external_phase_rate_support_xi':o['phase_rate_xi'],
    'material_RHS10_speed':o['speed10'],'Q5_tangent_speed':o['speed5'],
    'RHS_sum_angles_rate':o['alpha_rate_xi'],'RHS_gauge_rate':o['gauge_rate_xi'],
    'drive_power_physical_time':o['W_theta'],'dissipation_total':o['D_rho']+o['D_x'],
    'x_source_work':o['W_x']})
 return rows,events,deriv_rows

def interp_10(xi,Y,x):
 return np.array([np.interp(float(x),xi,Y[:,j]) for j in range(10)],float)

def full_one_turn_relative_map(xi,series,phase):
 rows=[]; anchor_n=int(C['recurrence_anchors']); ntwopi=2*math.pi
 for name,Y in series.items():
  seed=name.split('_',1)[0]
  if seed not in SEEDS:raise RuntimeError('BAD_NATIVE_SERIES '+name)
  p=phase[seed]
  if not (len(Y)==len(xi)==len(p) and np.all(np.diff(p)<0)):
   raise RuntimeError('LATE_PHASE_NOT_MONOTONIC '+name)
  total=float(p[0]-p[-1]);
  if total<ntwopi:continue
  pa=np.linspace(float(p[0]),float(p[-1]+ntwopi),anchor_n)
  for j,phase_a in enumerate(pa):
   xa=float(np.interp(phase_a,p[::-1],xi[::-1])); xb=float(np.interp(phase_a-ntwopi,p[::-1],xi[::-1]))
   a=interp_10(xi,Y,xa);b=interp_10(xi,Y,xb)
   eight=(0,1,3,4,5,7,8,9)
   d8=np.array([b[t]-a[t] for t in eight]);err_alpha=unwrap_angle((b[2]+b[6])-(a[2]+a[6]))
   quot=float(np.hypot(np.linalg.norm(d8),err_alpha));scaled=quot/max(1.,np.linalg.norm(a[list(eight)]))
   da=unwrap_angle(b[2]-a[2]);db=unwrap_angle(b[6]-a[6]);gauge_shift=unwrap_angle(.5*(da-db))
   diffs={'logr1':b[0]-a[0],'v1':b[1]-a[1],'omega1':b[3]-a[3],
     'logr2':b[4]-a[4],'v2':b[5]-a[5],'omega2':b[7]-a[7],
     'z':b[8]-a[8],'vz':b[9]-a[9]}
   pwr1=math.exp(2*a[0])*a[3]-math.exp(2*a[4])*a[7]
   pwr2=math.exp(2*b[0])*b[3]-math.exp(2*b[4])*b[7]
   rows.append({'series':name,'seed':seed,'anchor_index':j,'xi_start':xa,'xi_after_phase_turn':xb,
     'support_spacing':xb-xa,'material_time_spacing':float(np.interp(xb,xi,Y[:,10])-np.interp(xa,xi,Y[:,10])),
     'phase_winding':-1,'phase_match_error_rad':float(np.interp(xb,xi,p)-np.interp(xa,xi,p)+ntwopi),
     'quotient_error_absolute_mixed_units':quot,'quotient_error_scaled_descriptive':scaled,
     'alpha_sum_phase_error_rad':err_alpha,'gauge_angle_shift_rad':gauge_shift,
     'Jdiff_change':pwr2-pwr1,'max_nonangular_component_absolute_difference':float(np.max(np.abs(d8))),
     **{'difference_'+k:float(v) for k,v in diffs.items()},
     'interpretation':'ONE_EXTERNAL_PHASE_TURN_QUOTIENT_RECURRENCE_SCOUT_NOT_POINCARE_STABILITY_CERTIFICATE'})
 return rows

def plot_data(deriv,cycle):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except Exception:return []
 paths=[]
 for seed in SEEDS:
  rows=[r for r in deriv if r['seed']==seed]
  if rows:
   t=np.array([r['xi_minus_zero'] for r in rows])
   fig,ax=plt.subplots(3,1,figsize=(9,8),sharex=True,constrained_layout=True)
   ax[0].plot(t,[r['source_clock_sign'] for r in rows],label='prescribed signed drive');ax[0].axhline(0,c='k',lw=.6)
   ax[1].plot(t,[r['Q5_tangent_speed'] for r in rows],label='physical 5D tangent speed')
   ax[1].plot(t,[r['material_RHS10_speed'] for r in rows],label='original 10D RHS norm',alpha=.7)
   ax[2].plot(t,[r['drive_power_physical_time'] for r in rows],label='phase work power')
   ax[2].plot(t,[r['x_source_work'] for r in rows],label='x input work',alpha=.7)
   for a in ax:a.axvline(0,c='gray',ls='--',lw=.8);a.grid(alpha=.25);a.legend(fontsize=8)
   ax[2].set_xlabel('Support xi relative to actual drive-zero event')
   fig.suptitle('Q550 v18.52w '+seed+' — drive zero vs physical material velocity')
   p=OUT/(seed+'_drive_zero_posture.png');fig.savefig(p,dpi=135);plt.close(fig);paths.append(p.name)
  rows=[r for r in cycle if r['seed']==seed and r['series']==seed+'__reference']
  if rows:
   fig,ax=plt.subplots(2,1,figsize=(9,7),sharex=True,constrained_layout=True)
   t=[r['xi_start'] for r in rows];ax[0].plot(t,[r['quotient_error_absolute_mixed_units'] for r in rows],'o-',label='full quotient return residual')
   ax[1].plot(t,[abs(r['alpha_sum_phase_error_rad']) for r in rows],'o-',label='|angle sum mismatch|')
   for a in ax:a.grid(alpha=.25);a.legend(fontsize=8)
   ax[1].set_xlabel('Start xi for one external-phase turn')
   fig.suptitle('Q550 v18.52w '+seed+' — one-turn scout, NOT cycle-limit certificate')
   p=OUT/(seed+'_relative_phase_turn_recurrence.png');fig.savefig(p,dpi=135);plt.close(fig);paths.append(p.name)
 return paths

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--native-v1852l',type=Path,default=BASE/L52/'results/extended_horizon_10d_states.npz')
 ap.add_argument('--require-parent',action='store_true');ap.add_argument('--require-native',action='store_true');args=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 par_candidates=[args.parent_v1851s,BASE/(S51+'_server_sendback.zip'),Path('/root/elab1/hysteresis')/(S51+'_server_sendback.zip'),Path('/root/Downloads')/(S51+'_server_sendback.zip')]
 parent=next((x for x in par_candidates if x is not None and Path(x).is_file()),None)
 report={'version':'v18.52w','title':C['title'],'historical_parent_available':parent is not None,
    'native_late_available':args.native_v1852l.is_file(),'finite_late_archive_starts_after_drive_zero':True,
    'new_native_IVP_runs':0,'observed_full_forward_return_orbit_pair':False,
    'limit_cycle_certified':False,'relative_periodic_orbit_certified':False,
    'conservative_full_dynamics_certified':False,'global_energy_bound_certified':False}
 if args.require_parent and parent is None:raise RuntimeError('REQUIRED_V18_51S_PARENT_MISSING')
 if args.require_native and not args.native_v1852l.is_file():raise RuntimeError('REQUIRED_V18_52L_NATIVE_MISSING')
 if parent is None:
  report['scientific_decision']='PARTIAL_LOCAL_PROTOCOL_TESTS_ONLY_ANCESTOR_MISSING';writejson('report.json',report);return
 temp=None
 try:
  early_xi,early_y,ks,k,origin=extract_early(Path(parent))
  temp,engine,_,native_k,ctx,sources=stage_engine(Path(parent),early_xi)
  if any(abs(v-k)>1e-12 for v in native_k.values()):raise RuntimeError('MIDPOINT_KAPPA_MISMATCH')
  # Recover source code directly into the genealogy, rather than trusting a verbal sign summary.
  source={'clock_cubic_sign':sources['cubic_clock_sign'],'clock_excursion':sources['cubic_phase_excursion'],
          'make_context':sources['make_context'],'rhs_from_original':sources['rhs_from_original'],
          'original_integrate':sources['original_integrate'],'parent_authentication':sources['authenticated_parent']}
  writejson('authenticated_source_genealogy.json',source)
  report['early_lineage']=origin
  events,evrep,curve=zero_crossing(early_xi,early_y,k,engine,ctx)
  csvwrite('actual_drive_zero_posture_events.csv',events);csvwrite('drive_zero_local_curves.csv',curve)
  report['actual_drive_zero_events']=evrep
  report['phase_derivative_matches_source']=all(q['source_phase_derivative_minus_omega_sign_abs']<2e-6 for q in evrep.values())
  # The native 52l states have no drive-zero crossing. They are ONLY used for relative recurrence.
  cycle=[];xi_range=None
  if args.native_v1852l.is_file():
   with args.native_v1852l.open('rb') as f:
    h=hashlib.sha256()
    for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
   if h.hexdigest()!=C['v1852l_10d_sha256']:raise RuntimeError('LATE_NATIVE_V18_52L_SHA256_MISMATCH')
   with np.load(args.native_v1852l,allow_pickle=False) as z:
    late_xi=np.array(z['common_xi'],float)
    names=[n for n in z.files if n.endswith('__Y_tau')]
    series={n.removesuffix('__Y_tau'):np.array(z[n],float) for n in names}
   if not(len(series)==20 and late_xi[0]>1 and late_xi[-1]>late_xi[0]):raise RuntimeError('LATE_NATIVE_SHAPE_OR_HORIZON_MISMATCH')
   late_phase={seed:np.array([ctx[seed]['theta'](ctx[seed]['event']+float(v)) for v in late_xi],float) for seed in SEEDS}
   cycle=full_one_turn_relative_map(late_xi,series,late_phase)
   csvwrite('one_turn_full_quotient_recurrence_scout.csv',cycle)
   xi_range=[float(late_xi[0]),float(late_xi[-1])]
  report['native_late_xi_range']=xi_range;report['one_turn_rows']=len(cycle)
  report['one_turn_min_scaled_quotient_residual']=float(min((v['quotient_error_scaled_descriptive'] for v in cycle),default=math.inf)) if cycle else None
  report['one_turn_max_scaled_quotient_residual']=float(max((v['quotient_error_scaled_descriptive'] for v in cycle),default=-math.inf)) if cycle else None
  report['plots']=plot_data(curve,cycle)
  # Gate source and physical continuity are descriptive; no arbitrary 'limit-cycle PASS'.
  report['physical_motion_stopped_at_drive_zero']=all(v['physical_Q5_speed_xi_at_zero']<1e-12 for v in evrep.values())
  report['maximum_source_vs_recorded_early_tangent_relative']=max(q['early_Y_derivative_vs_source_rhs_relative'] for q in evrep.values())
  report['scientific_decision']='AUTHENTIC_ZERO_DRIVE_POSTURE_AND_ONE_TURN_QUOTIENT_RECURRENCE_METROLOGY__NO_LIMIT_CYCLE_CERTIFICATE' if cycle else 'AUTHENTIC_ZERO_DRIVE_POSTURE_METROLOGY__LATE_RECURRENCE_NOT_AVAILABLE'
  writejson('report.json',report)
  text=['# Q550 v18.52w — authenticated server result','',
     '## Status',report['scientific_decision'],'',
     '## What is authenticated','- SHA256-pinned original v18.51s ancestry and nested v18.51r early midpoint 10D trajectories.',
     '- Event-specific signed clock derivative and 10D physical RHS at the actual drive-zero point.',
     '- The distinct late 52l native state horizon for one-turn quotient return scouting.',
     '', '## Drive zero events']
  for seed,v in evrep.items():text.extend([f'- {seed}: xi_zero={v["early_xi_crossing"]:.12g}, Q5 speed={v["physical_Q5_speed_xi_at_zero"]:.12g}, RHS10 norm={v["physical_RHS10_norm_xi_at_zero"]:.12g}, source/trace derivative relative={v["early_Y_derivative_vs_source_rhs_relative"]:.6g}.'])
  text.extend(['','## No hidden promotion to theorem',
     'A nonzero material tangent at one drive-zero event does not establish an attracting limit cycle.',
     'A one-turn phase-matched near return is not a full return map, Floquet stability, proof of dissipative balance, or global energy conservation.',
     'The 52l horizon starts at xi=12.8 and must not be described as an observation of the early drive zero.'])
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(text)+'\n')
 finally:
  if temp is not None:temp.cleanup()

if __name__=='__main__':
 try:main()
 except Exception as e:
  OUT.mkdir(parents=True,exist_ok=True)
  (OUT/'FAILURE_DIAGNOSTIC.txt').write_text(traceback.format_exc())
  print('v18.52w FAILURE',repr(e),flush=True)
  raise
