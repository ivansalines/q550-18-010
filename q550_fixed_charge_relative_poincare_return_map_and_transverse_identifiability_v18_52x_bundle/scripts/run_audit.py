#!/usr/bin/env python3
"""Q550 v18.52x — authenticated one-step relative Poincare map, fixed-charge identifiability.

The 52l library has 20 long native 10D return states, with a common external phase
per seed. Each orbit conserves its own J, but the kicks in this library change J.
Therefore the library does NOT support an 8D fixed-J transverse Floquet spectrum.
This audit computes authentic, one-step section pairs, source/data checks, charge
mismatch, and 1D family directional response only. It never adjusts stored kicks
to pretend that they were integrated at J_reference.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,sys,traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;OUT=ROOT/'results'
C=json.loads((ROOT/'inputs_config.json').read_text()); SEEDS=('original','autonomous'); TWOPI=2*math.pi
NATIVE_STEM='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
PARENT_STEM='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'


def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def modangle(a):return math.atan2(math.sin(float(a)),math.cos(float(a)))

def J(y):
 y=np.asarray(y,float)
 return float(np.exp(2*y[0])*y[3]-np.exp(2*y[4])*y[7])

def reduce8(y):
 """8 coordinates of R10 / counterrotation at fixed J.

 u=(logr1,v1,omega1,logr2,v2,z,vz,alpha=(theta1+theta2) mod 2pi).
 omega2 can be reconstructed from J and u. This is a CHART, not a Euclidean
 kinetic metric; scales are descriptive and never define a Floquet certificate.
 """
 y=np.asarray(y,float)
 return np.array([y[0],y[1],y[3],y[4],y[5],y[8],y[9],modangle(y[2]+y[6])],float)

def lift8(u,j):
 u=np.asarray(u,float)
 y=np.zeros(10,float)
 y[0],y[1],y[3],y[4],y[5],y[8],y[9]=u[:7]
 y[2]=u[7]/2;y[6]=u[7]/2
 y[7]=(math.exp(2*y[0])*y[3]-j)/math.exp(2*y[4])
 return y

def diff8(a,b):
 d=np.asarray(b)-np.asarray(a);d=d.copy();d[7]=modangle(d[7]);return d

def interpolate_state(x,y,t,method='linear'):
 x=np.asarray(x,float);idx=int(np.searchsorted(x,t))
 if idx==0 or idx==len(x):
  if abs(t-x[min(max(idx-1,0),len(x)-1)])>1e-8:raise ValueError('SECTION_OUTSIDE_NATIVE_INTERVAL')
  return np.asarray(y[min(max(idx-1,0),len(x)-1),:10],float)
 if method=='linear':
  a=(t-x[idx-1])/(x[idx]-x[idx-1]);return np.asarray((1-a)*y[idx-1,:10]+a*y[idx,:10],float)
 if method=='quadratic':
  ii=np.array([max(0,idx-2),idx-1,min(len(x)-1,idx)])
  if len(np.unique(ii))<3:ii=np.array([idx-1,idx,min(len(x)-1,idx+1)])
  if len(np.unique(ii))<3:return interpolate_state(x,y,t,'linear')
  xx=x[ii];yy=y[ii,:10];w=[]
  for k in range(3):
   ids=[i for i in range(3) if i!=k];w.append(np.prod([(t-xx[i])/(xx[k]-xx[i]) for i in ids]))
  return np.asarray(np.array(w)@yy,float)
 raise ValueError(method)

def source_context(parent):
 mod_path=ROOT/'reference/v18_52w_run_audit.py'
 spec=importlib.util.spec_from_file_location('q550_v1852x_authenticated_w',mod_path)
 mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
 xi,states,ks,mid,origin=mod.extract_early(parent)
 temp,eng,_,matched,contexts,sources=mod.stage_engine(parent,xi)
 if abs(mid-C['v1852l_fixed_model_kappa'])>1e-12 or any(abs(v-mid)>1e-12 for v in matched.values()):
  temp.cleanup();raise RuntimeError('AUTHENTIC_FIXED_KAPPA_MISMATCH')
 return temp,eng,contexts,origin,sources

def phase_data(eng,ctx,xi):
 des=eng.original_ramp.Design('v18_50p','beta:2:2',.24)
 s=np.array([eng.clock.cubic_clock_sign(float(ctx['event']+v),des) for v in [xi[0],0.5*(xi[0]+xi[-1]),xi[-1]]],float)
 if not np.all(np.abs(s+1)<1e-10):raise RuntimeError('NON_CONSTANT_LATE_SIGN__NOT_STROBOSCOPIC')
 # The historical ctx theta source, not native column 10 (material tau).
 ph=np.array([ctx['theta'](ctx['event']+float(v)) for v in xi],float)
 if not np.all(np.diff(ph)<0):raise RuntimeError('LATE_PHASE_NOT_MONOTONE')
 return ph,s

def phase_gates(xi,ph,n):
 if ph[0]-ph[-1]<=TWOPI:raise RuntimeError('ONE_TURN_NOT_AVAILABLE')
 gates=[]
 for a in np.linspace(float(ph[0]),float(ph[-1]+TWOPI),int(n)):
  x0=float(np.interp(a,ph[::-1],xi[::-1]))
  x1=float(np.interp(a-TWOPI,ph[::-1],xi[::-1]))
  gates.append((x0,x1,float(a)))
 return gates

def native_names(files):
 allnames=sorted(k[:-7] for k in files if k.endswith('__Y_tau'))
 if len(allnames)!=20 or any(not n.startswith(SEEDS) for n in allnames):raise RuntimeError('EXPECTED_20_NATIVE_SERIES')
 return allnames

def bare(name):return name.replace('__reference_refined','__reference').replace('__kick_refined','__kick')

def short_csv(name,rows):
 if not rows:return
 with (OUT/name).open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def save(name,data): (OUT/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')

def row_one(series,seed,anchor,xa,xb,ph_target,a,b):
 ua=reduce8(a);ub=reduce8(b);du=diff8(ua,ub)
 ja,jb=J(a),J(b)
 # No mixed-unit claim: separately report physical channel components and chart norm.
 return {'seed':seed,'series':series,'anchor_index':anchor,'xi_start':xa,'xi_one_turn':xb,
  'delta_xi':xb-xa,'phase_gate_radians':modangle(ph_target),
  'J_start':ja,'J_after':jb,'J_absolute_drift':jb-ja,
  'relative_8D_chart_norm':float(np.linalg.norm(du)),
  'relative_8D_chart_normalized_descriptive':float(np.linalg.norm(du)/max(1.,np.linalg.norm(ua))),
  'alpha_wrap_difference_rad':float(du[7]),
  **{f'delta_{lab}':float(v) for lab,v in zip(['logr1','v1','omega1','logr2','v2','z','vz','alpha'],du)},
  'full_recurrence_certified':False,
  'note':'ONE_STEP_OBSERVED_SECTION_PAIR__NOT_REPEAT_RETURN_OR_FLOQUET_SPECTRUM'}

def rank_audit(seed,anchor,records):
 """Rank uses four non-refined kicks; 8D chart alone is not charge-conditioned."""
 ref=records[f'{seed}__reference'];u0=reduce8(ref['start']);v0=reduce8(ref['end']);j0=J(ref['start'])
 pair=[]
 for amp in (1e-4,3e-4):
  def pick(sign):
   q=f'{seed}_{"positive" if sign>0 else "negative"}_{amp:.1e}__kick'
   if q not in records:raise RuntimeError('MISSING_KICK '+q)
   return records[q]
  rp=pick(+1);rm=pick(-1)
  plus=(diff8(u0,reduce8(rp['start']))-diff8(u0,reduce8(rm['start'])))/(2*amp)
  plus_out=(diff8(v0,reduce8(rp['end']))-diff8(v0,reduce8(rm['end'])))/(2*amp)
  cj=(J(rp['start'])-J(rm['start']))/(2*amp)
  pair.append({'seed':seed,'anchor_index':anchor,'amplitude':amp,
   'input_8D_directional_norm':float(np.linalg.norm(plus)),
   'output_8D_directional_norm':float(np.linalg.norm(plus_out)),
   'directional_gain_mixed_J':float(np.linalg.norm(plus_out)/max(1e-25,np.linalg.norm(plus))),
   'charge_derivative_per_unit_kick':float(cj),
   'kick_charge_mixing_nonzero':abs(cj)>1e-8,
   'is_fixed_charge_floquet_multiplier':False,
   'note':'CENTRAL_SECANT_OF_MIXED_CHARGE_KICK_FAMILY__NOT_J_RESTRICTED_JACOBIAN'})
 # Input span: four kicks, at MOST rank four in 8D and all with delta J!=0.
 mat=[];charges=[]
 for amp in (1e-4,3e-4):
  for sign in (-1,1):
   q=f'{seed}_{"positive" if sign>0 else "negative"}_{amp:.1e}__kick';rec=records[q]
   mat.append(diff8(u0,reduce8(rec['start']))/(sign*amp))
   charges.append(J(rec['start'])-j0)
 sig=np.linalg.svd(np.array(mat).T,compute_uv=False)
 eff=int(np.sum(sig>sig[0]*1e-2)) if sig.size and sig[0]>0 else 0
 info={'seed':seed,'anchor_index':anchor,'input_displacement_count':4,'chart_dimension_after_fixed_J':8,
  'input_span_max_possible_rank':4,'observed_rank_threshold_1percent':eff,
  'singular_values_scaled_by_nominal_kick':[float(x) for x in sig],
  'J_reference':float(j0),'charge_offsets_for_kicks':[float(x) for x in charges],
  'max_abs_charge_offset':float(max(abs(x) for x in charges)),
  'all_kicks_on_reference_charge':all(abs(x)<1e-10 for x in charges),
  'fixed_charge_8D_Jacobian_identifiable':False,
  'full_floquet_spectrum_identifiable':False,
  'reason':'FOUR_KICK_POINTS_SPAN_AT_MOST_FOUR_DIRECTIONS_AND_EACH_HAS_DIFFERENT_CHARGE__ALSO_ONLY_ONE_RETURN'}
 return pair,info

def plot(sections,dirgain,charge):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return []
 saved=[]
 for seed in SEEDS:
  d=[r for r in sections if r['seed']==seed and r['series']==seed+'__reference']
  g=[r for r in dirgain if r['seed']==seed and r['amplitude']==1e-4]
  if not d:continue
  fig,ax=plt.subplots(3,1,figsize=(9,8.5),sharex=True,constrained_layout=True)
  x=[r['xi_start'] for r in d]
  ax[0].plot(x,[r['relative_8D_chart_norm'] for r in d],'o-',label='one-step reduced state drift')
  ax[1].plot(x,[abs(r['alpha_wrap_difference_rad']) for r in d],'o-',label='sum-angle phase mismatch')
  ax[1].plot(x,[abs(r['delta_omega1']) for r in d],'s-',label='omega1 difference')
  ax[1].plot(x,[abs(r['delta_v2']) for r in d],'^-',label='v2 difference')
  ax[2].plot(x,[r['directional_gain_mixed_J'] for r in g],'o-',label='observed kick-family secant gain (mixed J)')
  for a in ax:a.grid(alpha=.3);a.legend(fontsize=8)
  ax[-1].set_xlabel('Section start xi; one negative external phase winding')
  fig.suptitle('Q550 v18.52x '+seed+' — one step only; no full Floquet spectrum')
  p=OUT/(seed+'_relative_return_map.png');fig.savefig(p,dpi=130);plt.close(fig);saved.append(p.name)
 return saved

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--native-v1852l',type=Path,default=BASE/NATIVE_STEM/'results/extended_horizon_10d_states.npz')
 ap.add_argument('--require-parent',action='store_true');ap.add_argument('--require-native',action='store_true')
 args=ap.parse_args();OUT.mkdir(exist_ok=True,parents=True)
 cand=[args.parent_v1851s,BASE/(PARENT_STEM+'_server_sendback.zip'),Path('/root/elab1/hysteresis')/(PARENT_STEM+'_server_sendback.zip'),Path('/root/Downloads')/(PARENT_STEM+'_server_sendback.zip')]
 parent=next((p for p in cand if p is not None and p.is_file()),None)
 report={'version':'v18.52x','title':C['title'],'parent_available':bool(parent),
  'native_available':args.native_v1852l.is_file(),'new_native_integrations':0,
  'return_map_observations_are_one_step_only':True,'global_limit_cycle_certified':False,
  'full_8D_fixed_J_transverse_floquet_identified':False,'conservative_full_dynamics_certified':False}
 if (args.require_parent and parent is None) or (args.require_native and not args.native_v1852l.is_file()):
  raise RuntimeError('REQUIRED_AUTHENTIC_ANCESTOR_MISSING parent='+str(bool(parent))+' native='+str(args.native_v1852l.is_file()))
 if parent is None or not args.native_v1852l.is_file():
  report['scientific_decision']='PARTIAL_LOCAL_SYNTHETIC_ONLY_MISSING_AUTHENTIC_INPUTS'
  save('report.json',report);return
 if digest(parent)!=C['v1851s_sha256']:raise RuntimeError('PARENT_SHA256_MISMATCH')
 if digest(args.native_v1852l)!=C['v1852l_10d_sha256']:raise RuntimeError('NATIVE_SHA256_MISMATCH')
 tmp=None
 try:
  tmp,eng,ctx,origin,sources=source_context(parent)
  report['original_source_authentication']=origin
  report['verified_parent_sha256']=digest(parent)
  report['verified_native_sha256']=digest(args.native_v1852l)
  report['fixed_native_kappa']=C['v1852l_fixed_model_kappa']
  with np.load(args.native_v1852l,allow_pickle=False) as z:
   xi=np.array(z['common_xi'],float);names=native_names(z.files)
   if not (len(xi)==65537 and abs(xi[0]-12.8)<1e-9 and abs(xi[-1]-25.6)<1e-9):
    raise RuntimeError('UNEXPECTED_NATIVE_LATE_HORIZON')
   ph={};gates={};clock_data={}
   for seed in SEEDS:
    ph[seed],sg=phase_data(eng,ctx[seed],xi)
    gates[seed]=phase_gates(xi,ph[seed],C['phase_section_anchors'])
    clock_data[seed]={'source_sign_samples':[float(v) for v in sg],
      'phase_decrease_radians':float(ph[seed][0]-ph[seed][-1]),
      'phase_turn_count':float((ph[seed][0]-ph[seed][-1])/TWOPI),
      'two_complete_returns_possible':bool((ph[seed][0]-ph[seed][-1])>=2*TWOPI)}
   allrows=[];recstore={seed:[] for seed in SEEDS};jrows=[];interpolation=[]
   for name in names:
    seed=name.split('_',1)[0]
    Y=np.asarray(z[name+'__Y_tau'],float)
    if Y.shape!=(len(xi),11):raise RuntimeError('NATIVE_DIMENSIONS_DRIFT '+name)
    jj=np.exp(2*Y[:,0])*Y[:,3]-np.exp(2*Y[:,4])*Y[:,7]
    delta=float(np.max(np.abs(jj-jj[0])));jrows.append({'series':name,'seed':seed,'J_initial':float(jj[0]),'J_terminal':float(jj[-1]),'max_abs_J_drift':delta})
    if delta>1e-9:raise RuntimeError('CHARGE_NOT_CONSERVED '+name)
    for h,(xa,xb,phase_target) in enumerate(gates[seed]):
     a=interpolate_state(xi,Y,xa);b=interpolate_state(xi,Y,xb)
     qrow=row_one(name,seed,h,xa,xb,phase_target,a,b);allrows.append(qrow)
     recstore[seed].append((name,h,a,b))
     if (h in (0,len(gates[seed])//2,len(gates[seed])-1)) and ('refined' not in name):
      a2=interpolate_state(xi,Y,xa,'quadratic');b2=interpolate_state(xi,Y,xb,'quadratic')
      interpolation.append({'series':name,'anchor_index':h,'start_method_discrepancy_8D':float(np.linalg.norm(diff8(reduce8(a),reduce8(a2)))),
       'end_method_discrepancy_8D':float(np.linalg.norm(diff8(reduce8(b),reduce8(b2))))})
  short_csv('one_step_fixed_phase_return_pairs.csv',allrows)
  short_csv('per_trajectory_charge_conservation.csv',jrows)
  short_csv('section_interpolation_precision.csv',interpolation)
  direction=[];ranks=[]
  for seed in SEEDS:
   for anchor in range(len(gates[seed])):
    records={name:{'start':a,'end':b} for name,h,a,b in recstore[seed] if h==anchor}
    pairs,info=rank_audit(seed,anchor,records)
    direction.extend(pairs);ranks.append(info)
  short_csv('mixed_charge_kick_family_directional_response.csv',direction)
  save('section_input_rank_and_identifiability.json',ranks)
  # Accuracy of refined traces: compare same-channel reduced returns, do not treat refined as new independent kicks.
  rindex={(r['series'],r['anchor_index']):r for r in allrows}
  refined=[]
  for name in names:
   if not ('refined' in name):continue
   base=bare(name)
   for anchor in range(len(gates[name.split('_',1)[0]])):
    a=rindex[(name,anchor)];b=rindex[(base,anchor)]
    refined.append({'series_refined':name,'series_base':base,'anchor_index':anchor,
     'difference_return_norm':float(a['relative_8D_chart_norm']-b['relative_8D_chart_norm']),
     'difference_charge_drift':float(a['J_absolute_drift']-b['J_absolute_drift'])})
  short_csv('refined_vs_base_return_precision.csv',refined)
  report['late_clock']=clock_data;report['native_series']=len(names);report['section_anchors_per_seed']=C['phase_section_anchors'];report['one_step_pairs']=len(allrows)
  report['one_step_reference_min_norm']={seed:min(r['relative_8D_chart_norm'] for r in allrows if r['series']==seed+'__reference') for seed in SEEDS}
  report['one_step_reference_max_norm']={seed:max(r['relative_8D_chart_norm'] for r in allrows if r['series']==seed+'__reference') for seed in SEEDS}
  report['max_abs_J_drift']=max(r['max_abs_J_drift'] for r in jrows)
  report['max_reference_charge_mismatch_from_kicks']=max(r['max_abs_charge_offset'] for r in ranks)
  report['max_input_rank_1percent']=max(r['observed_rank_threshold_1percent'] for r in ranks)
  report['min_input_rank_1percent']=min(r['observed_rank_threshold_1percent'] for r in ranks)
  report['max_quadratic_vs_linear_interpolation_residual']=max(max(r['start_method_discrepancy_8D'],r['end_method_discrepancy_8D']) for r in interpolation)
  report['observed_at_least_two_full_returns']=all(v['two_complete_returns_possible'] for v in clock_data.values())
  report['one_step_relative_fixed_J_pairs_meaning']='Individual trajectory carries constant J; different kicks have different J; pairs are not samples of a single fixed-J manifold map.'
  report['eight_dimensional_chart']='(logr1,v1,omega1,logr2,v2,z,vz,theta1+theta2 mod 2pi); omega2 reconstructed from J'
  report['floquet_multipliers']=None;report['relative_limit_cycle_certified']=False
  report['scientific_decision']='AUTHENTIC_ONE_STEP_PHASE_SECTION_MAP_WITH_CHARGE_STRATIFICATION__FIXED_J_8D_FLOQUET_UNIDENTIFIABLE'
  report['plots']=plot(allrows,direction,jrows)
  save('report.json',report)
  lines=['# Q550 v18.52x — authenticated server verdict','',report['scientific_decision'],'',
   '## Scope','- External phase section mod 2pi on the constant signed-drive late return branch.',
   '- Relative symmetry quotient reduces 10D physical state to a nine-dimensional chart; fixing J gives eight coordinates.',
   '- The 20 trajectories give one genuine return per selected phase gate; they do NOT give two iterates at the same gate.','',
   '## Key numerical facts',f'- {len(allrows)} native single-step return pairs, {len(names)} states, {len(ranks)} gate identifiability checks.',
   f'- Maximum native conservation drift in J: {report["max_abs_J_drift"]:.9g}.',
   f'- Maximum charge mismatch of kicked orbits relative to seed reference: {report["max_reference_charge_mismatch_from_kicks"]:.9g}.',
   f'- Observed input-rank 1% range: {report["min_input_rank_1percent"]}..{report["max_input_rank_1percent"]}, rank 8 required.',
   f'- Full two-returns observed: {report["observed_at_least_two_full_returns"]}.',
   '','## Scientific classification',
   'No fixed-J 8D Jacobian, Floquet spectrum, limit cycle, or attraction can be certified from only four mixed-charge kicked orbits per seed.',
   'Secant gain is defined only for the historical mixed-J kick family, not as a Floquet multiplier.',
   'A future independent fixed-J tangent-basis run and at least two genuine section crossings would be needed for transverse stability.']
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
 finally:
  if tmp is not None:tmp.cleanup()

if __name__=='__main__':
 try:main()
 except Exception as exc:
  OUT.mkdir(exist_ok=True,parents=True)
  (OUT/'FAILURE_DIAGNOSTIC.txt').write_text(traceback.format_exc())
  print('Q550 v18.52x FAILURE',repr(exc),flush=True)
  raise
