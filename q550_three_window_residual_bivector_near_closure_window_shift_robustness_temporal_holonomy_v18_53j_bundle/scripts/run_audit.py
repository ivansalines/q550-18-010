#!/usr/bin/env python3
"""Q550 v18.53j: predeclared three-window robustness and polygon-plane transport.

Historical states and forcing are recovered ONLY from SHA authenticated v18.51s
and v18.52l ancestors.  Material-time skew torques are integrated along true
parameter-distinct IVPs.  Beyond xi=25.6 the forcing clock is extrapolated.
No claim of a physical orbit return, Floquet theorem, or intrinsic holonomy.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REF=ROOT/'reference';OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'scripts'))
from bivector_geometry import principal_skew,cancel_metric
from three_window_geometry import three_window,accumulators,interval,independent_direct_window_integral
CFG=json.loads((ROOT/'inputs_config.json').read_text());T=CFG['three_window'];K=CFG['kappas']

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def put(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False,ensure_ascii=False)+'\n')

def csv_out(name,rows):
 if not rows:return
 fields=list(dict.fromkeys(x for r in rows for x in r))
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def hist_module():
 p=REF/'v18_53h_fix1_run_audit.py'
 spec=importlib.util.spec_from_file_location('q550_53j_authentic_h_module',p)
 module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
 return module

def validate_genealogy(h):
 h.provenance() # Includes all 53g embedded source and output hash checks.
 for file,digest in CFG['v1853h_references_sha256'].items():
  if sha(REF/file)!=digest:raise RuntimeError('H_BASELINE_INPUT_DRIFT '+file)
 for file,digest in CFG['v1853i_references_sha256'].items():
  if sha(REF/file)!=digest:raise RuntimeError('I_BASELINE_INPUT_DRIFT '+file)
 ledger=json.loads((REF/'v18_53i_OUTPUT_MANIFEST.json').read_text())['files']
 needed={'v18_53i_authenticated_report.json':'results/report.json',
     'v18_53i_residual_bivector_windows.csv':'results/residual_bivector_windows.csv',
     'v18_53i_integrated_kappa_runs.csv':'results/integrated_kappa_runs.csv',
     'v18_53i_principal_plane_time_transport.csv':'results/principal_plane_time_transport.csv',
     'v18_53i_kappa_principal_plane_separation.csv':'results/kappa_principal_plane_separation.csv',
     'v18_53i_original_run_audit.py':'scripts/run_audit.py',
     'v18_53i_original_test_synthetic.py':'tests/test_v1853i_synthetic.py'}
 for local,old in needed.items():
  if ledger.get(old)!=sha(REF/local):raise RuntimeError('I_OUTPUT_MANIFEST_SHA256_FAIL '+local)
 old=json.loads((REF/'v18_53i_authenticated_report.json').read_text())
 if old.get('scientific_decision')!='AUTHENTIC_DUAL_KAPPA_PARAMETER_SENSITIVITY_AND_EXTENDED_MODEL_CLOCK_IVPS__SAMPLED_GEOMETRY_PERSISTENT' or old.get('actual_completed_Q550_ivps')!=26:
  raise RuntimeError('I_BASELINE_UNAUTHENTIC_OR_INCOMPLETE')
 if old.get('parent_status')!='SHA256_AUTHENTICATED' or old.get('native_status')!='SHA256_AUTHENTICATED':
  raise RuntimeError('I_SOURCE_PARENT_MISSING')
 if dict((k,v) for k,v in old['kappa_parameter_grid'] if k in ('original','autonomous','shared'))!=K:
  raise RuntimeError('I_HISTORIC_KAPPA_MISMATCH')
 with (REF/'v18_53i_residual_bivector_windows.csv').open() as f:
  old_windows={(r['seed'],r['kappa_rule'],r['window']):r for r in csv.DictReader(f)}
 if len(old_windows)!=132:raise RuntimeError('I_EXPECTED_132_WINDOW_REFERENCE_ROWS')
 return old,old_windows

def xi_grid(full_native_xi):
 orig=np.asarray(full_native_xi,float)[::CFG['native_stride']]
 if orig.size<100 or abs(orig[0]-T['start_xi'])>1e-9 or abs(orig[-1]-T['native_end_xi'])>1e-9:
  raise RuntimeError('NATIVE_GRID_NOT_EXPECTED')
 n=int(round((T['continuation_xi_max']-T['native_end_xi'])/T['grid_step_extension_xi']))
 extra=np.linspace(T['native_end_xi'],T['continuation_xi_max'],n+1)[1:]
 return np.r_[orig,extra],len(orig)

def window_triples():
 for o in T['origins_xi']:
  for w in T['widths_xi']:
   if o+3*w>T['continuation_xi_max']+1e-8:continue
   yield ('origin_width',f'origin={o:g}|width={w:g}',[o,o+w,o+2*w,o+3*w],dict(origin_shift_xi=o-T['start_xi'],width_xi=w))
 st=T['start_xi'];width=T['historical_width_xi']
 for delta1 in T['boundary_jitter_xi']:
  for delta2 in T['boundary_jitter_xi']:
   yield ('internal_boundary_jitter',f'boundary1_delta={delta1:g}|boundary2_delta={delta2:g}',
      [st,st+width+delta1,st+2*width+delta2,T['fixed_boundary_end_xi']],
      dict(internal_boundary_delta1_xi=delta1,internal_boundary_delta2_xi=delta2))

def simple_row(sc,group,name,seed,kn,kappa,extra):
 ho1=sc['plane1_geodesic_polygon_holonomy'];ho2=sc['plane2_geodesic_polygon_holonomy']
 return dict(seed=seed,kappa_rule=kn,kappa=float(kappa),scan_family=group,scan_id=name,
     xi_start=sc['xi_start'],xi_end=sc['xi_end'],boundary1_xi=sc['left_boundary_xi'],boundary2_xi=sc['right_boundary_xi'],
     segment1_xi=sc['window_lengths_xi'][0],segment2_xi=sc['window_lengths_xi'][1],segment3_xi=sc['window_lengths_xi'][2],
     R1_norm=sc['norms'][0],R2_norm=sc['norms'][1],R3_norm=sc['norms'][2],
     R_total_norm=sc['net_norm'],cancellation_three_window=sc['normalized_residual_closure_score'],
     residual_magnitude_cv=sc['equal_window_residual_norm_cv'],
     angle12_deg=sc['pairwise_bivector_angles_deg'][0],angle23_deg=sc['pairwise_bivector_angles_deg'][1],
     angle13_deg=sc['pairwise_bivector_angles_deg'][2],
     max_120_deviation_deg=sc['max_deviation_from_120_deg'],
     integrated_radial_phase_cancellation=sc['whole_interval_radial_phase_cancellation'],
     gram_min_eigen_fraction=sc['gram_smallest_eigenvalue_fraction'],
     holo1_polygon_deg=ho1['angle_deg'],holo2_polygon_deg=ho2['angle_deg'],
     holo1_defined=ho1['defined'],holo2_defined=ho2['defined'],
     holo1_reason=ho1['reason'],holo2_reason=ho2['reason'],
     holo1_min_overlap=ho1['min_adjacent_overlap_cosine'],holo2_min_overlap=ho2['min_adjacent_overlap_cosine'],
     additive_error=sc['additive_bivector_error'],
     strong_near_closure=sc['normalized_residual_closure_score']>=T['strong_near_closure_threshold'],
     near_120_and_balanced=(sc['max_deviation_from_120_deg']<=T['near_equiangular_tolerance_deg'] and
          sc['equal_window_residual_norm_cv']<=T['equal_magnitude_cv_threshold']),
     **extra)

def tensor_coords(R):
 return {f'R{i}{j}':float(R[i,j]) for i in range(5) for j in range(i+1,5)}

def check_old_windows(acc,seed,kn,old_windows):
 rows=[]
 for ix in range(1,4):
  a=T['start_xi']+(ix-1)*T['historical_width_xi'];b=a+T['historical_width_xi']
  U,C=interval(acc,a,b);R=U+C;pre=old_windows[(seed,kn,f'segment_{ix}')]
  dn=abs(float(np.linalg.norm(R))-float(pre['residual_norm']))
  dc=abs(float(cancel_metric(U,C)['invariant_cancellation'])-float(pre['cancellation']))
  if dn>T['v1853i_residual_norm_abs_tolerance'] or dc>T['v1853i_cancellation_abs_tolerance']:
   raise RuntimeError(f'V1853I_SEPARATE_SEGMENT_REPLAY_FAIL {seed} {kn} {ix} norm={dn} cancellation={dc}')
  rows.append(dict(seed=seed,kappa_rule=kn,segment=ix,old_norm=float(pre['residual_norm']),
      new_norm=float(np.linalg.norm(R)),norm_abs_diff=dn,old_cancellation=float(pre['cancellation']),
      new_cancellation=float(cancel_metric(U,C)['invariant_cancellation']),cancellation_abs_diff=dc))
 return rows

def plot(rows):
 try:
  import matplotlib;matplotlib.use('Agg')
  import matplotlib.pyplot as plt
 except Exception:return
 for seed in ('original','autonomous'):
  data=[r for r in rows if r['seed']==seed and r['kappa_rule']=='shared' and r['scan_family']=='origin_width']
  fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
  for width in T['widths_xi']:
   part=sorted([r for r in data if abs(r['width_xi']-width)<1e-6],key=lambda v:v['origin_shift_xi'])
   ax.plot([r['origin_shift_xi'] for r in part],[r['cancellation_three_window'] for r in part],label=f'width {width:g}',marker='o')
  ax.axhline(T['strong_near_closure_threshold'],ls='--',color='grey',label='predeclared 90%')
  ax.set_xlabel('shift of common starting xi');ax.set_ylabel('Three-window near-closure score')
  ax.grid(alpha=.2);ax.legend(fontsize=8,ncol=2)
  ax.set_title(f'v18.53j — whole-interval shift and width ({seed} seed, shared κ₀)')
  fig.savefig(OUT/f'three_window_shift_width_{seed}.png',dpi=155);plt.close(fig)
  jitter=[r for r in rows if r['seed']==seed and r['kappa_rule']=='shared' and r['scan_family']=='internal_boundary_jitter']
  fig,ax=plt.subplots(figsize=(8,5.8),layout='constrained')
  mat=np.full((len(T['boundary_jitter_xi']),len(T['boundary_jitter_xi'])),np.nan)
  for r in jitter:
   i=T['boundary_jitter_xi'].index(r['internal_boundary_delta1_xi']);j=T['boundary_jitter_xi'].index(r['internal_boundary_delta2_xi'])
   mat[i,j]=r['max_120_deviation_deg']
  im=ax.imshow(mat,origin='lower',aspect='auto')
  ax.set_xticks(range(5),[f'{v:g}' for v in T['boundary_jitter_xi']]);ax.set_yticks(range(5),[f'{v:g}' for v in T['boundary_jitter_xi']])
  ax.set_xlabel('shift of 2nd internal boundary (xi)');ax.set_ylabel('shift of 1st internal boundary (xi)')
  ax.set_title(f'v18.53j — angular distortion at FIXED full endpoints ({seed} seed)')
  fig.colorbar(im,ax=ax,label='max deviation from 120 deg')
  fig.savefig(OUT/f'fixed_endpoint_boundary_jitter_{seed}.png',dpi=155);plt.close(fig)

def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--base',type=Path,default=ROOT.parent)
 parser.add_argument('--native-npz',type=Path);parser.add_argument('--parent-v1851s',type=Path)
 parser.add_argument('--require-parent',action='store_true');parser.add_argument('--require-native',action='store_true')
 parser.add_argument('--require-dynamic',action='store_true')
 args=parser.parse_args()
 status=dict(version='v18.53j',scientific_decision='UNCLASSIFIED',parent_status='MISSING',
     native_status='MISSING',actual_completed_Q550_ivps=0,
     actual_completed_refinement_ivps=0,original_historical_trace_ends_xi=25.6,
     continued_external_clock_to_xi=T['continuation_xi_max'],
     not_a_physical_orbit_holonomy_or_full_ramp_certificate=True,
     fixed_endpoints_algebraic_additivity_not_robustness_evidence=True)
 temp=None;exit_code=0
 try:
  h=hist_module();old,ref_windows=validate_genealogy(h)
  status['v1853i_reference_status']='SHA256_AND_OUTPUT_MANIFEST_AUTHENTICATED'
  L=h.loader();parent=L.find_parent(args.parent_v1851s);native=args.native_npz or L.find_npz('10d',args.base)
  if parent and Path(parent).is_file():
   if sha(parent)!=CFG['parent_v1851s_sha256']:raise RuntimeError('V1851S_PARENT_SHA_DRIFT')
   status['parent_status']='SHA256_AUTHENTICATED'
  if native and Path(native).is_file():
   if sha(native)!=CFG['native_10d_sha256']:raise RuntimeError('V1852L_NATIVE_SHA_DRIFT')
   status['native_status']='SHA256_AUTHENTICATED'
  if status['parent_status']=='SHA256_AUTHENTICATED' and status['native_status']=='SHA256_AUTHENTICATED':
   temp,eng,contexts,meta=L.stage_original(parent);status['engine_provenance']=meta
   scans=[];components=[];reconcile=[];runs=[];refinements=[];baseline=[];quadrature=[]
   with np.load(native,allow_pickle=False) as f:
    xi,hcount=xi_grid(f['common_xi'])
    for seed in ('original','autonomous'):
     ctx=contexts[seed]
     b,clock_error,archived=h.native_53g_reconciliation(f,xi[:hcount],seed,ctx,eng.omega)
     if clock_error>CFG['long_horizon']['material_clock_rate_rel_tolerance']:
      raise RuntimeError('HISTORICAL_MATERIAL_CLOCK_FAILURE '+seed)
     reconcile.append(dict(seed=seed,v1853g_reference_error=b,material_clock_rate_relative_error=clock_error))
     z0=np.asarray(archived[0],float)
     design=eng.original_ramp.Design('v18_50p','beta:2:2',.24)
     clock_probe=np.linspace(T['start_xi'],T['continuation_xi_max'],33)
     theta=np.asarray([ctx['theta'](ctx['event']+float(u)) for u in clock_probe])
     clock=np.asarray([eng.clock.cubic_clock_sign(ctx['event']+float(u),design) for u in clock_probe])
     if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(clock)) or np.max(np.abs(clock))>1.00001:
      raise RuntimeError('AUTHENTIC_CLOCK_CONTINUATION_NOT_DEFINED')
     status.setdefault('clock',{})[seed]=dict(end_theta=float(theta[-1]),end_drive_sign=float(clock[-1]),
        historical_clock_extrapolated_beyond_xi25_6=True)
     for kn in ('shared','original','autonomous'):
      k=K[kn]
      z,info=h.solve(eng,ctx,k,xi,z0)
      status['actual_completed_Q550_ivps']+=1
      if not np.all(np.isfinite(z)) or np.max(np.abs(z[:,[0,4]]))>40 or np.max(np.abs(z[:,8]))>700:
       raise RuntimeError('BAD_INTEGRATED_Q550_STATE '+seed+' '+kn)
      src_error=h.validated_rhs(eng,ctx,k,xi,z)
      if src_error>T['source_field_rel_tolerance']:raise RuntimeError('AUTHENTIC_SOURCE_FIELD_ERROR')
      d=h.primary_metrics(xi,z,ctx,k)
      jd=d['metric']['charge_range']
      if jd>T['charge_abs_drift_tolerance']:raise RuntimeError('CHARGE_CONSERVATION_ERROR')
      original_replay=None
      if kn=='shared':
       original_replay=float(np.max(np.abs(z[:hcount]-archived)/np.maximum(1,np.abs(archived))))
       if original_replay>T['historical_replay_max_scaled_error']:
        raise RuntimeError('NATIVE_52L_SHARED_REPLAY_FAIL')
      if kn in ('original','autonomous'):
       z2,ii=h.solve(eng,ctx,k,xi,z0,tight=True)
       status['actual_completed_Q550_ivps']+=1
       status['actual_completed_refinement_ivps']+=1
       diff=float(np.max(np.abs(z-z2)/np.maximum(1,np.abs(z2))))
       refinements.append(dict(seed=seed,kappa_rule=kn,max_scaled_error=diff,nfev=ii['nfev']))
       if diff>T['integrator_refinement_scaled_max']:raise RuntimeError('LONG_IVP_REFINEMENT_FAILURE')
      acc=accumulators(xi,d['t'],d['B_U'],d['B_C'])
      checked=check_old_windows(acc,seed,kn,ref_windows);baseline.extend(checked)
      runs.append(dict(seed=seed,kappa_rule=kn,kappa=k,nfev=info['nfev'],true_integrated_IVP=True,
        max_scaled_native_replay=original_replay,charge_range=jd,source_field_error=src_error,
        material_clock_not_a_historical_extension=True))
      for fam,label,bounds,extras in window_triples():
       sc,segments=three_window(acc,bounds,T['holonomy_min_adjacent_cosine'],T['principal_gap_warning_relative'])
       rec=simple_row(sc,fam,label,seed,kn,k,extras)
       # The three residuals sum exactly at fixed full endpoints; it is the
       # PAIRWISE geometry and normalization that may depend on internal cuts.
       scans.append(rec)
       for j,v in enumerate(segments,1):
        components.append(dict(seed=seed,kappa_rule=kn,kappa=k,scan_family=fam,scan_id=label,
          segment=j,xi_start=bounds[j-1],xi_end=bounds[j],
          radial_fro=float(np.linalg.norm(v['U'])),phase_fro=float(np.linalg.norm(v['C'])),
          **tensor_coords(v['R'])))
      # Independently recompute a small predeclared set of OFF-grid quadrature
      # integrals, not merely verify the telescoping cumulative formula.
      for a,bx in [(13.17,26.31),(24.89,41.63),(40.21,55.77)]:
       u,c=interval(acc,a,bx)
       ui=independent_direct_window_integral(xi,d['t'],d['B_U'],a,bx)
       ci=independent_direct_window_integral(xi,d['t'],d['B_C'],a,bx)
       mismatch=float(np.linalg.norm((u+c)-(ui+ci))/max(1,np.linalg.norm(u+c)))
       quadrature.append(dict(seed=seed,kappa_rule=kn,xi_start=a,xi_end=bx,
           independent_interpolated_quadrature_relative_error=mismatch))
       if mismatch>T['max_quadrature_halfgrid_relative']:
        raise RuntimeError('INDEPENDENT_OFFGRID_QUADRATURE_FAILURE '+str(mismatch))
      print('V18_53J_IVP_DONE',seed,kn,'completed',status['actual_completed_Q550_ivps'],flush=True)
   if status['actual_completed_Q550_ivps']!=T['expected_primary_ivps']+T['expected_refined_ivps']:
    raise RuntimeError('EXPECTED_IVP_COUNT_NOT_MATCHED')
   if len(scans)!=2*3*(len(T['origins_xi'])*len(T['widths_xi'])+len(T['boundary_jitter_xi'])**2):
    raise RuntimeError('EXPECTED_PREDECLARED_TRIPLE_COUNT_NOT_MATCHED')
   # Explicitly prove fixed-endpoint net invariance and refuse to count it as a
   # physical robustness check. The denominator and angles are free to change.
   for seed in ('original','autonomous'):
    for kn in ('shared','original','autonomous'):
     jit=[r for r in scans if r['seed']==seed and r['kappa_rule']==kn and r['scan_family']=='internal_boundary_jitter']
     changes=np.ptp([r['R_total_norm'] for r in jit])
     if changes>1e-11:raise RuntimeError('FIXED_ENDPOINT_RESIDUAL_INVARIANCE_FAIL')
   baseline_cases=[r for r in scans if r['scan_family']=='origin_width' and abs(r['xi_start']-12.8)<1e-8 and abs(r['segment1_xi']-12.8)<1e-8]
   if len(baseline_cases)!=6:raise RuntimeError('BASELINE_THREE_WINDOW_53I_RECONCILIATION_MISSING')
   shift=[r for r in scans if r['scan_family']=='origin_width']
   jitter=[r for r in scans if r['scan_family']=='internal_boundary_jitter']
   min_shift=min(r['cancellation_three_window'] for r in shift)
   max_shift=max(r['cancellation_three_window'] for r in shift)
   max_dev=max(r['max_120_deviation_deg'] for r in shift)
   num_holo=sum(int(r['holo1_defined'])+int(r['holo2_defined']) for r in scans)
   status['scientific_decision']=('AUTHENTIC_THREE_WINDOW_SHIFT_WIDTH_SCAN__WINDOW_SENSITIVE_OR_NONCLOSING'
       if min_shift<T['strong_near_closure_threshold'] else
       'AUTHENTIC_THREE_WINDOW_SHIFT_WIDTH_SCAN__NEAR_CLOSURE_SAMPLED_PERSISTENT')
   csv_out('window_shift_width_triplets.csv',shift)
   csv_out('fixed_endpoint_internal_boundary_jitter.csv',jitter)
   csv_out('all_residual_bivector_components.csv',components)
   csv_out('v1853i_three_window_reconciliation.csv',baseline)
   csv_out('integrated_ivp_replay.csv',runs)
   csv_out('offgrid_quadrature_reconciliation.csv',quadrature)
   put('integrator_refinements.json',refinements)
   put('v1853g_clock_reconciliations.json',reconcile)
   put('baseline_three_window_geometry.json',baseline_cases)
   status['summary']=dict(n_completed_ivps=status['actual_completed_Q550_ivps'],
      n_unique_dynamical_traces=len(runs),n_window_triplets=len(scans),
      n_whole_interval_origin_width_triplets=len(shift),n_fixed_endpoints_boundary_jitter_triplets=len(jitter),
      n_principal_plane_geodesic_polygon_diagnostics_defined=num_holo,
      n_strong_near_closure_origin_width=sum(r['strong_near_closure'] for r in shift),
      n_near_equiangular_and_balanced_origin_width=sum(r['near_120_and_balanced'] for r in shift),
      min_shift_width_closure_score=min_shift,max_shift_width_closure_score=max_shift,
      max_shift_width_120_degree_deviation=max_dev,
      baseline_shared_original_seed_closure_score=next(r['cancellation_three_window'] for r in baseline_cases if r['seed']=='original' and r['kappa_rule']=='shared'),
      maximum_native_replay_error=max(r['max_scaled_native_replay'] or 0 for r in runs),
      maximum_authentic_source_error=max(r['source_field_error'] for r in runs),
      maximum_charge_drift=max(r['charge_range'] for r in runs),
      maximum_refinement_error=max(r['max_scaled_error'] for r in refinements),
      maximum_53i_residual_reconciliation_error=max(max(r['norm_abs_diff'],r['cancellation_abs_diff']) for r in baseline),
      maximum_independent_offgrid_integral_error=max(r['independent_interpolated_quadrature_relative_error'] for r in quadrature),
      first_12_8_to_25_6_history_replayed=True,
      beyond_xi_25_6_is_model_clock_continuation=True)
   plot(scans)
  else:
   status['scientific_decision']='PARTIAL_AUTHENTIC_V1853I_REFERENCE__NEW_IVPS_NOT_RUN'
 except Exception as exc:
  status['scientific_decision']='AUDIT_FAIL_NOT_CERTIFIED';status['failure']=str(exc)
  status['traceback']=traceback.format_exc();exit_code=2
 finally:
  if temp is not None:temp.cleanup()
  put('report.json',status)
  (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join([
    '# Q550 v18.53j — Three-Window Near-Closure, Shift Robustness & Temporal Holonomy Audit','',
    '## Scientific decision',f"`{status['scientific_decision']}`",'',
    '## Completed numerical integrations',str(status['actual_completed_Q550_ivps']),'',
    '## Source authentication',f"v18.51s: {status['parent_status']}; v18.52l: {status['native_status']}; v18.53i ancestor: {status.get('v1853i_reference_status','NOT_CHECKED')}",'',
    '## Critical interpretation',
    'Internal-boundary jitter at FIXED external endpoints cannot change the exact summed bivector. This is algebraic additivity, NOT persistence evidence. Angular and relative-strength changes under this jitter are informative.',
    'Whole-origin translation and three-window width changes MODIFY the total interval and can genuinely falsify an apparently fine-tuned near-closure.',
    'The optional geodesic polygon of THREE principal 2-planes is closed by construction in the Grassmannian. Its angle is a gauge-invariant OBSERVER-SPACE POLYGON diagnostic, NOT an independently established temporal holonomy of a physically closed orbit.',
    'The two seed families, historical original/autonomous kappa values and shared kappa are all independently integrated from the same initial state within each seed family.',
    'The source history ends at xi=25.6. For xi>25.6 the authenticated source clock is continued formally; no unseen historical orbit or complete ramp is claimed.',
    'No global periodic orbit, physical hinge Noether charge, or three-step dynamical law is certified by integrated residual near-closure.','',
    '## Numeric summary',json.dumps(status.get('summary',{'not_computed':status.get('failure','missing inputs')}),indent=2)
  ])+'\n')
  print('V18_53J_DECISION',status['scientific_decision'],'ACTUAL_Q550_IVPS',status['actual_completed_Q550_ivps'],flush=True)
  if args.require_parent and status['parent_status']!='SHA256_AUTHENTICATED':exit_code=max(exit_code,3)
  if args.require_native and status['native_status']!='SHA256_AUTHENTICATED':exit_code=max(exit_code,3)
  if args.require_dynamic and not status['scientific_decision'].startswith('AUTHENTIC_THREE_WINDOW_SHIFT_WIDTH_SCAN'):exit_code=max(exit_code,3)
 return exit_code
if __name__=='__main__':sys.exit(main())
