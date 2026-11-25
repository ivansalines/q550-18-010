#!/usr/bin/env python3
"""Q550 v18.53p: retrospective convergence audit from authenticated v18.53o data.
No new physical IVPs are run; all output is derived numerical geometry.
"""
import argparse,csv,hashlib,io,json,sys,zipfile
from pathlib import Path
import numpy as np
from valley_convergence_math import Primitive,optimum_width,eval_curve,arcs
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs/v18_53o_authenticated_server_sendback.zip'

def sha_bytes(x):return hashlib.sha256(x).hexdigest()
def sha_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as fd:
        for x in iter(lambda:fd.read(1<<20),b''):h.update(x)
    return h.hexdigest()
def authenticated_o():
    if sha_file(PARENT)!=CFG['parent_v18_53o_sha256']:
        raise RuntimeError('PARENT_O_SHA256_MISMATCH')
    with zipfile.ZipFile(PARENT) as z:
        if z.testzip() is not None:raise RuntimeError('PARENT_O_ZIP_CRC_FAILURE')
        names=z.namelist()
        if any(Path(x).is_absolute() or '..' in Path(x).parts for x in names):raise RuntimeError('UNSAFE_ARCHIVE_MEMBER')
        roots={x.split('/')[0] for x in names if '/' in x}
        if len(roots)!=1:raise RuntimeError('PARENT_O_MULTIPLE_ROOTS')
        prefix=next(iter(roots))+'/'
        counts={}
        for manifest in ['SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json']:
            obj=json.loads(z.read(prefix+manifest))
            for relative,h in obj['files'].items():
                if sha_bytes(z.read(prefix+relative))!=h:
                    raise RuntimeError('PARENT_HASH_MISMATCH '+relative)
            counts[manifest]=len(obj['files'])
        report=json.loads(z.read(prefix+'results/report.json'))
        if (report.get('scientific_decision')!=CFG['required_parent_decision'] or
                report.get('actual_Q550_ivps')!=10 or report.get('refined_Q550_ivps')!=4):
            raise RuntimeError('PARENT_O_NOT_AUTHENTIC_DYNAMICAL_RUN')
        source=np.load(io.BytesIO(z.read(prefix+'results/derived_primitives.npz')),allow_pickle=False)
        def rows(file):return list(csv.DictReader(io.StringIO(z.read(prefix+'results/'+file).decode())))
        fixed=rows('fixed_width_oriented_transport.csv')
        valley=rows('local_valley_oriented_transport.csv')
        old_summary=rows('valley_curve_summary.csv')
        if len(fixed)!=846 or len(valley)!=846 or len(old_summary)!=12:
            raise RuntimeError('PARENT_O_INCOMPLETE_CURVES')
        data={key:np.asarray(source[key]).copy() for key in source.files}
        source.close()
    return data,fixed,valley,old_summary,counts,report

def csv_write(filename,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r.keys()))
    with (OUT/filename).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)

def clean_old(rows,seed,rule):
    rr=[r for r in rows if r['seed']==seed and r['kappa_rule']==rule]
    return sorted(rr,key=lambda r:int(r['index']))

def run(deep=False):
    data,old_fixed,old_valley,old_summary,pins,parent_report=authenticated_o()
    xi=data['common_xi']
    points=CFG['grid_points']+(CFG['deep_grid_points'] if deep else [])
    summary=[]; samples=[]; parent_checks=[]; vectors={}; benchmark=[]
    for seed in ['original','autonomous']:
      for rule in ['shared','original','autonomous']:
        label=seed+'__'+rule
        c=data[label+'__primitive']
        r_fixed=clean_old(old_fixed,seed,rule)
        r_valley=clean_old(old_valley,seed,rule)
        if len(r_fixed)!=141 or len(r_valley)!=141:raise RuntimeError('MISSING_HISTORICAL_CASE '+label)
        mid=r_fixed[70];o0=float(mid['origin_xi']);w0=float(mid['width_xi'])
        if abs(float(r_valley[70]['origin_xi'])-o0)>1e-12:raise RuntimeError('PARENT_CENTER_MISMATCH')
        low=w0-CFG['width_halfband'];high=w0+CFG['width_halfband']
        for mode in CFG['interpolants']:
          p=Primitive(xi,c,mode)
          for npoints in points:
            xx=np.linspace(o0-CFG['origin_halfspan'],o0+CFG['origin_halfspan'],npoints)
            if xx[0]<xi[0] or xx[-1]+3*high>xi[-1]:raise RuntimeError('OUTSIDE_TRACE')
            ww=np.array([optimum_width(p,float(o),low,high)[0] for o in xx])
            fixed_width=np.full_like(xx,w0)
            for curve,w in [('fixed_width',fixed_width),('optimized_valley',ww)]:
              stats,series=eval_curve(p,xx,w)
              if min(series['r'])<CFG['min_norm_threshold']:
                  raise RuntimeError('UNEXPECTED_NEAR_ZERO_DIRECTION_UNDEFINED')
              is_bound=int(np.sum((abs(w-low)<1e-7)|(abs(w-high)<1e-7))) if curve=='optimized_valley' else 0
              row=dict(seed=seed,kappa_rule=rule,interpolant=mode,curve=curve,
                n_points=npoints,origin_center=o0,original_width=w0,
                n_width_bounds=is_bound,**stats)
              summary.append(row)
              for i in range(npoints):
                samples.append(dict(seed=seed,kappa_rule=rule,interpolant=mode,curve=curve,
                 n_points=npoints,i=i,origin_xi=float(xx[i]),width_xi=float(w[i]),end_xi=float(series['end'][i]),
                 norm_10d=float(series['r'][i]),width_slope=float(series['wdot'][i]),
                 width_slope_pchip=float(series['wdot_pchip'][i]),
                 directional_speed=float(series['speed'][i]),
                 nearest_endpoint_knot_distance=float(p.clearance(series['end'][i])),
                 arc_to_next=float(series['arc_steps'][i]) if i<npoints-1 else 0.))
              # Retain unit vectors in binary, allowing exact independent reconstruction.
              key=label+'__'+mode+'__'+str(npoints)+'__'+curve
              vectors[key+'__unit10d']=series['u']
              vectors[key+'__origin']=xx
              vectors[key+'__width']=w
              vectors[key+'__oriented_polygon45d']=series['polygon_area']
              if mode=='linear' and npoints==141:
                historic=r_fixed if curve=='fixed_width' else r_valley
                old_r=np.array([float(r['norm_10d']) for r in historic])
                old_w=np.array([float(r['width_xi']) for r in historic])
                old_u=np.stack([np.fromstring(r['unit_10d_components'],sep=' ') for r in historic])
                old_speed=np.array([float(r['directional_angular_speed_per_xi']) for r in historic])
                old_arc=float(sum(arcs(old_u)))
                old_chain=float(np.trapezoid(old_speed,xx))
                hrow=dict(seed=seed,kappa_rule=rule,curve=curve,
                  historical_direct_geodesic_arc=old_arc,historical_chain_integral=old_chain,
                  historical_chord_minus_chain=old_arc-old_chain,
                  replay_max_absolute_norm_difference=float(np.max(abs(old_r-series['r']))),
                  replay_max_absolute_width_difference=float(np.max(abs(old_w-w))),
                  replay_max_unit_vector_angle_deg=float(np.rad2deg(np.max(arcs(np.stack([old_u[70],series['u'][70]]))))),
                  new_exact_global_optimizer_length=stats['arc_geodesic'])
                parent_checks.append(hrow)
                if curve=='fixed_width' and hrow['replay_max_absolute_norm_difference']>1e-9:
                    raise RuntimeError('HISTORICAL_FIXED_WIDTH_REPLAY_MISMATCH')
                if curve=='optimized_valley' and hrow['replay_max_absolute_norm_difference']>0.01:
                    raise RuntimeError('HISTORICAL_OPTIMIZED_CURVE_DRIFT_EXCESS')
                print('CASE',label,curve,'historical_chord_minus_chain',format(old_arc-old_chain,'.9g'),
                      'new_geodesic_length',format(stats['arc_geodesic'],'.9g'),flush=True)
          print('MODE_DONE',label,mode,flush=True)
    csv_write('arclength_convergence_summary.csv',summary)
    csv_write('width_and_direction_samples.csv',samples)
    csv_write('historical_v18_53o_replay.csv',parent_checks)
    np.savez_compressed(OUT/'unit_direction_and_oriented_area.npz',**vectors)
    # Compare all resolutions while preserving interpolant, case and curve.
    ref=[]
    for s in summary:
       if s['n_points']!=points[0]:continue
       peers=[r for r in summary if all(r[k]==s[k] for k in ['seed','kappa_rule','interpolant','curve'])]
       peers.sort(key=lambda x:x['n_points'])
       for a,b in zip(peers,peers[1:]):
           ref.append(dict(seed=s['seed'],kappa_rule=s['kappa_rule'],interpolant=s['interpolant'],curve=s['curve'],
             n_low=a['n_points'],n_high=b['n_points'],
             geodesic_length_change=float(b['arc_geodesic']-a['arc_geodesic']),
             chain_discrepancy_low=a['arc_geodesic_minus_chain'],
             chain_discrepancy_high=b['arc_geodesic_minus_chain'],
             oriented_polygon_norm_change=float(b['area_polygon_norm']-a['area_polygon_norm']),
             signed_arc_endpoints_change=float(b['geodesic_endpoint']-a['geodesic_endpoint']),
             min_norm_change=float(b['min_norm']-a['min_norm'])))
    csv_write('resolution_pair_comparison.csv',ref)
    # Explicitly call out curves where the coarse-grid discrepancy persists.
    high=max(points)
    highest=[r for r in summary if r['n_points']==high]
    report=dict(version='v18.53p',scientific_decision='AUTHENTIC_V18_53O_DERIVED_GEOMETRIC_CONVERGENCE__NO_NEW_IVPS__NO_MATERIAL_HOLONOMY',
      source_parent_sha256=CFG['parent_v18_53o_sha256'],parent_manifests=pins,
      inherited_physical_ivps=10,new_physical_ivps=0,physical_refinement_reexecuted=0,
      num_case_combinations=6,num_curve_summaries=len(summary),num_sample_rows=len(samples),
      grid_resolutions=points,interpolants=CFG['interpolants'],
      historical_valley_direct_minus_chain_rad_range=[float(min(x['historical_chord_minus_chain'] for x in parent_checks if x['curve']=='optimized_valley')),
         float(max(x['historical_chord_minus_chain'] for x in parent_checks if x['curve']=='optimized_valley'))],
      finest_pl_valley_chord_minus_chain_range=[float(min(r['arc_geodesic_minus_chain'] for r in highest if r['interpolant']=='linear' and r['curve']=='optimized_valley')),
         float(max(r['arc_geodesic_minus_chain'] for r in highest if r['interpolant']=='linear' and r['curve']=='optimized_valley'))],
      finest_pchip_valley_chord_minus_chain_range=[float(min(r['arc_geodesic_minus_chain'] for r in highest if r['interpolant']=='pchip' and r['curve']=='optimized_valley')),
         float(max(r['arc_geodesic_minus_chain'] for r in highest if r['interpolant']=='pchip' and r['curve']=='optimized_valley'))],
      all_endpoints_in_formal_clock_continuation=True,
      caution='PCHIP and piecewise-linear changes are interpolation choices of the SAME archived physical data. Non-smooth PL knots and a discrete argmin path can bias chain-rule tangent quadrature. No material holonomy or full 11D closure is established.')
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    with (OUT/'SCIENTIFIC_REPORT_SERVER.md').open('w') as f:
      f.write('# Q550 v18.53p — numerical post-processing audit\n\n')
      f.write('**' + report['scientific_decision'] + '**\n\n')
      f.write('Authenticated v18.53o physical history reused; **0 new physical IVPs**.\n\n')
      f.write('Historical o valley direct-minus-chain arc: '+str(report['historical_valley_direct_minus_chain_rad_range'])+' rad.\n\n')
      f.write('At finest resolution '+str(high)+', piecewise linear direct-minus-chain: '+str(report['finest_pl_valley_chord_minus_chain_range'])+' rad.\n\n')
      f.write('At finest resolution '+str(high)+', PCHIP direct-minus-chain: '+str(report['finest_pchip_valley_chord_minus_chain_range'])+' rad.\n\n')
      f.write('See resolution and width/knot tables before interpreting any change as dynamics. No physical state closure, new IVP, or holonomy certificate.\n')
    try:
      import matplotlib
      matplotlib.use('Agg')
      import matplotlib.pyplot as plt
      for seed in ['original','autonomous']:
        fig,ax=plt.subplots(2,1,figsize=(8,7))
        for mode in CFG['interpolants']:
          d=[r for r in summary if r['seed']==seed and r['kappa_rule']=='original' and r['interpolant']==mode and r['curve']=='optimized_valley']
          ax[0].plot([r['n_points'] for r in d],[r['arc_geodesic_minus_chain'] for r in d],marker='o',label=mode)
          high_series=[r for r in samples if r['seed']==seed and r['kappa_rule']=='original' and r['interpolant']==mode and r['curve']=='optimized_valley' and r['n_points']==high]
          ax[1].plot([r['origin_xi'] for r in high_series],[r['width_xi'] for r in high_series],label=mode)
        ax[0].set(xlabel='Number of origin samples',ylabel='Geodesic − chain-rule length (rad)',title='v18.53p derived curve consistency — '+seed)
        ax[1].set(xlabel='Origin ξ',ylabel='Optimized three-window width w')
        for a in ax:a.legend();a.grid(alpha=.25)
        fig.tight_layout();fig.savefig(OUT/f'valley_arclength_convergence_{seed}.png',dpi=130);plt.close(fig)
    except Exception as e:
      (OUT/'PLOT_WARNING.txt').write_text(repr(e)+'\n')
    print('V1853P_POSTPROCESS_COMPLETE summaries',len(summary),'samples',len(samples),'new_ivps=0',flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--preflight',action='store_true');ap.add_argument('--deep',action='store_true')
    a=ap.parse_args()
    if a.preflight:
        data,f,v,s,pins,report=authenticated_o()
        assert len(data['common_xi'])>1000 and len(f)==len(v)==846
        print('V1853P_GENEALOGY_PREFLIGHT_PASS',pins,'parent_ivps=10',flush=True)
        return
    run(a.deep)

if __name__=='__main__':
    try:main()
    except Exception as e:
        import traceback
        traceback.print_exc();sys.exit(2)
