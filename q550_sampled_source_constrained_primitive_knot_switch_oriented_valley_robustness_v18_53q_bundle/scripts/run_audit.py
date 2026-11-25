#!/usr/bin/env python3
"""v18.53q source-constrained valley analysis using authenticated archived Q550 samples."""
import argparse,csv,hashlib,io,json,sys,zipfile
from pathlib import Path
import numpy as np
from source_valley_math import Primitive,optimized_width,curve
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs'/CFG['parent_filename']
sha=lambda b:hashlib.sha256(b).hexdigest()
def file_sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for t in iter(lambda:f.read(1048576),b''):h.update(t)
    return h.hexdigest()
def authenticate(z,expected,report_name):
    if z.testzip() is not None:raise ValueError('CRC_FAIL '+report_name)
    names=z.namelist()
    if any(Path(p).is_absolute() or '..' in Path(p).parts for p in names):raise ValueError('UNSAFE_ZIP_MEMBER')
    roots={p.split('/')[0] for p in names if '/' in p}
    if len(roots)!=1:raise ValueError('MULTIPLE_ROOTS')
    prefix=next(iter(roots))+'/'
    counts={}
    for label in ['SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json']:
        man=json.loads(z.read(prefix+label)); files=man['files']
        for path,dig in files.items():
            if sha(z.read(prefix+path))!=dig:raise ValueError(f'{report_name}_{label}_SHA_MISMATCH {path}')
        counts[label]=len(files)
    report=json.loads(z.read(prefix+'results/report.json'))
    if report.get('scientific_decision')!=expected:raise ValueError(report_name+'_DECISION_MISMATCH')
    return prefix,counts,report

def authenticated():
    if file_sha(PARENT)!=CFG['parent_sha256']:raise ValueError('PARENT_P_ARCHIVE_SHA_MISMATCH')
    with zipfile.ZipFile(PARENT) as z:
        pp,pc,pr=authenticate(z,CFG['parent_required_decision'],'P')
        if pr.get('new_physical_ivps')!=0 or pr.get('inherited_physical_ivps')!=10:raise ValueError('P_LINEAGE_WRONG_IVP_COUNT')
        inner_name=pp+'inputs/v18_53o_authenticated_server_sendback.zip'
        inner_bytes=z.read(inner_name)
        pins=json.loads(z.read(pp+'inputs_config.json'))
        if sha(inner_bytes)!=pins['parent_v18_53o_sha256']:raise ValueError('O_SHA_MISMATCH_FROM_P')
    with zipfile.ZipFile(io.BytesIO(inner_bytes)) as o:
        op,oc,orr=authenticate(o,CFG['parent_o_required_decision'],'O')
        if orr.get('actual_Q550_ivps')!=10 or orr.get('refined_Q550_ivps')!=4:raise ValueError('O_NOT_AUTHENTIC_DYNAMICS')
        with np.load(io.BytesIO(o.read(op+'results/derived_primitives.npz')),allow_pickle=False) as npz:
            data={k:np.asarray(npz[k]).copy() for k in npz.files}
        orig_rows=list(csv.DictReader(io.StringIO(o.read(op+'results/fixed_width_oriented_transport.csv').decode())))
    return data,orig_rows,{'p':pc,'o':oc,'sha_p':CFG['parent_sha256'],'sha_o':sha(inner_bytes)}

def save_csv(name,rows):
    if not rows:raise ValueError('EMPTY_RESULTS_'+name)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def run(maxcases=6,points=None,deep=False,save=True):
    data,old,pins=authenticated()
    xi=data['common_xi'];keys=[(s,k) for s in ('original','autonomous') for k in ('shared','original','autonomous')][:maxcases]
    resolutions=points or CFG['grid_points']; nsearch=CFG['deep_search_grid_points'] if deep else CFG['search_grid_points']
    summary=[];raw=[];source=[];vectors={};comparisons=[]
    for seed,rule in keys:
        case=seed+'__'+rule
        cr=sorted([x for x in old if x['seed']==seed and x['kappa_rule']==rule],key=lambda x:int(x['index']))
        if len(cr)!=141:raise ValueError('PARENT_O_CURVE_MISSING '+case)
        center_o=float(cr[70]['origin_xi']);center_w=float(cr[70]['width_xi']);lo=center_w-CFG['width_halfband'];hi=center_w+CFG['width_halfband']
        C=data[case+'__primitive'];g=data[case+'__sampled_source_10d']*data[case+'__clock_rate'][:,None]
        if not np.all(np.isfinite(g)):raise ValueError('SOURCE_NAN')
        for mode in ('pchip','hermite','source_integral'):
            p=Primitive(xi,C,g,mode)
            source.append(dict(seed=seed,kappa_rule=rule,mode=mode,
                max_source_knot_derivative_mismatch=float(np.max(p.node_source_derivative_mismatch())),
                max_nodal_primitive_deviation=float(np.max(p.node_mismatch())),
                endpoint_nodal_primitive_deviation=float(p.node_mismatch()[-1]),
                primitive_samples=len(xi),grid_step_min=float(np.min(np.diff(xi))),grid_step_max=float(np.max(np.diff(xi)))))
            for count in resolutions:
                ox=np.linspace(center_o-CFG['origin_halfspan'],center_o+CFG['origin_halfspan'],count)
                if ox[0]<xi[0] or ox[-1]+3*hi>xi[-1]:raise ValueError('EVALUATION_OUTSIDE_TRACE')
                mins=[optimized_width(p,float(o),lo,hi,nsearch=nsearch) for o in ox]
                w=np.array([x['w'] for x in mins]); wp=np.array([x['implicit_wprime'] for x in mins]);
                curve_result=curve(p,ox,w,wp)
                if curve_result['min_norm']<CFG['min_norm_floor']:raise ValueError('UNEXPECTED_NEAR_ZERO')
                nmask=int(np.sum(~np.isfinite(wp)));bounds=int(sum(x['bound'] for x in mins))
                switches=int(sum(mins[j]['label']!=mins[j-1]['label'] for j in range(1,len(mins))))
                base=dict(seed=seed,kappa_rule=rule,mode=mode,n_points=count,origin_center=center_o,
                          width_center=center_w,origin_start=float(ox[0]),origin_end=float(ox[-1]),
                          width_bounds_count=bounds,mask_implicit_derivative_count=nmask,
                          candidate_label_switches=switches,
                          min_stationary_curvature=float(min(x['curvature'] for x in mins)),
                          min_candidate_score_gap=float(min(x['score_margin'] for x in mins)),
                          min_endpoint_knot_clearance=float(min(x['clearance'] for x in mins)),
                          max_absolute_stationarity_in_interior=float(max(abs(x['stationarity']) for x in mins if not x['bound']) if bounds<len(mins) else 0),
                          all_after_formal_clock_boundary=bool(np.all(ox+3*w>CFG['formal_continuation_boundary_xi'])))
                stats={k:v for k,v in curve_result.items() if k not in ['u','norm','wd','directional_speed','oriented_polygon','oriented_tangent']}
                summary.append({**base,**stats})
                for i,(o,m) in enumerate(zip(ox,mins)):
                    raw.append(dict(seed=seed,kappa_rule=rule,mode=mode,n_points=count,index=i,
                       origin_xi=float(o),width_xi=m['w'],end_xi=float(o+3*m['w']),
                       net_norm=m['norm'],minimum_kind=m['label'],num_stationary_candidates=m['num_candidates'],
                       score_margin=m['score_margin'],curvature=m['curvature'],
                       implicit_width_derivative=m['implicit_wprime'],stationarity=m['stationarity'],
                       endpoint_knot_clearance=m['clearance'],directional_speed_10d=float(curve_result['directional_speed'][i]),width_slope_used=float(curve_result['wd'][i])))
                tag=f'{case}__{mode}__{count}'
                vectors[tag+'__unit10d']=curve_result['u'];vectors[tag+'__origin']=ox;vectors[tag+'__width']=w
                vectors[tag+'__oriented_polygon45d']=curve_result['oriented_polygon']
                vectors[tag+'__oriented_tangent45d']=curve_result['oriented_tangent']
                print('V1853Q_CASE',case,mode,count,'direct',curve_result['direct'],'chain',curve_result['chain'],flush=True)
        for count in resolutions:
            subset=[x for x in summary if x['seed']==seed and x['kappa_rule']==rule and x['n_points']==count]
            par={x['mode']:x for x in subset}
            for a,b in [('hermite','pchip'),('source_integral','pchip'),('source_integral','hermite')]:
                ua=vectors[f'{case}__{a}__{count}__unit10d'];ub=vectors[f'{case}__{b}__{count}__unit10d']
                aa=vectors[f'{case}__{a}__{count}__oriented_polygon45d'];ab=vectors[f'{case}__{b}__{count}__oriented_polygon45d']
                wa=vectors[f'{case}__{a}__{count}__width'];wb=vectors[f'{case}__{b}__{count}__width']
                dot=np.clip(np.sum(ua*ub,axis=1),-1,1)
                comparisons.append(dict(seed=seed,kappa_rule=rule,n_points=count,mode_a=a,mode_b=b,
                    max_width_abs_difference=float(np.max(abs(wa-wb))),
                    max_unit_direction_angle_rad=float(np.max(np.arccos(dot))),
                    oriented_area_angle_deg=float(np.rad2deg(np.arccos(np.clip(aa@ab/max(np.linalg.norm(aa)*np.linalg.norm(ab),1e-100),-1,1)))),
                    geodesic_length_difference=float(par[a]['direct']-par[b]['direct'])))
    if not save:return dict(num_cases=len(keys),summary=summary,source=source,comparisons=comparisons,pins=pins)
    save_csv('convergence_summary.csv',summary);save_csv('valley_samples.csv',raw)
    save_csv('primitive_source_consistency.csv',source);save_csv('cross_interpolant_comparison.csv',comparisons)
    np.savez_compressed(OUT/'oriented_forms.npz',**vectors)
    report=dict(version='v18.53q',scientific_decision='AUTHENTIC_SAMPLED_SOURCE_CONSTRAINED_INTERPOLANT_DISCRIMINATION__NO_NEW_IVPS__NO_MATERIAL_HOLONOMY',
        inherited_physical_ivps=10,new_physical_ivps=0,physical_refinement_reexecuted=0,
        parent_genealogy=pins,cases=len(keys),num_summaries=len(summary),num_samples=len(raw),
        methods=['pchip_cumulative','cubic_hermite_matched_sampled_source','integral_of_piecewise_linear_sampled_source'],
        max_source_integrated_nodal_drift=float(max(x['max_nodal_primitive_deviation'] for x in source if x['mode']=='source_integral')),
        max_cross_interpolant_width_difference=float(max(x['max_width_abs_difference'] for x in comparisons)),
        maximum_cross_interpolant_oriented_area_angle_deg=float(max(x['oriented_area_angle_deg'] for x in comparisons)),
        max_direct_chain_discrepancy=float(max(abs(x['arc_diff']) for x in summary)),
        all_curves_after_formal_clock_boundary=all(x['all_after_formal_clock_boundary'] for x in summary),
        limitations=['No genuinely new sampling of the physical ODE occurs.','Sampled source and cumulant are from the same archived trajectory.','All near-return endpoints rely on formal xi>25.6 clock continuation.','No claim of 11D closure or material holonomy.','Width search brackets the observed stationary candidates; missed narrow double roots remain possible.'])
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 v18.53q — sampled-source-constrained audit\n\n**'+report['scientific_decision']+'**\n\n'
        +'Reused ten inherited physical integrations, performed zero new IVPs.\n\n'
        +f"Scanned {len(raw)} width-optimized origin samples over {len(summary)} case/method/resolution rows.\n\n"
        +f"Maximum sampled-source-integrated nodal mismatch against archived C: {report['max_source_integrated_nodal_drift']:.9g}.\n\n"
        +f"Maximum cross-interpolant width difference: {report['max_cross_interpolant_width_difference']:.9g}.\n\n"
        +f"Maximum cross-interpolant oriented-area angle: {report['maximum_cross_interpolant_oriented_area_angle_deg']:.9g} deg.\n\n"
        +'No material spin/holonomy certificate. Any remaining difference requires independent denser 11D source evolution on the user server.\n')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for seed in ('original','autonomous'):
            fig,ax=plt.subplots(2,1,figsize=(9,7),sharex=True)
            for mode in ('pchip','hermite','source_integral'):
                rows=[v for v in raw if v['seed']==seed and v['kappa_rule']=='original' and v['mode']==mode and v['n_points']==561]
                if not rows:continue
                xx=np.array([r['origin_xi'] for r in rows])
                ax[0].plot(xx,[r['width_xi'] for r in rows],label=mode)
                ax[1].plot(xx,[r['directional_speed_10d'] for r in rows],label=mode)
            ax[0].set_ylabel('Optimized width w(xi)')
            ax[1].set_ylabel('10D directional speed / xi')
            ax[1].set_xlabel('Window origin xi')
            for a in ax:a.grid(alpha=.25);a.legend()
            fig.suptitle('v18.53q — sampled-source-constrained valley: '+seed)
            fig.tight_layout();fig.savefig(OUT/f'source_constrained_valley_{seed}.png',dpi=140);plt.close(fig)
    except Exception as ex:
        (OUT/'PLOT_WARNING.txt').write_text(str(ex)+'\n')
    return report

def main():
    a=argparse.ArgumentParser();a.add_argument('--preflight',action='store_true');a.add_argument('--smoke',action='store_true');a.add_argument('--deep-search',action='store_true')
    args=a.parse_args()
    if args.preflight:
        d,o,p=authenticated();print('V1853Q_PREFLIGHT_PASS',p,'cases',len([x for x in d if x.endswith('__primitive')]));return
    if args.smoke:
        r=run(maxcases=1,points=[11],deep=False,save=False)
        print('V1853Q_SMOKE_PASS',len(r['summary']),'rows');return
    r=run(deep=args.deep_search)
    print('V1853Q_AUDIT_PASS',r['num_summaries'],'rows',r['max_source_integrated_nodal_drift'])
if __name__=='__main__':main()
