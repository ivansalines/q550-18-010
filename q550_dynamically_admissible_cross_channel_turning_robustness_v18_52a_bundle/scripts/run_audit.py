#!/usr/bin/env python3
"""Q550 v18.52a — Physically admissible kappa-branch turning robustness.

This script reuses six independently integrated historical nonlinear trajectories,
then evaluates the ORIGINAL structural jet/Jacobian and partial-kappa source on
EACH trajectory. It solves a NEW forced tangent IVP for each branch. No channel
mixing, arbitrary forcing, or fabricated model perturbation is used.

The parameter band is the original-vs-autonomous historical kappa0 interval,
not a proof of robustness for other dynamical perturbations or longer times.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, io, json, platform, sys, tempfile, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results'
VERSION = 'v18.52a'
TITLE = 'Dynamically Admissible Cross-Channel Turning Robustness Audit'
SNAME = 'q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip'
ZNAME = 'q550_posture_elegance_observer_induced_catwalk_v18_51z_bundle_server_sendback.zip'
SROOT = 'q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
RROOT = 'q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
ZROOT = 'q550_posture_elegance_observer_induced_catwalk_v18_51z_bundle'
SEEDS = ('original','autonomous')
BRANCHES = ('autonomous','midpoint','original')
CHANNELS = ('v2','omega2','other','total')
MASK_START = 0.01


def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def load_import(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod)
    return mod


def wedge(a,b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def angular_budget(xi,v):
    """For a physical 2D plane: v has (n, 2, 4), ordered channels."""
    if v.ndim!=3 or v.shape[1:]!=(2,4):raise ValueError('Expected (N,2,4)')
    d=np.gradient(v,xi,axis=0,edge_order=2)
    total=v[:,:,3]
    den=np.sum(total*total,axis=1)
    cross=lambda i,j:wedge(v[:,:,i],d[:,:,j]) # shape (N,) 
    self_a=cross(0,0);self_b=cross(1,1);self_o=cross(2,2)
    mixed=cross(0,1)+cross(1,0)
    other_mixed=cross(0,2)+cross(2,0)+cross(1,2)+cross(2,1)
    full=wedge(total,d[:,:,3])
    valid=den>1e-28
    f=lambda x:np.divide(x,den,out=np.full(len(xi),np.nan),where=valid)
    self_ab=f(self_a+self_b);mix=f(mixed)
    total_rate=f(full)
    remainder=f(self_o+other_mixed)
    return {
        'rate_total':total_rate, 'rate_mixed':mix,
        'rate_self_ab':self_ab, 'rate_remaining':remainder,
        'rate_self_v2':f(self_a),'rate_self_omega2':f(self_b),
        'wedge_reconstruction_max_abs':float(np.nanmax(np.abs(total_rate-(mix+self_ab+remainder)))),
        'oriented_area':wedge(v[:,:,0],v[:,:,1]),
        'dot':np.sum(v[:,:,0]*v[:,:,1],axis=1),
        'channel_norms':np.linalg.norm(v[:,:,:3],axis=1),
    }


def find_input(name):
    candidates=[ROOT/'inputs'/name, ROOT.parent/name, Path('/root/elab1/hysteresis')/name]
    matches=[p for p in candidates if p.is_file()]
    if not matches:raise FileNotFoundError(f'Missing authenticated parent {name}; place in inputs/ or the parent workspace')
    return matches[0]


def authenticate_zip(path, expected_sha, root, version):
    if sha256_file(path)!=expected_sha:raise RuntimeError(f'Parent bytes drifted: {path.name}')
    with zipfile.ZipFile(path) as z:
        if z.testzip():raise RuntimeError(f'Parent CRC failure: {path.name}')
        report=json.loads(z.read(root+'/results/report.json'))
        assert report['version']==version and report['all_passed'], (root,report.get('version'))
    return report


def load_parents():
    expected=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())['parents']
    s_path=find_input(SNAME); z_path=find_input(ZNAME)
    authenticate_zip(s_path,expected['v18_51s_sha256'],SROOT,'v18.51s')
    authenticate_zip(z_path,expected['v18_51z_sha256'],ZROOT,'v18.51z')
    with zipfile.ZipFile(s_path) as z:
        rblob=z.read(SROOT+'/inputs/'+RROOT+'_server_sendback.zip')
        srun=z.read(SROOT+'/scripts/run_audit.py')
    with zipfile.ZipFile(io.BytesIO(rblob)) as z:
        rr=json.loads(z.read(RROOT+'/results/report.json'))
        assert rr['version']=='v18.51r' and all(rr['checks'].values())
        raw=z.read(RROOT+'/results/extended_paired_kappa_nonlinear_states.npz')
    with np.load(io.BytesIO(raw),allow_pickle=False) as obj:
        states={k:obj[k].copy() for k in obj.files}
    with zipfile.ZipFile(z_path) as z:
        rawz=z.read(ZROOT+'/results/posture_observer_turning_arrays.npz')
    with np.load(io.BytesIO(rawz),allow_pickle=False) as obj:
        baseline={k:obj[k].copy() for k in obj.files}
    xi=states['xi']
    assert len(xi)==513 and abs(float(xi[-1])-.2)<1e-13
    assert np.array_equal(xi,baseline['xi'])
    return s_path,z_path,rblob,srun,states,baseline


def sample_jet(xi,Y,kappa,ctx,engine,ids,h=1e-5):
    """Structural tangent operator and physical parameter derivative at fixed state."""
    J=[];B=[];source_check=[];rhs_max=0.
    for i in ids:
        yy=np.asarray(Y[i],float)
        tau_event=float(ctx['event']+xi[i]);zz=np.r_[yy,float(ctx['theta'](tau_event))]
        ff,jj,_=engine.jet.structural_jet(zz,kappa,engine.params,engine.slope)
        j=np.asarray(jj[:10,:10],float)
        def dk(hh):
            fp=np.asarray(engine.jet.structural_jet(zz,kappa+hh,engine.params,engine.slope)[0][:10],float)
            fm=np.asarray(engine.jet.structural_jet(zz,kappa-hh,engine.params,engine.slope)[0][:10],float)
            return (fp-fm)/(2*hh)
        bb=dk(h); b_coarse=dk(2*h)
        check=float(np.linalg.norm(bb-b_coarse)/max(1e-6,np.linalg.norm(bb)))
        if not np.all(np.isfinite(j)) or not np.all(np.isfinite(bb)):
            raise RuntimeError('Nonfinite authenticated structural jet')
        source_check.append(check);J.append(j);B.append(bb)
        rhs_max=max(rhs_max,float(np.linalg.norm(np.asarray(ff[:10]))))
    return np.asarray(J),np.asarray(B),max(source_check),rhs_max


def integrate_sources(xi,ids,J,b,engine_reference):
    jj=CubicSpline(xi[ids],J,axis=0)
    bb=CubicSpline(xi[ids],b,axis=0)
    def rhs(x,flat):
        s=flat.reshape(10,4);src=np.asarray(bb(x))
        B=np.column_stack([engine_reference.source_channel(src,s) for s in ('v2','omega2','other')]+[src])
        return (jj(x)@s+B).ravel()
    sol=solve_ivp(rhs,(xi[0],xi[-1]),np.zeros(40),method='DOP853',t_eval=xi,
                  rtol=3e-11,atol=3e-13,max_step=.001)
    if not sol.success:raise RuntimeError('Forced tangent IVP failed: '+sol.message)
    return sol.y.T.reshape(-1,10,4),int(sol.nfev)


def analyze_branch(seed,branch,xi,Y,Q,transport,bs,ref):
    s,nfev=transport
    c=np.stack([ref.j_immersion(y,q)@ss for y,q,ss in zip(Y,Q,s)],axis=0)
    c2=c[:,2:4,:]
    bud=angular_budget(xi,c2)
    valid=(xi>=MASK_START) & np.isfinite(bud['rate_total'])
    self_term=bud['rate_self_ab']; mixed=bud['rate_mixed']; rest=bud['rate_remaining']
    # The exact signed margin for v2/omega2 self-term dominance is mixed + self_ab.
    margin=mixed+self_term
    total=bud['rate_total']
    angle=np.arctan2(bud['oriented_area'],bud['dot'])
    norm=max(1e-20,np.linalg.norm(c[-1,2:4,3]))
    endpoint=dict(
        mixed=float(mixed[-1]),self_ab=float(self_term[-1]),
        other=float(rest[-1]),total=float(total[-1]),
        dominance_margin=float(margin[-1]),
        mutual_angle_deg=float(np.rad2deg(angle[-1])),
        response_norm=float(norm),
        relative_mixed_strength=float(mixed[-1]/max(1e-20,abs(self_term[-1]))),
    )
    summ={
       'seed':seed,'kappa_branch':branch,'n_qualified_sections':int(np.count_nonzero(valid)),
       'kappa_value':None,'source_max_off_position_absolute':float(np.max(np.linalg.norm(bs[:,[0,2,4,6,8]],axis=1))),
       'max_response_channel_additivity':float(np.max(np.linalg.norm(s[:,:,:3].sum(axis=2)-s[:,:,3],axis=1))),
       'max_turning_budget_residual':bud['wedge_reconstruction_max_abs'],
       'min_mixed_positive_margin':float(np.min(margin[valid])),
       'min_total_turning_rate':float(np.min(total[valid])),
       'max_total_turning_rate':float(np.max(total[valid])),
       'n_sections_mixed_margin_nonpositive':int(np.count_nonzero(margin[valid]<=0)),
       'n_sections_total_turning_nonpositive':int(np.count_nonzero(total[valid]<=0)),
       'n_sections_self_ab_nonnegative':int(np.count_nonzero(self_term[valid]>=0)),
       'fraction_positive_channel_area':float(np.mean(bud['oriented_area'][valid]>0)),
       'endpoint':endpoint,'tangent_evaluations':nfev,
    }
    return summ,bud,c2,s


def plot(all_data,ks,xi):
    fig,axs=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
    style={'autonomous':'--','midpoint':'-','original':':'}
    col={'original':'#2b70aa','autonomous':'#d16843'}
    valid=xi>=MASK_START
    for seed in SEEDS:
        for br in BRANCHES:
            d=all_data[(seed,br)];b=d['budget'];l=f'{seed} {br}'
            opts=dict(color=col[seed],linestyle=style[br],label=l,linewidth=1.45)
            axs[0,0].plot(xi[valid],b['rate_mixed'][valid],**opts)
            axs[0,1].plot(xi[valid],(b['rate_mixed']+b['rate_self_ab'])[valid],**opts)
            axs[1,0].plot(xi[valid],b['rate_total'][valid],**opts)
            ang=np.degrees(np.arctan2(b['oriented_area'],b['dot']))
            axs[1,1].plot(xi[valid],ang[valid],**opts)
    names=['Cross-channel turning contribution','Signed mixed-vs-self dominance margin',
           'Total response turning rate','Mutual response-vector angle']
    ys=['rad / xi','rad / xi','rad / xi','degrees']
    for ax,n,y in zip(axs.flat,names,ys):
        ax.set(xlabel='xi',ylabel=y,title=n);ax.grid(alpha=.23);ax.axhline(0,color='gray',lw=.6)
        ax.legend(fontsize=7,ncol=2,loc='best')
    fig.suptitle('Q550 v18.52a — Actual kappa branch dynamics and cross-channel turning')
    fig.savefig(OUT/'physical_kappa_branch_cross_channel_turning_robustness.png',dpi=160)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=TITLE)
    parser.add_argument('--skip-green',action='store_true',help='Skip optional baseline Green quadrature (not required for branch survey)')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    s_path,z_path,rblob,srun,states,baseline=load_parents()
    xi=states['xi']; k_orig=float(states['kappa_original']);k_auto=float(states['kappa_autonomous'])
    ks=dict(autonomous=k_auto,midpoint=.5*(k_orig+k_auto),original=k_orig)
    print('Historical physical kappa band:',ks,flush=True)
    ref=load_import(ROOT/'reference'/'parent_v18_51t_run_audit.py','q550_v1852a_parent_t')
    records=[];out={};all_checks={};data={};refinements={};source_consistency={}
    with tempfile.TemporaryDirectory(prefix='q550_v1852a_historical_model_') as td:
        print('AUTHENTICATED PHYSICAL FIELD: recovering historical 19-ancestor engine',flush=True)
        _,_,seed_api,observer,engine,unit,drpt,creport,_ = ref.stage_engine(td,rblob,srun)[4]
        for seed in SEEDS:
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            for branch in BRANCHES:
                k=ks[branch]
                Y=states[f'{seed}_{branch}_Y'];Q=states[f'{seed}_{branch}_Q']
                ids=np.arange(0,len(xi),2)
                print(f'  physical branch {seed}/{branch}: kappa={k:.12f}, jets={len(ids)}',flush=True)
                J,b,ds,rhs_max=sample_jet(xi,Y,k,ctx,engine,ids)
                forced,nfev=integrate_sources(xi,ids,J,b,ref)
                summ,bud,plane,s=analyze_branch(seed,branch,xi,Y,Q,(forced,nfev),b,ref)
                summ['kappa_value']=k;summ['max_source_h_refinement_relative']=ds
                summ['max_nonlin_rhs_norm_at_jet_nodes']=rhs_max
                # Model-defined fixed-support channels are not arbitrary mixtures.
                all_checks[f'{seed}_{branch}_linear_recomposition']=summ['max_response_channel_additivity']<1e-10
                all_checks[f'{seed}_{branch}_angular_budget_recomposition']=summ['max_turning_budget_residual']<1e-7
                all_checks[f'{seed}_{branch}_physical_source_position_support']=summ['source_max_off_position_absolute']<1e-7
                all_checks[f'{seed}_{branch}_parameter_derivative_step_consistency']=ds<1e-3
                out[(seed,branch)]={'summary':summ,'budget':bud,'plane':plane}
                source_consistency[(seed,branch)]=ds
                for i in np.where(xi>=MASK_START)[0]:
                    records.append({'seed':seed,'branch':branch,'kappa':k,'index':int(i),'xi':float(xi[i]),
                         'turn_mixed':float(bud['rate_mixed'][i]),'turn_self_v2':float(bud['rate_self_v2'][i]),
                         'turn_self_omega2':float(bud['rate_self_omega2'][i]),
                         'turn_self_ab':float(bud['rate_self_ab'][i]),
                         'turn_residual_channels':float(bud['rate_remaining'][i]),
                         'turn_total':float(bud['rate_total'][i]),
                         'signed_cross_dominance_margin':float((bud['rate_mixed']+bud['rate_self_ab'])[i]),
                         'mutual_angle_deg':float(np.rad2deg(np.arctan2(bud['oriented_area'][i],bud['dot'][i]))),
                         'oriented_area':float(bud['oriented_area'][i])})
                if branch=='midpoint':
                    b0=baseline[f'{seed}_lab_channel_vectors']
                    disagreement=float(np.max(np.abs(plane-b0)))
                    all_checks[f'{seed}_independent_Z_midpoint_replay']=disagreement<1e-7
                    summ['max_midpoint_Z_cartesian_channel_difference']=disagreement
                    # The same physical J and b, sampled more sparsely, form an independent grid test.
                    coarse=integrate_sources(xi,ids[::2],J[::2],b[::2],ref)[0]
                    refinement=float(np.max(np.linalg.norm(coarse[20:,:,3]-forced[20:,:,3],axis=1)
                         /np.maximum(1e-4,np.linalg.norm(forced[20:,:,3],axis=1))))
                    refinements[seed]=refinement
                    all_checks[f'{seed}_independent_jet_grid_refinement']=refinement<2e-3
                    summ['max_stride_four_vs_two_relative_difference']=refinement
                    if not args.skip_green:
                        # Independent propagator/Green quadrature only for the canonical midpoint.
                        g=ref.transport(xi,ids,J,b)
                        good=float(np.max(np.linalg.norm(g['green'][2::2,:,3]-forced[2::2,:,3],axis=1)
                                   /np.maximum(1e-4,np.linalg.norm(forced[2::2,:,3],axis=1))))
                        all_checks[f'{seed}_independent_midpoint_Green']=good<1e-3
                        summ['max_midpoint_green_quadrature_relative_difference']=good
                    print(f'    midpoint replay vs Z {disagreement:.3e}; grid-refine {refinement:.3e}',flush=True)
                print(f'    endpoint rates mixed {summ["endpoint"]["mixed"]:+.6f} self {summ["endpoint"]["self_ab"]:+.6f} total {summ["endpoint"]["total"]:+.6f}',flush=True)
        # End of physical engine staging.
    all_checks['three_authentic_distinct_kappa_branches']=k_auto<ks['midpoint']<k_orig
    all_checks['two_original_seed_preparations']=float(np.linalg.norm(states['original_midpoint_Y'][0]-states['autonomous_midpoint_Y'][0]))>1e-10
    all_checks['same_seed_same_initial_data_across_kappa']=all(
        np.linalg.norm(states[f'{seed}_original_Y'][0]-states[f'{seed}_autonomous_Y'][0])<1e-12 for seed in SEEDS)
    all_checks['fixed_horizon_and_sample_grid']=len(xi)==513 and abs(float(xi[-1])-.2)<1e-13
    all_checks['all_outputs_finite']=all(np.all(np.isfinite([r['turn_total'],r['turn_mixed']])) for r in records)
    all_checks['two_lineage_hashes_authenticated']=True
    print('WRITING RESULTS',flush=True)
    summary={f'{seed}_{br}':out[(seed,br)]['summary'] for seed in SEEDS for br in BRANCHES}
    # Scientific decisions are empirical RESULTS, not pass/fail criteria for the runner.
    physical_robust={
      'mixed_beats_negative_self_all_branches':all(x['n_sections_mixed_margin_nonpositive']==0 for x in summary.values()),
      'total_response_counterclockwise_all_branches':all(x['n_sections_total_turning_nonpositive']==0 for x in summary.values()),
      'both_self_channels_negative_everywhere_all_branches':all(x['n_sections_self_ab_nonnegative']==0 for x in summary.values()),
      'positive_oriented_channel_area_all_branches':all(x['fraction_positive_channel_area']==1. for x in summary.values())}
    decision=('SUPPORTED_IN_AUTHENTIC_THREE_KAPPA_BRANCH_BAND' if all(physical_robust.values())
         else 'NOT_UNIFORMLY_SUPPORTED_IN_AUTHENTIC_THREE_KAPPA_BRANCH_BAND')
    with (OUT/'physical_parameter_robustness_sections.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    arrays={'xi':xi,'kappa_original':np.array(k_orig),'kappa_autonomous':np.array(k_auto)}
    for seed in SEEDS:
        for branch in BRANCHES:
            obj=out[(seed,branch)]
            arrays[f'{seed}_{branch}_plane2_channels']=obj['plane']
            for key in ('rate_total','rate_mixed','rate_self_ab','rate_remaining','oriented_area','dot'):
                arrays[f'{seed}_{branch}_{key}']=obj['budget'][key]
    np.savez_compressed(OUT/'physical_branch_turning_arrays.npz',**arrays)
    report={
       'version':VERSION,'title':TITLE,'all_passed':all(all_checks.values()),
       'numerical_checks':all_checks,'n_numerical_checks':len(all_checks),
       'scientific_decision':decision,'scientific_robustness_tests':physical_robust,
       'classification':'MODEL_PARAMETER_BAND_CONDITIONAL_LOCAL_CROSS_CHANNEL_TURNING_ROBUSTNESS_DIAGNOSTIC',
       'initial_conditions':'Same initial state across kappa within each seed, original and autonomous independently prepared seeds.',
       'parameter_protocol':{'autonomous':k_auto,'midpoint':ks['midpoint'],'original':k_orig,
            'total_kappa_band':k_orig-k_auto,'xi_interval':[0.,.2],
            'qualified_turning_interval':[MASK_START,.2],
            'physical_perturbation':'Historical parameter kappa0 in original nonlinear 10D field; J and d_kappa f recomputed on saved nonlinear paths.',
            'source_channels':'Original fixed supports velocity state indices 5,7,with unaltered other residual; no source-channel recombination.'},
       'parent_z_sha256':sha256_file(z_path),'parent_s_sha256':sha256_file(s_path),
       'branches':summary,'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__},
       'limitations':[
          'Only the historical original-to-autonomous kappa0 band and the existing two seed preparations are tested; no general admissible-perturbation ball is certified.',
          'The previously integrated nonlinear trajectories are reused; only their Jacobian, parameter source and forced channel responses are recomputed.',
          'The signed turning decomposition is coordinate/channel dependent; source supports are anchored to the original model equations.',
          'The model does not establish a physical material-return section on xi in [0,0.2].',
          'Discrete sample signs and ordinary numerical convergence are not rigorous interval bounds on sign preservation between samples.',
          'The scientific robust/nonrobust decision is recorded separately from numerical audit validity.'
       ]}
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    plot(out,ks,xi)
    report_lines=[f'# Q550 {VERSION}: {TITLE}','',f'**Numerical audit:** {sum(all_checks.values())}/{len(all_checks)} PASS',
                  f'**Scientific decision:** `{decision}`','', '## Model-defined physical perturbation','',
                  f'- Three already integrated kappa0 branches: {k_auto:.15f}, {ks["midpoint"]:.15f}, {k_orig:.15f}.',
                  '- Each branch uses its own archived nonlinear Y(t), current structural Jacobian J(t), and fresh parameter source partial_kappa f(t).',
                  '- Independent channel responses are forced tangent IVP solutions; no channel mixing is used.','',
                  '## Endpoint signed turning budget (rad / xi)','',
                  '| Seed | Branch | Mixed | Self v2+omega2 | Other terms | Total | Min margin |','|---|---|---:|---:|---:|---:|---:|']
    for seed in SEEDS:
        for br in BRANCHES:
            x=summary[f'{seed}_{br}'];e=x['endpoint']
            report_lines.append(f'| {seed} | {br} | {e["mixed"]:+.8g} | {e["self_ab"]:+.8g} | {e["other"]:+.8g} | {e["total"]:+.8g} | {x["min_mixed_positive_margin"]:+.8g} |')
    report_lines.extend(['','## Scientific scope','', *('- '+x for x in report['limitations']),''])
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(report_lines))
    print('NUMERICAL CHECKS',sum(all_checks.values()),'/',len(all_checks),'PASS' if report['all_passed'] else 'FAIL',flush=True)
    print('SCIENTIFIC DECISION:',decision,flush=True)
    if not report['all_passed']:
        print('FAILED CHECKS',[k for k,v in all_checks.items() if not v]);raise SystemExit(2)


if __name__=='__main__':main()
