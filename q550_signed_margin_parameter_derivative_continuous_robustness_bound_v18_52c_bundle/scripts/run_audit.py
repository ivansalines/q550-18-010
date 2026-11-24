#!/usr/bin/env python3
"""Q550 v18.52c: physical signed-margin parameter slope and conditional continuum bound.

Numerical secants and PCHIP derivatives do NOT prove a bound on the exact
nonlinear model's parameter derivative. This distinction is enforced in all
reports, including after a successful physical midpoint replay.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, io, json, os, shutil, sys, tempfile, zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
VERSION='v18.52c'
TITLE='Signed-Margin Parameter Derivative & Continuous-Robustness Bound Audit'
BROOT='q550_signed_margin_continuation_physical_parameter_robustness_radius_v18_52b_bundle'
BNAME=BROOT+'_server_sendback.zip'
SNAME='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip'
SEEDS=('original','autonomous')
FIELDS=('margin','total','mixed','self_ab','oriented_area')


def sha_bytes(v):return hashlib.sha256(v).hexdigest()
def sha_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as fd:
        for block in iter(lambda:fd.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def source_manifest():return json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())

def find_input(name):
    places=[ROOT/'inputs'/name,ROOT.parent/name,Path('/root/elab1/hysteresis')/name]
    for p in places:
        if p.is_file():return p
    raise FileNotFoundError('Missing input: '+name+'; place in inputs/ or /root/elab1/hysteresis')


def verify_source_manifest():
    src=source_manifest()
    for rel,expected in src['files'].items():
        path=ROOT/rel
        if not path.is_file() or sha_file(path)!=expected:
            raise RuntimeError('Source manifest mismatch: '+rel)
    path=find_input(BNAME)
    if sha_file(path)!=src['parents']['v18_52b_sha256']:
        raise RuntimeError('Historical v18.52b sendback bytes changed')
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise RuntimeError('B sendback CRC failure')
        r=BROOT+'/'
        b_report=json.loads(z.read(r+'results/report.json'))
        if b_report['version']!='v18.52b' or not b_report['all_passed']:
            raise RuntimeError('Historical B did not pass')
        b_source=json.loads(z.read(r+'SOURCE_MANIFEST.json'))
        b_output=json.loads(z.read(r+'SHA256SUMS_OUTPUT.json'))
        for rel,expected in b_source['files'].items():
            if sha_bytes(z.read(r+rel))!=expected:raise RuntimeError('B internal source drift: '+rel)
        for rel,expected in b_output['file_sha256'].items():
            if sha_bytes(z.read(r+rel))!=expected:raise RuntimeError('B recorded output drift: '+rel)
        # Parent B output manifest is historically incomplete; explicitly
        # check all five genuine server result files and record their digests.
        required=['results/report.json','results/physical_parameter_continuation_scan.csv',
                  'results/physical_parameter_continuation_curves.npz',
                  'results/physical_kappa_signed_margin_continuation.png',
                  'results/SCIENTIFIC_REPORT_SERVER.md']
        extra_hash={rel:sha_bytes(z.read(r+rel)) for rel in required}
        with np.load(io.BytesIO(z.read(r+'results/physical_parameter_continuation_curves.npz')),
                     allow_pickle=False) as packed:
            arrays={key:packed[key].copy() for key in packed.files}
        rows=list(csv.DictReader(io.StringIO(z.read(r+'results/physical_parameter_continuation_scan.csv').decode('utf-8'))))
    if len(rows)!=26 or b_report['total_branch_records']!=26:
        raise RuntimeError('Expected 26 authenticated physical B branch records')
    return path,b_report,rows,arrays,extra_hash


def load_surface(rows,a,min_xi=.01):
    per={}
    for seed in SEEDS:
        rr=sorted((r for r in rows if r['seed']==seed),key=lambda r:float(r['kappa']))
        ks=np.array([float(r['kappa']) for r in rr],float)
        if len(rr)!=13 or not np.all(np.diff(ks)>0):raise RuntimeError('Parameter grid incomplete')
        xi=np.asarray(a[f'{seed}__{rr[0]["label"]}__xi'],float)
        if xi.shape!=(513,) or not np.all(np.diff(xi)>0):raise RuntimeError('Invalid physical time grid')
        mask=xi>=min_xi
        if np.count_nonzero(mask)!=487:raise RuntimeError('Unexpected qualified grid')
        data={}
        for field in FIELDS:
            names={'total':('total','rate_total'),
                   'mixed':('mixed','rate_mixed'),
                   'self_ab':('self_ab','rate_self_ab')}
            def getfield(r):
                pre=f'{seed}__{r["label"]}__'
                for n in names.get(field,(field,)):
                    if pre+n in a:return a[pre+n]
                raise KeyError('Missing both legacy and new naming for '+pre+field)
            vals=np.stack([getfield(r) for r in rr])
            if vals.shape!=(13,513):raise RuntimeError('Malformed field '+field)
            if not np.isfinite(vals[:,mask]).all():raise RuntimeError('Nonfinite historical '+field)
            data[field]=vals[:,mask]
        per[seed]={'k':ks,'xi':xi[mask],'data':data,'rows':rr}
    if not np.array_equal(per['original']['k'],per['autonomous']['k']):
        raise RuntimeError('Seed parameter grids are different')
    return per


def extrema_pchip(poly,knots):
    """Exact real extrema of each PCHIP polynomial piece in the kappa direction.

    PPoly coefficients are in ascending powers of local dx when reversed.
    Return min and max of value, and max absolute derivative.
    """
    cc=poly.c # (4,ninterval,nxi) descending powers
    mi=np.inf;ma=-np.inf;max_deriv=0.
    for seg in range(cc.shape[1]):
        dx=float(knots[seg+1]-knots[seg])
        for j in range(cc.shape[2]):
            a,b,c,d=[float(v) for v in cc[:,seg,j]]
            critical=[0.,dx]
            if abs(a)>1e-26:
                discr=4*b*b-12*a*c
                if discr>=0:
                    rt=np.sqrt(discr)
                    critical.extend([t for t in ((-2*b-rt)/(6*a),(-2*b+rt)/(6*a)) if 0<t<dx])
            elif abs(b)>1e-26:
                t=-c/(2*b)
                if 0<t<dx:critical.append(t)
            vals=[((a*t+b)*t+c)*t+d for t in critical]
            mi=min(mi,*vals);ma=max(ma,*vals)
            deriv_checks=[0.,dx]
            if abs(a)>1e-26:
                t=-b/(3*a)
                if 0<t<dx:deriv_checks.append(t)
            max_deriv=max(max_deriv,*[abs((3*a*t+2*b)*t+c) for t in deriv_checks])
    return float(mi),float(ma),float(max_deriv)


def derivative_diagnostics(per, safety_factor=64.,model_error_allowance=1e-7):
    if not (1<=safety_factor<=1e6) or model_error_allowance<0:raise ValueError('Invalid margins')
    out={};pieces=[]
    for seed,pp in per.items():
        ks=pp['k'];xi=pp['xi'];m=pp['data']['margin']
        dk=np.diff(ks);dt=np.diff(xi)
        sec_k=np.diff(m,axis=0)/dk[:,None]
        sec_xi=np.diff(m,axis=1)/dt[None,:]
        cubic=PchipInterpolator(ks,m,axis=0,extrapolate=False)
        min_poly,max_poly,Lk_pchip=extrema_pchip(cubic,ks)
        Lk_sec=float(np.max(abs(sec_k)));Lxi_sec=float(np.max(abs(sec_xi)))
        Lk_obs=max(Lk_sec,Lk_pchip)
        Lxi_obs=Lxi_sec
        center_min=float(m.min());gap_k=float(np.max(dk));gap_xi=float(np.max(dt))
        penalty_k=Lk_obs*gap_k*.5
        penalty_xi=Lxi_obs*gap_xi*.5
        lo=center_min-safety_factor*(penalty_k+penalty_xi)-model_error_allowance
        critical_factor=(center_min-model_error_allowance)/(penalty_k+penalty_xi)
        diag={'kappa_count':len(ks),'xi_count':len(xi),'kappa_interval':[float(ks[0]),float(ks[-1])],
              'xi_grid_interval':[float(xi[0]),float(xi[-1])],
              'min_sampled_margin':center_min,'max_sampled_margin':float(m.max()),
              'min_pchip_parameter_surrogate':min_poly,'max_pchip_parameter_surrogate':max_poly,
              'max_abs_kappa_secant':Lk_sec,'max_abs_parameter_pchip_derivative':Lk_pchip,
              'max_abs_xi_secant':Lxi_sec,
              'max_kappa_spacing':gap_k,'max_xi_spacing':gap_xi,
              'assumed_derivative_inflation_factor':safety_factor,
              'assumed_parameter_lipschitz':safety_factor*Lk_obs,
              'assumed_time_lipschitz':safety_factor*Lxi_obs,
              'assumed_absolute_numerical_allowance':model_error_allowance,
              'conditional_2d_lower_bound':lo,
              'empirical_headroom_factor':float(critical_factor),
              'is_2d_conditional_bound_positive':bool(lo>0),
              'derivative_bound_is_rigorous':False,
              'continuous_physical_positivity_certified':False}
        out[seed]=diag
        for j,(lo_k,hi_k) in enumerate(zip(ks[:-1],ks[1:])):
            pr=cubic(.5*(lo_k+hi_k))
            pieces.append({'seed':seed,'gap_index':j,'kappa_lo':lo_k,'kappa_hi':hi_k,
                           'kappa_mid':.5*(lo_k+hi_k),'gap_width':hi_k-lo_k,
                           'max_abs_kappa_secant':float(np.max(abs(sec_k[j]))),
                           'min_endpoint_margin':float(min(m[j].min(),m[j+1].min())),
                           'min_surrogate_midpoint':float(np.min(pr)),
                           'surrogate_midpoint_at_endpoint_time':float(pr[-1])})
        pp['pchip']=cubic
        pp['dsec_k']=sec_k
    return out,pieces


def module_from(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj
    spec.loader.exec_module(obj)
    return obj


def safe_extract_parent(path,area):
    with zipfile.ZipFile(path) as z:
        for m in z.infolist():
            if m.is_dir():continue
            p=Path(m.filename)
            if p.is_absolute() or '..' in p.parts or p.parts[0]!=BROOT:
                raise RuntimeError('Unsafe authenticated B archive path')
            dest=area/p;dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(z.read(m.filename))
    return area/BROOT


def probe_physical_midpoints(per,pieces,source_b_path,source_s_path,mode):
    """Run true nonlinear+forced tangent at new kappa values. Never interpolate the model.

    Requires exact historical 51s genealogy to stage the physical generator.
    """
    if mode=='none':return [],[],{}
    s_expected=source_manifest()['parents']['v18_51s_sha256']
    if sha_file(source_s_path)!=s_expected:raise RuntimeError('v18.51s source changed')
    selection=[r for r in pieces if mode=='all' or r['gap_index'] in (0,11)]
    found=[];failures=[];probe_curves={}
    with tempfile.TemporaryDirectory(prefix='q550_v1852c_phys_') as td:
        stage=Path(td)
        b_root=safe_extract_parent(source_b_path,stage)
        b=module_from(b_root/'scripts/run_audit.py','q550_v1852c_true_52b')
        a_root=b.stage_parent(b_root/'inputs'/b.ANAME,stage/'ancestor',source_s_path)
        a=module_from(a_root/'scripts/run_audit.py','q550_v1852c_true_52a')
        t=module_from(a_root/'reference/parent_v18_51t_run_audit.py','q550_v1852c_true_51t')
        print('Physical model genealogy authenticated. Staging original 10D generator...',flush=True)
        _,_,rblob,srun,_,_=a.load_parents()
        q,p,v,vr,parts,kappas,lineage,_=t.stage_engine(stage/'engine',rblob,srun)
        _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
        if len(lineage)!=19:raise RuntimeError('Unexpected physical generator ancestry')
        if abs(kappas['midpoint']-float(np.mean(per['original']['k'][5:8])))>1e-11:
            raise RuntimeError('Kappa historical midpoint drift')
        mx=float(engine.params.m_x)
        for seed in SEEDS:
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            for r in selection:
                if r['seed']!=seed:continue
                kval=float(r['kappa_mid'])
                # B inspect_real reads the actual parent 52a grid from its arguments.
                # We never reuse an interpolation as a physical trajectory.
                physical_xi=np.linspace(0.,.2,513)
                series={}
                label='midpoint_gap'+str(r['gap_index'])
                print(f'Q550 v18.52c physical probe {seed} gap {r["gap_index"]}/11 kappa={kval:.12f}',flush=True)
                try:
                    rr=b.inspect_real(seed,label,kval,ctx,physical_xi,a,t,q,p,v,engine,mx,
                                      float(kappas['midpoint']),float((per[seed]['k'][-1]-per[seed]['k'][0])/64),
                                      series,refine=False)
                    curve=series[seed+':'+label]['margin']
                    mask=physical_xi>=.01
                    baseline=per[seed];src=baseline['data']['margin']
                    j=r['gap_index'];expected=baseline['pchip'](kval)
                    actual=curve[mask]
                    max_defect=float(np.max(abs(actual-expected)))
                    mid_sec0=(actual-src[j])/(kval-r['kappa_lo'])
                    mid_sec1=(src[j+1]-actual)/(r['kappa_hi']-kval)
                    data={**r,'minimum_physical_midpoint_margin':float(np.min(actual)),
                          'max_abs_pchip_midpoint_defect':max_defect,
                          'max_abs_midpoint_left_slope':float(np.max(abs(mid_sec0))),
                          'max_abs_midpoint_right_slope':float(np.max(abs(mid_sec1))),
                          'model_solver_nfev':rr['solver_nfev'],
                          'tangent_nfev':rr['tangent_nfev'],
                          'max_recomposition':rr['max_recomposition'],
                          'max_budget_residual':rr['max_budget_residual'],
                          'sampled_physical_status':rr['status']}
                    found.append(data)
                    probe_curves[f'{seed}__gap{j}__margin']=curve
                    if rr['status'] not in ('SAMPLED_PASS','MARGIN_SIGN_LOSS',
                                            'TOTAL_TURNING_SIGN_LOSS','ORIENTED_AREA_SIGN_LOSS',
                                            'SELF_CHANNEL_SIGN_LOSS'):
                        failures.append({'seed':seed,'gap':j,'kappa':kval,
                                         'category':'NUMERIC_OR_NEAR_ZERO_UNRESOLVED','details':rr['status']})
                except Exception as exc:
                    failures.append({'seed':seed,'gap':r['gap_index'],'kappa':kval,
                                     'category':'MODEL_OR_SOLVER_UNRESOLVED','details':str(exc)})
                    print('Q550 v18.52c unresolved probe:',str(exc),flush=True)
    return found,failures,probe_curves


def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='',encoding='utf-8') as stream:
        names=list(rows[0]); w=csv.DictWriter(stream,fieldnames=names,extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def plot(per,diag,probes):
    fig,axs=plt.subplots(2,2,figsize=(13,8),constrained_layout=True)
    for seed,color in [('original','tab:blue'),('autonomous','tab:orange')]:
        pp=per[seed];k=pp['k'];xi=pp['xi'];m=pp['data']['margin']
        axs[0,0].plot(k,np.min(m,axis=1),'.-',color=color,label=seed)
        axs[0,1].plot(xi,np.max(np.abs(pp['dsec_k']),axis=0),color=color,label=seed)
        q=diag[seed]
        axs[1,0].bar(seed,q['conditional_2d_lower_bound'],alpha=.6)
        axs[1,1].bar(seed,q['empirical_headroom_factor'],alpha=.6)
    axs[0,0].set_title('Minimum sampled signed margin vs kappa')
    axs[0,0].set_xlabel('Physical kappa0');axs[0,0].set_ylabel('Signed margin (rad / xi)')
    axs[0,1].set_title('Observed max |parameter secant| by xi')
    axs[0,1].set_xlabel('xi');axs[0,1].set_ylabel('|delta margin / delta kappa0|')
    axs[1,0].set_title('Conditional 2D lower bound (derivative assumption)')
    axs[1,1].set_title('Factor before conditional margin reaches zero')
    for ax in axs.ravel():ax.grid(alpha=.25);ax.axhline(0,color='gray',lw=.8)
    axs[0,0].legend();axs[0,1].legend()
    fig.suptitle('Q550 v18.52c — Measured response vs CONDITIONAL continuum reasoning',fontsize=13)
    fig.savefig(OUT/'signed_margin_parameter_derivative_conditional_bound.png',dpi=155)
    plt.close(fig)
    if probes:
        fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
        for seed,color in [('original','tab:blue'),('autonomous','tab:orange')]:
            rr=[p for p in probes if p['seed']==seed]
            if rr:ax.plot([r['kappa_mid'] for r in rr],
                          [r['max_abs_pchip_midpoint_defect'] for r in rr],'.-',color=color,label=seed)
        ax.set_title('Independent nonlinear midpoint vs parameter PCHIP interpolant')
        ax.set_xlabel('kappa0');ax.set_ylabel('max absolute midpoint interpolation discrepancy')
        ax.grid(alpha=.3);ax.legend()
        fig.savefig(OUT/'physical_midpoint_interpolation_defects.png',dpi=155);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=TITLE)
    parser.add_argument('--preflight-only',action='store_true',help='Analyze genuine B outputs without new physical model IVPs')
    parser.add_argument('--probes',choices=('all','outer','none'),default='all',
                        help='Run 24, 4, or 0 genuine nonlinear+tangent interval-midpoint probes')
    parser.add_argument('--safety-factor',type=float,default=64.,
                        help='Hypothetical upper bound multiplier; NOT rigorously established')
    parser.add_argument('--absolute-error-allowance',type=float,default=1e-7,
                        help='Assumed uniform numeric allowance for CONDITIONAL inequality')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    b_path,b_report,rows,arrays,hashes=verify_source_manifest()
    per=load_surface(rows,arrays)
    diag,pieces=derivative_diagnostics(per,args.safety_factor,args.absolute_error_allowance)
    first={'version':VERSION,'mode':'preflight','parent_version':'v18.52b',
           'parent_sha256':sha_file(b_path),'parent_output_sha256_direct':hashes,
           'diagnostics':diag,'derivative_bounds_independently_verified':False,
           'physical_continuous_certificate':False}
    (OUT/'LOCAL_PREFLIGHT.json').write_text(json.dumps(first,indent=2)+'\n')
    write_csv(OUT/'parameter_gap_derivative_ledger.csv',pieces)
    plot(per,diag,[])
    print('Q550 v18.52c: authenticated 26 physical branches; surrogate and conditional bounds computed.',flush=True)
    for seed in SEEDS:
        d=diag[seed]
        print(f'  {seed}: min margin={d["min_sampled_margin"]:.9f}, '
              f'max secant kappa={d["max_abs_kappa_secant"]:.6f}, '
              f'conditional lower={d["conditional_2d_lower_bound"]:.6f}, '
              f'critical inflation={d["empirical_headroom_factor"]:.1f}',flush=True)
    if args.preflight_only:
        print('PREFLIGHT ONLY — no new physical midpoint branches executed',flush=True)
        return
    mode=args.probes
    s_path=find_input(SNAME) if mode!='none' else None
    probes,failures,curves=probe_physical_midpoints(per,pieces,b_path,s_path,mode)
    write_csv(OUT/'physical_midpoint_probe_ledger.csv',probes)
    write_csv(OUT/'physical_probe_failures.csv',failures)
    if curves:np.savez_compressed(OUT/'physical_midpoint_margin_curves.npz',**curves)
    plot(per,diag,probes)
    expected_count={'all':24,'outer':4,'none':0}[mode]
    for seed in SEEDS:
        related=[p for p in probes if p['seed']==seed]
        d=diag[seed]
        d['physical_midpoint_probes_completed']=len(related)
        d['minimum_physical_midpoint_margin']=min((p['minimum_physical_midpoint_margin'] for p in related),default=None)
        d['maximum_physical_midpoint_pchip_defect']=max((p['max_abs_pchip_midpoint_defect'] for p in related),default=None)
        d['maximum_observed_midpoint_subsecant']=max((max(p['max_abs_midpoint_left_slope'],
                p['max_abs_midpoint_right_slope']) for p in related),default=None)
        d['hypothetical_Lk_exceeded_by_midpoint_secants']=(
            d['maximum_observed_midpoint_subsecant'] is not None and
            d['maximum_observed_midpoint_subsecant']>d['assumed_parameter_lipschitz'])
    checks={
        'authenticated_52b_26_physical_branches':len(rows)==26,
        'historical_52b_computation_passed':b_report['all_passed'],
        'valid_two_seed_parameter_grid':all(len(per[s]['k'])==13 for s in SEEDS),
        'parameter_interpolant_positive_on_qualified_grid':all(d['min_pchip_parameter_surrogate']>0 for d in diag.values()),
        'sampled_margins_all_positive':all(d['min_sampled_margin']>0 for d in diag.values()),
        'interpolant_has_finite_derivative':all(np.isfinite(d['max_abs_parameter_pchip_derivative']) for d in diag.values()),
        'derivative_assumption_is_explicit_not_certified':all(not d['derivative_bound_is_rigorous'] for d in diag.values()),
        'hypothetical_2d_lower_positive':all(d['is_2d_conditional_bound_positive'] for d in diag.values()),
        'new_midpoint_count_completed':len(probes)==expected_count,
        'new_midpoint_recompositions_small':all(float(p['max_recomposition'])<1e-8 for p in probes),
        'new_midpoint_angular_budget_small':all(float(p['max_budget_residual'])<1e-7 for p in probes),
        'hypothetical_lipschitz_not_falsified_by_new_secants':all(
            not d['hypothetical_Lk_exceeded_by_midpoint_secants'] for d in diag.values()),
    }
    if mode=='all' and not failures and len(probes)==24:
        science='PHYSICAL_MIDPOINTS_POSITIVE__CONDITIONAL_CONTINUUM_BOUND_ONLY'
    elif failures:
        science='PHYSICAL_MIDPOINT_PROBES_PARTIAL_OR_UNRESOLVED'
    else:
        science='PHYSICAL_MIDPOINT_COVERAGE_INCOMPLETE__CONDITIONAL_BOUND_ONLY'
    if any(p['minimum_physical_midpoint_margin']<=0 for p in probes):
        science='SAMPLED_PHYSICAL_SIGN_LOSS_FOUND'
    if any(d['hypothetical_Lk_exceeded_by_midpoint_secants'] for d in diag.values()):
        science='ASSUMED_DERIVATIVE_BOUND_REFUTED_BY_MIDPOINTS'
    report={'version':VERSION,'title':TITLE,'mode':'server_physical_midpoint_probe',
            'all_numerical_checks_passed':bool(all(checks.values()) and not failures),
            'numerical_checks':checks,
            'scientific_decision':science,
            'classification':'CONDITIONAL_DERIVATIVE_BOUND_NUMERICAL_DIAGNOSTIC_ONLY',
            'historical_52b_branches':26,'new_physical_midpoints_completed':len(probes),
            'requested_new_physical_midpoints':expected_count,
            'new_physical_midpoint_failures':failures,
            'differentiation':diag,
            'max_physical_vs_surrogate_midpoint_defect':float(max((p['max_abs_pchip_midpoint_defect'] for p in probes),default=float('nan'))) if probes else None,
            'physical_continuous_certificate':False,
            'sufficiency_condition':'If abs(dm/dkappa)<=Lk, abs(dm/dxi)<=Lxi everywhere in the rectangle, '
                                    'and uniform absolute numerical uncertainty <=epsilon, '
                                    'then m >= min_grid_margin - Lk*max_gap_kappa/2 - Lxi*max_gap_xi/2 - epsilon.',
            'limitations':['No rigorous enclosure of the exact model parameter derivative between kappa samples.',
              'No rigorous enclosure of the exact model time derivative between xi samples.',
              'A safety-factor inflated secant is a hypothesis, not a mathematical derivative bound.',
              'PCHIP positivity is a property of an interpolant, not of the original nonlinear dynamics.',
              'New physical midpoint probes cannot exclude narrow sub-grid negative pockets.',
              'The signed angular margin depends on the original model channel attribution.',
              'Parameter interval is limited to the tested 52b band and xi to the qualified grid.']}
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# Q550 v18.52c — Server numerical report','',
           f'- Scientific classification: **{science}**',
           f'- New genuine nonlinear/tangent midpoint branches: **{len(probes)}/{expected_count}**',
           f'- Numerical checks: **{sum(checks.values())}/{len(checks)}**','']
    for seed in SEEDS:
        d=diag[seed]
        lines.extend([f'## {seed}',f'- Minimum sampled margin: {d["min_sampled_margin"]:.12g}',
            f'- Max observed |kappa secant|: {d["max_abs_kappa_secant"]:.12g}',
            f'- Max observed |xi secant|: {d["max_abs_xi_secant"]:.12g}',
            f'- Hypothetical conditional 2D lower bound: {d["conditional_2d_lower_bound"]:.12g}',
            f'- Inflation factor at conditional zero: {d["empirical_headroom_factor"]:.8g}',''])
    lines.extend(['## Epistemic boundary','','**No unconditional continuous-robustness certificate has been established.**',
                  'Physical midpoint probes are independent nonlinear/tangent evaluations but remain a finite parameter sample.',
                  'A true continuous certificate requires validated derivative enclosures over the kappa–xi rectangle.',''])
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines))
    print('Q550 v18.52c:',science,'checks',sum(checks.values()),'/',len(checks),flush=True)
    if not report['all_numerical_checks_passed']:raise SystemExit(3)


if __name__=='__main__':main()
