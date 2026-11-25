#!/usr/bin/env python3
"""Q550 v18.53t: genuinely coupled on-the-fly 11D+10D accumulator.

Only the default full run performs physical IVPs. Preflight authenticates the complete
historical numerical genealogy. Scientific documentation stays in English.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,tempfile,traceback,zipfile
from contextlib import ExitStack
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
from physics_21d import rhs_21d,source_array,source_10d,angular_polygon

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'inputs_config.json').read_text())
OUT=ROOT/'results'; PARENT=ROOT/'inputs'/CFG['parent_filename']
SHA=lambda b:hashlib.sha256(b).hexdigest()


def add_module(path,name):
    path=Path(path)
    if str(path.parent) not in sys.path:sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod)
    return mod


def authenticate_zip(z):
    assert z.testzip() is None,'S_PARENT_CRC_FAILED'
    names=z.namelist()
    assert names and not any(Path(p).is_absolute() or '..' in Path(p).parts for p in names),'S_PARENT_PATHS_UNSAFE'
    roots={p.split('/')[0] for p in names if '/' in p}
    assert len(roots)==1,'S_PARENT_ROOT_AMBIGUOUS'
    prefix=next(iter(roots))+'/'
    nums={}
    for key in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
        files=json.loads(z.read(prefix+key))['files']
        for rel,digest in files.items():
            if SHA(z.read(prefix+rel))!=digest:raise RuntimeError('PARENT_HASH_MISMATCH '+rel)
        nums[key]=len(files)
    report=json.loads(z.read(prefix+'results/report.json'))
    if report['scientific_decision']!=CFG['parent_decision'] or report['new_physical_ivps']!=6:
        raise RuntimeError('S_PARENT_NOT_AUTHENTIC_PHYSICAL_CHECKPOINT')
    return prefix,nums,report


def load_ancestors(stack):
    if not PARENT.is_file() or SHA(PARENT.read_bytes())!=CFG['parent_sha256']:
        raise RuntimeError('S_PARENT_SHA_MISMATCH')
    with zipfile.ZipFile(PARENT) as z:
        prefix,counts,report=authenticate_zip(z)
        folder=Path(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853t_s_')))
        z.extractall(folder)
    sroot=folder/prefix.rstrip('/')
    s=add_module(sroot/'scripts/run_audit.py','q550_v1853t_s')
    r,rr,old,centers,pins,rcounts,rreport,archive=s.authenticated_parent(stack)
    o,n,m,l,k,j,h,ncount,mcount=r.require_ancestry(archive,stack)
    from collections import Counter
    rep=dict(s=counts,r=rcounts,n=ncount,m=mcount,ancestors='s->r->q->p->o->n->m->l->k->j->h')
    return s,r,h,m,k,o,j,centers,sroot,rr,rep


def csvrows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))


def csvwrite(name,rows):
    if not rows:raise RuntimeError('EMPTY_CSV '+name)
    keys=list(dict.fromkeys(key for row in rows for key in row))
    with (OUT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def maxrow(a):
    return float(np.max(np.linalg.norm(np.asarray(a,float),axis=-1)))


def solve_ivp_dense(fun,t0,t1,y0,rtol,atol,label):
    out=solve_ivp(fun,(t0,t1),y0,method='DOP853',dense_output=True,
                  rtol=rtol,atol=atol,max_step=CFG['max_step_xi'])
    if not out.success or out.sol is None or not np.all(np.isfinite(out.y)):
        raise RuntimeError(label+'_IVP_FAILED '+str(out.message))
    return out


def accepted_step_diag(sol):
    dt=np.diff(sol.t)
    return dict(nfev=int(sol.nfev),accepted_steps=int(len(dt)),
                median_accepted_step=float(np.median(dt)),
                max_accepted_step=float(np.max(dt)),
                fraction_max_step_limited=float(np.mean(dt>=.999*CFG['max_step_xi'])))


def symmetric_accepted_mesh_gap(a,b):
    """Maximum nearest-neighbour xi distance between the accepted-step meshes."""
    def directed(u,v):
        ix=np.searchsorted(v,u)
        left=np.abs(u-v[np.clip(ix-1,0,len(v)-1)])
        right=np.abs(u-v[np.clip(ix,0,len(v)-1)])
        return float(np.max(np.minimum(left,right)))
    return max(directed(np.asarray(a,float),np.asarray(b,float)),
               directed(np.asarray(b,float),np.asarray(a,float)))


def geometry_case(r,seed,rule,variant,xi,C,g,center,stored,records):
    # C is integrated in xi-clock; its additive integration constant cancels.
    proxy=r.Primitive(xi,C,g,'hermite')
    oc,wc=center
    origins=np.linspace(oc-CFG['origin_halfspan'],oc+CFG['origin_halfspan'],CFG['curve_points'])
    wlo,wup=wc-CFG['width_halfband'],wc+CFG['width_halfband']
    if origins[0]<xi[0] or origins[-1]+3*wup>xi[-1]:raise RuntimeError('GEOMETRY_OUTSIDE_TRACE')
    mins=[r.optimized_width(proxy,float(a),wlo,wup,nsearch=CFG['search_brackets']) for a in origins]
    ww=np.array([v['w'] for v in mins]);dw=np.array([v['implicit_wprime'] for v in mins])
    geom=r.curve(proxy,origins,ww,dw)
    if np.min(geom['norm'])<CFG['min_norm_screen']:raise RuntimeError('NEAR_ZERO_UNDEFINED_DIRECTION')
    row=dict(seed=seed,kappa_rule=rule,variant=variant,origin_center=float(oc),
             width_center=float(ww[len(ww)//2]),norm_center=float(geom['norm'][len(ww)//2]),
             min_residual_10d=float(geom['min_norm']),
             direction_speed_center=float(geom['directional_speed'][len(ww)//2]),
             arclength_direct=float(geom['direct']),
             arclength_tangent=float(geom['chain']),
             oriented_polygon_norm=float(np.linalg.norm(geom['oriented_polygon'])),
             width_sign_changes=int(geom['width_sign_changes']),
             width_total_variation=float(geom['width_total_variation']),
             bound_count=int(sum(v['bound'] for v in mins)),
             all_far_endpoints_formal=bool(np.all(origins+3*ww>CFG['formal_clock_boundary'])))
    records.append(row)
    key=seed+'__'+rule+'__'+variant
    stored[key+'__origins']=origins
    stored[key+'__widths']=ww
    stored[key+'__u10d']=geom['u']
    stored[key+'__polygon45d']=geom['oriented_polygon']
    stored[key+'__direction_speed']=geom['directional_speed']
    return origins,ww,geom,row


def physical_cases(args,s,r,h,m,k,o,j,centers,sroot,rr):
    parent_report=csvrows(sroot/'results/physical_step_reference.csv')
    expected={(a['seed'],a['kappa_rule']):a for a in parent_report}
    with np.load(sroot/'results/oriented_valley_quadrature_comparison.npz',allow_pickle=False) as data:
        archived={name:np.array(data[name],copy=True) for name in data.files
                  if name.endswith(('__tight_xi','__tight_g10d','__tight_z11d'))}
    loader=h.loader()
    ancestor=loader.find_parent(args.parent_v1851s)
    native=args.native_npz or loader.find_npz('10d',args.base)
    if not ancestor or SHA(Path(ancestor).read_bytes())!=o.CFG['physical_parent_source_sha256']:
        raise RuntimeError('HISTORICAL_V1851S_SOURCE_NOT_FOUND_OR_HASH_WRONG')
    if not native or SHA(Path(native).read_bytes())!=o.CFG['physical_native_source_sha256']:
        raise RuntimeError('HISTORICAL_10D_NATIVE_NOT_FOUND_OR_HASH_WRONG')
    stage,engine,contexts,meta=loader.stage_original(ancestor)
    rows=[];geo=[];store={};warnings=[]
    try:
        with np.load(native,allow_pickle=False) as f:
            for seed in CFG['seeds']:
                ctx=contexts[seed]
                for rule in CFG['kappa_rules']:
                    tag=seed+'__'+rule
                    xi=archived[tag+'__tight_xi']
                    z_s=archived[tag+'__tight_z11d']
                    g_s=archived[tag+'__tight_g10d']
                    if len(xi)!=8961:raise RuntimeError('UNEXPECTED_4X_GRID_LENGTH')
                    if not np.all(np.diff(xi)>0):raise RuntimeError('GRID_NOT_STRICTLY_MONOTONE')
                    xbase,nh=j.xi_grid(f['common_xi'])
                    if not np.allclose(xi[::4],xbase,rtol=0,atol=2e-12):
                        raise RuntimeError('CLOCK_BASE_GRID_MISMATCH')
                    _,_,hist=h.native_53g_reconciliation(f,xbase[:nh],seed,ctx,engine.omega)
                    z0=np.asarray(hist[0],float)
                    kp=float(k.K[rule])
                    rtol=CFG['rtol'];atol=CFG['atol_physical'];max_step=CFG['max_step_xi']
                    sol11=solve_ivp_dense(lambda t,z:h.rhs(engine,ctx,kp,t,z),
                                          float(xi[0]),float(xi[-1]),z0,rtol,atol,'11D_'+tag)
                    z11=sol11.sol(xi).T
                    g11=source_array(h,m,ctx,kp,xi,z11,engine.omega)
                    # Validate vectorized source against direct source formula at original xi nodes.
                    scalar_probe=[]
                    for i in np.linspace(0,len(xi)-1,9,dtype=int):
                        src=source_10d(h,m,ctx,kp,float(xi[i]),z11[i],
                                      tau_rate=float(engine.omega)/float(expit(z11[i,8])))
                        scalar_probe.append(np.linalg.norm(src-g11[i]))
                    source_scalar_error=float(max(scalar_probe))
                    # Validate original physical invariants separately from accumulator.
                    rhs_error_11=h.validated_rhs(engine,ctx,kp,xi,z11)
                    # Integrated 10D on physical trajectory frozen to the 11D dense-output path.
                    # The 10D solution does not feed back on the physical integrator.
                    sol10=solve_ivp_dense(
                        lambda t,c:source_10d(h,m,ctx,kp,t,sol11.sol(t),
                          tau_rate=float(engine.omega)/float(expit(sol11.sol(t)[8]))),
                        float(xi[0]),float(xi[-1]),np.zeros(10),rtol,
                        CFG['atol_accumulator'],'10D_LOCKED_'+tag)
                    C_locked=sol10.sol(xi).T
                    # Fully coupled numerical integration; physical RHS remains 11D-autonomous.
                    yinit=np.r_[z0,np.zeros(10)]
                    atvec=np.r_[np.full(11,atol),np.full(10,CFG['atol_accumulator'])]
                    sol21=solve_ivp_dense(
                        lambda t,y:rhs_21d(h,m,engine,ctx,kp,t,y),
                        float(xi[0]),float(xi[-1]),yinit,rtol,atvec,'21D_'+tag)
                    state21=sol21.sol(xi).T
                    z21=state21[:,:11];C21=state21[:,11:]
                    g21=source_array(h,m,ctx,kp,xi,z21,engine.omega)
                    rhs_error_21=h.validated_rhs(engine,ctx,kp,xi,z21)
                    _,_,charge,_=h.embed(z21[:,:10],h.P['m_x'])
                    charge_range=float(np.ptp(charge))
                    if charge_range>float(k.CFG['three_window']['charge_abs_drift_tolerance']):
                        raise RuntimeError('AUGMENTED21_CHARGE_DRIFT '+tag+' '+str(charge_range))
                    # Exactly the sampled source as in v18.53s: integration is post hoc.
                    sourceC=s.source_pchip_integral(xi,g_s)
                    trapezoidC=s.cumulative_trapezoid(xi,g11)
                    direct_compare=maxrow(g11-g_s)
                    coupling_state=maxrow(z21-z11)
                    coupling_primitive=maxrow(C21-C_locked)
                    lock_vs_sample=maxrow(C_locked-sourceC)
                    row=dict(seed=seed,kappa_rule=rule,reference_s_nfev=int(expected[seed,rule]['nfev']),
                             same_evaluation_nodes=len(xi),max_step_xi=max_step,rtol=rtol,atol11=atol,
                             max_state_11d_vs_parent_s=maxrow(z11-z_s),
                             source_scalar_vectorized_probe_error=source_scalar_error,
                             validated_rhs_11d_rel_residual=float(rhs_error_11),
                             validated_rhs_21d_rel_residual=float(rhs_error_21),
                             augmented21_charge_range=charge_range,
                             accepted_mesh_symmetric_gap=symmetric_accepted_mesh_gap(sol11.t,sol21.t),
                             max_source_10d_vs_parent_s=direct_compare,
                             max_11d_change_when_accumulator_included=coupling_state,
                             max_source_change_when_accumulator_included=maxrow(g21-g11),
                             max_integrated_primitive_coupled_vs_locked=coupling_primitive,
                             max_locked_integral_vs_s_source_pchip=lock_vs_sample,
                             max_coupled_integral_vs_s_source_pchip=maxrow(C21-sourceC),
                             max_locked_integral_vs_trapezoid=maxrow(C_locked-trapezoidC),
                             max_s_integral_vs_trapezoid=maxrow(sourceC-trapezoidC),
                             **{'physical11_'+kk:vv for kk,vv in accepted_step_diag(sol11).items()},
                             **{'locked10_'+kk:vv for kk,vv in accepted_step_diag(sol10).items()},
                             **{'augmented21_'+kk:vv for kk,vv in accepted_step_diag(sol21).items()})
                    rows.append(row)
                    for variant,Cval,gval in (
                        ('same_grid_s_source_PCHIP',sourceC,g_s),
                        ('standalone11_locked_integral',C_locked,g11),
                        ('on_the_fly_coupled21',C21,g21),
                        ('standalone11_trapezoid',trapezoidC,g11),
                    ):
                        result=geometry_case(r,seed,rule,variant,xi,Cval,gval,
                                             centers[seed,rule],store,geo)
                        if variant=='on_the_fly_coupled21':
                            origins,widths,geom,rec=result
                            # Direct dense-output C at optimised endpoints: independent
                            # check of cubic Hermite interpolation between stored knots.
                            ep=origins+3*widths
                            Nexact=(sol21.sol(ep)[11:,:]-sol21.sol(origins)[11:,:]).T
                            rrnorm=np.linalg.norm(Nexact,axis=1)
                            u=Nexact/rrnorm[:,None]
                            from source_valley_math import wedge
                            direct,poly=angular_polygon(u,wedge)
                            proxy=r.Primitive(xi,C21,g21,'hermite')
                            Nh=(proxy.eval(ep)-proxy.eval(origins))
                            rec['dense_evaluation_max_residual_interp_error']=maxrow(Nexact-Nh)
                            rec['dense_evaluation_arclength']=float(direct)
                            rec['dense_evaluation_arclength_minus_Hermite']=float(direct-geom['direct'])
                            rec['dense_evaluation_oriented45d_angle_deg']=float(np.degrees(np.arccos(
                                np.clip(poly@geom['oriented_polygon']/
                                max(np.linalg.norm(poly)*np.linalg.norm(geom['oriented_polygon']),1e-99),-1,1))))
                            rec['dense_evaluation_min_residual_10d']=float(np.min(rrnorm))
                            store[tag+'__dense_u10d']=u
                            store[tag+'__dense_oriented45d']=poly
                    # Preserve normalized physical and integral checks, not full dense sol.
                    store[tag+'__xi']=xi
                    store[tag+'__C21']=C21
                    store[tag+'__Clocked']=C_locked
                    store[tag+'__CsourcePCHIP']=sourceC
                    store[tag+'__g21']=g21
                    store[tag+'__z21']=z21
                    print('V1853T_PHYSICAL_CASE',seed,rule,'SOLVES',len(rows)*3,
                          '11D/21D_NFEV',sol11.nfev,sol21.nfev,flush=True)
                    if coupling_state>CFG['max_physical_drift_warning']:
                        warnings.append(tag+': augmentation changed 11D beyond warning threshold')
                    if direct_compare>CFG['max_source_reconstruction_warning']:
                        warnings.append(tag+': reproduction of s source above warning threshold')
                    if coupling_primitive>CFG['max_accumulator_discrepancy_warning']:
                        warnings.append(tag+': locked/coupled accumulator mismatch above warning threshold')
    finally:
        stage.cleanup()
    return rows,geo,store,warnings


def make_plot(geo):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for seed in CFG['seeds']:
            rows=[q for q in geo if q['seed']==seed and q['kappa_rule']=='original']
            fig,axs=plt.subplots(1,2,figsize=(13,4.5))
            labels=[q['variant'].replace('standalone11_','').replace('same_grid_s_','s:') for q in rows]
            axs[0].bar(labels,[q['min_residual_10d'] for q in rows]);axs[0].set_ylabel('Minimum 10D residual')
            axs[1].bar(labels,[q['arclength_direct'] for q in rows]);axs[1].set_ylabel('Direct angular length (rad)')
            for a in axs:a.tick_params(axis='x',rotation=55);a.grid(alpha=.2,axis='y')
            fig.suptitle('Q550 v18.53t on-the-fly versus sampled integrals — '+seed)
            fig.tight_layout()
            fig.savefig(OUT/f'augmented21_valley_{seed}.png',dpi=140)
            plt.close(fig)
    except Exception as e:
        (OUT/'PLOT_WARNING.txt').write_text(str(e)+'\n')


def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument('--preflight',action='store_true',help='Verify genealogical sources, no IVPs')
    ap.add_argument('--base',type=Path,default=ROOT.parent)
    ap.add_argument('--parent-v1851s',type=Path)
    ap.add_argument('--native-npz',type=Path)
    ap.add_argument('--require-parent',action='store_true')
    ap.add_argument('--require-native',action='store_true')
    ap.add_argument('--require-dynamic',action='store_true')
    args=ap.parse_args(argv)
    OUT.mkdir(exist_ok=True)
    report=dict(version='v18.53t',scientific_decision='NOT_EVALUATED',
                completed_physical11_ivps=0,completed_coupled21_ivps=0,
                completed_locked10_ivps=0)
    try:
        if args.preflight and (args.require_parent or args.require_native or args.require_dynamic):
            raise RuntimeError('PREFLIGHT_IS_NONPHYSICAL')
        with ExitStack() as stack:
            s,r,h,m,k,o,j,centers,sroot,rr,ancestry=load_ancestors(stack)
            report['authenticated_ancestry']=ancestry
            if args.preflight:
                print('V1853T_PREFLIGHT_PASS',ancestry,flush=True)
                return 0
            physical,geo,store,warnings=physical_cases(args,s,r,h,m,k,o,j,centers,sroot,rr)
        report['completed_physical11_ivps']=len(physical)
        report['completed_coupled21_ivps']=len(physical)
        report['completed_locked10_ivps']=len(physical)
        if len(physical)!=6 or len(geo)!=24:raise RuntimeError('INCOMPLETE_PHYSICAL_AND_GEOMETRIC_CASES')
        csvwrite('physical_21d_vs_11d.csv',physical)
        csvwrite('valley_geometry_21d.csv',geo)
        np.savez_compressed(OUT/'augmented21_primitives_and_oriented_valleys.npz',**store)
        make_plot(geo)
        report.update(
            sample_grid_nodes=8961,
            bivector_components=10,
            max_physical_state_drift=float(max(q['max_11d_change_when_accumulator_included'] for q in physical)),
            max_coupled_locked_accumulator_difference=float(max(q['max_integrated_primitive_coupled_vs_locked'] for q in physical)),
            max_coupled_vs_s_pchip_difference=float(max(q['max_coupled_integral_vs_s_source_pchip'] for q in physical)),
            max_dense_vs_hermite_residual_difference=float(max(q['dense_evaluation_max_residual_interp_error'] for q in geo if q['variant']=='on_the_fly_coupled21')),
            max_direct_arclength_delta_dense_vs_hermite=float(max(abs(q['dense_evaluation_arclength_minus_Hermite']) for q in geo if q['variant']=='on_the_fly_coupled21')),
            warnings=warnings,
            formal_clock_continuation_beyond_25p6=True,
            scientific_decision=(
              'AUTHENTIC_21D_ON_THE_FLY_PRIMITIVE__SENSITIVITY_FLAGGED__NO_MATERIAL_HOLONOMY'
              if warnings else
              'AUTHENTIC_21D_ON_THE_FLY_PRIMITIVE__CONTROLLED_11D_STEP_EFFECT__NO_MATERIAL_HOLONOMY'),
            limitations=[
                'Accumulator is mathematically passive, but it can affect adaptive-step acceptance through the solver error norm.',
                'The 10D Euclidean bivector residual is not the full 11D physical-state distance.',
                'Dense output is an interpolant of the numerical ODE solution, not a mathematically exact continuum certificate.',
                'Both the locked trajectory and augmented solve use the same authenticated original Q550 model.',
                'The endpoints beyond xi=25.6 are a formal-model-clock continuation, not independent observation.',
                'No global minimum, exact recurrence, material spin, or material holonomy is certified.'
            ]
        )
        (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text(
          '# Q550 v18.53t — On-the-fly 21D primitive\n\n'
          +'Decision: `'+report['scientific_decision']+'`\n\n'
          +'Physical 11D-only: 6; coupled 21D: 6; locked 10D quadratures: 6.\n\n'
          +'Maximum physical difference: '+repr(report['max_physical_state_drift'])+'\n\n'
          +'Maximum coupled vs locked primitive: '+repr(report['max_coupled_locked_accumulator_difference'])+'\n\n'
          +'Maximum dense vs Hermite 10D residual: '+repr(report['max_dense_vs_hermite_residual_difference'])+'\n\n'
          +'Limitations are detailed in `report.json`; no material closure claim.\n')
        print('V1853T',report['scientific_decision'],
              'PHYSICAL11',report['completed_physical11_ivps'],
              'COUPLED21',report['completed_coupled21_ivps'],
              'LOCKED10',report['completed_locked10_ivps'],flush=True)
        return 0
    except Exception as e:
        report.update(scientific_decision='FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE',
                      error=str(e),traceback=traceback.format_exc())
        (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print('V1853T_FAILED',str(e),flush=True)
        return 2

if __name__=='__main__':sys.exit(main())
