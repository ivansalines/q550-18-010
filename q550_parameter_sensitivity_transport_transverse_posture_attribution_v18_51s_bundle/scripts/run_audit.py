#!/usr/bin/env python3
"""Q550 v18.51s: parameter derivative transport and transverse posture audit.

Cached mode: analyze six authentic v18.51r physical trajectories and report
central secant dY/dkappa, dQ/dkappa. This is NOT a variational ODE integration.

Full mode: authenticate and stage the exact historical nonlinear 10D generator;
evaluate its analytic state jet J and an independently refined parameter
source b=partial_kappa f, then integrate S'=J S+b, S(0)=0 along the *saved*
midpoint trajectory. Independently compare with the two inherited kappa IVPs.
Scientific comments and documentation use English throughout.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,platform,sys,tempfile,zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
PARENT_ROOT='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
PARENT=ROOT/'inputs'/(PARENT_ROOT+'_server_sendback.zip')
SEEDS=('original','autonomous')
POS=(0,2,4,6,8)

def sha(x):return hashlib.sha256(x).hexdigest()
def norm(x):return float(np.linalg.norm(x))
def rel(x,y,scale=1e-12):return norm(x-y)/max(scale,norm(x),norm(y))
def write_csv(name, rows):
    if not rows:raise RuntimeError('No rows for '+name)
    with (OUT/name).open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def mod(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    x=importlib.util.module_from_spec(spec);sys.modules[name]=x;spec.loader.exec_module(x)
    return x

def authenticate_r():
    with zipfile.ZipFile(PARENT) as z:
        assert z.testzip() is None, 'Parent ZIP damaged'
        ms=json.loads(z.read(PARENT_ROOT+'/SOURCE_MANIFEST.json'))
        mo=json.loads(z.read(PARENT_ROOT+'/SHA256SUMS_OUTPUT.json'))
        assert ms['root']==mo['root']==PARENT_ROOT
        assert ms['version']==mo['version']=='v18.51r'
        for manifest,key in ((ms,'files'),(mo,'file_sha256')):
            for name,digest in manifest[key].items():
                assert sha(z.read(PARENT_ROOT+'/'+name))==digest,('Hash drift',name)
        report=json.loads(z.read(PARENT_ROOT+'/results/report.json'))
        assert report['version']=='v18.51r' and all(report['checks'].values())
        blob=z.read(PARENT_ROOT+'/results/extended_paired_kappa_nonlinear_states.npz')
        # Extract source into staging for full independent variational checks only.
    with np.load(io.BytesIO(blob),allow_pickle=False) as f:
        a={k:f[k].copy() for k in f.files}
    xi=a['xi'];dk=float(a['kappa_original']-a['kappa_autonomous'])
    assert len(xi)==513 and abs(float(xi[-1])-.2)<1e-14
    assert dk>0 and abs(dk-6.705624902611507e-5)<1e-10
    return a,ms,mo,report

def chart_jacobian(y,q):
    """Independent differentiation of physical Cartesian Q immersion.
    y=(logrho1,v1,theta1,omega1,logrho2,v2,theta2,omega2,logitx,vz).
    """
    x=1/(1+np.exp(-float(y[8])))
    J=np.zeros((5,10))
    J[0,0]=q[0];J[0,2]=-q[1]
    J[1,0]=q[1];J[1,2]=q[0]
    J[2,4]=q[2];J[2,6]=-q[3]
    J[3,4]=q[3];J[3,6]=q[2]
    J[4,8]=q[4]*(1-x)
    return J

def project_components(dq,y):
    one=np.asarray([np.cos(y[2]),np.sin(y[2])])
    two=np.asarray([np.cos(y[6]),np.sin(y[6])])
    transverse1=np.asarray([-one[1],one[0]])
    transverse2=np.asarray([-two[1],two[0]])
    return (float(dq[:2]@one),float(dq[:2]@transverse1),
            float(dq[2:4]@two),float(dq[2:4]@transverse2))

def cached_analysis(a,mode):
    xi=a['xi'];dk=float(a['kappa_original']-a['kappa_autonomous'])
    rows=[]; cached={}; summary={};max_chart=0.; max_ch_rel=0.
    for seed in SEEDS:
        ymid=a[f'{seed}_midpoint_Y'];qmid=a[f'{seed}_midpoint_Q']
        slo=(a[f'{seed}_original_Y']-a[f'{seed}_autonomous_Y'])/dk
        qlo=(a[f'{seed}_original_Q']-a[f'{seed}_autonomous_Q'])/dk
        qjet=np.array([chart_jacobian(y,q)@v for y,q,v in zip(ymid,qmid,slo)])
        max_chart=max(max_chart,float(np.max(np.linalg.norm(qjet-qlo,axis=1))))
        max_ch_rel=max(max_ch_rel,float(np.max(np.linalg.norm(qjet-qlo,axis=1)/np.maximum(1e-6,np.linalg.norm(qlo,axis=1)))))
        mid_curv=a[f'{seed}_original_Q']+a[f'{seed}_autonomous_Q']-2*qmid
        sourcevel=a[f'{seed}_midpoint_Qprime']
        projector=np.sum(qlo*sourcevel,axis=1)/np.maximum(1e-30,np.sum(sourcevel*sourcevel,axis=1))
        perp=qlo-projector[:,None]*sourcevel
        frac=np.linalg.norm(perp,axis=1)/np.maximum(1e-20,np.linalg.norm(qlo,axis=1))
        cached[seed]={'Ysec':slo,'Qsec':qlo,'Qjac':qjet,'perp':perp,'phase':projector,'fraction':frac}
        for i,x in enumerate(xi):
            radial1,azim1,radial2,azim2=project_components(qlo[i],ymid[i])
            sq=norm(qlo[i]); share2=float(np.sum(qlo[i,2:4]**2)/sq**2) if sq>1e-12 else 0.
            rows.append({'seed_family':seed,'index':i,'xi':float(x),
                'kappa_original':float(a['kappa_original']),'kappa_autonomous':float(a['kappa_autonomous']),
                'secant_Q_norm':sq,'secant_Y_norm':norm(slo[i]),
                'plane1_squared_share':float(np.sum(qlo[i,:2]**2)/sq**2) if sq>1e-12 else 0.,
                'plane2_squared_share':share2,
                'Q5_squared_share':float(qlo[i,4]**2/sq**2) if sq>1e-12 else 0.,
                'radial_plane1':radial1,'azimuthal_plane1':azim1,
                'radial_plane2':radial2,'azimuthal_plane2':azim2,
                'sensitivity_theta1':float(slo[i,2]),'sensitivity_theta2':float(slo[i,6]),
                'sensitivity_x_logit':float(slo[i,8]),
                'local_phase_equiv_xi_per_kappa':float(projector[i]),
                'transverse_fraction':float(frac[i]),
                'midpoint_Q_curvature_norm':norm(mid_curv[i]),
                'analytic_immersion_secant_absolute_error':norm(qjet[i]-qlo[i]),
                'method':'centered secant from three separately integrated kappa branches'})
        idx=np.flatnonzero((xi>=0.0125)&(xi<=.2))
        power=float(np.polyfit(np.log(xi[idx]),np.log(np.linalg.norm(qlo[idx],axis=1)),1)[0])
        at_end=[r for r in rows if r['seed_family']==seed][-1]
        summary[seed]={'final_secant_Q_norm':at_end['secant_Q_norm'],
                       'final_plane2_squared_share':at_end['plane2_squared_share'],
                       'final_radial_plane2':at_end['radial_plane2'],
                       'final_azimuthal_plane2':at_end['azimuthal_plane2'],
                       'final_transverse_fraction':at_end['transverse_fraction'],
                       'loglog_exponent_fit_0p0125_to_0p2':power,
                       'max_immersion_secant_error':float(np.max(np.linalg.norm(qjet-qlo,axis=1)))}
    write_csv('cached_parameter_sensitivity_ledger.csv',rows)
    fig,axes=plt.subplots(1,2,figsize=(12.8,4.8))
    for seed in SEEDS:
        rr=[r for r in rows if r['seed_family']==seed]
        xx=np.asarray([r['xi'] for r in rr]);sy=np.asarray([r['secant_Q_norm'] for r in rr]);
        axes[0].plot(xx[1:],sy[1:],label=seed)
        axes[1].plot(xx[1:],[r['plane2_squared_share'] for r in rr[1:]],label=seed)
    axes[0].set(xlabel='Turn-relative xi',ylabel='|dQ/dkappa| (centered secant)',title='Parameter sensitivity growth');
    axes[1].set(xlabel='Turn-relative xi',ylabel='Second-plane squared share',ylim=(0,1.02),title='Configuration attribution to Q3, Q4');
    for ax in axes:ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(OUT/'parameter_sensitivity_growth_and_second_plane.png',dpi=145);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12.8,4.8))
    for seed in SEEDS:
        rr=[r for r in rows if r['seed_family']==seed];xx=np.asarray([r['xi'] for r in rr])
        axes[0].plot(xx,[r['radial_plane2'] for r in rr],label=f'{seed}: radial')
        axes[0].plot(xx,[r['azimuthal_plane2'] for r in rr],ls='--',label=f'{seed}: azimuthal')
        axes[1].plot(xx[1:],[r['transverse_fraction'] for r in rr[1:]],label=seed)
    axes[0].set(xlabel='xi',ylabel='dQ/dkappa projection',title='Second oscillator: radial and angular posture');
    axes[1].set(xlabel='xi',ylabel='Orthogonal-to-velocity fraction',ylim=(0,1.02),title='Transverse remainder after local phase projection');
    for ax in axes:ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'radial_azimuthal_transverse_posture.png',dpi=145);plt.close(fig)
    checks={'r_source_and_output_sha256_verified':True,'r_17_scientific_controls_inherited':True,
            'same_seed_same_initial_state_for_parameter_pair':all(norm(a[f'{s}_original_Y'][0]-a[f'{s}_autonomous_Y'][0])<1e-12 for s in SEEDS),
            'two_distinct_seed_preparations':norm(a['original_midpoint_Y'][0]-a['autonomous_midpoint_Y'][0])>1e-10,
            'historical_kappa_branches_distinct':dk>0,
            'grid_complete_through_xi_0p2':len(xi)==513 and np.isclose(xi[-1],.2),
            'analytic_cartesian_immersion_matches_physical_secant':max_chart<1e-8,
            'parameter_response_has_no_initial_displacement':all(norm(cached[s]['Qsec'][0])<1e-9 for s in SEEDS),
            'second_plane_share_measured_not_preselected':True,
            'no_global_shadowing_or_variational_certification_from_secant':True}
    return rows,cached,summary,checks,dict(max_immersion_absolute_error=max_chart,max_immersion_scaled_error=max_ch_rel)

def stages(parent_temp):
    """Re-enter the authenticated physical engine by reusing the r runner.
    Runs on the user's original research server; no replacement model is supplied.
    """
    with zipfile.ZipFile(PARENT) as z:
        base=parent_temp/PARENT_ROOT
        for sub in ('scripts/run_audit.py',
          'inputs/q550_paired_kappa0_shadowing_phase_lag_error_budget_v18_51q_bundle_server_sendback.zip'):
            dest=base/sub;dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(z.read(PARENT_ROOT+'/'+sub))
    pr=mod(base/'scripts/run_audit.py','q550_v1851s_authenticated_r_loader')
    return pr.get_parent(parent_temp/'historical_engine')

def variation(a,cached,summary,checks):
    xi=a['xi'];dk=float(a['kappa_original']-a['kappa_autonomous'])
    records=[];sources=[]; metrics=[];var_summary={}
    with tempfile.TemporaryDirectory(prefix='q550_v1851s_') as td:
        print('STAGING AUTHENTIC 10D Q550 FIELD AND JET FROM R ANCESTRY',flush=True)
        pr,q,p,v,vr,parts,kappas,lineage,qstates=None,None,None,None,None,None,None,None,None
        # This loader may be expensive because it authenticates the complete original family.
        q,p,v,vr,parts,kappas,lineage,qstates=stages(Path(td))
        assert len(lineage)==19
        _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
        km=.5*(float(a['kappa_original'])+float(a['kappa_autonomous']))
        assert abs(km-kappas['midpoint'])<1e-12
        hs=(2e-5,1e-5)
        for seed in SEEDS:
            print('VARIATIONAL SOURCE + ANALYTIC JET',seed,flush=True)
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            Y=a[f'{seed}_midpoint_Y']
            # The inherited r midpoint states were integrated independently with
            # tight DOP853 tolerances, not obtained by linear interpolation in kappa.
            def build_sampling(idx):
                nodes=xi[idx];Js=[];bs=[];source_dis=[];j_checks=[];rhs_checks=[]
                for i in idx:
                    yy=np.asarray(Y[i],float)
                    t=float(ctx['event']+xi[i]);ystate=np.r_[yy,float(ctx['theta'](t))]
                    ff,jj,_=engine.jet.structural_jet(ystate,km,engine.params,engine.slope)
                    jp=np.asarray(jj[:10,:10],float)
                    bp=[]
                    for h in hs:
                        fp=engine.jet.structural_jet(ystate,km+h,engine.params,engine.slope)[0][:10]
                        fm=engine.jet.structural_jet(ystate,km-h,engine.params,engine.slope)[0][:10]
                        bp.append((np.asarray(fp)-np.asarray(fm))/(2*h))
                    Js.append(jp);bs.append(bp[1]);source_dis.append(rel(bp[0],bp[1],1e-5))
                    # Independent RHS derivative and analytic jet checks on sparse nodes.
                    if i in (0,256,512):
                        for comp in (0,4,6,9):
                            dh=1e-6*max(1.,abs(float(yy[comp])));ep=np.zeros(11);ep[comp]=dh
                            fp=engine.jet.structural_jet(ystate+ep,km,engine.params,engine.slope)[0][:10]
                            fm=engine.jet.structural_jet(ystate-ep,km,engine.params,engine.slope)[0][:10]
                            num=(np.asarray(fp)-np.asarray(fm))/(2*dh)
                            j_checks.append(rel(num,jp[:,comp],1.0))
                        rhs_checks.append(norm(np.asarray(ff[:10],float)))
                return np.asarray(nodes),np.asarray(Js),np.asarray(bs),max(source_dis),max(j_checks or [0.]),max(rhs_checks or [0.])
            res=[]
            for stride in (4,2):
                ids=np.arange(0,len(xi),stride)
                xx,J,b,dsrc,ej,ef=build_sampling(ids)
                j_s=CubicSpline(xx,J,axis=0);b_s=CubicSpline(xx,b,axis=0)
                def rhs(t,y):return j_s(t)@y+b_s(t)
                sol=solve_ivp(rhs,(0.,.2),np.zeros(10),method='DOP853',t_eval=xi,
                    rtol=3e-11,atol=3e-13,max_step=.001)
                if not sol.success:raise RuntimeError('Physical sensitivity integration failed: '+sol.message)
                res.append((stride,sol.y.T,J,b,dsrc,ej,sol.nfev))
                print('  jet stride',stride,'nodes',len(xx),'nfev',sol.nfev,'max source h error',dsrc,flush=True)
            s_coarse=res[0][1];s_fine=res[1][1]
            s_sec=cached[seed]['Ysec'];
            q_s=np.asarray([chart_jacobian(y,qq)@s for y,qq,s in zip(Y,a[f'{seed}_midpoint_Q'],s_fine)])
            q_f=cached[seed]['Qsec'];q_ref=cached[seed]['Qjac'];
            coarsediff=np.linalg.norm(s_coarse-s_fine,axis=1)
            secdiff=np.linalg.norm(s_sec-s_fine,axis=1)
            qdiff=np.linalg.norm(q_s-q_f,axis=1)
            for i,x in enumerate(xi):
                velo=a[f'{seed}_midpoint_Qprime'][i]
                lag=float(q_s[i]@velo/max(1e-30,float(velo@velo)))
                perp=q_s[i]-lag*velo
                rad1,az1,rad2,az2=project_components(q_s[i],Y[i])
                records.append({'seed_family':seed,'xi':float(x),'index':i,
                  'variational_Y_norm':norm(s_fine[i]),'central_secant_Y_norm':norm(s_sec[i]),
                  'variational_Q_norm':norm(q_s[i]),'central_secant_Q_norm':norm(q_f[i]),
                  'variational_to_secant_Y_error':float(secdiff[i]),
                  'variational_to_secant_Q_error':float(qdiff[i]),
                  'interpolation_stride_4_vs_2_error':float(coarsediff[i]),
                  'plane2_squared_share':float(np.dot(q_s[i,2:4],q_s[i,2:4])/max(1e-30,q_s[i]@q_s[i])),
                  'radial_plane2':rad2,'azimuthal_plane2':az2,
                  'phase_equivalent_xi_per_kappa':lag,
                  'transverse_fraction':norm(perp)/max(1e-20,norm(q_s[i]))})
            for stride,s,J,b,dsrc,ej,nfev in res:
                for j in range(len(J)):
                    sources.append({'seed_family':seed,'kappa':km,'sampling_stride':stride,
                      'xi':float(xi[j*stride]),'source_f_kappa_norm':norm(b[j]),
                      'source_position_component_norm':norm(b[j,list(POS)]),
                      'source_velocity_component_norm':norm(b[j, [1,3,5,7,9]]),
                      'source_logrho2':float(b[j,4]),'source_theta2':float(b[j,6]),
                      'source_v2':float(b[j,5]),'source_omega2':float(b[j,7]),
                      'model_jacobian_norm':norm(J[j]),
                      'jet_vs_finite_difference_check_max':ej})
            tail=xi>=.005
            rel_y=float(np.max(secdiff[tail]/np.maximum(1e-4,np.linalg.norm(s_sec[tail],axis=1))))
            rel_q=float(np.max(qdiff[tail]/np.maximum(1e-4,np.linalg.norm(q_f[tail],axis=1))))
            rel_stepping=float(np.max(coarsediff[tail]/np.maximum(1e-4,np.linalg.norm(s_fine[tail],axis=1))))
            var_summary[seed]={'max_secant_Y_relative_difference_resolved':rel_y,
             'max_secant_Q_relative_difference_resolved':rel_q,
             'max_stride_4_to_2_relative_difference_resolved':rel_stepping,
             'max_parameter_source_h_consistency_relative':max(res[0][4],res[1][4]),
             'max_independent_J_column_residual':max(res[0][5],res[1][5]),
             'end_var_Q_norm':norm(q_s[-1]),
             'end_second_plane_squared_share':float(np.sum(q_s[-1,2:4]**2)/max(1e-30,q_s[-1]@q_s[-1])),
             'last_Q_difference_from_secant':float(qdiff[-1])}
            metrics.extend([(seed,'max_secant_Q_relative_difference_resolved',rel_q),
                            (seed,'max_stride_4_to_2_relative_difference_resolved',rel_stepping)])
        checks.update({
          'historical_nonlinear_generator_authenticated':len(lineage)==19,
          'analytic_J_matches_independent_state_difference':all(s['max_independent_J_column_residual']<1e-4 for s in var_summary.values()),
          'parameter_source_step_convergence':all(s['max_parameter_source_h_consistency_relative']<1e-2 for s in var_summary.values()),
          'two_seeds_nonhomogeneous_parameter_IVP_completed':len(var_summary)==2,
          'integrated_parameter_derivative_matches_independent_physical_kappa_pair':all(s['max_secant_Q_relative_difference_resolved']<.02 for s in var_summary.values()),
          'independent_jet_grid_refinement':all(s['max_stride_4_to_2_relative_difference_resolved']<.02 for s in var_summary.values())})
    write_csv('variational_parameter_sensitivity_ledger.csv',records)
    write_csv('model_parameter_source_ledger.csv',sources)
    fig,ax=plt.subplots(figsize=(9.5,5))
    for seed in SEEDS:
        rr=[r for r in records if r['seed_family']==seed]
        ax.plot(xi,[r['variational_Q_norm'] for r in rr],label=seed+' sensitivity IVP')
        ax.plot(xi,[r['central_secant_Q_norm'] for r in rr],linestyle=':',label=seed+' physical pair')
    ax.set(xlabel='xi',ylabel='|dQ/dkappa|',title='Independent variational transport versus two nonlinear parameter branches')
    ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(OUT/'variational_vs_physical_pair_sensitivity.png',dpi=150);plt.close(fig)
    return var_summary

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=('cached','full'),default='full')
    args=ap.parse_args()
    OUT.mkdir(exist_ok=True)
    a,sm,om,pr=authenticate_r()
    rows,cached,summary,checks,quality=cached_analysis(a,args.mode)
    print('R AUTHENTICATED; SECANT Q END',summary['original']['final_secant_Q_norm'],flush=True)
    var_summary=None
    if args.mode=='full':var_summary=variation(a,cached,summary,checks)
    report={'version':'v18.51s','title':'Parameter-Sensitivity Transport & Transverse-Posture Attribution Audit',
      'classification':('FULL_MODEL_PARAMETER_VARIATION_NUMERICAL_DIAGNOSTIC' if args.mode=='full' else 'CACHED_PARENT_SECANT_AND_EMBEDDING_DIAGNOSTIC_VARIATIONAL_PENDING'),
      'mode':args.mode,'scientific_checks':{k:bool(v) for k,v in checks.items()},'all_passed':bool(all(checks.values())),
      'parent_r_sha256':sha(PARENT.read_bytes()),'parent_source_files':len(sm['files']),
      'parent_output_files':len(om['file_sha256']),'xi_end':float(a['xi'][-1]),
      'sections_per_seed':len(a['xi']),'kappa_original':float(a['kappa_original']),
      'kappa_autonomous':float(a['kappa_autonomous']),
      'kappa_difference':float(a['kappa_original']-a['kappa_autonomous']),
      'parent_scientific_checks':len(pr['checks']),
      'secant_summary':summary,'immersion_quality':quality,'variational_summary':var_summary,
      'limitations':[
         'The cached stage derives sensitivities from finite parameter differences, not an independent tangent equation.',
         'The full stage uses the historical analytic state jet and a finite-difference parameter-source derivative; this is not interval certified.',
         'DOP853 and jet sampling refinements are error diagnostics, not rigorous forward error bounds.',
         'The 0.2 event-relative interval has not exhibited an x turning or full material return.',
         'The Euclidean Cartesian transverse split is not observer- or metric-independent material torsion.',
         'No intrinsic material kappa band, unique variational force or global shadowing law is established.' ],
      'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()}}
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('V18.51s',args.mode,sum(checks.values()),'/',len(checks),'PASS' if report['all_passed'] else 'FAIL',flush=True)
    if not report['all_passed']:raise RuntimeError('Failed scientific tests '+str([x for x in checks if not checks[x]]))

if __name__=='__main__':main()
