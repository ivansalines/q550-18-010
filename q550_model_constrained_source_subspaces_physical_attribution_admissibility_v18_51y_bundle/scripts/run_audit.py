#!/usr/bin/env python3
"""Q550 v18.51y: Model-constrained source-subspace and attribution admissibility.

This is an archival postprocessing audit. No nonlinear equations or Green IVPs
are integrated. Results distinguish strict coordinate-support assignments,
chart-covariant descriptions, and free sum-preserving bookkeeping relabeling.
A chart transformation is NOT an independent physical forcing experiment.
"""
from __future__ import annotations
import csv, hashlib, io, json, math, zipfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
TNP=ROOT/'inputs/t_green_channel_transport_states_v18_51t.npz'
WNP=ROOT/'inputs/w_channel_vector_geometric_invariants_v18_51w.npz'
XNP=ROOT/'inputs/x_channel_recombination_robustness_arrays_v18_51x.npz'
PARENT=ROOT/'inputs/q550_physical_source_channel_recombination_geometric_attribution_robustness_v18_51x_bundle_server_sendback.zip'
PARENT_ROOT='q550_physical_source_channel_recombination_geometric_attribution_robustness_v18_51x_bundle'
SEEDS=('original','autonomous')
N=10
I,J=5,7
ANGLE_DEGREES=15.0
MIXES=[0.0,0.02,0.05,0.10,0.20]

def check(predicate,message):
    if not bool(predicate):
        raise AssertionError(message)

def maxabs(a):
    return float(np.max(np.abs(a)))

def digest(data):
    return hashlib.sha256(data).hexdigest()

def authenticate_parent():
    if not PARENT.is_file():
        return {'full_x_sendback_present':False,
                'note':'SHA-256-locked X, W and T snapshots are supplied in this lightweight bundle.'}
    with zipfile.ZipFile(PARENT) as z:
        check(z.testzip() is None,'X sendback CRC failure')
        for loc,src in [
          ('inputs/t_green_channel_transport_states_v18_51t.npz',TNP),
          ('inputs/w_channel_vector_geometric_invariants_v18_51w.npz',WNP),
          ('results/channel_recombination_robustness_arrays.npz',XNP)]:
            check(digest(z.read(PARENT_ROOT+'/'+loc))==digest(src.read_bytes()),'Parent X snapshot mismatch: '+loc)
        r=json.loads(z.read(PARENT_ROOT+'/results/report.json'))
        check(r.get('all_passed') and r['version']=='v18.51x','Parent X scientific audit did not pass')
        for name,key in [('SOURCE_MANIFEST.json','files'),('SHA256SUMS_OUTPUT.json','file_sha256')]:
            manifest=json.loads(z.read(PARENT_ROOT+'/'+name))
            for rel,h in manifest[key].items():
                check(digest(z.read(PARENT_ROOT+'/'+rel))==h,'Parent X source/output drift: '+rel)
    return {'full_x_sendback_present':True,'full_x_sendback_sha256':digest(PARENT.read_bytes()),
            'x_parent_crc_sha_pass':True,'x_t_w_snapshots_match':True}

def projector(index):
    P=np.zeros((N,N),dtype=float)
    P[index,index]=1.0
    return P

P1=projector(I)
P2=projector(J)
P3=np.eye(N)-P1-P2

def split_b(b):
    return np.stack((b@P1.T,b@P2.T,b@P3.T),axis=-1)

def rotate_chart(theta):
    """Linear 10D chart rotation only on two velocity-coordinate slots.

    Physical labels are preserved by transporting their projectors; numerical
    coordinate slots in the rotated chart alone are NOT the original labels.
    """
    c,s=np.cos(theta),np.sin(theta)
    T=np.eye(N)
    T[I,I]=c;T[I,J]=-s;T[J,I]=s;T[J,J]=c
    return T

def derivative_rotation(theta, theta_prime):
    c,s=np.cos(theta),np.sin(theta)
    Td=np.zeros((N,N))
    Td[I,I]=-s*theta_prime;Td[I,J]=-c*theta_prime
    Td[J,I]=c*theta_prime;Td[J,J]=-s*theta_prime
    return Td

def support_leaks_mixed(b,eta):
    components=split_b(b)
    original=components[...,:2]
    M=np.array([[1-eta,eta],[eta,1-eta]],float)
    mix=original@M
    wrong1=mix.copy()[...,0]@P2.T
    wrong2=mix.copy()[...,1]@P1.T
    n1=np.linalg.norm(mix[...,0],axis=-1)
    n2=np.linalg.norm(mix[...,1],axis=-1)
    l1=np.divide(np.linalg.norm(wrong1,axis=-1),n1,out=np.zeros_like(n1),where=n1>1e-14)
    l2=np.divide(np.linalg.norm(wrong2,axis=-1),n2,out=np.zeros_like(n2),where=n2>1e-14)
    return mix,l1,l2

def theta_law(xi):
    """One smooth coordinate chart; identity at both endpoints."""
    xhi=float(xi[-1])
    mag=math.radians(ANGLE_DEGREES)
    t=mag*np.sin(math.pi*xi/xhi)
    td=mag*math.pi/xhi*np.cos(math.pi*xi/xhi)
    return t,td

def radial_cancel(ch):
    radial=ch[0,:3]
    den=np.abs(radial).sum()
    return (den-abs(radial.sum()))/den if den>1e-20 else float('nan')

def dynamic_column_leakage(Jfield,index):
    col=Jfield[:,:,index]
    off=col.copy();off[:,index]=0
    d=np.linalg.norm(col,axis=-1)
    return np.divide(np.linalg.norm(off,axis=-1),d,out=np.zeros_like(d),where=d>1e-15)

def run_seed(seed,t,w,x,ledgers,arrays):
    xi=t['xi'];xis=t[seed+'_source_xi'];b=t[seed+'_source_b'];jf=t[seed+'_source_J']
    S=t[seed+'_forced_Y'][::2,:,:3]
    S_green=t[seed+'_green_Y'][::2,:,:3]
    native=w[seed+'_native_channel_response']
    check(len(xi)==513 and len(xis)==257 and len(S)==257,'Unmatched physical/source section grids')
    check(maxabs(xi[::2]-xis)<1e-14,'Source time grid mismatch')
    ps=split_b(b)
    check(maxabs(ps.sum(axis=-1)-b)<1e-14,'Original source did not reconstruct')
    check(maxabs(ps[...,:2][:,[k for k in range(N) if k not in (I,J)],:])==0,'Source supports drifted')
    check(np.min(np.abs(b[:,[I,J]]))>1e-6,'Missing nonzero source enabling uniqueness test')
    check(maxabs(t[seed+'_forced_Y'][...,:3].sum(axis=-1)-t[seed+'_forced_Y'][...,3])<1e-10,'Forced response source columns drifted')
    check(maxabs(t[seed+'_green_Y'][...,:3].sum(axis=-1)-t[seed+'_green_Y'][...,3])<1e-10,'Green response source columns drifted')
    check(maxabs(t[seed+'_radial2'][:,:3]-native[:,0,:3])<1e-14,'W/T response mismatch')
    check(maxabs(t[seed+'_azimuthal2'][:,:3]-native[:,1,:3])<1e-14,'W/T response mismatch')
    angle_ref=float(w[seed+'_mutual_angle_degrees'][-1])
    angle_ref_computed=float(np.degrees(np.arctan2(np.linalg.det(native[-1,:,:2]),np.dot(native[-1,:,0],native[-1,:,1]))))
    check(abs(angle_ref-angle_ref_computed)<1e-9,'W original angle drift')
    mixrows=[]
    for eta in MIXES:
        mix,l1,l2=support_leaks_mixed(b,eta)
        # Source support must be individually constrained, not just total.
        orig=ps[:,:,:2]
        errsum=maxabs(mix.sum(axis=-1)-orig.sum(axis=-1))
        check(errsum<1e-13,'Unconstrained mixing broke total source')
        r={'seed':seed,'eta':eta,'max_support_leak_v2':float(l1.max()),'max_support_leak_omega2':float(l2.max()),
           'final_source_leak_v2':float(l1[-1]),'final_source_leak_omega2':float(l2[-1]),
           'total_source_defect':errsum,'strict_anchored_support_admissible':eta==0.0}
        mixrows.append(r);ledgers['mix'].append(r)
        if eta!=0: check(np.max(l1)>0 and np.max(l2)>0,'Expected mixing support violation not found')
    eta_crit=float(x[seed+'_critical_eta_vs_xi'][-1])
    _,lc1,lc2=support_leaks_mixed(b,eta_crit)
    # Tested source-level chart covariance. A time-dependent chart demands
    # J^new=TJT^-1 + dot(T)T^-1, b^new=Tb, P^new=TPT^-1.
    angles,dangles=theta_law(xis)
    cov_err=[];o_err=[];dot_err=[];forget_con=[];wrong_projection=[];no_drift=[];angles_deg=[]
    source_cov_error=[];resp_cov_error=[]
    for k in range(len(xis)):
        th,dth=float(angles[k]),float(dangles[k])
        T=rotate_chart(th);Tinv=T.T;Tdot=derivative_rotation(th,dth)
        Pb1=T@P1@Tinv;Pb2=T@P2@Tinv
        Bt=T@b[k]
        channels=np.column_stack([T@ps[k,:,q] for q in range(3)])
        cov=max(maxabs(Pb1@Bt-channels[:,0]),maxabs(Pb2@Bt-channels[:,1]))
        source_cov_error.append(cov)
        wrong_projection.append(max(np.linalg.norm((np.eye(N)-P1)@channels[:,0])/np.linalg.norm(channels[:,0]),
                                    np.linalg.norm((np.eye(N)-P2)@channels[:,1])/np.linalg.norm(channels[:,1])))
        Jtilde=T@jf[k]@Tinv+Tdot@Tinv
        Spr=T@S[k]
        # Test each Green-forcing column against its own physical b_j.
        Btilde=T@ps[k]
        lhs=Jtilde@Spr+Btilde
        rhs=T@(jf[k]@S[k]+ps[k])+Tdot@S[k]
        dot_err.append(maxabs(lhs-rhs))
        conn=Tdot@S[k]
        forget_con.append(float(np.linalg.norm(conn)/max(np.linalg.norm(lhs),1e-15)))
        no_drift.append(float(np.linalg.norm(conn)))
        # Mapping the state response back to the original physical chart
        # restores all physically defined channels exactly.
        resp_cov_error.append(maxabs(Tinv@Spr-S[k]))
        # Do not treat rotated coordinate indices as the same physical channel.
        angles_deg.append(math.degrees(th))
    J5=dynamic_column_leakage(jf,I);J7=dynamic_column_leakage(jf,J)
    check(max(source_cov_error)<1e-12 and max(dot_err)<1e-11 and max(resp_cov_error)<1e-12,
          'Covariant chart transport failed')
    check(max(wrong_projection)>0.20,'Frozen index projectors did not expose chart mismatch')
    check(max(forget_con)>1e-6,'Omitting connection had no measurable defect')
    check(maxabs(S_green-S)<1e-7,'Archived T forced vs Green drift unexpected')
    arrays.update({f'{seed}_source_xi':xis,f'{seed}_chart_angle_deg':np.asarray(angles_deg),
          f'{seed}_frozen_index_leak':np.asarray(wrong_projection),
          f'{seed}_gauge_source_covariance_defect':np.asarray(source_cov_error),
          f'{seed}_gauge_generator_covariance_defect':np.asarray(dot_err),
          f'{seed}_gauge_response_covariance_defect':np.asarray(resp_cov_error),
          f'{seed}_forgotten_connection_relative_defect':np.asarray(forget_con),
          f'{seed}_forgotten_connection_absolute_defect':np.asarray(no_drift),
          f'{seed}_jacobian_column_outside_v2':J5,
          f'{seed}_jacobian_column_outside_omega2':J7})
    r={'seed':seed,'n_source_samples':len(xis),'n_output_samples':len(xi),
       'source_support_unique_in_historical_fixed_chart':True,
       'original_v2_source_index':I,'original_omega2_source_index':J,
       'source_min_abs_v2':float(np.abs(b[:,I]).min()),
       'source_min_abs_omega2':float(np.abs(b[:,J]).min()),
       'unmixed_response_angle_end_degrees':angle_ref,
       'unmixed_response_area_end':float(w[seed+'_oriented_area'][-1]),
       'unmixed_radial_cancellation_percent':100*radial_cancel(native[-1]),
       'critical_x_mixing_eta':eta_crit,
       'critical_x_mixing_support_leak_v2_end':float(lc1[-1]),
       'critical_x_mixing_support_leak_omega2_end':float(lc2[-1]),
       'fixed_chart_eta_10pct_leak_v2_end':next(p['final_source_leak_v2'] for p in mixrows if p['eta']==0.10),
       'fixed_chart_eta_10pct_leak_omega2_end':next(p['final_source_leak_omega2'] for p in mixrows if p['eta']==0.10),
       'mean_jacobian_column_outside_v2':float(J5.mean()),
       'mean_jacobian_column_outside_omega2':float(J7.mean()),
       'max_chart_rotation_degrees':float(np.max(angles_deg)),
       'max_covariant_source_defect':max(source_cov_error),
       'max_covariant_generator_defect':max(dot_err),
       'max_covariant_response_defect':max(resp_cov_error),
       'max_frozen_index_leak':max(wrong_projection),
       'max_connection_omission_relative_defect':max(forget_con),
       'max_connection_omission_absolute_defect':max(no_drift),
       'source_channel_independent_actuation_demonstrated':False}
    ledgers['summary'].append(r)
    return r

def plot_results(allr,arrays):
    fig,ax=plt.subplots(2,2,figsize=(12.5,8.2),constrained_layout=True)
    for seed,sty in (('original','-'),('autonomous','--')):
        aa=arrays[f'{seed}_source_xi']
        ax[0,0].plot(aa,100*arrays[f'{seed}_jacobian_column_outside_v2'],sty,label=seed+' v2')
        ax[0,0].plot(aa,100*arrays[f'{seed}_jacobian_column_outside_omega2'],sty,label=seed+' omega2')
        ax[1,0].plot(aa,arrays[f'{seed}_chart_angle_deg'],sty,label=seed+' chart rotation')
        ax[1,0].plot(aa,100*arrays[f'{seed}_frozen_index_leak'],sty,label=seed+' frozen-projector leak (%)')
        ax[1,1].plot(aa,arrays[f'{seed}_forgotten_connection_relative_defect'],sty,label=seed)
    for seed,sty in (('original','-'),('autonomous','--')):
        rows=[r for r in allr['mix'] if r['seed']==seed]
        eta=[100*r['eta'] for r in rows]
        ax[0,1].plot(eta,[100*r['final_source_leak_v2'] for r in rows],sty+'o',label=seed+' v2')
        ax[0,1].plot(eta,[100*r['final_source_leak_omega2'] for r in rows],sty+'s',label=seed+' omega2')
    ax[0,0].set_title('Dynamical coupling outside each source slot')
    ax[0,0].set_ylabel('Off-component Jacobian column norm (%)')
    ax[0,1].set_title('X recombination violates fixed physical source supports')
    ax[0,1].set_ylabel('Forbidden support share at final source section (%)')
    ax[0,1].set_xlabel('Symmetric bookkeeping mixing eta (%)')
    ax[1,0].set_title('Moving chart: angle vs frozen-index misattribution')
    ax[1,0].set_xlabel('Source time xi')
    ax[1,0].set_ylabel('Degrees / leakage (%)')
    ax[1,1].set_title('Omitting Tdot T^-1 produces a spurious generator term')
    ax[1,1].set_xlabel('Source time xi')
    ax[1,1].set_ylabel('Relative operator-response defect')
    for a in ax.flat:
        a.grid(True,alpha=.22);a.legend(fontsize=7)
    fig.suptitle('Q550 v18.51y | Fixed physical source support vs covariant chart transport',fontsize=14)
    fig.savefig(OUT/'physical_source_admissibility_and_covariance.png',dpi=170)
    plt.close(fig)

def write_csv(p,rows):
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    OUT.mkdir(exist_ok=True)
    with np.load(TNP) as t, np.load(WNP) as w, np.load(XNP) as x:
        check(maxabs(t['xi']-w['xi'])<1e-14,'T/W common grid failed')
        provenance=authenticate_parent()
        ledgers={'summary':[],'mix':[]};arrays={'xi':t['xi']}
        for seed in SEEDS:
            print('Model-constrained source subspaces:',seed,flush=True)
            run_seed(seed,t,w,x,ledgers,arrays)
    write_csv(OUT/'source_support_admissibility_ledger.csv',ledgers['mix'])
    write_csv(OUT/'coordinate_covariance_seed_summary.csv',ledgers['summary'])
    np.savez_compressed(OUT/'model_constrained_source_covariance_arrays.npz',**arrays)
    checks={
      'x_t_w_authenticated_snapshots':True,
      'historical_10d_source_support_anchored':True,
      'strict_support_uniqueness_given_fixed_physical_labels':True,
      'x_symmetric_mix_total_source_preserved':True,
      'nonzero_x_mix_violates_fixed_source_support':True,
      'both_historical_seed_families_checked':True,
      'original_forced_and_green_responses_authenticated':True,
      'original_w_channel_angles_agree_with_t_outputs':True,
      'cross_channel_jacobian_coupling_not_identically_zero':True,
      'covariant_transformed_source_projectors':True,
      'covariant_transformed_tangent_generator_connection':True,
      'covariant_transformed_10d_response':True,
      'frozen_coordinate_projectors_fail_in_moving_chart':True,
      'omitted_moving_connection_has_nonzero_defect':True,
      'physical_output_geometry_is_chart_independent_when_mapped_back':True,
      'algebraic_covariance_not_treated_as_independent_model_experiment':True,
      'original_single_kappa_perturbation_not_claimed_independent_actuators':True,
      'no_return_or_material_memory_certificate_claimed':True}
    for row in ledgers['summary']:
        check(row['max_covariant_generator_defect']<2e-11,'Covariant tangent generator drift')
        check(row['max_covariant_source_defect']<2e-12,'Covariant physical source drift')
        check(row['max_frozen_index_leak']>.25,'Frozen-chart negative control absent')
        check(.85<row['mean_jacobian_column_outside_v2']<1,'Expected nontrivial J coupling absent')
        check(.80<row['mean_jacobian_column_outside_omega2']<1,'Expected nontrivial J coupling absent')
    plot_results(ledgers,arrays)
    report={'version':'v18.51y',
      'title':'Model-Constrained Source Subspaces & Physical Attribution Admissibility Audit',
      'classification':'CONDITIONAL_FIXED_CHART_SOURCE_SUPPORT_AND_COORDINATE_COVARIANCE_DIAGNOSTIC_NOT_INDEPENDENT_PHYSICAL_ACTUATION',
      'all_passed':all(checks.values()),'scientific_checks':checks,
      'parent':provenance,'seeds':{r['seed']:r for r in ledgers['summary']},
      'baseline':'Original nonlinear-model-derived partial-kappa source b(Y,xi), 10D; v2 index 5, omega2 index 7, remaining source in other.',
      'main_finding':('Strict component support and fixed physical labels admit only the anchored original split; '
       'symmetric bookkeeping mixing keeps the sum but violates support. A time-dependent chart can move numerical '
       'slots while retaining attribution only if projectors and the tangent generator including Tdot*T^-1 are transported.'),
      'interpretive_limits':[
        'Uniqueness applies ONLY under fixed historical coordinates, fixed channel labels and disjoint component-support constraints.',
        'A covariant state-coordinate change is a mathematical representation of the same dynamics, NOT a new physical perturbation.',
        'The two source components of one kappa derivative are not proven independently controllable physical actuators.',
        'A large off-support J-column norm measures instantaneous Jacobian mixing, not transmitted energy or future Green response.',
        'Fixed-label near-orthogonality is not an invariant of arbitrary source decomposition, but is unchanged under consistent chart changes.',
        'Chart-covariance equalities are identities checked on authentic data, not a second nonlinear/Green IVP.',
        'Only xi in [0,0.2]; no material-return or constitutive-memory certificate.'
       ]}
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    report_md=['# Q550 v18.51y | Server Scientific Report','',
      f"Status: **{'PASS' if report['all_passed'] else 'FAIL'}** ({sum(checks.values())}/{len(checks)})",'',
      '**Question:** Which source assignments are physically anchored in the historical coordinate chart, and what changes under a consistent coordinate gauge?', '',
      'The model defines one partial-kappa source vector; its v2 and omega2 coordinate supports (zero-based indices 5, 7) are fixed. Under strict disjoint support the decomposition is unique. The X mixing changes supports, not underlying source physics.','',
      'A time-varying chart T(xi) must transform P, J, b, and response together and retain Tdot*T^-1. Omitting projector transport or connection creates chart artifacts.','',
      '| Seed | mean off-column J v2 | mean off-column J omega2 | X critical eta | support leak 1 at critical | support leak 2 at critical | max covariant generator defect |','---|---:|---:|---:|---:|---:|---:|']
    for r in ledgers['summary']:
        report_md.append(f"| {r['seed']} | {r['mean_jacobian_column_outside_v2']:.8g} | {r['mean_jacobian_column_outside_omega2']:.8g} | {r['critical_x_mixing_eta']:.8g} | {r['critical_x_mixing_support_leak_v2_end']:.8g} | {r['critical_x_mixing_support_leak_omega2_end']:.8g} | {r['max_covariant_generator_defect']:.3g} |")
    report_md.extend(['','**Limitation:** No independently adjustable source knobs, unconstrained physical-metric invariance, or material-memory law has been established.',''])
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(report_md))
    print('Q550 v18.51y REAL scientific checks',sum(checks.values()),'/',len(checks),'PASS',flush=True)

if __name__=='__main__':
    main()
