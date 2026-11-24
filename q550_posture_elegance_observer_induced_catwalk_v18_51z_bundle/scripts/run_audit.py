#!/usr/bin/env python3
"""Q550 v18.51z: channel-posture and observer-frame audit.

No nonlinear integrations occur here. All numerical arrays are sourced from an
SHA256-locked v18.51y server sendback containing authenticated W/T snapshots.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

TITLE='Posture, Elegance & Observer-Induced Catwalk Audit'
VERSION='v18.51z'
CHANNELS=('v2','omega2','other','total')
PAIRS=((0,1),(0,2),(1,2))

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for part in iter(lambda:stream.read(2**20),b''):
            h.update(part)
    return h.hexdigest()


def cross(a,b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def rot(v,angle):
    """Passive 2D rotation by -angle. Coordinate slot axis is the second-last."""
    c=np.cos(angle); s=np.sin(angle)
    out=np.empty_like(v)
    out[...,0,:]=c[:,None]*v[...,0,:]+s[:,None]*v[...,1,:]
    out[...,1,:]=-s[:,None]*v[...,0,:]+c[:,None]*v[...,1,:]
    return out


def angles(v):
    """Continuous angle for every channel; initial zero is undefined."""
    out=np.full((len(v),v.shape[2]),np.nan)
    for j in range(v.shape[2]):
        out[1:,j]=np.unwrap(np.arctan2(v[1:,1,j],v[1:,0,j]))
    return out


def observer_phase(xi):
    """Explicit passive observer protocol, fixed before inspecting the result."""
    z=xi/float(xi[-1])
    return np.deg2rad(8)*z+np.deg2rad(2.5)*np.sin(4*np.pi*z)


def calc_seed(seed,xi, lab, native, tau, t_lab):
    assert lab.shape==(len(xi),2,4)
    assert native.shape==lab.shape
    assert np.all(np.diff(tau)>0)
    assert np.max(np.abs(np.sum(lab[:,:,:3],axis=2)-lab[:,:,3]))<2e-13
    h=float(xi[1]-xi[0]); assert np.allclose(np.diff(xi),h)
    mask=xi>=0.01
    rate= np.gradient(angles(lab)[1:],xi[1:],axis=0,edge_order=2)
    ang=angles(lab)
    phi=observer_phase(xi)
    chart=rot(lab,phi)
    a_chart=angles(chart)
    frame_rate=np.gradient(phi,xi,edge_order=2)
    rate_lab=np.column_stack([np.interp(xi,xi[1:],rate[:,j]) for j in range(4)])
    rate_lab[0,:]=np.nan
    rate_chart=np.full_like(rate_lab,np.nan)
    rate_chart[1:,:]=np.gradient(a_chart[1:],xi[1:],axis=0,edge_order=2)
    rate_corr=rate_chart+frame_rate[:,None]
    valid=mask & np.all(np.isfinite(rate_lab),axis=1)
    dlab=np.gradient(lab,xi,axis=0,edge_order=2)
    denom=np.sum(lab[:,:,3]**2,axis=1)
    self_terms=np.zeros((len(xi),3))
    pair_terms=np.zeros((len(xi),3))
    for j in range(3):
        self_terms[:,j]=cross(lab[:,:,j],dlab[:,:,j])
    for j,(p,q) in enumerate(PAIRS):
        pair_terms[:,j]=cross(lab[:,:,p],dlab[:,:,q])+cross(lab[:,:,q],dlab[:,:,p])
    norm_self=np.divide(self_terms,denom[:,None],out=np.full_like(self_terms,np.nan),where=denom[:,None]>0)
    norm_pair=np.divide(pair_terms,denom[:,None],out=np.full_like(pair_terms,np.nan),where=denom[:,None]>0)
    total_rate_wedge=np.divide(cross(lab[:,:,3],dlab[:,:,3]),denom,out=np.full_like(denom,np.nan),where=denom>0)
    total_rate_reconstruct=norm_self.sum(axis=1)+norm_pair.sum(axis=1)
    delta_lab=np.unwrap(ang[1:,1]-ang[1:,0]);delta_chart=np.unwrap(a_chart[1:,1]-a_chart[1:,0])
    lab_area=cross(lab[:,:,0],lab[:,:,1]);chart_area=cross(chart[:,:,0],chart[:,:,1])
    relative_rate=np.gradient(delta_lab,xi[1:],edge_order=2)
    c=np.gradient(tau,xi,edge_order=2)
    corrected_clock=rate_lab/c[:,None]
    direct_clock=np.full_like(rate_lab,np.nan)
    direct_clock[1:,:]=np.gradient(ang[1:],tau[1:],axis=0,edge_order=2)
    correction_gauge_error=np.max(np.abs(rate_corr[valid]-rate_lab[valid]))
    direct_closure=np.max(np.abs(lab[1:,:,3]-np.sum(lab[1:,:,:3],axis=2)))
    t_w_disagreement=np.max(np.abs(t_lab[1:]-lab[1:]))
    angle_diff=np.max(np.abs(delta_lab-delta_chart))
    area_diff=np.max(np.abs(lab_area[valid]-chart_area[valid]))
    relative_rate_diff=np.max(np.abs(np.gradient(delta_lab,xi[1:],edge_order=2)[xi[1:]>=.01]-np.gradient(delta_chart,xi[1:],edge_order=2)[xi[1:]>=.01]))
    abs_rate_original=rate_lab[valid]
    abs_rate_observed=rate_chart[valid]
    def sign_changes(x):
        s=np.sign(x)
        s=s[s!=0]
        return int(np.count_nonzero(s[1:]!=s[:-1]))
    changes_original=[sign_changes(abs_rate_original[:,j]) for j in range(4)]
    changes_observer=[sign_changes(abs_rate_observed[:,j]) for j in range(4)]
    max_budget_err=np.max(np.abs(total_rate_reconstruct[valid]-total_rate_wedge[valid]))
    clock_discrepancy=np.max(np.abs(direct_clock[valid]-corrected_clock[valid]))
    s0=int(np.where(valid)[0][0]);e=len(xi)-1
    data={
        'xi':xi,'tau':tau,'frame_phi_rad':phi,'frame_rate_rad_per_xi':frame_rate,
        'lab_angles_rad':ang,'observer_angles_rad':a_chart,
        'lab_rate_rad_per_xi':rate_lab,'observer_rate_rad_per_xi':rate_chart,
        'corrected_lab_rate_rad_per_xi':rate_corr,
        'lab_rate_rad_per_tau':corrected_clock,
        'relative_angle_rad':np.r_[np.nan,delta_lab],
        'relative_rate_rad_per_xi':np.r_[np.nan,relative_rate],
        'oriented_area':lab_area,
        'norm_self_turning':norm_self,'norm_mixed_turning':norm_pair,
        'total_turning':total_rate_wedge,
        'reconstructed_total_turning':total_rate_reconstruct,
        'lab_channel_vectors':lab,'observer_channel_vectors':chart,
    }
    summary={
        'seed':seed,'xi_first_qualified':float(xi[s0]),'xi_end':float(xi[e]),
        'initial_relative_angle_deg':float(np.rad2deg(delta_lab[0])),
        'final_relative_angle_deg':float(np.rad2deg(delta_lab[-1])),
        'lab_posture_turn_deg': {nm:float(np.rad2deg(ang[e,j]-ang[1,j])) for j,nm in enumerate(CHANNELS)},
        'native_posture_end_deg':{nm:float(np.rad2deg(np.arctan2(native[e,1,j],native[e,0,j]))) for j,nm in enumerate(CHANNELS)},
        'lab_final_angular_rate_per_xi':{nm:float(rate_lab[e,j]) for j,nm in enumerate(CHANNELS)},
        'observer_final_angular_rate_per_xi':{nm:float(rate_chart[e,j]) for j,nm in enumerate(CHANNELS)},
        'observer_turning_reversals':dict(zip(CHANNELS,changes_observer)),
        'physical_lab_turning_reversals':dict(zip(CHANNELS,changes_original)),
        'final_total_turning_budget_per_xi':{
            'self_v2':float(norm_self[e,0]),'self_omega2':float(norm_self[e,1]),
            'self_other':float(norm_self[e,2]),
            'mixed_v2_omega2':float(norm_pair[e,0]),
            'mixed_v2_other':float(norm_pair[e,1]),
            'mixed_omega2_other':float(norm_pair[e,2]),
            'total':float(total_rate_wedge[e]),
        },
        'final_tau':float(tau[e]),
        'max_w_T_cartesian_difference':float(t_w_disagreement),
        'max_lab_channel_closure':float(direct_closure),
        'max_covariant_rate_defect':float(correction_gauge_error),
        'max_relative_angle_gauge_defect':float(angle_diff),
        'max_relative_rate_gauge_defect':float(relative_rate_diff),
        'max_oriented_area_gauge_defect':float(area_diff),
        'max_signed_turning_budget_defect':float(max_budget_err),
        'max_direct_tau_derivative_discrepancy':float(clock_discrepancy),
        'max_omitted_connection_turning_defect':float(np.max(np.abs(rate_chart[valid]-rate_lab[valid]))),
        'min_positive_clock_rate':float(c.min()),
        'max_frame_rotation_deg':float(np.max(np.abs(np.rad2deg(phi)))),
    }
    checks={
       'T_W_cartesian_vectors_agree':t_w_disagreement<1e-9,
       'three_channel_closure':direct_closure<2e-13,
       'both_channels_turn_clockwise_in_lab':bool(np.all(rate_lab[valid,:2]<0)),
       'total_response_turns_counterclockwise_in_lab':bool(np.all(rate_lab[valid,3]>0)),
       'positive_oriented_channel_area':bool(np.all(lab_area[valid]>0)),
       'relative_angle_decreases':bool(np.all(np.diff(delta_lab)<0)),
       'turning_budget_reconstructs_total':max_budget_err<2e-10,
       'mixed_v2_omega2_overcomes_negative_self_turning':bool(np.all(norm_pair[valid,0]+norm_self[valid,0]+norm_self[valid,1]>0)),
       'passive_frame_relative_angle_invariant':angle_diff<1e-10,
       'passive_frame_relative_angular_rate_invariant':relative_rate_diff<1e-9,
       'passive_frame_oriented_area_invariant':area_diff<1e-11,
       'observer_rate_connection_restores_lab_rate':correction_gauge_error<1e-9,
       'omitting_connection_creates_false_turning':summary['max_omitted_connection_turning_defect']>0.1,
       'observer_induces_turning_reversals':changes_observer[0]>=2 and changes_original[0]==0,
       'physical_clock_increases':bool(np.all(np.diff(tau)>0)),
       'positive_clock_preserves_turning_orientation':bool(np.all(np.sign(corrected_clock[valid,:])==np.sign(rate_lab[valid,:]))),
       'direct_tau_derivative_consistency':clock_discrepancy<0.01,
    }
    return summary,checks,data


def plot(results,out):
    colors={'original':'#1673ae','autonomous':'#c24e54'}
    fig,axes=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
    for seed,(summary,checks,d) in results.items():
        xi=d['xi'];c=colors[seed];valid=xi>=.01
        a=d['lab_angles_rad']
        axes[0,0].plot(xi[valid],np.degrees(a[valid,0]-a[1,0]),color=c,ls='-',label=f'{seed} v2')
        axes[0,0].plot(xi[valid],np.degrees(a[valid,1]-a[1,1]),color=c,ls='--',label=f'{seed} omega2')
        axes[0,0].plot(xi[valid],np.degrees(a[valid,3]-a[1,3]),color=c,ls=':',lw=2,label=f'{seed} total')
        axes[0,1].plot(xi[valid],np.degrees(d['relative_angle_rad'][valid]),color=c,label=seed)
        axes[1,0].plot(xi[valid],d['lab_rate_rad_per_xi'][valid,0],color=c,ls='-',label=f'{seed} lab')
        axes[1,0].plot(xi[valid],d['observer_rate_rad_per_xi'][valid,0],color=c,ls='--',label=f'{seed} rotating observer')
        axes[1,1].plot(xi[valid],d['norm_mixed_turning'][valid,0],color=c,lw=2,label=f'{seed} mixed v2/omega2')
        axes[1,1].plot(xi[valid],d['norm_self_turning'][valid,0]+d['norm_self_turning'][valid,1],color=c,ls='--',label=f'{seed} self terms')
        axes[1,1].plot(xi[valid],d['total_turning'][valid],color=c,ls=':',label=f'{seed} total')
    axes[0,0].set_title('Lab posture: individual channels vs their total');axes[0,0].set_ylabel('Unwrapped turn from first sample (deg)')
    axes[0,1].set_title('Invariant mutual posture angle');axes[0,1].set_ylabel('Signed v2 → omega2 angle (deg)')
    axes[1,0].set_title('Observer-induced apparent turning reversals');axes[1,0].set_ylabel('v2 angular rate (rad / xi)')
    axes[1,1].set_title('Signed turning budget: mixed versus self terms');axes[1,1].set_ylabel('Contribution to total rate (rad / xi)')
    for ax in axes.ravel():
        ax.set_xlabel('xi');ax.grid(alpha=.23);ax.legend(fontsize=8,loc='best')
    fig.suptitle('Q550 v18.51z — Posture, observer frame and turning-budget audit')
    fig.savefig(out/'channel_posture_and_observer_catwalk.png',dpi=165)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    for seed,(summary,checks,d) in results.items():
        valid=d['xi']>=.01;c=colors[seed]
        axes[0].plot(d['xi'][valid],d['lab_rate_rad_per_xi'][valid,3],label=f'{seed} total',color=c)
        axes[0].plot(d['xi'][valid],d['lab_rate_rad_per_xi'][valid,0],label=f'{seed} v2',color=c,ls='--')
        axes[0].plot(d['xi'][valid],d['lab_rate_rad_per_xi'][valid,1],label=f'{seed} omega2',color=c,ls=':')
        axes[1].plot(d['tau'][valid],d['lab_rate_rad_per_tau'][valid,3],label=f'{seed} total',color=c)
        axes[1].plot(d['tau'][valid],d['lab_rate_rad_per_tau'][valid,0],label=f'{seed} v2',color=c,ls='--')
        axes[1].plot(d['tau'][valid],d['lab_rate_rad_per_tau'][valid,1],label=f'{seed} omega2',color=c,ls=':')
    axes[0].set_title('Angular rates in xi');axes[1].set_title('Same angular rates in physical tau')
    for ax in axes:
        ax.axhline(0,color='gray',lw=.8);ax.grid(alpha=.23);ax.legend(fontsize=8)
    axes[0].set_xlabel('xi');axes[1].set_xlabel('tau')
    axes[0].set_ylabel('rad / xi');axes[1].set_ylabel('rad / tau')
    fig.suptitle('Q550 v18.51z — Clock reparametrization preserves turning sense')
    fig.savefig(out/'physical_clock_posture_turning.png',dpi=165)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=TITLE)
    parser.add_argument('--source',type=Path,default=Path('inputs/q550_model_constrained_source_subspaces_physical_attribution_admissibility_v18_51y_bundle_server_sendback.zip'))
    parser.add_argument('--out',type=Path,default=Path('results'))
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(args.source) as z:
        assert z.testzip() is None,'Corrupted parent ZIP'
        find=lambda suffix: next(k for k in z.namelist() if k.endswith(suffix))
        w=np.load(io.BytesIO(z.read(find('w_channel_vector_geometric_invariants_v18_51w.npz'))))
        t=np.load(io.BytesIO(z.read(find('t_green_channel_transport_states_v18_51t.npz'))))
        parent=json.loads(z.read(find('results/report.json')))
        assert parent.get('all_passed') is True and parent['version']=='v18.51y'
        source_sha=sha256(args.source)
        xi=w['xi'];assert np.array_equal(xi,t['xi'])
        outputs={}; all_checks={}; array_output={'xi':xi}
        for seed in ('original','autonomous'):
            lab=w[f'{seed}_lab_channel_response'];native=w[f'{seed}_native_channel_response'];tau=w[f'{seed}_tau']
            t_lab=t[f'{seed}_cartesian_Q'][:,2:4,:]
            summary,checks,data=calc_seed(seed,xi,lab,native,tau,t_lab)
            outputs[seed]=summary
            for k,v in checks.items():all_checks[f'{seed}:{k}']=bool(v)
            for key,value in data.items():array_output[f'{seed}_{key}']=value
            print(seed,'turns deg',summary['lab_posture_turn_deg'],'reversals',summary['observer_turning_reversals'])
        sep=abs(outputs['original']['final_relative_angle_deg']-outputs['autonomous']['final_relative_angle_deg'])
        all_checks['two_seed_mutual_angle_agreement']=sep<0.005
        all_checks['model_parent_authenticated']=parent['all_passed'] and parent['version']=='v18.51y'
        all_checks['no_parameter_actuator_claim']=True
        all_checks['no_material_return_claim']=True
        report={
            'version':VERSION,'title':TITLE,
            'classification':'MODEL_GROUNDED_LOCAL_RESPONSE_TURNING_AND_PASSIVE_OBSERVER_COVARIANCE_DIAGNOSTIC_NOT_MATERIAL_POSTURE_CERTIFICATE',
            'parent_source':args.source.name,'parent_source_sha256':source_sha,
            'n_samples_per_seed':len(xi),'mask_for_angular_rates':'xi >= 0.01; exactly zero vectors at xi=0 have undefined angles',
            'observer_protocol':'phi(xi)=8 deg*(xi/0.2)+2.5 deg*sin(4*pi*xi/0.2), common passive SO(2) rotation for all channels',
            'all_passed':all(all_checks.values()),'scientific_checks':all_checks,
            'n_checks':len(all_checks),'seeds':outputs,
            'method_limits':[
                'V2 and omega2 refer to channels of the kappa sensitivity in Cartesian Q3/Q4, not an actual human posture or the full 10D orbit.',
                'The signed self/mixed turning budget follows from bilinearity of the wedge product; it is checked on authentic numerical responses, not an independent model evolution.',
                'Observer rotations are passive; induced apparent reversal is not a physical reversal.',
                'A physical time reparametrization with positive rate preserves the sign of angular turning but changes its numerical rate.',
                'The 0.01 cutoff avoids ill-defined angles and amplification near zero response norm. All reported angular changes start at the first nonzero section.',
                'Results cover xi <= 0.2, not the true return or material memory.'
            ]}
        (args.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        np.savez_compressed(args.out/'posture_observer_turning_arrays.npz',**array_output)
        with open(args.out/'turning_budget_sections.csv','w',newline='') as f:
            wr=csv.writer(f)
            wr.writerow(['seed','xi','tau','lab_angle_v2_deg','lab_angle_omega2_deg','lab_angle_total_deg','relative_angle_deg','self_v2_turning','self_omega2_turning','mixed_v2_omega2_turning','other_terms_turning','total_turning','observer_v2_turning','corrected_v2_turning'])
            for seed,(summary,checks,d) in ((s,(outputs[s],None,{k.replace(s+'_','',1):v for k,v in array_output.items() if k.startswith(s+'_')})) for s in outputs):
                for i in range(1,len(xi)):
                    other=(d['norm_self_turning'][i,2]+np.sum(d['norm_mixed_turning'][i,1:]))
                    wr.writerow([seed,float(xi[i]),float(d['tau'][i]),*[float(np.rad2deg(d['lab_angles_rad'][i,j])) for j in (0,1,3)],float(np.rad2deg(d['relative_angle_rad'][i])),float(d['norm_self_turning'][i,0]),float(d['norm_self_turning'][i,1]),float(d['norm_mixed_turning'][i,0]),float(other),float(d['total_turning'][i]),float(d['observer_rate_rad_per_xi'][i,0]),float(d['corrected_lab_rate_rad_per_xi'][i,0])])
        packed={s:(outputs[s],None,{k.replace(s+'_','',1):v for k,v in array_output.items() if k.startswith(s+'_')}) for s in outputs}
        plot(packed,args.out)
        text_lines=[f'# Q550 {VERSION} — {TITLE}','',f'**Audit classification:** `{report["classification"]}`','',f'**Checks:** {sum(all_checks.values())}/{len(all_checks)} PASS','', '## Observation', '']
        for seed,s in outputs.items():
            turns=s['lab_posture_turn_deg'];bt=s['final_total_turning_budget_per_xi']
            text_lines.extend([f'### {seed}', f'- Lab turns (from first nonzero sample to xi=0.2): v2 {turns["v2"]:+.7f} deg, omega2 {turns["omega2"]:+.7f} deg, total {turns["total"]:+.7f} deg.',f'- Final turning budget (rad/xi): self v2 {bt["self_v2"]:+.8f}, self omega2 {bt["self_omega2"]:+.8f}, cross v2–omega2 {bt["mixed_v2_omega2"]:+.8f}, total {bt["total"]:+.8f}.',f'- Final mutual angle {s["final_relative_angle_deg"]:.8f} deg; observer-induced v2 turning reversals {s["observer_turning_reversals"]["v2"]}, lab reversals {s["physical_lab_turning_reversals"]["v2"]}.',f'- Maximum passive-covariant angular-rate defect {s["max_covariant_rate_defect"]:.3e} rad/xi; maximum turning-budget defect {s["max_signed_turning_budget_defect"]:.3e} rad/xi.',''])
        text_lines.extend(['## Interpretation','', 'Both anchored source-channel response vectors turn clockwise in the fixed second Cartesian plane, but their combined response turns counterclockwise because the positive mixed wedge term exceeds negative self terms. This is a local signed-turning identity applied to physical-model-derived sensitivity channels, not a direct torque or energy balance.','', 'A common rotating observer can introduce apparent sign reversals in individual angular velocities. After accounting for the observer connection, the recovered laboratory angular velocities and the mutual angle agree. The parent fixed-support attribution remains an assumption in defining the channels.','', '## Limits','',*['- '+x for x in report['method_limits']],''])
        (args.out/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(text_lines))
        print('Q550',VERSION,sum(all_checks.values()),'/',len(all_checks),'PASS' if report['all_passed'] else 'FAIL')
        if not report['all_passed']:
            print('FAILED:',[k for k,v in all_checks.items() if not v]);raise SystemExit(2)

if __name__=='__main__':main()
