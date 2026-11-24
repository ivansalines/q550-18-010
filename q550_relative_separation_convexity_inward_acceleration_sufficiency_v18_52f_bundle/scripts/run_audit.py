#!/usr/bin/env python3
"""Q550 v18.52f: sampled separation convexity and inward acceleration sufficiency."""
import argparse
import csv
import io
import json
import os
from pathlib import Path
import zipfile

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PARENT = ROOT / 'inputs' / 'q550_lateral_curvature_relative_acceleration_work_restoring_tendency_v18_52e_bundle_server_sendback.zip'


def analyze_vectors(t, q, v, a):
    """Recompute scalar observables without trusting archived scalars."""
    t, q, v, a = map(np.asarray, (t, q, v, a))
    assert t.ndim == 1 and q.shape == v.shape == a.shape == (len(t), 5)
    assert len(t) >= 16 and np.all(np.diff(t) > 0)
    assert np.isfinite(q).all() and np.isfinite(v).all() and np.isfinite(a).all()
    pos = np.einsum('ij,ij->i', v, v)
    inward = -np.einsum('ij,ij->i', q, a)
    rate = np.einsum('ij,ij->i', q, v)
    D = 0.5 * np.einsum('ij,ij->i', q, q)
    C = pos - inward
    R = np.divide(inward, pos, out=np.full_like(pos, np.nan), where=pos > 0)
    F = np.divide(pos, inward, out=np.full_like(pos, np.nan), where=inward > 0)
    return dict(tau=t,positive=pos,inward=inward,rate=rate,separation=D,
                convexity=C,ratio=R,needed_factor=F)


def sample_summary(name, x):
    t=x['tau'];P=x['positive'];I=x['inward'];C=x['convexity'];G=x['rate'];R=x['ratio'];F=x['needed_factor']
    qualified = P>0
    after_initial = np.arange(len(t))>0
    ffinite = F[np.isfinite(F)]
    cs=CubicSpline(t,G)
    # Exclude endpoint derivatives to avoid introducing endpoint-specific differentiation bias.
    deriv_rel=float(np.max(np.abs(cs.derivative()(t)[2:-2]-C[2:-2]))/max(np.max(np.abs(C)),1e-300))
    C_integral=float(trapezoid(C,t)); Gchange=float(G[-1]-G[0]);
    D_integral=float(trapezoid(G,t)); Dchange=float(x['separation'][-1]-x['separation'][0])
    def first_idx(mask):
        inds=np.flatnonzero(mask)
        return int(inds[0]) if len(inds) else None
    # PCHIP is an interpolant-only diagnostic, NOT a rigorous bound of the true flow.
    dense=np.linspace(t[0],t[-1],4*(len(t)-1)+1)
    pchip_C=PchipInterpolator(t,C)(dense)
    pchip_G=PchipInterpolator(t,G)(dense)
    idxc=first_idx(qualified & (C<=0))
    idxg=first_idx(after_initial & (G<=0))
    return dict(
      label=name,n_samples=int(len(t)),tau_end=float(t[-1]),
      ratio_initial=float(R[0]),ratio_final=float(R[-1]),ratio_max_sampled=float(np.nanmax(R)),
      ratio_min_sampled=float(np.nanmin(R)),
      min_convexity=float(np.min(C)),min_normalized_convexity=float(np.min(C[qualified]/P[qualified])),
      final_normalized_convexity=float(C[-1]/P[-1]),
      final_sufficiency_factor_shortfall=float(F[-1]),
      min_positive_needed_factor=float(np.min(ffinite)) if len(ffinite) else None,
      first_sampled_convexity_reversal_tau=float(t[idxc]) if idxc is not None else None,
      first_sampled_separation_rate_reversal_tau=float(t[idxg]) if idxg is not None else None,
      fraction_sampled_strict_convexity=float(np.mean(C>0)),
      fraction_sampled_outward_rate_after_initial=float(np.mean(G[1:]>0)),
      first_sampled_ratio_above_25pct_tau=float(t[first_idx(R>=0.25)]) if np.any(R>=0.25) else None,
      first_sampled_ratio_above_50pct_tau=float(t[first_idx(R>=0.50)]) if np.any(R>=0.50) else None,
      inward_ratio_non_decreasing_on_samples=bool(np.all(np.diff(R)>=-1e-11)),
      sampled_convexity_all_positive=bool(np.all(C>0)),
      sampled_position_rate_all_positive_after_initial=bool(np.all(G[1:]>0)),
      pchip_interpolant_min_convexity=float(np.min(pchip_C)),
      pchip_interpolant_min_rate_after_initial=float(np.min(pchip_G[1:])),
      spline_rate_derivative_vs_C_rel=float(deriv_rel),
      integral_C_predicts_rate_abs=float(abs(C_integral-Gchange)),
      integral_C_predicts_rate_rel=float(abs(C_integral-Gchange)/max(abs(Gchange),1e-300)),
      integral_rate_predicts_D_abs=float(abs(D_integral-Dchange)),
      integral_rate_predicts_D_rel=float(abs(D_integral-Dchange)/max(abs(Dchange),1e-300)),
      scientific_classification=(
         'SAMPLED_CONVEXITY_REVERSAL_PRESENT' if idxc is not None else
         'SAMPLED_INWARD_ACCELERATION_INSUFFICIENT_FOR_CONVEXITY_REVERSAL'
      ),
      continuous_time_flow_certified=False,
      physical_position_return_certified=False,
    )


def load_archive(path):
    with zipfile.ZipFile(path) as z:
        bad=z.testzip()
        if bad is not None: raise RuntimeError('Corrupt parent ZIP member: '+bad)
        matches=[n for n in z.namelist() if n.endswith('/results/physical_relative_vector_diagnostics.npz')]
        if len(matches)!=1: raise RuntimeError('Expected exactly one physical vectors NPZ in parent sendback')
        # Validate the parent archive's own source/output manifest entries.
        import hashlib
        for manifest_name in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            mnames=[n for n in z.namelist() if n.endswith('/'+manifest_name)]
            if len(mnames)!=1:raise RuntimeError('Missing or duplicate parent '+manifest_name)
            data=json.loads(z.read(mnames[0]))
            prefix=mnames[0].rsplit('/',1)[0]+'/'
            for rel,expected in data['file_sha256'].items():
                full=prefix+rel
                got=hashlib.sha256(z.read(full)).hexdigest()
                if got!=expected:raise RuntimeError('Parent manifest hash mismatch '+rel)
        dat=np.load(io.BytesIO(z.read(matches[0])),allow_pickle=False)
        return {k:np.asarray(dat[k]) for k in dat.files}


def plot_all(series, summaries, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors={'original':'C0','autonomous':'C1'}
    fig,axs=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
    for label,x in series.items():
        seed=label.split('_')[0]; alpha=.45 if '3.0e-04' in label else .75
        style='-' if '_positive_' in label else '--'
        t=x['tau']; r=x['ratio']
        axs[0,0].plot(t,100*r,color=colors[seed],ls=style,alpha=alpha,lw=1.35)
        axs[0,1].plot(t,100*x['convexity']/x['positive'],color=colors[seed],ls=style,alpha=alpha,lw=1.35)
        axs[1,0].plot(t,x['rate']/x['positive'][0],color=colors[seed],ls=style,alpha=alpha,lw=1.35)
        axs[1,1].plot(t,x['needed_factor'],color=colors[seed],ls=style,alpha=alpha,lw=1.35)
    axs[0,0].axhline(100,color='C3',lw=1.0,ls=':',label='Convexity reversal (100%)')
    axs[0,0].set(title='Inward-to-kinematic ratio',ylabel='100 × (-dQ·da)/|dV|² (%)',xlabel='Material time τ')
    axs[0,1].axhline(0,color='C3',ls=':',lw=1.0)
    axs[0,1].set(title='Normalized separation convexity',ylabel='D" / |dV|² (%)',xlabel='Material time τ')
    axs[1,0].axhline(0,color='C3',ls=':',lw=1.0)
    axs[1,0].set(title='Separation rate remains outward',ylabel='(dQ·dV) / |dV(0)|² (τ units)',xlabel='Material time τ')
    axs[1,1].set_ylim(bottom=0,top=30)
    axs[1,1].set(title='Frozen-state inward amplification needed',ylabel='|dV|² / (-dQ·da), clipped at 30',xlabel='Material time τ')
    for ax in axs.flat:
       ax.grid(alpha=.18)
       ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    from matplotlib.lines import Line2D
    fig.legend([Line2D([0],[0],color=colors[s],lw=2) for s in colors],list(colors),loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.5,1.035))
    fig.suptitle('Q550 v18.52f | Separation convexity and inward-acceleration sufficiency',fontweight='bold')
    target=out_dir/'relative_separation_convexity_inward_sufficiency.png'
    fig.savefig(target,dpi=170,bbox_inches='tight')
    plt.close(fig)
    return target


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--parent',type=Path,default=DEFAULT_PARENT)
    args=ap.parse_args()
    raw=load_archive(args.parent)
    labels=sorted({k.split('__')[0] for k in raw})
    if len(labels)!=8:raise RuntimeError('Need all eight authentic physical kick traces')
    out=ROOT/'results';out.mkdir(exist_ok=True)
    summaries=[];series={};rows=[];checks=[]
    for label in labels:
        def field(key):return raw[f'{label}__{key}']
        x=analyze_vectors(field('tau'),field('dQ'),field('dV'),field('da'))
        s=sample_summary(label,x)
        for our,parent in [('rate','relative_position_rate'),('convexity','second_position_separation_derivative')]:
            diff=float(np.max(np.abs(x[our]-field(parent))))
            tol=1e-13*max(1.0,np.max(np.abs(x[our])))
            checks.append((f'{label}.{our}_matches_parent',diff<=tol,diff))
        diff=float(np.max(np.abs(-x['inward']-field('position_acceleration_dot'))))
        checks.append((f'{label}.inward_matches_parent',diff<=1e-13,diff))
        checks.append((f'{label}.spline_derivative_closure',s['spline_rate_derivative_vs_C_rel']<5e-6,s['spline_rate_derivative_vs_C_rel']))
        checks.append((f'{label}.quadrature_rate_closure',s['integral_C_predicts_rate_rel']<1e-4,s['integral_C_predicts_rate_rel']))
        checks.append((f'{label}.quadrature_distance_closure',s['integral_rate_predicts_D_rel']<1e-4,s['integral_rate_predicts_D_rel']))
        summaries.append(s);series[label]=x
        for i,t in enumerate(x['tau']):
           rows.append(dict(label=label,index=i,tau=t,relative_speed_term=x['positive'][i],inward_acceleration_term=x['inward'][i],separation_convexity=x['convexity'][i],inward_ratio=x['ratio'][i],required_inward_amplification=x['needed_factor'][i],position_separation_rate=x['rate'][i],position_distance_squared_half=x['separation'][i]))
    with (out/'separation_convexity_ledger.csv').open('w',newline='') as fp:
       w=csv.DictWriter(fp,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    image=plot_all(series,summaries,out)
    pass_count=sum(k[1] for k in checks)
    report=dict(version='v18.52f',title='Relative-Separation Convexity & Inward-Acceleration Sufficiency Audit',
          stage='AUTHENTIC_PARENT_VECTOR_DIAGNOSTIC_NO_NEW_NONLINEAR_INTEGRATIONS',
          parent_archive=args.parent.name,all_numerical_checks_passed=pass_count==len(checks),
          check_count=len(checks),checks_passed=pass_count,
          checks=[{'name':name,'passed':bool(ok),'residual':float(res)} for name,ok,res in checks],
          physical_kicks=summaries,
          global_summary=dict(kick_count=len(summaries),samples_total=len(rows),
            maximum_inward_fraction=float(max(s['ratio_max_sampled'] for s in summaries)),
            terminal_inward_ratio_range=[float(min(s['ratio_final'] for s in summaries)),float(max(s['ratio_final'] for s in summaries))],
            terminal_required_inward_factor_range=[float(min(s['final_sufficiency_factor_shortfall'] for s in summaries)),float(max(s['final_sufficiency_factor_shortfall'] for s in summaries))],
            min_normalized_convexity=float(min(s['min_normalized_convexity'] for s in summaries)),
            sampled_convexity_reversals=sum(s['first_sampled_convexity_reversal_tau'] is not None for s in summaries),
            sampled_position_rate_reversals=sum(s['first_sampled_separation_rate_reversal_tau'] is not None for s in summaries),
            scientific_classification=('NO_SAMPLED_CONVEXITY_REVERSAL_OR_POSITION_RETURN' if all(s['sampled_convexity_all_positive'] and s['sampled_position_rate_all_positive_after_initial'] for s in summaries) else 'MIXED_SAMPLED_SIGNS_REQUIRING_EVENT_REVIEW'),
            certified_continuous_bound=False,physical_recovery_certified=False),
          caveats=['True ODE behavior between samples is not certified by PCHIP or cubic spline.',
                   'The amplification factor is an instantaneous arithmetic ratio, not an admissible physical force change.',
                   'No extension beyond the authenticated time horizon has been performed.'])
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    g=report['global_summary']
    lines=['# Q550 v18.52f — Server replay scientific report','',f'Numerical checks: **{pass_count}/{len(checks)}**.',
     f'Authentic physical trajectories: **{len(summaries)}**, sampled sections **{len(rows)}**.',
     f'Maximum sampled inward-to-positive ratio: **{g["maximum_inward_fraction"]:.9f}** (threshold=1).',
     f'Final required inward factor range (frozen-state diagnostic only): **{g["terminal_required_inward_factor_range"][0]:.6f}–{g["terminal_required_inward_factor_range"][1]:.6f}**.',
     f'Minimum normalized convexity across all sampled states: **{g["min_normalized_convexity"]:.9f}**.',
     f'Scientific classification: **{g["scientific_classification"]}**.','',
     '| Kick | R at final | Minimum C/P | Required factor at final | First sampled C<=0 | First sampled G<=0 |',
     '|---|---:|---:|---:|---|---|']
    for s in summaries:
      lines.append(f'| {s["label"]} | {s["ratio_final"]:.7f} | {s["min_normalized_convexity"]:.7f} | {s["final_sufficiency_factor_shortfall"]:.5f} | {s["first_sampled_convexity_reversal_tau"]} | {s["first_sampled_separation_rate_reversal_tau"]} |')
    lines.extend(['','A loss of separation convexity (C<=0) would **not by itself** imply separation return (G<0).',
      'No rigorous continuous-time sign enclosure, global dynamics, or asymptotic recovery certificate is claimed.',
      'The parent nonlinear trajectories were not reintegrated in this audit.',''])
    (out/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines))
    print('Q550 v18.52f checks',pass_count,'/',len(checks))
    print('Inward max ratio',g['maximum_inward_fraction'])
    print('Terminal frozen-state factor range',g['terminal_required_inward_factor_range'])
    print('Scientific result',g['scientific_classification'])
    if pass_count!=len(checks):raise SystemExit('Numerical checks failed; see results/report.json')

if __name__=='__main__':main()
