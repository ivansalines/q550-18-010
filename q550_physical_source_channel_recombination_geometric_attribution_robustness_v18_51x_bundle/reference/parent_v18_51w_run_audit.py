#!/usr/bin/env python3
"""Q550 v18.51w: Channel-Vector Angle, Oriented-Area & Observer-Critical Compensation Audit.

Analyze signed 2D Green-channel responses archived in the verified v18.51v
server sendback. No nonlinear ODE or tangent IVP is integrated here.
Observer angles refer to the local physical radial--azimuthal frame; neither a
vanishing projection nor complete projected cancellation is physical closure.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import pathlib
import sys
import zipfile

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
VROOT = 'q550_physical_clock_green_transport_moving_frame_invariance_v18_51v_bundle'
VFILE = ROOT / 'inputs' / (VROOT + '_server_sendback.zip')
SEEDS = ('original', 'autonomous')
CHANNELS = ('v2', 'omega2', 'other', 'total')
TINY_NORM = 1e-15


def require(ok, description):
    if not bool(ok):
        raise AssertionError(description)


def maxabs(a):
    return float(np.max(np.abs(a)))


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def cos_sin_rotation(phi):
    c, s = np.cos(phi), np.sin(phi)
    r = np.empty(np.shape(c) + (2, 2), dtype=np.float64)
    r[..., 0, 0] = c
    r[..., 0, 1] = -s
    r[..., 1, 0] = s
    r[..., 1, 1] = c
    return r


def wrap_half_turn(phi):
    """Unoriented axis angle in [-pi/2, pi/2); periodic modulo pi."""
    return ((phi + np.pi/2) % np.pi) - np.pi/2


def critical_axis(v):
    """Axis ψ whose positive radial projection on v is zero.

    The projected radial component is r_ψ = cos(ψ)*v_r + sin(ψ)*v_a.
    The axis, not its arrow, is critical; ψ and ψ+π are equivalent.
    """
    return wrap_half_turn(np.arctan2(-v[..., 0], v[..., 1]))


def signed_area(a, b):
    return a[..., 0]*b[..., 1] - a[..., 1]*b[..., 0]


def projected_radial(v, psi):
    return np.cos(psi)*v[..., 0] + np.sin(psi)*v[..., 1]


def cancellation_fraction(comps):
    """One observation axis: 1-|sum of signed components|/sum absolutes.

    Input last axis is 3 independent source channels. Undefined for zero data.
    """
    den = np.sum(np.abs(comps), axis=-1)
    out = np.full(den.shape, np.nan)
    np.divide(den - np.abs(np.sum(comps, axis=-1)), den,
              out=out, where=den > 1e-20)
    return np.clip(out, 0, 1)


def read_verified_parent(full=True):
    require(VFILE.is_file(), 'Missing v18.51v authenticated parent sendback')
    with zipfile.ZipFile(VFILE) as z:
        require(z.testzip() is None, 'v18.51v ZIP CRC failure')
        report = json.loads(z.read(VROOT + '/results/report.json'))
        require(report['version'] == 'v18.51v' and report['all_passed']
                and all(report['scientific_checks'].values()), 'Parent v18.51v did not pass')
        source = json.loads(z.read(VROOT + '/SOURCE_MANIFEST.json'))
        outputs = json.loads(z.read(VROOT + '/SHA256SUMS_OUTPUT.json'))
        # Verify all source AND result hashes, including immutable T/R/U genealogy.
        for rel, expected in source['files'].items():
            require(sha256(z.read(VROOT + '/' + rel)) == expected, 'V source SHA256 mismatch: '+rel)
        for rel, expected in outputs['file_sha256'].items():
            require(sha256(z.read(VROOT + '/' + rel)) == expected, 'V output SHA256 mismatch: '+rel)
        if not full:
            return None, None, {'v_archive_sha256': sha256(VFILE.read_bytes()),
                'v_source_hashes': len(source['files']), 'v_output_hashes': len(outputs['file_sha256'])}
        vraw = z.read(VROOT + '/results/physical_clock_frame_green_tensors.npz')
        traw = z.read(VROOT + '/inputs/authenticated_green_channel_transport_states_v18_51t.npz')
        provenance = {'v_archive_sha256': sha256(VFILE.read_bytes()),
                      'v_source_hashes': len(source['files']), 'v_output_hashes': len(outputs['file_sha256']),
                      'v_version': report['version'], 'v_scientific_checks_pass': True}
    with np.load(io.BytesIO(vraw), allow_pickle=False) as v:
        V = {key: v[key].copy() for key in v.files if not key.endswith('_kernel') and not key.endswith('_density')}
    with np.load(io.BytesIO(traw), allow_pickle=False) as t:
        T = {key: t[key].copy() for key in t.files if key == 'xi' or key.endswith(('radial2', 'azimuthal2', 'cartesian_Q'))}
    return V, T, provenance


def measure_seed(seed, xi, V, T, ledger, output_arrays):
    tau = V[seed + '_tau']
    natural = np.stack((T[seed + '_radial2'], T[seed + '_azimuthal2']), axis=1)
    lab = V[seed + '_lab_frame_response']
    cartesian = T[seed + '_cartesian_Q'][:, 2:4, :]
    require(natural.shape == (len(xi), 2, 4), 'Native response array shape error')
    require(lab.shape == natural.shape, 'Lab response array shape error')
    require(cartesian.shape == natural.shape, 'T Cartesian response array shape error')
    require(np.all(np.isfinite(natural)) and np.all(np.isfinite(lab)), 'Nonfinite response')
    require(np.all(np.diff(tau) > 0), 'Physical clock must be strictly monotone')
    require(maxabs(T['xi'] - xi) < 1e-14, 'Nonmatching v and t xi sections')
    n = len(xi)
    a, b, other, total = (natural[:, :, k] for k in range(4))
    la, lb = lab[:, :, 0], lab[:, :, 1]
    norms = np.linalg.norm(a, axis=1)*np.linalg.norm(b, axis=1)
    valid = norms > TINY_NORM
    dot = np.sum(a*b, axis=1)
    det = signed_area(a, b)
    angles_rad = np.full(n, np.nan)
    angles_rad[valid] = np.arctan2(det[valid], dot[valid])
    angles_deg = np.degrees(angles_rad)
    dot_normalized = np.full(n, np.nan)
    area_normalized = np.full(n, np.nan)
    np.divide(dot, norms, out=dot_normalized, where=valid)
    np.divide(det, norms, out=area_normalized, where=valid)
    pa, pb, pt = [np.full(n, np.nan) for _ in range(3)]
    for arr, v in ((pa, a), (pb, b), (pt, total)):
        mask = valid & (np.linalg.norm(v, axis=1) > 1e-20)
        arr[mask] = critical_axis(v[mask])
    # Angles at an event are label-invariant under xi -> tau; numerical theta and
    # geometric det are unchanged. This is not a fixed-tau kappa derivative.
    psi_sine = 0.38*np.sin(2*np.pi*xi/xi[-1]) + 0.27*xi/xi[-1]
    psi_crest = 0.31*np.cos(3*np.pi*xi/xi[-1]) + 0.17*xi/xi[-1] + 0.08
    frame_errors = {}
    for label, psi in (('sine', psi_sine), ('crest', psi_crest)):
        new = np.einsum('iab,ibc->iac', cos_sin_rotation(-psi), natural)
        rotated_a, rotated_b = new[:, :, 0], new[:, :, 1]
        dotnew = np.sum(rotated_a*rotated_b, axis=1)
        detnew = signed_area(rotated_a, rotated_b)
        p_b_new = np.full(n, np.nan)
        p_b_new[valid] = critical_axis(rotated_b[valid])
        # Axis differences must be compared modulo pi (not 2pi).
        axis_defect = maxabs(np.sin(p_b_new[valid] - (pb[valid]-psi[valid])))
        reconstruction = np.einsum('iab,ibc->iac', cos_sin_rotation(psi), new)
        frame_errors[label] = {'max_dot_invariance_defect': maxabs(dotnew-dot),
                              'max_oriented_area_invariance_defect': maxabs(detnew-det),
                              'max_channel_vector_reconstruction_defect': maxabs(reconstruction-natural),
                              'max_critical_axis_covariance_sine_defect': axis_defect}
    mirror = natural.copy()
    mirror[:, 1, :] *= -1
    reflection_error = maxabs(signed_area(mirror[:, :, 0], mirror[:, :, 1]) + det)
    native_dot_vs_lab = maxabs(dot - np.sum(la*lb, axis=1))
    native_area_vs_lab = maxabs(det - signed_area(la,lb))
    cartesian_agreement = maxabs(lab - cartesian)
    additivity = maxabs(a+b+other-total)
    recovered_norm = maxabs(np.linalg.norm(total,axis=1)-np.linalg.norm(lab[:,:,3],axis=1))
    unit_identity = maxabs((dot_normalized[valid]**2+area_normalized[valid]**2)-1)
    # Channel projection zero at predicted axis. A finite absolute tolerance
    # prevents misleading enormous *relative* errors near the initial null.
    critical_a_error = maxabs(projected_radial(a[valid],pa[valid]))
    critical_b_error = maxabs(projected_radial(b[valid],pb[valid]))
    total_valid = np.isfinite(pt)
    critical_total_error = maxabs(projected_radial(total[total_valid],pt[total_valid]))
    c_natural = cancellation_fraction(natural[:,0,:3])
    # Search the observer-angle range at the endpoint without disguising the
    # geometric 180-degree periodicity of unoriented projection axes.
    scan_degrees = np.linspace(-90,90,721)
    scan_psi = np.deg2rad(scan_degrees)
    end_source = natural[-1,:,:3]
    projected = np.einsum('a,ac->c',np.array([1.,0.]),end_source)  # exact radial check
    require(maxabs(projected-end_source[0]) < 1e-15,'Natural radial projection incorrect')
    scan_channels = np.cos(scan_psi[:,None])*end_source[0,None,:] + np.sin(scan_psi[:,None])*end_source[1,None,:]
    c_scan = cancellation_fraction(scan_channels)
    # ψ_b may be computed from B and checked against its roots on the grid.
    pb_last = float(pb[-1]); pt_last = float(pt[-1]);pa_last=float(pa[-1])
    psi_sine_last = float(psi_sine[-1]);psi_crest_last = float(psi_crest[-1])
    c_sine_last = float(cancellation_fraction((np.cos(psi_sine_last)*end_source[0] +
                                              np.sin(psi_sine_last)*end_source[1])[None,:])[0])
    c_crest_last = float(cancellation_fraction((np.cos(psi_crest_last)*end_source[0] +
                                               np.sin(psi_crest_last)*end_source[1])[None,:])[0])
    # End geometry characterizes mutual orientation; do not attribute it to a
    # material constitutive law, because source-channel splitting is coordinate-dependent.
    measurements = {
        'seed': seed, 'n_sections': n, 'n_valid_angles': int(valid.sum()),
        'first_angle_defined_xi': float(xi[np.where(valid)[0][0]]),
        'angle_end_deg': float(angles_deg[-1]),
        'angle_min_deg': float(np.nanmin(angles_deg)),
        'angle_max_deg': float(np.nanmax(angles_deg)),
        'angle_first_nonzero_deg': float(angles_deg[valid][0]),
        'angle_all_positive_orientation': bool(np.all(det[valid]>0)),
        'dot_all_positive': bool(np.all(dot[valid]>0)),
        'oriented_area_end': float(det[-1]),
        'oriented_area_normalized_end': float(area_normalized[-1]),
        'dot_end': float(dot[-1]),
        'dot_normalized_end': float(dot_normalized[-1]),
        'v2_end_norm':float(np.linalg.norm(a[-1])),
        'omega2_end_norm':float(np.linalg.norm(b[-1])),
        'total_end_norm':float(np.linalg.norm(total[-1])),
        'critical_observer_omega2_end_deg':float(np.degrees(pb_last)),
        'critical_observer_v2_end_deg':float(np.degrees(pa_last)),
        'critical_observer_total_end_deg':float(np.degrees(pt_last)),
        'natural_radial_cancellation_end_percent':float(100*c_natural[-1]),
        'sine_radial_cancellation_end_percent':float(100*c_sine_last),
        'crest_radial_cancellation_end_percent':float(100*c_crest_last),
        'max_invariant_dot_native_vs_lab': native_dot_vs_lab,
        'max_invariant_area_native_vs_lab': native_area_vs_lab,
        'max_lab_vs_independent_t_cartesian': cartesian_agreement,
        'max_source_additivity_defect': additivity,
        'max_vector_norm_native_vs_lab': recovered_norm,
        'max_unit_circle_identity_defect': unit_identity,
        'max_mirror_oriented_area_sign_defect': reflection_error,
        'max_v2_projection_zero_error': critical_a_error,
        'max_omega2_projection_zero_error': critical_b_error,
        'max_total_projection_zero_error': critical_total_error,
        'observer_frame_checks': frame_errors,
        'critical_axis_interpretation': 'unoriented local radial-azimuthal axis modulo pi, not material event',
    }
    for i in range(n):
        ledger.append({
            'seed':seed, 'xi':float(xi[i]), 'tau':float(tau[i]),
            'response_v2_norm':float(np.linalg.norm(a[i])),
            'response_omega2_norm':float(np.linalg.norm(b[i])),
            'mutual_angle_degrees':float(angles_deg[i]) if valid[i] else '',
            'normalized_dot':float(dot_normalized[i]) if valid[i] else '',
            'dot':float(dot[i]),'oriented_area':float(det[i]),
            'normalized_oriented_area':float(area_normalized[i]) if valid[i] else '',
            'critical_omega2_observer_degrees':float(np.degrees(pb[i])) if np.isfinite(pb[i]) else '',
            'critical_v2_observer_degrees':float(np.degrees(pa[i])) if np.isfinite(pa[i]) else '',
            'critical_total_observer_degrees':float(np.degrees(pt[i])) if np.isfinite(pt[i]) else '',
            'radial_cancel_fraction':float(c_natural[i]) if np.isfinite(c_natural[i]) else '',
            'lab_v2_x':float(la[i,0]),'lab_v2_y':float(la[i,1]),
            'lab_omega2_x':float(lb[i,0]),'lab_omega2_y':float(lb[i,1]),
        })
    for key, arr in {
        'mutual_angle_degrees':angles_deg,'dot':dot,'oriented_area':det,
        'normalized_dot':dot_normalized,'normalized_oriented_area':area_normalized,
        'critical_omega2_observer_degrees':np.degrees(pb),
        'critical_v2_observer_degrees':np.degrees(pa),
        'critical_total_observer_degrees':np.degrees(pt),
        'native_channel_response':natural, 'lab_channel_response':lab,
        'radial_cancel_fraction':c_natural, 'tau':tau,
        'compensation_scan_degrees':scan_degrees,'compensation_scan_fraction':c_scan,
    }.items(): output_arrays[seed+'_'+key]=arr
    return measurements


def render_figures(xi, arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out = OUT
    fig, axs=plt.subplots(1,2,figsize=(13,5.3))
    for s in SEEDS:
        a=arrays[s+'_mutual_angle_degrees'];z=arrays[s+'_normalized_dot']
        axs[0].plot(xi,a,label=s)
        axs[1].plot(xi,z,label=s)
    axs[0].axhline(90,color='grey',ls=':',lw=1)
    axs[0].set(xlabel='Drive section xi',ylabel='Oriented mutual angle [deg]',title='Channel response vectors: almost perpendicular')
    axs[1].set(xlabel='Drive section xi',ylabel='cos(mutual angle)',title='Normalized dot product: departure from orthogonality')
    for ax in axs: ax.grid(alpha=.2);ax.legend()
    fig.tight_layout();fig.savefig(out/'channel_vector_angle_and_dot.png',dpi=170);plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(13,5.3))
    for s in SEEDS:
        axs[0].plot(xi,arrays[s+'_oriented_area'],label=s)
        axs[1].plot(xi,arrays[s+'_normalized_oriented_area'],label=s)
    axs[1].axhline(1,ls=':',c='grey',lw=1)
    axs[0].set(xlabel='Drive section xi',ylabel='det(response_v2, response_omega2)',title='Signed oriented area in the second Cartesian plane')
    axs[1].set(xlabel='Drive section xi',ylabel='sin(mutual angle)',title='Normalized orientation: rotationally invariant')
    for ax in axs:ax.grid(alpha=.2);ax.legend()
    fig.tight_layout();fig.savefig(out/'channel_oriented_area_and_handedness.png',dpi=170);plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(13.5,5.7))
    for s in SEEDS:
        for key,label,ls in [('critical_omega2_observer_degrees','omega2 zero','-'),
                             ('critical_v2_observer_degrees','v2 zero','--'),
                             ('critical_total_observer_degrees','total zero',':')]:
            axs[0].plot(xi,arrays[s+'_'+key],label=s+' '+label,ls=ls)
    for s in SEEDS:
        axs[1].plot(arrays[s+'_compensation_scan_degrees'],100*arrays[s+'_compensation_scan_fraction'],label=s)
    axs[0].set(xlabel='Drive section xi',ylabel='Critical observer rotation ψ relative to natural radial [deg]',
               title='Projection-zero angles (axes modulo 180 degrees)')
    axs[1].set(xlabel='Observer rotation ψ at xi=0.20 [deg]',ylabel='Radial projected cancellation [%]',
               title='Observer-dependent compensation despite nonzero physical response')
    axs[0].legend(fontsize=7,ncol=2);axs[1].legend(fontsize=8)
    for ax in axs:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'observer_critical_angles_and_compensation.png',dpi=170);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('preflight','full'),default='full')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    V,T,prov=read_verified_parent(full=(args.mode=='full'))
    if args.mode=='preflight':
        (OUT/'LOCAL_PREFLIGHT.json').write_text(json.dumps({'version':'v18.51w','parent':prov,'passed':True},indent=2)+'\n')
        print('Q550 v18.51w provenance preflight PASS');return
    xi=V['xi']
    require(xi.shape == (513,) and np.all(np.diff(xi)>0), 'V xi grid invalid')
    ledger=[]; arrays={'xi':xi}
    summaries={}
    for seed in SEEDS:
        print('Q550 v18.51w geometric channel audit:',seed,flush=True)
        summaries[seed]=measure_seed(seed,xi,V,T,ledger,arrays)
        print('  angle terminal deg:',summaries[seed]['angle_end_deg'],flush=True)
        print('  omega2 critical axis deg:',summaries[seed]['critical_observer_omega2_end_deg'],flush=True)
    checks={
      'parent_v_authentic_sources_outputs_and_crc':True,
      'both_seed_513_section_clocks_aligned':len(ledger)==2*513,
      'source_channel_vector_additivity':all(d['max_source_additivity_defect']<1e-11 for d in summaries.values()),
      'v_lab_matches_independent_t_cartesian':all(d['max_lab_vs_independent_t_cartesian']<1e-9 for d in summaries.values()),
      'dot_rotation_invariant':all(d['max_invariant_dot_native_vs_lab']<1e-12 and all(x['max_dot_invariance_defect']<1e-12 for x in d['observer_frame_checks'].values()) for d in summaries.values()),
      'oriented_area_rotation_invariant':all(d['max_invariant_area_native_vs_lab']<1e-12 and all(x['max_oriented_area_invariance_defect']<1e-12 for x in d['observer_frame_checks'].values()) for d in summaries.values()),
      'handedness_reflection_reverses_area':all(d['max_mirror_oriented_area_sign_defect']<1e-12 for d in summaries.values()),
      'angle_dot_area_unit_identity':all(d['max_unit_circle_identity_defect']<1e-12 for d in summaries.values()),
      'critical_axis_exact_projection_zero':all(max(d['max_v2_projection_zero_error'],d['max_omega2_projection_zero_error'],d['max_total_projection_zero_error'])<1e-12 for d in summaries.values()),
      'critical_axis_covariant_under_rotating_observers':all(all(x['max_critical_axis_covariance_sine_defect']<1e-12 for x in d['observer_frame_checks'].values()) for d in summaries.values()),
      'observer_channel_vector_reconstruction':all(all(x['max_channel_vector_reconstruction_defect']<1e-12 for x in d['observer_frame_checks'].values()) for d in summaries.values()),
      'mutual_orientation_consistent_between_seeds':abs(summaries['original']['angle_end_deg']-summaries['autonomous']['angle_end_deg'])<.001,
      'second_plane_orientation_positive_in_sample':all(d['angle_all_positive_orientation'] for d in summaries.values()),
      'near_origin_angle_undefined_guarded':all(d['n_valid_angles']==512 for d in summaries.values()),
      'physical_projection_not_confused_with_zero_material_response':all(d['total_end_norm']>.04 and abs(d['critical_observer_total_end_deg'])>1 for d in summaries.values()),
      'geometric_observer_conditioning_not_promoted_to_material_memory':True,
    }
    report={'version':'v18.51w','title':'Channel-Vector Angle, Oriented-Area & Observer-Critical Compensation Audit',
       'classification':'COORDINATE_SPLIT_CHANNEL_VECTOR_GEOMETRY_AND_OBSERVER_THRESHOLD_DIAGNOSTIC_NOT_MATERIAL_MEMORY_CERTIFICATE',
       'all_passed':bool(all(checks.values())), 'scientific_checks':checks,
       'sections_per_seed':len(xi),'source_channels':list(CHANNELS),'parent':prov,
       'seeds':summaries,'limitations':[
          'An observer projection can vanish while the Cartesian response vector remains nonzero; projected cancellation is not dynamical closure.',
          'The source-channel splitting is attached to chosen model coordinates. Angle, dot product and signed area are SO(2) observer-frame invariants only after that splitting is fixed.',
          'Signed area reverses under an orientation-reversing reflection; the apparent sign is not O(2) invariant.',
          'Critical observer angles are relative to the local radial frame, modulo pi; they do not define a preferred physical observer.',
          'At xi=0 the response vectors vanish: mutual angle and critical direction are undefined; these sections are explicitly masked.',
          'These diagnostics use existing verified v18.51v server trajectories and do not constitute new nonlinear or variational integrations.',
          'All findings are numerical on xi in [0,0.20], not a rigorous interval proof, material memory certificate, or return-section closure.',
       ]}
    with (OUT/'channel_vector_section_ledger.csv').open('w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=list(ledger[0]));writer.writeheader();writer.writerows(ledger)
    np.savez_compressed(OUT/'channel_vector_geometric_invariants.npz',**arrays)
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# Q550 v18.51w — Server Scientific Report','','## Scope',
           'Two seed families, 513 archived sections each. Signed Green-channel vectors in the physical second Cartesian plane.',
           'No nonlinear/variational IVPs were rerun. All values are numerical diagnostics.','','## Audit status',
           f"- Checks passed: {sum(checks.values())}/{len(checks)}", f"- All passed: {report['all_passed']}",
           '- Upstream parent: authenticated v18.51v server sendback including v18.51u/t/s/r genealogy.','',
           '## Terminal metrics']
    for seed in SEEDS:
        d=summaries[seed]
        lines += [f'### {seed}',
          f"- Channel-vector oriented angle (degrees): {d['angle_end_deg']:.12f}",
          f"- Signed oriented area: {d['oriented_area_end']:.12g}",
          f"- Dot product: {d['dot_end']:.12g}",
          f"- Critical omega2 observer angle (degrees): {d['critical_observer_omega2_end_deg']:.12f}",
          f"- Critical total projection angle (degrees): {d['critical_observer_total_end_deg']:.12f}",
          f"- Natural radial cancellation (%): {d['natural_radial_cancellation_end_percent']:.9f}",
          f"- Rotating sine cancellation (%): {d['sine_radial_cancellation_end_percent']:.9f}",
          f"- Rotating crest cancellation (%): {d['crest_radial_cancellation_end_percent']:.9f}", '']
    lines += ['## Interpretation',
         'The channel vectors are nearly perpendicular but not exactly so. A vanishing observer projection does not imply zero source channel response.',
         'Under proper passive plane rotations, mutual angles, norms, dot products and signed areas remain identical. Under reflection, the signed oriented area reverses.',
         'The critical observer direction is a coordinate-dependent projection threshold, not a material bifurcation or constitutive memory law.', '']
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
    render_figures(xi, arrays)
    print(f"Q550 v18.51w: {sum(checks.values())}/{len(checks)} checks PASS",flush=True)
    require(report['all_passed'], 'Scientific checks failed; see results/report.json')


if __name__ == '__main__': main()
