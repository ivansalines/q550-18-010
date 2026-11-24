#!/usr/bin/env python3
"""Independent lightweight output sanity checks, run after server audit."""
from pathlib import Path
import json,sys,numpy as np
root=Path(__file__).resolve().parents[1]
report=json.loads((root/'results/report.json').read_text())
assert report['version']=='v18.51v' and report['mode']=='full', 'Incorrect output version'
assert report['all_passed'] and all(report['scientific_checks'].values()),'Scientific checks failed'
with np.load(root/'results/physical_clock_frame_green_tensors.npz',allow_pickle=False) as a:
    xi=a['xi'];assert len(xi)==513 and abs(xi[-1]-.2)<1e-14
    for s in ('original','autonomous'):
        c=a[s+'_clock'];tau=a[s+'_tau']
        assert np.all(c>0) and np.all(np.diff(tau)>0)
        K=a[s+'_kernel_physical_time_density']
        assert K.shape==(513,513,2,4)
        assert np.max(abs(K[np.triu_indices(513,1)]))==0
        I=a[s+'_lab_frame_response']
        W=a[s+'_rotating_sine_response']
        assert I.shape==(513,2,4) and W.shape==(513,2,4)
        assert np.max(abs(I[:,:,:3].sum(axis=2)-I[:,:,3]))<1e-8
        assert np.max(abs(W[:,:,:3].sum(axis=2)-W[:,:,3]))<1e-8
        sc=a[s+'_clock_shift_cartesian'];st=a[s+'_fixed_tau_matched_cartesian_sensitivity']
        assert sc.shape==st.shape==(513,5)
        assert np.max(abs(sc[0]))<1e-12
for name in ('physical_clock_vs_fixed_tau_sensitivity.png',
             'moving_frame_components_vs_physical_vectors.png',
             'physical_clock_causal_green_density.png',
             'SCIENTIFIC_REPORT_SERVER.md','physical_clock_frame_section_ledger.csv'):
    assert (root/'results'/name).exists(),name
print('Q550 v18.51v independent server output tests PASS')
