#!/usr/bin/env python3
"""Check full server-produced two-time outputs, no model integrations."""
import csv,json
from pathlib import Path
import numpy as np
r=Path(__file__).resolve().parents[1]/'results'
a=json.loads((r/'report.json').read_text())
assert a['version']=='v18.51u' and a['mode']=='full', 'Server full audit missing'
assert a['all_passed'] and all(a['scientific_checks'].values())
with (r/'compensation_birth_section_ledger.csv').open() as f:
    ledger=list(csv.DictReader(f))
with (r/'causal_source_window_contributions.csv').open() as f:
    windows=list(csv.DictReader(f))
assert len(ledger)==1026 and len(windows)>0
with np.load(r/'green_causal_two_time_kernels.npz',allow_pickle=False) as z:
    assert z['xi'].shape==(513,)
    for seed in ('original','autonomous'):
        K=z[seed+'_kernel'];I=z[seed+'_integrated_radial_azimuthal']
        C=z[seed+'_cancellation_fraction']
        assert K.shape==(513,513,2,4) and I.shape==(513,2,4)
        assert np.all(np.isfinite(K)) and np.all(np.isfinite(I))
        assert np.all(K[np.triu_indices(513,k=1)]==0)
        assert np.max(np.abs(I[:,:,:3].sum(axis=2)-I[:,:,3]))<5e-10
        assert abs(C[-1]-a['seeds'][seed]['end_cancellation_fraction_radial'])<1e-12
        # All 8 exact Simpson source windows partition the terminal integral.
        for p,obs in enumerate(('radial','azimuthal')):
            for ch,c in enumerate(('v2','omega2','other','total')):
                rows=[w for w in windows if w['seed_family']==seed
                     and abs(float(w['xi_observation'])-.2)<1e-13
                     and w['projection']==obs and w['source_channel']==c]
                assert len(rows)==8
                assert abs(sum(float(w['signed_window_response']) for w in rows)-I[-1,p,ch])<1e-10
for seed in ('original','autonomous'):
    assert (r/f'green_kernel_causal_map_{seed}.png').stat().st_size>1000
assert (r/'compensation_birth_and_channel_history.png').stat().st_size>1000
assert (r/'SCIENTIFIC_REPORT_SERVER.md').stat().st_size>100
print('Q550 v18.51u full-result verification PASS')
