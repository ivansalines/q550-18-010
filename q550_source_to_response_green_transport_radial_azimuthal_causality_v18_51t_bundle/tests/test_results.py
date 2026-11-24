#!/usr/bin/env python3
"""Post-run validation of Q550 v18.51t results, independently from main."""
import csv,json
from pathlib import Path
import numpy as np
r=Path(__file__).resolve().parents[1]/'results'
a=json.loads((r/'report.json').read_text())
assert a['version']=='v18.51t' and a['mode']=='full'
assert a['all_passed'] and all(a['scientific_checks'].values())
assert len(a['green_channel_results'])==2
with (r/'green_transport_attribution_ledger.csv').open() as f: rows=list(csv.DictReader(f))
assert len(rows)==1026
with np.load(r/'green_channel_transport_states.npz') as f:
    assert len(f['xi'])==513
    for seed in ('original','autonomous'):
        x=f[seed+'_forced_Y'];q=f[seed+'_cartesian_Q'];g=f[seed+'_green_Y']
        assert x.shape==(513,10,4) and q.shape==(513,5,4) and g.shape==(513,10,4)
        assert np.all(np.isfinite(x)) and np.all(np.isfinite(g))
        assert np.linalg.norm(x[0])<1e-14
        assert np.max(np.linalg.norm(x[:,:,:3].sum(axis=2)-x[:,:,3],axis=1))<1e-8
        assert np.max(np.linalg.norm(q[:,:,:3].sum(axis=2)-q[:,:,3],axis=1))<1e-8
        m=a['green_channel_results'][seed]
        assert abs(np.linalg.norm(q[-1,:,3])-m['S_Q_total_norm_end'])<1e-12
        assert abs(m['S_Q_radial2_from_v2_end']+m['S_Q_radial2_from_omega2_end'] + float(f[seed+'_radial2'][-1,2]) - m['S_Q_radial2_total_end'])<1e-10
        assert abs(m['S_Q_azimuthal2_from_v2_end']+m['S_Q_azimuthal2_from_omega2_end'] + float(f[seed+'_azimuthal2'][-1,2]) - m['S_Q_azimuthal2_total_end'])<1e-10
assert (r/'source_to_response_green_radial_azimuthal_transport.png').stat().st_size>1000
print('Q550 v18.51t result verification PASS')
