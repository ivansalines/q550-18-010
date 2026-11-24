#!/usr/bin/env python3
"""Post-run output checks, no reinterpretation of closure."""
import json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
r=json.loads((root/'results/report.json').read_text())
assert r['version']=='v18.52j' and r['stage']=='SERVER_AUTHENTIC_POST_I_10D_CONTINUATION'
assert r['all_numerical_checks_passed'] and r['refinement_performed']
assert len(r['results'])==8 and len(r['checks'])>=48 and all(r['checks'].values())
with np.load(root/'results/long_horizon_10d_states.npz',allow_pickle=False) as f:
    assert abs(float(f['common_xi'][0])-3.2)<1e-13
    assert abs(float(f['common_xi'][-1])-r['settings']['xi_end'])<1e-13
    for key in f.files:
        assert np.isfinite(f[key]).all(),key
        if key!='common_xi':assert f[key].shape[1]==11,key
with np.load(root/'results/long_horizon_5d_vectors_and_traces.npz',allow_pickle=False) as f:
    for label in r['results']:
        for suffix in ('dQ','dV','da'):
            assert f[label+'__'+suffix].shape[1]==5
        for suffix in ('angle_from_h_seam_deg','angle_from_seam_deg'):
            assert f[label+'__'+suffix].shape==f[label+'__D'].shape
        assert np.isfinite(f[label+'__D']).all()
print('SERVER RESULTS TEST PASS',len(r['results']),'kicks')
