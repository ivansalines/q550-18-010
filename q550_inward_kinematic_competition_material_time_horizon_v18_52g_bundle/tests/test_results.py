#!/usr/bin/env python3
"""Validate a completed v18.52g server run; preserve scientific non-findings."""
import json
from pathlib import Path
import numpy as np
r=Path(__file__).resolve().parents[1]/'results'
p=json.loads((r/'report.json').read_text())
assert p['version']=='v18.52g'
assert p['stage']=='SERVER_AUTHENTIC_NONLINEAR_EXTENSION'
assert len(p['results'])==8
assert p['all_numerical_checks_passed'],p['numerical_checks']
assert (r/'material_horizon_events.csv').is_file()
assert (r/'material_horizon_ledger.csv').is_file()
assert (r/'material_horizon_traces.npz').is_file()
assert (r/'extended_10d_nonlinear_states.npz').is_file()
with np.load(r/'material_horizon_traces.npz',allow_pickle=False) as traces:
    for label,res in p['results'].items():
        tau=traces[label+'__tau'];d=traces[label+'__D'];g=traces[label+'__G']
        c=traces[label+'__C'];pos=traces[label+'__P'];inward=traces[label+'__I']
        assert np.isfinite(tau).all() and np.all(np.diff(tau)>0)
        assert len(tau)==p['settings']['samples']
        assert np.max(abs(c-(pos-inward)))<1e-14
        assert d[0]<1e-19 and np.max(d)>d[0]
        assert res['tau_end']>0.33
        assert len(p['events'][label]['position_separation_maxima']) <= len(p['events'][label]['rate_zero'])
print('PASS: v18.52g authenticated-source server results and event/report consistency')
