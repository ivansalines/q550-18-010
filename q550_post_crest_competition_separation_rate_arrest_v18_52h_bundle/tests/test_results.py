#!/usr/bin/env python3
"""Post-run independent QA. No requirement for a physically positive finding."""
import json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
r=root/'results'
p=json.loads((r/'report.json').read_text())
assert p['version']=='v18.52h'
assert p['stage']=='SERVER_AUTHENTIC_POST_CREST_10D_CONTINUATION'
assert p['all_numerical_checks_passed'],p['numerical_checks']
assert len(p['results'])==8
with np.load(r/'post_crest_traces.npz',allow_pickle=False) as a:
    for label,rec in p['results'].items():
        t=a[label+'__tau'];D=a[label+'__D'];G=a[label+'__G'];C=a[label+'__C']
        P=a[label+'__P'];I=a[label+'__I']
        assert len(t)==p['settings']['samples']
        assert (np.diff(t)>0).all() and np.isfinite(t).all()
        assert np.max(abs(C-(P-I)))<1e-14
        assert D[0]>0 and np.isfinite(D).all()
        assert rec['tau_end']>rec['tau_start']>1.3
        ev=p['events'][label]
        assert rec['n_D_maxima']==sum(e['crossing']=='positive_to_negative' for e in ev['rate_zero'])
        assert rec['n_C_recoveries']==sum(e['crossing']=='negative_to_positive' for e in ev['convexity_zero'])
        assert rec['n_G_zeros']==len(ev['rate_zero'])
with np.load(r/'post_crest_continued_10d_states.npz',allow_pickle=False) as a:
    assert np.isclose(a['common_xi'][0],.8) and np.isclose(a['common_xi'][-1],p['settings']['xi_end'])
    for label,rec in p['results'].items():
        assert a[label+'__kick__Y_tau'].shape==(p['settings']['samples'],11)
assert (r/'post_crest_events.csv').exists()
assert (r/'post_crest_competition_rate_arrest.png').exists()
print('RESULTS PASS: 8 continuations, independent checks, event consistency')
