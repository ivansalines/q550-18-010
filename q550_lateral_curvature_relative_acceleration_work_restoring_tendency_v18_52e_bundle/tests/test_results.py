#!/usr/bin/env python3
"""Check numerical consistency without requiring a preferred scientific sign."""
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'results/report.json'
if not p.is_file():raise SystemExit('Missing server report.json')
r=json.loads(p.read_text())
assert r['version']=='v18.52e'
assert r['stage']=='SERVER_AUTHENTIC_NONLINEAR_SIGNED_VECTOR_DIAGNOSTIC'
assert r['all_numerical_checks_passed'] and all(r['checks'].values())
assert len(r['physical_kicks'])==8
assert len({x['seed'] for x in r['physical_kicks']})==2
for x in r['physical_kicks']:
    assert 0<=x['negative_signed_relative_power_fraction']<=1
    assert 0<=x['positive_signed_position_separation_fraction']<=1
    assert x['max_abs_norm_trace_replay_difference']<2e-7
assert not r['asymptotic_stability_or_dissipation_certified']
print('Q550 v18.52e server report consistency PASS')
