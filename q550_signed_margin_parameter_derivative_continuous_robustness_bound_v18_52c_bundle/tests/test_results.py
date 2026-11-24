#!/usr/bin/env python3
"""Validate scientific honesty and output structure, not force a positive result."""
import json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=R/'results/report.json'
if not p.is_file():raise SystemExit('Missing generated physical audit report')
a=json.loads(p.read_text())
checks={
 'correct_version':a['version']=='v18.52c',
 'twenty_six_parent':a['historical_52b_branches']==26,
 'no_unconditional_certificate':a['physical_continuous_certificate'] is False,
 'explicit_derivative_limitation':any('rigorous' in v for v in a['limitations']),
 'distinct_scientific_decision':isinstance(a['scientific_decision'],str),
 'two_families':len(a['differentiation'])==2,
 'no_fake_completed_probes':a['new_physical_midpoints_completed']<=a['requested_new_physical_midpoints'],
 'real_report':p.stat().st_size>1000,
}
bad=[k for k,v in checks.items() if not v]
print('Q550 v18.52c result structure:',sum(checks.values()),'/',len(checks),flush=True)
if bad:raise SystemExit('Invalid result: '+str(bad))
