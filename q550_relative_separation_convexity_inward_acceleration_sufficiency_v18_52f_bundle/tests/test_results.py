import json
from pathlib import Path
r=Path(__file__).resolve().parents[1]/'results'
d=json.loads((r/'report.json').read_text())
assert d['version']=='v18.52f'
assert d['check_count']==d['checks_passed']==48
assert d['all_numerical_checks_passed']
assert len(d['physical_kicks'])==8
assert d['global_summary']['samples_total']==8*513
assert not d['global_summary']['certified_continuous_bound']
for s in d['physical_kicks']:
 assert s['ratio_final']<1 and s['final_sufficiency_factor_shortfall']>1
 assert s['first_sampled_convexity_reversal_tau'] is None
 assert s['first_sampled_separation_rate_reversal_tau'] is None
 assert s['sampled_convexity_all_positive']
 assert s['sampled_position_rate_all_positive_after_initial']
print('Q550 v18.52f result test 8/8 PASS')
