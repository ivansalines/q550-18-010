#!/usr/bin/env python3
"""Check generated audit files without elevating signed margin to a guaranteed theorem."""
import csv,json,sys
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
out=root/'results'
r=json.loads((out/'report.json').read_text())
with (out/'physical_parameter_continuation_scan.csv').open(newline='') as f:
 rows=list(csv.DictReader(f))
checks={
 'right_version':r.get('version')=='v18.52b',
 'both_seeds_in_rows':set(row['seed'] for row in rows)=={'original','autonomous'},
 'six_historical_anchors':sum(row['kind']=='authenticated_parent_branch' for row in rows)==6,
 'parent_sha_present':len(r['sources']['v18_52a_parent_sha256'])==64,
 'signed_margin_finite_all_reported':all(np.isfinite(float(row['min_margin'])) for row in rows),
 'scientific_decision_separate_from_numeric':r['scientific_decision']!='' and 'numerical_audit_status' in r,
 'no_spurious_interval_claim':'SAMPLED' in r['classification'],
 'scan_figure_exists':(out/'physical_kappa_signed_margin_continuation.png').is_file(),
 'scan_npz_exists':(out/'physical_parameter_continuation_curves.npz').is_file(),
 'boundary_directions':len(r['boundary_search'])==4,
 'at_least_six_anchors':len(rows)>=6,
 'branch_count_consistent':len(rows)==r['total_branch_records'],
}
print('Q550 v18.52b output validation',sum(checks.values()),'/',len(checks),flush=True)
for k,v in checks.items():
 if not v:print('FAIL',k)
if not all(checks.values()):sys.exit(2)
