#!/usr/bin/env python3
"""Write output hashes, excluding the manifest itself."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1]
names=['results/report.json','results/physical_parameter_robustness_sections.csv',
  'results/physical_branch_turning_arrays.npz',
  'results/physical_kappa_branch_cross_channel_turning_robustness.png',
  'results/SCIENTIFIC_REPORT_SERVER.md']
m={'root':R.name,'version':'v18.52a','file_sha256':{
    n:hashlib.sha256((R/n).read_bytes()).hexdigest() for n in names}}
(R/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(m,indent=2)+'\n')
print('Q550 v18.52a wrote',len(names),'output digests')
