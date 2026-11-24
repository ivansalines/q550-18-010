#!/usr/bin/env python3
"""Exact file hashes of all finished non-preflight scientific outputs."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[1]
results=root/'results'
expected=['report.json','SCIENTIFIC_REPORT_SERVER.md','paired_hinge_views_single_trajectory.png',
          'single_trajectory_reduced_3d_and_hinge_radius.png','paired_observer_gate_windows.png',
          'angle_and_hidden_direction_sensitivity.csv','tilted_hinge_section_crossings_5d.csv',
          'cross_view_sequential_posture_pairs_5d.csv','single_observer_consecutive_section_return_map_5d.csv',
          'single_trajectory_hinge_observer_metrics.csv']
for f in expected:
    if not (results/f).is_file():raise RuntimeError('Missing output '+f)
manifest={'version':'v18.52m','file_sha256':{
    'results/'+p.name:hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted(results.iterdir()) if p.is_file() and p.name!='PREFLIGHT.json'}}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('OUTPUT MANIFEST WRITTEN',len(manifest['file_sha256']),'files')
