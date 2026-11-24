#!/usr/bin/env python3
"""Lightweight genuine v18.52a historic parameter-band visualization."""
import importlib.util,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('audit52b',root/'scripts/run_audit.py')
audit=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=audit;spec.loader.exec_module(audit)
_,anchor=audit.authenticate_parent(audit.find_input(audit.ANAME))
center,half,rows,series=audit.anchor_rows(anchor)
fig,ax=plt.subplots(1,2,figsize=(11.8,4.8),constrained_layout=True)
for seed in audit.SEEDS:
 data=sorted([r for r in rows if r['seed']==seed],key=lambda r:r['kappa'])
 xx=[(r['kappa']-center)/half for r in data]
 ax[0].plot(xx,[r['min_margin'] for r in data],marker='o',label=seed)
 ax[1].plot(xx,[r['min_total'] for r in data],marker='o',label=seed)
for aa,title in zip(ax,['Signed mixed + self minimum','Full turning minimum']):
 aa.set(xlabel='(κ₀ − historical midpoint) / historical half-width',ylabel='rad / ξ',title=title)
 aa.axhline(0,color='gray',lw=.7);aa.grid(alpha=.25);aa.legend()
fig.suptitle('Q550 v18.52a: authenticated six-branch anchor for continuation (not new radius)')
(root/'results').mkdir(exist_ok=True)
fig.savefig(root/'results/historical_52a_signed_margin_anchor.png',dpi=160)
plt.close(fig)
print('Q550 v18.52b preflight figure from six authenticated v18.52a branches',flush=True)
