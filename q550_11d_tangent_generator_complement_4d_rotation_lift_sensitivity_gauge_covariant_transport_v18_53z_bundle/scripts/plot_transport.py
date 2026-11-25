#!/usr/bin/env python3
"""Scientific plots from authenticated-model transport rows."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plots(rows,out):
    for seed in ('original','autonomous'):
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
        for epoch,linestyle in [('preformal','-'),('full_formal','--')]:
            for lift,marker in [('fixed_cartesian_velocity','o'),('fixed_generalized_velocity','s')]:
                rs=sorted([r for r in rows if r['seed']==seed and r['kappa_rule']=='shared' and r['epoch']==epoch
                           and r['temporal_lift']=='midpoint' and r['horizontal_lift']==lift],key=lambda r:r['lambda_value'])
                L=np.array([r['lambda_value'] for r in rs]);delta=np.array([r['model_vs_minimal_same_lambda_polygon_angle_deg'] for r in rs]);arc=np.array([r['model_minus_minimal_arc_rad'] for r in rs])
                label=epoch+' / '+lift.replace('fixed_','')
                axes[0].plot(L,delta,linestyle,marker=marker,markevery=4,markersize=3,label=label)
                axes[1].plot(L,arc,linestyle,marker=marker,markevery=4,markersize=3,label=label)
        axes[0].set(ylabel='Same-path polygon difference (degrees)',xlabel='Metric lambda',title='Complement-induced orientation difference')
        axes[1].set(ylabel='Additional covariant arclength (radians)',xlabel='Metric lambda',title='Path-dependent angular transport')
        for a in axes:a.grid(alpha=.3);a.legend(fontsize=7)
        fig.suptitle('Q550 v18.53z · seed '+seed+' · midpoint temporal lift')
        fig.savefig(out/('model_tangent_complement_'+seed+'.png'),dpi=140)
        plt.close(fig)
