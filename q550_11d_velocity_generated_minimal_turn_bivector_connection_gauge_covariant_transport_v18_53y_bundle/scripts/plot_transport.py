"""Lightweight overview plots for v18.53y, using only analysis CSV rows."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make_plots(rows,out):
    for seed in ('original','autonomous'):
        fig,ax=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
        for epoch in ('preformal','full_formal'):
            for lift in ('start','midpoint','end'):
                tab=sorted((r for r in rows if r['seed']==seed and r['kappa_rule']=='shared' and r['epoch']==epoch and r['lift']==lift),key=lambda r:r['lambda_value'])
                ax[0].plot([r['lambda_value'] for r in tab],[r['dynamic_frame_angle_deg'] for r in tab],label=f'{epoch}/{lift}')
                ax[1].plot([r['lambda_value'] for r in tab],[r['dynamic_minus_canonical_arc_rad'] for r in tab],label=f'{epoch}/{lift}')
        for a in ax:a.set_xlabel('Metric strength lambda');a.grid(alpha=.22);a.legend(fontsize=7)
        ax[0].set_ylabel('Angle relative to lambda=0 [degrees]');ax[0].set_title('11D velocity-induced frame: oriented polygon')
        ax[1].set_ylabel('Covariant minus fixed arc [rad]');ax[1].set_title('Different temporal lifts of 10D integrated residual')
        fig.suptitle('Q550 v18.53y — Kinematic SO(5) minimum-turn induced 10D geometry')
        fig.savefig(out/f'dynamic_connection_oriented_transport_{seed}.png',dpi=140)
        plt.close(fig)
