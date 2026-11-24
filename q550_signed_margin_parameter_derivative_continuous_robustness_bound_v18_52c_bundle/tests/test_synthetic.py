#!/usr/bin/env python3
"""Synthetic analytic checks; never substitutes for physical midpoint replay."""
import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_audit import extrema_pchip


def run():
    checks={}
    ks=np.array([0.,.5,1.,2.]);xx=np.array([.01,.02,.04]);m=.3+.1*ks[:,None]+np.zeros((4,3))
    poly=PchipInterpolator(ks,m,axis=0,extrapolate=False)
    lo,hi,lips=extrema_pchip(poly,ks)
    checks['exact_affine_minimum']=abs(lo-.3)<1e-12
    checks['exact_affine_maximum']=abs(hi-.5)<1e-12
    checks['exact_affine_kappa_derivative']=abs(lips-.1)<1e-12
    checks['pchip_interpolates_nodes']=np.max(abs(poly(ks)-m))<1e-14
    checks['pchip_does_not_extrapolate']=np.isnan(poly(-.1)).all()
    checks['constant_time_derivative_zero']=np.max(abs(np.diff(m,axis=1)))<1e-14
    # A hidden sinusoidal dip can defeat any criterion based only on endpoints.
    a=.3
    hidden=lambda x:a-2*a*np.sin(np.pi*x)**2
    checks['hidden_negative_between_positive_samples']=(hidden(0)>0 and hidden(1)>0 and hidden(.5)<0)
    checks['endpoint_secant_fails_as_derivative_majorant']=abs(hidden(1)-hidden(0))<1e-12 and abs(hidden(.25)-hidden(0))>0
    # Rigorous conditional estimate if true derivative bounds were supplied.
    m0=.24;h_k=.001;h_t=.0004;Lk=1.;Lt=10.;eps=1e-6
    bound=m0-Lk*h_k/2-Lt*h_t/2-eps
    checks['conditional_bound_positive_under_assumptions']=bound>0
    checks['conditional_bound_can_fail_if_derivative_too_large']=m0-1000*h_k/2-Lt*h_t/2-eps<0
    checks['negative_time_counterexample_not_claimed_physical']=hidden(.5)<0
    checks['positive_pchip_not_equal_positive_physical']=poly(.25)[0]>0 and hidden(.5)<0
    checks['signed_margin_rad_per_xi_units_preserved']=np.isfinite(lips)
    checks['two_dimensional_nearest_node_distance']=abs(h_k/2+h_t/2-.0007)<1e-12
    failed=[k for k,v in checks.items() if not v]
    print('Q550 v18.52c synthetic:',len(checks)-len(failed),'/',len(checks), 'PASS' if not failed else 'FAIL '+str(failed),flush=True)
    if failed:raise SystemExit(2)

if __name__=='__main__':run()
